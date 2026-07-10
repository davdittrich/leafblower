"""benchmarks/study/common/run_matrix.py -- WU-8 (leafblower-2ouc.9).

Pure-logic helpers shared by python/run_arm.py (and mirrored by
R/run_arm.R's common/run_matrix.R) for the registry-driven run-matrix
driver: applicability gating, run-cell enumeration + randomization,
hardware-isolation capture, the pre-run frozen-tag SHA-check hook, and
timing statistics (median/percentile/bootstrap CI).

This module intentionally does NOT import any adapter, pandas, numpy, or
leafblower -- it is pure stdlib so both the driver's orchestrator and its
unit tests can exercise gating/environment/timing logic without paying
subprocess or package-import cost (DESIGN.md Sec5/Sec8).

STRICT SEPARATION (user constraint 2026-07-08): imports nothing from
leafblower's own src/, r_bridge.cpp, R/, python/leafblower/. No registry
"applicability pin" field (arm/families/bounds/K_max/builds) is ever
mutated here -- only `applicable_problems`, which contract.md/registry_
schema.json document as the WU-10-installed-version-gated, cross-checkable
derived field this WU is explicitly permitted to extend/validate.
"""

from __future__ import annotations

import json
import math
import random
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional

_THIS_DIR = Path(__file__).resolve().parent              # benchmarks/study/common
STUDY_DIR = _THIS_DIR.parent                               # benchmarks/study
REPO_ROOT = STUDY_DIR.parent.parent                          # repo root

FROZEN_TAG = "benchmark-config-freeze-v9"

# Paths (relative to REPO_ROOT) whose git tree must match the frozen tag's
# tree before any timed cell runs. registry.json is DELIBERATELY excluded --
# contract.md/this WU's own ticket authorize extending/validating its
# `applicable_problems` field, so it is not part of the frozen perimeter.
FROZEN_PATHS = [
    "benchmarks/study/spec",
    "benchmarks/study/common/problem_io.R",
    "benchmarks/study/common/problem_io.py",
    "benchmarks/study/common/metrics.R",
    "benchmarks/study/common/metrics.py",
    "benchmarks/study/common/instance_family.R",
    "benchmarks/study/common/instance_family.py",
    "benchmarks/study/R/competitors.R",
    "benchmarks/study/R/leafblower_adapter.R",
    "benchmarks/study/python/competitors.py",
    "benchmarks/study/python/leafblower_adapter.py",
]

THREAD_SWEEP = (1, 4)

# ott-jax is excluded from ranked {1,4}-thread timing (DESIGN.md Sec2/Sec5,
# WU-OTT leafblower-2ouc.18): XLA ignores OMP_/OPENBLAS_/MKL_NUM_THREADS, so
# a thread-sweep row for it would misrepresent the isolation contract.
RANKED_TIMING_EXCLUDED = {"ott_jax_sinkhorn"}


# ---------------------------------------------------------------------------
# Registry / spec loading
# ---------------------------------------------------------------------------

def load_registry(path: Path | str = STUDY_DIR / "registry.json") -> dict[str, Any]:
    with open(path) as f:
        return json.load(f)


def load_problem_specs(spec_dir: Path | str = STUDY_DIR / "spec") -> dict[str, dict[str, Any]]:
    """Raw spec/*.json dicts keyed by `id` (schema.json shape) -- NOT the
    resolved problem_io object. Only the scalar fields (objective_families,
    bounds, K) needed for applicability gating are used; bulk `data`/
    `data_ref` are left untouched (never loaded here)."""
    out: dict[str, dict[str, Any]] = {}
    for p in sorted(Path(spec_dir).glob("*.json")):
        if p.name in ("schema.json", "runs_schema.json", "registry_schema.json",
                      "status_enum.json", "hyperparams.json", "tol_mapping.json",
                      "instance_family.json"):  # WU-3 manifest (not a spec; per-instance specs live in instance_family/)
            continue
        with open(p) as f:
            spec = json.load(f)
        out[spec["id"]] = spec
    return out


def load_instance_family_specs(
    spec_dir: Path | str = STUDY_DIR / "spec", regenerate: bool = True,
) -> dict[str, dict[str, Any]]:
    """Materializes (WU-3 instance_family.generate_instance_family_specs) and
    loads the 32-instance synthetic family from spec/instance_family/<id>.json.
    Lazy-imports instance_family (which pulls pandas) so this module's top-level
    import surface stays pure-stdlib for the unit tests / static-only callers.
    Generation is cheap regardless of n -- the gen: data is resolved lazily at
    solve time, never here."""
    import instance_family  # lazy: keep run_matrix's module import surface stdlib-only
    spec_dir = Path(spec_dir)
    inst_dir = spec_dir / "instance_family"
    if regenerate or not inst_dir.exists():
        instance_family.generate_instance_family_specs(spec_dir)
    out: dict[str, dict[str, Any]] = {}
    for p in sorted(inst_dir.glob("*.json")):
        with open(p) as f:
            spec = json.load(f)
        out[spec["id"]] = spec
    return out


def load_all_problem_specs(
    spec_dir: Path | str = STUDY_DIR / "spec",
    include_instance_family: bool = True,
    regenerate: bool = True,
) -> dict[str, dict[str, Any]]:
    """The MATRIX's full problem set: the 4 canonical static specs UNION the
    WU-3 instance family (32 instances). load_problem_specs() is retained,
    static-only, for callers/tests that must see ONLY the canonical 4 -- the
    instance family is the sole arena declaring `ot`/`newton_kl`
    objective_families (giving the OT / Newton-KL solvers their applicable
    problems), so folding it into the static loader would change what the
    static-only consumers observe."""
    specs = load_problem_specs(spec_dir)
    if include_instance_family:
        specs.update(load_instance_family_specs(spec_dir, regenerate=regenerate))
    return specs


def instance_family_n(spec: dict[str, Any]) -> Optional[int]:
    """Row count `n` parsed from a `gen:instance_family?n=...` data_ref; None
    for any spec that is NOT an instance-family spec (static file:/pkg:/inline).
    Lets --smoke bound the family to small instances (a huge-n synthetic
    instance would spend minutes in the pure-Python generator only to be
    subsampled away)."""
    data_ref = spec.get("data_ref", "") or ""
    if not data_ref.startswith("gen:instance_family"):
        return None
    query = data_ref.partition("?")[2]
    for part in query.split("&"):
        key, _, val = part.partition("=")
        if key == "n":
            try:
                return int(val)
            except ValueError:
                return None
    return None


def resolve_spec_path(problem_id: str, spec_dir: Path | str = STUDY_DIR / "spec") -> Path:
    """Locates a scheduled cell's spec JSON: canonical spec/<id>.json first,
    then the generated spec/instance_family/<id>.json. Raises FileNotFoundError
    if neither exists (a scheduled cell whose spec vanished is a hard error, not
    a silent skip)."""
    spec_dir = Path(spec_dir)
    direct = spec_dir / f"{problem_id}.json"
    if direct.exists():
        return direct
    nested = spec_dir / "instance_family" / f"{problem_id}.json"
    if nested.exists():
        return nested
    raise FileNotFoundError(
        f"no spec JSON for problem id {problem_id!r} under {spec_dir} "
        f"or {spec_dir / 'instance_family'}"
    )


def load_hyperparams(path: Path | str = STUDY_DIR / "spec" / "hyperparams.json") -> dict[str, Any]:
    with open(path) as f:
        return json.load(f)


def solver_package(solver_id: str, hyperparams: dict[str, Any]) -> str:
    """(language,package) baseline-RSS key. Falls back to solver_id for
    leafblower's own methods (package == "leafblower", shared across all 9)
    and for any competitor whose hyperparams.json entry omits `package`."""
    if solver_id.startswith("leafblower_"):
        return "leafblower"
    for table in ("R_competitors", "Python_competitors"):
        entry = hyperparams.get(table, {}).get(solver_id)
        if entry and entry.get("package"):
            return str(entry["package"])
    return solver_id


# ---------------------------------------------------------------------------
# Applicability gating (contract.md Sec4 / registry_schema.json)
# ---------------------------------------------------------------------------

def problem_bounds_kind(problem_bounds: dict[str, Any]) -> str:
    """'unbounded' iff max is null/None or +inf; else 'bounded'."""
    mx = problem_bounds.get("max")
    if mx is None:
        return "unbounded"
    try:
        if math.isinf(float(mx)):
            return "unbounded"
    except (TypeError, ValueError):
        pass
    return "bounded"


def is_applicable(entry: dict[str, Any], problem_spec: dict[str, Any]) -> bool:
    """True iff `entry` (a registry.json SolverEntry) may legally be
    scheduled against `problem_spec` (a raw spec/*.json dict): family
    intersection, bounds compatibility, K_max cap (registry_schema.json)."""
    if not (set(entry["families"]) & set(problem_spec["objective_families"])):
        return False
    kind = problem_bounds_kind(problem_spec["bounds"])
    if entry["bounds"] not in ("both", kind):
        return False
    k_max = entry.get("K_max")
    if k_max is not None and int(problem_spec["K"]) > int(k_max):
        return False
    data_ref = problem_spec.get("data_ref", "")
    if entry["arm"] == "python" and data_ref.startswith("pkg:"):
        # problem_io.py's _PKG_LOADERS is empty by design (DESIGN.md Sec3):
        # pkg: data_refs are R-package datasets (survey/sampling/anesrake/
        # ebal/optweight/balance), resolvable generically by problem_io.R's
        # requireNamespace()-based loader but NOT by Python (no pkg loader
        # registered). Gate python-arm out rather than schedule a cell that
        # is guaranteed to fail with NotImplementedError at solve time.
        return False
    return True


def compute_applicable_problems(
    registry: dict[str, Any], problem_specs: dict[str, dict[str, Any]],
    available: Optional[set[str]] = None, arm: Optional[str] = None,
) -> dict[str, Any]:
    """Returns a COPY of `registry` with every solver's `applicable_problems`
    recomputed from is_applicable() against `problem_specs`. Every other
    SolverEntry field (arm/families/bounds/K_max/builds -- the frozen
    "applicability pins") is passed through unchanged.

    `available`+`arm` (when given) gate out THIS-ARM registry solvers that do NOT
    resolve to a runnable adapter (e.g. `samplics`, installed but deliberately
    unwrapped -- svy is its successor): such a solver gets an empty
    applicable_problems so the matrix never schedules a cell that would KeyError
    at adapter-resolution time. Solvers of OTHER arms are NEVER gated by this
    arm's `available` set -- their applicable_problems are recomputed normally
    (is_applicable is arm-agnostic), so a single-arm sync never clobbers the
    other arm's applicability."""
    out = json.loads(json.dumps(registry))  # deep copy, dict order preserved
    for solver_id, entry in out["solvers"].items():
        if arm is not None and entry.get("arm") != arm:
            continue  # other arm: leave applicable_problems untouched (its own driver syncs it)
        if available is not None and solver_id not in available:
            entry["applicable_problems"] = []
            continue
        applicable = [
            pid for pid, spec in sorted(problem_specs.items())
            if is_applicable(entry, spec)
        ]
        entry["applicable_problems"] = applicable
    return out


# ---------------------------------------------------------------------------
# Run-cell enumeration
# ---------------------------------------------------------------------------

def build_matrix(
    registry: dict[str, Any],
    problem_specs: dict[str, dict[str, Any]],
    arm: str,
    threads: tuple[int, ...] = THREAD_SWEEP,
    rng_seed: Optional[int] = None,
    exclude_ranked_timing: bool = True,
    available: Optional[set[str]] = None,
) -> list[dict[str, Any]]:
    """One dict per (solver,problem,thread,build) cell gated by
    registry[solver]['arm']==arm and is_applicable(); order RANDOMIZED
    (DESIGN.md Sec5 hardware-isolation requirement). `rng_seed` makes the
    shuffle reproducible for tests; omit (None) for a real run (uses
    `random.SystemRandom`-quality default seeding via `random.Random()`)."""
    cells: list[dict[str, Any]] = []
    for solver_id, entry in registry["solvers"].items():
        if entry["arm"] != arm:
            continue
        if available is not None and solver_id not in available:
            continue  # registry solver with no runnable adapter this arm (e.g. samplics)
        if exclude_ranked_timing and solver_id in RANKED_TIMING_EXCLUDED:
            continue
        for pid in entry.get("applicable_problems") or []:
            spec = problem_specs.get(pid)
            if spec is None or not is_applicable(entry, spec):
                continue  # defensive re-check; registry may be stale
            for build in entry["builds"]:
                for thread in threads:
                    cells.append(dict(solver=solver_id, problem=pid,
                                       thread=thread, build=build))
    random.Random(rng_seed).shuffle(cells)
    return cells


def weights_store_filename(solver: str, problem: str, thread: int, build: str) -> str:
    """DESIGN.md Sec8 / contract.md Sec2.1 weights store key convention."""
    return f"{solver}__{problem}__t{thread}__{build}.parquet"


# ---------------------------------------------------------------------------
# Pre-run frozen-tag SHA-check hook (WU-9 wires enforcement; this WU commits
# to calling it before any timed cell executes -- present-but-graceful when
# the tag has not been re-cut since a legitimate post-freeze amendment).
# ---------------------------------------------------------------------------

def _git(args: list[str], repo_root: Path) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=repo_root, capture_output=True, text=True)


def assert_frozen_tag(
    tag: str = FROZEN_TAG,
    paths: Optional[list[str]] = None,
    repo_root: Path = REPO_ROOT,
) -> dict[str, Any]:
    """Halts (raises RuntimeError) if `tag` exists but HEAD has drifted from
    it under any FROZEN_PATHS entry -- CLAUDE.md Compliance Ultimatum
    Anti-Pivot rule (no silent workaround on a detected spec/adapter drift).
    Warns (does not raise) if the tag does not exist yet."""
    paths = paths if paths is not None else FROZEN_PATHS
    exists = _git(["rev-parse", "-q", "--verify", f"refs/tags/{tag}"], repo_root).returncode == 0
    if not exists:
        msg = (f"WARN: frozen tag {tag!r} not found -- freeze not yet cut "
               f"(WU-9 wires enforcement); proceeding unchecked.")
        print(msg, file=sys.stderr)
        return dict(tag_exists=False, ok=True, mismatches=[], message=msg)

    mismatches = []
    for p in paths:
        tag_obj = _git(["rev-parse", f"{tag}:{p}"], repo_root)
        head_obj = _git(["rev-parse", f"HEAD:{p}"], repo_root)
        if tag_obj.returncode != 0 or head_obj.returncode != 0:
            mismatches.append(p)  # path missing on one side -- also drift
            continue
        if tag_obj.stdout.strip() != head_obj.stdout.strip():
            mismatches.append(p)
    if mismatches:
        raise RuntimeError(
            f"SPEC_FAILURE: frozen tag {tag!r} exists but HEAD has drifted "
            f"from it under: {mismatches} -- halting per CLAUDE.md Anti-Pivot "
            f"rule (no workaround; re-tag+re-rehearse per contract.md Sec6)."
        )
    return dict(tag_exists=True, ok=True, mismatches=[])


def assert_runnable_frozen(tag: str, repo_root: Path = REPO_ROOT) -> dict[str, Any]:
    """WU-9 pre-run runnable-tree SHA gate (mirrors run_matrix.R). Broader than
    assert_frozen_tag: a SCORED/timed run must execute against the EXACT frozen
    runnable tree. Raises RuntimeError (SPEC_FAILURE, hard stop) if the
    runnable-tree tag is absent, the benchmarks/study working tree has any
    uncommitted/untracked non-ignored change (run outputs are gitignored, so
    excluded), or the benchmarks/study tree at HEAD has drifted from the tag.
    Call ONLY for scored runs (WU-11), NOT the rehearsal (unfrozen by design)."""
    if _git(["rev-parse", "-q", "--verify", f"refs/tags/{tag}"], repo_root).returncode != 0:
        raise RuntimeError(
            f"SPEC_FAILURE: runnable-tree freeze tag {tag!r} not found -- refusing "
            f"to execute a scored/timed run against an unfrozen tree (WU-9). Cut the signed tag first.")
    st = _git(["status", "--porcelain", "--", "benchmarks/study"], repo_root)
    dirty = [ln for ln in st.stdout.split("\n") if ln.strip()]
    drifted = _git(["diff", "--quiet", tag, "HEAD", "--", "benchmarks/study"], repo_root).returncode != 0
    if dirty or drifted:
        parts = ""
        if dirty:
            parts += f"; uncommitted under benchmarks/study: {' | '.join(dirty[:5])}"
        if drifted:
            parts += "; benchmarks/study tree at HEAD has DRIFTED from the tag"
        raise RuntimeError(
            f"SPEC_FAILURE: refusing to run a scored/timed matrix against a dirty/unfrozen "
            f"runnable tree (WU-9). tag={tag!r}{parts}. Commit + re-tag (amendment protocol, "
            f"contract.md Sec6) before running.")
    return dict(tag_exists=True, ok=True, clean=True, tag=tag)


# ---------------------------------------------------------------------------
# Hardware isolation capture (DESIGN.md Sec5 Blocker F)
# ---------------------------------------------------------------------------

def _read_text(path: Path) -> Optional[str]:
    try:
        return path.read_text().strip()
    except OSError:
        return None


def _cpu_governors() -> list[str]:
    vals = []
    for p in sorted(Path("/sys/devices/system/cpu").glob("cpu[0-9]*/cpufreq/scaling_governor")):
        v = _read_text(p)
        if v is not None:
            vals.append(v)
    return vals


def _turbo_boost_enabled() -> Optional[bool]:
    """True/False if determinable, else None (unknown platform). AMD/Arm
    generic cpufreq 'boost' file: '1' = boost ENABLED. Intel P-State
    'no_turbo': '1' = turbo DISABLED (inverted polarity)."""
    boost = _read_text(Path("/sys/devices/system/cpu/cpufreq/boost"))
    if boost is not None:
        return boost.strip() == "1"
    no_turbo = _read_text(Path("/sys/devices/system/cpu/intel_pstate/no_turbo"))
    if no_turbo is not None:
        return no_turbo.strip() == "0"
    return None


def capture_environment(pin_core: Optional[int] = None) -> dict[str, Any]:
    """Captures + asserts hardware-isolation state (governor/turbo/pinning)
    per DESIGN.md Sec5 Blocker F. Does NOT attempt to *set* governor/turbo
    (requires root) -- records the observed state and WARNS to stderr (does
    not raise) when unisolated, so the driver still runs but the gap is
    always visible in environment.json (matches WU-10's own env/environment.json
    'hardware_state_for_timed_runs' precedent)."""
    import shutil

    governors = _cpu_governors()
    governor = governors[0] if governors and len(set(governors)) == 1 else (
        "mixed" if governors else None
    )
    turbo = _turbo_boost_enabled()
    taskset = shutil.which("taskset") is not None
    numactl = shutil.which("numactl") is not None

    governor_ok = governor == "performance"
    turbo_ok = turbo is False
    pinning_active = pin_core is not None

    state = dict(
        cpu_governor_all_cores=governor,
        cpu_governor_required_for_timed_runs="performance",
        cpu_governor_conformant=governor_ok,
        turbo_boost_enabled=turbo,
        turbo_boost_required_for_timed_runs="disabled",
        turbo_boost_conformant=turbo_ok,
        core_pinning_active=pinning_active,
        core_pinning_pin_core=pin_core,
        core_pinning_mechanism_available=(taskset or numactl),
        core_pinning_mechanism=("taskset" if taskset else ("numactl" if numactl else None)),
    )
    isolated = governor_ok and turbo_ok and pinning_active
    if not isolated:
        gaps = [k for k, ok in (("governor", governor_ok), ("turbo", turbo_ok),
                                 ("pinning", pinning_active)) if not ok]
        print(f"WARN: hardware not fully isolated for timed runs (gaps: {gaps}); "
              f"see environment.json 'hardware_state_for_timed_runs' -- proceeding "
              f"(warn-not-fail per WU-8 scope; root-level fix is WU-11's job).",
              file=sys.stderr)

    import datetime
    return dict(
        captured_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        captured_by="WU-8 run_arm.py (leafblower-2ouc.9)",
        arm="python",
        hardware_state_for_timed_runs=state,
        isolated=isolated,
    )


# ---------------------------------------------------------------------------
# Timing statistics: median + 5/95 percentile + percentile-bootstrap CI
# ---------------------------------------------------------------------------

def _percentile(sorted_vals: list[float], q: float) -> float:
    """Linear-interpolation percentile (numpy 'linear' method equivalent),
    q in [0, 1]. `sorted_vals` must already be sorted ascending."""
    n = len(sorted_vals)
    if n == 1:
        return sorted_vals[0]
    pos = q * (n - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return sorted_vals[lo]
    frac = pos - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def summarize_timing(
    reps: list[float],
    n_boot: int = 2000,
    ci: float = 0.90,
    rng_seed: Optional[int] = None,
) -> dict[str, float]:
    """Median, 5th/95th percentile, and a percentile-bootstrap CI on the
    median (DESIGN.md Sec5). `reps` are the timed (post-warmup) wall_time_s
    values for one (solver,problem,thread,build) cell."""
    if not reps:
        raise ValueError("summarize_timing: reps must be non-empty")
    reps_sorted = sorted(reps)
    median = statistics.median(reps)
    if len(reps) == 1:
        return dict(median=median, p05=reps[0], p95=reps[0],
                     boot_ci_lo=reps[0], boot_ci_hi=reps[0], n_reps=1)
    p05 = _percentile(reps_sorted, 0.05)
    p95 = _percentile(reps_sorted, 0.95)
    rng = random.Random(rng_seed)
    boot_medians = sorted(
        statistics.median(rng.choices(reps, k=len(reps))) for _ in range(n_boot)
    )
    alpha = (1.0 - ci) / 2.0
    return dict(
        median=median, p05=p05, p95=p95,
        boot_ci_lo=_percentile(boot_medians, alpha),
        boot_ci_hi=_percentile(boot_medians, 1.0 - alpha),
        n_reps=len(reps),
    )
