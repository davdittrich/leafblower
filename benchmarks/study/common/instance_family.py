"""instance_family.py -- parametric synthetic instance-family generator.

WU-3 (leafblower-2ouc.4). Design of record: docs/benchmark/DESIGN.md §3C
(instance-family axes), §7 (Dolan-More profiles need a problem
DISTRIBUTION), §11 DoD. Emits a frozen 30-100-instance sweep over
n / K / cell-cardinality / target-skew / margin-collinearity /
infeasibility-slack, each materialized lazily through a `gen:instance_family`
data_ref resolved by the WU-2 loader (common/problem_io.py).

STRICT SEPARATION (user constraint 2026-07-08): this file calls NOTHING in
leafblower. It only emits standardized problem specs (DataFrame + targets +
bounds), consumed downstream by adapters that themselves call leafblower's
harvest() -- out of this file's scope entirely.

Non-home-turf rationale, determinism, and R<->Python exact-parity strategy
(hand-rolled Park-Miller LCG, modulus 2^31-1, multiplier 48271 -- IEEE-754
double arithmetic on operands that stay exactly representable in a 53-bit
mantissa, so R and Python compute BIT-IDENTICAL uniforms): see
instance_family.R's module docstring (mirrored line-for-line; this file is
the Python port of that exact algorithm -- do not let the two drift).

Loader wiring (no problem_io.py edit -- out of WU-3 file scope):
`install_gen_resolver()` monkey-patches problem_io._resolve_data_ref at
runtime. Python resolves the free name `_resolve_data_ref` used inside
`load_problem_spec` from the module's global namespace at CALL time (not
definition time), so reassigning `problem_io._resolve_data_ref` after import
changes what `load_problem_spec()` invokes. The patch intercepts ONLY
`gen:instance_family...` data_refs it recognizes; every other `gen:<id>`
(including WU-2's own frozen `gen:toy_recipe` regression guard in
test_problem_io.py) falls through unchanged to the original resolver, which
still raises its WU-3-not-yet-implemented error.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import pandas as pd

import problem_io

# ---- LCG core --------------------------------------------------------

_LCG_M = 2147483647  # 2^31 - 1 (Mersenne prime, Park-Miller modulus)
_LCG_A = 48271  # Park-Miller 1993 revised multiplier


def _lcg_next(state: int) -> int:
    return (_LCG_A * state) % _LCG_M


def _lcg_seed0(master_seed: int, local_seed: int) -> int:
    # EXACTNESS NOTE: master_seed (~2.03e7) * 2654435761 (~2.65e9) ~= 5.38e16,
    # which EXCEEDS a double's 53-bit-mantissa exact range (2^53 ~= 9.007e15)
    # -- R's `*` on doubles rounds that product while Python's arbitrary-
    # precision int does not, silently breaking R<->Python bit-parity. Fix
    # (mirrors instance_family.R): SMALL multiplier (40503) on the large,
    # run-constant master_seed (product ~8.2e11, exact); LARGE multiplier
    # (2654435761) on the small, per-instance local_seed (product <= 2.65e11
    # for local_seed <= 100 as used by build_instance_grid(), exact).
    raw = (master_seed * 40503 + local_seed * 2654435761 + 1) % (_LCG_M - 1)
    return raw + 1  # avoid the absorbing state 0


# ---- Zipf target distribution -----------------------------------------


def _zipf_probs(c: int, s: float) -> tuple[list[float], list[float]]:
    raw = [1.0] * c if s == 0 else [float(j) ** (-s) for j in range(1, c + 1)]
    total = sum(raw)
    probs = [x / total for x in raw]
    cum = []
    running = 0.0
    for p in probs:
        running += p
        cum.append(running)
    cum[-1] = 1.0
    return probs, cum


def _inv_cdf(u: float, cum: list[float]) -> int:
    for j, edge in enumerate(cum, start=1):
        if u <= edge:
            return j
    return len(cum)


# ---- Row-level category generation -------------------------------------


def _generate_categories(
    n: int, k: int, c: int, s: float, rho: float, seed: int
) -> tuple[list[list[int]], tuple[list[float], list[float]]]:
    zp = _zipf_probs(c, s)
    cum = zp[1]
    state = _lcg_seed0(20260708, seed)
    mat = [[0] * k for _ in range(n)]
    for i in range(n):
        state = _lcg_next(state)
        cat1 = _inv_cdf(state / _LCG_M, cum)
        mat[i][0] = cat1
        for kk in range(2, k + 1):
            state = _lcg_next(state)
            u_dep = state / _LCG_M
            if u_dep <= rho:
                mat[i][kk - 1] = ((cat1 - 1 + (kk - 1)) % c) + 1
            else:
                state = _lcg_next(state)
                mat[i][kk - 1] = _inv_cdf(state / _LCG_M, cum)
    return mat, zp


# ---- Axis level tables (frozen; must mirror instance_family.R exactly) --

IF_CARD_LEVELS = {"low": 4, "medium": 8, "high": 16}
IF_SKEW_LEVELS = {"none": 0, "moderate": 1, "extreme": 3}
IF_COND_LEVELS = {"well": 0.0, "moderate": 0.5, "ill": 0.9}
IF_INFEAS_MAX: dict[str, float | None] = {"loose": None, "moderate": 3.0, "tight": 1.05}
IF_N_LEVELS = [1000, 10000, 100000, 1580000]
IF_K_LEVELS = [2, 4, 9]
IF_MASTER_SEED = 20260708

# ---- Instance grid (frozen recipe; deterministic, no RNG) --------------


def build_instance_grid() -> list[dict[str, Any]]:
    g: list[dict[str, Any]] = []

    def add(n: int, k: int, card: str, skew: str, cond: str, infeas: str, seed: int, tag: str) -> None:
        g.append(dict(n=n, K=k, card=card, skew=skew, cond=cond, infeas=infeas, seed=seed, tag=tag))

    # A: primary n x K grid at baseline card/skew/cond/infeas. 1.58M capped
    # to K=9 only ("only a handful at 1.58M", DESIGN.md §4).
    for n in IF_N_LEVELS:
        for k in IF_K_LEVELS:
            if n == 1580000 and k != 9:
                continue
            add(n, k, "medium", "moderate", "well", "loose", 0, "primary")
    # B: OFAT stress on each difficulty axis at baseline n=10000, K=4.
    for lvl in ("low", "high"):
        add(10000, 4, lvl, "moderate", "well", "loose", 0, "ofat_card")
    for lvl in ("none", "extreme"):
        add(10000, 4, "medium", lvl, "well", "loose", 0, "ofat_skew")
    for lvl in ("moderate", "ill"):
        add(10000, 4, "medium", "moderate", lvl, "loose", 0, "ofat_cond")
    for lvl in ("moderate", "tight"):
        add(10000, 4, "medium", "moderate", "well", lvl, 0, "ofat_infeas")
    # C: combined worst-case stress (all axes extreme) at a few n levels.
    for n in (10000, 100000, 1580000):
        add(n, 9, "high", "extreme", "ill", "tight", 0, "stress")
    add(1000, 2, "high", "extreme", "ill", "tight", 0, "stress")
    # D: pairwise axis interactions at baseline n=10000, K=4.
    add(10000, 4, "high", "moderate", "ill", "loose", 0, "pairwise")
    add(10000, 4, "medium", "extreme", "ill", "loose", 0, "pairwise")
    add(10000, 4, "high", "extreme", "well", "loose", 0, "pairwise")
    add(10000, 4, "high", "moderate", "well", "tight", 0, "pairwise")
    # E: seed replicates for profile density (same structural params, fresh draws).
    for seed in (1, 2):
        add(10000, 4, "medium", "moderate", "well", "loose", seed, "replicate")
    for seed in (1, 2):
        add(1000, 4, "medium", "moderate", "well", "loose", seed, "replicate")
    for seed in (1, 2):
        add(10000, 4, "high", "moderate", "well", "loose", seed, "replicate")
    return g


def _instance_id(inst: dict[str, Any]) -> str:
    return (
        f"if_n{inst['n']}_K{inst['K']}_{inst['card']}_{inst['skew']}_"
        f"{inst['cond']}_{inst['infeas']}_s{inst['seed']}"
    )


def _data_ref(inst: dict[str, Any]) -> str:
    return (
        f"gen:instance_family?n={inst['n']}&K={inst['K']}&card={inst['card']}"
        f"&skew={inst['skew']}&cond={inst['cond']}&infeas={inst['infeas']}&seed={inst['seed']}"
    )


# ---- Spec materialization (schema.json-conformant) ----------------------


def _build_spec(inst: dict[str, Any]) -> dict[str, Any]:
    k = inst["K"]
    c = IF_CARD_LEVELS[inst["card"]]
    s = IF_SKEW_LEVELS[inst["skew"]]
    probs, _ = _zipf_probs(c, s)
    cats = [f"c{j}" for j in range(1, c + 1)]
    target_one = dict(zip(cats, probs))
    margins = [f"m{j}" for j in range(1, k + 1)]
    targets = {m: dict(target_one) for m in margins}
    bmax = IF_INFEAS_MAX[inst["infeas"]]
    families = ["kl", "chi2", "logit", "newton_kl", "minimax"]
    if k == 2:
        families = families + ["ot"]
    return dict(
        id=_instance_id(inst),
        data_ref=_data_ref(inst),
        margins=margins,
        targets=targets,
        bounds=dict(min=0, max=bmax),
        tol=1e-8,
        objective_families=families,
        K=k,
    )


def generate_instance_family_specs(out_dir: str | Path = "benchmarks/study/spec") -> dict[str, Any]:
    """Writes one spec/instance_family/<id>.json per instance plus the grid
    manifest spec/instance_family.json (axis levels + instance roster).
    Cheap regardless of n: data is NOT materialized here (gen: is lazy)."""
    out_dir = Path(out_dir)
    inst_dir = out_dir / "instance_family"
    inst_dir.mkdir(parents=True, exist_ok=True)
    grid = build_instance_grid()
    ids = []
    for inst in grid:
        spec = _build_spec(inst)
        ids.append(spec["id"])
        with open(inst_dir / f"{spec['id']}.json", "w") as f:
            json.dump(spec, f, indent=2)
    manifest = dict(
        recipe="instance_family",
        master_seed=IF_MASTER_SEED,
        n_instances=len(grid),
        axis_levels=dict(
            n=IF_N_LEVELS,
            K=IF_K_LEVELS,
            card=IF_CARD_LEVELS,
            skew=IF_SKEW_LEVELS,
            cond=IF_COND_LEVELS,
            infeas={"loose": "Inf", "moderate": 3.0, "tight": 1.05},
        ),
        instances=ids,
        frozen=True,
    )
    with open(out_dir / "instance_family.json", "w") as f:
        json.dump(manifest, f, indent=2)
    return dict(ids=ids, out_dir=str(inst_dir))


# ---- gen: origin resolution + loader wiring -----------------------------


def _parse_query(qs: str) -> dict[str, str]:
    out = {}
    for part in qs.split("&"):
        key, _, val = part.partition("=")
        out[key] = val
    return out


def _resolve_gen_data(data_ref: str) -> pd.DataFrame | None:
    """Resolves a `gen:instance_family?...` data_ref into a DataFrame with
    columns m1..mK, matching problem_io's expected data-frame shape. Returns
    None if the recipe id isn't "instance_family" (caller falls through to
    the original resolver for unrecognized recipes)."""
    recipe_and_query = data_ref[len("gen:"):]
    recipe_id, _, query = recipe_and_query.partition("?")
    if recipe_id != "instance_family":
        return None
    q = _parse_query(query)
    n = int(q["n"])
    k = int(q["K"])
    c = IF_CARD_LEVELS[q["card"]]
    s = IF_SKEW_LEVELS[q["skew"]]
    rho = IF_COND_LEVELS[q["cond"]]
    seed = int(q["seed"])
    mat, _ = _generate_categories(n, k, c, s, rho, seed)
    cols = {f"m{j}": [f"c{row[j - 1]}" for row in mat] for j in range(1, k + 1)}
    return pd.DataFrame(cols)


_installed = False


def install_gen_resolver() -> None:
    """Monkey-patches problem_io._resolve_data_ref to dispatch
    `gen:instance_family?...` to this module's generator, delegating every
    other data_ref (incl. unrecognized `gen:<id>`) to the original function.
    Idempotent; safe to call more than once."""
    global _installed
    if _installed:
        return
    original = problem_io._resolve_data_ref

    def patched(data_ref: str, inline_data: list[dict[str, Any]] | None) -> pd.DataFrame:
        if data_ref.startswith("gen:instance_family"):
            out = _resolve_gen_data(data_ref)
            if out is not None:
                return out
        return original(data_ref, inline_data)

    problem_io._resolve_data_ref = patched
    _installed = True
