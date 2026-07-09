"""benchmarks/study/python/adapter_ott.py -- WU-OTT (leafblower-2ouc.18).

Dedicated adapter for the ott-jax (JAX/XLA) Sinkhorn competitor, carved out
of competitors.py (WU-6, leafblower-2ouc.7) because JAX/XLA is a footgun
under this repo's single-thread-BLAS + float64 protocol (CLAUDE.md):

1. x64 (MANDATORY, module top): `jnp.float64` silently downcasts to
   float32 unless `jax.config.update("jax_enable_x64", True)` runs before
   any array op in the process (mirrors the OMP-before-numpy rule). Done
   immediately after `import jax` below. Verified in
   test_adapter_ott.py::test_x64_actually_in_effect (a test that would FAIL
   under jax's float32 default).

2. XLA threads (MANDATORY, module top): XLA reads XLA_FLAGS at process
   init, before `import jax`, so `os.environ.setdefault("XLA_FLAGS", ...)`
   runs before the `import jax` line below.

   VERIFIED THIS WU (jaxlib 0.10.2, 2026-07-09, 32-core host) -- do not
   silently trust the commonly-copied FAQ snippet:
     - `XLA_FLAGS="--intra_op_parallelism_threads=1"` -> FATAL
       "Unknown flag in XLA_FLAGS" (not a registered flag in this build).
     - `XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1"`
       (the textbook snippet, second token missing "--") does NOT error,
       but the second token is silently dropped as a non-flag positional
       arg -- it does nothing.
     - Neither `--xla_cpu_multi_thread_eigen=false` alone, NOR that flag
       plus OMP_NUM_THREADS=1/OPENBLAS_NUM_THREADS=1/MKL_NUM_THREADS=1,
       measurably reduce XLA's CPU intra-op parallelism: measured
       (ru_utime+ru_stime)/wall-clock ratio for a jitted float64 matmul was
       ~27x with or without these settings (identical, within noise).
     - Only OS-level CPU affinity (`os.sched_setaffinity`, applied BEFORE
       `import jax` or as a temporary restriction -- XLA sizes its
       CPU thread pool from the process's affinity mask at backend init)
       verifiably forces single-thread execution: ratio -> 1.00 measured.
   This CONFIRMS -- rather than contradicts -- DESIGN.md Section 2's
   stated reason ott-jax is excluded from ranked {1,4}-thread timing by
   the WU-8 run matrix ("XLA ignores OMP/BLAS env"): no env-var-only
   mechanism reproduces the other arms' OMP/BLAS thread-count protocol in
   this jaxlib build. The XLA_FLAGS value + x64 state actually in effect
   is disclosed via `xla_config_snapshot()` in every vignette row, per
   DESIGN.md's disclosure requirement. `ott_scaling_vignette(...,
   pin_single_core=True)` additionally offers a genuinely single-threaded
   row using the verified `os.sched_setaffinity` mechanism -- opt-in only
   (default False), and the prior affinity mask is restored afterwards, so
   importing this module has no global process side effect.

Shared OT-reduction plumbing (`_adapt` decorator, `_HP`, `_OT_EPSILON`,
`_ot_seed_cost`, `_ot_weights_from_coupling`) is imported FROM
competitors.py (single source, pragmatic DRY) rather than duplicated.
competitors.py's own re-export of `run_ott_jax_sinkhorn` defers ITS import
of this module to call time (see competitors.py comment) specifically to
avoid a competitors<->adapter_ott circular import, since this module needs
`_adapt` et al. from competitors.py eagerly at decoration time below.
"""

from __future__ import annotations

import os

os.environ.setdefault("XLA_FLAGS", "--xla_cpu_multi_thread_eigen=false")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import jax  # noqa: E402

jax.config.update("jax_enable_x64", True)

import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from ott.geometry import geometry  # noqa: E402
from ott.problems.linear import linear_problem  # noqa: E402
from ott.solvers.linear import sinkhorn  # noqa: E402

_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

import competitors as _c  # noqa: E402  (single source for shared OT/adapter plumbing)

_HP = _c._HP
_OT_EPSILON = _c._OT_EPSILON
_ot_seed_cost = _c._ot_seed_cost
_ot_weights_from_coupling = _c._ot_weights_from_coupling
_adapt = _c._adapt


def xla_config_snapshot() -> dict:
    """DESIGN.md Section 2 disclosure requirement: the XLA/x64 process
    state actually in effect, recorded alongside every ranked-timing-
    EXCLUDED ott-jax result (contract.md's harmonized 7-key contract has no
    dedicated field for this -- it lives in vignette/diagnostic output)."""
    return dict(
        xla_flags=os.environ.get("XLA_FLAGS", ""),
        jax_enable_x64=bool(jax.config.jax_enable_x64),
        jax_local_device_count=int(jax.local_device_count()),
        jax_default_backend=jax.default_backend(),
        intra_op_threads_env_pinning_verified_effective=False,
        note=(
            "jaxlib 0.10.2 (verified this WU): 'intra_op_parallelism_threads' "
            "is not a registered XLA_FLAGS token (fatal 'Unknown flag' with "
            "'--', silently dropped without it); xla_cpu_multi_thread_eigen="
            "false + OMP/OPENBLAS/MKL_NUM_THREADS=1 do NOT reduce measured "
            "CPU intra-op parallelism (cpu-time/wall-time ratio ~27x on a "
            "32-core host, unchanged by these settings). This confirms "
            "DESIGN.md Section 2's stated rationale for excluding ott-jax "
            "from ranked {1,4}-thread timing. Only OS-level CPU affinity "
            "(os.sched_setaffinity) verifiably forces single-thread "
            "execution (ratio -> 1.00); see "
            "ott_scaling_vignette(pin_single_core=True)."
        ),
    )


# ---------------------------------------------------------------------------
# ott-jax sinkhorn (registry id "ott_jax_sinkhorn"; OT family, K_max=2).
# Carved from competitors.py (WU-OTT). Per DESIGN.md Section 2 / contract
# DoD: wrapped + home-turf golden'd like every other competitor, but
# EXCLUDED FROM RANKED TIMING by the WU-8 run matrix (XLA/JAX thread-
# pinning under {1,4}-thread OMP/BLAS caps is not apples-to-apples with the
# other arms, see module docstring) -- that exclusion is a run-matrix
# (WU-8) concern; this adapter itself is fully functional.
# Home-turf golden: ott-jax's own pointcloud/Sinkhorn quickstart shape,
# adapted here to the same support-cost 2-margin reduction as POT so both
# OT arms are pinned against the identical IPF golden.
# ---------------------------------------------------------------------------

@_adapt("ott_jax_sinkhorn")
def run_ott_jax_sinkhorn(problem: dict) -> tuple:
    margins = problem["margins"]
    if len(margins) != 2:
        return None, None, "bad_arg", (
            f"ott_jax_sinkhorn: OT reduction requires K==2 (registry K_max=2); got K={len(margins)}"
        )
    m1, m2 = margins
    df = problem["data"]
    levels1 = list(problem["targets"][m1].keys())
    levels2 = list(problem["targets"][m2].keys())
    a = jnp.asarray([problem["targets"][m1][k] for k in levels1], dtype=jnp.float64)
    b = jnp.asarray([problem["targets"][m2][k] for k in levels2], dtype=jnp.float64)
    M = _ot_seed_cost(df[m1].to_numpy(), df[m2].to_numpy(), levels1, levels2, problem["design_weights"])

    hp = _HP["ott_jax_sinkhorn"]
    geom = geometry.Geometry(cost_matrix=jnp.asarray(M), epsilon=_OT_EPSILON)
    lp = linear_problem.LinearProblem(geom, a=a, b=b)
    solver = sinkhorn.Sinkhorn(threshold=hp["threshold"], min_iterations=hp["min_iterations"],
                                max_iterations=hp["max_iterations"], inner_iterations=hp["inner_iterations"])
    out = solver(lp)

    u = np.asarray(out.scalings[0], dtype=np.float64)
    v = np.asarray(out.scalings[1], dtype=np.float64)
    weights = _ot_weights_from_coupling(u, v, df[m1].to_numpy(), df[m2].to_numpy(),
                                         levels1, levels2, problem["design_weights"])
    status = "converged" if bool(out.converged) else "no_conv"
    return weights, int(out.n_iters), status, None


# ---------------------------------------------------------------------------
# Scaling vignette (n-sweep). NOT part of ranked {1,4}-thread timing (WU-8
# run matrix exclusion) -- a standalone diagnostic with disclosed XLA/x64
# config, per DESIGN.md Section 2's disclosure requirement.
# ---------------------------------------------------------------------------

def _synthetic_2margin_problem(n: int, problem_id: str) -> dict:
    """Scalable synthetic 2-margin problem for the n-sweep: n rows drawn
    from the same interacting (A/B x X/Y) cell structure as the hand-solved
    IPF golden in test_competitors.py (odds ratio != 1, genuine 2-D
    adjustment required), just resized -- every n in the sweep exercises
    the identical OT-reduction shape."""
    rng = np.random.default_rng(0)
    m1 = rng.choice(["A", "B"], size=n, p=[0.5, 0.5])
    interact_p = np.where(m1 == "A", 0.8, 0.2)
    m2 = np.where(rng.random(n) < interact_p, "X", "Y")
    df = pd.DataFrame({"m1": m1, "m2": m2})
    return dict(
        id=problem_id,
        data=df,
        design_weights=np.ones(n, dtype=np.float64),
        margins=["m1", "m2"],
        targets={"m1": {"A": 0.7, "B": 0.3}, "m2": {"X": 0.6, "Y": 0.4}},
        bounds={"min": 0.0, "max": float("inf")},
        tol=1e-3,
        objective_families=["kl"],
        K=2,
    )


def ott_scaling_vignette(n_values: list[int], pin_single_core: bool = False) -> dict:
    """Standalone n-sweep timing vignette for ott-jax (excluded from ranked
    {1,4}-thread timing -- WU-8 run matrix, DESIGN.md Section 2). Every
    entry is run under, and the return value discloses, the XLA/x64
    process config in effect (`xla_config_snapshot()`).

    `pin_single_core=True` additionally applies the verified OS-level CPU
    affinity mechanism (module docstring) for the duration of the sweep
    only, restoring the prior affinity mask afterwards -- opt-in, so the
    default call leaves no global process side effect.
    """
    config = xla_config_snapshot()
    prior_affinity = None
    if pin_single_core and hasattr(os, "sched_getaffinity"):
        prior_affinity = os.sched_getaffinity(0)
        os.sched_setaffinity(0, {next(iter(prior_affinity))})
        config = dict(config, pinned_single_core=True,
                      pinned_affinity_mask=sorted(os.sched_getaffinity(0)))
    try:
        rows = []
        for n in n_values:
            problem = _synthetic_2margin_problem(n, f"ott-vignette-n{n}")
            res = run_ott_jax_sinkhorn(problem)
            rows.append(dict(n=n, wall_time_s=res["wall_time_s"], status=res["status"],
                              converged=res["converged"], iterations=res["iterations"]))
        return dict(config=config, rows=rows, excluded_from_ranked_timing=True)
    finally:
        if prior_affinity is not None:
            os.sched_setaffinity(0, prior_affinity)
