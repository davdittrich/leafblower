# STUDY-BRANCH-ONLY-DO-NOT-MERGE
# WU-12c (leafblower-2ouc.42) trajectory-capture runner: one (solver, problem)
# harvest with LBW_TRAJECTORY_AT / LBW_TRAJECTORY_OUT set by the caller, so the
# solver writes its per-iteration RQ3 convergence CSV. Cold path (byte-identical
# weights, proven). Mirrors python/leafblower_adapter.py's harvest arg mapping.
import os
import sys

_COMMON = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common")
sys.path.insert(0, _COMMON)

import instance_family  # noqa: E402
instance_family.install_gen_resolver()  # gen: instance specs (harmless for file:)
import problem_io  # noqa: E402
import leafblower as lb  # noqa: E402


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit("usage: capture_one.py <method> <spec_path>")
    method, spec_path = sys.argv[1], sys.argv[2]
    p = problem_io.load_problem_spec(spec_path)
    # LBW_TRAJECTORY_AT / LBW_TRAJECTORY_OUT set in the env by the caller.
    try:
        lb.harvest(
            data=p["data"], targets=p["targets"], method=method,
            min_weight=p["bounds"]["min"], max_weight=p["bounds"]["max"],
            design_weights=p["design_weights"],
            convergence={"absolute": p["tol"]}, attach_weights=False,
        )
    except Exception as e:  # infeasible/error cells have no convergence curve -- skip cleanly
        sys.stderr.write(f"skip {method}/{os.path.basename(spec_path)}: {e}\n")


if __name__ == "__main__":
    main()
