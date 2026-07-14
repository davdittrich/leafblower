import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import pandas as pd  # noqa: E402


def test_trajectory_csv_smoke(tmp_path):
    """leafblower-9nuo: regression smoke for the former SIGSEGV crash path.

    Before leafblower-9nuo (WU-1), `write_trajectory_csv` in
    src/oris_trajectory.cpp used std::ofstream/iostream, which SIGSEGVs in
    std::codecvt do_unshift when driven from the module-side _leafblower.so
    (it static-links libstdc++). This test drives that exact code path (an
    ORIS solve with LBW_TRAJECTORY_AT/LBW_TRAJECTORY_OUT set) and asserts the
    process does NOT segfault, and that the writer actually produced a
    non-vacuous CSV (header + >=1 data row).
    """
    from leafblower import harvest

    out_path = tmp_path / "trajectory.csv"
    os.environ["LBW_TRAJECTORY_AT"] = "1,2,3"
    os.environ["LBW_TRAJECTORY_OUT"] = str(out_path)
    try:
        df = pd.DataFrame({
            "a": [["1", "2", "3"][i % 3] for i in range(600)],
            "b": [["x", "y"][i % 2] for i in range(600)],
        })
        targets = {
            "a": {"1": 0.5, "2": 0.3, "3": 0.2},
            "b": {"x": 0.6, "y": 0.4},
        }
        # Reaching this call and returning proves no SIGSEGV in the writer.
        harvest(df, targets, method="oris", max_iterations=50,
                attach_weights=False)
    finally:
        os.environ.pop("LBW_TRAJECTORY_AT", None)
        os.environ.pop("LBW_TRAJECTORY_OUT", None)

    assert out_path.exists(), "trajectory CSV was not written"
    lines = out_path.read_text().splitlines()
    assert lines[0] == "iter,errRp"
    assert len(lines) >= 2, "trajectory CSV has no data rows (vacuous smoke)"
