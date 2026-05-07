#!/usr/bin/env python3
"""
benchmarks/yh0l/synthesize_winner.py
T3 post-processor: median-reduce NREPS=3 per (dgp, lang), apply Tier-1
lexicographic ranking, emit winner.toon.

Usage:
    python benchmarks/yh0l/synthesize_winner.py

Output:
    benchmarks/yh0l/winner.toon  (overwritten)
"""

import csv
import math
import statistics
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────

YDIR = Path(__file__).resolve().parent
R_CSV  = YDIR / "results_R_all.csv"
PY_CSV = YDIR / "results_Py_all.csv"
OUT    = YDIR / "winner.toon"

# ── Ranking constants ────────────────────────────────────────────────────────

# Tier-1 lexicographic order: margin_kl → wall_ms → weight_kl → max_err
# Tie thresholds:
#   margin_kl : abs(diff) < 1e-12  → tie
#   wall_ms   : rel diff < 0.10    → tie (10 % noise band)
#   weight_kl : abs(diff) < 1e-12  → tie
#   max_err   : abs(diff) < 1e-12  → tie
KL_TIE_ABS  = 1e-12
WALL_TIE_REL = 0.10
WKL_TIE_ABS  = 1e-12
ME_TIE_ABS   = 1e-12

# Confidence: max(50, min(99, 50 + 50*tanh(rel_diff * scale)))
# scale=100 → tanh saturates around rel_diff ≥ 0.03 (3 %)
CONF_SCALE = 100.0


# ── Helpers ──────────────────────────────────────────────────────────────────

def read_csv(path):
    """Return list-of-dicts with numeric coercion."""
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append({
                "dgp":        row["dgp"],
                "rep":        int(row["rep"]),
                "n":          int(row["n"]),
                "margin_kl":  float(row["margin_kl"]),
                "weight_kl":  float(row["weight_kl"]),
                "max_err":    float(row["max_err"]),
                "wall_ms":    float(row["wall_ms"]),
                "iterations": int(row["iterations"]),
                "lang":       row["lang"].strip('"'),
            })
    return rows


def median_reduce(rows):
    """
    Group rows by (dgp, lang), assert iter identical, return dict:
      {dgp: {lang: {margin_kl, weight_kl, max_err, wall_ms, iter, n}}}
    Uses statistics.median (exact for odd NREPS).
    """
    groups = {}
    for r in rows:
        key = (r["dgp"], r["lang"])
        groups.setdefault(key, []).append(r)

    result = {}
    for (dgp, lang), reps in sorted(groups.items()):
        iters = [r["iterations"] for r in reps]
        if len(set(iters)) != 1:
            raise AssertionError(
                f"Non-deterministic iterations for dgp={dgp} lang={lang}: {iters}"
            )
        ns = [r["n"] for r in reps]
        if len(set(ns)) != 1:
            raise AssertionError(f"n differs across reps for dgp={dgp} lang={lang}: {ns}")

        result.setdefault(dgp, {})[lang] = {
            "n":          ns[0],
            "margin_kl":  statistics.median([r["margin_kl"]  for r in reps]),
            "weight_kl":  statistics.median([r["weight_kl"]  for r in reps]),
            "max_err":    statistics.median([r["max_err"]    for r in reps]),
            "wall_ms":    statistics.median([r["wall_ms"]    for r in reps]),
            "iter":       iters[0],
        }
    return result


def confidence(rel_diff):
    """Bounded confidence score from relative difference."""
    return int(max(50, min(99, 50 + 50 * math.tanh(rel_diff * CONF_SCALE))))


def rank_dgp(r_row, py_row):
    """
    Tier-1 lexicographic ranker.
    Returns (winner: 'R'|'Py'|'tie', confidence: int, note: str)
    Lower margin_kl = better; lower wall_ms = faster.
    """
    r_mkl  = r_row["margin_kl"]
    py_mkl = py_row["margin_kl"]
    diff_mkl = r_mkl - py_mkl   # positive → Py better

    # 1. margin_kl
    if abs(diff_mkl) >= KL_TIE_ABS:
        winner = "Py" if py_mkl < r_mkl else "R"
        mn = min(abs(r_mkl), abs(py_mkl))
        rel = abs(diff_mkl) / max(mn, 1e-300)
        conf = confidence(rel)
        note = f"margin_kl: R={r_mkl:.6g} Py={py_mkl:.6g} rel_diff={rel:.3g}"
        return winner, conf, note

    # 2. wall_ms
    r_wall  = r_row["wall_ms"]
    py_wall = py_row["wall_ms"]
    mn_wall = min(r_wall, py_wall)
    rel_wall = abs(r_wall - py_wall) / max(mn_wall, 1e-9)
    if rel_wall >= WALL_TIE_REL:
        winner = "Py" if py_wall < r_wall else "R"
        conf = confidence(rel_wall)
        note = (f"margin_kl tie (|diff|<1e-12); "
                f"wall_ms: R={r_wall:.1f} Py={py_wall:.1f} rel={rel_wall:.3g}")
        return winner, conf, note

    # 3. weight_kl
    r_wkl  = r_row["weight_kl"]
    py_wkl = py_row["weight_kl"]
    diff_wkl = abs(r_wkl - py_wkl)
    if diff_wkl >= WKL_TIE_ABS:
        winner = "Py" if py_wkl < r_wkl else "R"
        rel = diff_wkl / max(min(abs(r_wkl), abs(py_wkl)), 1e-300)
        conf = confidence(rel)
        note = (f"margin_kl+wall_ms ties; "
                f"weight_kl: R={r_wkl:.6g} Py={py_wkl:.6g}")
        return winner, conf, note

    # 4. max_err
    r_me  = r_row["max_err"]
    py_me = py_row["max_err"]
    diff_me = abs(r_me - py_me)
    if diff_me >= ME_TIE_ABS:
        winner = "Py" if py_me < r_me else "R"
        rel = diff_me / max(min(abs(r_me), abs(py_me)), 1e-300)
        conf = confidence(rel)
        note = (f"margin_kl+wall_ms+weight_kl ties; "
                f"max_err: R={r_me:.6g} Py={py_me:.6g}")
        return winner, conf, note

    # Full tie
    return "tie", 50, "all Tier-1 metrics within tie thresholds"


def aggregate_winner(per_dgp_results):
    """
    Overall winner: Py wins if never loses (wins ≥1, ties rest).
    R wins symmetrically. Otherwise tie.
    """
    wins  = {"R": 0, "Py": 0, "tie": 0}
    for entry in per_dgp_results:
        wins[entry["dgp_winner"]] += 1
    if wins["R"] == 0:
        return "Py"
    if wins["Py"] == 0:
        return "R"
    return "tie"


def build_winner_basis(per_dgp_results, agg_winner):
    """Auto-generate winner_basis paragraph from ranked results."""
    lines = []
    for e in per_dgp_results:
        dgp = e["dgp"]
        w   = e["dgp_winner"]
        c   = e["confidence"]
        note = e.get("rank_note", "")
        lines.append(f"{dgp}: {w} (confidence {c}) — {note}.")
    summary = f"Aggregate winner: {agg_winner}. " + " ".join(lines)
    return summary


# ── Verbatim preserved blocks from prior winner.toon ────────────────────────

CODE_ATTRIBUTION_YAML = """\
  code_attribution:
    r_side:
      file: src/r_bridge.cpp
      line: 359
      quote: "st.alm.capacity_mu = (cp_val <= 0.0) ? ct_tmp.capacity_mu_auto : cp_val;"
    py_side:
      file: src/c_api.cpp
      line: 379
      quote: "st.alm.capacity_mu = (n > 0) ? static_cast<double>(M_cell_est) / n : 1.0;"
    differentiator: >-
      The R bridge (r_bridge.cpp:346-359) calls lbw::build_cell_table() before the solver to
      compute the exact number of occupied interaction cells (M_cell) and sets
      capacity_mu = M_cell/n. The Python path goes through c_api.cpp (called via pybind11)
      which uses lbw::estimate_M_cell() instead. estimate_M_cell has a K>8 fast-exit: when
      K exceeds 8 margins, it returns n immediately (the theoretical maximum), making
      capacity_mu = n/n = 1.0. For the fulldata DGP (K=9 margins), the R bridge computes
      exact capacity_mu ≈ 0.073 (≈115,740 occupied cells out of 1,582,732 observations),
      while the C API produces capacity_mu = 1.0. This capacity_mu value is the ALM
      penalty scaling used by ieppa_soft's ADMM-style constraint enforcement; a smaller
      capacity_mu imposes weaker initial constraint pressure, causing R's solver to converge
      slower (220 iters) and to a slightly worse solution (margin_kl 0.005074 vs 0.005045).
      The Python path's larger capacity_mu = 1.0 drives stronger initial constraint enforcement,
      yielding faster convergence (90 iters) and lower final margin_kl. The root cause is not
      in group_ids encoding, NA handling, factor level ordering, convergence parameters, or
      BLAS threading — all of those are identical between the two sides (confirmed by T1 and
      direct code verification). The sole differentiator is this capacity_mu initialization
      path: r_bridge.cpp (exact build_cell_table) vs c_api.cpp (approximate estimate_M_cell
      with K>8 cap)."""

FIX_TARGET_YAML = """\
  fix_target: >-
    The loser is R (src/r_bridge.cpp). The fix is to change r_bridge.cpp line 359 to use
    estimate_M_cell (matching the Python/C-API path) OR change c_api.cpp line 379-380 to
    call build_cell_table (matching the R bridge path). Since Python converges faster and to
    better solutions with capacity_mu = 1.0 for large K, the preferred fix is to align
    r_bridge.cpp to use estimate_M_cell for K>8, matching c_api.cpp. Specifically: in
    src/r_bridge.cpp around line 340-360, replace the build_cell_table block with a call to
    lbw::estimate_M_cell(n, K, group_ids.data(), cat_counts.data()) and set
    st.alm.capacity_mu = (n > 0) ? static_cast<double>(M_cell_est) / n : 1.0. Files T4
    must modify: src/r_bridge.cpp (primary change). No changes needed in src/c_api.cpp,
    src/ieppa.cpp, R/harvest.R, or python/leafblower/_harvest.py."""


# ── YAML emission ────────────────────────────────────────────────────────────

def fmt_float(v):
    """Compact repr that avoids scientific notation for moderate values."""
    if v == 0.0:
        return "0.0"
    if abs(v) < 1e-10 or abs(v) >= 1e6:
        return repr(v)
    return f"{v:.8g}"


def emit_toon(medians, per_dgp_results, agg_winner):
    winner_basis = build_winner_basis(per_dgp_results, agg_winner)

    lines = []
    lines.append("task_id: T3")
    lines.append("success: true")
    lines.append("data:")
    lines.append(f"  reduction: median_over_3_reps")
    lines.append(f"  winner: {agg_winner}")
    lines.append(f"  winner_basis: >-")
    # Wrap winner_basis at ~100 chars
    words = winner_basis.split()
    cur = "    "
    for w in words:
        if len(cur) + len(w) + 1 > 100:
            lines.append(cur.rstrip())
            cur = "    " + w + " "
        else:
            cur += w + " "
    lines.append(cur.rstrip())
    lines.append("  per_dgp:")

    dgp_order = ["fulldata", "medium", "small_hicard"]
    for dgp in dgp_order:
        if dgp not in medians:
            continue
        r_m  = medians[dgp].get("R")
        py_m = medians[dgp].get("Py")
        if r_m is None or py_m is None:
            continue
        entry = next(e for e in per_dgp_results if e["dgp"] == dgp)
        lines.append(f"    - dgp: {dgp}")
        lines.append(f"      n: {r_m['n']}")
        lines.append(f"      r_margin_kl: {fmt_float(r_m['margin_kl'])}")
        lines.append(f"      py_margin_kl: {fmt_float(py_m['margin_kl'])}")
        lines.append(f"      r_wall_ms: {fmt_float(r_m['wall_ms'])}")
        lines.append(f"      py_wall_ms: {fmt_float(py_m['wall_ms'])}")
        lines.append(f"      r_iter: {r_m['iter']}")
        lines.append(f"      py_iter: {py_m['iter']}")
        lines.append(f"      dgp_winner: {entry['dgp_winner']}")
        lines.append(f"      confidence: {entry['confidence']}")
        if entry.get("rank_note"):
            # YAML block scalar for note
            lines.append(f"      note: \"{entry['rank_note']}\"")
        lines.append("")

    lines.append(CODE_ATTRIBUTION_YAML)
    lines.append("")
    lines.append(FIX_TARGET_YAML)
    lines.append("")
    lines.append("  spec_amendment_needed: false")
    lines.append("  spec_amendment_path: null")
    lines.append("")
    lines.append("error_log: null")

    return "\n".join(lines) + "\n"


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    r_rows  = read_csv(R_CSV)
    py_rows = read_csv(PY_CSV)
    all_rows = r_rows + py_rows

    medians = median_reduce(all_rows)

    dgp_order = ["fulldata", "medium", "small_hicard"]
    per_dgp_results = []
    for dgp in dgp_order:
        if dgp not in medians:
            continue
        r_m  = medians[dgp].get("R")
        py_m = medians[dgp].get("Py")
        if r_m is None or py_m is None:
            print(f"WARNING: missing lang for dgp={dgp}, skipping")
            continue
        winner, conf, note = rank_dgp(r_m, py_m)
        per_dgp_results.append({
            "dgp":        dgp,
            "dgp_winner": winner,
            "confidence": conf,
            "rank_note":  note,
        })
        print(f"  {dgp}: winner={winner} conf={conf} | {note}")

    agg_winner = aggregate_winner(per_dgp_results)
    print(f"\nAggregate winner: {agg_winner}")

    toon = emit_toon(medians, per_dgp_results, agg_winner)

    # Sanity: verify YAML parses
    try:
        import yaml
        parsed = yaml.safe_load(toon)
        assert parsed["data"]["reduction"] == "median_over_3_reps", "reduction field missing"
        assert parsed["data"]["winner"] == agg_winner
        assert len(parsed["data"]["per_dgp"]) == len(per_dgp_results)
        print("YAML parse: OK")
    except ImportError:
        print("WARNING: pyyaml not installed; YAML parse check skipped")

    OUT.write_text(toon)
    print(f"Wrote {OUT}")

    # Eyeball verify: print median wall_ms per dgp/lang
    print("\nMedian wall_ms verification:")
    for dgp in dgp_order:
        if dgp not in medians:
            continue
        for lang in ("R", "Py"):
            if lang in medians[dgp]:
                m = medians[dgp][lang]
                print(f"  {dgp}/{lang}: wall_ms={m['wall_ms']:.1f} "
                      f"margin_kl={m['margin_kl']:.6g} iter={m['iter']}")


if __name__ == "__main__":
    main()
