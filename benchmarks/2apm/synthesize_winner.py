#!/usr/bin/env python3
"""
benchmarks/2apm/synthesize_winner.py
T3 (leafblower-2apm.3) post-processor: median-reduce NREPS=3 per (dgp, lang),
apply Tier-1 lexicographic ranking, emit winner.toon.

Adapted from benchmarks/yh0l/synthesize_winner.py with chebyshev-specific
code-attribution paragraph (see CODE_ATTRIBUTION_YAML / FIX_TARGET_YAML).

Output: benchmarks/2apm/winner.toon  (overwritten)
"""

import csv
import math
import statistics
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────

APMDIR = Path(__file__).resolve().parent
R_CSV  = APMDIR / "results_R_all.csv"
PY_CSV = APMDIR / "results_Py_all.csv"
OUT    = APMDIR / "winner.toon"

# ── Ranking constants ────────────────────────────────────────────────────────
# Tier-1 lexicographic order: margin_kl → wall_ms → weight_kl → max_err
KL_TIE_ABS   = 1e-12
WALL_TIE_REL = 0.10
WKL_TIE_ABS  = 1e-12
ME_TIE_ABS   = 1e-12

CONF_SCALE = 100.0


def read_csv(path):
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
    return int(max(50, min(99, 50 + 50 * math.tanh(rel_diff * CONF_SCALE))))


def rank_dgp(r_row, py_row):
    r_mkl  = r_row["margin_kl"]
    py_mkl = py_row["margin_kl"]
    diff_mkl = r_mkl - py_mkl

    if abs(diff_mkl) >= KL_TIE_ABS:
        winner = "Py" if py_mkl < r_mkl else "R"
        mn = min(abs(r_mkl), abs(py_mkl))
        rel = abs(diff_mkl) / max(mn, 1e-300)
        conf = confidence(rel)
        note = f"margin_kl: R={r_mkl:.6g} Py={py_mkl:.6g} rel_diff={rel:.3g}"
        return winner, conf, note

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

    return "tie", 50, "all Tier-1 metrics within tie thresholds"


def aggregate_winner(per_dgp_results):
    wins  = {"R": 0, "Py": 0, "tie": 0}
    for entry in per_dgp_results:
        wins[entry["dgp_winner"]] += 1
    if wins["R"] == 0 and wins["Py"] > 0:
        return "Py"
    if wins["Py"] == 0 and wins["R"] > 0:
        return "R"
    if wins["R"] == 0 and wins["Py"] == 0:
        return "tie"
    return "tie"


def build_winner_basis(per_dgp_results, agg_winner):
    lines = []
    for e in per_dgp_results:
        dgp = e["dgp"]
        w   = e["dgp_winner"]
        c   = e["confidence"]
        note = e.get("rank_note", "")
        lines.append(f"{dgp}: {w} (confidence {c}) — {note}.")
    summary = f"Aggregate winner: {agg_winner}. " + " ".join(lines)
    return summary


# ── Code attribution (chebyshev) ─────────────────────────────────────────────
# Forensic finding: ieppa warm-start `inner_max_iter` lower bound differs
# between r_bridge.cpp and c_api.cpp chebyshev dispatch blocks.

CODE_ATTRIBUTION_YAML = """\
  code_attribution:
    r_side:
      file: src/r_bridge.cpp
      line: 648
      quote: "st_warm.inner_max_iter = std::max(5, std::min(100, st.inner_max_iter / 10));"
    py_side:
      file: src/c_api.cpp
      line: 359
      quote: "st_warm.inner_max_iter = std::max(50, std::min(100, st.inner_max_iter / 10));"
    differentiator: >-
      Both R (r_bridge.cpp:643-657) and Python (c_api.cpp:351-367) run the same ieppa warm-start
      sequence before invoking lbw::chebyshev_ipm(st, w_warm, delta_warm). The blocks are
      structurally identical EXCEPT for one constant in the inner_max_iter clamp expression:
      the R bypass uses std::max(5, ...) while c_api.cpp uses std::max(50, ...). For default
      max_iterations=3000 the inner cap is 100/10 = 10 inner Newton iters per ieppa outer iter
      on both sides (clamped above by 100, below by 5/50 respectively); the divergence kicks in
      when st.inner_max_iter / 10 is below 50. Additionally, c_api.cpp:352 carries the comment
      "ieppa warm-start; mirrors r_bridge.cpp:628-657" — confirming this block was hand-copied
      from the R bypass but drifted by one literal (50 vs 5 in the std::max floor). The warm-start `w_warm`
      and `delta_warm` propagate into chebyshev_ipm's IPM Newton initial state via
      lambda_A/lambda_B initialization derived from w_warm. Different warm states → different
      Newton trajectories → different fixed points within the 500-iter Chebyshev cap. T2
      already ruled out BLAS thread nondet and confirmed shared-code numerical drift (cause c).
      This file:line pair is the sole textual divergence in the chebyshev dispatch path between
      r_bridge.cpp and c_api.cpp."""

FIX_TARGET_YAML = """\
  fix_target: >-
    Files T4 must modify: src/c_api.cpp (single-line change at line 359). Replace
    `std::max(50, std::min(100, st.inner_max_iter / 10))` with
    `std::max(5, std::min(100, st.inner_max_iter / 10))` to align the Python path with the
    R bridge convention (matching the comment "mirrors r_bridge.cpp:628-657" at c_api.cpp:352).
    Direction (port Py→R) is chosen because (a) R won on fulldata margin_kl (0.7410 vs 0.7456),
    (b) the comment explicitly cites r_bridge as the template, and (c) the R-side bound (5) is
    the safer floor for users who pass small max_iterations values where Py would silently
    inflate the warm-start budget. Verification: post-fix, run benchmarks/2apm/winner_bench.{R,py}
    + synthesize_winner.py and confirm (margin_kl, weight_kl, max_err, iter) tuples match across
    R and Py to rtol≤1e-12 on all three DGPs. No changes needed in src/r_bridge.cpp,
    src/chebyshev.cpp, src/calib_dispatch.hpp, R/harvest.R, or python/leafblower/_harvest.py."""


def fmt_float(v):
    if v == 0.0:
        return "0.0"
    if abs(v) < 1e-10 or abs(v) >= 1e6:
        return repr(v)
    return f"{v:.8g}"


def emit_toon(medians, per_dgp_results, agg_winner, has_winner):
    winner_basis = build_winner_basis(per_dgp_results, agg_winner)

    branch = "cause_c_with_winner" if has_winner else "cause_c_no_winner"
    winner_field = agg_winner if has_winner else "none"

    lines = []
    lines.append("task_id: T3_2apm")
    lines.append("success: true")
    lines.append("data:")
    lines.append(f"  branch: {branch}")
    lines.append(f"  reduction: median_over_3_reps")
    lines.append(f"  winner: {winner_field}")
    lines.append(f"  winner_basis: >-")
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
        lines.append(f"      r_weight_kl: {fmt_float(r_m['weight_kl'])}")
        lines.append(f"      py_weight_kl: {fmt_float(py_m['weight_kl'])}")
        lines.append(f"      r_max_abs_diff_stability: {fmt_float(r_m['max_err'])}")
        lines.append(f"      py_max_abs_diff_stability: {fmt_float(py_m['max_err'])}")
        lines.append(f"      dgp_winner: {entry['dgp_winner']}")
        lines.append(f"      confidence: {entry['confidence']}")
        if entry.get("rank_note"):
            lines.append(f"      note: \"{entry['rank_note']}\"")
        lines.append("")

    lines.append(CODE_ATTRIBUTION_YAML)
    lines.append("")
    if has_winner:
        lines.append(FIX_TARGET_YAML)
    else:
        lines.append("  fix_target: null")
    lines.append("")
    lines.append(f"  no_winner_followup_ticket_required: {'false' if has_winner else 'true'}")
    lines.append("  spec_amendment_needed: false")
    lines.append("  spec_amendment_path: null")
    lines.append("")
    lines.append("error_log: null")

    return "\n".join(lines) + "\n"


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
    has_winner = agg_winner in ("R", "Py")
    print(f"\nAggregate winner: {agg_winner} (branch={'cause_c_with_winner' if has_winner else 'cause_c_no_winner'})")

    out_name = "winner.toon" if has_winner else "no_winner.toon"
    out_path = APMDIR / out_name
    toon = emit_toon(medians, per_dgp_results, agg_winner, has_winner)

    try:
        import yaml
        parsed = yaml.safe_load(toon)
        assert parsed["data"]["reduction"] == "median_over_3_reps"
        assert parsed["data"]["winner"] == (agg_winner if has_winner else "none")
        assert len(parsed["data"]["per_dgp"]) == len(per_dgp_results)
        print("YAML parse: OK")
    except ImportError:
        print("WARNING: pyyaml not installed; YAML parse check skipped")

    out_path.write_text(toon)
    print(f"Wrote {out_path}")

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
