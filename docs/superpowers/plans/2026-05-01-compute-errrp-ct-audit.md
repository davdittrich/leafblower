# compute_errRp_ct Dead-Code Audit Implementation Plan

**Goal:** Verify whether `compute_errRp_ct` (raking.cpp:34-49) has any callers in the codebase, and delete it if not.
**Architecture:** `compute_errRp_ct` is a static helper computing per-iter errRp over the cell table. Its only definition site lives in raking.cpp. Recent refactors may have replaced inline call sites without removing the helper.
**Tech Stack:** C++17, R CMD INSTALL --preclean ., grep.

**Mechanism:** grep audit; if dead, surgical deletion.
**Forbidden:** removing other unused statics in raking.cpp not flagged by this ticket; refactoring the inline replacement at the call sites; changing function signatures of any retained code.

---

## Task T1: Caller audit

Steps:

1. Grep across all C++ source and headers:

```bash
grep -n "compute_errRp_ct" src/*.cpp src/*.hpp src/*.h
```

2. Initial finding (already verified at plan time): only matches are at `src/raking.cpp:34` (definition) — no callers.

3. **Halt and decide** based on the audit:
   - If 0 non-definition matches: proceed to T2 (delete).
   - If ≥1 caller: close the ticket as "still in use" with a comment showing the call sites.

Confidence: 95 — already grepped during plan write-up; only the definition appears.

---

## Task T2: Delete

Edit `src/raking.cpp`, remove lines 33-50 (the helper plus the trailing blank line):

```cpp
// Cell-table errRp: O(K * M_cell). bucket pre-allocated to max_cats.
static double compute_errRp_ct(const CalibState& st,
                                const CellTable& ct,
                                const std::vector<double>& X,
                                std::vector<double>& bucket) {
    double W = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W += X[c];
    double err = 0.0;
    for (int k = 0; k < st.K; k++) {
        lbw::aggregate_to_margin(ct, X, k, st.cat_counts[k], bucket.data());
        for (int j = 0; j < st.cat_counts[k]; j++) {
            double e = std::fabs(bucket[j] / W - st.targets[k][j]);
            if (e > err) err = e;
        }
    }
    return err;
}
```

Confidence: 95.

---

## Task T3: Rebuild + test

```bash
R CMD INSTALL --preclean .
Rscript -e 'devtools::test(filter = "raking")'
Rscript -e 'devtools::test(filter = "calibration-solvers")'
```

Pass criteria: all raking and calibration-solver tests pass; no compiler warning about an unused function (-Wunused-function).

---

## Task T4: Commit

`refactor(raking): remove dead compute_errRp_ct helper`

Body: cites grep showing zero callers (leafblower-9jmj).
