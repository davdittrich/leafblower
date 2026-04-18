# anatomy.md

> Auto-maintained by OpenWolf. Last scanned: 2026-04-18T12:46:12.093Z
> Files: 34 tracked | Anatomy hits: 0 | Misses: 0

## ./

- `.gitignore` — Git ignore rules (~20 tok)
- `.Rbuildignore` (~27 tok)
- `AGENTS.md` — Agent Instructions (~726 tok)
- `CLAUDE.md` — OpenWolf (~564 tok)
- `configure` — Detect C++17 support. Note: -c is required; linking /dev/null fails without main(). (~100 tok)
- `DESCRIPTION` (~172 tok)
- `LICENSE` (~290 tok)
- `NAMESPACE` (~54 tok)

## .beads/

- `.gitignore` — Git ignore rules (~444 tok)
- `.local_version` (~2 tok)
- `config.yaml` — Beads Configuration File (~597 tok)
- `interactions.jsonl` (~0 tok)
- `metadata.json` (~46 tok)
- `README.md` — Project documentation (~562 tok)

## .beads/embeddeddolt/

- `.lock` (~0 tok)

## .beads/embeddeddolt/leafblower/.dolt/

- `config.json` (~1 tok)
- `repo_state.json` (~24 tok)

## .beads/embeddeddolt/leafblower/.dolt/noms/

- `journal.idx` (~3967 tok)
- `LOCK` (~0 tok)
- `manifest` (~39 tok)
- `vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv` (~62644 tok)

## .beads/hooks/

- `post-checkout` — --- BEGIN BEADS INTEGRATION v1.0.0 --- (~213 tok)
- `post-merge` — --- BEGIN BEADS INTEGRATION v1.0.0 --- (~210 tok)
- `pre-commit` — --- BEGIN BEADS INTEGRATION v1.0.0 --- (~210 tok)
- `pre-push` — --- BEGIN BEADS INTEGRATION v1.0.0 --- (~208 tok)
- `prepare-commit-msg` — --- BEGIN BEADS INTEGRATION v1.0.0 --- (~218 tok)

## .claude/

- `settings.json` (~539 tok)

## .claude/rules/

- `openwolf.md` (~313 tok)

## R/

- `zzz.R` — R_init_leafblower() in r_bridge.cpp is called automatically by R when the (~86 tok)

## docs/superpowers/plans/

- `2026-04-18-leafblower-core.md` — Leafblower Core — Implementation Plan (~24756 tok)

## src/

- `leafblower.h` — ifndef LEAFBLOWER_H (~767 tok)
- `Makevars.in` (~36 tok)
- `types.hpp` — lbw::CalibState struct; internal state for calibration algorithms (~80 tok)
- `c_api.cpp` — rk_params_init defaults, validate_inputs, rk_calibrate stub; extern "C" (~220 tok)

## tasks/

- `prd-leafblower-core.md` — Leafblower Core — Product Requirements Document (~11700 tok)

## tests/

- `testthat.R` — This file is part of the R package leafblower. (~49 tok)

## tests/testthat/

- `test-harvest.R` — Placeholder; real BADARG RED tests deferred to Task 7 (~20 tok)
