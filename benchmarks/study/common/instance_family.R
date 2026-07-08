# instance_family.R -- parametric synthetic instance-family generator.
#
# WU-3 (leafblower-2ouc.4). Design of record: docs/benchmark/DESIGN.md §3C
# (instance-family axes), §7 (Dolan-More profiles need a problem
# DISTRIBUTION), §11 DoD. Emits a frozen 30-100-instance sweep over
# n / K / cell-cardinality / target-skew / margin-collinearity /
# infeasibility-slack, each materialized lazily through a `gen:instance_family`
# data_ref resolved by the WU-2 loader (common/problem_io.R).
#
# STRICT SEPARATION (user constraint 2026-07-08): this file calls NOTHING in
# leafblower. It only emits standardized problem specs (data.frame + targets
# + bounds), consumed downstream by adapters that themselves call
# leafblower::harvest() -- out of this file's scope entirely.
#
# Non-home-turf rationale (DESIGN.md §3C "designed independently of
# leafblower's known convergence behavior"): the recipe below is a generic
# categorical-margin generator with two knobs borrowed from standard,
# citable techniques outside the calibration literature --
#   (1) Zipf/power-law category-frequency skew (theta exponent), the
#       standard technique for tunable-skew synthetic benchmarks (e.g.
#       sysbench's Zipfian distribution for workload generation) -- NOT
#       derived from any leafblower fixture or observed solver behavior;
#   (2) a copy-vs-independent-draw collinearity knob (rho) that induces
#       near-rank-deficient margin combinations, a standard synthetic
#       ill-conditioning construction. A targeted web search (2026-07-09,
#       WU-3 authoring) found no single published benchmark suite combining
#       raking-style margin calibration + tunable Zipf skew + a reproducible
#       instance generator (closest: sysbench's theta-parameterized Zipfian
#       distribution for the skew mechanism; no survey-calibration-specific
#       analog) -- confirming DESIGN.md's in-house fallback applies here.
#
# DETERMINISM / R<->Python PARITY (DoD: "identical ... exactly if
# deterministic construction"): row-level category draws use a hand-rolled
# Lehmer/Park-Miller multiplicative LCG (modulus m = 2^31-1, multiplier
# a = 48271) instead of R's/numpy's native RNG streams, which do NOT agree
# bit-for-bit across languages. Every LCG operation (multiply, %%, /,
# comparison) is IEEE-754 double arithmetic on operands that stay exactly
# representable in a double's 53-bit mantissa (max product
# 48271 * (2^31-2) ~= 1.037e14 << 2^53 ~= 9.007e15), so R and Python compute
# BIT-IDENTICAL uniforms for identical seeds -- exact parity, not rtol-based.
#
# Two-step loader wiring (no problem_io.R edit -- out of WU-3 file scope):
# `install_gen_resolver()` monkey-patches the *already-sourced* problem_io.R
# script-level function `.pio_resolve_data_ref` in the shared sourcing
# environment (both files are `source()`d as plain scripts into the same
# environment; R resolves free variables in a function body by searching its
# *enclosing* environment at CALL time, not at definition time, so
# reassigning `.pio_resolve_data_ref` after sourcing problem_io.R changes
# what `load_problem_spec()` invokes). The patch intercepts ONLY
# `gen:instance_family...` data_refs it recognizes; every other `gen:<id>`
# (including WU-2's own frozen `gen:toy_recipe` regression guard in
# test_problem_io.R) falls through unchanged to the original resolver, which
# still raises its WU-3-not-yet-implemented error.

suppressPackageStartupMessages({
  library(jsonlite)
})

# ---- LCG core --------------------------------------------------------

.if_lcg_m <- 2147483647  # 2^31 - 1 (Mersenne prime, Park-Miller modulus)
.if_lcg_a <- 48271       # Park-Miller 1993 revised multiplier

.if_lcg_next <- function(state) {
  (.if_lcg_a * state) %% .if_lcg_m
}

# Deterministic per-instance seed from (master_seed, instance-local seed).
# EXACTNESS NOTE: master_seed (~2.03e7) * 2654435761 (~2.65e9) ~= 5.38e16,
# which EXCEEDS a double's 53-bit-mantissa exact range (2^53 ~= 9.007e15) --
# R's `*` on doubles would round that product, while Python's arbitrary-
# precision `int` would not, silently breaking R<->Python bit-parity (caught
# by the roundtrip parity test). Fix: assign the SMALL multiplier (40503) to
# the large, run-constant master_seed (product ~8.2e11, exact) and the LARGE
# multiplier (2654435761) to the small, per-instance local_seed (product
# <= 2.65e11 for local_seed <= 100 as used by build_instance_grid(), exact)
# -- both operands' products individually stay << 2^53.
.if_lcg_seed0 <- function(master_seed, local_seed) {
  raw <- (master_seed * 40503 + local_seed * 2654435761 + 1) %% (.if_lcg_m - 1)
  raw + 1  # avoid the absorbing state 0
}

# ---- Zipf target distribution -----------------------------------------

# Category j=1..C weight j^(-s); s=0 -> uniform. Cumulative probs, last
# element forced to exactly 1 (avoids float residue leaving cat C
# unreachable by inverse-CDF comparisons).
.if_zipf_probs <- function(C, s) {
  j <- seq_len(C)
  raw <- if (s == 0) rep(1, C) else j^(-s)
  probs <- raw / sum(raw)
  cum <- cumsum(probs)
  cum[C] <- 1
  list(probs = probs, cum = cum)
}

.if_inv_cdf <- function(u, cum) {
  C <- length(cum)
  for (j in seq_len(C)) {
    if (u <= cum[j]) return(j)
  }
  C
}

# ---- Row-level category generation -------------------------------------

# Generates an n x K integer category matrix (values 1..C) for one instance.
# Margin 1 is drawn from the Zipf(C,s) distribution every row. Each margin
# k=2..K is, with probability rho, a deterministic shift of margin 1's
# category (collinearity knob); otherwise an independent Zipf(C,s) draw.
.if_generate_categories <- function(n, K, C, s, rho, seed) {
  zp <- .if_zipf_probs(C, s)
  cum <- zp$cum
  state <- .if_lcg_seed0(20260708L, seed)
  mat <- matrix(0L, nrow = n, ncol = K)
  for (i in seq_len(n)) {
    state <- .if_lcg_next(state)
    cat1 <- .if_inv_cdf(state / .if_lcg_m, cum)
    mat[i, 1L] <- cat1
    if (K >= 2L) {
      for (k in 2:K) {
        state <- .if_lcg_next(state)
        u_dep <- state / .if_lcg_m
        if (u_dep <= rho) {
          mat[i, k] <- ((cat1 - 1L + (k - 1L)) %% C) + 1L
        } else {
          state <- .if_lcg_next(state)
          mat[i, k] <- .if_inv_cdf(state / .if_lcg_m, cum)
        }
      }
    }
  }
  list(mat = mat, zipf = zp)
}

# ---- Axis level tables (frozen) ----------------------------------------

.IF_CARD_LEVELS <- c(low = 4L, medium = 8L, high = 16L)
.IF_SKEW_LEVELS <- c(none = 0, moderate = 1, extreme = 3)
.IF_COND_LEVELS <- c(well = 0.0, moderate = 0.5, ill = 0.9)
.IF_INFEAS_MAX <- list(loose = NULL, moderate = 3.0, tight = 1.05)
.IF_N_LEVELS <- c(1000L, 10000L, 100000L, 1580000L)
.IF_K_LEVELS <- c(2L, 4L, 9L)
.IF_MASTER_SEED <- 20260708L

# ---- Instance grid (frozen recipe; deterministic, no RNG) --------------

.if_add <- function(lst, n, K, card, skew, cond, infeas, seed, tag) {
  lst[[length(lst) + 1L]] <- list(
    n = n, K = K, card = card, skew = skew, cond = cond, infeas = infeas,
    seed = seed, tag = tag
  )
  lst
}

build_instance_grid <- function() {
  g <- list()
  # A: primary n x K grid at baseline card/skew/cond/infeas. 1.58M capped
  # to K=9 only ("only a handful at 1.58M", DESIGN.md §4).
  for (n in .IF_N_LEVELS) {
    for (K in .IF_K_LEVELS) {
      if (n == 1580000L && K != 9L) next
      g <- .if_add(g, n, K, "medium", "moderate", "well", "loose", 0L, "primary")
    }
  }
  # B: OFAT stress on each difficulty axis at baseline n=10000, K=4.
  for (lvl in c("low", "high")) g <- .if_add(g, 10000L, 4L, lvl, "moderate", "well", "loose", 0L, "ofat_card")
  for (lvl in c("none", "extreme")) g <- .if_add(g, 10000L, 4L, "medium", lvl, "well", "loose", 0L, "ofat_skew")
  for (lvl in c("moderate", "ill")) g <- .if_add(g, 10000L, 4L, "medium", "moderate", lvl, "loose", 0L, "ofat_cond")
  for (lvl in c("moderate", "tight")) g <- .if_add(g, 10000L, 4L, "medium", "moderate", "well", lvl, 0L, "ofat_infeas")
  # C: combined worst-case stress (all axes extreme) at a few n levels.
  for (n in c(10000L, 100000L, 1580000L)) g <- .if_add(g, n, 9L, "high", "extreme", "ill", "tight", 0L, "stress")
  g <- .if_add(g, 1000L, 2L, "high", "extreme", "ill", "tight", 0L, "stress")
  # D: pairwise axis interactions at baseline n=10000, K=4.
  g <- .if_add(g, 10000L, 4L, "high", "moderate", "ill", "loose", 0L, "pairwise")
  g <- .if_add(g, 10000L, 4L, "medium", "extreme", "ill", "loose", 0L, "pairwise")
  g <- .if_add(g, 10000L, 4L, "high", "extreme", "well", "loose", 0L, "pairwise")
  g <- .if_add(g, 10000L, 4L, "high", "moderate", "well", "tight", 0L, "pairwise")
  # E: seed replicates for profile density (same structural params, fresh draws).
  for (seed in c(1L, 2L)) g <- .if_add(g, 10000L, 4L, "medium", "moderate", "well", "loose", seed, "replicate")
  for (seed in c(1L, 2L)) g <- .if_add(g, 1000L, 4L, "medium", "moderate", "well", "loose", seed, "replicate")
  for (seed in c(1L, 2L)) g <- .if_add(g, 10000L, 4L, "high", "moderate", "well", "loose", seed, "replicate")
  g
}

.if_instance_id <- function(inst) {
  sprintf("if_n%d_K%d_%s_%s_%s_%s_s%d", inst$n, inst$K, inst$card, inst$skew,
          inst$cond, inst$infeas, inst$seed)
}

.if_data_ref <- function(inst) {
  sprintf("gen:instance_family?n=%d&K=%d&card=%s&skew=%s&cond=%s&infeas=%s&seed=%d",
          inst$n, inst$K, inst$card, inst$skew, inst$cond, inst$infeas, inst$seed)
}

# ---- Spec materialization (schema.json-conformant) ----------------------

.if_build_spec <- function(inst) {
  K <- inst$K
  C <- unname(.IF_CARD_LEVELS[[inst$card]])
  s <- unname(.IF_SKEW_LEVELS[[inst$skew]])
  zp <- .if_zipf_probs(C, s)
  cats <- paste0("c", seq_len(C))
  target_one <- as.list(setNames(zp$probs, cats))
  margins <- paste0("m", seq_len(K))
  targets <- setNames(rep(list(target_one), K), margins)
  bmax <- .IF_INFEAS_MAX[[inst$infeas]]
  families <- c("kl", "chi2", "logit", "newton_kl", "minimax")
  if (K == 2L) families <- c(families, "ot")
  list(
    id = .if_instance_id(inst),
    data_ref = .if_data_ref(inst),
    margins = I(margins),
    targets = targets,
    bounds = list(min = 0, max = bmax),
    tol = 1e-8,
    objective_families = I(families),
    K = K
  )
}

# Writes one spec/instance_family/<id>.json per instance plus the grid
# manifest spec/instance_family.json (axis levels + instance roster).
# Cheap regardless of n: data is NOT materialized here (gen: is lazy).
generate_instance_family_specs <- function(out_dir = file.path("benchmarks", "study", "spec")) {
  grid <- build_instance_grid()
  inst_dir <- file.path(out_dir, "instance_family")
  if (!dir.exists(inst_dir)) dir.create(inst_dir, recursive = TRUE)
  ids <- character(length(grid))
  for (i in seq_along(grid)) {
    spec <- .if_build_spec(grid[[i]])
    ids[i] <- spec$id
    path <- file.path(inst_dir, paste0(spec$id, ".json"))
    writeLines(jsonlite::toJSON(spec, auto_unbox = TRUE, digits = 12, null = "null"), path)
  }
  manifest <- list(
    recipe = "instance_family",
    master_seed = .IF_MASTER_SEED,
    n_instances = length(grid),
    axis_levels = list(
      n = I(as.list(.IF_N_LEVELS)),
      K = I(as.list(.IF_K_LEVELS)),
      card = as.list(.IF_CARD_LEVELS),
      skew = as.list(.IF_SKEW_LEVELS),
      cond = as.list(.IF_COND_LEVELS),
      infeas = list(loose = "Inf", moderate = 3.0, tight = 1.05)
    ),
    instances = I(ids),
    frozen = TRUE
  )
  writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, digits = 12, null = "null"),
             file.path(out_dir, "instance_family.json"))
  list(ids = ids, out_dir = inst_dir)
}

# ---- gen: origin resolution + loader wiring -----------------------------

.if_parse_query <- function(qs) {
  parts <- strsplit(qs, "&", fixed = TRUE)[[1]]
  kv <- lapply(parts, function(p) strsplit(p, "=", fixed = TRUE)[[1]])
  vals <- setNames(vapply(kv, `[`, character(1), 2), vapply(kv, `[`, character(1), 1))
  as.list(vals)
}

# Resolves a `gen:instance_family?...` data_ref into a data.frame with
# columns m1..mK (factor-ready character), matching problem_io's expected
# data-frame shape. Returns NULL if the recipe id isn't "instance_family"
# (caller falls through to the original resolver for unrecognized recipes).
.if_resolve_gen_data <- function(data_ref) {
  recipe_and_query <- sub("^gen:", "", data_ref)
  parts <- strsplit(recipe_and_query, "?", fixed = TRUE)[[1]]
  recipe_id <- parts[1]
  if (!identical(recipe_id, "instance_family")) return(NULL)
  q <- .if_parse_query(parts[2])
  n <- as.integer(q$n)
  K <- as.integer(q$K)
  C <- unname(.IF_CARD_LEVELS[[q$card]])
  s <- unname(.IF_SKEW_LEVELS[[q$skew]])
  rho <- unname(.IF_COND_LEVELS[[q$cond]])
  seed <- as.integer(q$seed)
  gen <- .if_generate_categories(n, K, C, s, rho, seed)
  df <- as.data.frame(gen$mat)
  names(df) <- paste0("m", seq_len(K))
  for (nm in names(df)) df[[nm]] <- paste0("c", df[[nm]])
  df
}

.if_installed <- FALSE

# Monkey-patches problem_io.R's script-level `.pio_resolve_data_ref` (must
# already be sourced) to dispatch `gen:instance_family?...` to this file's
# generator, delegating every other data_ref (incl. unrecognized `gen:<id>`)
# to the original function. Idempotent; safe to call more than once.
install_gen_resolver <- function() {
  if (.if_installed) return(invisible(NULL))
  if (!exists(".pio_resolve_data_ref", mode = "function")) {
    stop("instance_family: problem_io.R must be source()d before install_gen_resolver()",
         call. = FALSE)
  }
  original <- get(".pio_resolve_data_ref", mode = "function")
  patched <- function(data_ref, inline_data) {
    if (grepl("^gen:instance_family", data_ref)) {
      out <- .if_resolve_gen_data(data_ref)
      if (!is.null(out)) return(out)
    }
    original(data_ref, inline_data)
  }
  # Install into the same environment problem_io.R's functions were defined
  # in (the caller's -- typically globalenv() when both are top-level
  # source()d), so load_problem_spec()'s free-variable lookup finds it.
  assign(".pio_resolve_data_ref", patched, envir = environment(original))
  .if_installed <<- TRUE
  invisible(NULL)
}
