library(leafblower)

# T-D: Strategy B orthogonality validator tests.
# Tests lbw::validate_orthogonal_split via:
#   (a) C_hier_orthogonal_probe — direct C++ unit tests
#   (b) harvest() end-to-end — Guard (13) fires on BADARG, passes on orthogonal

# ── Probe helper ─────────────────────────────────────────────────────────────
# group_ids: list of K integer vectors (0-indexed, -1=NA), each length n
# coarse_mask: integer vector length K (1=coarse, 0=fine)
# cat_counts: integer vector length K
.orthogonal_probe <- function(group_ids, n, K, coarse_mask, cat_counts) {
  .Call("C_hier_orthogonal_probe",
        group_ids,
        as.integer(n),
        as.integer(K),
        as.integer(coarse_mask),
        as.integer(cat_counts),
        PACKAGE = "leafblower")
}

RK_OK        <- 0L
RK_ERR_BADARG <- 3L

# ── Direct unit tests via probe ───────────────────────────────────────────────

# (1a) Orthogonal split: age_bracket (coarse) vs age_year (fine) nested.
# age_bracket 0 contains age_year 0,1,2  (young: 0-2)
# age_bracket 1 contains age_year 3,4,5  (old:   3-5)
# No fine level spans two coarse cells -> RK_OK.
test_that("orthogonal nested split: fine levels nested inside coarse -> RK_OK", {
  n <- 120L
  set.seed(42)
  # age_bracket: 0=young (obs 0..59), 1=old (obs 60..119)
  # age_year:    0,1,2 within young; 3,4,5 within old — strict nesting
  age_bracket <- c(rep(0L, 60L), rep(1L, 60L))
  age_year    <- c(rep(c(0L, 1L, 2L), 20L), rep(c(3L, 4L, 5L), 20L))

  res <- .orthogonal_probe(
    group_ids   = list(age_bracket, age_year),
    n           = n,
    K           = 2L,
    coarse_mask = c(1L, 0L),
    cat_counts  = c(2L, 6L)
  )
  expect_equal(res$rc, RK_OK)
  expect_identical(res$diagnostic, "")
})

# (1b) Non-orthogonal split: fine margin 'region' (fine) spans two coarse cells.
# coarse = age_group (0=young, 1=old); fine = region (0=north, 1=south).
# Some north obs are in both age groups -> region level 0 spans cells {0, 1}.
test_that("non-orthogonal split: fine margin level spans 2 coarse cells -> RK_ERR_BADARG with diagnostic", {
  n <- 100L
  # Mix: first 50 are young (age=0), last 50 are old (age=1).
  # Region interleaved: 0,1,0,1,... -> every region level appears in both coarse cells.
  age    <- c(rep(0L, 50L), rep(1L, 50L))
  region <- rep(c(0L, 1L), 50L)

  res <- .orthogonal_probe(
    group_ids   = list(age, region),
    n           = n,
    K           = 2L,
    coarse_mask = c(1L, 0L),
    cat_counts  = c(2L, 2L)
  )
  expect_equal(res$rc, RK_ERR_BADARG)
  # Diagnostic must name margin index, level, and coarse cells.
  expect_match(res$diagnostic, "Fine margin")
  expect_match(res$diagnostic, "coarse cells")
  expect_match(res$diagnostic, "orthogonal split")
})

# (1c) K=3 margins: two coarse, one fine — fine perfectly nested.
test_that("K=3 with two coarse margins, one fine nested margin -> RK_OK", {
  n <- 60L
  # coarse: gender (0/1), region (0/1/2)
  # fine:   age_year nested: gender=0,region=0 -> age=0; gender=0,region=1 -> age=1; etc.
  gender  <- c(rep(0L, 30L), rep(1L, 30L))
  region  <- c(rep(0L, 10L), rep(1L, 10L), rep(2L, 10L),
               rep(0L, 10L), rep(1L, 10L), rep(2L, 10L))
  age_fine <- c(rep(0L, 10L), rep(1L, 10L), rep(2L, 10L),
                rep(3L, 10L), rep(4L, 10L), rep(5L, 10L))

  res <- .orthogonal_probe(
    group_ids   = list(gender, region, age_fine),
    n           = n,
    K           = 3L,
    coarse_mask = c(1L, 1L, 0L),
    cat_counts  = c(2L, 3L, 6L)
  )
  expect_equal(res$rc, RK_OK)
  expect_identical(res$diagnostic, "")
})

# (1d) Diagnostic format check: names first offending margin, level, cell IDs.
test_that("diagnostic names offending fine margin index, level value, and coarse cells", {
  n <- 60L
  # Margin 0 is coarse (2 levels). Margin 1 is fine but level 0 spans BOTH coarse cells.
  coarse <- c(rep(0L, 30L), rep(1L, 30L))
  fine   <- rep(0L, 60L)  # single level 0 in every obs -> spans coarse cells {0, 1}

  res <- .orthogonal_probe(
    group_ids   = list(coarse, fine),
    n           = n,
    K           = 2L,
    coarse_mask = c(1L, 0L),
    cat_counts  = c(2L, 1L)
  )
  expect_equal(res$rc, RK_ERR_BADARG)
  # Must mention margin 1 (fine) and level 0 and multiple cell IDs.
  expect_match(res$diagnostic, "margin 1")
  expect_match(res$diagnostic, "level 0")
  expect_match(res$diagnostic, "\\{")   # opening brace of cell set
})

# ── End-to-end via harvest() ──────────────────────────────────────────────────
# Shared helper: two-margin setup K=2.
# Orthogonal: age_bracket (coarse, 2 levels), age_year (fine, 4 levels nested).
#   young -> years 0,1; old -> years 2,3. No overlap.
make_orthogonal_hier <- function(n = 200L, method = "raking") {
  set.seed(7)
  age_bracket_vec <- rep(c("young", "old"), each = n / 2L)
  age_year_vec    <- c(
    sample(c("y0", "y1"), n / 2L, replace = TRUE),  # young uses y0, y1
    sample(c("y2", "y3"), n / 2L, replace = TRUE)   # old   uses y2, y3
  )
  df  <- data.frame(
    age_bracket = factor(age_bracket_vec, levels = c("young", "old")),
    age_year    = factor(age_year_vec,    levels = c("y0", "y1", "y2", "y3"))
  )
  tgt <- list(
    age_bracket = c(young = 0.5, old = 0.5),
    age_year    = c(y0 = 0.25, y1 = 0.25, y2 = 0.25, y3 = 0.25)
  )
  hier <- list(
    coarse_mask       = c(1L, 0L),   # age_bracket coarse, age_year fine
    min_cell_n        = 5L,
    mode              = 1L,           # exact — Strategy B
    outer_tol         = 1e-3,
    outer_iterations  = 10L
  )
  list(df = df, tgt = tgt, hier = hier, method = method)
}

# Non-orthogonal: region (fine) spans BOTH age groups.
make_nonorthogonal_hier <- function(n = 200L, method = "raking") {
  set.seed(8)
  age_bracket_vec <- rep(c("young", "old"), each = n / 2L)
  region_vec      <- sample(c("north", "south"), n, replace = TRUE)  # crosses coarse cells
  df  <- data.frame(
    age_bracket = factor(age_bracket_vec, levels = c("young", "old")),
    region      = factor(region_vec,      levels = c("north", "south"))
  )
  tgt <- list(
    age_bracket = c(young = 0.5, old = 0.5),
    region      = c(north = 0.5, south = 0.5)
  )
  hier <- list(
    coarse_mask       = c(1L, 0L),   # age_bracket coarse, region fine — region crosses!
    min_cell_n        = 5L,
    mode              = 1L,           # exact
    outer_tol         = 1e-3,
    outer_iterations  = 10L
  )
  list(df = df, tgt = tgt, hier = hier, method = method)
}

# (2a) mode="exact" with orthogonal split AND P1 method (raking) -> RK_OK.
test_that("mode=exact + orthogonal split + raking (P1) -> succeeds", {
  s <- make_orthogonal_hier(method = "raking")
  expect_no_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier)
  )
})

# (2b) mode="exact" with non-orthogonal split -> BADARG via harvest() stop().
test_that("mode=exact + non-orthogonal split -> error naming orthogonal split", {
  s <- make_nonorthogonal_hier(method = "raking")
  expect_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier),
    regexp = "orthogonal split"
  )
})

# (2c) mode="refine" with same non-orthogonal split -> NO orthogonality error
# (Guard 13 fires only when mode==1). LOGIT+refine is already BADARG from Guard 10;
# use raking+refine to confirm Guard 13 is silent.
test_that("mode=refine + non-orthogonal split + raking -> no Guard-13 error", {
  s <- make_nonorthogonal_hier(method = "raking")
  s$hier$mode <- 0L   # refine — Guard 13 must NOT fire
  expect_no_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier)
  )
})

# (2d) mode="exact" + orthogonal split + LOGIT (P2) -> succeeds (Guard 10 only
# rejects LOGIT+refine, not LOGIT+exact).
test_that("mode=exact + orthogonal split + logit (P2) -> succeeds", {
  s <- make_orthogonal_hier(method = "logit")
  expect_no_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier)
  )
})
