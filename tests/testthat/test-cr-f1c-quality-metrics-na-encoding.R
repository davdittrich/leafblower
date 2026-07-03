# CR-F1c (dtkn.13): compute_quality_metrics recoded NA observations into an "NA"
# bin whenever the target merely CONTAINED an "NA" level name — even for a HAND-BUILT
# "NA" target with add_na_proportion=FALSE, where the solver maps real NAs to gid -1
# (EXCLUDED). That reported a KL for an encoding the solver never used. The metric now
# keys on na_margins (the actual injected-NA set the solver used, harvest.R:352/554),
# so it mirrors the solver: only injected margins recode; a hand-built "NA" target
# falls through to the valid-mask path (NA excluded) and, being unmatched, surfaces
# as Inf — the honest signal that the solver did not fit it.

cqm <- leafblower:::compute_quality_metrics

test_that("injected NA margin (in na_margins) recodes NAs into the 'NA' bin → finite margin_kl (CR-F1c)", {
  set.seed(1); n <- 60
  a <- sample(c("x", "y"), n, TRUE)
  a[1:12] <- NA  # 20% NA — the add_na_proportion injection scenario
  df <- data.frame(a = a, b = factor(sample(c("p", "q"), n, TRUE)))
  na_frac <- mean(is.na(df$a))
  tg <- list(a = c(x = 0.5 * (1 - na_frac), y = 0.5 * (1 - na_frac), "NA" = na_frac),
             b = c(p = 0.5, q = 0.5))
  w <- rep(1, n)
  # Solver recoded 'a' (v %in% .na_margins) → NAs go to the "NA" bin. Metric must match.
  mk <- cqm(w, tg, df, na_margins = "a")$margin_kl
  expect_true(is.finite(mk))
})

test_that("hand-built 'NA' target NOT in na_margins excludes NAs (matches solver gid -1) → Inf (CR-F1c)", {
  set.seed(1); n <- 60
  a <- sample(c("x", "y"), n, TRUE)
  a[1:12] <- NA
  df <- data.frame(a = a, b = factor(sample(c("p", "q"), n, TRUE)))
  # Same target shape, but the user hand-built "NA" WITHOUT add_na_proportion, so the
  # solver never injected it (na_margins empty) and maps NAs to gid -1 (excluded).
  tg <- list(a = c(x = 0.3, y = 0.3, "NA" = 0.4), b = c(p = 0.5, q = 0.5))
  w <- rep(1, n)
  mk <- cqm(w, tg, df, na_margins = character(0))$margin_kl
  # NAs excluded → the positive-mass "NA" target level is unmatched → Inf (honest),
  # NOT the previously-reported finite KL from an encoding the solver never used.
  expect_true(is.infinite(mk))
})

test_that("no 'NA' target level → metric unchanged and finite regardless of na_margins (CR-F1c)", {
  set.seed(3); n <- 60
  df <- data.frame(a = factor(sample(c("x", "y"), n, TRUE)),
                   b = factor(sample(c("p", "q"), n, TRUE)))
  tg <- list(a = c(x = 0.5, y = 0.5), b = c(p = 0.5, q = 0.5))
  w <- rep(1, n)
  expect_equal(cqm(w, tg, df, na_margins = character(0))$margin_kl,
               cqm(w, tg, df, na_margins = "a")$margin_kl)
})
