# benchmarks/plot_helpers.R
# Shared plotting helpers for 2D and 3D Bayesian benchmarks.

# ── make_plots (2D) ───────────────────────────────────────────────────────────
# Generates contour and uncertainty plots for 2D design space.
# Args:
#   state: list with elements $gp (km object), $design (n×2 matrix), $y (numeric)
#   candidates: m×2 matrix of candidate grid points
#   threshold: scalar decision boundary for contour line
#   out_dir: directory for output PDFs (default "benchmarks")
# Returns: invisibly, list of ggplot2 objects (contour, uncertainty)
make_plots <- function(state, candidates, threshold, out_dir = "benchmarks") {
  pred    <- DiceKriging::predict(state$gp,
                                  newdata    = as.data.frame(candidates),
                                  type       = "UK",
                                  checkNames = FALSE)
  cand_df <- as.data.frame(candidates)
  names(cand_df) <- c("log_complexity", "log_tol")
  cand_df$mean <- pred$mean
  cand_df$sd   <- pred$sd

  # Mark evaluated points: circle = LHC initial, cross = adaptive
  n_lhc    <- 8L
  n_total  <- nrow(state$design)
  pt_df    <- as.data.frame(state$design)
  names(pt_df) <- c("log_complexity", "log_tol")
  pt_df$type <- c(rep("LHC", min(n_lhc, n_total)),
                  rep("Adaptive", max(0L, n_total - n_lhc)))

  # ── Contour plot ───────────────────────────────────────────────────────────
  p1 <- ggplot2::ggplot(cand_df, ggplot2::aes(log_complexity, log_tol)) +
    ggplot2::geom_tile(ggplot2::aes(fill = mean)) +
    ggplot2::geom_contour(ggplot2::aes(z = mean),
                          breaks    = threshold,
                          colour    = "red",
                          linewidth = 1.2) +
    ggplot2::geom_point(data  = pt_df,
                        ggplot2::aes(shape = type),
                        colour = "white", size = 2) +
    ggplot2::scale_fill_viridis_c(name = "log(t_ORIS/baseline)") +
    ggplot2::scale_shape_manual(values = c(LHC = 16L, Adaptive = 4L)) +
    ggplot2::labs(
      title   = "GP posterior mean - red line = 1.2x contour",
      x       = "log10(complexity = n x sum(cat_counts))",
      y       = "log10(tol_abs)",
      shape   = "Design point"
    ) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(out_dir, "algo_selection_contour.pdf"), p1, width = 9, height = 6)

  # ── Uncertainty plot ───────────────────────────────────────────────────────
  p2 <- ggplot2::ggplot(cand_df, ggplot2::aes(log_complexity, log_tol)) +
    ggplot2::geom_tile(ggplot2::aes(fill = sd)) +
    ggplot2::scale_fill_viridis_c(name = "posterior sd", option = "magma") +
    ggplot2::labs(
      title = "GP posterior uncertainty - high sd = unreliable region",
      x     = "log10(complexity = n x sum(cat_counts))",
      y     = "log10(tol_abs)"
    ) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(out_dir, "algo_selection_uncertainty.pdf"), p2, width = 9, height = 6)

  cat(sprintf("Plots saved to %s/\n", out_dir))
  invisible(list(contour = p1, uncertainty = p2))
}

# ── make_plots_3d ─────────────────────────────────────────────────────────────
# Generates slice contour plots for 3D design space via 2D projections.
# Args:
#   state: list with elements $gp (km object), $design (n×3 matrix), $y (numeric),
#          $bounds (list of 3 numeric ranges [c(min, max), ...])
#   candidates: (unused, for API compatibility)
#   threshold: scalar decision boundary for contour line
#   out_dir: directory for output PDFs (default "benchmarks")
#   slice_dim: which dimension to slice (default 3)
#   slice_values: which values of slice_dim to use (default NULL → c(0, 1, 2, 3))
# Returns: invisibly, list of ggplot2 objects (one per slice)
make_plots_3d <- function(state, candidates, threshold, out_dir = "benchmarks",
                          slice_dim = 3, slice_values = NULL) {
  if (is.null(slice_values)) {
    slice_values <- c(0, 1, 2, 3)
  }

  # Map slice dimension to kept dimensions
  # slice_dim=1 → fix x1, vary x2,x3
  # slice_dim=2 → fix x2, vary x1,x3
  # slice_dim=3 → fix x3, vary x1,x2
  if (!(slice_dim %in% c(1, 2, 3))) {
    stop("slice_dim must be 1, 2, or 3")
  }
  kept_dims <- setdiff(1:3, slice_dim)

  plots <- list()

  for (val in slice_values) {
    # Construct 50×50 grid in the kept dimensions
    rng1 <- state$bounds[[kept_dims[1]]]
    rng2 <- state$bounds[[kept_dims[2]]]

    grid_1d_1 <- seq(rng1[1], rng1[2], length.out = 50L)
    grid_1d_2 <- seq(rng2[1], rng2[2], length.out = 50L)
    grid_base <- expand.grid(g1 = grid_1d_1, g2 = grid_1d_2)

    # Expand to 3D by inserting slice_value at slice_dim
    grid_3d <- matrix(NA_real_, nrow = nrow(grid_base), ncol = 3)
    grid_3d[, kept_dims[1]] <- grid_base$g1
    grid_3d[, kept_dims[2]] <- grid_base$g2
    grid_3d[, slice_dim] <- val

    # Predict on grid
    pred <- DiceKriging::predict(state$gp,
                                 newdata    = as.data.frame(grid_3d),
                                 type       = "UK",
                                 checkNames = FALSE)

    # Build data frame for plotting
    cand_df <- as.data.frame(grid_base)
    names(cand_df) <- c("x1", "x2")
    cand_df$mean <- pred$mean
    cand_df$sd   <- pred$sd

    # Dimension labels
    dim_labels <- c("log10_complexity", "log10_tol", "log10_compression")
    x1_label <- dim_labels[kept_dims[1]]
    x2_label <- dim_labels[kept_dims[2]]
    slice_label <- dim_labels[slice_dim]

    # Create contour plot
    p <- ggplot2::ggplot(cand_df, ggplot2::aes(x1, x2)) +
      ggplot2::geom_tile(ggplot2::aes(fill = mean)) +
      ggplot2::geom_contour(ggplot2::aes(z = mean),
                            breaks    = threshold,
                            colour    = "red",
                            linewidth = 1.2) +
      ggplot2::scale_fill_viridis_c(name = "log(t_ORIS/t_raking)") +
      ggplot2::labs(
        title   = sprintf("GP posterior mean: %s = %.1f", slice_label, val),
        x       = x1_label,
        y       = x2_label
      ) +
      ggplot2::theme_minimal()

    # Save PDF with underscore-replaced value (dots → underscores for file safety)
    val_str <- sprintf("%.1f", val)
    val_str <- gsub("\\.", "_", val_str)
    filename <- sprintf("oris_vs_raking_3d_slice_%s.pdf", val_str)
    ggplot2::ggsave(file.path(out_dir, filename), p, width = 9, height = 6)

    plots[[length(plots) + 1]] <- p
    cat(sprintf("Saved: %s\n", file.path(out_dir, filename)))
  }

  invisible(plots)
}
