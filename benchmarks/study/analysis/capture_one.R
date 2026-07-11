# STUDY-BRANCH-ONLY-DO-NOT-MERGE
# WU-12c (leafblower-2ouc.42) trajectory-capture runner (R side): one
# (solver, problem) harvest with LBW_TRAJECTORY_AT / LBW_TRAJECTORY_OUT set by
# the caller, writing the per-iteration RQ3 convergence CSV. Cold path
# (byte-identical weights, proven). Mirrors R/leafblower_adapter.R's mapping.
suppressMessages(library(leafblower))

.here <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                  error = function(e) ".")
.common <- file.path(.here, "..", "common")
source(file.path(.common, "problem_io.R"))       # load_problem_spec
source(file.path(.common, "instance_family.R"))  # install_gen_resolver
install_gen_resolver()                            # gen: instance specs (harmless for file:)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: capture_one.R <method> <spec_path>")
method <- args[[1]]; spec_path <- args[[2]]

p <- load_problem_spec(spec_path)
# LBW_TRAJECTORY_AT / LBW_TRAJECTORY_OUT are set in the environment by the caller.
invisible(tryCatch(
  harvest(
    data = p$data, target = p$targets, method = method,
    min_weight = p$bounds$min, max_weight = p$bounds$max,
    design_weights = p$design_weights,
    convergence = list(absolute = p$tol), attach_weights = FALSE
  ),
  error = function(e)  # infeasible/error cells have no convergence curve -- skip cleanly
    message(sprintf("skip %s/%s: %s", method, basename(spec_path), conditionMessage(e)))
))
