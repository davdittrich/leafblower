# problem_io.R — leafblower benchmark study problem-spec loader.
#
# Resolves spec/*.json problem specs (data_ref origins file:/pkg:/gen:/inline)
# into a standardized problem object. Design of record: docs/benchmark/DESIGN.md
# §3; schema: benchmarks/study/spec/schema.json; ticket leafblower-2ouc.3 (WU-2).
#
# All file paths in specs (data_ref 'file:', targets_ref) are relative to the
# repo root, matching the existing benchmarks/*.R convention
# (arrow::read_parquet("benchmarks/...")) and CLAUDE.md's Build & Test
# invocation (run from repo root).
#
# Standardized problem object (list), returned by load_problem_spec():
#   id                  character(1)
#   data                data.frame, nrow = n; margin columns coerced to factor
#   design_weights      numeric(n), d_i (defaults to rep(1, n) when omitted)
#   margins             character(K), column names in `data`
#   targets             named list of K named numeric vectors T_kj, each
#                       margin's category values normalized to sum to 1
#   bounds              list(min = <double>, max = <double>); NULL sides
#                       resolved to min=0 / max=Inf (harvest.R:454 allows +Inf)
#   tol                 numeric(1)
#   objective_families  character
#   K                   integer(1), == length(margins)

# `[[<name>, exact = TRUE]]` -- NOT `$` -- is mandatory for every field read
# on a parsed spec list. R's `$` does *partial* prefix matching on list names
# when there is no exact match; "data"/"targets" are literal prefixes of
# "data_ref"/"targets_ref" in this exact schema, so `spec$data` on a
# data_ref='file:...' spec silently returns the data_ref STRING instead of
# NULL. `exact = TRUE` disables partial matching entirely (returns NULL).
.pio_get <- function(spec, name) spec[[name, exact = TRUE]]

load_problem_spec <- function(spec_path) {
  if (!file.exists(spec_path)) {
    stop("problem_io: spec file not found: ", spec_path, call. = FALSE)
  }
  spec <- jsonlite::fromJSON(spec_path, simplifyVector = FALSE)
  .pio_validate_spec(spec)

  data <- .pio_resolve_data_ref(.pio_get(spec, "data_ref"), .pio_get(spec, "data"))
  n <- nrow(data)

  design_weights <- .pio_resolve_design_weights(.pio_get(spec, "design_weights"), data, n)

  targets_raw <- .pio_resolve_targets(.pio_get(spec, "targets_ref"), .pio_get(spec, "targets"))
  targets <- lapply(targets_raw, function(t) {
    v <- unlist(t)
    v / sum(v)
  })

  margins <- unlist(.pio_get(spec, "margins"))
  missing_cols <- setdiff(margins, names(data))
  if (length(missing_cols)) {
    stop("problem_io: margin column(s) not found in resolved data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!setequal(names(targets), margins)) {
    stop("problem_io: targets_ref/targets keys must match `margins` exactly. ",
         "margins=[", paste(margins, collapse = ","), "] targets=[",
         paste(names(targets), collapse = ","), "]", call. = FALSE)
  }
  for (nm in margins) data[[nm]] <- factor(data[[nm]])

  K <- as.integer(.pio_get(spec, "K"))
  if (K != length(margins)) {
    stop("problem_io: spec$K (", K, ") != length(margins) (", length(margins),
         ")", call. = FALSE)
  }

  spec_bounds <- .pio_get(spec, "bounds")
  bounds <- list(
    min = if (is.null(.pio_get(spec_bounds, "min"))) 0 else as.numeric(.pio_get(spec_bounds, "min")),
    max = if (is.null(.pio_get(spec_bounds, "max"))) Inf else as.numeric(.pio_get(spec_bounds, "max"))
  )

  list(
    id = .pio_get(spec, "id"),
    data = data,
    design_weights = design_weights,
    margins = margins,
    targets = targets,
    bounds = bounds,
    tol = as.numeric(.pio_get(spec, "tol")),
    objective_families = unlist(.pio_get(spec, "objective_families")),
    K = K
  )
}

.pio_validate_spec <- function(spec) {
  required <- c("id", "data_ref", "margins", "bounds", "tol",
                "objective_families", "K")
  missing <- setdiff(required, names(spec))
  if (length(missing)) {
    stop("problem_io: spec missing required field(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  data_ref <- .pio_get(spec, "data_ref")
  valid_origin <- grepl("^file:", data_ref) || grepl("^pkg:", data_ref) ||
    grepl("^gen:", data_ref) || identical(data_ref, "inline")
  if (!valid_origin) {
    stop("problem_io: invalid data_ref origin: '", data_ref,
         "' (expected file:/pkg:/gen:/inline)", call. = FALSE)
  }
  if (identical(data_ref, "inline") && is.null(.pio_get(spec, "data"))) {
    stop("problem_io: data_ref='inline' requires a `data` field in the spec",
         call. = FALSE)
  }
  has_ref <- !is.null(.pio_get(spec, "targets_ref"))
  has_inline <- !is.null(.pio_get(spec, "targets"))
  if (has_ref == has_inline) {
    stop("problem_io: spec must set exactly one of targets_ref / targets",
         call. = FALSE)
  }
}

.pio_resolve_data_ref <- function(data_ref, inline_data) {
  if (grepl("^file:", data_ref)) {
    path <- sub("^file:", "", data_ref)
    if (!file.exists(path)) {
      stop("problem_io: file not found: ", path, call. = FALSE)
    }
    ext <- tolower(tools::file_ext(path))
    df <- if (ext == "parquet") {
      arrow::read_parquet(path)
    } else if (ext == "csv") {
      utils::read.csv(path, stringsAsFactors = FALSE)
    } else {
      stop("problem_io: unsupported file: extension '", ext, "' (path=",
           path, ")", call. = FALSE)
    }
    df <- as.data.frame(df)
    if ("uuid" %in% names(df)) df$uuid <- NULL
    df
  } else if (grepl("^pkg:", data_ref)) {
    pkg_spec <- sub("^pkg:", "", data_ref)
    parts <- strsplit(pkg_spec, "::", fixed = TRUE)[[1]]
    if (length(parts) != 2) {
      stop("problem_io: pkg: data_ref must be 'pkg:<package>::<dataset>', got: ",
           data_ref, call. = FALSE)
    }
    pkg <- parts[1]
    dataset <- parts[2]
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("problem_io: package '", pkg, "' not installed (required for ",
           "data_ref='", data_ref, "')", call. = FALSE)
    }
    # Some packages bundle several datasets in one lazy-load group under a
    # different group name than the individual object (e.g. survey::apistrat
    # ships inside the "api" group: data(package="survey")$results lists
    # Item="apistrat (api)"). Resolve the group name from the package's own
    # data index rather than assuming dataset == group.
    idx <- utils::data(package = pkg)$results
    item_col <- idx[, "Item"]
    hit <- which(item_col == dataset |
                   sub(" \\([^)]*\\)$", "", item_col) == dataset)
    if (length(hit) == 0) {
      stop("problem_io: dataset '", dataset, "' not found in package '", pkg,
           "' (data_ref=", data_ref, ")", call. = FALSE)
    }
    load_name <- item_col[hit[1]]
    group <- sub("^.*\\(([^)]*)\\)$", "\\1", load_name)
    if (identical(group, load_name)) group <- dataset  # no "(group)" suffix
    env <- new.env()
    utils::data(list = group, package = pkg, envir = env)
    if (!exists(dataset, envir = env, inherits = FALSE)) {
      stop("problem_io: dataset '", dataset, "' not found after loading ",
           "group '", group, "' from package '", pkg, "' (data_ref=",
           data_ref, ")", call. = FALSE)
    }
    as.data.frame(get(dataset, envir = env))
  } else if (grepl("^gen:", data_ref)) {
    recipe_id <- sub("^gen:", "", data_ref)
    stop("problem_io: gen: data_ref origin (parametric synthetic instance ",
         "generator) is WU-3 scope, not yet implemented (recipe_id='",
         recipe_id, "')", call. = FALSE)
  } else if (identical(data_ref, "inline")) {
    if (is.null(inline_data)) {
      stop("problem_io: data_ref='inline' requires a `data` field in the spec",
           call. = FALSE)
    }
    rows <- lapply(inline_data, function(row) {
      as.data.frame(row, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  } else {
    stop("problem_io: unrecognised data_ref origin: '", data_ref,
         "' (expected file:/pkg:/gen:/inline)", call. = FALSE)
  }
}

.pio_resolve_design_weights <- function(design_weights, data, n) {
  if (is.null(design_weights)) {
    return(rep(1, n))
  }
  if (is.character(design_weights) && length(design_weights) == 1L) {
    if (!grepl("^column:", design_weights)) {
      stop("problem_io: unrecognised design_weights string: '",
           design_weights, "' (expected 'column:<name>')", call. = FALSE)
    }
    col <- sub("^column:", "", design_weights)
    if (!col %in% names(data)) {
      stop("problem_io: design_weights column '", col,
           "' not found in data", call. = FALSE)
    }
    return(as.numeric(data[[col]]))
  }
  arr <- as.numeric(unlist(design_weights))
  if (length(arr) != n) {
    stop("problem_io: inline design_weights length ", length(arr),
         " != n ", n, call. = FALSE)
  }
  arr
}

.pio_resolve_targets <- function(targets_ref, inline_targets) {
  has_ref <- !is.null(targets_ref)
  has_inline <- !is.null(inline_targets)
  if (has_ref == has_inline) {
    stop("problem_io: spec must set exactly one of targets_ref / targets",
         call. = FALSE)
  }
  if (has_ref) {
    if (!file.exists(targets_ref)) {
      stop("problem_io: targets_ref file not found: ", targets_ref,
           call. = FALSE)
    }
    jsonlite::fromJSON(targets_ref, simplifyVector = FALSE)
  } else {
    inline_targets
  }
}
