#' Gower dissimilarity for mixed-response clustering validation.
#'
#' Uses equal variable weights and cluster::daisy() with Gower's coefficient.
#' Continuous and count features are numeric and range-normalised by daisy;
#' binary and nominal features are treated as symmetric categorical variables;
#' ordinal variables retain ordered ranks. The external response label should
#' NEVER be included in `data`.
#'
#' @param data A complete predictor data frame, without external class labels.
#' @param types Named vector of response types matching columns in `data`.
#' @return A `dissimilarity` / `dist` object.
#' @export
mixeddeepmm_gower_distance <- function(data, types) {
  if (!requireNamespace("cluster", quietly=TRUE))
    stop("Install the 'cluster' package for Gower dissimilarity and silhouette.")
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (nrow(data) < 3L || ncol(data) < 1L || anyNA(data))
    stop("data must have at least 3 complete observations and 1 predictor.")
  if (is.null(names(types)) || anyDuplicated(names(types)) ||
      length(types) != ncol(data) || !setequal(names(types),names(data)))
    stop("types must name exactly the columns in data.")
  x <- data[, names(types), drop=FALSE]
  for (nm in names(types)) {
    kind <- tolower(types[[nm]])
    v <- x[[nm]]
    if (kind %in% c("binary","nominal")) {
      x[[nm]] <- factor(v)  # symmetric matches, including both binary levels
    } else if (kind=="ordinal") {
      if (is.ordered(v)) {
        x[[nm]] <- v
      } else if (is.numeric(v)) {
        x[[nm]] <- ordered(v, levels=sort(unique(v)))
      } else {
        stop("Ordinal predictor '", nm, "' must be an ordered factor or numeric with ordered codes.")
      }
    } else if (kind %in% c("continuous","count")) {
      if (!is.numeric(v) || any(!is.finite(v)))
        stop("Predictor '", nm, "' must be finite and numeric for Gower.")
      x[[nm]] <- as.numeric(v)
    } else {
      stop("Unknown response type for '", nm, "': ", kind)
    }
  }
  cluster::daisy(x, metric="gower")
}

#' Mean Gower silhouette width for a fixed clustering.
#'
#' This is a label-free *external-to-likelihood* clustering criterion.
#' Returns `NA_real_` for degenerate partitions (<2 or >=n occupied clusters).
#' Silhouette can be computed once per completed restart and used for model
#' selection without examining ARI or MR. Passing `distance` avoids repeated
#' calculation of the O(n^2) Gower matrix.
#'
#' @param clusters Cluster labels, one per observation.
#' @param distance Pre-computed Gower distance from mixeddeepmm_gower_distance.
#' @return Mean silhouette width (numeric scalar), or NA for invalid partitions.
#' @export
mixeddeepmm_gower_silhouette <- function(clusters, distance) {
  if (!requireNamespace("cluster", quietly=TRUE))
    stop("Install 'cluster' to evaluate Gower silhouette.")
  if (!inherits(distance, "dist")) stop("distance must be a dist/dissimilarity object.")
  n <- as.integer(attr(distance,"Size"))
  if (length(clusters)!=n || anyNA(clusters))
    stop("clusters must contain one non-missing label per observation.")
  k <- length(unique(clusters))
  if (k < 2L || k >= n) return(NA_real_)
  sil <- cluster::silhouette(as.integer(factor(clusters)), distance)
  val <- mean(sil[, "sil_width"])
  if (is.finite(val)) as.numeric(val) else NA_real_
}
