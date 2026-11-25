#' Misclassification rate (MR) after optimal relabelling of clusters.
#'
#' Computes the smallest fraction of incorrectly classified observations after
#' one-to-one matching the predicted clusters to the reference classes. Matching
#' is solved by the Hungarian algorithm; neither class nor cluster IDs need to
#' have the same numerical values or order. For unequal numbers of classes and
#' clusters, a zero-padded assignment matrix is used: unmatched groups count as
#' misclassified. This is an external validation metric: never use reference
#' labels to fit models or to choose the representative seed.
#'
#' @param truth Known class labels, with no missing values.
#' @param predicted Estimated cluster labels, with no missing values.
#' @return A numeric scalar between 0 and 1 (0 is perfect agreement).
#' @examples
#' mixeddeepmm_mr(c("A", "A", "B", "B"), c(9, 9, 4, 4))
#' mixeddeepmm_mr(c("A", "A", "B", "B"), c(9, 9, 4, 9))
#' @export
mixeddeepmm_mr <- function(truth, predicted) {
  if (!is.atomic(truth) || !is.atomic(predicted) ||
      !is.null(dim(truth)) || !is.null(dim(predicted)) ||
      !length(truth) || length(truth) != length(predicted) ||
      anyNA(truth) || anyNA(predicted)) {
    stop("truth and predicted must be nonempty, equal-length, non-missing vectors.")
  }

  tab <- unclass(table(as.character(truth), as.character(predicted)))
  nr <- nrow(tab)
  nc <- ncol(tab)
  size <- max(nr, nc)
  weight <- matrix(0, nrow = size, ncol = size)
  weight[seq_len(nr), seq_len(nc)] <- tab
  # Minimise cost == maximise correctly matched observations.
  cost <- max(weight) - weight

  # Hungarian shortest-augmenting-path implementation using one-based indices.
  # Column 1 is the sentinel (the mathematical 'column 0').
  u <- numeric(size + 1L)
  v <- numeric(size + 1L)
  p <- integer(size + 1L)
  way <- integer(size + 1L)
  for (i in seq_len(size)) {
    p[1L] <- i
    j0 <- 1L
    minv <- rep(Inf, size + 1L)
    used <- rep(FALSE, size + 1L)
    repeat {
      used[j0] <- TRUE
      i0 <- p[j0]
      delta <- Inf
      j1 <- 0L
      for (j in seq.int(2L, size + 1L)) {
        if (!used[j]) {
          cur <- cost[i0, j - 1L] - u[i0 + 1L] - v[j]
          if (cur < minv[j]) {
            minv[j] <- cur
            way[j] <- j0
          }
          if (minv[j] < delta) {
            delta <- minv[j]
            j1 <- j
          }
        }
      }
      for (j in seq_len(size + 1L)) {
        if (used[j]) {
          u[p[j] + 1L] <- u[p[j] + 1L] + delta
          v[j] <- v[j] - delta
        } else {
          minv[j] <- minv[j] - delta
        }
      }
      j0 <- j1
      if (p[j0] == 0L) break
    }
    repeat {
      j1 <- way[j0]
      p[j0] <- p[j1]
      j0 <- j1
      if (j0 == 1L) break
    }
  }
  assigned_col <- integer(size)
  for (j in seq.int(2L, size + 1L)) {
    if (p[j] > 0L) assigned_col[p[j]] <- j - 1L
  }
  matched <- sum(weight[cbind(seq_len(size), assigned_col)])
  max(0, min(1, (length(truth) - matched) / length(truth)))
}
