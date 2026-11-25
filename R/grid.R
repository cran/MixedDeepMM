#' Compare architecture candidates with an explicit candidate list.
#'
#' Each candidate is a named list with family, heads, layers. Results use the
#' model's independent MC log likelihood check; BIC remains approximate for
#' non-Gaussian mixed fits. Smaller BIC is preferred, conditional on MC error.
#'
#' @param data Input data frame.
#' @param types Named response type vector.
#' @param candidates List of candidate specification lists.
#' @param control Shared fitting controls.
#' @return List with `table` and `fits` (failed fits are NULL).
#' @export
mixeddeepmm_grid <- function(data,types=NULL,candidates,control=list()) {
  if (!is.list(candidates) || !length(candidates)) stop("candidates must be a nonempty list.")
  fits <- vector("list",length(candidates)); tab <- vector("list",length(candidates))
  for (i in seq_along(candidates)) {
    s <- candidates[[i]]
    ans <- try(mixeddeepmm(data,types=types,family=s$family,heads=s$heads,
                           layers=s$layers,control=control),silent=TRUE)
    if (inherits(ans,"try-error")) {
      tab[[i]] <- data.frame(candidate=i,family=as.character(s$family),
                            heads=as.character(s$heads),logLik=NA_real_,
                            BIC=NA_real_,mc_se=NA_real_,status=as.character(ans)[1])
    } else {
      fits[[i]] <- ans
      tab[[i]] <- data.frame(candidate=i,family=ans$family,heads=ans$heads,
                            logLik=ans$loglik,BIC=ans$bic,
                            mc_se=ans$mc_se,status="fitted")
    }
  }
  table <- do.call(rbind,tab)
  table <- table[order(table$BIC,na.last=TRUE),,drop=FALSE]
  rownames(table) <- NULL
  list(table=table,fits=fits)
}
