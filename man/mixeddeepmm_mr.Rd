\name{mixeddeepmm_mr}
\alias{mixeddeepmm_mr}
\title{Compute a label-invariant misclassification rate}
\description{Compute the minimum misclassification rate (MR) after optimally matching
  predicted cluster labels to known class labels, using one-to-one assignment.}
\usage{mixeddeepmm_mr(truth, predicted)}
\arguments{
  \item{truth}{Nonempty vector of reference class labels, without missing values.}
  \item{predicted}{Nonempty vector of estimated cluster labels, without missing values.}
}
\details{Cluster labels are arbitrary. The Hungarian algorithm finds the
  one-to-one correspondence maximising the diagonal contingency-table count,
  with zero padding for unequal numbers of reference classes and predicted
  clusters. Unmatched groups therefore contribute errors. Reference labels
  must be kept out of unsupervised fitting and model selection.}
\value{A numeric scalar between 0 and 1: 0 is perfect agreement and smaller is
  better.}
\examples{
mixeddeepmm_mr(c("A", "A", "B", "B"), c(9, 9, 4, 4))
mixeddeepmm_mr(c("A", "A", "B", "B"), c(9, 9, 4, 9))
}
