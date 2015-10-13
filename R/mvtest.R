#' Permutation test of F == G based on first two moments
#' To be replaced by Rcpp version
#' 
#' @param x (n1 x p) data matrix
#' @param y (n2 x p) data matrix
#' @param n.reps number of replications
#' @examples
#' n <- 20
#' p <- 3
#' x <- MASS::mvrnorm(n, rep(0, p), diag(rep(1, p)))
#' y <- MASS::mvrnorm(n, rep(0.1, p), diag(rep(1, p)))
#' mvtest(x, y, 1e3)
#' @export
mvtest <- function(x, y, n.reps = 1e3, seed = NA, fullOut = FALSE) {
  if (!is.na(seed)) {
    set.seed(seed)
  }
  n1 <- dim(x)[1]
  n2 <- dim(y)[1]
  p <- dim(x)[2]
  xy <- rbind(x, y)
  t_obs <- mvtest_stat(x, y)
  t_perm <- sapply(1:n.reps, function(i) {
    o <- sample(n1 + n2, n1)
    mvtest_stat(xy[o, ], xy[-o, ])
  })
  if (fullOut) {
    return(list(t_obs, t_perm))
  }
  sum(t_obs < t_perm)/n.reps
}

#' Permutation test of F == G based on first two moments
#' To be replaced by Rcpp version
#' 
#' @param x (n1 x p) data matrix
#' @param y (n2 x p) data matrix
#' @export
mvtest_stat <- function(x, y) {
  mu1 <- colMeans(x)
  mu2 <- colMeans(y)
  Sigma1 <- cov(x)
  Sigma2 <- cov(y)
  diff <- mu1 - mu2
  sum(diag(solve(Sigma1, Sigma2))) + sum(diag(solve(Sigma2, Sigma1))) +
    (t(diff) %*% (solve(Sigma1, diff) + solve(Sigma2, diff)))[1]
}

