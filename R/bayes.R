#' Tests accuracy of posterior distribution
#' 
#' Given black-box functions \code{prior_draw}, \code{data_draw}, and \code{posterior_draw},
#' checks that the joint distribution of the prior and data match the joint distribution
#' of the posterior and data.
#' Draws a number of (parameter, data, posterior) triplets 
#' and then partitions data into clusters.
#' Checks that the first two moments of the parameter are consistent with the 
#' first two moments of the posterior within each cluster.
#' @param prior_draw function() which produces a parameter vector \code{theta}
#' @param data_draw function(theta) which produces data vector \code{x}
#' @param posterior_draw function(x) which produces \code{theta_hat}
#' from the posterior conditional on x and the supplied prior
#' @param n Number of (theta, x) samples to draw
#' @param k Number of clusters to split the data
#' @param mc.cores Number of parallel cores to use
#' @return List with global p-value for equality of posterior
#' and prior in each of k clusters, as well as generated data and clusters
#' @import parallel
#' @export
#' @examples
#' p <- 3
#' prior_draw <- function() rnorm(p)
#' data_draw <- function(theta) theta + rnorm(p)
#' posterior_draw <- function(x) (x + rnorm(p))/2
#' res <- posterior_check(prior_draw, data_draw, posterior_draw)
posterior_check <- function(prior_draw, data_draw, posterior_draw,
                            n = 1000, k = 10, mc.cores = 1) {
  generate_triple <- function(i) {
    theta <- prior_draw()
    x <- data_draw(theta)
    theta_hat <- posterior_draw(x)
    list(theta = theta, x = x, theta_hat = theta_hat)
  }
  res <- lclapply(1:n, generate_triple, mc.cores)
  priors <- do.call(rbind, res$theta)
  datas <- do.call(rbind, res$x)
  posteriors <- do.call(rbind, res$theta_hat)
  Sigma <- cov(datas)
  wht <- isqrtm(Sigma)
  wht_datas <- datas %*% wht
  cl <- kmeans(wht_datas, k)$cluster
  prior_moments <- mclapply0(1:k, 
                             function(i) {
                               c(colMeans(priors[cl ==i, ]),
                                 as.numeric(cov(priors[cl ==i, ])))
                             }, mc.cores)
  post_moments <- mclapply0(1:k, 
                             function(i) {
                               c(colMeans(posteriors[cl ==i, ]),
                                 as.numeric(cov(posteriors[cl ==i, ])))
                             }, mc.cores)
  list(prior_moments = prior_moments, post_moments = post_moments) 
}