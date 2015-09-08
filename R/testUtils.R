#' Converts list(list(fields)) into list(fields(list))
#' 
#' @param dots list of lists
listcomb <- function(dots) {
  nms <- names(dots[[1]])
  ans <- as.list(nms)
  names(ans) <- nms
  for (nm in nms) {
    ans[[nm]] <- lapply(1:length(dots), function(i) dots[[i]][[nm]])
  }
  ans
}

#' Inverse square root of a matrix
#' 
#' @param m matrix
isqrtm <- function(m) {
  res <- eigen(m)
  d <- res$values
  if (min(d) < -1e-5) warning("Negative eigenvalues in isqrtm")
  d[d < 0] <- 0
  d[d > 0] <- 1/sqrt(d[d > 0])
  v <- res$vectors
  return (v %*% diag(d) %*% t(v))
}

#' Calls mclapply and combines the end result using listcomb
#' 
#' @import parallel
#' @examples
#' lclapply(1:10, function(i) list(a = i, b = i^2), mc.cores = 1)
lclapply <- function(x, f, mc.cores = 0) {
  if (mc.cores == 0) {
    return(listcomb(lapply(x, f)))
  } else {
    return(listcomb(mclapply(x, f, mc.cores = mc.cores)))    
  }
}

#' Either uses lapply or mclapply
#' 
#' @import parallel
#' @examples
#' lclapply(1:10, function(i) list(a = i, b = i^2), mc.cores = 1)
mclapply0 <- function(x, f, mc.cores = 0) {
  if (mc.cores == 0) {
    return(lapply(x, f))
  } else {
    return(mclapply(x, f, mc.cores = mc.cores))    
  }
}

