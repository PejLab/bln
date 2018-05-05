#' @include utilities.R
#'
NULL

# binom variables
# x, q <- quantiles (d, p)
# p <- probabilities (q)
# n <- number of observations (r)
# size <- number of trials (all)
# prob <- probabilitiy of success (all)
# lower.tail <- P[X <= x], else P[X > x] (p, q)

# norm variables
# x, q <- quantiles (d, p)
# p <- probabilities (q)
# n <- number of observations (r)
# mean <- means (all)
# sd <- stdevs (all)
# lower.tail <- P[X <= x], else P[X > x] (p, q)

#' The Binomial-Logit-Normal Distribution
#'
#' @param x,q A vector of quantiles
#' @param n ...
#' @param size Number of trials
#' @param prob ...
#' @param mean A vector of means
#' @param sd A vector of standard deviations
#'
#' @rdname bln
#' @name Binomial-logit-normal
#'
NULL

# Density (PDF)
#' @rdname bln
#' @aliases dbln
#' @return \code{dbln} gives the density
#' @export
#'
dbln <- function(x, size, mean = 0, sd = 1) {
  invisible(x = NULL)
}

# Probability function (CDF)
#' @rdname bln
#' @aliases pbln
#' @return \code{pbln} gives the distribution function
#' @importFrom stats integrate
#' @export
#'
pbln <- function(q, size, mean = 0, sd = 1) {
  length.use <- max(sapply(X = list(q, size, mean), FUN = length))
  q <- rep_len(x = q, length.out = length.use)
  size <- rep_len(x = size, length.out = length.use)
  mean <- rep_len(x = size, length.out = length.use)
  probs <- mapply(
    FUN = integrate,
    upper = q,
    size = size,
    mean = mean,
    MoreArgs = c(
      f = blnratio,
      lower = -Inf,
      sd = sd
    )
  )
  probs <- unlist(x = probs[1, ])
  names(x = probs) <- NULL
  return(binom(k = q, n = size) * normfactor(sd = sd) * probs)
  # binom(k = q, n = size)
  invisible(x = NULL)
}

# Quantile function
#' @rdname bln
#' @aliases qbln
#' @return \code{qbln} gives the quantile function
#' @export
#'
qbln <- function(...) {
  invisible(x = NULL)
}

# Random number generation
#' @rdname bln
#' @aliases rbln
#' @return \code{rbln} generates random deviates
#' @importFrom stats rnorm rbinom
#' @export
#'
rbln <- function(n, size, mean = 0, sd = 1) {
  x <- rnorm(n = n, mean = mean, sd = sd)
  y <- logistic(x = x)
  z <- rbinom(n = n, size = size, prob = y)
  return(z)
}

# @param x Binomally distributed variable
# @param xc \code{x + xc} is the total number of trials in the binomial
# @param mu mean
# @param v variance
# @param dbc DropBinomialCoeff, I'm not sure if this is used
#
bln <- function(x, xc, mu, v, dbc) {
  nn <- 100
  if (v < 1e-3) {
    prob <- dbinom(x = x, size = x + xc, prob = logistic(x = mu))
    if (is.infinite(x = prob)) {
      stop("Numerical failure!")
    }
    return(prob)
  }
  # Z <- gammaln(x = x + xc + 1) - gammaln(x = x + 1) - gammaln(x = xc + 1) -
  #   log(x = 2 * pi * v) * 0.5
  Z <- lgamma(x = x + xc + 1) - lgamma(x = x + 1) - lgamma(x = xc + 1) -
    log(x = 2 * pi * v) * 0.5
  rd <- seq.int(from = 0, to = 1, length.out = nn + 1)
  rd2 <- rd + ((rd[2] - rd[1]) / 2)
  rd2 <- rd[1:nn]
  invisible(x = NULL)
}

# @param x Binomiallly distributed variable
# @param n ...
# @param mean ...
# @param sd ...
#
ibln <- function(x, n, mean = 0, sd = 1) {
  # Some simple things
  variance <- sd ^ 2
  tau <- 2 * pi
  j <- binom(k = x, n = n) * (1 / sqrt(x = tau * variance))
  # Prepare integration
  length.use <- max(length(x = x), length(x = n))
  x <- rep_len(x = x, length.out = length.use)
  n <- rep_len(x = n, length.out = length.use)
  f <- function(r, x.use, n.use) {
    a <- r ^ (x.use - 1)
    b <- (1 - r) ^ (n.use - x.use - 1)
    c <- exp(x = -(((logit(x = r) - mean) ^ 2) / (2 * variance)))
    return(a * b * c)
  }
  # Do the integration
  to.return <- mapply(
    FUN = integrate,
    x.use = x,
    n.use = n,
    MoreArgs = c(
      f = f,
      lower = 0,
      upper = 1
    )
  )
  to.return <- unlist(x = to.return[1, ])
  return(j * to.return)
}
