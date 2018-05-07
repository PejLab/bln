#' @include internal.R
#'
NULL

#' The Binomial-Logit-Normal Distribution
#'
#' @param x,q A vector of quantiles
#' @param p A vector of probabilities
#' @param n Number of observations, if \code{length(x = n) > 1},
#' the length is taken to be the number required
#' @param size Number of trials
#' @param mean A vector of means
#' @param sd A vector of standard deviations
#'
#' @rdname bln
#' @name Binomial-logit-normal
#' @aliases bln
#'
NULL

# Density (PDF)
#' @rdname bln
#' @aliases dbln
#' @return \code{dbln} gives the density
#' @importFrom stats dbinom
#' @export
#'
dbln <- function(x, size, mean = 0, sd = 1) {
  # Equalize lengths
  length.use <- max(sapply(X = list(x, size, mean, sd), FUN = length))
  x <- rep_len(x = x, length.out = length.use)
  size <- rep_len(x = size, length.out = length.use)
  mean <- rep_len(x = mean, length.out = length.use)
  sd <- rep_len(x = sd, length.out = length.use)
  # Some constants
  num.integrate <- 100
  variance <- sd ^ 2
  tau <- 2 * pi
  if (any(variance < 1e-3)) {
    warning("Variance is less than 1e-3, giving binomial probability")
    return(dbinom(x = x, size = size, prob = logistic(x = mean)))
  }
  xc <- size - x
  z <- lgamma(x = size + 1) - lgamma(x = x + 1) - lgamma(x = xc + 1)
  z <- z - (log(x = tau * variance) * 0.5)
  rd <- seq.int(from = 0, to = 1, length.out = num.integrate + 1)
  rd <- rd + ((rd[2] - rd[1]) / 2)
  rd <- rd[1:num.integrate]
  return(mapply(
    FUN = function(r, x.use, xc.use, mean.use, variance.use, z.use) {
      f <- exp(
        x = (log(x = r) * (x.use - 1)) +
          (log(x = 1 - r) * (xc.use - 1)) -
          (((logit(x = r) - mean.use) ^ 2) / (2 * variance.use)) +
          z.use
      )
      return(exp(x = -log(x = num.integrate) + log(x = sum(f))))
    },
    x.use = x,
    xc.use = xc,
    mean.use = mean,
    variance.use = variance,
    z.use = z,
    MoreArgs = list(r = rd)
  ))
}

# Probability function (CDF)
#' @rdname bln
#' @aliases pbln
#' @return \code{pbln} gives the distribution function
#' @importFrom stats integrate
#' @export
#'
pbln <- function(q, size, mean = 0, sd = 1) {
  # Equalize lengths
  length.use <- max(sapply(X = list(q, size, mean, sd), FUN = length))
  q <- rep_len(x = q, length.out = length.use)
  size <- rep_len(x = size, length.out = length.use)
  mean <- rep_len(x = mean, length.out = length.use)
  sd <- rep_len(x = sd, length.out = length.use)
  # Calculate the integral
  probs <- mapply(
    FUN = integrate,
    x = q,
    size = size,
    mean = mean,
    sd = sd,
    MoreArgs = list(
      f = blnratio,
      lower = 0,
      upper = 1
    )
  )
  probs <- unlist(x = probs[1, ])
  names(x = probs) <- NULL
  #
  return(choose(n = size, k = q) * normfactor(sd = sd) * probs)
}

# Quantile function
#' @rdname bln
#' @aliases qbln
#' @return \code{qbln} gives the quantile function
#' @export
#'
qbln <- function(...) {
  E <- expression(
    r ^ (x - 1) *
      (1 - r) ^ (size - x - 1) *
      exp(x = -((x - mean) ^ 2) / 2 * (sd ^ 2))
  )
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

# # @param x Binomally distributed variable
# # @param xc \code{x + xc} is the total number of trials in the binomial
# # @param mu mean
# # @param v variance
# # @param dbc DropBinomialCoeff, I'm not sure if this is used
# #
# bln <- function(x, xc, mu, v, dbc) {
#   nn <- 100
#   if (v < 1e-3) {
#     prob <- dbinom(x = x, size = x + xc, prob = logistic(x = mu))
#     if (is.infinite(x = prob)) {
#       stop("Numerical failure!")
#     }
#     return(prob)
#   }
#   # Z <- gammaln(x = x + xc + 1) - gammaln(x = x + 1) - gammaln(x = xc + 1) -
#   #   log(x = 2 * pi * v) * 0.5
#   Z <- lgamma(x = x + xc + 1) - lgamma(x = x + 1) - lgamma(x = xc + 1) -
#     log(x = 2 * pi * v) * 0.5
#   rd <- seq.int(from = 0, to = 1, length.out = nn + 1)
#   rd2 <- rd + ((rd[2] - rd[1]) / 2)
#   rd2 <- rd[1:nn]
#   invisible(x = NULL)
# }

# # @param x Binomiallly distributed variable
# # @param n ...
# # @param mean ...
# # @param sd ...
# #
# ibln <- function(x, n, mean = 0, sd = 1) {
#   # Some simple things
#   variance <- sd ^ 2
#   tau <- 2 * pi
#   j <- choose(n = n, k = x) * (1 / sqrt(x = tau * variance))
#   # Prepare integration
#   length.use <- max(length(x = x), length(x = n))
#   x <- rep_len(x = x, length.out = length.use)
#   n <- rep_len(x = n, length.out = length.use)
#   f <- function(r, x.use, n.use) {
#     a <- r ^ (x.use - 1)
#     b <- (1 - r) ^ (n.use - x.use - 1)
#     c <- exp(x = -(((logit(x = r) - mean) ^ 2) / (2 * variance)))
#     return(a * b * c)
#   }
#   # Do the integration
#   to.return <- mapply(
#     FUN = integrate,
#     x.use = x,
#     n.use = n,
#     MoreArgs = list(
#       f = f,
#       lower = 0,
#       upper = 1
#     )
#   )
#   to.return <- unlist(x = to.return[1, ])
#   return(j * to.return)
# }
