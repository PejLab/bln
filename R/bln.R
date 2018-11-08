#' @include internal.R
#'
NULL

#' The Binomial-Logit-Normal Distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the binomial-logit-normal distribution with mean equal to \code{mean} and
#' standard deviation equal to \code{sd}, as well as distribution parameter \code{size}.
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
  # length.use <- max(sapply(X = list(x, size, mean, sd), FUN = length))
  length.use <- max.length(x, size, mean, sd)
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
  # length.use <- max(sapply(X = list(q, size, mean, sd), FUN = length))
  length.use <- max.length(q, size, mean, sd)
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
  probs <- unlist(x = probs[1, ], use.names = FALSE)
  return(choose(n = size, k = q) * normfactor(sd = sd) * probs)
}

# Quantile function
#' @rdname bln
#' @aliases qbln
#' @return \code{qbln} gives the quantile function
#' @export
#'
qbln <- function(p, size, mean = 0, sd = 1) {
  .NotYetImplemented()
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
