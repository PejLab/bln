#' @include internal.R
#' @useDynLib bln
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
#' @param drop Drop the binomial coefficient? If \code{FALSE}, calculates the
#' normalized value. If \code{NULL}, calculates the approximated value.
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
#' @importFrom pracma eps
#' @importFrom stats dbinom
#' @export
#'
dbln <- function(x, size, mean = 0, sd = 1, drop = NULL) {
  # Equalize lengths
  length.use <- max.length(x, size, mean, sd)
  x <- rep_len(x = x, length.out = length.use)
  size <- rep_len(x = size, length.out = length.use)
  mean <- rep_len(x = mean, length.out = length.use)
  sd <- rep_len(x = sd, length.out = length.use)
  results <- vector(mode = 'numeric', length = length.use)
  # Some constants
  num.integrate <- 250
  variance <- sd ^ 2
  min.var <- eps()
  xc <- size - x
  z <- if (isTRUE(x = drop)) {
    rep_len(x = 0, length.out = length.use)
  } else {
    lgamma(x = size + 1) - lgamma(x = x + 1) - lgamma(x = xc + 1)
  }
  dx <- seq.int(
    from = 0,
    to = 1,
    length.out = ifelse(test = is.null(x = drop), yes = num.integrate, no = 1000)
  )
  too.small <- which(x = variance < min.var)
  # Do shit
  if (length(x = too.small) > 0) {
    warning("One or more variance values is less than ", min.var, ", giving binomial probability")
  }
  for (i in 1:length.use) {
    results[i] <- if (i %in% too.small) {
      dbinom(x = x[i], size = size[i], prob = logistic(x = mean[i]))
    } else {
      tpx <- pln(q = dx, mean = mean[i], sd = sd[i])
      lower <- max(dx[tpx < eps()])
      upper <- min(dx[tpx > (1 - eps())])
      if (is.null(x = drop)) {
        rd <- seq.int(from = lower, to = upper, length.out = num.integrate + 1)
        rd <- rd + ((rd[2] - rd[1]) / 2)
        rd <- rd[1:(length(x = rd) - 1)]
        f <- fxpdf(
          ratio = rd,
          x = x[i],
          xc = xc[i],
          mean = mean[i],
          variance = variance[i],
          z = z[i] - (log(x = 2 * pi * variance[i]) * 0.5)
        )
        exp(x = -log(x = num.integrate / (upper - lower)) + log(x = sum(f)))
      } else {
        probs <- integrate(
          f = fxpdf,
          lower = lower,
          upper = upper,
          x = x[i],
          xc = xc[i],
          mean = mean[i],
          variance = variance[i],
          z = z[i],
          rel.tol = 1e-4,
          abs.tol = 0
        )
        1 / sqrt(x = 2 * pi * variance[i]) * probs$value
      }
    }
  }
  return(results)
}

# @rdname bln
# @aliases dblnx
# @references \code{dblnx} gives the exact density
# @export
#
# dblnx <- function(x, size, mean = 0, sd = 1, drop = FALSE) {
#   # Equalize lengths
#   length.use <- max.length(x, size, mean, sd)
#   x <- rep_len(x = x, length.out = length.use)
#   size <- rep_len(x = size, length.out = length.use)
#   mean <- rep_len(x = mean, length.out = length.use)
#   sd <- rep_len(x = sd, length.out = length.use)
#   xc <- size - x
#   # Some constants
#   variance <- sd ^ 2
#   min.var <- eps()
#   # Do shit
#   too.small <- which(x = variance < min.var)
#   if (length(x = too.small) > 0) {
#     warning("One or more variance values is less than ", min.var, ", giving binomial probability")
#   }
#   z <- if (drop) {
#     0
#   } else {
#     lgamma(x = size + 1) - lgamma(x + 1) - lgamma(x = xc + 1)
#   }
#   probs <- mapply(
#     FUN = integrate,
#     x = x,
#     xc = xc,
#     mean = mean,
#     variance = variance,
#     z = z,
#     MoreArgs = list(
#       f = fxpdf,
#       lower = 0,
#       upper = 1,
#       rel.tol = 1e-4,
#       abs.tol = 0
#     )
#   )
#   probs <- unlist(x = probs[1, ], use.names = FALSE)
#   return(1 / sqrt(x = 2 * pi * variance) * probs)
# }

#' @rdname bln
#' @aliases dblnpp
#' @return \code{dblnpp} gives the density using C++
#' @export
#'
dblnpp <- function(x, size, mean = 0, sd = 1) {
  # Equalize lengths
  length.use <- max.length(x, size, mean, sd)
  x <- rep_len(x = x, length.out = length.use)
  size <- rep_len(x = size, length.out = length.use)
  mean <- rep_len(x = mean, length.out = length.use)
  sd <- rep_len(x = sd, length.out = length.use)
  min.var <- eps()
  if (any(sd ^ 2 < min.var)) {
    warning("One or more variance values is less than ", min.var, ", giving binomial probability")
  }
  return(blnpdf(x = x, size = size, mean = mean, sd = sd))
}

# @rdname bln
# @aliases dblnxpp
# @return \code{dblnxpp} gives the exact density using C++
# @export
#
# dblnxpp <- function(x, size, mean = 0, sd = 1) {
#   invisible(x = NULL)
# }

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
