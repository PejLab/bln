# Calculate the integral part of the BLN distribution formula
#
# @param r The binomial ratio for x
# @param x A binomial distributed variable
# @param size Number of trials (\code{n})
# @param mean Mean of the normal distribution
# @param sd Standard deviation of the normal distribution
#
# @return
#
blnratio <- function(r, x, size, mean, sd) {
  return(
    (r ^ (x - 1)) *
      ((1 - r) ^ (size - x - 1)) *
      exp(x = -((logit(x = r) - mean) ^ 2) / (2 * (sd ^ 2)))
  )
}

# Perform logistic transformation
#
# @param x A numeric vector
#
# @return A vector of the same length as \code{x} containing transformed values
#
logistic <- function(x) {
  return(1 / (1 + exp(x = -x)))
}

# Perform logit transformation
#
# @param x A numeric vector
#
# @return A vector of the same length as \code{x} containing transformed values
#
logit <- function(x) {
  return(log(x = x / (1 - x)))
}

# Determine the largest length out of many vectors
#
# @param ... One or more vectors
#
# @return The maximal length of all vectors passed
#
max.length <- function(...) {
  return(max(vapply(X = list(...), FUN = length, FUN.VALUE = numeric(length = 1L))))
}

### Proof of concept stuff
# Not used for anything other than
# sanity checks and testing things

# Calculate the expoent part of a normal distribution
#
# @param x The value where the normal distribution is as
# @param mean The mean of the normal distribution
# @param The value of the distribution
#
# @return The calculated exponent part of the normal distribution
#
normexp <- function(x, mean, sd) {
  return(exp(x = -((x - mean) ^ 2) / (2 * (sd ^ 2))))
}

# Calculate the area factor for a normal distribution
#
# @param sd A standard deviation
#
# @return The 1/sqrt(2 * pi * (sd ^ 2)) value
#
normfactor <- function(sd) {
  return(1 / sqrt(x = 2 * pi * (sd ^ 2)))
}

# Normal PDF function
#
# @param x A vector of quantiles
# @param mean A vector of means
# @param sd A vector of standard deviations
#
normpdf <- function(x, mean = 0, sd = 1) {
  length.use <- max(sapply(X = list(x, mean, sd), FUN = length))
  x <- rep_len(x = x, length.out = length.use)
  mean <- rep_len(x = mean, length.out = length.use)
  sd <- rep_len(x = sd, length.out = length.use)
  return(normfactor(sd = sd) * normexp(x = x, mean = mean, sd = sd))
}

# Normal CDF function
#
# @param q A vector of quantiles
# @param mean A vector of means
# @param sd A vector of standard deviations
#' @importFrom stats integrate
#
normcdf <- function(q, mean = 0, sd = 1) {
  length.use <- max(sapply(X = list(q, mean, sd), FUN = length))
  q <- rep_len(x = q, length.out = length.use)
  mean <- rep_len(x = mean, length.out = length.use)
  sd <- rep_len(x = sd, length.out = length.use)
  probs <- mapply(
    FUN = integrate,
    upper = q,
    mean = mean,
    sd = sd,
    MoreArgs = c(
      f = normexp,
      lower = -Inf
    )
  )
  probs <- unlist(x = probs[1, ], use.names = FALSE)
  return(normfactor(sd = sd) * probs)
}

erf <- function(x) {
  probs <- mapply(
    FUN = integrate,
    upper = x,
    MoreArgs = c(
      f = function(t) {
        return(exp(x = -(t ^ 2)))
      },
      lower = 0
    )
  )
  probs <- unlist(x = probs[1, ], use.names = FALSE)
  return((2 / sqrt(x = pi)) * probs)
}

normqf <- function(p, mean = 0, sd = 1) {
  expr <- expression(exp(x = -((x - mean) ^ 2) / (2 * (sd ^ 2))))
  invisible(x = NULL)
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
