# Calculate the integral part of the BLN distribution formula
#
# @param r The binomial ratio for x
# @param x A binomial distributed variable
# @param size Number of trials (\code{n})
# @param mean Mean of the normal distribution
# @param sd Standard deviation of the normal distribution
#
# @return ...
#
blnratio <- function(r, x, size, mean, sd) {
  return(
    (r ^ (x - 1)) *
      ((1 - r) ^ (size - x - 1)) *
      exp(x = -((logit(x = r) - mean) ^ 2) / (2 * (sd ^ 2)))
  )
}

# fx, used in calculation of PDF
#
# @param ratio The binomial ratio for x
# @param x A binomially distributed variable
# @param xc ...
# @param mean Mean of the normal distribution
# @param sd Standard deviation of the normal distribution
# @param z ...
#
# @return ...
#
fxpdf <- function(ratio, x, xc, mean, variance, z) {
  return(exp(
    x = (log(x = ratio) * (x - 1)) +
      (log(x = 1 - ratio) * (xc - 1)) -
      (((logit(x = ratio) - mean) ^ 2) / (2 * variance)) +
      z
  ))
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
