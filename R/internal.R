# Perform logit transformation
#
# @param x A numeric vector
#
# @return A vector of the same length as \code{x} containing transformed values
#
logit <- function(x) {
  return(log(x = x / (1 - x)))
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

# Calculate the area factor for a normal distribution
#
# @param sd A standard deviation
#
# @return The 1/sqrt(2 * pi * (sd ^ 2)) value
#
normfactor <- function(sd) {
  return(1 / sqrt(x = 2 * pi * (sd ^ 2)))
}

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
      normexp(x = logit(x = r), mean = mean, sd = sd)
  )
}

# Proof of concept stuff
# Not used for anything other than
# sanity checks and testing things

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
