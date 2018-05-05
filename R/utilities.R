# Calculate the binomial coefficent
#
# @param k Number of observations
# @param n Number of trials
#
# @return Binomial coefficients, any values of NaN will be changed to a 1
#
binom <- function(k, n) {
  coefs <- mapply(
    FUN = function(k.use, n.use) {
      return(factorial(x = n.use) / (factorial(x = k.use) * factorial(x = n.use - k.use)))
    },
    k.use = k,
    n.use = n
  )
  coefs[is.nan(x = coefs)] <- 1
  return(coefs)
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

blnratio <- function(r, x, size, mean, sd) {
  a <- r ^ (x - 1)
  b <- (1 - r) ^ (size - x - 1)
  c <- normexp(x = logit(x = r), mean = mean, sd = sd)
  return(a * b * c)
}

# betalogitnormal <- function(x, mean = 0, sd = 1, size) {
#   num.integrate <- 100
#   variance <- sd ^ 2
#   tau <- 2 * pi
#   if (variance < 1e-3) {
#     warning("Variance is less than 1e-3, giving binomial probability")
#     return(dbinom(x = x, size = size, prob = logistic(x = mean)))
#   }
#   xc <- size - x
#   z <- lgamma(x = size + 1) - lgamma(x = x + 1) - lgamma(x = xc + 1) - log(x = tau * variance) * 0.5
#   rd <- seq.int(from = 0, to = 1, by = num.integrate + 1)
#   rd <- rd + ((rd[2] - rd[1]) / 2)
#   rd <- rd[1:num.integrate]
# }

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
  probs <- unlist(x = probs[1, ])
  names(x = probs) <- NULL
  return(normfactor(sd = sd) * probs)
}
