# Tests for BLN distributions

# Load in the data
blnfile <- system.file('extdata', 'matlab_simluations.txt', package = 'bln')
matlab <- read.table(
  file = blnfile,
  header = TRUE,
  comment.char = '#',
  as.is = TRUE
)

# Test dbln
context(desc = "dbln")

test_that(
  desc = "dbln calculates approximated PDF correctly",
  code = {
    r.approx <- dbln(
      x = matlab$x,
      size = matlab$x + matlab$xc,
      mean = matlab$mu,
      sd = sqrt(x = matlab$v)
    )
    expect_length(object = r.approx, n = nrow(x = matlab))
    expect_equal(
      object = round(x = cor(x = r.approx, y = matlab$px)),
      expected = 1
    )
  }
)

test_that(
  desc = "dbln calculates accurate PDF correctly",
  code = {
    r.exact <- dbln(
      x = matlab$x,
      size = matlab$x + matlab$xc,
      mean = matlab$mu,
      sd = sqrt(x = matlab$v),
      drop = FALSE
    )
    expect_length(object = r.exact, n = nrow(x = matlab))
    expect_equal(
      object = round(x = cor(x = r.exact, y = matlab$px)),
      expected = 1
    )
  }
)
