# BLN

An implementation of the binomial-logit-normal distribution in R

## Installation

You can install `bln` from GitHub using [devtools](https://www.rstudio.com/products/rpackages/devtools/). Note that to install from a private repository, you need to generate a [personal access token](https://github.com/settings/tokens). The token should have all of the `repo` scopes, no others are needed. You can set this as an environement variable (`GITHUB_PAT`) to save it and not to enter it every time.

``` r
# install.packages("devtools")
# leave out the auth_token argument if set as an environement variable
devtools::install_github(repo = "mojaveazure/bln", auth_token = 'abc')
```
