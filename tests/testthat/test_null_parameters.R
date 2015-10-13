# unit testing for the CNV classification functions

library(cifer)
library(testthat)

context("CNV classification")
set.seed(2)

test_that("get_maxima output is correct", {
    # Note that these simple, explicitly normally distributed mixture models
    # do not completely represent the complexity of real distributions of
    # log2 ratio data.
    
    # true for the easiest case, where there is a single, normally distributed
    # population
    z_scores = rnorm(1000, mean=0, sd=0.5)
    expect_true(length(get_maxima(z_scores)) == 1)
    
    # try when there is another peak, with sufficient bulk and distant from the
    # first peak
    z_scores = c(z_scores, rnorm(100, mean=5, sd=0.1))
    expect_true(length(get_maxima(z_scores)) == 2)
    
    # check that adding in miniscule peaks only returns the major peaks
    z_scores = c(z_scores, rnorm(10, mean=10, sd=0.5))
    expect_true(length(get_maxima(z_scores)) == 2)
    
    # and add in another smaller peak near the secondary peak, so there are
    # three peaks which exceed the significance criteria, but since two of the
    # small peaks are close to each other, they get collapsed into a single
    # maxima.
    z_scores = c(z_scores, rnorm(100, mean=5.4, sd=0.1))
    expect_true(length(get_maxima(z_scores)) == 2)
})

test_that("get_null_parameters is correct", {
    
    z_scores = list()
    
    # check the null parameters for a normally distributed with mean of zero, and sd of 1
    z_scores = rnorm(4000, mean=0, sd=1)
    parameters = get_null_parameters(z_scores)
    expect_true(parameters$null_mean < 0.40 & parameters$null_mean > -0.40)
    expect_true(parameters$null_sd < 1.50 & parameters$null_sd > 0.50)
    
    # check the null parameters for a population with mean of zero, and sd of 1
    z_scores = rnorm(4000, mean=1, sd=2)
    parameters = get_null_parameters(z_scores)
    expect_true(parameters$null_mean < 1.40 & parameters$null_mean > 0.60)
    expect_true(parameters$null_sd < 3.00 & parameters$null_sd > 1.00)
    
    # check a distribution where we have a mixture of two normally distributed
    # models
    z_scores = c(rnorm(4000, mean=0, sd=1), rnorm(800, mean=4, sd=1))
    parameters = get_null_parameters(z_scores)
    expect_true(parameters$null_mean < 0.40 & parameters$null_mean > -0.40)
    expect_true(parameters$null_sd < 1.50 & parameters$null_sd > 0.50)
    
    # check a distribution where we have a mixture of multiple normally
    # distributed models
    z_scores = c(rnorm(800, mean=0, sd=1), rnorm(4000, mean=4, sd=1),
        rnorm(800, mean=8, sd=1))
    parameters = get_null_parameters(z_scores)
    expect_true(parameters$null_mean < 4.50 & parameters$null_mean > 3.50)
    expect_true(parameters$null_sd < 1.50 & parameters$null_sd > 0.50)
    
    # check that we get an error with distributions where the mixture model
    # estimation is unable to cope
    # NOTE: I still haven't tracked down the precise circumstances where mixture
    # NOTE: models fail.
    z_scores = c(runif(4000, -2, 2))
    parameters = get_null_parameters(z_scores)
})
