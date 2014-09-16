# unit testing for the CNV classification functions

library(cifer)
library(testthat)


context("CNV classification")

test_that("check categorise_p_value output", {
    
    # check the threshold boundaries of the p-value categorisation
    expect_identical(categorize_p_value(0.000099), "reject")
    expect_identical(categorize_p_value(0.0001), "reject")
    expect_identical(categorize_p_value(0.00011), "uncertain")
    expect_identical(categorize_p_value(0.00499), "uncertain")
    expect_identical(categorize_p_value(0.005), "uncertain")
    expect_identical(categorize_p_value(0.0051), "null")
    
    # and also allow for user-specified cutoffs
    expect_identical(categorize_p_value(0.01, null_cutoff=0.1, 
        uncertain_cutoff=0.02), "reject")
    
    # and catch potential error with user-specified cutoffs
    expect_error(categorize_p_value(0.01, null_cutoff=0.001, 
        uncertain_cutoff=0.02))
})

test_that("check classify_inheritance output", {
    
    # check the different inheritance classifications, given the trio p-values
    expect_identical(classify_inheritance(0.1, 0.1, 0.00001), "not_inherited")
    expect_identical(classify_inheritance(0.1, 0.1, 0.001), "uncertain")
    expect_identical(classify_inheritance(0.1, 0.1, 0.1), "false_positive")
    
    # check the different uncertain categories
    expect_identical(classify_inheritance(0.001, 0.1, 0.0001), "uncertain")
    expect_identical(classify_inheritance(0.1, 0.001, 0.0001), "uncertain")
    expect_identical(classify_inheritance(0.001, 0.001, 0.0001), "uncertain")
    
    # check the different inheritance categories
    expect_identical(classify_inheritance(0.0001, 0.001, 0.0001), "maternal_inh")
    expect_identical(classify_inheritance(0.0001, 0.01, 0.0001), "maternal_inh")
    expect_identical(classify_inheritance(0.001, 0.0001, 0.0001), "paternal_inh")
    expect_identical(classify_inheritance(0.01, 0.0001, 0.0001), "paternal_inh")
    expect_identical(classify_inheritance(0.0001, 0.0001, 0.0001), "biparental_inh")
})

test_that("check get_maxima output", {
    # true for the easiest case, where there is a single, normally distributed
    # population
    z_scores = rnorm(1000, mean=0, sd=1)
    expect_true(length(get_maxima(z_scores)) == 1)
    
    # try whe there is another peak, with sufficient bulk and distant from the
    # first peak
    z_scores = c(z_scores, rnorm(100, mean=5, sd=1))
    expect_true(length(get_maxima(z_scores)) == 2)
    
    # check that adding in miniscule peaks only returns the major peaks
    z_scores = c(z_scores, rnorm(10, mean=10, sd=1))
    expect_true(length(get_maxima(z_scores)) == 2)
    
    # and add in another smaller peak near the secondary peak, so there are
    # three peaks which exceed the significance criteria, but since two of the
    # small peaks are close to each other, they get collapsed into a single
    # maxima.
    z_scores = c(z_scores, rnorm(100, mean=5.5, sd=1))
    expect_true(length(get_maxima(z_scores)) == 2)
})

test_that("check get_null_parameters output", {
    
    z_scores = list()
    z_scores$population = rnorm(1000, mean=0, sd=1)
    
    parameters = get_null_parameters(z_scores)
    expect_true(parameters$null_mean < 0.05 & parameters$null_mean > -0.05)
    expect_true(parameters$null_sd < 1.05 & parameters$null_sd > 0.95)
    
    z_scores$population = rnorm(1000, mean=1, sd=2)
    parameters = get_null_parameters(z_scores)
    expect_true(parameters$null_mean < 1.05 & parameters$null_mean > 0.95)
    expect_true(parameters$null_sd < 2.1 & parameters$null_sd > 1.9)
    
    
    
})