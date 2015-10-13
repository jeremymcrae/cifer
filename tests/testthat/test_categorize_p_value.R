
library(cifer)
library(testthat)

context("p-value classification")

test_that("categorize_p_value is correct", {
    
    # check the threshold boundaries of the p-value categorisation
    expect_identical(categorize_p_value(0.000099), "reject")
    expect_identical(categorize_p_value(0.0005), "reject")
    expect_identical(categorize_p_value(0.00051), "uncertain")
    expect_identical(categorize_p_value(0.00499), "uncertain")
    expect_identical(categorize_p_value(0.005), "uncertain")
    expect_identical(categorize_p_value(0.0051), "null")
    
    # and also allow for user-specified cutoffs
    expect_identical(categorize_p_value(0.01, null_cutoff=0.1,
        uncertain_cutoff=0.02), "reject")
    
    # and catch potential error with user-specified cutoffs
    expect_error(categorize_p_value(0.01, null_cutoff=0.001,
        uncertain_cutoff=0.02))
    
    # and make sure we return "uncertain" when the parents lack data, and have
    # NA p-values
    expect_equal(categorize_p_value(NA), "uncertain")
})
