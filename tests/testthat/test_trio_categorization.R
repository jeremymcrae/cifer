
library(cifer)
library(testthat)

context("inheritance classification")

test_that("classify_classify_from_trio_p_values is correct", {
    
    # check the different inheritance classifications, given the trio p-values
    expect_identical(classify_from_trio_p_values(0.1, 0.1, 0.00001), "not_inherited")
    expect_identical(classify_from_trio_p_values(0.1, 0.1, 0.001), "uncertain")
    expect_identical(classify_from_trio_p_values(0.1, 0.1, 0.1), "false_positive")
})

test_that("classify_from_trio_p_values is correct for uncertain categories", {
    # check the different uncertain categories
    expect_identical(classify_from_trio_p_values(0.001, 0.1, 0.0001), "uncertain")
    expect_identical(classify_from_trio_p_values(0.1, 0.001, 0.0001), "uncertain")
    expect_identical(classify_from_trio_p_values(0.001, 0.001, 0.0001), "uncertain")
})

test_that("classify_from_trio_p_values is correct for inherited categories", {
    # check the different inheritance categories
    expect_identical(classify_from_trio_p_values(0.0001, 0.001, 0.0001), "maternal_inh")
    expect_identical(classify_from_trio_p_values(0.0001, 0.01, 0.0001), "maternal_inh")
    expect_identical(classify_from_trio_p_values(0.001, 0.0001, 0.0001), "paternal_inh")
    expect_identical(classify_from_trio_p_values(0.01, 0.0001, 0.0001), "paternal_inh")
    expect_identical(classify_from_trio_p_values(0.0001, 0.0001, 0.0001), "biparental_inh")
})
