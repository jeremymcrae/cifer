
library(cifer)
library(testthat)

context("inheritance prediction")

test_that("predict_inheritance is correct", {
    
    # construct a null distribution with mean of 0, and sd of 1. Then define
    # trio values such that the mom and proband are outliers
    population = rnorm(500, mean=0, sd=1)
    family = list("mom"=10, "dad"=0, "child"=10)
    
    # since the mother is an outlier, we expect maternal inheritance
    expect_identical(predict_inheritance(population, family)$inheritance, "maternal_inh")
})

test_that("predict_inheritance is correct for biparental inheritance", {
    
    # make both parents outliers, for biparental inheritance
    population = rnorm(500, mean=0, sd=1)
    family = list("mom"=10, "dad"=10, "child"=10)
    
    expect_identical(predict_inheritance(population, family)$inheritance, "biparental_inh")
})

test_that("predict_inheritance is correct for paternal inheritance", {
    # shift the mother back to the null, leaving only the father as outlier, for
    # paternal inheritance
    population = rnorm(500, mean=0, sd=1)
    family = list("mom"=0, "dad"=10, "child"=10)
    
    expect_identical(predict_inheritance(population, family)$inheritance, "paternal_inh")
})

test_that("predict_inheritance is correct for de novo mutation", {
    # if both parents are within the null, then if the probands value is an
    # outlier, the variant is not inherited
    population = rnorm(500, mean=0, sd=1)
    family = list("mom"=0, "dad"=0, "child"=10)
    
    expect_identical(predict_inheritance(population, family)$inheritance, "not_inherited")
})

test_that("predict_inheritance is correct for false positives", {
    # if both parents are within the null, as well as the probands value, then
    # we have a false_positive
    population = rnorm(500, mean=0, sd=1)
    family = list("mom"=0, "dad"=0, "child"=0)
    
    expect_identical(predict_inheritance(population, family)$inheritance, "false_positive")
})

test_that("predict_inheritance is correct for mixture distribution", {
    # construct a mixture distribution, and where only the father is an outlier
    # to the null, but within the secondary peak
    population = c(rnorm(500, mean=0, sd=1), rnorm(100, mean=5, sd=1))
    family = list("mom"=0, "dad"=5, "child"=5)
    
    expect_identical(predict_inheritance(population, family)$inheritance, "paternal_inh")
})

test_that("predict_inheritance is correct for single NA parent", {
    # check that it can cope with NA parents
    population = rnorm(500, mean=0, sd=1)
    family = list("mom"=NA, "dad"=5, "child"=5)
    
    expect_identical(predict_inheritance(population, family)$inheritance, "paternal_inh")
})

test_that("predict_inheritance is correct for both parents NA", {
    # check it can cope when both parents are NA
    population = rnorm(500, mean=0, sd=1)
    family = list("mom"=NA, "dad"=NA, "child"=5)
    
    expect_identical(predict_inheritance(population, family)$inheritance, "uncertain")
})

test_that("predict_inheritance is correct for both parents NA and child in null", {
    # check correctness when both parents are NA, and the proband's value lies
    # within the null distribution
    population = rnorm(500, mean=0, sd=1)
    family = list("mom"=NA, "dad"=NA, "child"=0)
    
    expect_identical(predict_inheritance(population, family)$inheritance, "false_positive")
})
