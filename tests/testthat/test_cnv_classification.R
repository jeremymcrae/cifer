# unit testing for the CNV classification functions

library(cifer)
library(testthat)

context("CNV classification")


test_that("get_l2r_z_scores output is correct", {
    # construct a dataframe for two samples, one of which will have a Z 
    # transformed mean close to -0.707106, the other close to 0.707106, and 
    # trio members with Z transformed values close to 2.12.
    pop = data.frame(a = c(1,2,3,4,5), b = c(3,4,5,6,7))
    mother = c(5,6,7,8,9)
    father = c(5,6,7,8,9)
    proband = c(5,6,7,8,9)
    z_scores = get_l2r_z_scores(mother, father, proband, pop)
    
    # check that the Z transformed values are close to their expected values.
    # Note that I don't include the full precision, since the values have 22
    # significant figures, so it's easier to check if they are close to the 
    # predicted values.
    expect_true(z_scores$population[["a"]] - -0.707 < 0.01)
    expect_true(z_scores$population[["b"]] - 0.707 < 0.01)
    expect_true(z_scores$dad - 2.12 < 0.01)
    expect_true(z_scores$mom - 2.12 < 0.01)
    expect_true(z_scores$proband - 2.12 < 0.01)
    
    # and check that the function copes with parents lacking probe data
    mother = NULL
    z_scores = get_l2r_z_scores(mother, father, proband, pop)
    expect_identical(z_scores$mom, NA)
})

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

test_that("predict_inheritance is correct", {
    
    z_scores = list()
    
    # construct a null distribution with mean of 0, and sd of 1. Then define
    # trio values such that the mom and proband are outliers
    z_scores$population = rnorm(500, mean=0, sd=1)
    z_scores$mom = 10
    z_scores$dad = 0
    z_scores$proband = 10
    
    # since the mother is an outlier, we expect maternal inheritance
    expect_identical(predict_inheritance(z_scores)$inheritance, "maternal_inh")
    
    # make both parents outliers, for biparental inheritance
    z_scores$dad = 10
    expect_identical(predict_inheritance(z_scores)$inheritance, "biparental_inh")
    
    # shift the mother back to the null, leaving only the father as outlier, for
    # paternal inheritance
    z_scores$mom = 0
    expect_identical(predict_inheritance(z_scores)$inheritance, "paternal_inh")
    
    # if both parents are within the null, then if the probands value is an
    # outlier, the variant is not inherited 
    z_scores$dad = 0
    expect_identical(predict_inheritance(z_scores)$inheritance, "not_inherited")
    
    # if both parents are within the null, as well as the probands value, then
    # we have a false_positive
    z_scores$proband = 0
    expect_identical(predict_inheritance(z_scores)$inheritance, "false_positive")
    
    # construct a mixture distribution, and where only the father is an outlier 
    # to the null, but within the secondary peak
    z_scores$population = c(rnorm(500, mean=0, sd=1), rnorm(100, mean=5, sd=1))
    z_scores$dad = 5
    z_scores$proband = 5
    expect_identical(predict_inheritance(z_scores)$inheritance, "paternal_inh")
    
    # check that it can cope with NA parents
    z_scores$population = c(rnorm(500, mean=0, sd=1))
    z_scores$mom = NA
    expect_identical(predict_inheritance(z_scores)$inheritance, "paternal_inh")
    
    # check it can cope when both parents are NA
    z_scores$dad = NA
    expect_identical(predict_inheritance(z_scores)$inheritance, "uncertain")
    
    # check correctness when both parents are NA, and the proband's valie lies 
    # within the null distribution
    z_scores$proband = 0
    expect_identical(predict_inheritance(z_scores)$inheritance, "false_positive")
})

test_that("categorise_p_value is correct", {
    
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
    
    # and make sure we return "uncertain" when the parents lack data, and have
    # NA p-values
    expect_equal(categorize_p_value(NA), "uncertain")
})

test_that("classify_inheritance output is correct", {
    
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





