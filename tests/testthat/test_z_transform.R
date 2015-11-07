
library(cifer)
library(testthat)

context("Z-transform data")

test_that("get_l2r_z_scores output is correct", {
    # construct a dataframe for two samples, one of which will have a Z
    # transformed mean close to -0.707106, the other close to 0.707106, and
    # trio members with Z transformed values close to 2.12.
    pop = data.frame(a=c(1,2,3,4,5), b=c(3,4,5,6,7))
    z_scores = get_l2r_z_scores(mom_data=c(5,6,7,8,9), dad_data=c(5,6,7,8,9), child_data=c(5,6,7,8,9), pop)
    
    # check that the Z transformed values are close to their expected values.
    # Note that I don't include the full precision, since the values have 22
    # significant figures, so it's easier to check if they are close to the
    # predicted values.
    expect_true(z_scores$population[["a"]] - -0.707 < 0.01)
    expect_true(z_scores$population[["b"]] - 0.707 < 0.01)
    expect_true(z_scores$dad - 2.12 < 0.01)
    expect_true(z_scores$mom - 2.12 < 0.01)
    expect_true(z_scores$child - 2.12 < 0.01)
    
    # and check that the function copes with parents lacking probe data
    z_scores = get_l2r_z_scores(mom_data=NULL, dad_data=c(5,6,7,8,9), child_data=c(5,6,7,8,9), pop)
    expect_identical(z_scores$mom, NA)
})
