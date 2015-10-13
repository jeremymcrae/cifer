

library(cifer)
library(testthat)

context("CNV call processing")

test_that("process_cnv_call is correct", {
    
    cohort_n = 500
    sample_ids = paste("sample", 1:cohort_n, sep="_")
    
    samples = data.frame("individual_id"=sample_ids,
        "is_proband"=c(rep(FALSE, length(sample_ids) - 1), TRUE))
    
    probes = data.frame(matrix(rnorm(length(sample_ids) * 5), nrow=5))
    names(probes) = sample_ids
    
    # set the child probe values to distant from the population values
    probes[sample_ids[length(sample_ids)]] = rnorm(5, mean=10, sd=1)
    
    # define the sample IDs for the trio members
    child_id = sample_ids[length(sample_ids)]
    mom_id = sample_ids[1]
    dad_id = sample_ids[2]
    
    values = process_cnv_call(samples, probes, child_id, mom_id, dad_id)
    
    expect_equal(values$inheritance, "not_inherited")
})
