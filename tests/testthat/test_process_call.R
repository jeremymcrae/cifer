

library(cifer)
library(testthat)

context("CNV call processing")

test_that("process_cnv_call is correct", {
    
    cohort_n = 500
    sample_ids = paste("sample", 1:cohort_n, sep="_")
    
    samples = data.frame("individual_id"=sample_ids,
        "is_proband"=c(rep(FALSE, length(sample_ids) - 1), TRUE))
    probes = read.table(header=TRUE, text="
            probe A B C D E F G H I J K
            p1    3 5 4 4 5 4 4 5 2 3 9
            p2    3 2 3 3 2 4 3 4 2 5 8
            p3    2 3 3 2 3 3 3 3 2 4 9
            p4    4 4 4 2 2 3 2 5 3 2 9")
    
    probes = data.frame(matrix(rnorm(length(sample_ids) * 5), nrow=5))
    names(probes) = sample_ids
    
    probes[sample_ids[length(sample_ids)]] = rnorm(5, mean=10, sd=1)
    
    child_id = sample_ids[length(sample_ids)]
    mom_id = sample_ids[1]
    dad_id = sample_ids[2]
    
    values = process_cnv_call(samples, probes, child_id, mom_id, dad_id)
    
    expect_equal(values$inheritance, "not_inherited")
})
