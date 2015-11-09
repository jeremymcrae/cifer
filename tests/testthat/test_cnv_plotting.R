

library(cifer)
library(testthat)

context("CNV plotting")

test_that("plot_cnv runs without errors", {
    
    population = rnorm(500, mean=0, sd=1)
    z_scores = list("population"=population, "mom"=10, "dad"=0, "child"=10)
    
    cnv = list("person_id"="sample_1", "CHROM"="1", "POS"=1000, "INFO.END"=2000,
        "Inheritance"="unknown")
    
    plot_cnv(z_scores, cnv)
})

test_that("process_cnv_call runs without errors when including CNV details", {
    
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
    
    cnv = list("person_id"="sample_1", "CHROM"="1", "POS"=1000, "INFO.END"=2000,
        "Inheritance"="unknown")
    
    values = process_cnv_call(samples, probes, child_id, mom_id, dad_id, cnv)
})
