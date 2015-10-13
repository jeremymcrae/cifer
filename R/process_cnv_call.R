#' predict the inheritance state of a CNV call
#'
#' @param samples dataframe listing sample IDs, and file paths
#' @param probes dataframe of log-2-ratio (or adm3 score) values for all
#'        participants, for the exome probes that lie within the CNV region.
#' @param child_id ID of the proband
#' @param mom_id sample ID for the proband's mother
#' @param dad_id sample ID for the proband's father
#' @param cnv row information for CNV
#' @export
#'
#' @return list of inheritance classification, the mother's p-value and the
#'        father's p-value.
#'
#' @examples
#' samples = read.table(header=TRUE, text="
#'        individual_id is_proband
#'        A             FALSE
#'        B             FALSE
#'        C             FALSE
#'        D             FALSE
#'        E             TRUE
#'        F             FALSE
#'        G             FALSE
#'        H             FALSE
#'        I             FALSE
#'        J             FALSE
#'        K             TRUE
#'        L             FALSE",
#'        colClasses=c("character", "logical"))
#' probes = read.table(header=TRUE, text="
#'        probe A B C D E F G H I J K
#'        p1    5 6 4 5 5 6 5 5 4 6 9
#'        p2    5 6 4 5 5 6 5 5 4 6 9
#'        p3    5 6 4 5 5 6 5 5 4 6 9
#'        p4    5 6 4 5 5 6 5 5 4 6 9")
#' child_id = "K"
#' mom_id = "J"
#' dad_id = "I"
#' process_cnv_call(samples, probes, child_id, mom_id, dad_id)
#'
#' # run another example, this time with a larger dataset
#' cohort_n = 500
#' sample_ids = paste("sample", 1:cohort_n, sep="_")
#' samples = data.frame("individual_id"=sample_ids,
#'     "is_proband"=c(rep(FALSE, length(sample_ids) - 1), TRUE))
#'
#' # define the population as having probes values cenetred around zero
#' probes = data.frame(matrix(rnorm(length(sample_ids) * 5), nrow=5))
#' names(probes) = sample_ids
#'
#' # set the child probe values to distant from the population values
#' probes[sample_ids[length(sample_ids)]] = rnorm(5, mean=10, sd=1)
#'
#' # define the sample IDs for the trio members
#' child_id = sample_ids[length(sample_ids)]
#' mom_id = sample_ids[1]
#' dad_id = sample_ids[2]
#'
#' process_cnv_call(samples, probes, child_id, mom_id, dad_id)
process_cnv_call <- function(samples, probes, child_id, mom_id, dad_id, cnv=NA) {
    
    # make sure all the parents are included in the probe dataset, so that we
    # can order the parental z scores and correlations correctly.
    parent_ids = unique(samples[!samples$is_proband, ]$individual_id)
    
    # drop out the parents for the current child, so we can be sure the population
    # parameters aren't biased by the child's parents.
    parent_ids = parent_ids[!(parent_ids %in% c(mom_id, dad_id))]
    missing_parents = parent_ids[!parent_ids %in% names(probes)]
    for (parent in missing_parents) {
        probes[[parent]] = NA
    }
    
    family = list("mom"=probes[[mom_id]], "dad"=probes[[dad_id]], "child"=probes[[child_id]])
    
    # find the Z scores for the population, and the family of interest
    z_scores = get_l2r_z_scores(family, probes[, parent_ids])
    population = z_scores$population
    family = z_scores$family
    
    # predict the inheritance state of the childs CNV
    prediction = predict_inheritance(population, family)
    
    if (!is.na(cnv)) {
        include_graphs(samples, probes, z_scores, child_id, mom_id, dad_id, cnv)
    }
    
    return(prediction)
}
