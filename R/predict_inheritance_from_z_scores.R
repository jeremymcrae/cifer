#' predict the inheritance of a child's CNV from population and family scores
#'
#' @param population vector of Z scores of L2R for the parental population
#' @param family list of Z scores for "mom", "dad" and "child"
#' @export
#'
#' @return list of inheritance classification, and P values for the mother and
#'     father, where the values indicate how likely each parents data belongs
#'     to the uninherited cluster.
#'
#' @examples
#' family = list("mom"=0, "dad"=0, "child"=5)
#' population = rnorm(100)
#' predict_inheritance(population, family)
#'
#' family = list("mom"=NA, "dad"=NA, "child"=1)
#' population = rnorm(100)
#' predict_inheritance(population, family)
#' z_scores = list("population"=population, "mom"=NA, "dad"=NA, "child"=1)
#' predict_inheritance(z_scores)
predict_inheritance <- function(population, family=NULL) {
    
    # Previously the function allowed putting all the variables in a single list
    # variable, extract out the population and family data if this is the case.
    if ("population" %in% names(population) & is.null(family)) {
        family = population
        population = population$population
    }
    
    parameters = try(get_null_parameters(population), silent=TRUE)
    
    # sometimes the mixture model is unable to get a good fix on the model,
    # and raises an error, we capture those errors, and simply return a
    # invalid classification.
    if (class(parameters) == "try-error") {
        return(list("inheritance"="unable_to_evaluate_probes", "mom_value"=NA,
            "dad_value"=NA, "proband_value"=NA, "proband_z_score"=NA))
    }
    
    null_mean = parameters$null_mean
    null_sd = parameters$null_sd
    
    # estimate the probability of getting the parental Z scores, given the
    # population null distribution
    paternal_z = (family$dad - null_mean)/null_sd
    dad_value = 2 * stats::pnorm(-abs(paternal_z))
    maternal_z = (family$mom - null_mean)/null_sd
    mom_value = 2 * stats::pnorm(-abs(maternal_z))
    child_z = (family$child - null_mean)/null_sd
    child_value = 2 * stats::pnorm(-abs(child_z))
    
    inheritance = classify_from_trio_p_values(mom_value, dad_value, child_value)
    family = list("mom_value"=mom_value, "dad_value"=dad_value,
        "proband_value"=child_value, "proband_z_score"=child_z)
    values = list("inheritance"=inheritance, "family"=family)
    
    return(values)
}
