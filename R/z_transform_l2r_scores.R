#' get Z transformed L2R scores in the population, and in the proband family
#'
#' @param mom_data vector of l2r values for the mother, or NULL if unavailable
#' @param dad_data vector of l2r values for tthe father, or NULL if unavailable
#' @param child_data vector of l2r values for the child
#' @param population l2r values for parents not in proband's family
#' @export
#'
#' @return list of population Z scores (as vector), and family Z scores (as list)
#'
#' @examples
#' population = data.frame(a=c(0.5, 0.1), b=c(0.3, 0.2))
#' get_l2r_z_scores(mom_data=c(1,2), dad_data=c(0,3), child_data=c(3,3), population)
get_l2r_z_scores <- function(mom_data, dad_data, child_data, population) {
    # determine the log2 ratio population distribution
    if (is.null(dim(population))) {
        population_l2r = population
    } else {
        population_l2r = colMeans(population, na.rm = TRUE)
    }
    pop_mean = mean(population_l2r, na.rm=TRUE)
    pop_sd = stats::sd(population_l2r, na.rm=TRUE)
    
    # get the z-scores of mean log2 ratios for the population, and for the
    # trio members
    population = (population_l2r - pop_mean)/pop_sd
    mom = NA
    dad = NA
    # make sure the parents data exists (they might not be present in the
    # population) before trying to Z transform
    if (!is.null(mom_data)) {mom = (mean(unlist(mom_data), na.rm=TRUE) - pop_mean)/pop_sd}
    if (!is.null(dad_data)) {dad = (mean(unlist(dad_data), na.rm=TRUE) - pop_mean)/pop_sd}
    child = (mean(unlist(child_data), na.rm=TRUE) - pop_mean)/pop_sd
    
    z_scores = list("population"=population, "mom"=mom, "dad"=dad, "child"=child)
    
    return(z_scores)
}
