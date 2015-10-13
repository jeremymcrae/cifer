#' get Z transformed L2R scores in the population, and in the proband family
#'
#' @param family list of l2r values for the current trio, ie values for mother,
#'     father and child. The parent values should be NULL if the parents do not
#'     have data available.
#' @param population l2r values for parents not in proband's family
#' @export
#'
#' @return list of population Z scores (as vector), and family Z scores (as list)
#'
#' @examples
#' family = list("mom"=c(1,2), "dad"=c(0,3), "child"=c(3,3))
#' population = data.frame(a=c(0.5, 0.1), b=c(0.3, 0.2))
#' get_l2r_z_scores(family, population)
get_l2r_z_scores <- function(family, population) {
    # determine the log2 ratio population distribution
    if(is.null(dim(population))) {
        population_l2r = population
    } else {
        population_l2r = colMeans(population, na.rm = TRUE)
    }
    pop_mean = mean(population_l2r, na.rm=TRUE)
    pop_sd = sd(population_l2r, na.rm=TRUE)
    
    # get the z-scores of mean log2 ratios for the population, and for the
    # trio members
    population = (population_l2r - pop_mean)/pop_sd
    mom = NA
    dad = NA
    # make sure the parents data exists (they might not be present in the
    # population) before trying to Z transform
    if (!is.null(family$mom)) {mom = (mean(family$mom, na.rm=TRUE) - pop_mean)/pop_sd}
    if (!is.null(family$dad)) {dad = (mean(family$dad, na.rm=TRUE) - pop_mean)/pop_sd}
    child = (mean(family$child, na.rm=TRUE) - pop_mean)/pop_sd
    
    family = list("mom"=mom, "dad"=dad, "child"=child)
    z_scores = list("population"=population, "family"=family)
    
    return(z_scores)
}
