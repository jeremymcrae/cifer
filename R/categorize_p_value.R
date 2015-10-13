#' categorize whether a p-value belongs to the null distribution
#'
#' @param p_value probability that sample fits within the null distribution
#' @param null_cutoff threshold at which to assign the sample as belonging to
#'     the null distribution
#' @param uncertain_cutoff threshold at which to assign the sample as belonging to
#'     the null distribution
#' @export
#' @return string of "reject" if the p-value belongs to the null distribution,
#'     "null" if the p-value belongs to the null distribution, and
#'     "uncertain" if it is unclear whether the p-value belongs to the null
#'     distribution or not.
#' @examples
#' categorize_p_value(0.000001)
#' categorize_p_value(0.1)
#' categorize_p_value(NA)
#' categorize_p_value(0.001, uncertain_cutoff=0.002)
categorize_p_value <- function(p_value, null_cutoff=0.005, uncertain_cutoff=0.0005) {
    
    # warn if the null and uncertain cutoff have been modified, and are
    # inconsistent
    stopifnot(null_cutoff > uncertain_cutoff)
    
    # return "uncertain" when the p-value is NA, for example when one of the parents
    # does not have data
    if (is.na(p_value)) { return("uncertain") }
    
    # figure out the parents CNV null model status
    category = "reject"
    if (p_value > uncertain_cutoff) {category = "uncertain"}
    if (p_value > null_cutoff) {category = "null"}
    
    return(category)
}
