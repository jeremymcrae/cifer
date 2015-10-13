#' classify the inheritance state of the child's CNV call
#'
#' @param mom_p_value p-value for mom's scores belonging to the null model
#' @param dad_p_value p-value for dad's scores belonging to the null model
#' @param proband_p_value p-value for proband's scores belonging to the null model
#' @export
#'
#' @return inheritance classification for the CNV eg("de_novo", "uncertain",
#'         "paternal_inh" etc)
#'
#' @examples
#' mom_p = 0.1
#' dad_p = 0.1
#' child_p = 0.0001
#' classify_from_trio_p_values(mom_p, dad_p, child_p)
classify_from_trio_p_values <- function(mom_p_value, dad_p_value, proband_p_value) {
    
    # categorise the trio members as "null", "uncertain" or "reject"
    mom = categorize_p_value(mom_p_value)
    dad = categorize_p_value(dad_p_value)
    proband = categorize_p_value(proband_p_value)
    
    # classify the child's CNV inheritance based on the parental statuses
    inh = NA
    if (dad == "null" & mom == "null") {inh = "not_inherited"}
    if (dad == "uncertain" & mom == "null") {inh = "uncertain"}
    if (dad == "null" & mom == "uncertain") {inh = "uncertain"}
    if (dad == "uncertain" & mom == "uncertain") {inh = "uncertain"}
    if (dad == "uncertain" & mom == "reject") {inh = "maternal_inh"}
    if (dad == "reject" & mom == "uncertain") {inh = "paternal_inh"}
    if (dad == "reject" & mom == "reject") {inh = "biparental_inh"}
    if (dad == "null" & mom == "reject") {inh = "maternal_inh"}
    if (dad == "reject" & mom == "null") {inh = "paternal_inh"}
    
    # if the variant is in the not inherited category, or when we lack data for
    # both parents, we check the probands status to see if we can discount the
    # CNV as a false positive, since many apparent de novo CNV calls are false
    # positives rather than true de novos
    if (inh == "not_inherited" | (is.na(mom_p_value) & is.na(dad_p_value))) {
        if (proband == "null") {
            inh = "false_positive"
        } else if (proband == "uncertain") {
            inh = "uncertain"
        }
    }
    
   return(inh)
}
