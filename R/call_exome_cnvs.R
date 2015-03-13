# find the inheritance status for CNVs called by CONVEX from exome data
# 

#' function to classify single CNV based on sample IDs, and CNV coordinates
#' 
#' @param proband_id sample ID for the proband (eg DDDP100001)
#' @param maternal_id sample ID for the proband's mother
#' @param paternal_id sample ID for the proband's father
#' @param chrom chromosome that the CNV is on (eg "1", "2", ..., "X")
#' @param start_pos start nucleotide of the CNV
#' @param stop_pos stop nucleotide of the CNV
#' @param cnv one row dataframe that contains the details of the CNV. Optional,
#'     only used if we need to plot the CNV, when we need to note any
#'     existing inheritance classification.
#' @param DATAFREEZE_DIR path to datafreeze directory containing files with family
#'     relationships and sanger IDs
#' @export
#' 
#' @return inheritance classification as a string eg "paternal_inh", "de_novo" etc
classify_exome_cnv <- function(proband_id, maternal_id, paternal_id, chrom, 
    start_pos, stop_pos, cnv=NA, DATAFREEZE_DIR="/nfs/ddd0/Data/datafreeze/1133trios_20131218") {
    
    ddd = get_ddd_individuals(DATAFREEZE_DIR)
    probes = try(get_probes_data(ddd, chrom, start_pos, stop_pos), silent = TRUE)
    
    # set the default return values, in case loading probe data gave an error
    prediction = list(inheritance="no_probe_data_for_CNV_region", mom_value=NA, dad_value=NA)
    
    # if we have loaded some probe scores, then try to predict inheritance
    if (class(probes) != "try-error") {
        prediction = process_cnv_call(ddd, probes, proband_id, maternal_id, paternal_id, cnv)
    } 
    
    return(prediction)
}

#' get a list of DDD participants, with their family relationships, and
#' sanger sample IDs
#' 
#' @param samples_dir path to datafreeze directory containing files with family
#'     relationships and sanger IDs
#' @export
#' 
#' @return dataframe with sample information, including sanger IDs
get_ddd_individuals <- function(samples_dir) {
    
    families_path = file.path(samples_dir, "family_relationships.shared.txt")
    families = read.table(families_path, sep = "\t", header = TRUE, stringsAsFactors=FALSE)
    
    sanger_ids_path = file.path(samples_dir, "person_sanger_decipher.private.txt")
    sanger_ids = read.table(sanger_ids_path, sep = "\t", header = TRUE, stringsAsFactors=FALSE)
    
    ddd = merge(families, sanger_ids, by.x = "individual_id", by.y = 
        "stable_id", all.x = TRUE)
    ddd$is_proband = ddd$dad_id != 0
    
    return(ddd)
}

#' get Z transformed L2R scores in the population, and in the proband family
#' 
#' @param mother l2r values for the mother of the current proband, or NULL if 
#'     the mother is not present
#' @param father l2r values for the father of the current proband, or NULL if
#'     the father is not present
#' @param proband l2r values for the current proband
#' @param population l2r values for parents not in proband's family
#' @export
#' 
#' @return list of z scores for population (as vector), and z scores for proband
#'     family (as list)
#' @examples
#' get_l2r_z_scores(c(1,2), c(0,3), c(3,3), data.frame(a=c(0.5, 0.1), b=c(0.3, 0.2)))
get_l2r_z_scores <- function(mother, father, proband, population) {
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
    if (!is.null(mother)) {mom = (mean(mother, na.rm = TRUE) - pop_mean)/pop_sd}
    if (!is.null(father)) {dad = (mean(father, na.rm = TRUE) - pop_mean)/pop_sd}
    proband = (mean(proband, na.rm = TRUE) - pop_mean)/pop_sd
    
    z_scores = list(population = population, mom = mom, dad = dad, proband = proband)
    
    return(z_scores)
}

#' estimate local maxima within a numeric distribution
#' 
#' Note that this function tends to overestimate the number of maxima in the
#' distribution, which is compensated for by only picking the biggest of the 
#' peaks later on.
#' 
#' @param z_scores numeric vector of Z-transformed log2-ratio data for parental
#'     population
#' @export
#' 
#' @return vector of maxima positions
#' @examples
#' get_maxima(rnorm(100, mean=0, sd=1))
#' get_maxima(c(rnorm(100, mean=0, sd=1), rnorm(20, mean=4, sd=1)))
get_maxima <- function(z_scores) {
    # figure out if we need a mixture model by examining the local maxima of the 
    # density
    dens = density(z_scores)
    maxima = which(diff(sign(diff(dens$y))) == -2) + 1
    
    # exclude maxima that are simply blips, and thus do not contribute greatly,
    # and drop out single maxima from pairs that are too close together
    maxima = maxima[dens$y[maxima] / max(dens$y[maxima]) > 0.05]
    if (length(maxima) > 1) {
        close_points = diff(dens$x[maxima]) < 0.7
        if (any(close_points)) {
            maxima = maxima[-c(which(close_points) + 1)]
        }
    }
    
    # convert the maxima back into their original population values
    maxima = dens$x[maxima]
    
    return(maxima)
}

#' get the parameters of the null distrubution
#' 
#' @param z_scores Z score transformed log2 ratio data for parents unrelated 
#'     to the proband currently being classified.
#' @export
#' 
#' @return a list containing the null model's mean and standard deviation, or
#'     raises an error if generating a mixture model and the code has too 
#'     many retries
#' @examples
#' get_null_parameters(rnorm(100, mean=0, sd=1))
#' get_null_parameters(c(rnorm(100, mean=0, sd=1), rnorm(30, mean=5, sd=1)))
get_null_parameters <-function(z_scores) {
    
    z_scores = z_scores[!is.na(z_scores)]
    
    # model the density, to figure out if we have multiple mixture models. 
    maxima = get_maxima(z_scores)
    
    # if there is only one maxima, then the vast majority of the parental 
    # population is tightly and normally distributed around a single mode, and 
    # the population mean and sd are good estimates for the parameters of the 
    # mode.
    if (length(maxima) == 1) {
        null_mean = mean(z_scores)
        null_sd = sd(z_scores)
    } else {
        # use a mixture distribution with > two local maxima, as this will allow
        # for non-rare CNVs that might be distributed differently from the null.
        # Note that this occurs when the blips are prevalent enough to affect
        # fitting a normal distribution, ie has local maxima greater than 5% of
        # the peak maxima.
        l2r_model = mixtools::normalmixEM(z_scores, k = length(maxima), maxrestarts = 50)
        
        # use the model closest to 0 as the null distribution, or could use the
        # model that forms the greatest proportion of the population
        # null_pos = which.min(abs(l2r_model$mu))
        null_pos = which.max(l2r_model$lambda)
        null_mean = l2r_model$mu[null_pos]
        null_sd = l2r_model$sigma[null_pos]
    }
    
    parameters = list(null_mean = null_mean, null_sd = null_sd)
    
    return(parameters)
}

#' quantify the probability that the parental CNVs belong to the null 
#' (uninherited) distribution
#' 
#' @param z_scores list of Z scores of L2R for the parental population, and 
#'     Z scores for the mother, father and current proband
#' @export
#' @return list of inheritance classification, and P values for the mother and 
#'     father, where the values indicate how likely each parents data belongs 
#'     to the uninherited cluster.
predict_inheritance <- function(z_scores) {
    
    parameters = try(get_null_parameters(z_scores$population), silent = TRUE)
    
    # sometimes the mixture model is unable to get a good fix on the model,
    # and raises an error, we capture those errors, and simply return a 
    # invalid classification.
    if (class(parameters) == "try-error") {
        return(list(inheritance = "unable_to_evaluate_probes", mom_value = NA, dad_value = NA))
    }
    
    null_mean = parameters$null_mean
    null_sd = parameters$null_sd
    
    # estimate the probability of getting the parental Z scores, given the 
    # population null distribution
    paternal_z = (z_scores$dad - null_mean)/null_sd
    dad_value = 2 * pnorm(-abs(paternal_z))
    maternal_z = (z_scores$mom - null_mean)/null_sd
    mom_value = 2 * pnorm(-abs(maternal_z))
    proband_z = (z_scores$proband - null_mean)/null_sd
    proband_value = 2 * pnorm(-abs(proband_z))
    
    inheritance = classify_inheritance(mom_value, dad_value, proband_value)
    family = list(mom_value = mom_value, dad_value = dad_value, 
        proband_value = proband_value, proband_z_score = proband_z)
    values = list(inheritance = inheritance, family = family)
    
    return(values)
}

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

#' classify the inheritance state of the childs CNV call
#' 
#' @param mom_p_value p-value for mom's scores belonging to the null model
#' @param dad_p_value p-value for dad's scores belonging to the null model
#' @param proband_p_value p-value for proband's scores belonging to the null model
#' @export
#' @return inheritance classification for the CNV eg("de_novo", "uncertain", 
#'     "paternal_inh" etc)
classify_inheritance <- function(mom_p_value, dad_p_value, proband_p_value) {
    
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

#' predict the inheritance state of the CNV call
#' 
#' @param ddd dataframe listing sample IDs, and file paths
#' @param probes dataframe of log-2-ratio (or adm3 score) values for all DDD
#'     participants, for the exome probes that lie within the CNV region.
#' @param proband_id ID of the proband
#' @param maternal_id DDD ID for the proband's mother
#' @param paternal_id DDD ID for the proband's father
#' @param cnv row information for CNV
#' @export
#' @return list of inheritance classification, the mother's p-value and the 
#'     father's p-value.
process_cnv_call <- function(ddd, probes, proband_id, maternal_id, paternal_id, cnv=NA) {
    
    # make sure all the parents are included in the probe dataset, so that we 
    # can order the parental z scores and correlations correctly.
    parent_ids = unique(ddd[!ddd$is_proband, ]$individual_id)
    parent_ids = parent_ids[!(parent_ids %in% c(maternal_id, paternal_id))]
    missing_parents = parent_ids[!parent_ids %in% names(probes)]
    for (parent in missing_parents) {
        probes[[parent]] = NA
    }
    
    # get a correctly ordered dataframe of parental L2R values
    parent_set = subset(probes, select=parent_ids)
    
    # find the Z scores for the population, and the family of interest
    z_scores = get_l2r_z_scores(probes[[maternal_id]], probes[[paternal_id]], 
        probes[[proband_id]], parent_set)
    
    # predict the inheritance state of the childs CNV
    prediction = predict_inheritance(z_scores)
    
    if (!is.na(cnv)) { include_graphs(ddd, probes, z_scores, proband_id, maternal_id, 
            paternal_id, cnv)
    }
    
    return(prediction)  
}

