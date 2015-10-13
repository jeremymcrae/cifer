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
    
    samples = get_ddd_individuals(DATAFREEZE_DIR)
    probes = try(get_probes_data(samples, chrom, start_pos, stop_pos), silent = TRUE)
    
    # set the default return values, in case loading probe data gave an error
    prediction = list(inheritance="no_probe_data_for_CNV_region", mom_value=NA, dad_value=NA)
    
    # if we have loaded some probe scores, then try to predict inheritance
    if (class(probes) != "try-error") {
        prediction = process_cnv_call(samples, probes, proband_id, maternal_id, paternal_id, cnv)
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
    
    samples = merge(families, sanger_ids, by.x = "individual_id", by.y =
        "stable_id", all.x = TRUE)
    samples$is_proband = samples$dad_id != 0
    
    return(samples)
}

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
predict_inheritance <- function(population, family) {
    
    parameters = try(get_null_parameters(population), silent=TRUE)
    
    # sometimes the mixture model is unable to get a good fix on the model,
    # and raises an error, we capture those errors, and simply return a
    # invalid classification.
    if (class(parameters) == "try-error") {
        return(list(inheritance = "unable_to_evaluate_probes", mom_value=NA, dad_value=NA))
    }
    
    null_mean = parameters$null_mean
    null_sd = parameters$null_sd
    
    # estimate the probability of getting the parental Z scores, given the
    # population null distribution
    paternal_z = (family$dad - null_mean)/null_sd
    dad_value = 2 * pnorm(-abs(paternal_z))
    maternal_z = (family$mom - null_mean)/null_sd
    mom_value = 2 * pnorm(-abs(maternal_z))
    child_z = (family$child - null_mean)/null_sd
    child_value = 2 * pnorm(-abs(child_z))
    
    inheritance = classify_inheritance(mom_value, dad_value, child_value)
    family = list("mom_value"=mom_value, "dad_value"=dad_value,
        "proband_value"=child_value, "proband_z_score"=child_z)
    values = list("inheritance"=inheritance, "family"=family)
    
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
#' mom = 0.1
#' dad = 0.1
#' child = 0.0001
#' classify_inheritance(mom, dad, child)
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
    
    if (!is.na(cnv)) { include_graphs(samples, probes, z_scores, child_id, mom_id,
            dad_id, cnv)
    }
    
    return(prediction)
}
