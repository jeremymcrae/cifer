# find the inheritance status for CNVs called by CONVEX from exome data
# 

library(mixtools)

source("load_convex_data.R")

DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/1133trios_20131218"

classify_exome_cnv <- function(proband_id, maternal_id, paternal_id, chrom, start, stop) {
    # function to classify single CNV based on sample IDs, and CNV coordinates
    # 
    # Args:
    #     proband_id: sample ID for the proband (eg DDDP100001)
    #     maternal_id: sample ID for the proband's mother
    #     paternal_id: sample ID for the proband's father
    #     chrom: chromosome that the CNV is on (eg "1", "2", ..., "X")
    #     start: start nucleotide of the CNV
    #     stop: stop nucleotide of the CNV
    # 
    # Returns:
    #     inheritance classification as a string eg "paternal_inh", "de_novo" etc
    
    ddd = get_ddd_individuals(DATAFREEZE_DIR)
    probes = try(get_probes_data(ddd, chrom, start, stop), silent = TRUE)
    if (class(probes) != "try-error") {
        prediction = process_cnv_call(ddd, probes, proband_id, maternal_id, paternal_id)
    } else {
        prediction = list(inheritance = "no_probe_data_for_CNV_region", mom_value = NA, dad_value = NA)
    }
    
    return(prediction)
}

get_ddd_individuals <- function(samples_dir) {
    # get a list of DDD participants, with their family relationships, and
    # sanger sample IDs
    
    families_path = file.path(samples_dir, "family_relationships.shared.txt")
    families = read.table(families_path, sep = "\t", header = TRUE, stringsAsFactors=FALSE)
    
    sanger_ids_path = file.path(samples_dir, "person_sanger_decipher.private.txt")
    sanger_ids = read.table(sanger_ids_path, sep = "\t", header = TRUE, stringsAsFactors=FALSE)
    
    ddd = merge(families, sanger_ids, by.x = "individual_id", by.y = 
        "stable_id", all.x = TRUE)
    ddd$is_proband = ddd$dad_id != 0
    
    return(ddd)
}

get_l2r_z_scores <- function(mother, father, population_l2r_set) {
    # get the L2R Z scores in the population, and in the proband family
    # 
    # Args:
    #     mother: l2r values for the mother of the current proband
    #     father: l2r values for the father of the current proband
    #     population_l2r_set: l2r values for parents not in proband's family
    # 
    # Returns:
    #     list of z scores for population (as vector), and z scores for proband
    #     family (as list)
    
    # determine the log2 ratio population distribution
    population_l2r = as.vector(colMeans(population_l2r_set, na.rm = TRUE))
    pop_mean = mean(population_l2r, na.rm=TRUE)
    pop_sd = sd(population_l2r, na.rm=TRUE)
    
    # get the z-scores of mean log2 ratios for the population, and for the 
    # parents of the proband
    population_z_scores = (population_l2r - pop_mean)/pop_sd
    mom = (mean(mother, na.rm = TRUE) - pop_mean)/pop_sd
    dad = (mean(father, na.rm = TRUE) - pop_mean)/pop_sd
    
    z_scores = list(population_z_scores = population_z_scores, mom = mom, dad = dad)
    
    return(z_scores)
}

predict_inheritance <- function(z_scores, maternal_z_score, paternal_z_score) {
    # quantify the probability that the parental CNVs belong to the null 
    # (uniherited) distribution
    # 
    # Args:
    #     z_scores: Z scores of L2R for the parental population
    #     maternal_z_score: Z score for the mothers L2R
    #     paternal_z_score: Z score for the fathers L2R
    # 
    # Returns:
    #     list of, inheritance classification, and P values for the mother and 
    #     father, where the values indicate how likely each parents data belongs 
    #     to the uninherited cluster.
    
    z_scores = z_scores[!is.na(z_scores)]
    
    # model the density, to figure out if we have multiple mixture models. 
    # Figure out if we need mixture model by examining the local maxima of the 
    # density
    dens = density(z_scores)
    maxima = which(diff(sign(diff(dens$y))) == -2) + 1
    
    # exclude maxima that are simply blips, and do not contribute greatly, and
    # drop out single maxima from pairs that are too close together
    maxima = maxima[dens$y[maxima] / max(dens$y[maxima]) > 0.05]
    maxima = maxima[maxima - c(-100000, maxima) > 0.7]
    
    if (length(maxima) == 1) {
        null_mean = mean(z_scores)
        null_sd = sd(z_scores)
    } else {
        # use a mixture distribution with > two local maxima, as this will allow
        # for non-rare CNVs that might be distributed differently from the null.
        # Note that this occurs when the blips are prevalent enough to affect
        # fitting a normal distribution, ie has local maxima greater than 5% of
        # the peak maxima.
        l2r_model = try(normalmixEM(z_scores, k = length(maxima), maxrestarts = 100))
        
        # sometimes the mixture model is unable to get a good fix on the model,
        # and raises an error, we capture those errors, and simply return a 
        # invalid classification.
        if (class(l2r_model) == "try-error") {
            return(list(inheritance = "unable_to_evaluate_probes", mom_value = NA, dad_value = NA))
        }
        
        # use the model closest to 0 as the null distribution
        null_pos = which.min(abs(l2r_model$mu))
        null_mean = l2r_model$mu[null_pos]
        null_sd = l2r_model$sigma[null_pos]
    }
    
    # estimate the probability of getting the parental Z scores, given the 
    # population null distribution 
    paternal_z = (paternal_z_score - null_mean)/null_sd
    dad_value = 2 * pnorm(-abs(paternal_z))
    maternal_z = (maternal_z_score - null_mean)/null_sd
    mom_value = 2 * pnorm(-abs(maternal_z))
    
    inheritance = classify_inheritance(mom_value, dad_value)
    values = list(inheritance = inheritance, mom_value = mom_value, dad_value = dad_value)
    
    return(values)
}

classify_inheritance <- function(mom_value, dad_value) {
    # classify the inheritance state of the childs CNV call
    # 
    # Args:
    #     mom_value: p-value for mom's scores belonging to the null model
    #     dad_value: p-value for dad's scores belonging to the null model
    # 
    # Returns:
    #     inheritance classification for the CNV eg("de_novo", "uncertain", 
    #     "paternal_inh" etc)
    
    uncertain_cutoff = 0.0001
    null_cutoff = 0.005
    
    # figure out the parents CNV null model status
    mom = "reject"
    if (mom_value > uncertain_cutoff) {mom = "uncertain"}
    if (mom_value > null_cutoff) {mom = "null"}
    
    # figure out the parents CNV null model status
    dad = "reject"
    if (dad_value > uncertain_cutoff) {dad = "uncertain"}
    if (dad_value > null_cutoff) {dad = "null"}
    
    # classify the child's CNV inheritance based on the parental statuses
    inh = NA
    if (dad == "null" & mom == "null") {inh = "de_novo"}
    if (dad == "uncertain" & mom == "null") {inh = "uncertain"}
    if (dad == "null" & mom == "uncertain") {inh = "uncertain"}
    if (dad == "uncertain" & mom == "uncertain") {inh = "uncertain"}
    if (dad == "uncertain" & mom == "reject") {inh = "maternal_inh"}
    if (dad == "reject" & mom == "uncertain") {inh = "paternal_inh"}
    if (dad == "reject" & mom == "reject") {inh = "biparental_inh"}
    if (dad == "null" & mom == "reject") {inh = "maternal_inh"}
    if (dad == "reject" & mom == "null") {inh = "paternal_inh"}
    
   return(inh)
}

process_cnv_call <- function(ddd, probes, person_id, maternal_id, paternal_id, row=NA, graphs=FALSE) {
    # predict the inheritance state of the CNV call
    # 
    # Args:
    #     probes: dataframe of log-2-ratio (or adm3 score) values for all DDD
    #         participants, for the exome probes that lie within the CNV region.
    #     person_id: ID of the proband
    #     maternal_id: DDD ID for the proband's mother
    #     paternal_id: DDD ID for the proband's father
    #     row: row information for CNV
    #     graphs: true/false for whether to plot graphs
    # 
    # Returns:
    #     list of inheritance classification, the mother's p-value and the 
    #     father's p-value.
    
    # make sure all the parents are included in the probe dataset, so that we 
    # can order the parental z scores and correlations correctly.
    parent_ids = unique(ddd[!ddd$is_proband, ]$individual_id)
    missing_parents = parent_ids[!parent_ids %in% names(probes)]
    for (parent in missing_parents) {
        probes[[parent]] = NA
    }
    
    # get a correctly ordered dataframe of parental L2R values
    parent_set = subset(probes, select = parent_ids)
    
    # find the Z scores for the population, and the family of interest
    z_scores = get_l2r_z_scores(probes[[maternal_id]], probes[[paternal_id]], parent_set)
    population_z_scores = z_scores$population_z_scores
    maternal_z_score = z_scores$mom
    paternal_z_score = z_scores$dad
    
    # predict the inheritance state of the childs CNV
    predictions = predict_inheritance(population_z_scores, maternal_z_score, paternal_z_score)
    
    if (graphs) {
        include_graphs(ddd, probes, population_z_scores, maternal_z_score,
            paternal_z_score, person_id, maternal_id, paternal_id, row)
    }
    
    return(predictions)  
}

