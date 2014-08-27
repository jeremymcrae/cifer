# functions to plot L2R scores against correlations

library(Cairo)

get_parental_correlations <- function(proband, mother, father) {
    # find the correlation in L2R values between a child and their parents
    # 
    # Args:
    #     proband: L2R values for proband (or NULL if not available)
    #     mother: L2R values for mother (or NULL if not available)
    #     father: L2R values for father (or NULL if not available)
    # 
    # Returns:
    #     list of correlations for mother and father's L2R scores with probands
    #     L2R values, or NA where L2R values are unavailable.
    
    mom = NA
    dad = NA
    # make sure that the values are available for checking correlation
    if (!is.null(proband)) {
        if (!is.null(mother)) { mom = cor(proband, mother) }
        if (!is.null(father)) { dad = cor(proband, father) }  
    }
    
    correlations = list(mom = mom, dad = dad)
    
    return(correlations)
}

include_graphs <- function(ddd, probes, person_id, maternal_id, paternal_id, row) {
    # set up the data for plotting a graph, then call the plotting function
    # 
    # Initially we expected that correlations between the proband's scores with
    # parental scores might be useful, which is why this function calculates the
    # correlation, before plotting.
    # 
    # Args:
    #     ddd: dataframe containing all the DDD sample information.
    #     probes: L2R or ADM3 scores for all individuals for probes within the CNV
    #     person_id: DDD ID of the proband
    #     maternal_id: DDD ID of the proband's mother
    #     paternal_id: DDD ID of the proband's father
    #     row: one row dataframe that contains the CNV's details (used to 
    #         identify the plot)
    # 
    # Returns:
    #     nothing
    
    # get the correlation values for the proband to their parents
    parental_corr = get_parental_correlations(probes[[person_id]], 
        probes[[maternal_id]], probes[[paternal_id]])
    paternal_corr = parental_corr$dad
    maternal_corr = parental_corr$mom
    
    # get the sample IDs for different subsets of the population
    proband_ids = unique(ddd[ddd$is_proband, ]$individual_id)
    proband_ids = proband_ids[proband_ids != person_id]
    
    # set up a vector of parental IDs, ordered by the proband ID, with mother ID
    # before father ID for each pair. This is so the vector of population Z 
    # scores and vector of population correlations are ordered identically.
    parent_ids = c()
    for (proband in proband_ids) {
        dad_id = unique(ddd[ddd$individual_id == proband, ]$dad_id)
        mom_id = unique(ddd[ddd$individual_id == proband, ]$mum_id)
        parent_ids = c(parent_ids, mom_id, dad_id)
    }
    
    # get the correlation between each proband and their parents L2R values
    population_corrs = c()
    for (proband in proband_ids) {
        dad_id = unique(ddd[ddd$individual_id == proband, ]$dad_id)
        mom_id = unique(ddd[ddd$individual_id == proband, ]$mum_id)
        fam_corrs = get_parental_correlations(probes[[proband]], 
            probes[[mom_id]], probes[[dad_id]])
        population_corrs = c(population_corrs, fam_corrs$mom, fam_corrs$dad)
    }
    
    # get the correlation between the values for the proband, and the values
    # for the different sample sets
    proband_values = probes[[person_id]]
    
    # ignore samples we cannot we get correlations from (perhaps should examine
    # the z scores alone)
    if (nrow(probes) == 1 | length(unique(proband_values)) == 1) {return}
    
    # # plot the clusters
    plot_cluster(population_z_scores, population_corrs, maternal_z_score, 
        paternal_z_score, maternal_corr, paternal_corr, row)
}


plot_cluster <- function(z_scores, correlations, maternal_z_score, paternal_z_score, maternal_corr, paternal_corr, cnv_row) {
    # plot the Z scores and correlations for the population, so that we can also
    # plot the parental data points, in order to quickly assess how far away
    # from the normal population that parental data points are.
    # 
    # Args:
    #     z_scores: Z scores of L2R for the parental population
    #     correlations: correlations of parental L2R to child L2R for the 
    #         parental population
    #     maternal_z_score: Z score for the mothers L2R
    #     paternal_z_score: Z score for the fathers L2R
    #     maternal_corr: correlations of mother's L2R to proband's L2R
    #     paternal_corr: correlations of father's L2R to proband's L2R
    #     cnv_row: row of data frame for CNV containing the proband ID, chrom, 
    #         start, stop, manually reviewed call, which we use in the title of 
    #         the plot.
    
    # combine the parental and population Z scores and correlations, so that we 
    # can construct axes that include all data points.
    all_z_scores = c(z_scores, maternal_z_score, paternal_z_score)
    all_z_scores = all_z_scores[!is.na(all_z_scores)]
    all_corrs = c(correlations, maternal_corr, paternal_corr)
    all_corrs = all_corrs[!is.na(all_corrs)]
    
    # set the title of the plot (useful for when we include multiple plots in 
    # a single file)
    title = paste(row$person_id, " chr", row$CHROM, ":", row$POS, "-", 
        row$INFO.END, " ", row$Inheritance, sep = "")
    
    # plot the population Z scores and correlations
    plot(z_scores, correlations, xlab = "L2R Z scores", 
        ylab = "correlation", xlim = c(min(all_z_scores), max(all_z_scores)), 
        ylim = c(min(all_corrs), max(all_corrs)), main = title, cex.main = 0.8)
    
    # add the parental data points, so we can see how far away from the main 
    # cluster they lie
    points(maternal_z_score, maternal_corr, col = "blue", pch = 16)
    points(paternal_z_score, paternal_corr, col = "green", pch = 16)
    legend("bottomleft", pch = 16, legend = c("mom", "dad"), col = c("blue", "green"), cex = 1)
}

