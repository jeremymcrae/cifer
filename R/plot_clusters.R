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

get_population_correlation <- function(ddd, probes, parent_ids) {
    # get the correlation between each proband and their parents L2R values
    # 
    # Args:
    #     ddd: dataframe containing all the DDD sample information.
    #     probes: L2R or ADM3 scores for all individuals for probes within the CNV
    #     parent_ids: vector of parental IDs, except for the trio under
    #         ainvestigation, sorted as for the vector of parental Z scores. 
    # 
    # Returns:
    #     numeric vector of correlation scores, (with sample IDs as names), 
    #     sorted according to the parent_ids vector.
    
    correlations = rep(NA, length(parent_ids))
    names(correlations) = parent_ids
    for (child in proband_ids) {
        proband_row = ddd[ddd$individual_id == child, ]
        dad_id = unique(proband_row$dad_id)
        mom_id = unique(proband_row$mum_id)
        
        # ignore parents that are not in the required set
        if (!(dad_id %in% parent_ids)) { next } 
        
        corrs = get_parental_correlations(probes[[child]], probes[[mom_id]], probes[[dad_id]])
        correlations[[dad_id]] = corrs$dad
        correlations[[mom_id]] = corrs$mom
    }
    
    return(correlations)
}

include_graphs <- function(ddd, probes, z_scores, proband_id, maternal_id, paternal_id, cnv) {
    # set up the data for plotting a graph, then call the plotting function
    # 
    # Initially we expected that correlations between the proband's scores with
    # parental scores might be useful, which is why this function calculates the
    # correlation, before plotting.
    # 
    # Args:
    #     ddd: dataframe containing all the DDD sample information.
    #     probes: L2R or ADM3 scores for all individuals for probes within the CNV
    #     z_scores: list of Z scores for unrelated parents, and Z scores for 
    #         members of the current trio
    #     proband_id: DDD ID of the proband
    #     maternal_id: DDD ID of the proband's mother
    #     paternal_id: DDD ID of the proband's father
    #     cnv: one row dataframe that contains the CNV's details (used to 
    #         identify the plot)
    # 
    # Returns:
    #     nothing
    
    # ignore samples we cannot we get correlations from 
    if (nrow(probes) == 1 | length(unique(probes[[proband_id]])) == 1) { return("no plot") }
    
    # get the correlation values for the proband to their parents
    trio_corr = get_parental_correlations(probes[[proband_id]], 
        probes[[maternal_id]], probes[[paternal_id]])
    trio_corr$proband = 1
    
    # get the proband IDs, except for the proband currently being analysed
    proband_ids = unique(ddd[ddd$is_proband, ]$individual_id)
    proband_ids = proband_ids[proband_ids != proband_id]
    
    # set up a vector of parental IDs, so that we can match the Z score and 
    # correlation order
    parent_ids = names(z_scores$population)
    
    # get the correlation between each proband and their parents L2R values
    correlations = get_population_correlation(ddd, probes, parent_ids)
    
    # # plot the clusters
    plot_cluster(z_scores, correlations, trio_corr, cnv)
}

plot_cluster <- function(z_scores, correlations, trio_corr, cnv) {
    # plot the Z scores and correlations for the population, so that we can also
    # plot the parental data points, in order to quickly assess how far away
    # from the normal population that parental data points are.
    # 
    # Args:
    #     z_scores: Z scores of L2R for the parental population
    #     correlations: correlations of parental L2R to child L2R for the 
    #         parental population
    #     trio_corr: correlations of mother's L2R to proband's L2R
    #     row: row of data frame for CNV containing the proband ID, chrom, 
    #         start, stop, manually reviewed call, which we use in the title of 
    #         the plot.
    
    # combine the parental and population Z scores and correlations, so that we 
    # can construct axes that include all data points.
    z = c(z_scores$population, z_scores$mom, z_scores$dad, z_scores$proband)
    corrs = c(correlations, trio_corr$mom, trio_corr$mom, trio_corr$proband)
    all = data.frame(z_scores = z, corrs)
    all = all[complete.cases(all), ]
    
    # set the title of the plot (useful for when we include multiple plots in 
    # a single file)
    title = paste(cnv$person_id, " chr", cnv$CHROM, ":", cnv$POS, "-", 
        cnv$INFO.END, " ", cnv$Inheritance, sep = "")
    
    # plot the population Z scores and correlations
    plot(all$z_scores, all$corrs, xlab = "L2R Z scores", 
        ylab = "correlation", xlim = c(min(all$z_scores), max(all$z_scores)), 
        ylim = c(min(all$corrs), max(all$corrs)), main = title, cex.main = 0.8)
    
    # add the parental data points, so we can see how far away from the main 
    # cluster they lie
    points(z_scores$mom, trio_corr$mom, col = "blue", pch = 16)
    points(z_scores$dad, trio_corr$dad, col = "green", pch = 16)
    points(z_scores$proband, trio_corr$proband, col = "red", pch = 16)
    legend("bottomleft", pch = 16, legend = c("mom", "dad", "proband"), 
        col = c("blue", "green", "red"), cex = 1)
}

