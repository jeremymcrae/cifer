# check the performance of predicting inheritance, by subsampling from
# large, probe-rich manually classified CNVs.
# 

library(cifer)

REVIEWED_CNV_PATH = "/nfs/users/nfs_j/jm33/apps/cifer/data/exome_only_denovo_cnvs.xlsx"
DATAFREEZE_DIR="/nfs/ddd0/Data/datafreeze/1133trios_20131218"

#' test the inheritance classification of subsets of probes in a CNV
#' 
#' We subset the probes, in order to check the performance of the classifier
#' with smaller numbers of probes. Use a sliding window of n probes, where n
#' varies between 1 and 10.
#' 
#' @param cnv_predictions blank data frame to put classifcations into
#' @param cnvs data frame containing CNV information
#' @param row_num row of the cnvs dataframe to examine
#' @param ddd dataframe containing all the DDD sample information.
#' 
#' @return data frame containing CNV information (with number of probes), and
#'     the inheritance classification.
test_cnv_subset <- function(cnv_predictions, cnvs, row_num, ddd) {
    
    row = cnvs[row_num, ]
        
    paternal_id = unique(ddd[ddd$individual_id == row$person_id, ]$dad_id)
    maternal_id = unique(ddd[ddd$individual_id == row$person_id, ]$mum_id)
    
    print(row$person_id)
    probes = get_probes_data(ddd, row$CHROM, row$POS, row$INFO.END)
    
    # for the probes that exist for a CNV, take subsets of those probes, so
    # that we can examine the performance under different stringencies
    for (num_probes in 1:10) {
        print(num_probes)
        start_probe = 1
        end_probe = start_probe + num_probes - 1
        
        # slide a window of probes across the set of full probes, and
        # predict the inheritance from that
        while (end_probe <= nrow(probes)) {
            probe_subset = probes[start_probe:end_probe, ]
            predictions = process_cnv_call(ddd, probe_subset, row$person_id, maternal_id, paternal_id)
            
            temp_row = data.frame("person_id" = row$person_id, "CHROM" = row$CHROM, 
                "POS" = row$POS, "INFO.END" = row$INFO.END, probe_count = num_probes, 
                "Inheritance" = row$Inheritance, "predicted_inheritance" = predictions$inheritance, 
                "mom_value" = predictions$mom_value, "dad_value" = predictions$dad_value)
            
            cnv_predictions = rbind(cnv_predictions, temp_row)
            
            start_probe = start_probe + 1
            end_probe = end_probe + 1
        }
    }
    
    return(cnv_predictions)
}

#' plot the number of correct predictions by number of probes selected, for 
#' each of the inheritance categories
#' 
#' @param cnv_predictions data frame of CNVs, including coordinates, number of
#'     probes, inheritance classification (manually reviewed and predicted)
plot_correct_ratios <- function(cnv_predictions) {
    
    # count the number of predictions that were correct, or incorrect, for each
    # inheritance type
    counts = cast(cnv_predictions, Inheritance + consistent ~ probe_count, value = "person_id", fun = length)
    
    ratios = data.frame(probe_count = 1:10)
    for (inh in c("Pat_INH", "Mat_INH", "DNM")){
        temp = melt(counts)
        
        false_n = temp[temp$Inheritance == inh & temp$consistent == FALSE, ]
        true_n = temp[temp$Inheritance == inh & temp$consistent == TRUE, ]
        
        # add the proportions to a data frame
        inh_ratio = true_n/(true_n + false_n)
        inh_ratio = data.frame(inh_ratio$value)
        names(inh_ratio) = c(inh)
        ratios = cbind(ratios, inh = inh_ratio)
    }
    
    # plot the ratios at each tested number of probes for the different
    # inheritance types
    par(cex.lab = 2, cex.axis = 2, cex = 2)
    Cairo(file = "../cnv_inh.large_subset_ratios.pdf", type = "pdf", width = 15, height = 15, units = "cm")
    plot(ratios$probe_count, ratios$DNM, type = "b", ylim = c(0,1), las = 1, 
        xlab = "Number of probes", ylab = "Proportion with correct inheritance",
        pch = 15, col = "black")
    points(ratios$probe_count, ratios$Pat_INH, pch = 16, col = "darkgreen", type = "b")
    points(ratios$probe_count, ratios$Mat_INH, pch = 17, col = "blue", type = "b")
    legend("bottomright", legend = c("de novo", "paternal", "maternal"), 
        pch = c(15, 16, 17), col = c("black", "darkgreen", "blue"), 
        lwd = 1, bty = "n")
    
    dev.off()
}

main <- function() {
    # load the datafiles
    ddd = get_ddd_individuals(DATAFREEZE_DIR)
    ddd = find_convex_files(ddd)
    
    allowed_inh = c("Pat_INH", "Mat_INH", "Biparental", "DNM")
    
    # open the list of manually reviewed CNVs, subset down to the ones with 
    # permitted inheritances, and with sufficient probes and length
    cnvs = read.xls(REVIEWED_CNV_PATH, stringsAsFactors = FALSE)
    cnvs = cnvs[cnvs$Inheritance %in% allowed_inh, ]
    cnvs = cnvs[cnvs$INFO.NUMBERPROBESCONVEX > 10 & cnvs$INFO.SVLEN > 10000, ]
    
    # order the CNVs by how many probes each CNV has
    cnvs = cnvs[order(cnvs$INFO.NUMBERPROBESCONVEX, decreasing = TRUE), ]
    
    cnv_predictions = data.frame("person_id" = character(0), "CHROM" = character(0), 
        "POS" = character(0), "INFO.END" = numeric(0), probe_count = numeric(), 
        "Inheritance" = character(0), "predicted_inheritance" = character(0), 
        "mom_value" = numeric(0), "dad_value" = numeric(0))
    
    for (row_num in 1:nrow(cnvs)) {
        cnv_predictions = test_cnv_subset(cnv_predictions, cnvs, row_num, ddd)
    }
    
    # convert the inheritance strings to those used in the CNV excel file
    cnv_predictions$fixed_inheritance = NA
    cnv_predictions$fixed_inheritance[cnv_predictions$predicted_inheritance == "paternal_inh"] = "Pat_INH"
    cnv_predictions$fixed_inheritance[cnv_predictions$predicted_inheritance == "maternal_inh"] = "Mat_INH"
    cnv_predictions$fixed_inheritance[cnv_predictions$predicted_inheritance == "biparental"] = "Biparental_INH"
    cnv_predictions$fixed_inheritance[cnv_predictions$predicted_inheritance == "de_novo"] = "DNM"
    
    # check which predictions match the original reviewed call for the full CNV
    cnv_predictions$consistent = as.character(cnv_predictions$Inheritance) == as.character(cnv_predictions$fixed_inheritance)
    
    # trim out one CNV that has consistently wrong predictions - it's unclear 
    # whether the manually reviewed status is incorrect
    cnv_predictions = cnv_predictions[cnv_predictions$person_id != "DDDP107519" | cnv_predictions$CHROM != "X" | cnv_predictions$POS != 53449333, ]
    write.table(cnv_predictions, file = "../cnv_inh.large_cnv_subset_predictions.txt", quote = FALSE, sep = "\t", row.names = FALSE)
    
    plot_correct_ratios(cnv_predictions)
}

main()
