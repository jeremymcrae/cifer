# classifies a set of CVNs that have been assigned inhertance states by VICAR,
# so we can compare the automated predictions vs the VICAR classifications.

library(cifer)

num = commandArgs(trailingOnly = TRUE)

PLOT_GRAPHS = FALSE
DATAFREEZE_DIR="/nfs/ddd0/Data/datafreeze/1133trios_20131218"
POSITIVE_CONTROL_CNV_PATH = paste("/nfs/users/nfs_j/jm33/apps/cifer/data/control_CNVs_inheritance_jeremy.txt", num, ".txt", sep = "")

#' open a dataset containing CNV regions by individual, with reviewed calls
#' 
#' @param path path to excel file listing manually reviewed CNV inheritances
#' @return dataframe containing CNVs with confident inheritance classifications
get_cnvs <- function(path) {
   
    # open the file, and select the columns of interest
    cnvs = read.table(path, stringsAsFactors = FALSE, header = TRUE)
    cnvs$chr = as.character(cnvs$chr)
    cnvs = subset(cnvs, select = c(proband_stable_id, chr, chr_start, chr_end, inherit_type))
    
    return(cnvs)
}

classify_control_cnvs <- function() {
    # classify a set of CNVs with VICAR inheritance states
    # 
    
    ddd = get_ddd_individuals(DATAFREEZE_DIR)
    cnvs = get_cnvs(POSITIVE_CONTROL_CNV_PATH)
    cnvs$predicted_inheritance = NA
    cnvs$mom_value = NA
    cnvs$dad_value = NA
    
    if (PLOT_GRAPHS) {
        plot_path = paste("../cnv_inh.adm3_vs_correlations.pdf", sep = "")
        Cairo(file = plot_path, type = "pdf", height = 30, width = 22, units = "cm")
        par(mfrow = c(4, 3))
    }
    
    for (row_num in 1:nrow(cnvs)) {
        
        row = cnvs[row_num, ]
        
        # define the parameters of the CNV
        proband_id = row$proband_stable_id
        chrom = row$chr
        start = row$chr_start
        stop = row$chr_end
        
        paternal_id = unique(ddd[ddd$individual_id == proband_id, ]$dad_id)
        maternal_id = unique(ddd[ddd$individual_id == proband_id, ]$mum_id)
        
        # divert probands who don't fit into the standard population
        if (!(proband_id %in% ddd$individual_id)) { 
            cnvs[row_num, ]$predicted_inheritance = "proband_not_in_datafreeze"
            next
        }
        
        # run the CNV inheritance classification
        print(c(proband_id, chrom, start))
        inh = classify_exome_cnv(proband_id, maternal_id, paternal_id, chrom, start, stop)
        cnvs[row_num, ]$predicted_inheritance = inh$inheritance
        cnvs[row_num, ]$mom_value = inh$mom_value
        cnvs[row_num, ]$dad_value = inh$dad_value
    }
    
    if (PLOT_GRAPHS) {
        dev.off()
    }
    
    write.table(cnvs, file = paste("../cnv_inh.control_predictions", num, ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
}

classify_control_cnvs()

