# classifies a set of CVNs that have been manually assigned inhertance states,
# so we can compare the automated predictions vs the manual classifications.

library(cifer)
library(gdata)
library(Cairo)

PLOT_GRAPHS = TRUE
REVIEWED_CNV_PATH = "/nfs/users/nfs_j/jm33/apps/cifer/control_data/exome_only_denovo_cnvs.xlsx"
DATAFREEZE_DIR="/nfs/ddd0/Data/datafreeze/1133trios_20131218"

#' open a dataset containing CNV regions by individual, with reviewed calls
#' 
#' @param path path to excel file listing manually reviewed CNV inheritances
#' 
#' @return dataframe containing CNVs with confident inheritance classifications
get_reviewed_cnvs <- function(path) {
    
    # open the file, and select the columns of interest
    cnvs = read.xls(path, stringsAsFactors = FALSE)
    cnvs = subset(cnvs, select = c(person_id, CHROM, POS, INFO.END, Inheritance))
    
    # remove the unsure calls, and unreviewed CNVs
    cnvs = cnvs[!grepl("\\?", cnvs$Inheritance), ]
    cnvs = cnvs[cnvs$Inheritance != "", ]
    
    return(cnvs)
}

classify_reviewed_cnvs <- function() {
    # classify a set of CNVs with manually reviwed inheritance states
    # 
    
    ddd = get_ddd_individuals(DATAFREEZE_DIR)
    cnvs = get_reviewed_cnvs(REVIEWED_CNV_PATH)
    cnvs$predicted_inheritance = NA
    cnvs$mom_value = NA
    cnvs$dad_value = NA
    cnvs$proband_value = NA
    cnvs$proband_z_score = NA
    
    if (PLOT_GRAPHS) {
        plot_path = paste("../cnv_inh.adm3_vs_correlations.pdf", sep = "")
        Cairo(file = plot_path, type = "pdf", height = 30, width = 22, units = "cm")
        par(mfrow = c(4, 3))
    }
    
    for (row_num in 1:nrow(cnvs)) {
        
        cnv = cnvs[row_num, ]
        
        # define the parameters of the CNV
        proband_id = cnv$person_id
        chrom = cnv$CHROM
        start = cnv$POS
        stop = cnv$INFO.END
        paternal_id = unique(ddd[ddd$individual_id == proband_id, ]$dad_id)
        maternal_id = unique(ddd[ddd$individual_id == proband_id, ]$mum_id)
        
        # run the CNV inheritance classification
        print(c(proband_id, chrom, start))
        inh = classify_exome_cnv(proband_id, maternal_id, paternal_id, chrom, start, stop, cnv)
        cnvs[row_num, ]$predicted_inheritance = inh$inheritance
        cnvs[row_num, ]$mom_value = inh$family$mom_value
        cnvs[row_num, ]$dad_value = inh$family$dad_value
        cnvs[row_num, ]$proband_value = inh$family$proband_value
        cnvs[row_num, ]$proband_z_score = inh$family$proband_z_score
        
    }
    
    if (PLOT_GRAPHS) {
        dev.off()
    }
    
    write.table(cnvs, file = "../cnv_inh.predictions.txt", quote = FALSE, sep = "\t", row.names = FALSE)
}

classify_reviewed_cnvs()

