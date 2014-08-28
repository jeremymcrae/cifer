# classifies a set of CVNs that have been manually assigned inhertance states,
# so we can compare the automated predictions vs the manual classifications.

library(gdata)

source("plot_clusters.R")
source("call_exome_cnvs.R")

PLOT_GRAPHS = FALSE
REVIEWED_CNV_PATH = "/nfs/users/nfs_j/jm33/apps/exome_cnv_inheritance/data/exome_only_denovo_cnvs.xlsx"

get_reviewed_cnvs <- function(path) {
    # open a dataset containing CNV regions by individual, with reviewed calls
    # 
    # Args:
    #     path: path to excel file listing manually reviewed CNV inheritances
    # 
    # Returns:
    #     dataframe containing CNVs with confident inheritance classifications
    
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
    
    if (PLOT_GRAPHS) {
        plot_path = paste("../cnv_inh.adm3_vs_correlations.pdf", sep = "")
        Cairo(file = plot_path, type = "pdf", height = 30, width = 22, units = "cm")
        par(mfrow = c(4, 3))
    }
    
    for (row_num in 1:nrow(cnvs)) {
        
        row = cnvs[row_num, ]
        
        # define the parameters of the CNV
        proband_id = row$person_id
        chrom = row$CHROM
        start = row$POS
        stop = row$INFO.END
        paternal_id = unique(ddd[ddd$individual_id == proband_id, ]$dad_id)
        maternal_id = unique(ddd[ddd$individual_id == proband_id, ]$mum_id)
        
        # run the CNV inheritance classification
        print(c(proband_id, chrom, start))
        inh = classify_exome_cnv(proband_id, maternal_id, paternal_id, chrom, start, stop)
        cnvs[row_num, ]$predicted_inheritance = inh
    }
    
    if (PLOT_GRAPHS) {
        dev.off()
    }
    
    write.table(cnvs, file = "../cnv_inh.predictions.txt", quote = FALSE, sep = "\t", row.names = FALSE)
}

classify_reviewed_cnvs()

