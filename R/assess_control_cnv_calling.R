# script to examine the output from predicting CNV inheritance for the control 
# CNVs provided by Tom Fitzgerald.
## NOTE: requires R-3.1 or higher for dplyr library

if (!(as.numeric(R.Version()$major) >= 3 & as.numeric(R.Version()$minor) > 0.1)) {
    stop("need R version 3.0.1 or higher for dplyr.")
}

library(dplyr)
library(reshape)
library(ggplot2)
library(Cairo)

EXOME_DIR = "/nfs/users/nfs_j/jm33/apps/exome_cnv_inheritance"
DATA_DIR = file.path(EXOME_DIR, "data")
RESULTS_DIR = file.path(EXOME_DIR, "results")
POSITIVE_CONTROL_CNV_PATH = file.path(DATA_DIR, "control_CNVs_inheritance_jeremy.txt")

open_inheritance_calls <- function(folder) {
    # open the inheritance calls predicted from the exome data
    # 
    # We might start off with the exome-based CNV cals split across multiple 
    # files. This function loads all the separate files (if necessary), and
    # generates a file with all the CNV calls from the separate files in it, so
    # that we can load that single file later.
    # 
    # Args:
    #     folder: folder containing files with the exome inheritance calls
    # 
    # Returns:
    #     dataframe containing CNV inheritance calls from the exome-based 
    #     classifier
    
    # if we have already merged the files together, don't worry about doing this
    # again
    merged_path = file.path(folder, "cnv_inh.control_predictions.txt")
    if (file.exists(merged_path)) {
        cnvs = read.table(merged_path, header = TRUE)
        return(cnvs)
    }
    
    # otherwise combine the data from the separate files, then make a merged file
    paths = Sys.glob(file.path(folder, "cnv_inh.control_predictions*"))
    files = vector("list", length = length(paths))
    
    for (pos in 1:length(paths)) {
        files[[pos]] = read.table(paths[pos], header = TRUE)
    }
    
    cnvs = rbind_all(files)
    write.table(cnvs, file = merged_path, sep = "\t", quote = FALSE, row.names = FALSE)
    
    return(cnvs)
}

open_vicar_calls <- function(path) {
    # open the inheritance calls predicted from the array data.
    # 
    # Args:
    #     path: path to the VICAR classification data file
    # 
    # Returns:
    #     dataframe containing VICAR CNVs
    
    cnvs = read.table(path, stringsAsFactors = FALSE, header = TRUE)
    cnvs$chr = as.character(cnvs$chr)
    
    # # trim some unnessary columns
    # cnvs$proband_file_path = NULL
    # cnvs$father_file_path = NULL
    # cnvs$mother_file_path = NULL
    # cnvs$father_stable_id = NULL
    # cnvs$mother_stable_id = NULL
    # cnvs$proband_sanger_id = NULL
    
    return(cnvs)
}

combine_cnv_calls <- function(exome_calls, array_calls) {
    # combines the VICAR and exome-based inheritance calls
    # 
    # Args:
    #     exome_calls: CNVs, with inheritance predictions from the exome 
    #         inheritance classifier
    #     array_calls: CNVs with inheritance predictions from VICAR
    # 
    # Returns:
    #     a dataframe where the VICAR and exome predictions are merged, and the
    #     prediction names are standardised between classifiers
    
    cnvs = merge(array_calls, exome_calls, by = c("chr", "chr_start", 
        "chr_end", "proband_stable_id", "inherit_type"))
    
    # exclude CNVs without matching exome data
    cnvs = cnvs[cnvs$predicted_inheritance != "proband_not_in_datafreeze", ]
    cnvs = cnvs[cnvs$predicted_inheritance != "no_probe_data_for_CNV_region", ]
    cnvs = cnvs[cnvs$predicted_inheritance != "unable_to_evaluate_probes", ]
    
    # standardise the inheritance state names to those given by VICAR
    cnvs$inheritance = NA
    cnvs$inheritance[cnvs$predicted_inheritance == "de_novo"] = "deNovo"
    cnvs$inheritance[cnvs$predicted_inheritance == "paternal_inh"] = "paternal"
    cnvs$inheritance[cnvs$predicted_inheritance == "maternal_inh"] = "maternal"
    cnvs$inheritance[cnvs$predicted_inheritance == "biparental_inh"] = "biparental"
    cnvs$inheritance[cnvs$predicted_inheritance == "uncertain"] = "uncertain"
    
    # and reorganise the column names
    cnvs$predicted_inheritance = cnvs$inheritance
    cnvs$inheritance = NULL
    
    # I assume inheritedDuo means the same thing as biparental. Why the difference?
    cnvs$inherit_type[cnvs$inherit_type == "inheritedDuo"] = "biparental"
    
    cnvs$consistent = as.character(cnvs$inherit_type) == as.character(cnvs$predicted_inheritance)
    
    return(cnvs)
}

get_certain_cnvs <- function(cnvs) {
    # restrict ourselves to the set of CNVs with VICAR calls that are certain,
    # ie drop all of the "insufficient", "inconclusive" etc calls
    # 
    # Args:
    #     cnvs: dataframe of CNVs
    # 
    # Returns:
    #     dataframe of CNVs, restrcicted to ones with certain calls
    
    uncertain_vicar_calls = c("inconclusive", "noCNVfound", "inconclusiveDuo", 
        "insufficient", "noCNVDuo", "noCNVfoundProband", "noCNVprobandDuo")
    certain_cnvs = cnvs[!(cnvs$inherit_type %in% uncertain_vicar_calls), ]
    
    return(certain_cnvs)
}

bin_cnvs <- function(cnvs, column, title) {
    # check how different metrics impact matching the inheritance state
    # 
    # Args:
    #     cnvs: dataframe of CNVs
    #     column: column from CNVs data frame, to be used to bin the CNVs
    #     title: string for the title of the plot
    # 
    # Returns:
    #    a qplot, to be include in a multi-plot later on
    
    probs = seq(0.0, 1, 0.2)
    breaks = quantile(column, probs, names=FALSE)
    quantile = try(cut(column, breaks, include.lowest=TRUE), silent = TRUE)
    
    # sometimes we have non-unique breakpoints, just default to very broad
    # quantiles in this case
    if (class(quantile) == "try-error") {
        breaks = quantile(column, c(0, 0.5, 1), names=FALSE)
        quantile = cut(column, breaks, include.lowest=TRUE)
    }
    
    cnvs$quantile = quantile
    ratios = cast(cnvs, predicted_inheritance + quantile ~ consistent, value = "proband_stable_id", length)
    
    # cacluate how many calls are correct at each qwuantile, for the different 
    # call types
    ratios$match_ratio = ratios[["TRUE"]]/(ratios[["TRUE"]] + ratios[["FALSE"]])
    
    # plot the ratio of correct calls by quantile, for the different call types
    plot = qplot(quantile, match_ratio, data = ratios) 
    plot = plot + aes(color = predicted_inheritance, group = predicted_inheritance) 
    plot = plot + geom_line() + theme_bw() + ggtitle(title)
    
    return(plot)
}

vp.setup <- function(x, y){
    # create a new layout with grid
    grid.newpage()
    # define viewports and assign it to grid layout
    pushViewport(viewport(layout = grid.layout(x, y)))
}

vp.layout <- function(x, y){
    # define function to easily access layout (row, col)
    viewport(layout.pos.row=x, layout.pos.col=y)
}

analyse_performance <- function(cnvs, pdf_filename) {
    # analyse how well the exome predictions match the VICAR predictions
    # 
    # Args:
    #     cnvs: dataframe of CNVs with exome-based and VIAR inheritance states
    #     pdf_filename: filename to export plots into
    # 
    # Returns:
    #     nothing
    
    # tabulate the differences between VICAR and exome calls
    print(cast(cnvs, predicted_inheritance ~ consistent, value = "proband_stable_id", length))
    print(cast(cnvs, predicted_inheritance ~ inherit_type, value = "proband_stable_id", length))
    
    # plot the ratio of correct matches across quantiles of different CNV scores
    Cairo(file = file.path(EXOME_DIR, pdf_filename), type = "pdf", height = 25, 
        width = 30, units = "cm")
    vp.setup(2,2)
    print(bin_cnvs(cnvs, cnvs$number_exome_probes, "number of exome probes"), vp=vp.layout(1,1))
    print(bin_cnvs(cnvs, cnvs$adjw_score, "adjusted w-score"), vp=vp.layout(1,2))
    print(bin_cnvs(cnvs, cnvs$ratio_mean, "ratio mean"), vp=vp.layout(2,1))
    print(bin_cnvs(cnvs, cnvs$mad_region, "mad region"), vp=vp.layout(2,2))
    dev.off()
}

main <- function() {
    vicar_calls = open_vicar_calls(POSITIVE_CONTROL_CNV_PATH)
    exome_calls = open_inheritance_calls(RESULTS_DIR)
    
    # combine the CNV calls, and export the file, so we can analyse them later
    cnvs = combine_cnv_calls(exome_calls, vicar_calls)
    write.table(cnvs, file = file.path(EXOME_DIR, "control_CNVs_VICAR_and_exome.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    
    # get different CNV subsets
    certain_cnvs = get_certain_cnvs(cnvs)
    rare_cnvs = cnvs[cnvs$rare == 1, ]
    common_cnvs = cnvs[cnvs$rare == 0, ]
    
    # analyse how well the inheritance states match between the VICAR calls, and
    # the exome-based calls
    analyse_performance(cnvs, "all_cnvs.pdf")
    analyse_performance(certain_cnvs, "cnvs_with_certain_vicar_classifications.pdf")
    analyse_performance(rare_cnvs, "rare_cnvs.pdf")
    analyse_performance(common_cnvs, "common_cnvs.pdf")
    
}

main()
