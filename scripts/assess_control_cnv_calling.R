# script to examine the output from predicting CNV inheritance for the control 
# CNVs provided by Tom Fitzgerald.

library(dplyr)
library(reshape)
library(Cairo)
library(ggplot2)
library(grid)

EXOME_DIR = "/nfs/users/nfs_j/jm33/apps/cifer"
DATA_DIR = file.path(EXOME_DIR, "control_data")
RESULTS_DIR = file.path(EXOME_DIR, "results")
POSITIVE_CONTROL_CNV_PATH = file.path(DATA_DIR, "control_CNVs_inheritance_jeremy.txt")
DATAFREEZE_DIR="/nfs/ddd0/Data/datafreeze/1133trios_20131218"

# seed the random number generator, so we get repeatable selections of CNVs
set.seed(1043290)

#' open the inheritance calls predicted from the exome data
#' 
#' We might start off with the exome-based CNV cals split across multiple 
#' files. This function loads all the separate files (if necessary), and
#' generates a file with all the CNV calls from the separate files in it, so
#' that we can load that single file later.
#' 
#' @param folder folder containing files with the exome inheritance calls
#' @return dataframe containing CNV inheritance calls from the exome-based 
#'     classifier
open_inheritance_calls <- function(folder) {
    
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

#' open the inheritance calls predicted from the array data.
#' 
#' @param path path to the VICAR classification data file
#' 
#' @return dataframe containing VICAR CNVs
open_vicar_calls <- function(path) {
    
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

#' combines the VICAR and exome-based inheritance calls
#' 
#' @param exome_calls CNVs, with inheritance predictions from the exome 
#'     inheritance classifier
#' @param array_calls CNVs with inheritance predictions from VICAR
#' 
#' @return a dataframe where the VICAR and exome predictions are merged, and the
#'     prediction names are standardised between classifiers
combine_cnv_calls <- function(exome_calls, array_calls) {
    
    cnvs = merge(array_calls, exome_calls, by = c("chr", "chr_start", 
        "chr_end", "proband_stable_id", "inherit_type"))
    
    # exclude cnvs with low MAD ratios, since they are more likely to be false
    # positives, and also filter to the rare CNVs (to match some analysis that
    # Tom ran)
    cnvs$mad_ratio = abs(cnvs$ratio_mean/cnvs$mad_region)
    cnvs = cnvs[cnvs$mad_ratio > 10, ]
    cnvs = cnvs[cnvs$rare == 1, ]
    
    # exclude CNVs without matching exome data
    cnvs = cnvs[cnvs$predicted_inheritance != "proband_not_in_datafreeze", ]
    cnvs = cnvs[cnvs$predicted_inheritance != "no_probe_data_for_CNV_region", ]
    cnvs = cnvs[cnvs$predicted_inheritance != "unable_to_evaluate_probes", ]
    
    # standardise the inheritance state names to those given by VICAR
    cnvs$inheritance = NA
    cnvs$inheritance[cnvs$predicted_inheritance == "deNovo"] = "deNovo"
    cnvs$inheritance[cnvs$predicted_inheritance == "de_novo"] = "deNovo"
    cnvs$inheritance[cnvs$predicted_inheritance == "not_inherited"] = "deNovo"
    cnvs$inheritance[cnvs$predicted_inheritance == "paternal_inh"] = "paternal"
    cnvs$inheritance[cnvs$predicted_inheritance == "maternal_inh"] = "maternal"
    cnvs$inheritance[cnvs$predicted_inheritance == "biparental_inh"] = "biparental"
    cnvs$inheritance[cnvs$predicted_inheritance == "uncertain"] = "uncertain"
    cnvs$inheritance[cnvs$predicted_inheritance == "false_positive"] = "false_positive"
    
    # and reorganise the column names
    cnvs$predicted_inheritance = cnvs$inheritance
    cnvs$inheritance = NULL
    
    # I assume inheritedDuo means the same thing as biparental. Why the difference?
    cnvs$inherit_type[cnvs$inherit_type == "inheritedDuo"] = "biparental"
    cnvs$consistent = as.character(cnvs$inherit_type) == as.character(cnvs$predicted_inheritance)
    
    # drop probands with excessive number of CNV calls
    cnv_counts = table(cnvs$proband_stable_id)
    keep_probands = names(cnv_counts[cnv_counts <= 10])
    cnvs = cnvs[cnvs$proband_stable_id %in% keep_probands, ]
    
    return(cnvs)
}

#' restrict ourselves to the set of CNVs with VICAR calls that are certain,
#' ie drop all of the "insufficient", "inconclusive" etc calls
#' 
#' @param cnvs dataframe of CNVs
#' 
#' @return dataframe of CNVs, restrcicted to ones with certain calls
get_certain_cnvs <- function(cnvs) {
    
    uncertain_vicar_calls = c("inconclusive", "noCNVfound", "inconclusiveDuo", 
        "insufficient", "noCNVDuo", "noCNVfoundProband", "noCNVprobandDuo")
    certain_cnvs = cnvs[!(cnvs$inherit_type %in% uncertain_vicar_calls), ]
    
    return(certain_cnvs)
}

#' check how different metrics impact matching the inheritance state
#' 
#' @param cnvs dataframe of CNVs
#' @param column column from CNVs data frame, to be used to bin the CNVs
#' @param title string for the title of the plot
#' 
#' @return a qplot, to be include in a multi-plot later on
bin_cnvs <- function(cnvs, column, title) {
    
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

#' analyse how well the exome predictions match the VICAR predictions
#' 
#' @param cnvs dataframe of CNVs with exome-based and VIAR inheritance states
#' @param pdf_filename filename to export plots into
analyse_performance <- function(cnvs, pdf_filename) {
    
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

#' pick CNVs with mismatching inheritance predictions, for manual review
#' 
#' @param cnvs dataframe of high-quality CNVs with cifer and VICAR inheritance
#'     predictions
#'
#' @return assess_cnvs dataframe of CNVs to be manually assessed
select_cnvs_to_assess <- function(cnvs) {
    
    # select only the CNVs that do not match, so we can examine a sample of them
    cnvs = cnvs[cnvs$consistent == FALSE, ]
    
    # split the cnvs into the different mismatch categories, so, for example, we
    # could look at the CNVs where cifer predicted paternal, and VICAR 
    # predicted inconclusive
    cnvs$mismatch_type = paste(cnvs$inherit_type, cnvs$predicted_inheritance, sep = "-")
    
    # construct a url from which we will be able to examine the array CGH and
    # exome data used to make the CNV call in the members of the trio
    cnvs$url = paste("http://ddd-vm3.internal.sanger.ac.uk:9090/ddd.html?", cnvs$proband_sanger_id, ";", cnvs$chr, ";", cnvs$chr_start, ";", cnvs$chr_end, sep = "")
    
    cnvs = subset(cnvs, select = c(mismatch_type, proband_stable_id, chr, chr_start, chr_end, url))
    
    assess_cnvs = cnvs[0, ]
    # for each mistmatch category, get a randome sample of up to ten CNV calls
    for (mismatch in unique(cnvs$mismatch_type)) {
        mismatches = cnvs[cnvs$mismatch_type == mismatch, ]
        
        selection = mismatches[sample(nrow(mismatches), 10, replace = TRUE), ]
        selection = unique(selection)
        assess_cnvs = rbind(assess_cnvs, selection)
    }
    
    return(assess_cnvs)
}

main <- function() {
    vicar_calls = open_vicar_calls(POSITIVE_CONTROL_CNV_PATH)
    exome_calls = open_inheritance_calls(RESULTS_DIR)
    
    # combine the CNV calls, and export the file, so we can analyse them later
    cnvs = combine_cnv_calls(exome_calls, vicar_calls)
    write.table(cnvs, file = file.path(EXOME_DIR, "control_CNVs_VICAR_and_exome.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    
    # get different CNV subsets
    certain_cnvs = get_certain_cnvs(cnvs)
    
    # analyse how well the inheritance states match between the VICAR calls, and
    # the exome-based calls
    analyse_performance(cnvs, "all_cnvs.pdf")
    analyse_performance(certain_cnvs, "cnvs_with_certain_vicar_classifications.pdf")
    
    cnvs_for_assessment = select_cnvs_to_assess(cnvs)
    write.table(cnvs_for_assessment, file = file.path(EXOME_DIR, "cnv_calls_to_review.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}

main()
