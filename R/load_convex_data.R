# functions to load probe data for the DDD population
# 

#' find the convex data files for each DDD individual
#' 
#' @param ddd data frame of participant informatyion, including sample IDs
#' @param convex_dir path to folder containing individual L2R datafiles
#' @export
#' @return dataframe with additional column of file paths
find_convex_files <- function(ddd, convex_dir) {
    
    # convex_dir = "/nfs/ddd2/ddd_convex_data"
    
    ddd$convex_path = NA
    individuals = unique(ddd$individual_id)
    
    for (person_id in individuals) {
        sanger_ids = ddd$sanger_id[ddd$individual_id == person_id]
        
        # many individuals have multiple sanger IDs, so we need to find which
        # one of those has convex L2R data associated with it. Currently only
        # one of the sanger IDs for an individual will have CONVEX L2R data, 
        # or where more than one of the sanger IDs have data files, both files
        # contain the same data.
        for (sanger_id in sanger_ids) {
            basename = paste("ProbeRD_", sanger_id, "_LR2_GAM.dat",  sep="")
            # split_path = paste(substring(person_id, seq(1, nchar(person_id), 2), seq(2, nchar(person_id), 2)), collapse="/")
            # path = file.path(convex_dir, split_path, person_id, "convex", sanger_id, basename)
            path = file.path(convex_dir, basename)
            
            # if the file exists, then we can use it
            if (file.exists(path)) {
                ddd$convex_path[ddd$sanger_id == sanger_id] = path
            }
        }
    }
    
    return(ddd)
}

#' convert a chrom string to a number (make sex chroms 23 and beyond)
#' 
#' @param chrom chromosome ID as string eg  "1", "2" ... "X"
#' @export
#' @return chromosome as number, where X and Y are 23 and 24
chrom_num <- function(chrom) {
    
    # construct a dictionary-like object
    chrom_dict = vector(mode = "list", length=24)
    chrom_dict[1:24] = 1:24
    names(chrom_dict) = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
        "X", "Y")
    
    chrom_num = chrom_dict[[chrom]]
    
    return(chrom_num)
}

#' find where the required data lies within the file
#' 
#' Seek in a file to the section containing the data that we want. We start
#' with a reasonable guess about where in the file the data is located, but 
#' we seek up or downstream of that (depending on whether the data at that 
#' position is before or after the data that we want), then bisect to pin 
#' down where to start loading data from.
#' 
#' @param path path to the datafile
#' @param start_pos reasonable guess about where the correct data is located in 
#'     the file
#' @param chrom chrom that we are searching for
#' @param start start position of the CNV
#' @param filesize size of L2R file in bytes
#' @export
#' @return the approximate byte offset in the file (to within a few hundred
#'    bytes) where the required data starts
find_file_pos <- function(path, start_pos, chrom, start, file_size) {
    
    # seek to the correct part of the file
    size = file.info(path)$size
    file_con = file(path, open = "r")
    prev_pos = seek(file_con, start_pos, origin = "start")
    
    read = scan(file_con, what = character(), skip = 1, sep = "\t", nmax = 2, quiet = TRUE)
    temp_chrom = read[1]
    temp_start = as.numeric(read[2])
    
    # get an initial end point to bisect within, make the search window 
    # proportional to the filesize, so that we check across a larger range with
    # larger files
    offset = -file_size/200
    if (temp_start <= start | (temp_chrom != chrom)) {offset = abs(offset)}
    
    # get the upper and lower bounds to bisect within
    lower = min(start_pos, start_pos + offset)
    upper = max(start_pos, start_pos + offset)
    
    # make sure we don't use negative seek positions, or go beyond the file size
    lower = max(lower, 0)
    upper = min(upper, size)
    
    while (upper - lower > 100) {
        mid_point = lower + (upper - lower)/2
        prev_pos = seek(file_con, mid_point, origin = "start")
        
        read = scan(file_con, what = character(), skip = 1, sep = "\t", nmax = 2, quiet = TRUE)
        temp_chrom = read[1]
        temp_start = as.numeric(read[2])
        
        # if we have gone to the next chromosome, adjust the mid point until 
        # we hit the correct chrom
        if (chrom_num(temp_chrom) < chrom_num(chrom)) {lower = mid_point; next}
        if (chrom_num(temp_chrom) > chrom_num(chrom)) {upper = mid_point; next}
        
        # narrow the boundaries to half the previous size
        if (temp_start > start) {upper = mid_point}
        if (temp_start <= start) {lower = mid_point}
    }
    
    close(file_con)
    
    # make sure we don't use negative seek positions
    seek_position = max(lower - 100, 0)
    
    return(seek_position)
}

#' quickly open a section of a file
#' 
#' @param path path to L2R dataset
#' @param file_pos position in the file to jump ahead to
#' @param lines_n number of lines to read in
#' @param chrom chromosome string for current CNV
#' @param start CNV start nucleotide position as integer
#' @param filesize size of L2R file in bytes
#' @param l2r_type string indicating the type of L2R scores to load, "adm3" or 
#'     "l2r"
#' @export
#' @return dataframe of probe information, with one row per probe
extract_lines_from_file <- function(path, file_pos, lines_n, chrom, start, filesize, l2r_type) {
    
    pos = find_file_pos(path, file_pos, chrom, start, filesize)
    
    # seek to the correct part of the file
    file_con = file(path, open = "r")
    prev_pos = seek(file_con, pos, origin = "start")
    
    # figure out which columns of the line to load, depending on whether we are
    # using L2R or ADM3 scores
    line_classes = c("character", "numeric", "numeric", "NULL", "NULL", "NULL",
        "NULL", "NULL")
    if (l2r_type == "adm3") {
        line_classes[8] = "numeric"
    } else if (l2r_type == "l2r") {
        line_classes[4] = "numeric"
    }
    
    # read sufficient lines from the file, to contain the correct positions. 
    # Skip the first line, since it is probably not a full line.
    probes = read.table(file_con, skip = 1, nrows = lines_n + 3, header = FALSE,
        sep = "\t", stringsAsFactors = FALSE, colClasses = line_classes)
    close(file_con)
    
    names(probes) = c("chrom", "start", "stop", "score")
    
    return(probes)
}

#' open convex data from an aggregated convex dataset
#' 
#' @param ddd dataframe listing sample IDs, and file paths
#' @param chrom chromosome of interest as string (eg '1', '2', 'X')
#' @param start chromosome position that the CNV region starts at
#' @param stop chromosome position that the CNV region ends at
#' @param l2r_type string indicating the type of L2R scores to load, "adm3" or 
#'     "l2r"
#' @param l2r_dir path to folder containing the aggregated L2R scores
#' @export
#' @return a dataframe containing all the data for the different probes in 
#'     the CNV region, for all the particpants in the ddd dataframe
open_probes_data_fast <- function(ddd, chrom, start, stop, l2r_type, l2r_dir) {
    
    regions_path = file.path(l2r_dir, paste(l2r_type, "_index.chr", chrom, 
        ".txt", sep = ""))
    index_data = read.table(regions_path, header = TRUE, sep = "\t",
        comment.char="", colClasses=c("character", "numeric", "numeric"))
    
    path = file.path(l2r_dir, paste(l2r_type, "_dataset.chr", chrom, ".txt", sep = ""))
    header = read.table(path, header = TRUE, sep = "\t", nrows = 1)
    
    # figure out which lines of the index file correspond to the correct region
    index = which(index_data$chrom == chrom & index_data$start <= stop & 
        index_data$stop >= start)
    
    # raise error, if there isn't any probe data for the CNV region
    if (length(index) == 0) { stop("no probe data for CNV") }
    
    start_index = min(index)
    lines_n = (max(index) - start_index) + 1
    
    # figure out the correct place to start reading from the data file
    file_size = file.info(path)$size
    file_pos = ((start_index/nrow(index_data)) * file_size) - ncol(header) * 150
    correct_pos = find_file_pos(path, file_pos, chrom, start, file_size)
    
    file_con = file(path, open = "r")
    prev_pos = seek(file_con, correct_pos, origin = "start")
    
    # read sufficient lines from the file, to contain the correct positions.
    probes = read.table(file_con, skip = 1, nrows = lines_n + 2, header = FALSE,
        sep = "\t", stringsAsFactors = FALSE)
    close(file_con)
    names(probes) = names(header)
    
    # select the correct rows, and drop the first three columns with the chrom,
    # start and stop
    probes = probes[probes$chrom == chrom & probes$start <= stop & 
        probes$stop >= start, ]
    probes = subset(probes, select = -c(chrom, start, stop))
    
    return(probes)
}

#' parse all the convex files for the individuals, to get the L2R values 
#' for every individual
#' 
#' @param ddd dataframe listing sample IDs, and file paths
#' @param chrom chromosome of interest as string (eg '1', '2', 'X')
#' @param start chromosome position that the CNV region starts at
#' @param stop chromosome position that the CNV region ends at
#' @param l2r_type string indicating the type of L2R scores to load, "adm3" or 
#'     "l2r"
#' @export
#' @return a dataframe containing all the L2R data for the different probes in 
#'     the CNV region, for all the particpants in the ddd dataframe
open_probes_data_slow <- function(ddd, chrom, start, stop, l2r_type) {
    
    convex = ddd[!is.na(ddd$convex_path), ]
    header = c("chrom", "start", "stop", "l2r", "V5", "V6", "V7", "adm3")
    
    regions_path = unique(convex$convex_path)[1]
    index_data = read.table(regions_path, header = TRUE, sep = "\t",
        comment.char="", colClasses=c("character", "numeric", "numeric", NULL,
        NULL, NULL, NULL, NULL))
    names(index_data) = header
    
    # figure out which lines of the file correspond to the correct region
    index = which(index_data$chrom == chrom & index_data$start <= stop & 
        index_data$stop >= start)
    
    # raise error, if there isn't any probe data for the CNV region
    if (length(index) == 0) { stop("no probe data for CNV") }
    
    start_index = min(index)
    lines_n = (max(index) - start_index) + 1
    
    all_samples = vector("list", nrow(convex))
    
    for (row in 1:nrow(convex)) {
        sample = convex[row, ]
        if (row %% 50 == 0) {print(paste("loaded", row, "samples of", nrow(convex), sep = " ")) }
        
        # estimate the correct place to start reading from the files
        file_size = file.info(sample$convex_path)$size
        file_pos = ((start_index/nrow(index_data)) * file_size) - 5000
        
        probes = extract_lines_from_file(sample$convex_path, file_pos, lines_n, chrom, start, file_size, l2r_type)
        
        # select the correct rows, and drop the first three columns with the 
        # chrom, start and stop
        probes = probes[probes$chrom == chrom & probes$start <= stop & 
            probes$stop >= start, ]
        
        all_samples[[row]] = probes$score
        
        if (length(probes$score) == 0) {stop("missed samples")}
    }
    
    probes = do.call(cbind, all_samples)
    probes = as.data.frame(probes)
    names(probes) = convex$individual_id
    
    return(probes)
}

#' get probe scores for all the individuals, so that we can pick the values 
#' for every individual at once
#' 
#' @param ddd dataframe listing sample IDs, and file paths
#' @param chrom chromosome of interest as string (eg '1', '2', 'X')
#' @param start chromosome position that the CNV region starts at
#' @param stop chromosome position that the CNV region ends at
#' @export
#' @return a dataframe containing all the L2R data for the different probes in 
#'     the CNV region, for all the particpants in the ddd dataframe
get_probes_data <- function(ddd, chrom, start, stop, l2r_type="adm3", convex_dir="/lustre/scratch113/projects/ddd/users/ddd/DDD_ANALYSIS/convex/data", l2r_dir="/lustre/scratch113/projects/ddd/users/jm33/aggregated_convex") {    
    regions_path = file.path(l2r_dir, paste(l2r_type, "_index.chr",
        chrom, ".txt", sep = ""))
    
    # we have two options, if we have pre-generated single data files per
    # chromosome (using python combine_l2r_data.py), we use a faster function 
    # than having to extract data from all the separate files one by one
    if (file.exists(regions_path)) {
        probes = open_probes_data_fast(ddd, chrom, start, stop, l2r_type, l2r_dir)
    } else {
        if (!("convex_path" %in% names(ddd))) {
           ddd = find_convex_files(ddd, convex_dir)
        }
        probes = open_probes_data_slow(ddd, chrom, start, stop, l2r_type)
    }
    
    return(probes)
}
