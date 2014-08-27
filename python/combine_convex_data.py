""" aggregates L2R data from different individuals into per chromosome files
"""

from __future__ import division
from __future__ import print_function

import os
import argparse

CHROMS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", \
    "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]

def get_options():
    """ set up parsing of command line options
    """
    
    # hard code a few default folder locations, which can be over-ridden with
    # the command line switches
    L2R_DIR = "/lustre/scratch113/projects/ddd/users/ddd/DDD_ANALYSIS/convex/data"
    DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/1133trios_20131218"
    EXPORT_DIR = "/lustre/scratch113/projects/ddd/users/jm33/aggregated_convex"
    
    parser = argparse.ArgumentParser(description="Aggregate L2R or ADM3 values by chromosome from DDD individuals.")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--chrom", dest="chrom", help="chromosome to examine")
    group.add_argument("--all-chroms", dest="all_chroms", action="store_true", default=False, help="whether to process all the chromosomes")
    
    parser.add_argument("--dir", dest="outdir", default=EXPORT_DIR, help="folder for exporting the datafiles, default = " + EXPORT_DIR)
    parser.add_argument("--type", dest="l2r_type", choices=["adm3", "l2r"], default="adm3", help="whether to extract ADM3 or L2R scores")
    parser.add_argument("--datafreeze", dest="datafreeze_dir", default=DATAFREEZE_DIR, help="directory containing participant info files, in particular \"person_sanger_decipher.private.txt\", default = " + DATAFREEZE_DIR)
    parser.add_argument("--l2r-dir", dest="l2r_dir", default=L2R_DIR, help="directory containing L2R datafiles, default = " + L2R_DIR)
    
    args = parser.parse_args()
    
    columns = {"adm3": 7, "lr2": 3}
    column_num = columns[args.l2r_type]
    
    return args.chrom, args.outdir, args.all_chroms, args.datafreeze_dir, \
        args.l2r_dir, column_num

def get_sanger_ids(datafreeze_dir):
    """get a list of DDD participants, with their family relationships, and
    sanger sample IDs
    
    Returns:
        a dictionary of sanger ID lists (eg [DDD_MAIN10000000], possibly more 
        than one per individual), indexed by their DDD study IDs (eg 
        DDDP1000001).
    """
    
    sanger_ids = {}
    sanger_ids_path = os.path.join(datafreeze_dir, "person_sanger_decipher.private.txt")
    
    with open(sanger_ids_path, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            ddd_id = line[0]
            sanger_id = line[2]
            
            if ddd_id not in sanger_ids:
                sanger_ids[ddd_id] = []
            
            sanger_ids[ddd_id].append(sanger_id)
    
    return sanger_ids

def get_regions_list(path):
    """ start a list of positions using a path to a L2R containing file
    
    Args:
        path: path to a file containing L2R data (tab separated)
    
    Returns:
        a list containing [chrom start, stop] lists for each line in the L2R
        file.
    """
    
    table = []
    
    with open(path, "r") as f:
        for line in f:
            line = line.split("\t")
            table.append(line[0:3])
    
    return table

def get_l2r_paths(l2r_dir, sanger_ids):
    """ get a dictionary of all the paths of L2R containing files
    
    Args:
        l2r_dir: path to directory containing L2R files
        sanger_ids: dictionary of sanger IDs indexed by DDD study ID
    
    Returns:
        dictionary of file names and filesizes, indexed by DDD study ID
    """
    
    # # the convex data is also located on an nfs drive, where the convex files 
    # # are separated into different folders for each individual. In that case we
    # # need to split the person_id into two-char substrings, and include the 
    # # sanger ID before using the ProbeRD[xxx] filename. Here is an example of
    # # a path, in case the L2R files in the lustre folder ever go away.
    # alt_convex_dir = "/nfs/ddd2/ddd_convex_data/DD/DP/10/00/09/DDDP100009/convex/DDD_MAIN5194302"
    
    file_names = {}
    for person_id in sanger_ids:
        alt_ids = sanger_ids[person_id]
        
        for alt_id in alt_ids:
            basename = "ProbeRD_{0}_LR2_GAM.dat".format(alt_id)
            path = os.path.join(l2r_dir, basename)
            
            if os.path.exists(path):
                file_names[person_id] = {}
                file_names[person_id]["path"] = path
                file_names[person_id]["filesize"] = os.path.getsize(path)
                break
    
    return file_names

def write_index(table, outdir, chrom):
    """ make a file from the generated table that we can use in R
    
    Args:
        table: list of [chrom, start stop] lists
        outdir: path to folder to place the index files in
        chrom: path to folder to place the index files in
        
    Returns:
        nothing
    """
    
    header = "chrom\tstart\tstop\n"
    
    basename = "adm3_index.chr" + chrom + ".txt"
    output_path = os.path.join(outdir, basename)      
    output = open(output_path, "w")
    output.write(header)
    
    current_chrom = None
    for line in table:
        current_chrom = line[0]
        
        if current_chrom != chrom:
            continue
        
        line = "\t".join(line) + "\n"
        output.write(line)
    
    output.close()

def add_l2r_for_sample(table, sample_values):
    """ run through the L2R values for a sample, and add them to a table
    
    Args:
        table: list containing [chrom start, stop, l2r_1, l2r_2 ... l2r_n] lists 
        sample_values: list of L2R values for person, in same order as table
    
    Returns:
        the L2R table, but with an additional data point added to each row
    """
    
    table_pos = 0
    for value in sample_values:
        table[table_pos].append(value)
        table_pos += 1
    
    return table

def chrom_num(chrom):
    """ convert a chrom string to a number (make sex chroms 23 and beyond)
    """
    
    chrom_dict = {"X": 23, "Y": 24}
    
    try:
        chrom_num = int(chrom)
    except ValueError:
        chrom_num = chrom_dict[chrom]
    
    return chrom_num

def find_file_pos(handle, start_pos, chrom, start, size):
    """ find where the required data lies within the file
    
    Seek in a file to the section containing the data that we want. We start
    with a reasonable guess about where in the file the data is located, but 
    we seek up or downstream of that (depending on whether the data at that 
    position is before or after the data that we want), then bisect to pin 
    down the region of the file to load.
    
    Args:
        path: path to the L2R datafile
        start_pos: reasonable guess about where the correct data is located in 
            the file
        chrom: chrom that we are searching for
        start: start position of the CNV
        size: size of the L2R file in bytes
    
    Returns:
        the approximate byte offset in the file (to within a few hundred
        bytes) where the required data starts
    """
    
    # seek to the correct part of the file
    cur_pos = handle.seek(start_pos)
    
    read = handle.readlines(200)
    line = read[1].split("\t")
    temp_chrom = line[0]
    temp_start = int(line[1])
    
    # get an initial end point to bisect within
    offset = -100000
    if temp_start <= start or temp_chrom != chrom:
        offset = abs(offset)
    
    # get the upper and lower bounds to bisect within
    lower = min(start_pos, start_pos + offset)
    upper = max(start_pos, start_pos + offset)
    
    # make sure we don't use negative seek positions, or go beyond the file size
    lower = max(lower, 0)
    upper = min(upper, size)
    
    while (upper - lower > 100):
        mid_point = int(lower + (upper - lower)/2)
        cur_pos = handle.seek(mid_point)
        
        read = handle.readlines(200)
        line = read[1].split("\t")
        temp_chrom = line[0]
        temp_start = int(line[1])
        
        # if we have gone to the next chromosome, adjust the mid point until 
        # we hit the correct chrom
        if chrom_num(temp_chrom) < chrom_num(chrom):
            lower = mid_point
            continue
        elif chrom_num(temp_chrom) > chrom_num(chrom):
            upper = mid_point
            continue
        
        # narrow the boundaries to half the previous size
        if temp_start > start:
            upper = mid_point
        elif temp_start <= start:
            lower = mid_point
    
    # make sure we don't use negative seek positions
    seek_position = max(lower - 50, 0)
    
    return seek_position

def get_values_from_file(handle, pos_estimate, chrom, start, keys, column_num, step_size, filesize):
    """ extract L2R values from a file for a given region
    
    Args:
        handle: file handle for current L2R file
        pos_estimate: byte offset guess for where the data is in the L2R file
        chrom: string for current chromosome
        start: nucleotide position for where the required segment starts
        keys: set of (chrom, position) keys for probes from the required lines
        column_num: position of the required value in the L2R line
        step_size: number of lines that we read from the L2R file
    
    Returns:
        list of L2R values, in  same order as the lines in the required region
    """
    
    correct_pos = find_file_pos(handle, pos_estimate, chrom, start, filesize)
    handle.seek(correct_pos)
    
    # after seeking to the given position, we could start in the middle of a
    # line, so drop the first line, except if we are at the very start of the
    # file
    if pos_estimate > 0:
        throwaway = handle.readline()
    
    l2r_values = []
    line_num = 0
    while line_num < step_size + 5:
        line = handle.readline()
        
        if not line: # don't go beyond the end of the file
            break
        
        line = line.strip().split("\t")
        
        # grab the coordinates from the line, which we use to check if the line
        # lies within our desired region
        cur_chrom = line[0]
        cur_start = line[1]
        
        key = (cur_chrom, cur_start)
        if key in keys:
            l2r_values.append(line[column_num])
        
        line_num += 1
    
    handle.close()
    
    return l2r_values

def write_output(outfile, chrom, table, header):
    """ write the output data to a file
    
    Args:
        outfile: file handle for writing data to, or None
        chrom: string for current chromosome
        table: list of rows for current segment of output data
        header: list of entries for header line
    
    Returns:
        nothing
    """
    
    for line in table:
        temp_chrom = line[0]
        
        # if the line contains a different chrom than expected, just ignore it 
        # (and any subsequent lines), since this will mark the transition to a
        # different chromosome
        if temp_chrom != chrom:
            break
        
        outfile.write("\t".join(line) + "\n")

def build_l2r_table_for_chrom(regions, l2r_paths, output_dir, column_num, chrom):
    """ construct L2R tables that aggregate L2R values from multiple individuals
    
    Rather than loading all lines from files at once, we process the L2R files 
    in segments, which reduces the RAM usage (otherwise loading all L2R files
    requires ~ 50 Gb of RAM). We open each file in turn, grab the L2R data from
    it, then add the data into a table, and once all the files have been opened,
    we write the data out, then move onto the next segment of the L2R files.
    
    Args:
        regions: list of [chrom, start, stop] regions in the L2R files
        l2r_paths: dict of L2R paths indexed by DDD sample IDs
        output_dir: folder for creating the aggregated data
        column_num: position of the L2R value in the L2R lines
        chrom: position of the L2R value in the L2R lines
    
    Returns:
        nothing
    """
    
    header = ["chrom", "start", "stop"] + sorted(l2r_paths)
    step_size = 1000
    
    basename = "adm3_dataset.chr" + chrom + ".txt"
    output_path = os.path.join(output_dir, basename)
    
    outfile = open(output_path, "w")
    outfile.write("\t".join(header) + "\n")
    
    pos = 0
    while pos < len(regions):
        current_chrom = regions[pos][0]
        start = int(regions[pos][1])
        
        if current_chrom != chrom:
            pos += 1
            continue
        
        # construct a set of (chrom, position) keys for the probes in the
        # current segment, we can use these to match lines from the L2R files.
        # Also add the probe coordinates to the table.
        keys = set()
        table = []
        for i in range(pos, pos + step_size):
            if i >= len(regions): # don't go past the end of the chroms
                break
            keys.add((regions[i][0], regions[i][1]))
            table.append(regions[i][0:3])
        
        # run through each L2R file and obtain the L2R values for the probes in
        # the curent chromosome segment
        for person_id in sorted(l2r_paths):
            path = l2r_paths[person_id]["path"]
            size = l2r_paths[person_id]["filesize"]
            handle = open(path, "r")
            
            # guess where in the file the data will lie (we get exact later)
            guess = int((pos/len(regions)) * size)
            l2r_values = get_values_from_file(handle, guess, current_chrom, start, keys, column_num, step_size, size)
            table = add_l2r_for_sample(table, l2r_values)
        
        write_output(outfile, current_chrom, table, header)
        pos += step_size
    
    outfile.close()

def main():
    
    chrom, outdir, all_chroms, datafreeze_dir, l2r_dir, column_num = get_options()
    
    sanger_ids = get_sanger_ids(datafreeze_dir)
    l2r_paths = get_l2r_paths(l2r_dir, sanger_ids)
    
    # load the first few columns from a file, which contain the chrom, start and
    # stop positions for each probe region. We'll use these to match positions
    # different files
    first_id = sorted(l2r_paths)[0]
    regions = get_regions_list(l2r_paths[first_id]["path"])
    
    # run through each individual in turn, adding their data to the ends of the
    # existing lists
    if all_chroms:
        for chrom in CHROMS:
            # generate index files to identify lines to load
            write_index(regions, outdir, chrom)
            build_l2r_table_for_chrom(regions, l2r_paths, outdir, column_num, chrom)
    else:
        write_index(regions, outdir, chrom)
        build_l2r_table_for_chrom(regions, l2r_paths, outdir, column_num, chrom)

if __name__ == '__main__':
    main()

