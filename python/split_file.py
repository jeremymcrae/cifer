
import os
import argparse

def get_options():
    
    parser = argparse.ArgumentParser(description="split a file")
    
    parser.add_argument("--in", dest="path", help="file to split")
    parser.add_argument("-n", "--num_lines", dest="num_lines", help="number of lines per file")
    args = parser.parse_args()
    
    return args.path, int(args.num_lines)

def split_file(path, num_lines):
    """ splits a file into  a number of files, all with the same number of lines
    while retaining the header line for all files.
    """
    
    file_num = 1
    with open(path, "r") as f:
        header = f.readline()
        
        lines_in_file = 0
        for line in f:
            if lines_in_file == 0:
                out_path = "{0}.{1}.txt".format(path, file_num)
                output = open(out_path, "w")
                output.write(header)
                file_num += 1
            
            output.write(line)
            lines_in_file += 1
            
            if lines_in_file > num_lines:
                lines_in_file = 0

def main():
    
    path, num_lines = get_options()
    split_file(path, num_lines)

if __name__ == '__main__':
    main()
