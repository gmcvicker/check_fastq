
import gzip
import sys
import os
import re

from fastq import read_fastq, write_fastq


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("usage: %s <left_fastq> <right_fastq>\n"
                         % sys.argv[0])
        exit(2)

    left_fastq = sys.argv[1]
    right_fastq = sys.argv[2]

    left_reader = read_fastq(left_fastq)
    right_reader = read_fastq(right_fastq)

    # base the output filenames on the input filenames
    left_out_prefix = re.sub(r".txt(.gz)?$", "", left_fastq)
    right_out_prefix = re.sub(r".txt(.gz)?$", "", right_fastq)

    left_out_filename = "%s.fixed.txt.gz" % (left_out_prefix)
    right_out_filename = "%s.fixed.txt.gz" % (right_out_prefix)
    unmatched_left_out_filename = "%s.unmatched.txt.gz" % (left_out_prefix)
    unmatched_right_out_filename = "%s.unmatched.txt.gz" % (right_out_prefix)

    sys.stderr.write("writing to files:\n  %s\n  %s\n" %
                     (left_out_filename, right_out_filename))

    left_outfile = gzip.open(left_out_filename, "wb")
    right_outfile = gzip.open(right_out_filename, "wb")
    unmatched_left_outfile = gzip.open(unmatched_left_out_filename, "wb")
    unmatched_right_outfile = gzip.open(unmatched_right_out_filename, "wb")

    right_unmatched = {}
    left_unmatched = {}

    line_num = 0
    
    n_matched = 0

    right_done = False
    left_done = False

    while (not left_done) or (not right_done):
        # try to get read from left lane
        if not left_done:
            try:
                left_lines = left_reader.next()
            except StopIteration:
                sys.stderr.write("left lane ended\n")
                left_done = True
        
        # try to get read from right lane
        if not right_done:
            try:
                right_lines = right_reader.next()
            except StopIteration:
                sys.stderr.write("right lane ended\n")
                right_done = True
        
        if line_num % 100000 == 0:
            sys.stderr.write(".")

        # strip off /1 or /2 from end of read names so they can be compared
        if not left_done:
            left_read_name = re.sub(r"/\d(\s)*$", "", 
                                    left_lines[0].split()[0])

        if not right_done:
            right_read_name = re.sub(r"/\d(\s)*$", "", 
                                     right_lines[0].split()[0])


        if ((not left_done) and (not right_done) and 
            (left_read_name == right_read_name)):
            # reads match
            n_matched += 1
            write_fastq(left_outfile, left_lines)
            write_fastq(right_outfile, right_lines)
        else:
            # reads do not match
            
            if not left_done:
                if left_read_name in right_unmatched:
                    write_fastq(left_outfile, left_lines)
                    write_fastq(right_outfile, right_unmatched[left_read_name])
                    n_matched += 1
                    del right_unmatched[left_read_name]
                else:
                    if left_read_name in left_unmatched:
                        raise ValueError("duplicate left read %s" %
                                         left_read_name)
                    left_unmatched[left_read_name] = left_lines
            
            if not right_done:
                if right_read_name in left_unmatched:
                    # sys.stderr.write("found match for right read\n")
                    write_fastq(left_outfile, left_unmatched[right_read_name])
                    write_fastq(right_outfile, right_lines)
                    n_matched += 1
                    del left_unmatched[left_read_name]
                else:
                    if right_read_name in right_unmatched:
                        raise ValueError("duplicate right read %s" %
                                         right_read_name)
                    right_unmatched[right_read_name] = right_lines

        line_num += 4

    sys.stderr.write("\n")

    # write out the unmatched records
    for left_read in left_unmatched.values():
      write_fastq(unmatched_left_outfile, left_read)
    for right_read in right_unmatched.values():
      write_fastq(unmatched_right_outfile, right_read)

    sys.stderr.write("done")
    sys.stderr.write("matched:%d, unmatched_left:%d, unmatched_right:%d\n"
                     % (n_matched, len(left_unmatched),
                        len(right_unmatched)))

    left_outfile.close()
    right_outfile.close()
    unmatched_left_outfile.close()
    unmatched_right_outfile.close()

                
                

main()
