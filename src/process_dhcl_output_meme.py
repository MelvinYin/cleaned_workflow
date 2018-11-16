import re
import os
import sys
from utils import read_cmd_args

def main(kwargs):
    dhcl_dir = kwargs['dhcl_dir']
    fasta_dir = kwargs['fasta_dir']
    output = kwargs['output']
    input_filenames = []
    for filename in os.listdir(dhcl_dir):
        input_filenames.append(filename)

    filename_loop_map = dict()

    for filename in input_filenames:
        if not filename.endswith("dhcl.txt"):
            continue
        with open(dhcl_dir + "/" + filename) as file:
            loops = []
            for line in file:
                if line.startswith("LOOPS"):
                    terms = re.findall("([0-9]+)\:A\>([0-9]+)\:", line)
                    for (a, b) in terms:
                        loops.append((int(a), int(b)))
        shortened_filename = filename.split(".", 2)[0]
        filename_loop_map[shortened_filename] = loops

    filename_seqs_map = dict()

    for raw_filename in filename_loop_map.keys():
        filename = raw_filename + ".fasta.txt"
        with open(fasta_dir + "/" + filename) as file:
            start = True
            seqs = ""
            for line in file:
                if start:
                    start = False
                    continue
                if line.startswith(">"):
                    break
                seqs += line.strip()  # Remove \n

        filename_seqs_map[raw_filename] = seqs

    # Splitting into segments
    filename_segments_map = dict()
    for filename, loops in filename_loop_map.items():
        seq = filename_seqs_map[filename]
        merged_seq = []
        for loop_i in loops:
            # Add main loop
            raw_start, raw_end = loop_i[0], loop_i[1]
            len_loop = raw_end - raw_start
            if len_loop >= 30:
                median_i = raw_start + int((raw_end - raw_start)/2)
                loop = seq[median_i-15:median_i+15]
                merged_seq.append(loop)
            else:
                median_i = raw_start + int((raw_end - raw_start)/2)
                if median_i - 15 < 0:
                    loop = seq[0:30]
                elif median_i + 15 >= len(seq):
                    loop = seq[-30:]
                else:
                    loop = seq[median_i-15:median_i+15]
                merged_seq.append(loop)
            # Add preceding loop
            if raw_start <= 15:
                loop = seq[:30]
                merged_seq.append(loop)
            elif len(seq) - raw_start <= 20:
                pass
            else:
                loop = seq[raw_start-15:raw_start+15]
                merged_seq.append(loop)

            # Add successive loop
            if len(seq) - raw_end <= 16:
                loop = seq[-30:]
                merged_seq.append(loop)
            elif raw_end < 20:
                pass
            else:
                loop = seq[raw_end-15:raw_end+15]
                merged_seq.append(loop)
        filename_segments_map[filename] = merged_seq

    # Output as fasta for seed seqs
    with open(output, "w") as file:
        for filename, seqs in filename_segments_map.items():
            for seq in seqs:
                assert len(seq) == 30
                file.write("{}\n".format(seq))

if __name__ == "__main__":
    # input_file_directory = "files/from_dhcl"
    # fasta_file_directory = "files/input_fasta"
    # output_filename = "files/consensus_seqs.txt"
    kwargs = read_cmd_args(sys.argv, 'dhcl_dir fasta_dir output')
    main(kwargs)
