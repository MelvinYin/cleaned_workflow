import re
import os

def get_dhcl_segments(kwargs):
    dhcl_dir = kwargs['dhcl_dir']
    fasta_dir = kwargs['fasta_dir']
    output_filename = kwargs['output']

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
        merged_seq = ""
        for loop_i in loops:
            raw_start, raw_end = loop_i[0], loop_i[1]
            len_loop = raw_end - raw_start
            if len_loop >= 30:
                median_i = raw_start + int((raw_end - raw_start)/2)
                loop = seq[median_i-15:median_i+15]
                merged_seq += loop
            else:
                median_i = raw_start + int((raw_end - raw_start)/2)
                if median_i - 15 < 0:
                    loop = seq[0:30]
                elif median_i + 15 >= len(seq):
                    loop = seq[-30:]
                else:
                    loop = seq[median_i-15:median_i+15]
                merged_seq += loop
        filename_segments_map[filename] = merged_seq

    # Output as fasta for seed seqs
    with open(output_filename, "w") as file:
        for filename, seq in filename_segments_map.items():
            file.write(">{}\n".format(filename))
            if len(seq) % 60 == 0:
                for i in range(int(len(seq) // 60)):
                    file.write(seq[i*60:(i+1)*60])
                    file.write("\n")
            else:
                for i in range(int(len(seq) // 60)):
                    file.write(seq[i*60:(i+1)*60])
                    file.write("\n")
                file.write(seq[int(len(seq) // 60) * 60:])
                file.write("\n")

def main(kwargs):
    # dhcl_filedir = "files/from_dhcl"
    # fasta_filedir = "files/input_fasta"
    # output = "files/init_seed_seqs.fasta"
    get_dhcl_segments(kwargs)
    return
