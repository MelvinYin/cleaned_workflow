import os
import re

def _get_loop_endpoints(midpoint, seq_len):
    if midpoint <= 15:
        loop = (0, 30)
    elif seq_len - midpoint <= 15:
        loop = (seq_len-31, seq_len-1)
    else:
        loop = (midpoint-15, midpoint+15)
    return loop

def extract_loops_from_dhcl(filename):
    loops = []
    with open(filename) as file:
        for line in file:
            if line.startswith("LOOPS"):
                terms = re.findall("([0-9]+)\:A\>([0-9]+)\:", line)
                for (start_i, end_i) in terms:
                    loops.append((int(start_i), int(end_i)))
    return loops

def extract_seq_from_fasta(filename):
    merged_seq = ""
    with open(filename) as file:
        next(file)
        for line in file:
            if line.startswith(">"):
                break
            merged_seq += line.strip()  # Remove \n
    return merged_seq

def build_all_loop_indices(seq_len, loops):
    loop_indices = []
    for loop in loops:
        assert loop[1] > loop[0]
        midpoint = int((loop[1] - loop[0]) / 2) + loop[0]
        main_loop = _get_loop_endpoints(midpoint, seq_len)
        loop_indices.append(main_loop)
        if midpoint > 30:
            preced_midpoint = midpoint - 15
            preced_loop = _get_loop_endpoints(preced_midpoint, seq_len)
            loop_indices.append(preced_loop)
        if seq_len - midpoint > 30:
            succ_midpoint = midpoint + 15
            succ_loop = _get_loop_endpoints(succ_midpoint, seq_len)
            loop_indices.append(succ_loop)
    return loop_indices

def match_indices_to_seq(loops, full_seq):
    loop_seqs = []
    for loop in loops:
        assert loop[1] > loop[0]
        seq = full_seq[loop[0]:loop[1]]
        loop_seqs.append(seq)
    return loop_seqs

def main(kwargs):
    # dhcl_dir = "files/from_dhcl"
    # fasta_dir = "files/input_fasta"
    # output = "files/input_seed_seqs.txt"
    dhcl_dir = kwargs['dhcl_dir']
    fasta_dir = kwargs['fasta_dir']
    output = kwargs['output']
    loops = []
    for filename in os.listdir(dhcl_dir):
        if not filename.endswith("dhcl.txt"):
            continue
        dhcl_filepath = f"{dhcl_dir}/{filename}"
        filename_no_suffix = filename.split(".", 2)[0]
        fasta_filepath = f"{fasta_dir}/{filename_no_suffix}.fasta.txt"
        raw_loops = extract_loops_from_dhcl(dhcl_filepath)
        full_seq = extract_seq_from_fasta(fasta_filepath)
        loop_indices = build_all_loop_indices(len(full_seq), raw_loops)
        loop_seqs = match_indices_to_seq(loop_indices, full_seq)
        loops += loop_seqs

    with open(output, "w") as file:
        for loop_seq in loops:
            file.write("{}\n".format(loop_seq))
