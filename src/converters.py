from collections import OrderedDict
import os
import re

# meme to minimal
def _parse_meme(fname):
    composition = ""
    pssms = []
    in_composition = False
    current_pssm = []
    with open(fname, 'r') as file:
        for line in file:
            if not in_composition \
                and line.startswith("Letter frequencies"):
                in_composition = True
                continue
            if in_composition and line.startswith("Background letter"):
                in_composition = False
                continue
            if in_composition:
                composition += line
                continue
            if line.startswith("letter-probability matrix"):
                current_pssm.append(line)
                continue
            if current_pssm and line.startswith("------------"):
                pssms.append(current_pssm)
                current_pssm = []
                continue
            if current_pssm:
                current_pssm.append(line[1:])   # remove an initial space
    return composition, pssms

def _format_minimal_output_meme(composition, pssms):
    output = []
    output.append("MEME version 4\n\n")
    output.append("ALPHABET= ACDEFGHIKLMNPQRSTVWY\n\n")
    output.append("Background letter frequencies\n")
    output += composition
    output.append("\n")
    for i, pssm in enumerate(pssms):
        output.append(f"MOTIF MEME-{i+1}\n")
        output += pssm
        output.append("\n")
    return output

def meme_to_minimal(kwargs):
    input_fname = kwargs['input']
    output = kwargs['output']
    composition, pssms = _parse_meme(input_fname)
    output_lines = _format_minimal_output_meme(composition, pssms)
    with open(output, 'w') as file:
        for line in output_lines:
            file.write(line)
    return

# Converge output to minimal
# Converts converge motif format to minimal meme format
# see http://meme-suite.org/doc/examples/sample-protein-motif.meme

def _parse_converge_output(filename):
    alphabets = ""
    length = 30
    matrices = OrderedDict()
    matrix = []
    nsite = 0
    matrix_count = 0
    with open(filename, "r") as file:
        for line in file:
            if line.startswith("BEGIN") and matrix_count != 0:
                assert len(matrix) == length, len(matrix)
                motif_name = "MEME-{}".format(matrix_count)
                matrices[motif_name] = (nsite, matrix)
                assert nsite != 0
                matrix = []
                nsite = 0
                continue
            if line.startswith("MATRIX"):
                matrix_count += 1
                match = re.search(r"K=([0-9]+)", line)
                if match is None:
                    raise AssertionError
                nsite = int(match[1])
                continue
            if (line.startswith("50") or line.startswith("30")):
                if not alphabets:
                    matched_alphabets = re.findall("[A-Z]", line)
                    alphabets = "".join(matched_alphabets)
                continue
            if re.match(" [0-9]", line) or re.match("[0-9]+", line):
                probs = re.findall(r"[0-1]\.[0-9]+", line)
                assert len(probs) == len(alphabets)
                matrix.append(probs)
                continue
    return alphabets, matrices

def _parse_converge_composition(filename):
    composition_map = dict()
    with open(filename, "r") as file:
        for line in file:
            if re.match("[A-Z]", line):
                alphabet = line[0]
                composition = line[2:]
                composition_map[alphabet] = float(composition)
                continue
    summed_composition = sum(composition_map.values())
    for key, value in composition_map.items():
        composition_map[key] = value / summed_composition
    return composition_map

def _format_minimal_output_conv(alphabets, composition_map,
                       matrices, output):
    m_to_write = list(range(len(matrices)))
    with open(output, 'w') as file:
        file.write("MEME version 4\n")
        file.write("\n")
        file.write("ALPHABET= " + alphabets + "\n")
        file.write("\n")
        file.write("Background letter frequencies\n")
        for i, alphabet in enumerate(alphabets):
            composition = composition_map[alphabet]
            file.write("{} {} ".format(alphabet, round(composition, 4)))
            if (i != 0) and (i % 9 == 0):
                file.write("\n")
        file.write("\n")
        file.write("\n")
        m_count = 0
        while matrices:
            motif_name, (nsite, matrix) = matrices.popitem(last=False)
            if m_count not in m_to_write:
                m_count += 1
                continue
            m_count += 1
            file.write("MOTIF {}".format(motif_name))
            file.write("\n")
            file.write("letter-probability matrix: alength= 20 w= 30 nsites= {} "
                       "E= 0.000".format(nsite))  # alength = len(alphabets)
            # E is just some random number for now
            # w = width of motif
            file.write("\n")
            for line in matrix:
                to_write = ""
                for prob in line:
                    to_write += prob + " "
                file.write(to_write)
                file.write("\n")
            file.write("\n")

def converge_to_minimal(kwargs):
    # input_pssm=''output.4.matrix.0''
    # composition='composition.txt'
    # output="meme_format.txt"
    input_pssm = kwargs['input_pssm']
    composition = kwargs['composition']
    output = kwargs['output']
    alphabets, matrices = _parse_converge_output(input_pssm)
    composition_map = _parse_converge_composition(composition)
    _format_minimal_output_conv(alphabets, composition_map, matrices, output)

# cons_to_conv_input
# Convert dhcl seed sequences to converge input seqs

def cons_to_conv_input(kwargs):
    seedseq_filename = kwargs['seed_seqs']
    output = kwargs['output']
    to_write = ""

    with open(seedseq_filename, 'r') as rfile:
        for line in rfile:
            to_write += line.strip()

    with open(output, 'w') as wfile:
        wfile.write(">RANDOM\n")
        for i in range(len(to_write) // 60):
            wfile.write(to_write[i*60:(i+1)*60] + "\n")
        if not (len(to_write) % 60 == 0):
            wfile.write(to_write[(len(to_write) // 60) * 60:])

# dhcl_to_cons
# Convert dhcl output to consensus seed sequences

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

def dhcl_to_cons(kwargs):
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



