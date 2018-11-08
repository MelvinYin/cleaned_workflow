# Converts converge motif format to minimal meme format
# see http://meme-suite.org/doc/examples/sample-protein-motif.meme

import re
from collections import OrderedDict

def read_converge_mat(filename):
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
                motif_name = filename + "_{}".format(matrix_count)
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

def read_composition(filename):
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

def meme_format_writer(alphabets, composition_map,
                       matrices, m_to_write, filename="meme_format.txt"):
    with open(filename, 'w') as file:
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
                       "E= 0.000001".format(nsite))  # alength = len(alphabets)
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









alphabets, matrices = read_converge_mat('output.4.matrix.0')
composition_map = read_composition('composition.txt')
m_to_write = list(range(len( matrices)))
# m_to_write = [3, 19, 24, 29, 32, 35, 38]
meme_format_writer(alphabets, composition_map, matrices, m_to_write)
