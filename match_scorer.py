import os
import re
import numpy as np
from leven import levenshtein
import subprocess
from collections import defaultdict



def cluster_descr_parser(src):
    with open(src, 'r') as file:
        raw_lines = file.readlines()
    combination_families = dict()
    curr_combination = None
    curr_groups = None
    for line in raw_lines:
        if line.startswith("Combination") and curr_groups is None:
            matched_values = re.findall("([0-9]+)", line)
            curr_combination = tuple([int(i) for i in matched_values])
            curr_groups = dict()
            continue
        if line.startswith("Combination"):
            combination_families[curr_combination] = curr_groups
            curr_groups = dict()
            matched_values = re.findall("([0-9]+)", line)
            curr_combination = tuple([int(i) for i in matched_values])
            continue
        if re.search(" : ", line) is not None:
            group_id = re.match("[A-Z]+", line).group(0)
            count = re.search("[0-9]+", line).group(0)
            curr_groups[group_id] = int(count)
    if curr_groups:
        combination_families[curr_combination] = curr_groups
    return combination_families

def _convert_to_percentage(combination_families):
    for key in combination_families.keys():
        new_curr_groups = dict()
        old_curr_groups = combination_families[key]
        total_count = sum(old_curr_groups.values())
        for group_id, count in old_curr_groups.items():
            new_curr_groups[group_id] = float(count) / total_count
        combination_families[key] = new_curr_groups
    return combination_families



def _determine_distribution(prov_combi, comb_fam):
    level_scores = dict()
    for comb in comb_fam.keys():
        if prov_combi == comb_fam:
            level_scores = dict()
            level_scores[comb] = 1.
            return level_scores
        unicode_str_a = ""
        for element in prov_combi:
            unicode_str_a += chr(element)
        unicode_str_b = ""
        for element in comb:
            unicode_str_b += chr(element)
        dist = levenshtein(unicode_str_a, unicode_str_b)
        level_scores[comb] = 1. / dist
    # normalise
    tot_score = sum(level_scores.values())
    for comb in level_scores.keys():
        level_scores[comb] = level_scores[comb] / tot_score
    return level_scores

def _assign_fam_probs(level_scores, combination_families):
    assigned_probs = defaultdict(float)
    for comb, comb_prob in level_scores.items():
        families_probs = combination_families[comb]
        for fam, fam_prob in families_probs.items():
            assigned_probs[fam] += fam_prob * comb_prob
    return assigned_probs



# given a sequence, first convert it to a tmp fasta file. Then, run it
# through mast. Before that, make a function that converts all the meme in
# clusters into a single meme file. Run mast using that meme file and the
# fasta file. Get the resultant combination from the mast.txt.

def meme_merger(meme_dir):
    output_fname = 'motifs.txt'
    output_dir = "."
    output_lines = ''
    added_motifs = set()
    for filename in os.listdir(meme_dir):
        to_copy = False
        with open(f"{meme_dir}/{filename}", 'r') as file:
            if not output_lines:
                for line in file:
                    if line.startswith('MOTIF'):
                        motif_no = int(re.search("[0-9]+", line).group(0))
                        added_motifs.add(motif_no)
                    output_lines += line
                continue
            for line in file:
                if line.startswith('MOTIF'):
                    motif_no = int(re.search("[0-9]+", line).group(0))
                    if motif_no not in added_motifs:
                        added_motifs.add(motif_no)
                        to_copy = True
                    else:
                        to_copy = False
                if to_copy:
                    output_lines += line
    with open(f"{output_dir}/{output_fname}", 'w') as file:
        file.writelines(output_lines)


def parse_mast_txt(input_fname):
    first_start = False
    second_start = False
    name_combi = defaultdict(str)
    with open(input_fname, 'r') as file:
        for line in file:
            if line.startswith("SECTION II"):
                first_start = True
                continue
            if first_start and line.startswith("-------------"):
                second_start = True
                continue
            if first_start and second_start:
                if line == "\n":
                    break
                name, remainder = line.split(" ", maxsplit=1)
                combi = re.findall("\[[0-9]+\]", remainder)
                combi_int = tuple([int(term[1:-1]) for term in combi])
                name_combi[combi_int] += name + " "
    name_combi = dict(name_combi)
    return name_combi


def _convert_to_fasta(seq, output_path):
    with open(output_path, 'w') as file:
        file.write(f">User\n")
        for i in range(len(seq) // 80):
            file.write(f"{seq[i*80:(i+1)*80]}\n")
        if len(seq) % 80 != 0:
            file.write(seq[(len(seq)//80)*80:])

# _convert_to_fasta("MVKITRLTTYRLPPRWMFLKVETDEGVTGWGEPVIEGRARTVEAAVHELSDYLIGQDPSR"
#                   "INDLWQTMYRAGFYRGGPILMSAIAGIDQALWDIKGKVLGVPVYELLGGLVRDKMRTYSW"
#                   "VGGDRPADVIAGMKALQAGGFDHFKLNGCEEMGIIDTSRAVDAAVARVAEIRSAFGNTVE"
#                   "FGLDFHGRVSAPMAKVLIKELEPYRPLFIEEPVLAEQAETYARLAAHTHLPIAAGERMFS"
#                   "RFDFKRVLEAGGVSILQPDLSHAGGITECVKIAAMAEAYDVALAPHCPLGPIALAACLHV"
#                   "DFVSWNATLQEQSMGIHYNKGAELLDYVRNKADFALEGGYIRPPRLPGLGVDIDEALVIE"
#                   "RSKEAPDWRNPVWRHADGSVAEWAENLYFQSHHHHHHWSHPQFEK", './input_fasta.fasta')


input_motif_dir = "./output/motifs"
input_cluster_descr = "./output/cluster_description.txt"
assert os.path.isdir(input_motif_dir)
assert os.path.isfile(input_cluster_descr)

combination_families = cluster_descr_parser(input_cluster_descr)
combination_families = _convert_to_percentage(combination_families)

meme_merger("./output/motifs")

mast_dir = 'external_scripts/meme/bin/mast'
input_seq = './input_fasta.fasta'
motif_file = "./motifs.txt"
output_mast_dir = "./mast"
# os.mkdir(output_mast_dir)
command = f'{mast_dir} {motif_file} {input_seq} -o {output_mast_dir}'
bash_exec = "/bin/bash"

subprocess.run(command, shell=True, executable=bash_exec)


name_combi = parse_mast_txt(f"{output_mast_dir}/mast.txt")
level_scores = _determine_distribution(list(name_combi.keys())[0],
                                       combination_families)
assigned_probs = _assign_fam_probs(level_scores, combination_families)
assert np.isclose(sum(assigned_probs.values()), 1.)
print(assigned_probs)
