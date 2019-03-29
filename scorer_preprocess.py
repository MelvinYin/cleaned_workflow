import os
import re
import pickle

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

def meme_merger(meme_dir, output_fname):
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
    with open(f"{output_fname}", 'w') as file:
        file.writelines(output_lines)
    return

def main():
    input_motif_dir = "./output/motifs"
    input_cluster_descr = "./output/cluster_description.txt"
    pkl_path = "./UI/combi_fam_data.pkl"
    assert os.path.isdir(input_motif_dir)
    assert os.path.isfile(input_cluster_descr)
    combination_families = cluster_descr_parser(input_cluster_descr)
    combination_families = _convert_to_percentage(combination_families)
    with open(pkl_path, 'wb') as file:
        pickle.dump(combination_families, file, -1)

    meme_merger("./output/motifs", "./UI/motifs.txt")
    return