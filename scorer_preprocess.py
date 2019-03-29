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

def _get_motif_mapping(motif_filepath):
    with open(motif_filepath, 'r') as file:
        raw_lines = file.readlines()
    motif_map = dict()
    curr_motif_count = 1
    for line in raw_lines:
        if line.startswith("MOTIF"):
            motif_no = int(re.search("[0-9]+", line).group(0))
            motif_map[motif_no] = curr_motif_count
            curr_motif_count += 1
    return motif_map

def _rewrite_motif_txt(motif_filepath, motif_map):
    with open(motif_filepath, 'r') as file:
        raw_lines = file.readlines()
    output_lines = ""
    for line in raw_lines:
        if line.startswith("MOTIF"):
            motif_no = int(re.search("[0-9]+", line).group(0))
            remapped_no = motif_map[motif_no]
            line = f"MOTIF MEME-{remapped_no}\n"
        output_lines += line
    with open(motif_filepath, 'w') as file:
        file.writelines(output_lines)

def _remap_comb_fam(comb_fam, motif_map):
    new_comb_fam = dict()
    for comb, comb_val in comb_fam.items():
        new_comb = []
        for element in comb:
            new_comb.append(motif_map[element])
        new_comb_fam[tuple(new_comb)] = comb_val
    return new_comb_fam

def main():
    input_motif_dir = "./output/motifs"
    input_cluster_descr = "./output/cluster_description.txt"
    pkl_path = "./src/UI/combi_fam_data.pkl"
    assert os.path.isdir(input_motif_dir)
    assert os.path.isfile(input_cluster_descr)
    meme_merger("./output/motifs", "./src/UI/motifs.txt")
    motif_map = _get_motif_mapping("./src/UI/motifs.txt")
    _rewrite_motif_txt("./src/UI/motifs.txt", motif_map)
    combination_families = cluster_descr_parser(input_cluster_descr)
    combination_families = _convert_to_percentage(combination_families)
    combination_families = _remap_comb_fam(combination_families, motif_map)
    with open(pkl_path, 'wb') as file:
        pickle.dump(combination_families, file, -1)
    return

main()