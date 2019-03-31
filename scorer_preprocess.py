import os
import re
import pickle

family_initial_to_name_map = dict(ENL='Enolase',
                                  GAD='Galactarate Dehydratase',
                                  GLD='Glucarate Dehydratase',
                                  MR='Mandelate Racemase',
                                  MD='Mannonate Dehydratase',
                                  MAL='Methylaspartate Ammonia-lyase',
                                  MC='Muconate Cycloisomerase')

def parse_cluster_descr(src):
    with open(src, 'r') as file:
        raw_lines = file.readlines()
    combi_famfreq = dict()
    combi = None
    family_freqs = None
    for line in raw_lines:
        if line.startswith("Combination") and family_freqs is None:
            combi_str = re.findall("([0-9]+)", line)
            combi = tuple([int(i) for i in combi_str])
            family_freqs = dict()
            continue
        if line.startswith("Combination"):
            combi_famfreq[combi] = family_freqs
            family_freqs = dict()
            combi_str = re.findall("([0-9]+)", line)
            combi = tuple([int(i) for i in combi_str])
            continue
        if re.search(" : ", line) is not None:
            family_initials = re.match("[A-Z]+", line).group(0)
            family = family_initial_to_name_map[family_initials]
            freq = re.search("[0-9]+", line).group(0)
            family_freqs[family] = int(freq)
    if family_freqs:
        combi_famfreq[combi] = family_freqs
    return combi_famfreq

def _convert_to_percent(combi_famfreq):
    combi_famprob = dict()
    for combi in combi_famfreq.keys():
        famprob = dict()
        famfreq = combi_famfreq[combi]
        total_freq = sum(famfreq.values())
        for fam, freq in famfreq.items():
            famprob[fam] = float(freq) / total_freq
        combi_famprob[combi] = famprob
    return combi_famprob

def merge_memefiles(meme_dir, output_fname):
    output_lines = ''
    added_motifs = set()
    to_copy = True
    # So initial lines of first file is copied, but not subsequent ones,
    # see to_copy=False at end of for loop
    for filename in os.listdir(meme_dir):
        with open(f"{meme_dir}/{filename}", 'r') as file:
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
        to_copy = False
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
    for comb, value in comb_fam.items():
        new_comb = []
        for element in comb:
            new_comb.append(motif_map[element])
        new_comb_fam[tuple(new_comb)] = value
    return new_comb_fam

def main():
    input_motif_dir = "./output/motifs"
    otuput_motif_file = "./src/UI/static/motifs.txt"
    input_cluster_descr = "./output/cluster_description.txt"
    pkl_path = "./src/UI/static/combi_fam_data.pkl"
    assert os.path.isdir(input_motif_dir)
    assert os.path.isfile(input_cluster_descr)
    merge_memefiles(input_motif_dir, otuput_motif_file)
    motif_map = _get_motif_mapping(otuput_motif_file)
    _rewrite_motif_txt(otuput_motif_file, motif_map)
    combi_famfreq = parse_cluster_descr(input_cluster_descr)
    combi_famprob = _convert_to_percent(combi_famfreq)
    combi_famprob = _remap_comb_fam(combi_famprob, motif_map)
    with open(pkl_path, 'wb') as file:
        pickle.dump(combi_famprob, file, -1)
    return

main()