from collections import defaultdict
from leven import levenshtein
import numpy as np
import operator
import os
import pickle
import re
import shutil
import subprocess

from ui_utils import _convert_url_to_bokeh

class Scorer:
    def __init__(self):
        self.classification_data_path = _convert_url_to_bokeh("combi_fam_data.pkl")
        self.class_data = self._load_classification_data()
        self.motif_path = _convert_url_to_bokeh("motifs.txt")
        self._tmp_seq_store = _convert_url_to_bokeh("tmp_input.fasta")
        self.bash_exec = "/bin/bash"
        self.mast_exec =  _convert_url_to_bokeh("mast")
        self.mast_dir = _convert_url_to_bokeh("mast_folder")
        return

    def _parse_mast_txt(self, input_fname):
        # copied from enolase code
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

    def _load_classification_data(self):
        with open(self.classification_data_path, 'rb') as file:
            class_data = pickle.load(file)
        return class_data

    def _store_to_fasta(self, seq):
        with open(self._tmp_seq_store, 'w') as file:
            file.write(f">User\n")
            for i in range(len(seq) // 80):
                file.write(f"{seq[i * 80:(i + 1) * 80]}\n")
            if len(seq) % 80 != 0:
                file.write(seq[(len(seq) // 80) * 80:])
        return

    def _assign_fam_probs(self, level_scores):
        assigned_probs = defaultdict(float)
        for comb, comb_prob in level_scores.items():
            families_probs = self.class_data[comb]
            for fam, fam_prob in families_probs.items():
                assigned_probs[fam] += 100. * fam_prob * comb_prob
        assert np.isclose(sum(assigned_probs.values()), 100.)

        return assigned_probs

    def _determine_distribution(self, prov_combi):
        level_scores = dict()
        for comb in self.class_data.keys():
            if prov_combi == comb:
                level_scores = dict()
                level_scores[comb] = 1.
                for comb in self.class_data.keys():
                    if comb != prov_combi:
                        level_scores[comb] = 0.
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

    def _score(self):
        if os.path.isdir(self.mast_dir):
            shutil.rmtree(self.mast_dir)
        command = f'{self.mast_exec} {self.motif_path} {self._tmp_seq_store} ' \
            f'-o {self.mast_dir}'
        subprocess.run(command, shell=True, executable=self.bash_exec)
        name_combi = self._parse_mast_txt(f"{self.mast_dir}/mast.txt")

        level_scores = self._determine_distribution(list(name_combi.keys())[0])
        assigned_probs = self._assign_fam_probs(level_scores)

        assigned_probs = sorted(assigned_probs.items(),
                                key=operator.itemgetter(1), reverse=True)
        for i, term in enumerate(assigned_probs):
            assigned_probs[i] = (term[0], str(term[1])[:5]+"%")
        return assigned_probs

    def score_fasta(self, seq_path):
        shutil.copy(seq_path, self._tmp_seq_store)
        assigned_probs = self._score()
        os.remove(self._tmp_seq_store)
        return assigned_probs

    def score_seq(self, seq):
        self._store_to_fasta(seq)
        assigned_probs = self._score()
        os.remove(self._tmp_seq_store)
        return assigned_probs
