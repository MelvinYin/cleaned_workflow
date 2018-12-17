import os
import re
import numpy as np
from scipy.stats import entropy
import pandas as pd

class PSSM:
    def __init__(self, **kwargs):
        self.fname = kwargs['filename']
        start, pssm_lines = self.extract_pssm_lines()
        self.composition = self.parse_composition(start)
        self.pssm_properties = self.extract_pssm_properties(pssm_lines)

        # self.pssms = self.parse_pssms()

    def extract_pssm_properties(self, pssms):
        full_df = pd.DataFrame(columns=('pssm_index', 'evalue', 'entropy',
                                        'pssm', 'nsites'))
        for i, pssm_lines in enumerate(pssms):
            label = None
            evalue = None
            nsites = None
            current = ""
            pssm = []
            for line in pssm_lines:
                if line.startswith("MOTIF"):
                    label = int(re.search("MEME-([0-9]+)", line).group(1))
                if line.startswith("letter-probability"):
                    evalue = re.search(
                        "E= ([0-9]+[\.]?[0-9]?e?[+-]?[0-9]*)", line)
                    evalue = float(evalue.group(1))
                    nsites = int(re.search("nsites= ([0-9]+)", line).group(1))
                if re.match("(1.[0]+)|(0.[0-9]+)", line):
                    current = []
                    probs_re = re.findall("(1.[0]+)|(0.[0-9]+)", line)
                    for (p1, p2) in probs_re:
                        current.append(float(p1 if p1 else p2))
                pssm.append(current)
            # evalue can be == 0 so avoid assert evalue
            assert evalue is not None
            assert label is not None
            assert nsites is not None
            assert pssm
            assert label > 0
            assert nsites > 0
            pssm = np.array(pssm, dtype=float)
            assert pssm.shape == (20, 30)
            _entropy = sum([entropy(line, self.composition) for line in pssm])
            full_df.loc[i, 'pssm_index'] = label
            full_df.loc[i, 'evalue'] = evalue
            full_df.loc[i, 'entropy'] = _entropy
            full_df.loc[i, 'nsites'] = nsites
            full_df.at[i, 'pssm'] = pssm
        full_df.set_index("pssm_index")
        return full_df

    def delete_pssm(self, to_delete):
        pass

    # def parse_pssms(self):
    #     for i, seq in full_df['seqs'].iteritems():
    #         full_df.loc[i, 'num_seqs'] = len(seq.split(" "))
    #
    #     # Drop small combi
    #     to_drop = []
    #     for i, num_seq in full_df['num_seqs'].iteritems():
    #         if num_seq < screen_threshold:
    #             to_drop.append(i)
    #     full_df = full_df.drop(to_drop, axis='index').reset_index()
    #
    #     dist_metric = get_dist_metric(full_df['combi'])
    #     for i in range(len(dist_metric)):
    #         full_df.at[i, 'dist_metric'] = dist_metric[i]
    #
    #     cluster_labels = cluster_metric(dist_metric, n_clusters=20)
    #     full_df['cluster_label'] = cluster_labels

    def relabel_pssms(self, pssm):
        for j, line in enumerate(pssm):
            if re.search("MEME-[0-9]+", line):
                # mast labels their meme starting from 1, so need to match
                line = re.sub("MEME-[0-9]+", f"MEME-{i+1}", line)
            pssm[j] = line
        return pssm

    def get_evalue(self):
        # todo
        evalues = []
        for pssm in self.pssms:
            evalue = None
            for line in pssm:
                if line.startswith("MOTIF"):
                    evalue = re.search(
                        "E= ([0-9]+[\.]?[0-9]?e?[+-]?[0-9]*)", line)
                    evalue = float(evalue.group(1))
                    break
            # evalue can be == 0 so avoid assert evalue
            assert evalue is not None
            evalues.append(evalue)
        return evalues

    def _convert_pssm_to_matrix(self):
        converted_pssms = []
        for pssm in self.pssms:
            current = []
            for line in pssm:
                converted_line = []
                if line.startswith("1.0") or line.startswith("0."):
                    probs_re = re.findall("(1.[0]+)|(0.[0-9]+)", line)
                    for i, (p1, p2) in enumerate(probs_re):
                        prob = p1 if p1 else p2
                        converted_line.append(prob)
                if converted_line:
                    current.append(converted_line)
            converted_pssms.append(current)
        converted_pssms = np.array(converted_pssms, dtype=float)
        assert converted_pssms.shape == (20, 30)
        for line in converted_pssms:
            assert np.isclose(sum(line), 1, atol=0.01)
        return converted_pssms

    def parse_composition(self, start):
        assert start
        composition = []
        at_composition = False
        for line in start:
            if line.startswith("Background"):
                at_composition = True
                continue
            if at_composition:
                freq_re = re.findall("(1.[0]+)|(0.[0-9]+)", line)
                if freq_re:  # skip any \n after frequencies
                    for (f1, f2) in freq_re:
                        freq = f1 if f1 else f2
                        composition.append(float(freq))
        assert len(composition) == 20  # num alphabets
        return composition


    # def get_entropy(self):
    #     converted_pssms = self._convert_pssm_to_matrix()
    #     kl_divs = []
    #     for pssm in converted_pssms:
    #         kl_summed = sum([entropy(line, self.composition) for line in pssm])
    #
    #
    #     def extract_pssm(segment):
    #         pssm = np.zeros((30, 20), dtype=float)
    #         line_no = 0
    #         for line in segment:
    #             if line.startswith("1.0") or line.startswith("0."):
    #                 probs_re = re.findall("(1.[0]+)|(0.[0-9]+)", line)
    #                 for i, (p1, p2) in enumerate(probs_re):
    #                     prob = p1 if p1 else p2
    #                     pssm[line_no][i] = float(prob)
    #                 line_no += 1
    #         assert line_no == 30
    #         for line in pssm:
    #             assert np.isclose(sum(line), 1, atol=0.01)
    #         return pssm
    #
    #     def measure_entropy(pssm, composition):
    #         kl_div = [entropy(line, composition) for line in pssm]
    #         kl_summed = sum(kl_div)
    #         return kl_summed
    #
    #     def main(kwargs):
    #         # This is for the one from converge
    #         pssm_fname = kwargs['pssm']
    #         lines_to_keep, pssms = conv_meme_parser(pssm_fname)
    #         composition = get_composition_conv(lines_to_keep)
    #         entrophies = []
    #         for i, pssm_lines in enumerate(pssms):
    #             pssm = extract_pssm(pssm_lines)
    #             entropy = measure_entropy(pssm, composition)
    #             entrophies.append(entropy)
    #         return entrophies

    # def get_entropy_bits(self):
    #     # TODO: need to change
    #     bits = []
    #     for pssm in self.pssms:
    #         entropy = None
    #         for line in pssm:
    #             if re.match("\(([0-9]+\.[0-9]+) bits\)", line):
    #                 entropy = float(re.match("\(([0-9]+\.[0-9]+) bits\)", line)
    #                                 .group(1))
    #                 break
    #         assert entropy is not None
    #         bits.append(entropy)
    #     return bits

    def get_composition(self):
        assert self.start
        composition = []
        at_composition = False
        for line in self.start:
            # Check if we are in composition
            if line.startswith("Background letter frequencies"):
                at_composition = True
                continue
            if line.startswith("MOTIF"):
                break
            if at_composition:
                freq_re = re.findall("(1.[0]+)|(0.[0-9]+)", line)
                if freq_re:
                    for (f1, f2) in freq_re:
                        # if f1=1.0, the rest == 0, hence f1 first otherwise f2
                        composition.append(float(f1) if f1 else float(f2))
        assert len(composition) == 20  # num alphabets
        return composition

    def extract_pssm_lines(self):
        assert os.path.isfile(self.fname)
        start = []
        pssms = []
        at_start = True
        at_pssm = False
        with open(self.fname, 'r') as file:
            current_pssm = []
            for line in file:
                # Check if we are in start
                if line.startswith("MOTIF") and at_start:
                    at_start = False
                    at_pssm = True
                if at_start:
                    start.append(line)

                # Split into pssms
                if line.startswith("MOTIF") and current_pssm:
                    assert len(current_pssm) > 30
                    pssms.append(current_pssm)
                    current_pssm = [line]
                elif at_pssm:
                    current_pssm.append(line)
        if current_pssm:
            pssms.append(current_pssm)
        assert start
        assert pssms
        return start, pssms

    # def delete_pssm(self, to_delete):
    #     # TODO: deletes and keeps should be done by pssm labels, unless relabel
    #     # is called
    #     # Deletion is only by order and not label
    #     for i in sorted(to_delete, reverse=True):
    #         del self.pssms[i]
    #     return
    #
    # def keep_pssm(self, to_keep):
    #     for i in range(len(self.pssms))[::-1]:
    #         if i not in to_keep:
    #             del self.pssms[i]
    #     return

    def output(self, fname=None):
        if not fname:
            fname = self.fname
        with open(fname, "w") as file:
            for line in self.start:
                file.write(line)
            for pssm in self.pssms:
                for line in pssm:
                    file.write(line)
        return

class PSSM_meme:
    def __init__(self, **kwargs):
        self.fname = kwargs['filename']
        self.start, self.pssms, self.end = self.meme_format_parser()

    def relabel_pssms(self):
        for i, pssm in enumerate(self.pssms):
            for j, line in enumerate(pssm):
                if re.search("MEME-[0-9]+", line):
                    line = re.sub("MEME-[0-9]+", f"MEME-{i}", line)
                pssm[j] = line
            self.pssms[i] = pssm
        return True

    def get_evalue(self):
        evalues = []
        for pssm in self.pssms:
            evalue = None
            for line in pssm:
                if line.startswith("MOTIF"):
                    evalue = re.search(
                        "E= ([0-9]+[\.]?[0-9]?e?[+-]?[0-9]*)", line)
                    evalue = float(evalue.group(1))
                    break
            # evalue can be == 0 so avoid assert evalue
            assert evalue is not None
            evalues.append(evalue)
        return evalues

    def get_entropy_bits(self):
        # TODO: need to change
        bits = []
        for pssm in self.pssms:
            entropy = None
            for line in pssm:
                if re.match("\(([0-9]+\.[0-9]+) bits\)", line):
                    entropy = float(re.match("\(([0-9]+\.[0-9]+) bits\)", line)
                                    .group(1))
                    break
            assert entropy is not None
            bits.append(entropy)
        return bits

    def get_composition(self):
        assert self.start
        composition = []
        at_composition = False
        for line in self.start:
            # Check if we are in composition
            if line.startswith("Letter frequencies in dataset"):
                at_composition = True
            if line.startswith("Background letter frequencies"):
                break
            if at_composition:
                freq_re = re.findall("(1.[0]+)|(0.[0-9]+)", line)
                for (f1, f2) in freq_re:
                    composition.append(float(f1) if f1 else float(f2))
        assert len(composition) == 20  # num alphabets
        return composition

    def meme_format_parser(self):
        assert os.path.isfile(self.fname)
        start = []
        pssms = []
        end = []
        at_start = True
        at_pssm = False
        at_end = False
        with open(self.fname, 'r') as file:
            current_pssm = []
            for line in file:
                # Check if we are in start
                if line.startswith("MOTIF") and at_start:
                    at_start = False
                    at_pssm = True
                if at_start:
                    start.append(line)

                # Split into pssms
                if line.startswith("MOTIF") and current_pssm:
                    assert len(current_pssm) > 30
                    pssms.append(current_pssm)
                    current_pssm = [line]
                elif line.startswith("Stopped"):
                    at_end = True
                    at_pssm = False
                    assert len(current_pssm) > 30
                    pssms.append(current_pssm)
                elif at_pssm:
                    current_pssm.append(line)

                # Check if we are at end
                if at_end:
                    end.append(line)
        assert start
        assert pssms
        assert end
        return start, pssms, end

    def delete_pssm(self, to_delete):
        # Deletion is only by order and not label
        for i in sorted(to_delete, reverse=True):
            del self.pssms[i]
        return

    def keep_pssm(self, to_keep):
        for i in range(len(self.pssms))[::-1]:
            if i not in to_keep:
                del self.pssms[i]
        return

    def output(self, fname=None):
        if not fname:
            fname = self.fname
        with open(fname, "w") as file:
            for line in self.start:
                file.write(line)
            for pssm in self.pssms:
                for line in pssm:
                    file.write(line)
            for line in self.end:
                file.write(line)
        return

class PSSM_conv:
    def __init__(self, kwargs):
        pass
