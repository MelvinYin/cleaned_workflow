import os
import re
import numpy as np
from scipy.stats import entropy
import pandas as pd

class PSSM:
    def __init__(self, **kwargs):
        self.fname = kwargs['filename']
        self.start, self.pssm_lines = self.extract_pssm_lines()
        self.composition = self.parse_composition(self.start)
        self.pssm_properties = self.extract_pssm_properties(self.pssm_lines)

    def extract_pssm_properties(self, pssms):
        # todo: converge2meme needs to be called for converge,
        #  and meme_to_minimal need to be called for meme.
        full_df = pd.DataFrame(columns=('pssm_index', 'evalue', 'entropy',
                                        'pssm', 'nsites'))
        for i, pssm_lines in enumerate(pssms):
            label = None
            evalue = None
            nsites = None
            pssm = []
            for line in pssm_lines:
                if line.startswith("MOTIF"):
                    label = int(re.search("MEME-([0-9]+)", line).group(1))
                if line.startswith("letter-probability"):
                    evalue = re.search("E= ([0-9]+[\.]?[0-9]?e?[+-]?[0-9]*)",
                                       line)
                    evalue = float(evalue.group(1))
                    nsites = int(re.search("nsites= ([0-9]+)", line).group(1))
                if re.match("(1\.[0]+)|(0\.[0-9]+)", line):
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
            assert pssm.shape == (30, 20)
            _entropy = sum([entropy(line, self.composition) for line in pssm])
            full_df.loc[i, 'pssm_index'] = label
            full_df.loc[i, 'evalue'] = evalue
            full_df.loc[i, 'entropy'] = _entropy
            full_df.loc[i, 'nsites'] = nsites
            full_df.at[i, 'pssm'] = pssm
        full_df = full_df.set_index("pssm_index")
        return full_df

    def relabel_pssms(self):
        length = len(self.pssm_properties)
        new_index = np.array(range(1, length+1))    # Fail if is list
        self.pssm_properties.set_index(new_index, inplace=True)
        return

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

    def get_output_lines(self):
        lines = []
        alphabets = 'ACDEFGHIKLMNPQRSTVWY'
        lines.append("MEME version 4\n\n")
        lines.append("ALPHABET= " + alphabets + "\n\n")
        lines.append("Background letter frequencies\n")
        for i, (alphabet, comp) in enumerate(zip(alphabets, self.composition)):
            comp = "{:.4f}".format(comp)
            lines.append(f"{alphabet} {comp} ")
            if (i != 0) and (i % 9 == 0):
                lines.append("\n")
        lines.append("\n\n")
        for motif_i, values in self.pssm_properties.iterrows():
            evalue = "{:.6f}".format(values.evalue)
            nsites = int(values.nsites)
            lines.append("MOTIF MEME-{}\n".format(motif_i))
            lines.append(f"letter-probability matrix: alength= 20 w= 30 "
                         f"nsites= {nsites} E= {evalue}\n")
            for line in values.pssm:
                line_to_write = ""
                for prob in line:
                    line_to_write += "{:.6f}".format(prob) + " "
                line_to_write += '\n'
                lines.append(line_to_write)
            lines.append("\n")
        return lines

    def output(self, fname=None):
        if not fname:
            fname = self.fname
        with open(fname, "w") as file:
            lines = self.get_output_lines()
            for line in lines:
                file.write(line)
        return

    def delete(self, i_to_delete):
        self.pssm_properties = self.pssm_properties.drop(i_to_delete)
        return

    def keep(self, i_to_keep):
        indices = list(self.pssm_properties.index.values)
        to_delete = []
        for index in indices[::-1]:
            if index not in i_to_keep:
                to_delete.append(index)
        self.delete(to_delete)
        return

    def get_evalue(self):
        return self.pssm_properties.evalue.values

    def get_entropy(self):
        return self.pssm_properties.entropy.values



# def get_evalue(self):
#     # todo
#     evalues = []
#     for pssm in self.pssm_properties.pssm:
#         evalue = None
#         for line in pssm:
#             if line.startswith("MOTIF"):
#                 evalue = re.search(
#                     "E= ([0-9]+[\.]?[0-9]?e?[+-]?[0-9]*)", line)
#                 evalue = float(evalue.group(1))
#                 break
#         # evalue can be == 0 so avoid assert evalue
#         assert evalue is not None
#         evalues.append(evalue)
#     return evalues

# def delete_pssm(self, to_delete):
#     self.pssm_properties = self.pssm_properties.drop(to_delete)
#     return
#
#
# def relabel_pssms(self, pssm):
#     for j, line in enumerate(pssm):
#         if re.search("MEME-[0-9]+", line):
#             # mast labels their meme starting from 1, so need to match
#             line = re.sub("MEME-[0-9]+", f"MEME-{i + 1}", line)
#         pssm[j] = line
#     return pssm


# def get_composition(self):
#     assert self.start
#     composition = []
#     at_composition = False
#     for line in self.start:
#         # Check if we are in composition
#         if line.startswith("Background letter frequencies"):
#             at_composition = True
#             continue
#         if line.startswith("MOTIF"):
#             break
#         if at_composition:
#             freq_re = re.findall("(1.[0]+)|(0.[0-9]+)", line)
#             if freq_re:
#                 for (f1, f2) in freq_re:
#                     # if f1=1.0, the rest == 0, hence f1 first otherwise f2
#                     composition.append(float(f1) if f1 else float(f2))
#     assert len(composition) == 20  # num alphabets
#     return composition