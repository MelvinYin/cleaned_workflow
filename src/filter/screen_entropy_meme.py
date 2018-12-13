import numpy as np
import re
from scipy.stats import entropy

def extract_pssm(segment):
    pssm = np.zeros((30, 20), dtype=float)
    line_no = 0
    for line in segment:
        if line.startswith("1.0") or line.startswith("0."):
            probs_re = re.findall("(1.[0]+)|(0.[0-9]+)", line)
            for i, (p1, p2) in enumerate(probs_re):
                prob = p1 if p1 else p2
                pssm[line_no][i] = float(prob)
            line_no += 1
    assert line_no == 30
    for line in pssm:
        assert np.isclose(sum(line), 1, atol=0.01)
    return pssm

def measure_entropy(pssm, composition):
    kl_div = [entropy(line, composition) for line in pssm]
    kl_summed = sum(kl_div)
    return kl_summed

def split_filelines(fname):
    start = []
    pssms = []
    end = []
    pssm_start = False
    end_start = False
    current_pssm = []
    with open(fname, 'r') as file:
        for line in file:
            if line.startswith("MOTIF"):
                pssm_start = True
            if not pssm_start:
                start.append(line)
                continue
            if line.startswith('MOTIF') and current_pssm:
                assert len(current_pssm) > 30, current_pssm
                pssms.append(current_pssm)
                current_pssm = [line]
                continue
            if line.startswith("Stopped"):
                assert len(current_pssm) > 30  # width of motif
                pssms.append(current_pssm)
                end_start = True
            if end_start:
                end.append(line)
                continue
            current_pssm.append(line)
    return start, pssms, end

def get_composition(lines):
    composition = []
    start_composition = False
    for line in lines:
        if line.startswith("Background"):
            start_composition = True
            continue
        if start_composition:
            frequencies = re.findall("(1.[0]+)|(0.[0-9]+)", line)
            if frequencies: # skip any \n after frequencies
                for (f1, f2) in frequencies:
                    freq = f1 if f1 else f2
                    composition.append(float(freq))
    assert len(composition) == 20   # num alphabets
    return composition

def main(kwargs):
    # This is for the one from converge
    pssm_fname = kwargs['pssm']
    entropy_threshold = kwargs['entropy_threshold']
    lines_to_keep, pssms = split_filelines(pssm_fname)
    composition = get_composition(lines_to_keep)
    for i, pssm_lines in enumerate(pssms):
        pssm = extract_pssm(pssm_lines)
        entropy = measure_entropy(pssm, composition)
        print("{} | {}".format(i, entropy))
        if entropy > entropy_threshold:
            lines_to_keep += pssm_lines
    with open(pssm_fname, 'w') as file:
        for line in lines_to_keep:
            file.write(line)
    return