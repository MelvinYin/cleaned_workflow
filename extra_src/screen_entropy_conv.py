import numpy as np
import re
from scipy.stats import entropy
from utils import get_composition_conv, conv_meme_parser

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

def main(kwargs):
    # This is for the one from converge
    pssm_fname = kwargs['pssm']
    entropy_threshold = kwargs['entropy_threshold']
    lines_to_keep, pssms = conv_meme_parser(pssm_fname)
    composition = get_composition_conv(lines_to_keep)
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