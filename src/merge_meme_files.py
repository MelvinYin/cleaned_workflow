import os

from utils import meme_format_parser, relabel_meme_pssm

def main(kwargs):
    folder = kwargs['meme_folder']
    fname = kwargs['memefile']
    pssm_label_count = 0

    start, pssms, end = meme_format_parser(fname)
    for i, pssm in enumerate(pssms):
        pssms[i] = relabel_meme_pssm(pssm, pssm_label_count)
        pssm_label_count += 1
    for _fname in os.listdir(folder):
        _start, _pssms, _end = meme_format_parser(f"{folder}/{_fname}")
        for i, _pssm in enumerate(_pssms):
            _pssms[i] = relabel_meme_pssm(_pssm, pssm_label_count)
            pssm_label_count += 1
        pssms += _pssms
    with open(fname, 'w') as file:
        for line in start:
            file.write(line)
        for pssm in pssms:
            for line in pssm:
                file.write(line)
        for line in end:
            file.write(line)
