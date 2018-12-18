import os

from utils import meme_format_parser, relabel_meme_pssm
from pssm_parser import PSSM

# def main(kwargs):
#     # todo: rewrite this once pssm class done
#     folder = kwargs['meme_folder']
#     fname = kwargs['memefile']
#     pssm_label_count = 0
#
#     # start, pssms, end = meme_format_parser(fname)
#     # for i, pssm in enumerate(pssms):
#     #     pssms[i] = relabel_meme_pssm(pssm, pssm_label_count)
#     #     pssm_label_count += 1
#     start = None
#     end = None
#     for _fname in os.listdir(folder):
#         _start, _pssms, _end = meme_format_parser(f"{folder}/{_fname}")
#         for i, _pssm in enumerate(_pssms):
#             _pssms[i] = relabel_meme_pssm(_pssm, pssm_label_count)
#             pssm_label_count += 1
#         if not start:
#             start = _start
#             end = _end
#             pssms = _pssms
#         pssms += _pssms
#     with open(fname, 'w') as file:
#         for line in start:
#             file.write(line)
#         for pssm in pssms:
#             for line in pssm:
#                 file.write(line)
#         for line in end:
#             file.write(line)

def main(kwargs):
    # todo: rewrite this once pssm class done
    folder = kwargs['meme_folder']
    fname = kwargs['memefile']
    main_instance = None
    for _fname in os.listdir(folder):
        pssm_instance = PSSM(filename=f"{folder}/{_fname}")
        if not main_instance:
            main_instance = pssm_instance
            continue
        main_instance.merge_with(pssm_instance)
    main_instance.relabel_pssms()
    main_instance.output(fname)
    return