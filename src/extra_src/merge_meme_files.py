import os

from pssm_parser import PSSM

def main(kwargs):
    folder = kwargs['pssm_folder']
    fname = kwargs['output']
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