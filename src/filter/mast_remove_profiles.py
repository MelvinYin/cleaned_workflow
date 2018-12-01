import re
import sys

from utils import meme_rewritter

def get_correlated_motifs(fname):
    to_remove = []
    with open(fname, 'r') as file:
        for line in file:
            if re.search("Removed motifs", line):
                motifs1 = re.findall("([0-9]+)\, ", line)
                motifs2 = list(re.findall("([0-9]+) and ([0-9]+)", line)[0])
                to_remove = list(int(i) for i in (motifs1 + motifs2))
                break
    return to_remove

def main(kwargs):
    # meme_in meme_out
    to_remove = get_correlated_motifs(kwargs['mast_in'])
    meme_rewritter(to_remove, kwargs['memefile'], to_keep=False)