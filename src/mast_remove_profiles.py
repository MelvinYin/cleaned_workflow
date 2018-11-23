import re
import sys

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

def meme_txt_rewritter(to_remove, r_fname, w_fname):
    motif_count = 0
    with open(r_fname, 'r') as rfile:
        with open(w_fname, 'w') as wfile:
            deleting = False
            for i, line in enumerate(rfile):
                if line.startswith("MOTIF") and deleting:
                    deleting = False
                if deleting:
                    continue
                if line.startswith("MOTIF"):
                    motif_count += 1
                    if motif_count in to_remove:
                        deleting = True
                        continue
                wfile.write(line)

def main(kwargs):
    # meme_in meme_out
    to_remove = get_correlated_motifs(kwargs['meme_in'])
    meme_txt_rewritter(to_remove, kwargs['meme_in'], kwargs['meme_out'])