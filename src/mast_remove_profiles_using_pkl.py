import re
import pickle

def get_correlated_motifs(fname='files/mast.txt'):
    to_remove = []
    with open(fname, 'r') as file:
        for line in file:
            if re.search("Removed motifs", line):
                motifs1 = re.findall("([0-9]+)\, ", line)
                motifs2 = list(re.findall("([0-9]+) and ([0-9]+)", line)[0])
                to_remove = list(int(i) for i in (motifs1 + motifs2))
                break
    return to_remove

def meme_txt_rewritter(to_remove, r_fname="meme_format2.txt",
                       w_fname="meme_format123.txt"):
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

# meme_txt_rewritter([3, 6, 7, 12, 14, 16, 25, 36, 38, 39, 46, 47, 54, 56, 58,
#                     59, 60, 61, 64, 69, 76, 77, 78, 82, 88, 89, 90, 91, 92,
#                     94, 98, 101, 102, 107, 108, 109])

with open("files/profiles_to_remove.pkl", 'rb') as file:
    to_remove = pickle.load(file)
meme_txt_rewritter(to_remove, r_fname="files/meme_format2.txt",
                       w_fname="files/meme_format3.txt")
# meme_txt_rewritter([], r_fname="meme_format123.txt",
#                        w_fname="meme_format_after_tomtom2.txt")