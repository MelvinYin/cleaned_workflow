import pickle
from utils import meme_rewritter
import re

def main(kwargs):
    centroid_pkl = kwargs['centroid_pkl']
    input_meme = kwargs['input_meme']
    motifs = kwargs['motifs']
    with open(centroid_pkl, 'rb') as file:
        cluster_centroids = pickle.load(file)

    for label, centroid in cluster_centroids['centroid'].items():
        output = f"{motifs}/motifs_in_cluster_{label}.txt"
        meme_rewritter(centroid, input_meme, to_keep=True)
        profiles = centroid
        fname = input_meme
        to_keep = True
        if output is None:
            output = fname
        motif_count = 0
        to_write = []
        with open(fname, 'r') as rfile:
            deleting = False
            for i, line in enumerate(rfile):
                if line.startswith("MOTIF") and deleting:
                    deleting = False
                if deleting:
                    continue
                if line.startswith("MOTIF"):
                    motif_count += 1
                    if to_keep and motif_count not in profiles:
                        deleting = True
                        continue
                    elif not to_keep and motif_count in profiles:
                        deleting = True
                        continue

                to_write.append(line)

        with open(output, 'w') as wfile:
            for line in to_write:
                wfile.write(line)