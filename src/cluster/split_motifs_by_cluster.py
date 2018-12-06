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
        meme_rewritter(centroid, input_meme, output=output, to_keep=True,
                       sub=True)