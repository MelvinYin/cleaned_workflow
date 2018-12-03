import pickle
from utils import meme_rewritter
import sys

def main(kwargs):
    cluster_df_pkl = kwargs['cluster_df_pkl']
    memefile = kwargs['memefile']
    with open(cluster_df_pkl, 'rb') as file:
        cluster_df = pickle.load(file)

    centroids = cluster_df['centroid']
    profiles_to_keep = set()
    for centroid in centroids:
        for profile in centroid:
            profiles_to_keep.add(profile)
    meme_rewritter(profiles_to_keep, memefile, to_keep=True)
    return True
