import re
import pickle

def meme_rewritter(to_keep, r_fname, w_fname):
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
                    if motif_count not in to_keep:
                        deleting = True
                        continue
                wfile.write(line)
    return True

def main(kwargs):
    # cluster_df_pkl=files/cluster_final.pkl
    # input=files/meme_format2.txt
    # output=files/meme_format3.txt
    cluster_df_pkl = kwargs['cluster_df_pkl']
    input_memefile = kwargs['input']
    output = kwargs['output']
    with open(cluster_df_pkl, 'rb') as file:
        cluster_df = pickle.load(file)

    # All profiles
    centroids = cluster_df['centroid']
    # pd.Series([(9, 6, 11, 2, 1, 7, 10, 5, 32, 3, 8),...])
    profiles_to_keep = set()
    for centroid in centroids:
        for profile in centroid:
            profiles_to_keep.add(profile)
    meme_rewritter(profiles_to_keep, input_memefile, output)
    return True

# with open("files/profiles_to_remove.pkl", 'rb') as file:
#     to_remove = pickle.load(file)
# meme_txt_rewritter(to_remove, r_fname="files/meme_format2.txt",
#                        w_fname="files/meme_format3.txt")
# meme_txt_rewritter([], r_fname="meme_format123.txt",
#                        w_fname="meme_format_after_tomtom2.txt")