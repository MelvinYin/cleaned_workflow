import pickle

def meme_txt_rewritter(to_keep, r_fname, w_fname):
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

def main(kwargs):
    centroid_pkl = kwargs['centroid_pkl']
    input_meme = kwargs['input_meme']
    motifs = kwargs['motifs']
    with open(centroid_pkl, 'rb') as file:
        cluster_centroids = pickle.load(file)

    for label, centroid in cluster_centroids['centroid'].items():
        output = f"{motifs}/motifs_in_cluster_{label}.txt"
        meme_txt_rewritter(centroid, input_meme, output)