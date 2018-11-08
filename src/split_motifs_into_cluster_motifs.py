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

def main():
    with open("files/cluster_centroids.pkl", 'rb') as file:
        cluster_centroids = pickle.load(file)

    r_fname = 'files/meme_format3.txt'
    for label, centroid in cluster_centroids.items():
        w_fname = 'files/motifs/motifs_in_cluster_{}.txt'.format(label)
        meme_txt_rewritter(centroid, r_fname, w_fname)

if __name__ == "__main__":
    main()