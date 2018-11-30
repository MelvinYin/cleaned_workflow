import pickle

def meme_rewritter(profiles, fname, to_keep = True):
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

    with open(fname, 'w') as wfile:
        for line in to_write:
            wfile.write(line)
    return

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
