import os
from collections import defaultdict
import re

def rename(kwargs):
    motif_filedir = kwargs['motif_filedir']
    output_logodir = kwargs['output_logodir']
    motif_filename_profile_no_map = defaultdict(list)
    for filename in os.listdir(motif_filedir):
        cluster_no = filename[18:-4]
        with open(motif_filedir + "/" + filename) as file:
            for line in file:
                if line.startswith("MOTIF"):
                    motif_no = re.search("MEME-([0-9]+)", line).group(1)
                    motif_filename_profile_no_map[cluster_no].append(motif_no)
    for cluster_no, logo_nos in motif_filename_profile_no_map.items():
        for i, logo_no in enumerate(logo_nos):
            folder_dir = output_logodir + "/cluster_{}/".format(cluster_no)
            orig_filename = folder_dir + "logo_{}.png".format(i+1)
            new_filename = folder_dir + "logos_{}.png".format(logo_no)
            os.rename(orig_filename, new_filename)
    return













