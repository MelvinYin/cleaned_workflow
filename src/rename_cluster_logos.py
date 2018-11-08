import os
from collections import defaultdict

MOTIF_FILE_DIRECTORY = "files/motifs"
LOGO_FILE_DIRECTORY = "output/logos"

def main():
    motif_filename_profile_no_map = defaultdict(list)
    for filename in os.listdir(MOTIF_FILE_DIRECTORY):
        with open(MOTIF_FILE_DIRECTORY + "/" + filename) as file:
            for line in file:
                if line.startswith("MOTIF"):
                    motif_no = line[30:].strip()
                    cluster_no = filename[18:-4]
                    motif_filename_profile_no_map[cluster_no].append(motif_no)
    for cluster_no, logo_nos in motif_filename_profile_no_map.items():
        for i, logo_no in enumerate(logo_nos):
            folder_dir = LOGO_FILE_DIRECTORY + "/cluster_{}/".format(cluster_no)
            orig_filename = folder_dir + "logo_{}.png".format(i+1)
            new_filename = folder_dir + "logos_{}.png".format(logo_no)
            os.rename(orig_filename, new_filename)

if __name__ == "__main__":
    main()