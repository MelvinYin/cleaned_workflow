from cluster import Cluster
from config import ClusterDir


file = "files/"
log = "files/log.txt"
trash = "files/_trash/"
input_seqs = "files/input_seqs.fasta"
input_pdb = "files/input_pdb_test"
fasta_for_pdb = 'files/input_fasta'
output = "output/"
p2_7_env = '/home/melvin/anaconda3/envs/p2.7/bin/python'
meme_dir = './external_scripts/meme/bin/'
dhcl_exec = './external_scripts/dhcl/executables/everything.py'
bash_exec = '/bin/bash'
single_seq = 'files/single_seq.fasta'
output_clusters = 'output/cluster_description.txt'
output_logos = 'output/logos/'

# cluster_dir = ClusterDir(
#     input_mast=None,
#     input_meme=None,
#     output_mast=output + "mast",
#     description=None,
#     logos=None,
#     cluster_pkl=None)
# Test 1: Only cluster_pkl
# Input: mast, meme
# Output: cluster_pkl


t1_dir = ClusterDir(
    file=None,
    log=None,
    trash=None,
    input_mast=None,
    input_meme=None,
    output_mast=None,
    description=None,
    logos=None,
    cluster_pkl=None)
t1 = Cluster()