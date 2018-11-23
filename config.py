from collections import namedtuple

ExecutorDir = namedtuple(
    'ExecutorDir',
    'file log trash input_seqs input_pdb fasta_for_pdb p2_7_env meme_dir '
    'dhcl_exec bash_exec output single_seq output_clusters output_logos '
    'output_mast')

FilterDir = namedtuple(
    'FilterDir',
    'file log trash meme_dir bash_exec single_seq orig cleaned input_seqs')

ClusterDir = namedtuple(
    "ClusterDir", "file log trash input_mast input_meme output_mast "
                  "description logos cluster_pkl meme_dir bash_exec")

class Dir:
    file = "files"
    log = "files/log.txt"
    trash = "files/_trash"
    input_pdb = "files/input_pdb_test"
    fasta_for_pdb = 'files/input_fasta'
    p2_7_env = '/home/melvin/anaconda3/envs/p2.7/bin/python'
    dhcl_exec = './external_scripts/dhcl/executables/everything.py'
    bash_exec = '/bin/bash'
    meme_dir = './external_scripts/meme/bin'
    input_seqs = "files/input_seqs.fasta"
    single_seq = 'files/single_seq.fasta'
    output = "output"
    output_clusters = 'output/cluster_description.txt'
    output_logos = 'output/logos'
    output_mast = 'output/mast'

    executor_dir = ExecutorDir(
        file=file,
        log=log,
        trash=trash,
        input_seqs=input_seqs,
        input_pdb=input_pdb,
        fasta_for_pdb=fasta_for_pdb,
        p2_7_env=p2_7_env,
        meme_dir=meme_dir,
        output=output,
        dhcl_exec=dhcl_exec,
        bash_exec=bash_exec,
        single_seq=single_seq,
        output_clusters=output_clusters,
        output_logos=output_logos,
        output_mast=output_mast)
