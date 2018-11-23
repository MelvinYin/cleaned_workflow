from collections import namedtuple

ExecutorDirectory = namedtuple(
    'ExecutorDirectory',
    "file log trash input_seqs input_pdb fasta_for_pdb p2_7_env meme_dir "
    "dhcl_exec bash_exec output single_seq output_clusters output_logos "
    "output_mast")

FilterDirectory = namedtuple(
    'FilterDirectory',
    "file log trash meme_dir input_seqs single_seq bash_exec orig cleaned")

ClusterDirectory = namedtuple(
    "ClusterDirectory",
    "file log trash meme_dir bash_exec input_mast input_meme output_mast "
    "description logos cluster_pkl")

class Directory:
    file = "files"
    output = "output"
    p2_7_env = "/home/melvin/anaconda3/envs/p2.7/bin/python"
    dhcl_exec = "./external_scripts/dhcl/executables/everything.py"
    bash_exec = "/bin/bash"
    meme_dir = "./external_scripts/meme/bin"
    log = f"{file}/log.txt"
    trash = f"{file}/_trash"
    input_seqs = f"{file}/input_seqs.fasta"
    single_seq = f"{file}/single_seq.fasta"
    input_pdb = f"{file}/input_pdb_test"
    fasta_for_pdb = f"{file}/input_fasta"
    output_clusters = f"{output}/cluster_description.txt"
    output_logos = f"{output}/logos"
    output_mast = f"{output}/mast"

    executor_dir = ExecutorDirectory(
        file=file,
        log=log,
        trash=trash,
        input_seqs=input_seqs,
        single_seq=single_seq,
        input_pdb=input_pdb,
        fasta_for_pdb=fasta_for_pdb,
        p2_7_env=p2_7_env,
        dhcl_exec=dhcl_exec,
        bash_exec=bash_exec,
        meme_dir=meme_dir,
        output=output,
        output_clusters=output_clusters,
        output_logos=output_logos,
        output_mast=output_mast)

    cluster_dir = ClusterDirectory(
        file=file,
        log=log,
        trash=trash,
        meme_dir=meme_dir,
        bash_exec=bash_exec,
        input_mast=None,    # Necessary
        input_meme=None,    # Necessary
        output_mast=None,   # Optional
        description=None,   # Optional
        logos=None,         # Optional
        cluster_pkl=None)   # Optional

    filter_dir = FilterDirectory(
        file=file,
        log=log,
        trash=trash,
        meme_dir=meme_dir,
        input_seqs=input_seqs,
        single_seq=single_seq,
        bash_exec=bash_exec,
        orig=None,
        cleaned=None)
