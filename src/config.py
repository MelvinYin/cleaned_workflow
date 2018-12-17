from collections import namedtuple

ExecutorDirectory = namedtuple(
    'ExecutorDirectory',
    "bash_exec converge_composition converge_dir converge_discard "
    "converge_exec converge_output dhcl_exec fasta_for_pdb file input_pdb "
    "input_seqdir input_seqs log meme_dir num_p output output_clusters "
    "output_logos output_mast p2_7_env trash seeds_divisor seq_divisor")

FilterDirectory = namedtuple(
    'FilterDirectory',
    "bash_exec file input_seqs log meme_dir memefile short_seq trash "
    "evalue_threshold entropy_bits_threshold")

ClusterDirectory = namedtuple(
    "ClusterDirectory",
    "bash_exec cluster_pkl description input_mast input_meme file log logos "
    "meme_dir trash")

class Directory:
    file = "files"
    output = "output"
    p2_7_env = "/home/melvin/anaconda3/envs/dhcl_p/bin/python"
    dhcl_exec = "external_scripts/dhcl/executables/everything.py"
    bash_exec = "/bin/bash"
    meme_dir = "external_scripts/meme/bin"
    converge_dir = "external_scripts/pipeline"
    converge_exec = f"{converge_dir}/converge"
    converge_composition = f"{converge_dir}/composition.txt"
    converge_output = f"{converge_dir}/output.4.matrix.0"
    converge_discard = f"{converge_dir}/output.1.matrix.0"
    input_seqdir = f"{file}/sfld_datasets"
    input_seqs = f"{file}/input_seqs.fasta"
    log = f"{file}/log.txt"
    trash = f"{file}/_trash"
    input_pdb = f"{file}/input_pdb"
    fasta_for_pdb = f"{file}/input_fasta"
    output_clusters = f"{output}/cluster_description.txt"
    output_logos = f"{output}/logos"
    output_mast = f"{output}/mast"
    num_processor = 7    # Memory use increases as well
    seeds_divisor = 10
    seq_divisor = 10
    evalue_threshold = 0.5
    entropy_bits_threshold = 30

    executor_dir = ExecutorDirectory(
        file=file,
        log=log,
        trash=trash,
        input_seqs=input_seqs,
        input_seqdir=input_seqdir,
        input_pdb=input_pdb,
        fasta_for_pdb=fasta_for_pdb,
        p2_7_env=p2_7_env,
        dhcl_exec=dhcl_exec,
        bash_exec=bash_exec,
        converge_exec=converge_exec,
        converge_dir=converge_dir,
        converge_composition=converge_composition,
        converge_output=converge_output,
        converge_discard=converge_discard,
        num_p=num_processor,
        meme_dir=meme_dir,
        output=output,
        output_clusters=output_clusters,
        output_logos=output_logos,
        output_mast=output_mast,
        seeds_divisor=seeds_divisor,
        seq_divisor=seq_divisor)

    cluster_dir = ClusterDirectory(
        file=file,
        log=log,
        trash=trash,
        meme_dir=meme_dir,
        bash_exec=bash_exec,
        input_mast=None,    # Necessary
        input_meme=None,    # Necessary
        description=None,   # Optional
        logos=None,         # Optional
        cluster_pkl=None)   # Optional

    filter_dir = FilterDirectory(
        file=file,
        log=log,
        trash=trash,
        meme_dir=meme_dir,
        input_seqs=input_seqs,
        bash_exec=bash_exec,
        short_seq=None,
        memefile=None,
        evalue_threshold=evalue_threshold,
        entropy_bits_threshold=entropy_bits_threshold)
