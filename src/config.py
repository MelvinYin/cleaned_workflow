from collections import namedtuple

ExecutorDirectory = namedtuple(
    'ExecutorDirectory',
    "bash_exec converge_composition converge_dir converge_discard "
    "converge_exec converge_output dhcl_exec fasta_for_pdb file input_pdb "
    "input_seqdir input_seqs log meme_dir num_p output output_clusters "
    "output_logos output_mast p2_7_env trash seeds_divisor seq_divisor "
    "num_cluster_final")

FilterDirectory = namedtuple(
    'FilterDirectory',
    "bash_exec file input_seqs log meme_dir memefile short_seq trash "
    "evalue_ceiling entropy_threshold kldiv_threshold num_cluster_screening")

ClusterDirectory = namedtuple(
    "ClusterDirectory",
    "bash_exec cluster_pkl combi_minsize cluster_minsize description "
    "input_mast input_meme num_cluster output file log logos "
    "meme_dir trash")

class Directory:
    bash_exec = "/bin/bash"
    converge_dir = "external_scripts/pipeline"
    dhcl_exec = "external_scripts/dhcl/executables/everything.py"
    file = "files"
    meme_dir = "external_scripts/meme/bin"
    output = "output"
    p2_7_env = "/home/melvin/anaconda3/envs/dhcl_p/bin/python"
    converge_composition = f"{converge_dir}/composition.txt"
    converge_discard = f"{converge_dir}/output.1.matrix.0"
    converge_exec = f"{converge_dir}/converge"
    converge_output = f"{converge_dir}/output.4.matrix.0"
    fasta_for_pdb = f"{file}/input_fasta"
    input_pdb = f"{file}/input_pdb"
    input_seqdir = f"{file}/sfld_datasets"
    input_seqs = f"{file}/input_seqs.fasta"
    log = f"{file}/log.txt"
    output_clusters = f"{output}/cluster_description.txt"
    output_logos = f"{output}/logos"
    output_mast = f"{output}/mast"
    trash = f"{file}/_trash"
    # KLdiv high for highly conserved, this removes flat signatures
    kldiv_threshold = 37
    # Entropy low for highly conserved, this removes all-conserved signatures
    entropy_threshold = 17
    evalue_ceiling = 0.5
    num_processor = 7
    seeds_divisor = 10
    seq_divisor = 10
    cluster_minsize = 70
    combi_minsize = 10
    num_cluster_screening = 12
    num_cluster_final = 10

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
        seq_divisor=seq_divisor,
        num_cluster_final=num_cluster_final)

    cluster_dir = ClusterDirectory(
        bash_exec=bash_exec,
        file=file,
        log=log,
        meme_dir=meme_dir,
        trash=trash,
        cluster_minsize=cluster_minsize,
        combi_minsize=combi_minsize,
        output=output,
        num_cluster=None,   # Necessary
        input_mast=None,    # Necessary
        input_meme=None,    # Necessary
        description=None,   # Optional
        logos=None,         # Optional
        cluster_pkl=None)   # Optional

    filter_dir = FilterDirectory(
        bash_exec=bash_exec,
        entropy_threshold=entropy_threshold,
        kldiv_threshold=kldiv_threshold,
        evalue_ceiling=evalue_ceiling,
        file=file,
        log=log,
        meme_dir=meme_dir,
        trash=trash,
        num_cluster_screening=num_cluster_screening,
        input_seqs=None,
        memefile=None,
        short_seq=None)
