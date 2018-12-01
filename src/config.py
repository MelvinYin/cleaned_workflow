from collections import namedtuple

ExecutorDirectory = namedtuple(
    'ExecutorDirectory',
    "file log trash input_seqs input_seqdir input_pdb fasta_for_pdb p2_7_env "
    "meme_dir dhcl_exec bash_exec converge_dir converge_exec "
    "converge_composition converge_output converge_discard num_p output "
    "output_clusters output_logos output_mast seeds_divisor seq_divisor")

FilterDirectory = namedtuple(
    'FilterDirectory',
    "file log trash meme_dir input_seqs short_seq bash_exec memefile")

ClusterDirectory = namedtuple(
    "ClusterDirectory",
    "file log trash meme_dir bash_exec input_mast input_meme description "
    "logos cluster_pkl")

class Directory:
    file = "files"
    output = "output"
    p2_7_env = "/home/melvin/anaconda3/envs/p2.7/bin/python"
    dhcl_exec = "./external_scripts/dhcl/executables/everything.py"
    bash_exec = "/bin/bash"
    meme_dir = "./external_scripts/meme/bin"

    converge_dir = "external_scripts/pipeline"
    converge_exec = f"{converge_dir}/converge"
    converge_composition = f"{converge_dir}/composition.txt"
    converge_output = f"{converge_dir}/output.4.matrix.0"
    converge_discard = f"{converge_dir}/output.1.matrix.0"
    # Either input_seqs is provided, or it is created using input_seqdir.
    # A path should still be set for input_seqs here, because it is used
    # outside of Executor (Filter, etc) and is set here in config rather than
    # in Executor.
    # See merge_input in Executor
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
        memefile=None)
