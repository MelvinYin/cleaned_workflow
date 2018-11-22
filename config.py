from collections import namedtuple

ExecutorDir = namedtuple(
    'ExecutorDir',
    'file log trash input_seqs input_pdb fasta_for_pdb p2_7_env meme_dir '
    'dhcl_exec bash_exec output single_seq output_clusters')

FilterDir = namedtuple(
    'FilterDir',
    'file log trash meme_dir bash_exec single_seq orig cleaned mast')


class Dir:
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
        output_clusters=output_clusters)

    filter_dir = FilterDir(
        file=file,
        log=log,
        trash=trash,
        meme_dir=meme_dir,
        bash_exec=bash_exec,
        single_seq=single_seq,
        cleaned=None,
        mast=None,
        orig=None)


# def set_dir(self):
#     dir = DirTemplate(
#         file="files/",
#         log="files/log.txt",
#         trash="files/_trash/",
#         input_seqs="files/input_seqs.fasta",
#         input_pdb="files/input_pdb_test",
#         fasta_for_pdb='files/input_fasta',
#         output="output/",
#         p2_7_env='/home/melvin/anaconda3/envs/p2.7/bin/python',
#         ext_meme='./external_scripts/meme/bin/',
#         ext_dhcl_exec='./external_scripts/dhcl/executables/everything.py',
#         bash_exec='/bin/bash',
#         single_seq='files/single_seq.fasta',
#         output_cluster_description='output/cluster_description.txt')
#     return dir

# DirTemplate = namedtuple(
#     'dir', 'file log trash input_seqs input_pdb fasta_for_pdb p2_7_env '
#            'ext_meme ext_dhcl_exec bash_exec output single_seq '
#            'output_cluster_description')

# DirTemplate = namedtuple(
#     'dir', 'file log trash meme_dir bash_exec single_seq orig cleaned mast')



# def set_dir(self):
#     dir = DirTemplate(
#         file="files/",
#         log="files/log.txt",
#         trash="files/_trash/",
#         meme_dir='./external_scripts/meme/bin/',
#         bash_exec='/bin/bash',
#         single_seq='files/single_seq.fasta',
#         orig="files/meme_merged.txt",
#         cleaned='files/meme_cleaned.txt',
#         mast='files/mast_onlycombi')
#     return dir