import shutil
import subprocess
import os
from collections import OrderedDict, namedtuple
from utils import move_replace
from filter import Filter

from config import Dir

IntDirTemplate = namedtuple(
    '_dir', 'dhcl_output consensus_seeds meme_full starter_meme meme_merged '
            'meme_cleaned noncombi_clusterparam noncombi_cluster output_mast '
            'cluster_description cluster_centroid')

class Executor:
    def __init__(self):
        self.switches = self.set_switches()
        self.dir = Dir.executor_dir
        self._dir = self.set_internal_dir()
        
    def screen_pssm(self):
        filter = Filter()
        filter.dir.orig = self._dir.meme_merged
        filter.dir.cleaned = self._dir.meme_cleaned
        filter.dir.mast = self._dir.output_mast
        filter.run()
        return True

    def run(self):
        for _switch, (to_run, function) in self.switches.items():
            if to_run:
                try:
                    function()
                except:
                    print("Error: {}".format(function.__name__))
                    raise
        return

    def set_internal_dir(self):
        # 'meme_evalue_screened meme_entropy_screened '
        # 'meme_correlated_screened meme_cleaned meme_cleaned')
        _dir = IntDirTemplate(
            dhcl_output="files/from_dhcl/",
            consensus_seeds="files/init_seed_seqs.txt",
            meme_full="files/meme_full",
            starter_meme="files/meme_starter.txt",
            meme_merged="files/meme_merged.txt",
            meme_cleaned='files/meme_cleaned.txt',
            noncombi_clusterparam='files/clustering_df.pkl',
            noncombi_cluster='files/cluster_final.pkl',
            output_mast='output/mast',
            cluster_centroid="./files/cluster_centroids.pkl",
            cluster_description='files/cluster_description.txt')
        return _dir

    def set_switches(self):
        switches = OrderedDict()
        switches['SHRINK_INPUT'] = (False, self.shrink_input)
        switches['RUN_DHCL'] = (False, self.run_dhcl)
        switches['EXTRACT_CONSENSUS'] = (False, self.extract_consensus)
        switches['REDUCE_DHCL'] = (False, self.reduce_dhcl)
        switches['BUILD_PSSM'] = (False, self.build_pssm)
        switches['BUILD_STARTER'] = (False, self.build_starter)
        switches['CLEAN_PSSM'] = (False, self.clean_pssm)
        switches['MERGE_PSSM'] = (False, self.merge_pssm)
        switches['SCREEN_PSSM'] = (False, self.screen_pssm)
        switches['ASSEMBLE_COMBI'] = (False, self.assemble_combi)
        switches['GET_CLUSTER_PARAMS'] = (False, self.get_cluster_params)
        switches['CLUSTER_COMBI'] = (False, self.cluster_combi)
        switches['CREATE_CLUSTER_MOTIFS'] = (False, self.create_cluster_motifs)
        switches['CREATE_CLUSTER_LOGOS'] = (False, self.create_cluster_logos)
        switches['DELETE_INTERMEDIATE'] = (True, self.delete_intermediate)
        return switches

    def to_trash(self, file):
        return move_replace(file, self.dir.trash)

    def run_dhcl(self):
        # Input: ./files/input_pdb
        # Output: ./files/from_dhcl
        print("RUN_DHCL:")
        assert os.path.isdir(self.dir.input_pdb)
        shutil.rmtree("./files/from_dhcl", ignore_errors=True)
        os.mkdir("./files/from_dhcl")
        command = f'{self.dir.p2_7_env} {self.dir.ext_dhcl_exec} -d ' \
                  f'{self.dir.input_pdb} --outdir '
        subprocess.run(command, shell=True)
        assert os.path.isdir('files/from_dhcl')
        print("Success!\n")
        return

    def extract_consensus(self):
        # Input: ./files/from_dhcl    ./files/input_fasta
        # Output: ./files/init_seed_seqs.txt
        print("extract_consensus:")
        assert os.path.isdir('files/from_dhcl')
        assert os.path.isdir('files/input_fasta')
        from process_dhcl_output_meme import main
        kwargs = dict(dhcl_dir='files/from_dhcl',
                      fasta_dir='files/input_fasta',
                      output='files/init_seed_seqs.txt')
        main(kwargs)
        assert os.path.isfile('files/init_seed_seqs.txt')
        print("Success!\n")
        return

    def build_pssm(self):
        # Input: ./files/init_seed_seqs.txt    ./files/input_seqs.fasta
        # Output: ./files/meme_full
        print("run_meme:")
        assert os.path.isfile('files/init_seed_seqs.txt')
        assert os.path.isfile(self.dir.input_seqs)
        shutil.rmtree("./files/meme_full", ignore_errors=True)
        os.mkdir("./files/meme_full")
        from run_meme_on_dhcl_loops import main
        kwargs = dict(consensus='./files/init_seed_seqs.txt',
                      seqs=self.dir.input_seqs,
                      output_folder='./files/meme_full')
        main(kwargs)
        assert os.path.isdir('files/meme_full')
        print("Success!\n")
        return

    def build_starter(self):
        # Input: ./files/input_seqs.fasta
        # Output: ./files/meme_starter.txt
        print("build_starter:")
        assert os.path.isfile('files/input_seqs.fasta')
        command = f'{self.dir.ext_meme}meme {self.dir.input_seqs} ' \
                  f'-text -protein -w 30 -p 7 -nmotifs 2 -nostatus &>> ' \
                  f'{self._dir.starter_meme}'
        subprocess.run(command, shell=True, executable=self.dir.bash_exec)
        assert os.path.isfile(self._dir.starter_meme)
        print("Success!\n")
        return

    def clean_pssm(self):
        # Input: ./files/meme_full    ./files/meme_starter.txt
        # Output: ./files/meme_full    ./files/meme_starter.txt
        print("clean_meme:")
        assert os.path.isdir(self._dir.meme_full)
        assert os.path.isfile(self._dir.starter_meme)
        from meme_cleaner import main
        kwargs = dict(input=self._dir.starter_meme,
                      output=self._dir.starter_meme)
        main(kwargs)
        for filename in os.listdir(self._dir.meme_full):
            kwargs = dict(input=self._dir.meme_full+filename,
                          output=self._dir.meme_full+filename)
            main(kwargs)
        assert os.path.isdir(self._dir.meme_full)
        assert os.path.isfile(self._dir.starter_meme)
        print("Success!\n")
        return

    def merge_pssm(self):
        # todo: remove SUMMARY OF MOTIFS for smaller input_seqs
        # Input: ./files/meme_full    ./files/meme_starter.txt
        # Output: ./files/meme_merged.txt
        print("merge_meme:")
        assert os.path.isdir(self._dir.meme_full)
        assert os.path.isfile(self._dir.starter_meme)
        from meme_merger import main
        kwargs = dict(meme_folder=self._dir.meme_full,
                      output=self._dir.meme_merged,
                      meme_starter=self._dir.starter_meme)
        main(kwargs)
        assert os.path.isfile(self._dir.meme_merged)
        print("Success!\n")
        return
    #
    # def screen_evalue(self):
    #     # Input: ./files/meme_merged.txt
    #     # Output: ./files/meme_evalue_screened.txt
    #     print("screen_evalue:")
    #     assert os.path.isfile(self._dir.meme_merged)
    #     from remove_motifs_with_low_evalue import main
    #     kwargs = dict(meme=self._dir.meme_merged,
    #                   output=self._dir.meme_evalue_screened)
    #     main(kwargs)
    #     assert os.path.isfile(self._dir.meme_evalue_screened)
    #     print("Success!\n")
    #     return
    #
    # def screen_entropy(self):
    #     # Input: ./files/meme_evalue_screened.txt
    #     # Output: ./files/meme_format.txt
    #     print("screen_entropy:")
    #     assert os.path.isfile(self._dir.meme_evalue_screened)
    #     from screen_motif_entropy import main
    #     kwargs = dict(meme=self._dir.meme_evalue_screened,
    #                   output=self._dir.meme_entropy_screened)
    #     main(kwargs)
    #     assert os.path.isfile(self._dir.meme_entropy_screened)
    #     print("Success!\n")
    #     return
    #
    # def screen_correlated(self):
    #     # Input: ./files/meme_format.txt    ./files/single_seq.fasta
    #     # Output: ./files/mast_single    ./files/meme_format2.txt
    #     print("screen_correlated:")
    #     assert os.path.isfile(self._dir.meme_entropy_screened)
    #     assert os.path.isfile(self.dir.single_seq)
    #     command = f'{self.dir.ext_meme}mast -remcorr ' \
    #               f'{self._dir.meme_entropy_screened} {self.dir.single_seq}'
    #     subprocess.run(command, shell=True, executable=self.dir.bash_exec)
    #     shutil.rmtree("", ignore_errors=True)
    #     os.mkdir(self._dir.mast_single)
    #     self.move_replace('mast_out', self.dir.file, 'mast_single')
    #     from mast_remove_profiles import main
    #     kwargs = dict(meme_in=self._dir.meme_entropy_screened,
    #                   meme_out=self._dir.meme_correlated_screened)
    #     main(kwargs)
    #     assert os.path.isdir(self._dir.mast_single)
    #     assert os.path.isfile(self._dir.meme_correlated_screened)
    #     print("Success!\n")
    #     return
    #
    # def screen_non_combi(self):
    #     # todo: split this up
    #     # Input: ./files/meme_format2.txt    ./files/single_seq.fasta
    #     # Output: ./files/mast_nocorr    ./files/meme_format3.txt
    #     print("screen_non_combi:")
    #     assert os.path.isfile(self._dir.meme_correlated_screened)
    #     assert os.path.isfile(self.dir.single_seq)
    #     command = f"{self.dir.ext_meme}mast -remcorr " \
    #               f"{self._dir.meme_correlated_screened} {self.dir.input_seqs}"
    #     subprocess.run(command, shell=True, executable=self.dir.bash_exec)
    #     shutil.rmtree(self._dir.mast_nocorr, ignore_errors=True)
    #     os.mkdir(self._dir.mast_nocorr)
    #     self.move_replace('mast_out', self.dir.file, 'mast_nocorr')
    #
    #     from generate_cluster_params import main
    #     kwargs = dict(input_fname=self._dir.mast_nocorr+"mast.txt",
    #                   screen_threshold=5,
    #                   pkl_path=self._dir.noncombi_clusterparam)
    #     main(kwargs)
    #
    #     from cluster_final import main
    #     kwargs = dict(cluster_threshold=50,
    #                   pkl_path=self._dir.noncombi_clusterparam,
    #                   output=None,
    #                   cluster_df_pkl="./files/cluster_final.pkl")
    #     main(kwargs)
    #
    #     from mast_remove_profiles_using_pkl import main
    #     kwargs = dict(cluster_df_pkl="./files/cluster_final.pkl",
    #                   input=self._dir.meme_correlated_screened,
    #                   output=self._dir.meme_cleaned)
    #     main(kwargs)
    #     self.to_trash(self._dir.meme_correlated_screened)
    #     self.to_trash(self._dir.noncombi_cluster)
    #     assert os.path.isfile(self._dir.meme_cleaned)
    #     assert os.path.isdir(self._dir.mast_nocorr)
    #     print("Success!\n")
    #     return

    def assemble_combi(self):
        # Input: .iles/meme_format3.txt    ./external_scripts/meme/input_seqs.fasta
        # Output: ./files/mast_onlycombi    ./output/mast
        print("mast_combi:")
        assert os.path.isfile(self._dir.meme_cleaned)
        command = f'{self.dir.ext_meme}mast -remcorr ' \
                  f'{self._dir.meme_cleaned} {self.dir.input_seqs}'
        subprocess.run(command, shell=True, executable=self.dir.bash_exec)
        self.to_trash(self._dir.mast_onlycombi)
        os.mkdir(self._dir.mast_onlycombi)
        self.move_replace('mast_out', self.dir.file, 'mast_onlycombi')
        shutil.copytree(self._dir.mast_onlycombi, self._dir.output_mast)
        assert os.path.isdir(self._dir.mast_onlycombi)
        assert os.path.isdir(self._dir.output_mast)
        print("Success!\n")
        return

    def get_cluster_params(self):
        # Input: ./files/mast_onlycombi/mast.txt
        # Output: ./files/clustering_df.pkl
        print("get_cluster_params:")
        assert os.path.isfile(self._dir.mast_onlycombi+'mast.txt')
        from generate_cluster_params import main
        kwargs = dict(input_fname=self._dir.mast_onlycombi+'mast.txt',
                      screen_threshold=5,
                      pkl_path=self._dir.noncombi_clusterparam)
        main(kwargs)
        assert os.path.isfile(self._dir.noncombi_clusterparam)
        print("Success!\n")
        return

    def cluster_combi(self):
        # Input: ./files/clustering_df.pkl
        # Output: output/cluster_description.txt
        #         (optional) files/cluster_centroids.pkl
        print("cluster_combi:")
        assert os.path.isfile(self._dir.noncombi_clusterparam)
        from cluster_final import main
        kwargs = dict(cluster_threshold=50,
                      pkl_path=self._dir.noncombi_clusterparam,
                      output=self._dir.cluster_description,
                      cluster_df_pkl=self._dir.cluster_centroid)
        main(kwargs)
        shutil.move(self._dir.cluster_description,
                    self.dir.output_cluster_description)
        assert os.path.isfile(self.dir.output_cluster_description)
        print("Success!\n")
        return

    def create_cluster_motifs(self):
        # Input: ./files/meme_format3.txt    ./files/cluster_centroids.pkl
        # Output: multiple ./files/motifs/motifs_in_cluster_{}.txt
        print("create_cluster_motifs:")
        assert os.path.isfile(self._dir.meme_cleaned)
        assert os.path.isfile(self._dir.cluster_centroid)
        from split_motifs_into_cluster_motifs import main
        self.to_trash('files/motifs')
        os.mkdir("files/motifs")
        kwargs = dict()
        main(kwargs)
        assert os.path.isdir('files/motifs')
        print("Success!\n")
        return

    def create_cluster_logos(self):
        # todo: check ordering of logo after renaming, see if it is a
        # file_folder options thing or my fault.

        # Input: multiple ./files/motifs/motifs_in_cluster_{}.txt
        # Output: multiple ./output/logos/cluster_{}/logo_{}.png
        print("create_cluster_logos:")
        assert os.path.isdir('files/motifs')
        self.to_trash('output/logos')
        os.mkdir("output/logos")
        for filename in os.listdir('files/motifs'):
            cluster_i = filename[18:-4]
            i=1
            os.mkdir("output/logos/cluster_{}".format(cluster_i))
            while True:
                command = f'{self.dir.ext_meme}ceqlogo -i{i} ' \
                      f'files/motifs/{filename} -o ' \
                      f'output/logos/cluster_{cluster_i}/logo_{i}.png -f PNG'
                i += 1
                returncode = subprocess.run(
                    command, shell=True, executable=self.dir.bash_exec,
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)\
                    .returncode
                if returncode != 0:
                    break
        from rename_cluster_logos import main
        kwargs = dict(motif_filedir="files/motifs",
                      output_logodir="output/logos")
        main(kwargs)
        assert os.path.isdir('output/logos')
        print("Success!\n")
        return

    def delete_intermediate(self):
        print("delete_intermediate:")
        self.to_trash('./files/init_seed_seqs.txt')
        self.to_trash('./files/meme_format.txt')
        self.to_trash('./files/mast_single')
        self.to_trash('./files/meme_merged.txt')
        self.to_trash('./files/meme_starter.txt')
        self.to_trash('./files/meme_format2.txt')
        self.to_trash('./files/meme_format3.txt')
        self.to_trash('./files/mast_onlycombi')
        self.to_trash('./files/cluster_centroids.pkl')
        self.to_trash('./files/clustering_df.pkl')
        self.to_trash('./files/meme_evalue_screened.txt')
        self.to_trash('./files/motifs')
        self.to_trash('./files/mast_nocorr')
        self.to_trash('./files/meme_full')
        self.to_trash('files/from_dhcl')
        self.to_trash(self.dir.log)
        print("Success!\n")
        return

    def reduce_dhcl(self):
        # Input: ./files/init_seed_seqs.txt
        # Output: ./files/init_seed_seqs.txt
        print("reduce_dhcl:")
        assert os.path.isfile('files/init_seed_seqs.txt')
        from reduce_dhcl_test import main
        kwargs = dict(input='./files/init_seed_seqs.txt',
                      denominator=5,
                      output='./files/init_seed_seqs.txt')
        main(kwargs)
        assert os.path.isfile('files/init_seed_seqs.txt')
        print("Success!\n")
        return

    def shrink_input(self):
        # Input: ./files/input_seqs.fasta
        # Output: ./files/input_seqs.fasta
        print("shrink_input:")
        assert os.path.isfile('files/input_seqs.fasta')
        from shrink_input_for_test import main
        kwargs = dict(seqs='./files/input_seqs.fasta',
                      output='./files/input_seqs.fasta',
                      divideby=20)
        main(kwargs)
        assert os.path.isfile('files/input_seqs.fasta')
        print("Success!\n")
        return


executor = Executor()
executor.run()