import shutil
import subprocess
import os
from collections import OrderedDict, namedtuple
from utils import move_replace
from filter import Filter
import sys

from config import Dir, ClusterDir
from cluster import Cluster

ExecIntDir = namedtuple(
    'ExecIntDir', 'dhcl_output consensus_seeds meme_full starter_meme '
                  'meme_merged meme_cleaned')

class Executor:
    def __init__(self):
        self.switches = self.set_switches()
        self.dir = Dir.executor_dir
        self._dir = self.set_internal_dir()

    def set_internal_dir(self):
        _dir = ExecIntDir(
            dhcl_output="files/from_dhcl",
            consensus_seeds="files/init_seed_seqs.txt",
            meme_full="files/meme_full",
            starter_meme="files/meme_starter.txt",
            meme_merged="files/meme_merged.txt",
            meme_cleaned='files/meme_cleaned.txt')
        return _dir

    def set_switches(self):
        switches = OrderedDict()
        switches['SHRINK_INPUT'] = (True, self.shrink_input)
        switches['RUN_DHCL'] = (True, self.run_dhcl)
        switches['EXTRACT_CONSENSUS'] = (True, self.extract_consensus)
        switches['REDUCE_CONSENSUS'] = (True, self.reduce_consensus)
        switches['BUILD_PSSM'] = (True, self.build_pssm)
        switches['BUILD_STARTER'] = (True, self.build_starter)
        switches['CLEAN_PSSM'] = (True, self.clean_pssm)
        switches['MERGE_PSSM'] = (True, self.merge_pssm)
        switches['SCREEN_PSSM'] = (True, self.screen_pssm)
        switches['ASSEMBLE_COMBI'] = (True, self.assemble_combi)
        switches['CLUSTER'] = (True, self.cluster)
        switches['DELETE_INTERMEDIATE'] = (False, self.delete_intermediate)
        return switches

    def shrink_input(self):
        # Input: ./files/input_seqs.fasta
        # Output: ./files/input_seqs.fasta
        print("shrink_input:")
        assert os.path.isfile(self.dir.input_seqs)
        from shrink_input_for_test import main
        kwargs = dict(seqs=self.dir.input_seqs, output=self.dir.input_seqs,
                      divideby=20)
        main(kwargs)
        assert os.path.isfile(self.dir.input_seqs)
        print("Success!\n")
        return

    def run_dhcl(self):
        # Input: ./files/input_pdb
        # Output: ./files/from_dhcl
        print("RUN_DHCL:")
        assert os.path.isdir(self.dir.input_pdb)
        shutil.rmtree(self._dir.dhcl_output, ignore_errors=True)
        os.mkdir(self._dir.dhcl_output)
        command = f'{self.dir.p2_7_env} {self.dir.dhcl_exec} -d ' \
                  f'{self.dir.input_pdb} --outdir {self._dir.dhcl_output}'
        subprocess.run(command, shell=True)
        assert os.path.isdir(self._dir.dhcl_output)
        print("Success!\n")
        return

    def extract_consensus(self):
        # Input: ./files/from_dhcl    ./files/input_fasta
        # Output: ./files/init_seed_seqs.txt
        print("extract_consensus:")
        assert os.path.isdir(self._dir.dhcl_output)
        assert os.path.isdir(self.dir.fasta_for_pdb)
        from process_dhcl_output_meme import main
        kwargs = dict(dhcl_dir=self._dir.dhcl_output,
                      fasta_dir=self.dir.fasta_for_pdb,
                      output=self._dir.consensus_seeds)
        main(kwargs)
        assert os.path.isfile(self._dir.consensus_seeds)
        print("Success!\n")
        return

    def reduce_consensus(self):
        # Input: ./files/init_seed_seqs.txt
        # Output: ./files/init_seed_seqs.txt
        print("reduce_consensus:")
        assert os.path.isfile(self._dir.consensus_seeds)
        from reduce_dhcl_test import main
        kwargs = dict(input=self._dir.consensus_seeds,
                      denominator=8,
                      output=self._dir.consensus_seeds)
        main(kwargs)
        assert os.path.isfile(self._dir.consensus_seeds)
        print("Success!\n")
        return

    def build_pssm(self):
        # Input: ./files/init_seed_seqs.txt    ./files/input_seqs.fasta
        # Output: ./files/meme_full
        print("build_pssm:")
        assert os.path.isfile(self._dir.consensus_seeds)
        assert os.path.isfile(self.dir.input_seqs)
        shutil.rmtree(self._dir.meme_full, ignore_errors=True)
        os.mkdir(self._dir.meme_full)
        from run_meme_on_dhcl_loops import main
        kwargs = dict(consensus=self._dir.consensus_seeds,
                      seqs=self.dir.input_seqs,
                      output_folder=self._dir.meme_full)
        main(kwargs)
        assert os.path.isdir(self._dir.meme_full)
        print("Success!\n")
        return

    def build_starter(self):
        # Input: ./files/input_seqs.fasta
        # Output: ./files/meme_starter.txt
        print("build_starter:")
        assert os.path.isfile(self.dir.input_seqs)
        command = f'{self.dir.meme_dir}/meme {self.dir.input_seqs} ' \
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
            kwargs = dict(input=f"{self._dir.meme_full}/{filename}",
                          output=f"{self._dir.meme_full}/{filename}")
            main(kwargs)
        assert os.path.isdir(self._dir.meme_full)
        assert os.path.isfile(self._dir.starter_meme)
        print("Success!\n")
        return

    def merge_pssm(self):
        # todo: remove SUMMARY OF MOTIFS for smaller input_seqs
        # Input: ./files/meme_full    ./files/meme_starter.txt
        # Output: ./files/meme_merged.txt
        print("merge_pssm:")
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

    def screen_pssm(self):
        # Input: self._dir.meme_merged    self.dir.input_seqs.fasta
        #        ./files/single_seq.fasta
        # Output: ./files/mast_single     ./files/mast_nocorr
        # self._dir.meme_cleaned
        from config import FilterDir
        filter_dir = FilterDir(
            file=self.dir.file,
            input_seqs=self.dir.input_seqs,
            log=self.dir.log,
            trash=self.dir.trash,
            meme_dir=self.dir.meme_dir,
            bash_exec=self.dir.bash_exec,
            single_seq=self.dir.single_seq,
            cleaned=self._dir.meme_cleaned,
            orig=self._dir.meme_merged)
        filter = Filter(filter_dir)
        filter.run()
        filter.delete_intermediate()
        return True

    def assemble_combi(self):
        # Input: .files/meme_format3.txt
        # ./external_scripts/meme/input_seqs.fasta
        # Output: ./files/mast_onlycombi    ./output/mast
        print("assemble_combi:")
        assert os.path.isfile(self._dir.meme_cleaned)
        self.to_trash(self.dir.output_mast)
        if not os.path.isdir('output'):
            os.mkdir('output')
        command = f'{self.dir.meme_dir}/mast -remcorr ' \
                  f'{self._dir.meme_cleaned} {self.dir.input_seqs} -o ' \
                  f'{self.dir.output_mast}'
        subprocess.run(command, shell=True, executable=self.dir.bash_exec)
        assert os.path.isdir(self.dir.output_mast)
        assert os.path.isfile(f"{self.dir.output_mast}/mast.txt")
        print("Success!\n")
        return

    def cluster(self):
        cluster_dir = ClusterDir(
            file=self.dir.file,
            log=self.dir.log,
            trash=self.dir.trash,
            input_mast=f"{self.dir.output_mast}/mast.txt",
            input_meme=self._dir.meme_cleaned,
            output_mast=f"{self.dir.output}/mast.txt",
            description=self.dir.output_clusters,
            logos=f"{self.dir.output}/logos",
            meme_dir=self.dir.meme_dir,
            cluster_pkl=None,
            bash_exec=self.dir.bash_exec)
        cluster = Cluster(cluster_dir)
        cluster.run()
        cluster.delete_intermediate()
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

    def to_trash(self, file):
        return move_replace(file, self.dir.trash)

    def delete_intermediate(self):
        print("delete_intermediate:")
        for file in self._dir:
            self.to_trash(file)
        print("Success!\n")
        return

executor = Executor()
executor.run()