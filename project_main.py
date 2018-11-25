from collections import OrderedDict, namedtuple
import os
import shutil
import subprocess
import sys

from cluster import Cluster
from config import Directory
from filter import Filter
from utils import move_replace

ExecIntDir = namedtuple(
    'ExecIntDir', 'dhcl_output consensus_seeds meme_full starter_meme '
                  'meme_merged meme_cleaned')

class Executor:
    def __init__(self):
        self.switches = self.set_switches()
        self.dir = Directory.executor_dir
        self._dir = self.set_internal_dir()

    def set_internal_dir(self):
        _dir = ExecIntDir(
            dhcl_output=f"{self.dir.file}/from_dhcl",
            consensus_seeds=f"{self.dir.file}/init_seed_seqs.txt",
            meme_full=f"{self.dir.file}/meme_full",
            starter_meme=f"{self.dir.file}/meme_starter.txt",
            meme_merged=f"{self.dir.file}/meme_merged.txt",
            meme_cleaned=f"{self.dir.file}/meme_cleaned.txt")
        return _dir

    def set_switches(self):
        switches = OrderedDict()
        switches['MERGE_INPUT'] = (False, self.merge_input)
        switches['SHRINK_INPUT'] = (False, self.shrink_input)
        switches['RUN_DHCL'] = (False, self.run_dhcl)
        switches['EXTRACT_CONSENSUS'] = (False, self.extract_consensus)
        switches['REDUCE_CONSENSUS'] = (False, self.reduce_consensus)
        switches['BUILD_PSSM'] = (False, self.build_pssm)
        switches['BUILD_STARTER'] = (False, self.build_starter)
        switches['CLEAN_PSSM'] = (False, self.clean_pssm)
        switches['MERGE_PSSM'] = (False, self.merge_pssm)
        switches['SCREEN_PSSM'] = (False, self.screen_pssm)
        switches['ASSEMBLE_COMBI'] = (False, self.assemble_combi)
        switches['CLUSTER'] = (False, self.cluster)
        switches['DELETE_INTERMEDIATE'] = (True, self.delete_intermediate)
        return switches

    def merge_input(self):
        # Input: self.dir.input_seqs
        # Output: self.dir.input_seqs
        assert os.path.isdir(self.dir.input_seqdir)
        from create_input_seqs import main
        kwargs = dict(input_dir=self.dir.input_seqdir,
                      output=self._dir.input_seqs)
        main(kwargs)
        assert os.path.isfile(self.dir.input_seqs)
        return

    def shrink_input(self):
        # Input: self.dir.input_seqs
        # Output: self.dir.input_seqs
        assert os.path.isfile(self.dir.input_seqs)
        from shrink_input_for_test import main
        kwargs = dict(seqs=self.dir.input_seqs, output=self.dir.input_seqs,
                      divideby=10)
        main(kwargs)
        assert os.path.isfile(self.dir.input_seqs)
        return

    def run_dhcl(self):
        # Input: self.dir.input_pdb
        # Output: self._dir.dhcl_output
        assert os.path.isdir(self.dir.input_pdb)
        shutil.rmtree(self._dir.dhcl_output, ignore_errors=True)
        os.mkdir(self._dir.dhcl_output)
        command = f'{self.dir.p2_7_env} {self.dir.dhcl_exec} -d ' \
                  f'{self.dir.input_pdb} --outdir {self._dir.dhcl_output}'
        subprocess.run(command, shell=True)
        assert os.path.isdir(self._dir.dhcl_output)
        return

    def extract_consensus(self):
        # Input: self._dir.dhcl_output | self.dir.fasta_for_pdb
        # Output: self._dir.consensus_seeds
        assert os.path.isdir(self._dir.dhcl_output)
        assert os.path.isdir(self.dir.fasta_for_pdb)
        from dhcl_output_to_consensus import main
        kwargs = dict(dhcl_dir=self._dir.dhcl_output,
                      fasta_dir=self.dir.fasta_for_pdb,
                      output=self._dir.consensus_seeds)
        main(kwargs)
        assert os.path.isfile(self._dir.consensus_seeds)
        return

    def reduce_consensus(self):
        # Input: self._dir.consensus_seeds
        # Output: self._dir.consensus_seeds
        assert os.path.isfile(self._dir.consensus_seeds)
        from reduce_dhcl_for_test import main
        kwargs = dict(input=self._dir.consensus_seeds,
                      divideby=5,
                      output=self._dir.consensus_seeds)
        main(kwargs)
        assert os.path.isfile(self._dir.consensus_seeds)
        return

    def build_pssm(self):
        # Input: self._dir.consensus_seeds | self.dir.input_seqs
        # Output: self._dir.meme_full
        assert os.path.isfile(self._dir.consensus_seeds)
        assert os.path.isfile(self.dir.input_seqs)
        shutil.rmtree(self._dir.meme_full, ignore_errors=True)
        os.mkdir(self._dir.meme_full)
        from run_meme_on_cons import main
        kwargs = dict(consensus=self._dir.consensus_seeds,
                      seqs=self.dir.input_seqs,
                      output_folder=self._dir.meme_full,
                      num_p=self.dir.num_p,
                      meme_dir=self.dir.meme_dir)
        main(kwargs)
        assert os.path.isdir(self._dir.meme_full)
        return

    def build_starter(self):
        # Input: self.dir.input_seqs
        # Output: self._dir.starter_meme
        assert os.path.isfile(self.dir.input_seqs)
        command = f'{self.dir.meme_dir}/meme {self.dir.input_seqs} ' \
                  f'-text -protein -w 30 -p {self.dir.num_p} -nmotifs 2 ' \
                  f'-nostatus &>> {self._dir.starter_meme}'
        subprocess.run(command, shell=True, executable=self.dir.bash_exec)
        assert os.path.isfile(self._dir.starter_meme)
        return

    def clean_pssm(self):
        # Input: self._dir.meme_full | self._dir.starter_meme
        # Output: self._dir.meme_full | self._dir.starter_meme
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
        return

    def merge_pssm(self):
        # Input: self._dir.meme_full | self._dir.starter_meme
        # Output: self._dir.meme_merged
        assert os.path.isdir(self._dir.meme_full)
        assert os.path.isfile(self._dir.starter_meme)
        from merge_meme_files import main
        kwargs = dict(meme_folder=self._dir.meme_full,
                      meme_starter=self._dir.starter_meme,
                      output=self._dir.meme_merged)
        main(kwargs)
        assert os.path.isfile(self._dir.meme_merged)
        return

    def screen_pssm(self):
        # Input: self._dir.meme_merged | self.dir.input_seqs |
        #        self.dir.single_seq
        # Output: self._dir.meme_cleaned
        filter_dir = Directory.filter_dir._replace(
            orig=self._dir.meme_merged,
            cleaned=self._dir.meme_cleaned)
        filter = Filter(filter_dir)
        filter.run()
        filter.delete_intermediate()
        return True

    def assemble_combi(self):
        # Input: self._dir.meme_cleaned | self.dir.input_seqs
        # Output: self.dir.output_mast
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
        return

    def cluster(self):
        # Input:self.dir.output_mast | self._dir.meme_cleaned
        # Output: self.dir.output_clusters | self.dir.output_logos
        cluster_dir = Directory.cluster_dir._replace(
            input_mast=f"{self.dir.output_mast}/mast.txt",
            input_meme=self._dir.meme_cleaned,
            description=self.dir.output_clusters,
            logos=self.dir.output_logos)
        cluster = Cluster(cluster_dir)
        cluster.run()
        cluster.delete_intermediate()
        return True

    def run(self):
        if not os.path.isdir(self.dir.trash):
            os.mkdir(self.dir.trash)
        if not os.path.isdir(self.dir.output):
            os.mkdir(self.dir.output)
        for (to_run, func) in self.switches.values():
            if to_run:
                try:
                    print(f"{func.__name__}:")
                    func()
                    print("Success!\n")
                except:
                    print("Error: {}".format(func.__name__))
                    raise
        return

    def to_trash(self, file):
        return move_replace(file, self.dir.trash)

    def delete_intermediate(self):
        for file in self._dir:
            self.to_trash(file)
        return

executor = Executor()
executor.run()