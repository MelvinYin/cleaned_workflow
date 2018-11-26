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
    'ExecIntDir', 'dhcl_output consensus_seeds converge_seeds meme_full '
                  'starter_meme converge_output converge_pssm '
                  'meme_merged meme_cleaned converge_meme converge_composition')

class Executor:
    def __init__(self):
        self.switches = self.set_switches()
        self.dir = Directory.executor_dir
        self._dir = self.set_internal_dir()

    def set_internal_dir(self):
        _dir = ExecIntDir(
            converge_output=f"{self.dir.file}/converge_output",
            converge_seeds=f"{self.dir.file}/converge_seeds.fasta",
            dhcl_output=f"{self.dir.file}/from_dhcl",
            consensus_seeds=f"{self.dir.file}/init_seed_seqs.txt",
            meme_full=f"{self.dir.file}/meme_full",
            starter_meme=f"{self.dir.file}/meme_starter.txt",
            meme_merged=f"{self.dir.file}/meme_merged.txt",
            meme_cleaned=f"{self.dir.file}/meme_cleaned.txt",
            converge_meme=f"{self.dir.file}/converge_meme.txt",
            converge_pssm=f"{self.dir.file}/converge_pssm",
            converge_composition=f"{self.dir.file}/converge_composition")
        return _dir

    def set_switches(self):
        switches = OrderedDict()
        switches['MERGE_INPUT'] = (True, self.merge_input)
        switches['SHRINK_INPUT'] = (True, self.shrink_input)
        switches['RUN_DHCL'] = (True, self.run_dhcl)
        switches['EXTRACT_CONSENSUS'] = (True, self.extract_consensus)
        switches['REDUCE_CONSENSUS'] = (True, self.reduce_consensus)
        # Either run this (using meme)
        switches['BUILD_PSSM'] = (False, self.build_pssm)
        switches['BUILD_STARTER'] = (False, self.build_starter)
        switches['CLEAN_PSSM'] = (False, self.clean_pssm)
        switches['MERGE_PSSM'] = (False, self.merge_pssm)
        switches['SCREEN_PSSM'] = (False, self.screen_pssm)
        # Or this (using converge)
        switches['BUILD_CONVERGE_SEEDS'] = (True, self.build_converge_seeds)
        switches['RUN_CONVERGE'] = (True, self.run_converge)
        switches['TO_MEME_FORMAT'] = (True, self.to_meme_format)
        switches['screen_converge_pssm'] = (True, self.screen_converge_pssm)

        switches['ASSEMBLE_COMBI'] = (True, self.assemble_combi)
        switches['CLUSTER'] = (True, self.cluster)
        switches['DELETE_INTERMEDIATE'] = (True, self.delete_intermediate)
        return switches

    def merge_input(self):
        # Input: self.dir.input_seqs
        # Output: self.dir.input_seqs
        assert os.path.isdir(self.dir.input_seqdir)
        from create_input_seqs import main
        kwargs = dict(input_dir=self.dir.input_seqdir,
                      output=self.dir.input_seqs)
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
                      divideby=10,
                      output=self._dir.consensus_seeds)
        main(kwargs)
        assert os.path.isfile(self._dir.consensus_seeds)
        return

    # Converge
    def build_converge_seeds(self):
        # Input: self._dir.consensus_seeds
        # Output: self._dir.converge_seeds
        assert os.path.isfile(self._dir.consensus_seeds)
        from seeds_to_converge_input import main
        kwargs = dict(seed_seqs=self._dir.consensus_seeds,
                      output=self._dir.converge_seeds)
        main(kwargs)
        assert os.path.isfile(self._dir.converge_seeds)
        return

    def run_converge(self):
        # Input: self._dir.converge_seeds
        # Output: self._dir.converge_pssm | self._dir.converge_composition
        assert os.path.isfile(self._dir.converge_seeds)

        seed_name = self._dir.converge_seeds.rsplit("/", maxsplit=1)[-1]
        input_seq_name = self.dir.input_seqs.rsplit("/", maxsplit=1)[-1]
        composition_name = self.dir.converge_composition.rsplit("/", maxsplit=1)[-1]
        exec_name = self.dir.converge_exec.rsplit("/", maxsplit=1)[-1]

        shutil.copy(self._dir.converge_seeds, self.dir.converge_dir)
        shutil.copy(self.dir.input_seqs, self.dir.converge_dir)

        command = f"cd {self.dir.converge_dir} && mpirun -np 7 " \
                f"./{exec_name} -B -E 1 -r 1 -f 1 -c {composition_name} " \
            f"-i {seed_name} -p {input_seq_name}"
        subprocess.run(command, shell=True, executable=self.dir.bash_exec,
                       stdout=subprocess.DEVNULL)
        if os.path.isfile(self._dir.converge_pssm):
            os.remove(self._dir.converge_pssm)
        if os.path.isfile(self._dir.converge_composition):
            os.remove(self._dir.converge_composition)
        shutil.move(self.dir.converge_output, self._dir.converge_pssm)
        shutil.move(self.dir.converge_composition,
                    self._dir.converge_composition)
        self.to_trash(self.dir.converge_discard)
        self.to_trash(f"{self.dir.converge_dir}/{seed_name}")
        self.to_trash(f"{self.dir.converge_dir}/{input_seq_name}")
        assert os.path.isfile(self._dir.converge_pssm)
        assert os.path.isfile(self._dir.converge_composition)
        return

    def to_meme_format(self):
        assert os.path.isfile(self._dir.converge_pssm)
        assert os.path.isfile(self._dir.converge_composition)
        from converge2meme import main
        kwargs = dict(input_pssm=self._dir.converge_pssm,
                      composition=self._dir.converge_composition,
                      output=self._dir.converge_meme)
        main(kwargs)
        assert os.path.isfile(self._dir.converge_meme)
        return

    def screen_converge_pssm(self):
        assert os.path.isfile(self._dir.converge_meme)
        from filter import Filter
        filter_dir = Directory.filter_dir._replace(
            memefile=self._dir.converge_meme)
        conv_filter = Filter(filter_dir)
        conv_filter.run(to_run=[conv_filter.screen_correlated,
                                conv_filter.screen_non_combi])
        shutil.copy(self._dir.converge_meme, self._dir.meme_cleaned)
        return True

    # Meme
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
        shutil.copy(self._dir.meme_cleaned, self._dir.meme_merged)
        filter_dir = Directory.filter_dir._replace(
            memefile=self._dir.meme_cleaned)
        Filter(filter_dir).run()
        return True

    # Resume
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
        Cluster(cluster_dir).run()
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
        Filter(Directory.filter_dir).delete_intermediate()
        Cluster(Directory.cluster_dir).delete_intermediate()
        return

executor = Executor()
executor.run()