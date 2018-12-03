from collections import OrderedDict, namedtuple
import os
import shutil
import subprocess
import sys
import numpy as np

# Add src folders to sys.path so import is straightforward
cwd = os.getcwd()
src = f"{cwd}/src"
cluster = f"{cwd}/src/cluster"
filter_ = f"{cwd}/src/filter"
sys.path.append(src)
sys.path.append(cluster)
sys.path.append(filter_)

from cluster import Cluster
from config import Directory
from filter import Filter
from utils import move_replace

ExecIntDir = namedtuple(
    'ExecIntDir', 'dhcl_output consensus_seeds converge_seeds meme_full '
                  'starter_meme converge_output converge_pssm pssm '
                  'meme_merged meme_cleaned converge_meme '
                  'converge_composition short_seq short_seq_len')

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
            converge_pssm=f"{self.dir.file}/converge_pssm.txt",
            pssm=f"{self.dir.file}/pssm.txt",
            converge_composition=f"{self.dir.file}/converge_composition.txt",
            short_seq=f"{self.dir.file}/short_seq.fasta",
            short_seq_len=3)
        return _dir

    def set_switches(self):
        switches = OrderedDict()
        # Get input_seqs
        switches['MERGE_INPUT'] = (False, self.merge_input)
        switches['SHRINK_INPUT'] = (False, self.shrink_input)
        switches['CREATE_SHORT_SEQS'] = (False, self.create_short_seqs)
        # Get consensus_loops
        switches['RUN_DHCL'] = (False, self.run_dhcl)
        switches['EXTRACT_CONSENSUS'] = (False, self.extract_consensus)
        switches['REDUCE_CONSENSUS'] = (False, self.reduce_consensus)
        # Get PSSM using:
        # Meme
        switches['BUILD_PSSM'] = (False, self.build_pssm)
        switches['BUILD_STARTER'] = (False, self.build_starter)
        switches['CLEAN_PSSM'] = (False, self.clean_pssm)
        switches['MERGE_PSSM'] = (False, self.merge_pssm)
        switches['SCREEN_PSSM'] = (False, self.screen_pssm)
        # Or converge
        switches['BUILD_CONVERGE_SEEDS'] = (True, self.build_converge_seeds)
        switches['RUN_CONVERGE'] = (True, self.run_converge)
        switches['TO_MEME_FORMAT'] = (True, self.to_meme_format)
        switches['SCREEN_CONVERGE_PSSM'] = (True, self.screen_converge_pssm)
        # Get combi
        switches['ASSEMBLE_COMBI'] = (True, self.assemble_combi)
        switches['CLUSTER'] = (True, self.cluster)
        switches['DELETE_INTERMEDIATE'] = (True, self.delete_intermediate)
        return switches

    def merge_input(self):
        # Input: self.dir.input_seqs
        # Output: self.dir.input_seqs
        if self.dir.input_seqdir:
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
        if self.dir.seq_divisor:
            from shrink_input_for_test import main
            kwargs = dict(seqs=self.dir.input_seqs, output=self.dir.input_seqs,
                          divisor=self.dir.seq_divisor)
            main(kwargs)
        assert os.path.isfile(self.dir.input_seqs)
        return

    def create_short_seqs(self):
        # Input: self.dir.input_seqs
        # Output: self._dir.short_seq
        assert os.path.isfile(self.dir.input_seqs)
        from create_short_seqs import main
        kwargs = dict(input=self.dir.input_seqs,
                      output=self._dir.short_seq,
                      length=self._dir.short_seq_len)
        main(kwargs)
        assert os.path.isfile(self._dir.short_seq)
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
        if self.dir.seeds_divisor:
            from reduce_dhcl_for_test import main
            kwargs = dict(input=self._dir.consensus_seeds,
                          divisor=self.dir.seeds_divisor,
                          output=self._dir.consensus_seeds)
            main(kwargs)
        assert os.path.isfile(self._dir.consensus_seeds)
        return

################################################
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
        # Output: self._dir.pssm | self._dir.converge_composition
        assert os.path.isfile(self._dir.converge_seeds)

        seeds = self._dir.converge_seeds.rsplit("/", maxsplit=1)[-1]
        input_seqs = self.dir.input_seqs.rsplit("/", maxsplit=1)[-1]
        composition = self.dir.converge_composition.rsplit("/", maxsplit=1)[-1]
        converge_exec = self.dir.converge_exec.rsplit("/", maxsplit=1)[-1]

        shutil.copy(self._dir.converge_seeds, self.dir.converge_dir)
        shutil.copy(self.dir.input_seqs, self.dir.converge_dir)

        command = f"cd {self.dir.converge_dir} && " \
            f"mpirun -np {self.dir.num_p} ./{converge_exec} -B -E 1 -r 1 -f 1 " \
            f"-c {composition} -i {seeds} -p {input_seqs}"
        subprocess.run(command, shell=True, executable=self.dir.bash_exec,
                       stdout=subprocess.DEVNULL)
        if os.path.isfile(self._dir.pssm):
            self.to_trash(self._dir.pssm)
        if os.path.isfile(self._dir.converge_composition):
            self.to_trash(self._dir.converge_composition)
        shutil.move(self.dir.converge_output, self._dir.pssm)
        shutil.move(self.dir.converge_composition,
                    self._dir.converge_composition)
        self.to_trash(self.dir.converge_discard)
        self.to_trash(f"{self.dir.converge_dir}/{seeds}")
        self.to_trash(f"{self.dir.converge_dir}/{input_seqs}")
        if __debug__:
            shutil.copy(self._dir.pssm, self._dir.converge_pssm)
        assert os.path.isfile(self._dir.pssm)
        assert os.path.isfile(self._dir.converge_composition)
        return

    def to_meme_format(self):
        # Input: self._dir.pssm | self._dir.converge_composition
        # Output: self._dir.pssm
        assert os.path.isfile(self._dir.pssm)
        assert os.path.isfile(self._dir.converge_composition)
        from converge2meme import main
        kwargs = dict(input_pssm=self._dir.pssm,
                      composition=self._dir.converge_composition,
                      output=self._dir.pssm)
        main(kwargs)
        if __debug__:
            shutil.copy(self._dir.pssm, self._dir.converge_meme)
        assert os.path.isfile(self._dir.pssm)
        return

    def screen_converge_pssm(self):
        # Input: self._dir.pssm
        # Output: self._dir.pssm
        assert os.path.isfile(self._dir.pssm)
        from filter import Filter
        filter_dir = Directory.filter_dir._replace(
            memefile=self._dir.pssm,
            short_seq=self._dir.short_seq)
        conv_filter = Filter(filter_dir)
        conv_filter.run(to_run=['correlated', 'noncombi'])
        if __debug__:
            shutil.copy(self._dir.pssm, self._dir.meme_cleaned)
        assert os.path.isfile(self._dir.pssm)
        return True

################################################
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
                  f'-nostatus &>> {self._dir.pssm}'
        subprocess.run(command, shell=True, executable=self.dir.bash_exec)
        if __debug__:
            shutil.copy(self._dir.pssm, self._dir.starter_meme)
        assert os.path.isfile(self._dir.pssm)
        return

    def clean_pssm(self):
        # Input: self._dir.meme_full | self._dir.starter_meme
        # Output: self._dir.meme_full | self._dir.starter_meme
        assert os.path.isdir(self._dir.meme_full)
        assert os.path.isfile(self._dir.pssm)
        from meme_cleaner import main
        kwargs = dict(input=self._dir.pssm,
                      output=self._dir.pssm)
        main(kwargs)
        for filename in os.listdir(self._dir.meme_full):
            kwargs = dict(input=f"{self._dir.meme_full}/{filename}",
                          output=f"{self._dir.meme_full}/{filename}")
            main(kwargs)
        if __debug__:
            shutil.copy(self._dir.pssm, self._dir.starter_meme)
        assert os.path.isdir(self._dir.meme_full)
        assert os.path.isfile(self._dir.pssm)
        return

    def merge_pssm(self):
        # Input: self._dir.meme_full | self._dir.starter_meme
        # Output: self._dir.meme_merged
        assert os.path.isdir(self._dir.meme_full)
        assert os.path.isfile(self._dir.pssm)
        from merge_meme_files import main
        kwargs = dict(meme_folder=self._dir.meme_full,
                      memefile=self._dir.pssm)
        main(kwargs)
        if __debug__:
            shutil.copy(self._dir.pssm, self._dir.meme_merged)
        assert os.path.isfile(self._dir.pssm)
        return

    def screen_pssm(self):
        # Input: self._dir.meme_merged | self.dir.input_seqs |
        #        self._dir.short_seq
        # Output: self._dir.meme_cleaned
        assert os.path.isfile(self._dir.pssm)
        filter_dir = Directory.filter_dir._replace(
            memefile=self._dir.pssm,
            short_seq=self._dir.short_seq)
        Filter(filter_dir).run()
        if __debug__:
            shutil.copy(self._dir.pssm, self._dir.meme_cleaned)
        assert os.path.isfile(self._dir.pssm)
        return True
    ###############################################

    # Resume
    def assemble_combi(self):
        # Input: self._dir.meme_cleaned | self.dir.input_seqs
        # Output: self.dir.output_mast
        assert os.path.isfile(self._dir.pssm)
        self.to_trash(self.dir.output_mast)
        if not os.path.isdir('output'):
            os.mkdir('output')
        command = f'{self.dir.meme_dir}/mast -remcorr ' \
                  f'{self._dir.pssm} {self.dir.input_seqs} -o ' \
                  f'{self.dir.output_mast}'
        subprocess.run(command, shell=True, executable=self.dir.bash_exec)
        assert os.path.isdir(self.dir.output_mast)
        assert os.path.isfile(f"{self.dir.output_mast}/mast.txt")
        return

    def cluster(self):
        # Input:self.dir.output_mast | self._dir.meme_cleaned
        # Output: self.dir.output_clusters | self.dir.output_logos
        assert os.path.isfile(self._dir.pssm)
        assert os.path.isfile(f"{self.dir.output_mast}/mast.txt")
        cluster_dir = Directory.cluster_dir._replace(
            input_mast=f"{self.dir.output_mast}/mast.txt",
            input_meme=self._dir.pssm,
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
                    print(f"Error: {func.__name__}")
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