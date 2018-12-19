from collections import OrderedDict, namedtuple
import os
import shutil
import subprocess
import sys

# Add src folders to sys.path so import is straightforward
cwd = os.getcwd()
src = f"{cwd}/src"
cluster = f"{cwd}/src/cluster"
filter_ = f"{cwd}/src/filter"
sys.path.append(src)
sys.path.append(cluster)
sys.path.append(filter_)

from cluster import Cluster
from filter import Filter
from config import Directory
from utils import move_replace, check_fasta_validity
from pssm_parser import PSSM

ExecIntDir = namedtuple(
    'ExecIntDir', 'dhcl_output consensus_seeds converge_seeds meme_full '
                  'converge_output pssm meme_merged meme_cleaned '
                  'converge_meme converge_composition short_seq short_seq_len')

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
            meme_merged=f"{self.dir.file}/meme_merged.txt",
            meme_cleaned=f"{self.dir.file}/meme_cleaned.txt",
            converge_meme=f"{self.dir.file}/converge_meme.txt",
            pssm=f"{self.dir.file}/pssm.txt",
            converge_composition=f"{self.dir.file}/converge_composition.txt",
            short_seq=f"{self.dir.file}/short_seq.fasta",
            short_seq_len=3)
        return _dir

    def set_switches(self):
        switches = OrderedDict()
        # Get input_seqs
        switches['MERGE_INPUT'] = (True, self.merge_input)
        switches['SHRINK_INPUT'] = (True, self.shrink_input)
        switches['CREATE_SHORT_SEQS'] = (True, self.create_short_seqs)
        # Get consensus_loops
        switches['RUN_DHCL'] = (True, self.run_dhcl)
        switches['EXTRACT_CONSENSUS'] = (True, self.extract_consensus)
        switches['REDUCE_CONSENSUS'] = (True, self.reduce_consensus)
        # Get PSSM using:
        # Meme
        switches['BUILD_PSSM'] = (True, self.build_pssm)
        switches['CLEAN_PSSM'] = (True, self.clean_pssm)
        switches['TO_MINIMAL'] = (True, self.to_minimal)
        switches['MERGE_PSSM'] = (True, self.merge_pssm)
        # Or converge
        switches['BUILD_CONVERGE_SEEDS'] = (False, self.build_converge_seeds)
        switches['RUN_CONVERGE'] = (False, self.run_converge)
        switches['CONV_TO_MINIMAL'] = (False, self.conv_to_minimal)
        # Get combi
        switches['SCREEN_PSSM'] = (False, self.screen_pssm)
        switches['ASSEMBLE_COMBI'] = (True, self.assemble_combi)
        switches['CLUSTER'] = (True, self.cluster)
        switches['DELETE_INTERMEDIATE'] = (True, self.delete_intermediate)
        return switches

    def merge_input(self):
        # Input: self.dir.input_seqdir
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
        if self.dir.seq_divisor:
            from shrink_input_for_test import main
            kwargs = dict(seqs=self.dir.input_seqs, output=self.dir.input_seqs,
                          divisor=self.dir.seq_divisor)
            main(kwargs)
        check_fasta_validity(self.dir.input_seqs)
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
        check_fasta_validity(self._dir.short_seq)
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
        from converters import dhcl_to_cons
        kwargs = dict(dhcl_dir=self._dir.dhcl_output,
                      fasta_dir=self.dir.fasta_for_pdb,
                      output=self._dir.consensus_seeds)
        dhcl_to_cons(kwargs)
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
        from converters import cons_to_conv_input
        kwargs = dict(seed_seqs=self._dir.consensus_seeds,
                      output=self._dir.converge_seeds)
        cons_to_conv_input(kwargs)
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
            f"mpirun -np {self.dir.num_p} ./{converge_exec} -B -E 1 -r 1 " \
            f"-f 1 -c {composition} -i {seeds} -p {input_seqs}"
        subprocess.run(command, shell=True, executable=self.dir.bash_exec,
                       stdout=subprocess.DEVNULL)
        if os.path.isfile(self._dir.pssm):
            self.to_trash(self._dir.pssm)
        if os.path.isfile(self._dir.converge_composition):
            self.to_trash(self._dir.converge_composition)
        # Moving this to files folder instead of pipeline folder
        shutil.move(self.dir.converge_output, self._dir.converge_output)
        shutil.move(self.dir.converge_composition,
                    self._dir.converge_composition)
        self.to_trash(self.dir.converge_discard)
        self.to_trash(f"{self.dir.converge_dir}/{seeds}")
        self.to_trash(f"{self.dir.converge_dir}/{input_seqs}")
        assert os.path.isfile(self._dir.converge_output)
        assert os.path.isfile(self._dir.converge_composition)
        return

    def conv_to_minimal(self):
        # Input: self._dir.converge_output | self._dir.converge_composition
        # Output: self._dir.pssm
        assert os.path.isfile(self._dir.converge_output)
        assert os.path.isfile(self._dir.converge_composition)
        from converters import converge_to_minimal
        kwargs = dict(input_pssm=self._dir.converge_output,
                      composition=self._dir.converge_composition,
                      output=self._dir.pssm)
        converge_to_minimal(kwargs)
        if __debug__:
            shutil.copy(self._dir.pssm, self._dir.converge_meme)
        assert os.path.isfile(self._dir.pssm)
        return

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

    def clean_pssm(self):
        # Input: self._dir.meme_full
        # Output: self._dir.meme_full
        assert os.path.isdir(self._dir.meme_full)
        from meme_cleaner import main
        for filename in os.listdir(self._dir.meme_full):
            kwargs = dict(input=f"{self._dir.meme_full}/{filename}",
                          output=f"{self._dir.meme_full}/{filename}")
            main(kwargs)
        assert os.path.isdir(self._dir.meme_full)
        return

    def to_minimal(self):
        # Input: self._dir.meme_full
        # Output: self._dir.meme_full
        assert os.path.isdir(self._dir.meme_full)
        from converters import meme_to_minimal
        for _fname in os.listdir(self._dir.meme_full):
            filename = f"{self._dir.meme_full}/{_fname}"
            kwargs = dict(input=filename,
                          output=filename)
            meme_to_minimal(kwargs)
        assert os.path.isdir(self._dir.meme_full)
        return

    def merge_pssm(self):
        # Input: self._dir.meme_full
        # Output: self._dir.meme_merged
        assert os.path.isdir(self._dir.meme_full)
        main_instance = None
        for _fname in os.listdir(self._dir.meme_full):
            pssm_instance = PSSM(filename=f"{self._dir.meme_full}/{_fname}")
            if not main_instance:
                main_instance = pssm_instance
                continue
            main_instance.merge_with(pssm_instance)
        main_instance.relabel_pssms()
        main_instance.output(self._dir.pssm)
        if __debug__:
            shutil.copy(self._dir.pssm, self._dir.meme_merged)
        assert os.path.isfile(self._dir.pssm)
        return

    ###############################################
    # Resume
    def screen_pssm(self):
        # Input: self._dir.pssm | self.dir.input_seqs |
        #        self._dir.short_seq
        # Output: self._dir.meme_cleaned
        assert os.path.isfile(self._dir.pssm)
        from filter import Filter
        filter_dir = Directory.filter_dir._replace(
            input_seqs=self.dir.input_seqs,
            memefile=self._dir.pssm,
            short_seq=self._dir.short_seq)
        Filter(filter_dir).run()
        if __debug__:
            shutil.copy(self._dir.pssm, self._dir.meme_cleaned)
        assert os.path.isfile(self._dir.pssm)
        return True

    def assemble_combi(self):
        # Input: self._dir.pssm | self.dir.input_seqs
        # Output: self.dir.output_mast
        assert os.path.isfile(self._dir.pssm)
        assert os.path.isfile(self.dir.input_seqs)
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
        Cluster(cluster_dir).run(make_logo=True)
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
        if self.switches['MERGE_INPUT'][0]:
            self.to_trash(self.dir.input_seqs)
        Filter(Directory.filter_dir).delete_intermediate()
        Cluster(Directory.cluster_dir).delete_intermediate()
        return

if __name__ == "__main__":
    executor = Executor()
    executor.run()