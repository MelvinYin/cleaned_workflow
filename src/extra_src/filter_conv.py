from collections import namedtuple, OrderedDict
import os
import shutil
import subprocess

from cluster import Cluster
from config import Directory
from utils import move_replace

FilterIntDir = namedtuple(
    'FilterIntDir', "post_evalue post_entropy post_corr mast_shortseq "
                    "mast_postcorr cluster_pkl")

class Filter:
    def __init__(self, dir):
        self.dir = dir
        self._dir = self.set_internal_dir()
        self.cluster = None
        self.switches = self.set_switches()

    def set_switches(self):
        switches = OrderedDict()
        # Get input_seqs
        switches['ENTROPY'] = (False, self.screen_entropy)
        switches['CORRELATED'] = (False, self.screen_correlated)
        switches['NON_COMBI'] = (False, self.screen_non_combi)
        return switches

    def run(self):
        for (to_run, func) in self.switches.values():
            if to_run:
                try:
                    print(f"Filter/{func.__name__}:")
                    func()
                    print("Filter/Success!\n")
                except:
                    print("Filter/Error: {}".format(func.__name__))
                    raise
        return

    def set_internal_dir(self):
        _dir = FilterIntDir(
            post_evalue=f"{self.dir.file}/meme_evalue_screened.txt",
            post_entropy=f"{self.dir.file}/meme_entropy_screened.txt",
            post_corr=f"{self.dir.file}/meme_correlated_screened.txt",
            mast_shortseq=f"{self.dir.file}/mast_single",
            mast_postcorr=f"{self.dir.file}/mast_nocorr",
            cluster_pkl=f"{self.dir.file}/cluster_final.pkl")
        return _dir

    def screen_entropy(self):
        # Input: self.dir.memefile
        # Output: self.dir.memefile
        assert os.path.isfile(self.dir.memefile)
        from screen_entropy_conv import main
        kwargs = dict(pssm=self.dir.memefile,
                      entropy_threshold=40)
        main(kwargs)
        if __debug__:
            shutil.copy(self.dir.memefile, self._dir.post_entropy)
        assert os.path.isfile(self.dir.memefile)
        return

    def screen_correlated(self):
        # Input: self.dir.memefile | self.dir.short_seq
        # Output: self.dir.memefile | self._dir.mast_shortseq
        assert os.path.isfile(self.dir.memefile)
        assert os.path.isfile(self.dir.short_seq)
        if os.path.isdir(self._dir.mast_shortseq):
            shutil.rmtree(self._dir.mast_shortseq)
        command = f"{self.dir.meme_dir}/mast -remcorr " \
                  f"{self.dir.memefile} {self.dir.short_seq} -o " \
                  f"{self._dir.mast_shortseq}"
        subprocess.run(command, shell=True, executable=self.dir.bash_exec)
        from mast_remove_profiles import main
        kwargs = dict(mast_in=f"{self._dir.mast_shortseq}/mast.txt",
                      memefile=self.dir.memefile)
        main(kwargs)
        if __debug__:
            shutil.copy(self.dir.memefile, self._dir.post_corr)
        assert os.path.isdir(self._dir.mast_shortseq)
        assert os.path.isfile(self.dir.memefile)
        return

    def screen_non_combi(self):
        # Input: self.dir.memefile | self.dir.input_seqs
        # Output: self.dir.memefile | self._dir.mast_postcorr
        assert os.path.isfile(self.dir.memefile)
        assert os.path.isfile(self.dir.input_seqs)
        if os.path.isdir(self._dir.mast_postcorr):
            shutil.rmtree(self._dir.mast_postcorr, ignore_errors=True)
        command = f"{self.dir.meme_dir}/mast -remcorr " \
                  f"{self.dir.memefile} {self.dir.input_seqs} -o " \
                  f"{self._dir.mast_postcorr}"
        subprocess.run(command, shell=True, executable=self.dir.bash_exec)
        cluster_dir = Directory.cluster_dir._replace(
            input_mast=f"{self._dir.mast_postcorr}/mast.txt",
            input_meme=self.dir.memefile,
            cluster_pkl=self._dir.cluster_pkl)
        Cluster(cluster_dir).run()

        from remove_noncentroid_profiles import main
        kwargs = dict(cluster_df_pkl=self._dir.cluster_pkl,
                      memefile=self.dir.memefile)
        main(kwargs)
        assert os.path.isfile(self.dir.memefile)
        return

    def delete_intermediate(self):
        for file in self._dir:
            self.to_trash(file)
        return

    def to_trash(self, file):
        return move_replace(file, self.dir.trash)