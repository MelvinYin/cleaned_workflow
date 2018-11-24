from collections import namedtuple
import os
import shutil
import subprocess

from cluster import Cluster
from config import Directory
from utils import move_replace


FilterIntDir = namedtuple(
    'FilterIntDir', "post_evalue post_entropy post_corr mast_singleseq "
                    "mast_postcorr cluster_pkl")

class Filter:
    def __init__(self, dir):
        self.dir = dir
        self._dir = self.set_internal_dir()

    def run(self):
        to_run = (self.screen_evalue, self.screen_entropy,
                        self.screen_correlated, self.screen_non_combi,
                        self.delete_intermediate)
        for func in to_run:
            try:
                print(f"Filter/{func.__name__}:")
                func()
                print("Filter/Success!\n")
            except:
                print("Filter/Error: {}".format(func.__name__))
                raise
        return

    def delete_intermediate(self):
        for file in self._dir:
            self.to_trash(file)
        return

    def to_trash(self, file):
        return move_replace(file, self.dir.trash)

    def set_internal_dir(self):
        _dir = FilterIntDir(
            post_evalue=f"{self.dir.file}/meme_evalue_screened.txt",
            post_entropy=f"{self.dir.file}/meme_entropy_screened.txt",
            post_corr=f"{self.dir.file}/meme_correlated_screened.txt",
            mast_singleseq=f"{self.dir.file}/mast_single",
            mast_postcorr=f"{self.dir.file}/mast_nocorr",
            cluster_pkl=f"{self.dir.file}/cluster_final.pkl")
        return _dir

    def screen_evalue(self):
        # Input: ./files/meme_merged.txt
        # Output: ./files/meme_evalue_screened.txt
        assert os.path.isfile(self.dir.orig), self.dir.orig
        from screen_motif_evalue import main
        kwargs = dict(meme=self.dir.orig,
                      output=self._dir.post_evalue)
        main(kwargs)
        assert os.path.isfile(self._dir.post_evalue)
        return

    def screen_entropy(self):
        # Input: ./files/meme_evalue_screened.txt
        # Output: ./files/meme_format.txt
        assert os.path.isfile(self._dir.post_evalue)
        from screen_motif_entropy import main
        kwargs = dict(meme=self._dir.post_evalue,
                      output=self._dir.post_entropy)
        main(kwargs)
        assert os.path.isfile(self._dir.post_entropy)
        return

    def screen_correlated(self):
        # Input: ./files/meme_format.txt    ./files/single_seq.fasta
        # Output: ./files/mast_single    ./files/meme_format2.txt
        assert os.path.isfile(self._dir.post_entropy)
        assert os.path.isfile(self.dir.single_seq)
        if os.path.isdir(self._dir.mast_singleseq):
            shutil.rmtree(self._dir.mast_singleseq)
        command = f"{self.dir.meme_dir}/mast -remcorr " \
                  f"{self._dir.post_entropy} {self.dir.single_seq} -o " \
                  f"{self._dir.mast_singleseq}"
        subprocess.run(command, shell=True, executable=self.dir.bash_exec)
        from mast_remove_profiles import main
        kwargs = dict(meme_in=self._dir.post_entropy,
                      meme_out=self._dir.post_corr)
        main(kwargs)
        assert os.path.isdir(self._dir.mast_singleseq)
        assert os.path.isfile(self._dir.post_corr)
        return

    def screen_non_combi(self):
        # Input: ./files/meme_format2.txt    ./files/input_seqs.fasta
        # Output: ./files/mast_nocorr    ./files/meme_format3.txt
        assert os.path.isfile(self._dir.post_corr)
        assert os.path.isfile(self.dir.input_seqs)
        if os.path.isdir(self._dir.mast_postcorr):
            shutil.rmtree(self._dir.mast_postcorr, ignore_errors=True)
        command = f"{self.dir.meme_dir}/mast -remcorr " \
                  f"{self._dir.post_corr} {self.dir.input_seqs} -o " \
                  f"{self._dir.mast_postcorr}"
        subprocess.run(command, shell=True, executable=self.dir.bash_exec)
        cluster_dir = Directory.cluster_dir._replace(
            input_mast=f"{self._dir.mast_postcorr}/mast.txt",
            input_meme=self._dir.post_corr,
            cluster_pkl=self._dir.cluster_pkl)
        cluster = Cluster(cluster_dir)
        cluster.run()
        cluster.delete_intermediate()

        from remove_noncentroid_profiles import main
        kwargs = dict(cluster_df_pkl=self._dir.cluster_pkl,
                      input=self._dir.post_corr,
                      output=self.dir.cleaned)
        main(kwargs)
        self.to_trash(self._dir.post_corr)
        self.to_trash(self._dir.cluster_pkl)
        assert os.path.isfile(self.dir.cleaned)
        assert os.path.isdir(self._dir.mast_postcorr)
        return