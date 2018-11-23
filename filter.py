import shutil
import subprocess
import os
from collections import OrderedDict, namedtuple
from utils import move_replace, rename
from config import ClusterDir
from cluster import Cluster

FilterIntDir = namedtuple(
    'FilterIntDir', 'post_evalue post_entropy post_corr mast_singleseq '
                    'mast_postcorr cluster_pkl')

class Filter:
    def __init__(self, dir):
        self.switches = self.set_switches()
        self.dir = dir
        self._dir = self.set_internal_dir()

    def run(self):
        for _switch, (to_run, function) in self.switches.items():
            if to_run:
                try:
                    function()
                except:
                    print("Error: {}".format(function.__name__))
                    raise
        return

    def delete_intermediate(self):
        print("delete_intermediate:")
        for file in self._dir:
            self.to_trash(file)
        print("Success!\n")
        return

    def to_trash(self, file):
        return move_replace(file, self.dir.trash)

    def set_switches(self):
        switches = OrderedDict()
        switches['SCREEN_EVALUE'] = (True, self.screen_evalue)
        switches['SCREEN_ENTROPY'] = (True, self.screen_entropy)
        switches['SCREEN_CORRELATED'] = (True, self.screen_correlated)
        switches['SCREEN_NON_COMBI'] = (True, self.screen_non_combi)
        return switches

    def set_internal_dir(self):
        _dir = FilterIntDir(
            post_evalue='files/meme_evalue_screened.txt',
            post_entropy='files/meme_entropy_screened.txt',
            post_corr='files/meme_correlated_screened.txt',
            mast_singleseq='files/mast_single',
            mast_postcorr='files/mast_nocorr',
            cluster_pkl='files/cluster_final.pkl')
        return _dir

    def screen_evalue(self):
        # Input: ./files/meme_merged.txt
        # Output: ./files/meme_evalue_screened.txt
        print("screen_evalue:")
        assert os.path.isfile(self.dir.orig), self.dir.orig
        from remove_motifs_with_low_evalue import main
        kwargs = dict(meme=self.dir.orig,
                      output=self._dir.post_evalue)
        main(kwargs)
        assert os.path.isfile(self._dir.post_evalue)
        print("Success!\n")
        return

    def screen_entropy(self):
        # Input: ./files/meme_evalue_screened.txt
        # Output: ./files/meme_format.txt
        print("screen_entropy:")
        assert os.path.isfile(self._dir.post_evalue)
        from screen_motif_entropy import main
        kwargs = dict(meme=self._dir.post_evalue,
                      output=self._dir.post_entropy)
        main(kwargs)
        assert os.path.isfile(self._dir.post_entropy)
        print("Success!\n")
        return

    def screen_correlated(self):
        # Input: ./files/meme_format.txt    ./files/single_seq.fasta
        # Output: ./files/mast_single    ./files/meme_format2.txt
        print("screen_correlated:")
        assert os.path.isfile(self._dir.post_entropy)
        assert os.path.isfile(self.dir.single_seq)
        if os.path.isdir(self._dir.mast_singleseq):
            shutil.rmtree(self._dir.mast_singleseq)
        command = f'{self.dir.meme_dir}/mast -remcorr ' \
                  f'{self._dir.post_entropy} {self.dir.single_seq} -o ' \
                  f'{self._dir.mast_singleseq}'
        subprocess.run(command, shell=True, executable=self.dir.bash_exec)
        from mast_remove_profiles import main
        kwargs = dict(meme_in=self._dir.post_entropy,
                      meme_out=self._dir.post_corr)
        main(kwargs)
        assert os.path.isdir(self._dir.mast_singleseq)
        assert os.path.isfile(self._dir.post_corr)
        print("Success!\n")
        return

    def screen_non_combi(self):
        # Input: ./files/meme_format2.txt    ./files/input_seqs.fasta
        # Output: ./files/mast_nocorr    ./files/meme_format3.txt
        print("screen_non_combi:")
        assert os.path.isfile(self._dir.post_corr)
        assert os.path.isfile(self.dir.input_seqs)
        if os.path.isdir(self._dir.mast_postcorr):
            shutil.rmtree(self._dir.mast_postcorr, ignore_errors=True)
        command = f"{self.dir.meme_dir}/mast -remcorr " \
                  f"{self._dir.post_corr} {self.dir.input_seqs} -o " \
                  f"{self._dir.mast_postcorr}"
        subprocess.run(command, shell=True, executable=self.dir.bash_exec)
        cluster_dir = ClusterDir(
            file=self.dir.file,
            log=self.dir.log,
            trash=self.dir.trash,
            input_mast=f"{self._dir.mast_postcorr}/mast.txt",
            input_meme=self._dir.post_corr,
            output_mast=None,
            description=None,
            logos=None,
            cluster_pkl=self._dir.cluster_pkl,
            meme_dir=self.dir.meme_dir,
            bash_exec=self.dir.bash_exec)
        cluster = Cluster(cluster_dir)
        cluster.run()
        cluster.delete_intermediate()

        from mast_remove_profiles_using_pkl import main
        kwargs = dict(cluster_df_pkl=self._dir.cluster_pkl,
                      input=self._dir.post_corr,
                      output=self.dir.cleaned)
        main(kwargs)
        self.to_trash(self._dir.post_corr)
        self.to_trash(self._dir.cluster_pkl)
        assert os.path.isfile(self.dir.cleaned)
        assert os.path.isdir(self._dir.mast_postcorr)
        print("Success!\n")
        return