import shutil
import subprocess
import os
from collections import OrderedDict, namedtuple
from utils import move_replace
from config import Dir

IntDirTemplate = namedtuple(
    '_dir', 'post_evalue post_entropy post_corr mast_singleseq mast_postcorr '
            'full_param_pkl cluster_param_pkl')

class Filter:
    def __init__(self):
        self.switches = self.set_switches()
        self.dir = Dir.filter_dir
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

    def to_trash(self, file):
        return move_replace(file, self.dir.trash)

    def set_switches(self):
        switches = OrderedDict()
        switches['SCREEN_EVALUE'] = (False, self.screen_evalue)
        switches['SCREEN_ENTROPY'] = (False, self.screen_entropy)
        switches['SCREEN_CORRELATED'] = (False, self.screen_correlated)
        switches['SCREEN_NON_COMBI'] = (False, self.screen_non_combi)
        return switches

    def set_internal_dir(self):
        _dir = IntDirTemplate(
            post_evalue='files/meme_evalue_screened.txt',
            post_entropy='files/meme_entropy_screened.txt',
            post_corr='files/meme_correlated_screened.txt',
            mast_singleseq='files/mast_single/',
            mast_postcorr='files/mast_nocorr/',
            full_param_pkl='files/clustering_df.pkl',
            cluster_param_pkl='files/cluster_final.pkl')
        return _dir

    def screen_evalue(self):
        # Input: ./files/meme_merged.txt
        # Output: ./files/meme_evalue_screened.txt
        print("screen_evalue:")
        assert os.path.isfile(self.dir.orig)
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
        command = f'{self.dir.meme_dir}mast -remcorr ' \
                  f'{self._dir.post_entropy} {self.dir.single_seq}'
        subprocess.run(command, shell=True, executable=self.dir.bash_exec)
        shutil.rmtree("", ignore_errors=True)
        os.mkdir(self._dir.mast_singleseq)
        move_replace('mast_out', self.dir.file, 'mast_single')
        from mast_remove_profiles import main
        kwargs = dict(meme_in=self._dir.post_entropy,
                      meme_out=self._dir.post_corr)
        main(kwargs)
        assert os.path.isdir(self._dir.mast_singleseq)
        assert os.path.isfile(self._dir.post_corr)
        print("Success!\n")
        return

    def screen_non_combi(self):
        # todo: split this up
        # Input: ./files/meme_format2.txt    ./files/single_seq.fasta
        # Output: ./files/mast_nocorr    ./files/meme_format3.txt
        print("screen_non_combi:")
        assert os.path.isfile(self._dir.post_corr)
        assert os.path.isfile(self.dir.single_seq)
        command = f"{self.dir.meme_dir}mast -remcorr " \
                  f"{self._dir.post_corr} {self.dir.orig}"
        subprocess.run(command, shell=True, executable=self.dir.bash_exec)
        shutil.rmtree(self._dir.mast_postcorr, ignore_errors=True)
        os.mkdir(self._dir.mast_postcorr)
        move_replace('mast_out', self.dir.file, 'mast_nocorr')

        from generate_cluster_params import main
        kwargs = dict(input_fname=self._dir.mast_postcorr+"mast.txt",
                      screen_threshold=5,
                      pkl_path=self._dir.full_param_pkl)
        main(kwargs)

        from cluster_final import main
        kwargs = dict(cluster_threshold=50,
                      pkl_path=self._dir.full_param_pkl,
                      output=None,
                      cluster_df_pkl="./files/cluster_final.pkl")
        main(kwargs)

        from mast_remove_profiles_using_pkl import main
        kwargs = dict(cluster_df_pkl="./files/cluster_final.pkl",
                      input=self._dir.post_corr,
                      output=self.dir.cleaned)
        main(kwargs)
        self.to_trash(self._dir.post_corr)
        self.to_trash(self._dir.cluster_param_pkl)
        assert os.path.isfile(self.dir.cleaned)
        assert os.path.isdir(self._dir.mast_postcorr)
        print("Success!\n")
        return