import os
import shutil
from config import Dir, namedtuple
import subprocess
from utils import move_replace
from collections import OrderedDict

ClusterIntDir = namedtuple(
    "ClusterIntDir", "full_param_pkl cluster_pkl motifs")

class Cluster:
    def __init__(self, dir):
        self.switches = self.set_switches()
        self.dir = dir
        self._dir = self.set_internal_dir()

    def set_internal_dir(self):
        # cluster_pkl only active when logo required but pkl not required.
        # This should then be deleted once cluster finish running.
        _dir = ClusterIntDir(
            full_param_pkl='files/full_param.pkl',
            cluster_pkl="files/_cluster.pkl",
            motifs="files/motifs/")
        return _dir

    def run(self):
        self.get_cluster_params()
        self.cluster_combi()
        if self.dir.logos:
            self.create_cluster_motifs()
            self.create_cluster_logos()
        return

    def reset(self):
        for name in self._dir:
            self.to_trash(name)

    def to_trash(self, file):
        return move_replace(file, self.dir.trash)

    def set_switches(self):
        switches = OrderedDict()
        switches['GET_CLUSTER_PARAMS'] = (True, self.get_cluster_params)
        switches['CLUSTER_COMBI'] = (True, self.cluster_combi)
        switches['CREATE_CLUSTER_MOTIFS'] = (True, self.create_cluster_motifs)
        switches['CREATE_CLUSTER_LOGOS'] = (True, self.create_cluster_logos)
        return switches

    def get_cluster_params(self):
        # Input: ./files/mast_onlycombi/mast.txt
        # Output: ./files/clustering_df.pkl
        print("get_cluster_params:")
        assert os.path.isfile(self.dir.input_mast)
        from generate_cluster_params import main
        kwargs = dict(input_fname=self.dir.input_mast,
                      screen_threshold=5,
                      pkl_path=self._dir.full_param_pkl)
        main(kwargs)
        assert os.path.isfile(self._dir.full_param_pkl)
        print("Success!\n")
        return

    def cluster_combi(self):
        # Input: ./files/clustering_df.pkl
        # Output: output/cluster_description.txt
        #         (optional) files/cluster_centroids.pkl
        print("cluster_combi:")
        assert os.path.isfile(self._dir.full_param_pkl)
        from cluster_final import main
        kwargs = dict(cluster_threshold=50,
                      pkl_path=self._dir.full_param_pkl,
                      output=self.dir.description,
                      cluster_df_pkl=self._dir.cluster_pkl)
        main(kwargs)
        if self.dir.description:
            assert os.path.isfile(self.dir.description)
        if self.dir.cluster_pkl:
            shutil.copy(self._dir.cluster_pkl, self.dir.cluster_pkl)
        print("Success!\n")
        return

    def create_cluster_motifs(self):
        # Input: ./files/meme_format3.txt    ./files/cluster_centroids.pkl
        # Output: multiple ./files/motifs/motifs_in_cluster_{}.txt
        print("create_cluster_motifs:")
        assert os.path.isfile(self.dir.input_meme)
        assert os.path.isfile(self._dir.cluster_pkl)
        from split_motifs_into_cluster_motifs import main
        self.to_trash(self._dir.motifs)
        os.mkdir(self._dir.motifs)
        kwargs = dict()
        kwargs['centroid_pkl'] = self._dir.cluster_pkl
        kwargs['input_meme'] = self.dir.input_meme
        kwargs['motifs'] = self._dir.motifs
        main(kwargs)
        assert os.path.isdir(self._dir.motifs)
        print("Success!\n")
        return

    def create_cluster_logos(self):
        # Input: multiple ./files/motifs/motifs_in_cluster_{}.txt
        # Output: multiple ./output/logos/cluster_{}/logo_{}.png
        print("create_cluster_logos:")
        assert os.path.isdir(self._dir.motifs)
        self.to_trash(self.dir.logos)
        os.mkdir(self.dir.logos)
        for filename in os.listdir(self._dir.motifs):
            cluster_i = filename[18:-4]
            i = 1
            path = self.dir.logos+"cluster_{}".format(cluster_i)
            os.mkdir(path)
            while True:
                command = f'{self.dir.meme_dir}ceqlogo -i{i} ' \
                      f'{self._dir.motifs}{filename} -o ' \
                      f'{path}/logo_{i}.png -f PNG'
                i += 1
                returncode = subprocess.run(
                    command, shell=True, executable=self.dir.bash_exec,
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)\
                    .returncode
                if returncode != 0:
                    break
        from rename_cluster_logos import main
        kwargs = dict(motif_filedir=self._dir.motifs,
                      output_logodir=self.dir.logos)
        main(kwargs)
        assert os.path.isdir(self.dir.logos)
        print("Success!\n")
        return


