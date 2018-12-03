from collections import namedtuple
import re
import os
import shutil
import subprocess

from utils import move_replace

ClusterIntDir = namedtuple(
    "ClusterIntDir", "full_param_pkl cluster_pkl motifs")

class Cluster:
    def __init__(self, dir):
        self.dir = dir
        self._dir = self.set_internal_dir()

    def set_internal_dir(self):
        # cluster_pkl only active when logo required but pkl not required.
        # This should then be deleted once cluster finish running.
        _dir = ClusterIntDir(
            full_param_pkl=f"{self.dir.file}/full_param.pkl",
            cluster_pkl=f"{self.dir.file}/_cluster.pkl",
            motifs=f"{self.dir.file}/motifs")
        return _dir

    def run(self):
        to_run = [self.get_cluster_params, self.cluster_combi]
        if self.dir.logos:
            to_run.append(self.create_cluster_motifs)
            to_run.append(self.create_cluster_logos)
            to_run.append(self.rename_logos)
        for func in to_run:
            try:
                print(f"Cluster/{func.__name__}:")
                func()
                print("Cluster/Success!\n")
            except:
                print("Cluster/Error: {}".format(func.__name__))
                raise
        return


    def get_cluster_params(self):
        # Input: ./files/mast_onlycombi/mast.txt
        # Output: ./files/clustering_df.pkl
        assert os.path.isfile(self.dir.input_mast), self.dir.input_mast
        from generate_cluster_params import main
        kwargs = dict(input_mast=self.dir.input_mast,
                      screen_threshold=5,
                      pkl_path=self._dir.full_param_pkl)
        main(kwargs)
        assert os.path.isfile(self._dir.full_param_pkl)
        return

    def cluster_combi(self):
        # Input: ./files/clustering_df.pkl
        # Output: output/cluster_description.txt
        #         (optional) files/cluster_centroids.pkl
        assert os.path.isfile(self._dir.full_param_pkl)
        from assemble_cluster_output import main
        kwargs = dict(cluster_threshold=50,
                      full_param_pkl=self._dir.full_param_pkl,
                      output=self.dir.description,
                      cluster_df_pkl=self._dir.cluster_pkl)
        main(kwargs)
        if self.dir.description:
            assert os.path.isfile(self.dir.description)
        if self.dir.cluster_pkl:
            shutil.copy(self._dir.cluster_pkl, self.dir.cluster_pkl)
        return

    def create_cluster_motifs(self):
        # Input: ./files/meme_format3.txt    ./files/cluster_centroids.pkl
        # Output: multiple ./files/motifs/motifs_in_cluster_{}.txt
        assert os.path.isfile(self.dir.input_meme)
        assert os.path.isfile(self._dir.cluster_pkl)
        from split_motifs_by_cluster import main
        self.to_trash(self._dir.motifs)
        os.mkdir(self._dir.motifs)
        kwargs = dict()
        kwargs['centroid_pkl'] = self._dir.cluster_pkl
        kwargs['input_meme'] = self.dir.input_meme
        kwargs['motifs'] = self._dir.motifs
        main(kwargs)
        assert os.path.isdir(self._dir.motifs)
        return

    def create_cluster_logos(self):
        # Input: multiple ./files/motifs/motifs_in_cluster_{}.txt
        # Output: multiple ./output/logos/cluster_{}/logo_{}.png
        assert os.path.isdir(self._dir.motifs)
        self.to_trash(self.dir.logos)
        os.mkdir(self.dir.logos)
        for filename in os.listdir(self._dir.motifs):
            cluster_i = re.match("motifs_in_cluster_([0-9]+)\.txt",
                                 filename).group(1)
            i = 1
            path = f"{self.dir.logos}/cluster_{cluster_i}"
            os.mkdir(path)
            while True:
                command = f'{self.dir.meme_dir}/ceqlogo -i{i} ' \
                      f'{self._dir.motifs}/{filename} -o ' \
                      f'{path}/logo_{i}.png -f PNG'
                i += 1
                returncode = subprocess.run(
                    command, shell=True, executable=self.dir.bash_exec,
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)\
                    .returncode
                if returncode != 0:
                    break
        assert os.path.isdir(self.dir.logos)
        return

    def rename_logos(self):
        assert os.path.isdir(self.dir.logos)
        assert os.path.isdir(self._dir.motifs)
        from rename_cluster_logos import main
        kwargs = dict(motif_filedir=self._dir.motifs,
                      output_logodir=self.dir.logos)
        main(kwargs)
        assert os.path.isdir(self.dir.logos)
        return

    def delete_intermediate(self):
        for file in self._dir:
            self.to_trash(file)
        return

    def to_trash(self, file):
        return move_replace(file, self.dir.trash)


