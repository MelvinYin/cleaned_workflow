from collections import namedtuple
import re
import os
import shutil
import subprocess

from utils import move_into
import pickle
from pssm_parser import PSSM

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

    def run(self, make_logo=False):
        if make_logo:
            assert self.dir.logos
            to_run = [self.get_cluster_params,
                      self.cluster_combi,
                      self.create_cluster_motifs,
                      self.create_cluster_logos,
                      self.rename_logos]
        else:
            to_run = [self.get_cluster_params, self.cluster_combi]
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
        assert os.path.isfile(self.dir.input_mast)
        from generate_cluster_params import main
        kwargs = dict(input_mast=self.dir.input_mast,
                      combi_minsize=self.dir.combi_minsize,
                      pkl_path=self._dir.full_param_pkl,
                      num_cluster=self.dir.num_cluster)
        main(kwargs)
        assert os.path.isfile(self._dir.full_param_pkl)
        return

    def cluster_combi(self):
        # Input: ./files/clustering_df.pkl
        # Output: output/cluster_description.txt
        #         (optional) files/cluster_centroids.pkl
        assert os.path.isfile(self._dir.full_param_pkl)
        from assemble_cluster_output import main
        kwargs = dict(cluster_threshold=self.dir.cluster_minsize,
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
        self.to_trash(self._dir.motifs)
        os.mkdir(self._dir.motifs)

        with open(self._dir.cluster_pkl, 'rb') as file:
            cluster_centroids = pickle.load(file)
        for label, centroid in cluster_centroids['centroid'].items():
            output = f"{self._dir.motifs}/motifs_in_cluster_{label}.txt"
            pssm_obj = PSSM(filename=self.dir.input_meme)
            pssm_obj.keep(centroid)
            pssm_obj.output(output)
        os.mkdir(f"{self.dir.output}/motifs")
        for file in os.listdir(self._dir.motifs):
            src_path = f"{self._dir.motifs}/{file}"
            dst_path = f"{self.dir.output}/motifs"
            shutil.copy(src_path, dst_path)
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
        from rename_cluster_logos import rename
        kwargs = dict(motif_filedir=self._dir.motifs,
                      output_logodir=self.dir.logos)
        rename(kwargs)
        assert os.path.isdir(self.dir.logos)
        return

    def delete_intermediate(self):
        for file in self._dir:
            self.to_trash(file)
        return

    def to_trash(self, file):
        if not (os.path.isdir(file) or os.path.isfile(file)):
            return
        return move_into(file, self.dir.trash)


