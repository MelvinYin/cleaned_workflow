import shutil
import subprocess
import os

class Executor:
    def __init__(self):
        self.switches = self.set_switches()
        self.switch_map = self.set_switch_map()
        self.CONST = self.set_CONST()

    def set_CONST(self):
        CONST = dict()
        CONST['log'] = "files/log.txt"
        CONST['trash'] = "files/_trash"
        CONST['input_seqs'] = "files/input_seqs_copy.fasta"
        CONST['input_pdb'] = "files/input_pdb_test"
        return CONST

    def set_switches(self):
        switches = dict()
        switches['COPY_INIT_SEQS'] = False
        switches['RUN_DHCL'] = False
        switches['EXTRACT_CONSENSUS'] = False
        switches['RUN_MEME'] = False
        switches['BUILD_STARTER'] = False
        switches['CLEAN_MEME'] = False
        switches['MERGE_MEME'] = False
        switches['SCREEN_EVALUE'] = False
        switches['SCREEN_ENTROPY'] = False
        switches['SCREEN_CORRELATED'] = False
        switches['SCREEN_NON_COMBI'] = False
        switches['MAST_COMBI'] = False
        switches['CLUSTER_COMBI'] = False
        switches['CREATE_CLUSTER_MOTIFS'] = False
        switches['CREATE_CLUSTER_LOGOS'] = False
        switches['DELETE_INIT_FASTA'] = False
        switches['DELETE_INTERMEDIATE'] = False

        # For testing purposes
        switches['REDUCE_DHCL'] = False
        switches['SHRINK_INPUT'] = False
        return switches

    def set_switch_map(self):
        switch_map = dict()
        switch_map['COPY_INIT_SEQS'] = self.copy_init_seqs
        switch_map['RUN_DHCL'] = self.run_dhcl
        switch_map['EXTRACT_CONSENSUS'] = self.extract_consensus
        switch_map['RUN_MEME'] = self.run_meme
        switch_map['BUILD_STARTER'] = self.build_starter
        switch_map['CLEAN_MEME'] = self.clean_meme
        switch_map['MERGE_MEME'] = self.merge_meme
        switch_map['SCREEN_EVALUE'] = self.screen_evalue
        switch_map['SCREEN_ENTROPY'] = self.screen_entropy
        switch_map['SCREEN_CORRELATED'] = self.screen_correlated
        switch_map['SCREEN_NON_COMBI'] = self.screen_non_combi
        switch_map['MAST_COMBI'] = self.mast_combi
        switch_map['CLUSTER_COMBI'] = self.cluster_combi
        switch_map['CREATE_CLUSTER_MOTIFS'] = self.create_cluster_motifs
        switch_map['CREATE_CLUSTER_LOGOS'] = self.create_cluster_logos
        switch_map['DELETE_INIT_FASTA'] = self.delete_init_fasta
        switch_map['DELETE_INTERMEDIATE'] = self.delete_intermediate

        # For testing purposes
        switch_map['REDUCE_DHCL'] = self.reduce_dhcl
        switch_map['SHRINK_INPUT'] = self.shrink_input
        return switch_map

    def copy_init_seqs(self):
        # Input: ./files/input_seqs.fasta    ./files/single_seq.fasta
        # Output: ./external_scripts/pipeline/input_seqs.fasta    ./external_scripts/meme/input_seqs.fasta
        #         ./external_scripts/meme/single_seq.fasta
        print("copy_init_seqs:")
        shutil.copy(self.CONST['input_seqs'], './external_scripts/pipeline')
        shutil.copy('./files/single_seq.fasta', './external_scripts/meme')
        shutil.copy(self.CONST['input_seqs'], './external_scripts/meme')
        print("Success!\n")
        return

    def run_dhcl(self):
        # Input: ./files/input_pdb
        # Output: ./files/from_dhcl
        print("RUN_DHCL:")
        shutil.rmtree("files/from_dhcl", ignore_errors=True)
        os.mkdir("files/from_dhcl")
        command = 'source activate p2.7 ' \
                  '&& python ./external_scripts/dhcl/executables/everything.py ' \
                  '-d {} --outdir files/from_dhcl ' \
                  '&& source deactivate'.format(self.CONST['input_pdb'])
        subprocess.run(command.split())
        print("Success!\n")
        return

    def extract_consensus(self):
        # Input: ./files/from_dhcl    ./files/input_fasta
        # Output: ./files/init_seed_seqs.fasta
        print("extract_consensus:")
        from process_dhcl_output_meme import main
        kwargs = dict(dhcl_dir='files/from_dhcl',
                      fasta_dir='files/input_fasta',
                      output='files/init_seed_seqs.fasta')
        main(kwargs)
        print("Success!\n")
        return

    def run_meme(self):
        # Input: ./files/consensus_seqs.txt    ./files/input_seqs.fasta
        # Output: ./files/meme_full
        print("run_meme:")
        from run_meme_on_dhcl_loops import main
        kwargs = dict(consensus='files/consensus_seqs.txt',
                      seqs=self.CONST['input_seqs'],
                      output_folder='files/meme_full')
        main(kwargs)
        print("Success!\n")
        return

    def build_starter(self):
        # Input: ./files/input_seqs.txt
        # Output: ./files/meme_starter.txt
        print("build_starter:")
        command = './external_scripts/meme/bin/meme -text -protein -w 30 ' \
                  '-p 7 -nmotifs 2 -nostatus {} &>>./files/meme_starter.txt' \
                  ''.format(self.CONST['input_seqs'])
        subprocess.run(command.split())
        print("Success!\n")
        return

    def clean_meme(self):
        # Input: ./files/meme_full    ./files/meme_starter.txt
        # Output: ./files/meme_full    ./files/meme_starter.txt
        print("clean_meme:")
        from meme_cleaner import main
        kwargs = dict(input=self.CONST['input_seqs'],
                      output=self.CONST['input_seqs'])
        main(kwargs)
        for filename in os.listdir('files/meme_full'):
            kwargs = dict(input='files/'+filename,
                          output='files/'+filename)
            main(kwargs)
        print("Success!\n")
        return

    def merge_meme(self):
        # Input: ./files/meme_full    ./files/meme_starter.txt
        # Output: ./files/meme_merged.txt
        print("merge_meme:")
        from meme_merger import main
        kwargs = dict(meme_folder='files/meme_full',
                      output='files/meme_merged.txt',
                      meme_starter='files/meme_starter.txt')
        main(kwargs)
        print("Success!\n")
        return

    def screen_evalue(self):
        # Input: ./files/meme_merged.txt
        # Output: ./files/meme_evalue_screened.txt
        print("screen_evalue:")
        shutil.copy('files/meme_merged.txt', './files/meme.txt')
        from remove_motifs_with_low_evalue import main
        kwargs = dict(meme='files/meme.txt',
                      output='files/meme_evalue_screened.txt')
        main(kwargs)
        shutil.move('files/meme.txt', self.CONST['trash'])
        print("Success!\n")
        return

    def screen_entropy(self):
        # Input: ./files/meme_evalue_screened.txt
        # Output: ./files/meme_format.txt
        print("screen_entropy:")
        shutil.copy('files/meme_evalue_screened.txt', './files/meme.txt')
        from screen_motif_entropy import main
        kwargs = dict(meme='files/meme.txt',
                      output='files/meme_format.txt')
        main(kwargs)
        shutil.move('files/meme.txt', self.CONST['trash'])
        print("Success!\n")
        return

    def screen_correlated(self):
        # Input: ./files/meme_format.txt    ./files/single_seq.fasta
        # Output: ./files/mast_single    ./files/meme_format2.txt
        print("screen_correlated:")
        command = './external_scripts/meme/bin/mast -remcorr ' \
                  'files/meme_format.txt single_seq.fasta'
        subprocess.run(command.split())
        shutil.rmtree("files/mast_single", ignore_errors=True)
        os.mkdir("files/mast_single")
        shutil.move('external_scripts/meme/mast_out', 'files/mast_single')
        shutil.copy('files/mast_single/mast.txt', './files')
        from mast_remove_profiles import main
        kwargs = dict(meme_in='files/meme_format.txt',
                      meme_out='files/meme_format2.txt')
        main(kwargs)
        shutil.move('files/mast.txt', self.CONST['trash'])
        print("Success!\n")
        return

    def screen_non_combi(self):
        # Input: files/meme_format2.txt    files/single_seq.fasta
        # Output: files/mast_nocorr    files/meme_format3.txt
        print("screen_non_combi:")
        command = './external_scripts/meme/bin/mast -remcorr ' \
                  'files/meme_format2.txt single_seq.fasta'
        subprocess.run(command.split())
        shutil.rmtree("files/mast_nocorr", ignore_errors=True)
        os.mkdir("files/mast_nocorr")
        shutil.move('external_scripts/meme/mast_out', 'files/mast_nocorr')
        shutil.copy('files/mast_nocorr/mast.txt', './files')
        from cluster_get_profiles import main
        kwargs = dict()
        main(kwargs)
        from mast_remove_profiles_using_pkl import main
        kwargs = dict()
        main(kwargs)
        shutil.move('files/mast.txt', self.CONST['trash'])
        print("Success!\n")
        return

    def mast_combi(self):
        # Input: .iles/meme_format3.txt    ./external_scripts/meme/consolidated.fasta
        # Output: files/mast_onlycombi    ./output/mast
        print("mast_combi:")
        command = 'external_scripts/meme/bin/mast -remcorr ' \
                  'files/meme_format3.txt {}'.format(self.CONST['input_seqs'])
        subprocess.run(command.split())
        shutil.move("files/mast_onlycombi", self.CONST['trash'])
        os.mkdir("files/mast_onlycombi")
        shutil.move('external_scripts/meme/mast_out', 'files/mast_onlycombi')
        shutil.copy('files/mast_onlycombi', 'output/mast')
        print("Success!\n")
        return

    def get_cluster_params(self):
        # Input: files/mast_onlycombi/mast.txt
        # Output: files/clustering_df.pkl
        print("get_cluster_params:")
        shutil.copy('files/mast_onlycombi/mast.txt', './files')
        from generate_cluster_params import main
        kwargs = dict(input_fname="../files/mast.txt",
                      screen_threshold=5,
                      pkl_path="../files/clustering_df.pkl")
        main(kwargs)
        shutil.move('files/mast.txt', self.CONST['trash'])
        print("Success!\n")
        return

    def cluster_combi(self):
        # Input: files/clustering_df.pkl
        # Output: files/cluster_description.txt  output/cluster_description.txt
        print("cluster_combi:")
        from cluster_final import main
        kwargs = dict(cluster_threshold=50,
                      pkl_path="files/clustering_df.pkl",
                      output="files/cluster_description.txt")
        main(kwargs)
        shutil.copy('files/cluster_description.txt',
                    'output/cluster_description.txt')
        print("Success!\n")
        return

    def create_cluster_motifs(self):
        # Input: ./files/meme_format3.txt    ./files/cluster_centroids.pkl
        # Output: multiple ./files/motifs/motifs_in_cluster_{}.txt
        print("create_cluster_motifs:")
        from split_motifs_into_cluster_motifs import main
        shutil.move("files/motifs", self.CONST['trash'])
        os.mkdir("files/motifs")
        kwargs = dict()
        main(kwargs)
        print("Success!\n")
        return

    def create_cluster_logos(self):
        # Input: multiple ./files/motifs/motifs_in_cluster_{}.txt
        # Output: multiple ./output/logos/cluster_{}/logo_{}.png
        print("create_cluster_logos:")
        shutil.move("files/output/logos", self.CONST['trash'])
        os.mkdir("files/output/logos")
        for filename in os.listdir('files/motifs'):
            filename = filename[15:]
            cluster_i = filename[18:-4]
            i=1
            os.mkdir("files/output/cluster_{}".format(cluster_i))
            command = 'external_scripts/meme/ceqlogo -i{} ' \
                      './files/motifs/{} -o ' \
                      './output/logos/cluster_{}/logo_{}.png -f ' \
                      'PNG'.format(i, filename, cluster_i, i)
            while True:
                try:
                    subprocess.run(command.split())
                    i += 1
                except:
                    break
        from rename_cluster_logos import main
        kwargs = dict()
        main(kwargs)
        print("Success!\n")
        return

    def delete_init_fasta(self):
        print("delete_init_fasta:")
        shutil.move('external_scripts/pipeline/{}'
                    .format(self.CONST['input_seqs']),
                    self.CONST['trash'])
        shutil.move('external_scripts/meme/single_seq.fasta',
                    self.CONST['trash'])
        shutil.move('external_scripts/meme/{}'.format(self.CONST['input_seqs']),
                    self.CONST['trash'])
        print("Success!\n")
        return

    def delete_intermediate(self):
        print("delete_intermediate:")
        shutil.move('files/init_seed_seqs.fasta', self.CONST['trash'])
        shutil.move('files/meme_format.fasta', self.CONST['trash'])
        shutil.move('files/mast_single', self.CONST['trash'])
        shutil.move('files/meme_format2.txt', self.CONST['trash'])
        shutil.move('files/mast_nocorr', self.CONST['trash'])
        shutil.move('files/meme_format3.txt', self.CONST['trash'])
        shutil.move('files/mast_onlycombi', self.CONST['trash'])
        shutil.move('files/cluster_centroids.pkl', self.CONST['trash'])
        shutil.move('files/motifs', self.CONST['trash'])
        shutil.move(self.CONST['log'], self.CONST['trash'])
        print("Success!\n")
        return

    def reduce_dhcl(self):
        # Input: ./files/consensus_seqs.txt
        # Output: ./files/consensus_seqs.txt
        from reduce_dhcl_test import main
        kwargs = dict(consensus='files/consensus_seqs.txt',
                      output='files/consensus_seqs.txt')
        main(kwargs)
        print("Success!\n")
        return

    def shrink_input(self):
        # Input: ./files/input_seqs.fasta
        # Output: ./files/input_seqs.fasta
        from shrink_input_for_test import main
        kwargs = dict(seqs='files/input_seqs.txt',
                      output='files/input_seqs.txt',
                      divideby=100)
        main(kwargs)
        print("Success!\n")
        return


