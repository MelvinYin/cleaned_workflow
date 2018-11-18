from collections import defaultdict, OrderedDict
import re
import numpy as np
import difflib
from sklearn.cluster import AgglomerativeClustering
import pickle
import pandas as pd
import operator
import sys
#
# def parse_mast_txt(input_fname="../files/mast.txt", screen_threshold=5):
#     first_start = False
#     second_start = False
#     name_combi = defaultdict(str)
#     with open(input_fname, 'r') as file:
#         for line in file:
#             if line.startswith("SECTION II"):
#                 first_start = True
#                 continue
#             if first_start and line.startswith("-------------"):
#                 second_start = True
#                 continue
#             if first_start and second_start:
#                 if line == "\n":
#                     break
#                 name, remainder = line.split(" ", maxsplit=1)
#                 combi = re.findall("\[[0-9]+\]", remainder)
#                 combi_int = tuple([int(term[1:-1]) for term in combi])
#
#                 name_combi[combi_int] += name + " "
#     name_combi = dict(name_combi)
#     # if screen_threshold:
#     #     name_combi = screen_combi(name_combi, screen_threshold)
#     return name_combi
#
# def screen_combi(name_combi, threshold):
#     to_del_combi = []
#     for combi, seqs in name_combi.items():
#         if len(seqs.split(" ")) < threshold:
#             to_del_combi.append(combi)
#             print(combi)
#             print(seqs)
#     for combi in to_del_combi:
#         del name_combi[combi]
#     return name_combi
#
# def lev_metric(x, y):
#     sm = difflib.SequenceMatcher(None, x, y)
#     return 1 - sm.ratio()
#
# def get_dist_metric(combinations):
#     metric = []
#     for curr in combinations:
#         line = []
#         for ref in combinations:
#             line.append(lev_metric(curr, ref))
#         metric.append(line)
#     metric = np.array(metric)
#     return np.array(metric)
#
# def cluster_metric(dist_metric, n_clusters):
#     # Clustered combi, showing labels that can be matched as indices to the
#     # combinations
#
#     # Remember that each combi may have many members, so don't judge
#     # based on number of combis per cluster, but members per cluster
#     class_ = AgglomerativeClustering(n_clusters=n_clusters,
#                                      affinity='precomputed',
#                                      linkage='average').fit(dist_metric)
#     cluster_labels = class_.labels_
#     return cluster_labels

def get_cluster_centroid(ind_df):
    min_score = None
    min_i = None
    for combi_i, full_dist in ind_df['dist_metric'].iteritems():
        reduced_dist = full_dist[ind_df.index]
        score = sum(reduced_dist * ind_df['num_seqs'])
        if min_score is None or score < min_score:
            min_score = score
            min_i = combi_i
    centroid = ind_df.loc[min_i, 'combi']
    return centroid

def get_seq_alloc_sorted(seqs):
    seq_alloc_unsorted = defaultdict(int)
    for merged_seq in seqs:
        for seq in merged_seq.split(" "):
            match_obj = re.match("([A-Z]+)Uniprot", seq)
            if match_obj is None:
                continue
            group = match_obj.group(1)
            seq_alloc_unsorted[group] += 1
    seq_alloc = OrderedDict()
    for key, value in sorted(seq_alloc_unsorted.items(),
                             key=operator.itemgetter(1), reverse=True):
        seq_alloc[key] = value
    return seq_alloc

def main(kwargs):
    cluster_threshold = kwargs['cluster_threshold']
    pkl_path = kwargs['pkl_path']
    output = kwargs['output']

    # Load cluster_params
    with open(pkl_path, 'rb') as file:
        full_df = pickle.load(file)

    # Assemble cluster_df
    cluster_df = pd.DataFrame(index=set(full_df['cluster_label']),
                              columns=['num_seqs', 'centroid', 'seq_alloc'])
    for cluster_label, ind_df in full_df.groupby('cluster_label'):
        # Get num_seqs
        num_seqs = sum(ind_df['num_seqs'])
        cluster_df.at[cluster_label, 'num_seqs'] = num_seqs

        # Get centroid
        centroid = get_cluster_centroid(ind_df)
        cluster_df.at[cluster_label, 'centroid'] = centroid

        # Get seq_alloc
        seq_alloc = get_seq_alloc_sorted(ind_df['seqs'])
        cluster_df.at[cluster_label, 'seq_alloc'] = seq_alloc

    # Drop small clusters
    to_drop = []
    for cluster_i, ind_df in cluster_df.iterrows():
        if ind_df['num_seqs'] < cluster_threshold:
            to_drop.append(cluster_i)
    cluster_df = cluster_df.drop(to_drop, axis='index')

    # Output clusters as txt
    with open(output, 'w') as file:
        for cluster_label, ind_df in cluster_df.iterrows():
            file.write("Cluster {}\n".format(cluster_label))
            file.write("Combination: {}\n".format(ind_df['centroid']))
            for group_name, seq_len in ind_df['seq_alloc'].items():
                file.write("{} : {}\n".format(group_name, seq_len))
            file.write("\n\n")
    return True

if __name__ == "__main__":
    kwargs = dict(cluster_threshold=50,
                  pkl_path="../files/clustering_df.pkl",
                  output="../files/cluster_description.txt")
    main(kwargs)
