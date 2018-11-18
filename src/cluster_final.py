from collections import Counter, defaultdict, OrderedDict
import re
import numpy as np
import difflib
from sklearn.cluster import AgglomerativeClustering
import pickle
import pandas as pd
import operator

def parse_mast_txt(screen=True, input_fname="../files/mast.txt"):
    # combi_count_map gives combi as key, number of sequences with that combi
    # as value
    # name_combi gives combi as key, list of name of sequences as value
    first_start = False
    second_start = False
    name_combi = defaultdict(str)
    with open(input_fname, 'r') as file:
        for line in file:
            if line.startswith("SECTION II"):
                first_start = True
                continue
            if first_start and line.startswith("-------------"):
                second_start = True
                continue
            if first_start and second_start:
                if line == "\n":
                    break
                name, remainder = line.split(" ", maxsplit=1)
                combi = re.findall("\[[0-9]+\]", remainder)
                combi_int = tuple([int(term[1:-1]) for term in combi])

                name_combi[combi_int] += name + " "
    name_combi = dict(name_combi)
    if screen:
        name_combi = screen_combi(name_combi, 5)
    return name_combi

def screen_combi(name_combi, threshold):
    to_del_combi = []
    for combi, seqs in name_combi.items():
        if len(seqs.split(" ")) < threshold:
            to_del_combi.append(combi)
    for combi in to_del_combi:
        del name_combi[combi]
    return name_combi

def lev_metric(x, y):
    sm = difflib.SequenceMatcher(None, x, y)
    return 1 - sm.ratio()

def get_dist_metric(combinations):
    # combinations = name_combi_map.keys()
    metric = []
    for curr in combinations:
        line = []
        for ref in combinations:
            line.append(lev_metric(curr, ref))
        metric.append(line)
    metric = np.array(metric)
    with open("dist_metric.pkl", "wb") as file:
        pickle.dump(metric, file, -1)
    return np.array(metric)

def cluster_metric(dist_metric, n_clusters):
    # Clustered combi, showing labels that can be matched as indices to the
    # combinations

    # Remember that each combi may have many members, so don't judge
    # based on number of combis per cluster, but members per cluster
    class_ = AgglomerativeClustering(n_clusters=n_clusters,
                                     affinity='precomputed',
                                     linkage='average').fit(dist_metric)
    cluster_labels = class_.labels_
    cluster_labels = np.array(cluster_labels)
    with open("cluster_labels.pkl", "wb") as file:
        pickle.dump(cluster_labels, file, -1)
    return cluster_labels

def get_index_clusters(labels):
#     Split labels into clusters, each with the index of the relevant combi
    _index_clusters = defaultdict(list)
    for i, label in enumerate(labels):
        _index_clusters[label].append(i)
    index_clusters = dict()
    for key, value in _index_clusters.items():
        index_clusters[key] = np.array(value)
    return index_clusters

# def get_cluster_combi(index_clusters, combinations):
#
#
#     return cluster_combi_map

def get_cluster_metrics(index_clusters, dist_metrics):

    return cluster_metric_map

def find_cluster_centroids(index_clusters, dist_metrics, combinations):

    cluster_combi_map = dict()
    # index_clusters has label: cluster_i, value: combi indices
    for label, index_cluster in index_clusters.items():
        cluster_combi_map[label] = combinations[index_cluster]
    # cluster_combi_map has key: cluster_i, value: [combis]
    # print(cluster_combi_map)
    cluster_metric_map = get_cluster_metrics(index_clusters, dist_metrics)

    cluster_centroids = dict()

    for label, metric in cluster_metric_map.items():
        min_score = 99999999999
        min_i = None
        for i, line in enumerate(metric):
            score = 0
            for j, element in enumerate(line):
                _combi = cluster_combi_map[label][j]
                count = len(name_combi_map[_combi])
                score += abs(element) * count
            if score < min_score:
                min_score = score
                min_i = i
        combi = cluster_combi_map[label][min_i]
        cluster_centroids[label] = combi
    with open("cluster_centroids.pkl", "wb") as file:
        pickle.dump(cluster_centroids, file, -1)
    return cluster_centroids

def get_name_clusters(index_clusters, name_combi_map):
    name_clusters = defaultdict(list)
    combinations = np.array(list(name_combi_map.keys()))
    for label, cluster_i in index_clusters.items():
        combi_in_cluster = combinations[cluster_i]
        for combi in combi_in_cluster:
            name_clusters[label].extend(name_combi_map[combi])
    return name_clusters

def get_cluster_allocation(name_clusters):
    cluster_allocation = dict()
    for label, names in name_clusters.items():
        _cluster_alloc = defaultdict(list)
        for name in names:
            group = re.match('[A-Z]+', name)[0][:-1]
            _cluster_alloc[group].append(name)
        _cluster_alloc = dict(_cluster_alloc)
        cluster_allocation[label] = _cluster_alloc
    return cluster_allocation

def find_clusters_too_few_members(cluster_allocation, threshold=30):
    labels_to_delete = []
    for label, clusters in cluster_allocation.items():
        tot_num = 0
        for i in clusters.values():
            tot_num += len(i)
        if tot_num < threshold:
            labels_to_delete.append(label)
    return labels_to_delete

# def assemble_clusters():
#     key: cluster centroid, value: seqs

import sys

def main_for_2nd_clustering():
    name_combi_map = parse_mast_txt()
    df = pd.DataFrame(index=range(len(name_combi_map)),
                      columns=('combi', 'seqs', 'dist_metric',
                               'cluster_label', 'num_seqs'))
    df['combi'] = list(name_combi_map.keys())
    df['seqs'] = list(name_combi_map.values())

    dist_metric = get_dist_metric(df['combi'])

    dist_metric_series = pd.Series()
    for i in range(len(dist_metric)):
        dist_metric_series.at[i] = dist_metric[i]
    df['dist_metric'] = dist_metric_series

    cluster_labels = cluster_metric(dist_metric, n_clusters=20)
    df['cluster_label'] = cluster_labels

    num_seqs = []
    for seq in df['seqs']:
        num_seqs.append(len(seq.split(" ")))
    df['num_seqs'] = num_seqs

    # Here, we spawn the new df
    cluster_df = pd.DataFrame(index=set(cluster_labels),
                              columns=['num_seqs', 'centroid', 'seq_alloc'])

    for cluster_label, ind_df in df.groupby('cluster_label'):
        num_seqs = sum(ind_df['num_seqs'])
        cluster_df.at[cluster_label, 'num_seqs'] = num_seqs
        min_score = None
        min_i = None
        for combi_i, full_dist in ind_df['dist_metric'].iteritems():
            reduced_dist = full_dist[ind_df.index]
            score = sum(reduced_dist * ind_df['num_seqs'])
            if min_score is None or score < min_score:
                min_score = score
                min_i = combi_i
        centroid = ind_df.loc[min_i, 'combi']
        cluster_df.at[cluster_label, 'centroid'] = centroid
        seq_alloc_unsorted = defaultdict(int)
        for merged_seq in ind_df['seqs']:
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
        cluster_df.at[cluster_label, 'seq_alloc'] = seq_alloc

    to_drop = []
    for cluster_i, ind_df in cluster_df.iterrows():
        if ind_df['num_seqs'] < 50:
            to_drop.append(cluster_i)

    cluster_df = cluster_df.drop(to_drop, axis='index')

    with open("../files/df.pkl", 'wb') as file:
        pickle.dump(df, file, -1)
    with open("../files/cluster_df.pkl", 'wb') as file:
        pickle.dump(cluster_df, file, -1)

    with open("../output/cluster_description.txt", 'w') as file:
        for cluster_label, ind_df in cluster_df.iterrows():
            file.write("Cluster {}\n".format(cluster_label))
            file.write("Combination: {}\n".format(ind_df['centroid']))
            for group_name, seq_len in ind_df['seq_alloc'].items():
                file.write("{} : {}\n".format(group_name, seq_len))
            file.write("\n\n")

if __name__ == "__main__":
    main_for_2nd_clustering()


# def sort_combi_data(combi_data):
#     # combi_data is now sorted by the number of sequences for that combi
#     to_sort = []
#     for key, value in combi_data.items():
#         to_sort.append((value, key))
#     to_sort.sort(key=lambda x: x[0], reverse=True)
#     to_sort = np.array([element[1] for element in to_sort])
#     return to_sort
#
# def create_set_of_conserved_profiles(cluster_centroids):
#     profiles = set()
#     for combi in cluster_centroids.values():
#         for i in combi:
#             profiles.add(i)
#     return profiles


    # to_trim = np.array([True] * len(df), dtype=bool)
    # for c_label, ind_df in df.groupby('cluster_label'):
    #     tot_seq_count = sum(ind_df['num_seqs'])
    #     if tot_seq_count < 50:
    #         # df = df[df.cluster_label != c_label]
    #         to_trim = to_trim & df.cluster_label != c_label
    #
    # df = df[to_trim]
    #
    # for i, ind_row in df.iterrows():
    #     df.at[i, 'dist_metric'] = df.loc[i]['dist_metric'][to_trim]
    #
    # assert len(df) == len(df['dist_metric'])
    # assert len(df['dist_metric']) == len(df['dist_metric'][0])
    #
    # df = df.assign(score=None)
    #
    # for cluster_label, cluster_df in df.groupby('cluster_label'):
    #     bool_mask = df['cluster_label'] == cluster_label
    #     for i, ind_row in cluster_df.iterrows():
    #         assert not np.any(ind_row['dist_metric'] < 0)
    #         # This should only be multiplied with those from the same cluster
    #         score = sum(ind_row['dist_metric'][bool_mask] * df['num_seqs'][bool_mask])
    #         df.loc[i, 'score'] = score

    # for i, ind_row in df.iterrows():
    #     assert not np.any(ind_row['dist_metric'] < 0)
    #     # This should only be multiplied with those from the same cluster
    #     score = sum(ind_row['dist_metric'] * df['num_seqs'])
    #     df.loc[i, 'score'] = score

    # centroid_series = pd.Series()
    # for cluster_label, cluster_df in df.groupby('cluster_label'):
    #     indices, scores = cluster_df['score'].index, cluster_df['score'].values
    #     min_i = indices[np.argmin(scores)]
    #     centroid = df.loc[min_i, 'combi']
    #     for i in indices:
    #         centroid_series.at[i] = centroid
    #
    # df['centroid'] = centroid_series
    # print(df.columns)
    # todo: it might make sense to spawn another df, specifically for
    # extracting info from the previous, larger df.

    # # cluster_centroid_table = list(cluster_centroids.values())
    # # profiles output
    #
    # name_clusters = get_name_clusters(index_clusters, name_combi_map)
    #
    #
    # cluster_allocation = get_cluster_allocation(name_clusters)
    # labels_to_delete = find_clusters_too_few_members(cluster_allocation)
    # for label in labels_to_delete:
    #     del cluster_centroids[label]
    #     del cluster_allocation[label]
    #
    # cluster_allocation = OrderedDict(
    #     sorted(cluster_allocation.items(), key=lambda t: t[0]))
    # cluster_centroids = OrderedDict(
    #     sorted(cluster_centroids.items(), key=lambda t: t[0]))