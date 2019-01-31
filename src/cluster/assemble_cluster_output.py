from collections import defaultdict, OrderedDict
import re
import pickle
import pandas as pd
import operator
import sys

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

def get_seq_alloc(seqs):
    seq_alloc_unsorted = defaultdict(int)
    for merged_seq in seqs:
        for seq in merged_seq.split(" "):
            match_obj = re.match("([A-Z]+)Uniprot", seq)
            if match_obj is None:
                continue
            group = match_obj.group(1)
            seq_alloc_unsorted[group] += 1
    seq_alloc = OrderedDict()
    for group, num_seq_in_group in \
            sorted(seq_alloc_unsorted.items(), key=operator.itemgetter(1),
                   reverse=True):
        seq_alloc[group] = num_seq_in_group
    return seq_alloc

def main(kwargs):
    cluster_threshold = kwargs['cluster_threshold']
    output = kwargs['output']
    full_param_pkl = kwargs['full_param_pkl']
    cluster_df_pkl = kwargs['cluster_df_pkl']

    # Load cluster_params
    with open(full_param_pkl, 'rb') as file:
        full_param_df = pickle.load(file)

    # Assemble cluster_df
    cluster_df = pd.DataFrame(index=set(full_param_df['cluster_label']),
                              columns=['num_seqs', 'centroid', 'seq_alloc'])
    for cluster_label, df_per_clust in full_param_df.groupby('cluster_label'):
        # Get num_seqs
        clust_num_seq = sum(df_per_clust['num_seqs'])
        cluster_df.at[cluster_label, 'num_seqs'] = clust_num_seq

        # Get centroid
        centroid = get_cluster_centroid(df_per_clust)
        cluster_df.at[cluster_label, 'centroid'] = centroid

        # Get seq_alloc
        # seq_alloc may be less than what num_seqs indicate, because some
        # lines are dropped, from the "Uniprot" screen.
        seq_allocation = get_seq_alloc(df_per_clust['seqs'])
        cluster_df.at[cluster_label, 'seq_alloc'] = seq_allocation

    # Drop small clusters
    to_drop = []
    for cluster_i, ind_df in cluster_df.iterrows():
        if ind_df['num_seqs'] < cluster_threshold:
            to_drop.append(cluster_i)
    cluster_df = cluster_df.drop(to_drop, axis='index')
    if cluster_df_pkl:
        with open(cluster_df_pkl, 'wb') as file:
            pickle.dump(cluster_df, file, -1)

    # Output clusters as txt
    if output:
        with open(output, 'w') as file:
            for cluster_label, ind_df in cluster_df.iterrows():
                file.write(f"Cluster {cluster_label}\n")
                file.write(f"Combination: {ind_df['centroid']}\n")
                for group_name, seq_len in ind_df['seq_alloc'].items():
                    file.write("{} : {}\n".format(group_name, seq_len))
                file.write("\n\n")
    return True
