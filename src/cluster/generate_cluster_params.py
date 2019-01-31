from collections import defaultdict
import difflib
import numpy as np
import pandas as pd
import pickle
import re
from sklearn.cluster import AgglomerativeClustering

def parse_mast_txt(input_fname):
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
    return name_combi

def _lev_metric(x, y):
    sm = difflib.SequenceMatcher(None, x, y)
    return 1 - sm.ratio()

def get_dist_metric(combinations):
    metric = []
    for curr in combinations:
        line = []
        for ref in combinations:
            line.append(_lev_metric(curr, ref))
        metric.append(line)
    metric = np.array(metric)
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
    return cluster_labels

def main(kwargs):
    input_mast = kwargs['input_mast']
    combi_minsize = kwargs['combi_minsize']
    pkl_path = kwargs['pkl_path']
    num_cluster = kwargs['num_cluster']
    # Assemble full_df
    full_df = pd.DataFrame(columns=('combi', 'seqs', 'dist_metric',
                               'cluster_label', 'num_seqs'))
    name_combi_map = parse_mast_txt(input_mast)
    full_df['combi'] = list(name_combi_map.keys())
    full_df['seqs'] = list(name_combi_map.values())

    for i, seq in full_df['seqs'].iteritems():
        full_df.loc[i, 'num_seqs'] = len(seq.split(" "))

    # Drop small combi
    to_drop = []
    for i, num_seq in full_df['num_seqs'].iteritems():
        if num_seq < combi_minsize:
            to_drop.append(i)
    full_df = full_df.drop(to_drop, axis='index').reset_index()

    dist_metric = get_dist_metric(full_df['combi'])
    for i in range(len(dist_metric)):
        full_df.at[i, 'dist_metric'] = dist_metric[i]

    cluster_labels = cluster_metric(dist_metric, n_clusters=num_cluster)
    full_df['cluster_label'] = cluster_labels

    # Dump
    with open(pkl_path, 'wb') as file:
        pickle.dump(full_df, file, -1)
    return True

if __name__ == "__main__":
    kwargs = dict(input_mast = "../files/mast.txt",
                  screen_threshold = 5,
                  pkl_path="../files/clustering_df.pkl",
                  num_cluster=20)
    main(kwargs)