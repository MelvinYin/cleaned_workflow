from collections import Counter, defaultdict, OrderedDict
import re
import numpy as np
import difflib
from sklearn.cluster import AgglomerativeClustering
import pickle

USE_SAVE = False

def parse_mast_txt(screen=True, input_fname="files/mast.txt"):
    # combi_count_map gives combi as key, number of sequences with that combi
    # as value
    # name_combi gives combi as key, list of name of sequences as value
    first_start = False
    second_start = False
    combi_data = Counter()
    name_combi = defaultdict(list)
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
                combi_data[combi_int] += 1
                name_combi[combi_int].append(name)
    name_combi = dict(name_combi)
    if screen:
        combi_data, name_combi = screen_combi(combi_data, name_combi)
    return combi_data, name_combi

def screen_combi(combi_data, name_combi):
    to_del_combi = []
    for combi, count in combi_data.items():
        if count < 5:
            to_del_combi.append(combi)
    for combi in to_del_combi:
        del combi_data[combi]
        del name_combi[combi]
    return combi_data, name_combi

def sort_combi_data(combi_data):
    # combi_data is now sorted by the number of sequences for that combi
    to_sort = []
    for key, value in combi_data.items():
        to_sort.append((value, key))
    to_sort.sort(key=lambda x: x[0], reverse=True)
    to_sort = np.array([element[1] for element in to_sort])
    return to_sort

def OUTPUT2_get_raw_combinations(combi_count_map):
    return list(combi_count_map.keys())

def lev_metric(x, y):
    sm = difflib.SequenceMatcher(None, x, y)
    return 1 - sm.ratio()

def get_dist_metric(name_combi_map, skip=True):
    try:
        if not USE_SAVE or not skip:
            raise Exception
        with open("dist_metric.pkl", "rb") as file:
            metric = pickle.load(file)
        return metric
    except:
        combinations = name_combi_map.keys()
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

def cluster_metric(dist_metric, n_clusters=20):
    # Clustered combi, showing labels that can be matched as indices to the
    # combinations
    try:
        if not USE_SAVE:
            raise Exception
        with open("cluster_labels.pkl", "rb") as file:
            cluster_labels = pickle.load(file)
        return cluster_labels
    except:
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

def get_cluster_combi(index_clusters, combinations):
    cluster_combi_map = dict()
    for label, index_cluster in index_clusters.items():
        cluster_combi_map[label] = combinations[index_cluster]
    return cluster_combi_map

def get_cluster_metrics(index_clusters, dist_metrics):
    cluster_metric_map = dict()
    for label, index_cluster in index_clusters.items():
        metric = np.array([line[index_cluster] for line in dist_metrics[index_cluster]])
        cluster_metric_map[label] = metric
    return cluster_metric_map

def find_cluster_centroids(index_clusters, combinations, dist_metrics, combi_count_map):
    try:
        if not USE_SAVE:
            raise Exception
        with open("cluster_centroids.pkl", "rb") as file:
            cluster_centroids = pickle.load(file)
        return cluster_centroids
    except:
        cluster_combi_map = get_cluster_combi(index_clusters, combinations)
        cluster_metric_map = get_cluster_metrics(index_clusters, dist_metrics)
        cluster_centroids = dict()

        for label, metric in cluster_metric_map.items():
            min_score = 99999999999
            min_i = None
            for i, line in enumerate(metric):
                score = 0
                for j, element in enumerate(line):
                    _combi = cluster_combi_map[label][j]
                    count = combi_count_map[_combi]
                    score += abs(element) * count
                if score < min_score:
                    min_score = score
                    min_i = i
            combi = cluster_combi_map[label][min_i]
            cluster_centroids[label] = combi
        with open("cluster_centroids.pkl", "wb") as file:
            pickle.dump(cluster_centroids, file, -1)
    return cluster_centroids

def create_set_of_conserved_profiles(cluster_centroids):
    profiles = set()
    for combi in cluster_centroids.values():
        for i in combi:
            profiles.add(i)
    return profiles

def OUTPUT_cluster_centroid_table(cluster_centroids):
    combi = list(cluster_centroids.values())
    return combi

def OUTPUT_profiles(profiles):
    return profiles

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

def main_for_2nd_clustering():
    combi_count_map, name_combi_map = parse_mast_txt()
    OUTPUT2_get_raw_combinations(combi_count_map)
    # combi_count_map = sort_combi_data(combi_count_map)
    # Next, we try to get the cluster centroids
    # From there, we can get the set of profiles, and then rebuild combi using
    # mast, maybe.
    dist_metric = get_dist_metric(name_combi_map, skip=False)
    cluster_labels = cluster_metric(dist_metric, n_clusters=20)
    # With set of cluster labels, we should form a name+combi+metric for each
    # cluster first. Metric is needed for centroid formation, combi to use to
    # find centroid profiles
    # print(cluster_labels)
    index_clusters = get_index_clusters(cluster_labels)
    combinations = np.array(list(combi_count_map.keys()))

    cluster_centroids = find_cluster_centroids(index_clusters, combinations,
                                               dist_metric, combi_count_map)
    profiles = create_set_of_conserved_profiles(cluster_centroids)
    OUTPUT_cluster_centroid_table(cluster_centroids)
    OUTPUT_profiles(profiles)
    name_clusters = get_name_clusters(index_clusters, name_combi_map)

    cluster_allocation = get_cluster_allocation(name_clusters)
    labels_to_delete = find_clusters_too_few_members(cluster_allocation)
    for label in labels_to_delete:
        del cluster_centroids[label]
        del cluster_allocation[label]

    cluster_allocation = OrderedDict(
        sorted(cluster_allocation.items(), key=lambda t: t[0]))
    cluster_centroids = OrderedDict(
        sorted(cluster_centroids.items(), key=lambda t: t[0]))

    with open("files/cluster_centroids.pkl", 'wb') as file:
        pickle.dump(cluster_centroids, file, -1)

    with open("output/cluster_description.txt", 'w') as file:
        for label, clusters in cluster_allocation.items():
            file.write("Cluster {}\n".format(label))
            file.write("Combination: [")
            for i in cluster_centroids[label]:
                file.write("{} ".format(i))
            file.write("]\n")
            name_cluster = []
            for group_name, seqs in clusters.items():
                name_cluster.append([group_name, len(seqs)])
            name_cluster.sort(key=lambda x: x[1])
            name_cluster = name_cluster[::-1]
            for group_name, seq_len in name_cluster:
                file.write("{} : {}\n".format(group_name, seq_len))
            file.write("\n\n")



    # profiles = create_set_of_conserved_profiles(cluster_centroids)
    # for label, clusters in cluster_allocation.items():
    #     print(clusters.keys())
    #     print(cluster_centroids[label])
    #     for i in clusters.values():
    #         print(len(i))
    #     print("")

# TODO: First, dump cluster_allocation and possibly profiles in separate pkl
# files. We then need a script that takes cluster_allocation, and build
# profile logo folders using ceqlogo for each cluster. Some form of labelling
#  (csv?) can be used to describe the clusters.

if __name__ == "__main__":
    main_for_2nd_clustering()
