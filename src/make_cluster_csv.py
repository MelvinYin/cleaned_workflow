import pickle


def main():
    with open("files/cluster_alloc.pkl", 'wb') as file:
        pickle.dump(cluster_allocation, file, -1)
