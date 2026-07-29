import argparse
import json

import pandas as pd
import matplotlib.pyplot as plt


# reads pacbio dictionary
def get_pacbio_dict(input_file):
    with open(input_file, "r") as f:
        pacbio_dict = json.load(f)
    return pacbio_dict

# groups by barcode in barcode: [[barcode_list], size of cluster]
def barcode_grouping(input_file, barcode_dict):
    df = pd.read_csv(input_file)

    cluster_ids = df["Cluster.ID"].astype(int)
    cluster_sizes = df["time_point_1"].astype(int)

    cluster_dict = {}

    for barcode, value in barcode_dict.items():

        for cluster_id in value[1]:
            cluster_id = int(cluster_id)

            if cluster_id not in cluster_dict:
                cluster_dict[cluster_id] = [[barcode]]
            else:
                cluster_dict[cluster_id][0].append(barcode)

    for i, cluster_id in enumerate(cluster_ids):
        cluster_id = int(cluster_id)

        if cluster_id in cluster_dict:
            cluster_dict[cluster_id].append(int(cluster_sizes.iloc[i]))

    return cluster_dict

def cluster_analysis(cluster_dict):
    cluster_size = []
    pacbio_in_cluster = []

    for cluster_id, value in cluster_dict.items():
        cluster_size.append(value[1])
        pacbio_in_cluster.append(len(value[0]))

    df = pd.DataFrame({
        "cluster_size": cluster_size,
        "unique_pacbio_barcodes": pacbio_in_cluster
    })

    summary = df.describe()

    plt.figure()
    plt.hist(df, bins=num_bins)

    plt.xlabel("Size of cluster/number of pacbio reads in clusters")
    plt.ylabel("Number of clusters")
    plt.title("Histogram of size per cluster with pacbio reads included")

    plt.tight_layout()

    return plt, summary

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cluster", required=True, help="bartender combined cluster .csv file")
    parser.add_argument("-p", "--pacbioCluster", required=True, help="pacbio barcode cluster association .json file")
    parser.add_argument("-o", "--output", required=True, help="output for cluster grouping of Pacbio barcodes in .json format")
    parser.add_argument("-s", "--summary", required=True, help="output for summary of the clustering process in .txt format")
    parser.add_argument("-hi", "--histogram", required=True, help="output for histogram of the pacbio cluster sizes in .png format")
    parser.add_argument("-b", "--bins", required=True, type=int, help="number of bins for histogram plot")

    args = parser.parse_args()

    cluster_file = args.cluster
    pacbio_file = args.pacbioCluster
    output_file = args.output
    summary_file = args.summary
    histogram_file = args.histogram
    num_bins = args.bins

    pacbio_dict = get_pacbio_dict(pacbio_file)

    clustered_barcodes = barcode_grouping(cluster_file, pacbio_dict)

    plt, summary = cluster_analysis(clustered_barcodes)

    plt.savefig(histogram_file, dpi=300)
    plt.close()

    with open(summary_file, 'w') as s:
        s.write((summary.to_string()))

    with open(output_file, 'w') as o:
        json.dump(clustered_barcodes, o)