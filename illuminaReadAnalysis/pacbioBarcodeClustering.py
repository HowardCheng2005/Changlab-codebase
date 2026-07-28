import argparse
import json

import pandas as pd

def get_pacbio_dict(input_file):
    with open(input_file, "r") as f:
        pacbio_dict = json.load(f)
    return pacbio_dict

def barcode_grouping(input_file, barcode_dict):
    df = pd.read_csv(input_file)

    cluster_ids = df["Cluster.ID"]
    cluster_sizes = df["time_point_1"]

    cluster_dict = {}

    for barcode, value in barcode_dict.items():
        barcode_id = value[1]

        if barcode_id not in cluster_dict:
            cluster_dict[barcode_id] = [[barcode]]
        else:
            cluster_dict[barcode_id][0].append(barcode)

    for i, id in enumerate(cluster_ids):
        if id in cluster_dict:
            cluster_dict[id].append(cluster_sizes[i])

    return cluster_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cluster", required=True, help="bartender combined cluster .csv file")
    parser.add_argument("-p", "--pacbioCluster", required=True, help="pacbio barcode cluster association .json file")
    parser.add_argument("-o", "--output", required=True, help="output for cluster grouping of Pacbio barcodes in .json format")

    args = parser.parse_args()

    cluster_file = args.cluster
    pacbio_file = args.pacbioCluster
    output_file = args.output
    summary_file = args.summary

    pacbio_dict = get_pacbio_dict(pacbio_file)

    clustered_barcodes = barcode_grouping(cluster_file, pacbio_dict)

    with open(output_file, 'w') as o:
        json.dump(clustered_barcodes, o)