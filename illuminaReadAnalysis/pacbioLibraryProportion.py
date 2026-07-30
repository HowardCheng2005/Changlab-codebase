from argparse import ArgumentParser
import pandas as pd
import json

def pacbio_dictionary_creator(pacbio_file):
    with open(pacbio_file) as f:
        pacbio_dictionary = json.load(f)
        print(f"dictionary length: {len(pacbio_dictionary.keys())}")

        return pacbio_dictionary

def pacbio_cluster_dictionary_creator(pacbio_dict):
    cluster_dictionary = {}

    for _, values in pacbio_dict.items():
        for value in values[1]:
            if value not in cluster_dictionary:
                cluster_dictionary[value] = 1

    return cluster_dictionary

def check_library_read_proportion(barcode_file, pacbio_dict, cluster_dict):
    df = pd.read_csv(barcode_file)

    num_unique_library_reads = 0
    num_unique_without_pacbio = 0
    num_total_library_reads = 0
    num_total_without_pacbio = 0

    reads = df['Unique.reads']
    frequency = df['Frequency']
    clusters = df['Cluster.ID']

    print(f"total reads: {print(sum(frequency))}")

    for i, read in enumerate(reads):
        if read not in pacbio_dict:
            num_unique_library_reads += 1
            num_total_library_reads += frequency[i]

            if clusters[i] not in cluster_dict:
                num_unique_without_pacbio += 1
                num_total_without_pacbio += frequency[i]

    return num_unique_library_reads, num_unique_without_pacbio, num_total_library_reads, num_total_without_pacbio


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True, help="input pacbio to cluster .json file")
    parser.add_argument("-b", "--barcode", dest="barcode", required=True, help="barcode file in .csv format containing pacbio + library reads")
    parser.add_argument("-o", "--output", dest="output", required=True, help="output file for clusters with no pacbio")


    args = parser.parse_args()

    input = args.input
    barcode = args.barcode
    output = args.output

    pacbio_dictionary = pacbio_dictionary_creator(input)
    cluster_dictionary = pacbio_cluster_dictionary_creator(pacbio_dictionary)

    unique_reads, unique_without_pacbio_reads, total_reads, total_without_pacbio_reads = check_library_read_proportion(barcode, pacbio_dictionary, cluster_dictionary)

    with open(output, "w") as f:
        f.write(f"Unique reads: {unique_reads}\n")
        f.write(f"Unique no pacbio: {unique_without_pacbio_reads}\n")
        f.write(f"Total reads: {total_reads}\n")
        f.write(f"Total without pacbio: {total_without_pacbio_reads}\n")
        f.close()
