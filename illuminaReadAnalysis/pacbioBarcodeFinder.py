import argparse
import json
import pandas as pd


def pacbio_barcode_reader(pacbio_file):
    with open(pacbio_file, 'r') as f:
        barcode_dict = {}

        for line in f:
            line = line.strip()

            if not line:
                continue

            parts = line.split()
            frequency = float(parts[0])
            sequence = parts[1]
            if sequence not in barcode_dict:
                barcode_dict[sequence] = [frequency, []]

    return barcode_dict

def pacbio_cluster_finder(pacbio_dict, input_file):
    df = pd.read_csv(input_file)

    reads = df["Unique.reads"]
    clusters = df["Cluster.ID"]

    for i, read in enumerate(reads):
        if read in pacbio_dict:
            pacbio_dict[read][1].append(int(clusters[i]))

def count_matched_barcodes(pacbio_dict):
    matched_barcodes = 0

    for read in pacbio_dict.values():
        if len(read[1]) != 0:
            matched_barcodes += 1

    return matched_barcodes


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bartender_input", required=True, help="input barcode .csv file from bartender")
    parser.add_argument("-p", "--pacbio_input", required=True, help="input Pacbio reads in .txt format")
    parser.add_argument("-o", "--output", required=True, help="output file in .json format")
    parser.add_argument("-s", "--summary", required=True, help="summary file for output in .txt format")

    args = parser.parse_args()

    bartender_file = args.bartender_input
    pacbio_file = args.pacbio_input
    output_file = args.output
    summary_file = args.summary

    pacbio_barcodes = pacbio_barcode_reader(pacbio_file)
    pacbio_cluster_finder(pacbio_barcodes, bartender_file)

    matched_barcodes = count_matched_barcodes(pacbio_barcodes)
    unmatched_barcodes = len(pacbio_barcodes) - matched_barcodes

    with open(output_file, 'w') as f:
        json.dump(pacbio_barcodes, f)

    with open(summary_file, 'w') as f:
        f.write(f"Total Pacbio barcodes: {len(pacbio_barcodes)}\n")
        f.write(f"Matched barcodes: {matched_barcodes}\n")
        f.write(f"Unmatched barcodes: {unmatched_barcodes}\n")

    print(f"Total PacBio barcodes: {len(pacbio_barcodes):,}")
    print(f"Matched barcodes: {matched_barcodes:,}")
    print(f"Unmatched barcodes: {unmatched_barcodes:,}")
    print(f"Output written to: {args.output}")









