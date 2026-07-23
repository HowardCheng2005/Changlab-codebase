import argparse
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="input cluster.csv file")
    parser.add_argument("-o", "--output", required=True, help="output for distribution in .png format histogram")
    parser.add_argument("-c", "--csv_output", required=True, help="output for singleton distribution in .csv format")
    parser.add_argument("-d", "--describe_output", required=True, help="output for a summary of the distribution in .txt format")
    parser.add_argument("-b", "--bins", required=True, type=int, help="number of bins to use")

    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    csv_output = args.csv_output
    summary_output = args.describe_output
    num_bins = args.bins

    df = pd.read_csv(input_file)
    cluster_sizes = df["time_point_1"]

    read_counts = (cluster_sizes.value_counts().rename_axis("reads_per_cluster").reset_index(name="number_of_clusters"))
    read_counts = read_counts.sort_values(by=["reads_per_cluster"], ascending=False)

    print(f"Describe the set of the cluster sizes:")
    summary = cluster_sizes.describe()
    print(summary)
    with open(summary_output, "w") as file:
        file.write(summary.to_string())

    plt.figure()
    plt.hist(cluster_sizes, bins=num_bins)

    plt.xlabel("Number of reads in cluster")
    plt.ylabel("Number of clusters")
    plt.title("Histogram of reads per cluster")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    read_counts.to_csv(csv_output, index=False)