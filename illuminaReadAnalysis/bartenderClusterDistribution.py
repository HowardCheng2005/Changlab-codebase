import argparse
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="input cluster .csv file")
    parser.add_argument("-o", "--output", required=True, help="output for distribution in .png format")
    parser.add_argument("-c", "--csv_output", required=True, help="output for distribution in .csv format")
    args = parser.parse_args()

    input_file = args.parse_args().input
    output_file = args.parse_args().output
    csv_output = args.parse_args().output

    df = pd.read_csv(input_file)

    read_counts = (df["time_point_1"].value_counts().rename_axis("reads_per_cluster").reset_index(name="number_of_clusters"))
    read_counts = read_counts.sort_values(by=["reads_per_cluster"], ascending=False)

    print(f"Distribution of top 50 largest clusters:")
    print(read_counts.head(50))

    plt.figure()
    plt.bar(read_counts["reads_per_cluster"], read_counts["number_of_clusters"], align='center')

    plt.xlabel("number of reads in cluster")
    plt.ylabel("number of clusters")
    plt.title("Distribution of reads in clusters")

    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

    read_counts.to_csv(csv_output, index=False)