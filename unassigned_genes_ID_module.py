import pandas as pd
import argparse
import os

def process_unassigned_genes(file_path, species_column, output_file):
    """
    Processes the Orthogroups_UnassignedGenes.tsv file to filter rows with unassigned genes.

    Args:
        file_path (str): Path to the Orthogroups_UnassignedGenes.tsv file.
        species_column (str): Name of the species column (e.g., "Melipona_bicolor").
        output_file (str): Path to save the filtered unassigned genes.

    Returns:
        None
    """
    # Load the dataset
    unassigned_genes = pd.read_csv(file_path, sep="\t", header=0)

    # Filter rows where the species column is not NA
    filtered_data = unassigned_genes[~unassigned_genes[species_column].isna()]

    # Extract only the 'Orthogroup' column and drop NA values
    orthogroup_ids = filtered_data[['Orthogroup']].dropna()

    # Append ".fa" to each entry in the Orthogroup column
    orthogroup_ids['Orthogroup'] = orthogroup_ids['Orthogroup'] + ".fa"

    # Save to a text file without headers or index
    orthogroup_ids.to_csv(output_file, index=False, header=False, sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter unassigned genes for a specific species")
    parser.add_argument("--file_path", type=str, required=True, help="Path to the Orthogroups_UnassignedGenes.tsv file")
    parser.add_argument("--species_column", type=str, required=True, help="Name of the species column")
    parser.add_argument("--output_file", type=str, required=True, help="Path to save the filtered unassigned genes")

    args = parser.parse_args()

    process_unassigned_genes(args.file_path, args.species_column, args.output_file)
    print(f"Output saved to: {args.output_file}")
    print(f"Current Working Directory: {os.getcwd()}")
