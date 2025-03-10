import pandas as pd
import argparse
import os

def process_orthogroups(file_path, species_column, output_file):
    """
    Processes the Orthogroups.GeneCount.tsv file to filter species-specific orthogroups.

    Args:
        file_path (str): Path to the Orthogroups.GeneCount.tsv file.
        species_column (str): Name of the species column (e.g., "Melipona_bicolor").
        output_file (str): Path to save the filtered orthogroups.

    Returns:
        None
    """
    # Load the dataset
    orthogroups_gene_count = pd.read_csv(file_path, sep="\t", header=0)

    # Filter rows where species-specific column equals the Total column
    species_specific_orthogroups = orthogroups_gene_count[
        orthogroups_gene_count[species_column] == orthogroups_gene_count['Total']
    ]

    # Extract only the 'Orthogroup' column and append ".fa" to each entry
    orthogroup_ids = species_specific_orthogroups[['Orthogroup']].copy()
    orthogroup_ids['Orthogroup'] = orthogroup_ids['Orthogroup'] + ".fa"

    # Save to a text file without headers or index
    orthogroup_ids.to_csv(output_file, index=False, header=False, sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter species-specific orthogroups")
    parser.add_argument("--file_path", type=str, required=True, help="Path to the Orthogroups.GeneCount.tsv file")
    parser.add_argument("--species_column", type=str, required=True, help="Name of the species column")
    parser.add_argument("--output_file", type=str, required=True, help="Path to save the filtered orthogroups")

    args = parser.parse_args()

    process_orthogroups(args.file_path, args.species_column, args.output_file)
    print(f"Output saved to: {args.output_file}")
    print(f"Current Working Directory: {os.getcwd()}")
