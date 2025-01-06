#!/usr/bin/env python3
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Concatenating the files into HPP_host_prediction.txt
# Date: 20250106
#!/usr/bin/env python3
import pandas as pd

def preprocess_file_inplace(file_path):
    # Read and process the file
    with open(file_path, 'r') as infile:
        lines = infile.readlines()
    
    with open(file_path, 'w') as outfile:
        for line in lines:
            # Split the line into columns
            parts = line.strip().split()
            # If there are more than 3 parts, concatenate the extra columns
            if len(parts) > 3:
                parts[2] = "-".join(parts[2:])
                parts = parts[:3]  # Keep only the first three columns
            # Write the fixed line back to the same file
            outfile.write(" ".join(parts) + "\n")

# Define file paths
phylum_file = "Homologs-based-host-prediction-phylum.txt"
family_file = "Homologs-based-host-prediction-family.txt"
genus_file = "Homologs-based-host-prediction-genus.txt"

# Preprocess files in place to fix any irregularities
preprocess_file_inplace(phylum_file)
preprocess_file_inplace(family_file)
preprocess_file_inplace(genus_file)

# Load the fixed files into DataFrames
phylum_df = pd.read_csv(phylum_file, sep=" ", header=None, names=["Virus_ID", "Homologs", "Predicted_phylum"])
family_df = pd.read_csv(family_file, sep=" ", header=None, names=["Virus_ID", "Homologs", "Predicted_family"])
genus_df = pd.read_csv(genus_file, sep=" ", header=None, names=["Virus_ID", "Homologs", "Predicted_genus"])

# Combine Homologs and Predicted columns with square brackets
phylum_df["Predicted_phylum"] = phylum_df["Predicted_phylum"] + " [" + phylum_df["Homologs"].astype(str) + "]"
family_df["Predicted_family"] = family_df["Predicted_family"] + " [" + family_df["Homologs"].astype(str) + "]"
genus_df["Predicted_genus"] = genus_df["Predicted_genus"] + " [" + genus_df["Homologs"].astype(str) + "]"

# Merge the DataFrames on Virus_ID, keeping all rows (outer join)
merged_df = pd.merge(phylum_df[["Virus_ID", "Predicted_phylum"]], 
                     family_df[["Virus_ID", "Predicted_family"]], on="Virus_ID", how="outer", suffixes=('_phylum', '_family'))
merged_df = pd.merge(merged_df, genus_df[["Virus_ID", "Predicted_genus"]], on="Virus_ID", how="outer", suffixes=('', '_genus'))

# Rearrange columns and handle missing values
merged_df = merged_df[["Virus_ID", "Predicted_phylum", "Predicted_family", "Predicted_genus"]]
merged_df.fillna("NA", inplace=True)

# Save the concatenated result to a TSV file
output_file = "HHP_host_prediction.tsv"
merged_df.to_csv(output_file, sep="\t", index=False)

print(f"Concatenated file saved as TSV: {output_file}")
print(f"Files cleaned and saved in place: {phylum_file}, {family_file}, {genus_file}")

