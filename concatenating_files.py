#!/usr/bin/env python3
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Concatenating the files into HPP_host_prediction.txt
# Date: 20241221
import pandas as pd

# Define file paths
phylum_file = "Homologs-based-host-prediction-phylum.txt"
family_file = "Homologs-based-host-prediction-family.txt"
genus_file = "Homologs-based-host-prediction-genus.txt"

# Load the files into DataFrames
phylum_df = pd.read_csv(phylum_file, sep=" ", header=None, names=["Virus_ID", "Homologs", "Predicted_phylum"])
family_df = pd.read_csv(family_file, sep=" ", header=None, names=["Virus_ID", "Homologs", "Predicted_family"])
genus_df = pd.read_csv(genus_file, sep=" ", header=None, names=["Virus_ID", "Homologs", "Predicted_genus"])

# Save each updated DataFrame back to the original files with the header added
phylum_df.to_csv(phylum_file, sep=" ", index=False, header=True)
family_df.to_csv(family_file, sep=" ", index=False, header=True)
genus_df.to_csv(genus_file, sep=" ", index=False, header=True)

# Merge the DataFrames on Virus_ID, keeping all rows (outer join)
merged_df = pd.merge(phylum_df, family_df, on="Virus_ID", how="outer", suffixes=('_phylum', '_family'))
merged_df = pd.merge(merged_df, genus_df, on="Virus_ID", how="outer", suffixes=('', '_genus'))

# Rearrange columns and handle missing values
merged_df = merged_df[["Virus_ID", "Predicted_phylum", "Predicted_family", "Predicted_genus"]]
merged_df.fillna("NA", inplace=True)

# Save the concatenated result to a new file
output_file = "HPP_host_prediction.txt"
merged_df.to_csv(output_file, sep="\t", index=False)

print(f"Concatenated file saved to {output_file}")
print(f"Headers added to {phylum_file}, {family_file}, and {genus_file}")
