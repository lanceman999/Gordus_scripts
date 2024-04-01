#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
from tqdm import tqdm

# Load in files
## columns I am naming "irr" are columns that I am not interested in (irrelevant)
pileup_column_names = ['chrom', 'pos', 'NT', 'depth', 'irr_1', 'irr_2']
pileup = pd.read_csv("/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/Spider/sequencing_bias/pileups/WholeFemale.scaffold_10.mpileup", sep='\t', header = None, index_col = False, names = pileup_column_names)

gff_column_names = ['chrom', 'irr_1', 'exon', 'start', 'end', 'irr_2', 'irr_3', 'irr_4', 'ID']
first_gff = pd.read_csv("/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/Spider/sequencing_bias/exons_gff_files/scaffold_10.exon1.gff3", sep='\s+', index_col = False, header = None, names = gff_column_names)
last_gff = pd.read_csv("/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/Spider/sequencing_bias/exons_gff_files/scaffold_10.last.exon.gff3", sep='\s+', index_col = False, header = None, names = gff_column_names)

#Sanity check!
#print(pileup.head())
#print(first_gff.head())
#print(last_gff.head())  

# Extract coverage depth at the first and last exon for every gene
Exon1_bed = pd.DataFrame(columns=['chrom', 'start', 'end', 'ID', 'pos', 'depth'])
ExonLast_bed = pd.DataFrame(columns=['chrom', 'start', 'end', 'ID', 'pos', 'depth'])

for i, row in tqdm(pileup.iterrows(), total=pileup.shape[0]):
    # Check if position is within the start and end of the first exon
    exon1_mask = (first_gff['start'] <= row['pos']) & (first_gff['end'] >= row['pos'])
    if exon1_mask.any():
        # Get the exon data
        exon1_data = first_gff[exon1_mask]
        # Append the data to the bed dataframe
        Exon1_bed = Exon1_bed.append({
            'chrom': exon1_data['chrom'].values[0],
            'start': exon1_data['start'].values[0],
            'end': exon1_data['end'].values[0],
            'ID': exon1_data['ID'].values[0],
            'pos': row['pos'],
            'depth': row['depth']
        }, ignore_index=True)

    # Check if position is within the start and end of the last exon
    exonLast_mask = (last_gff['start'] <= row['pos']) & (last_gff['end'] >= row['pos'])
    if exonLast_mask.any():
        # Get the exon data
        exonLast_data = last_gff[exonLast_mask]
        # Append the data to the bed dataframe
        ExonLast_bed = ExonLast_bed.append({
            'chrom': exonLast_data['chrom'].values[0],
            'start': exonLast_data['start'].values[0],
            'end': exonLast_data['end'].values[0],
            'ID': exonLast_data['ID'].values[0],
            'pos': row['pos'],
            'depth': row['depth']
        }, ignore_index=True)

# Average coverage depth at the first and last exon for every gene
Exon1_bed = Exon1_bed.groupby('ID').agg({'chrom':'first', 'start':'first', 'end':'first', 'pos':'first', 'depth':'mean'}).reset_index()
Exon1_bed.rename(columns={'depth': 'average_depth'}, inplace=True)
Exon1_bed = Exon1_bed.drop(columns=['pos'])

ExonLast_bed = ExonLast_bed.groupby('ID').agg({'chrom':'first', 'start':'first', 'end':'first', 'pos':'first', 'depth':'mean'}).reset_index()
ExonLast_bed.rename(columns={'depth': 'average_depth'}, inplace=True)
ExonLast_bed = ExonLast_bed.drop(columns=['pos'])

Exon1_bed.to_csv('Exon1.bed', sep='\t', index=False)
ExonLast_bed.to_csv('ExonLast.bed', sep='\t', index=False)

# Plot the coverage depth at the first and last exon for every gene
plt.hist(Exon1_bed['average_depth'], bins=50, alpha=0.5, label='First Exon', color='violet')
plt.hist(ExonLast_bed['average_depth'], bins=50, alpha=0.5, label='Last Exon', color='green')
plt.axvline(np.mean(Exon1_bed['average_depth']), color='red', linestyle='dashed', linewidth=1)
plt.axvline(np.mean(ExonLast_bed['average_depth']), color='red', linestyle='dashed', linewidth=1)

# Perform Mann-Whitney U test
stat, p = mannwhitneyu(Exon1_bed['average_depth'], ExonLast_bed['average_depth'])
plt.text(0.1, 0.9, f'Mann-Whitney p-value: {p:.3f}', transform=plt.gca().transAxes)

plt.xlabel('Average Coverage Depth for Each Gene')
plt.ylabel('Frequency')
plt.title('Coverage Depth at First and Last Exon for every gene')

plt.legend()
plt.show()

# Save the plot
#plt.savefig('coverageDepthHistogram.png', dpi = 900)