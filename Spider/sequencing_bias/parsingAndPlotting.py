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

# Extract coverage depth at the first and last exon for every gene
Exon1_depth = []
ExonLast_depth = []

for i, row in tqdm(pileup.iterrows(), total=pileup.shape[0]):
    if ((first_gff['start'] <= row['pos']) & (first_gff['end'] >= row['pos'])).any(): 
        Exon1_depth.append(row['depth'])
    if ((last_gff['start'] <= row['pos']) & (last_gff['end'] >= row['pos'])).any():
        ExonLast_depth.append(row['depth'])

# Plot the coverage depth at the first and last exon for every gene
plt.hist(Exon1_depth, bins=50, alpha=0.5, label='First Exon', color='pink')
plt.hist(ExonLast_depth, bins=50, alpha=0.5, label='Last Exon', color='cyan')
plt.axvline(np.mean(Exon1_depth), color='red', linestyle='dashed', linewidth=1, label = 'Mean Coverage')
plt.axvline(np.mean(ExonLast_depth), color='blue', linestyle='dashed', linewidth=1)

# Perform Mann-Whitney U test (non-parametric)
stat, p = mannwhitneyu(Exon1_depth, ExonLast_depth)
plt.text(0.1, 0.9, f'Mann-Whitney p-value: {p:.3e}', transform=plt.gca().transAxes) # the e specifies scientific notation

plt.xlabel('Coverage Depth for a NT')
plt.ylabel('Frequency')
plt.yscale("log")
plt.title('Coverage Depth at First and Last Exon for Every NT')
plt.legend()

# Save the plot
plt.savefig('/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/spider/sequencing_bias/plots/coverageDepthHistogram.png', dpi = 900)
plt.show()
