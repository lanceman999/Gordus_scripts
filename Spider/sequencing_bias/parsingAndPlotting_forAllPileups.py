#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
import glob
import os

class Coverage:
    
    def __init__(self, pileup_glooblist, gff_FirstE_gloobList, gff_LastE_gloobList):

        self.pileup_column_names = ['chrom', 'pos', 'NT', 'depth', 'irr_1', 'irr_2']
        self.gff_column_names = ['chrom', 'irr_1', 'exon', 'start', 'end', 'irr_2', 'irr_3', 'irr_4', 'ID']

        self.pileup = pd.read_csv(pileup_glooblist, sep='\t', header = None, index_col = False, names = self.pileup_column_names)
        self.gff_FirstE = pd.read_csv(gff_FirstE_gloobList, sep='\s+', index_col = False, header = None, names = self.gff_column_names)
        self.gff_LastE = pd.read_csv(gff_LastE_gloobList, sep='\s+', index_col = False, header = None, names = self.gff_column_names)

    # Extract coverage depth at the first and last exon for every gene
    def extract_coverage(self):
        Exon1_depth = []
        ExonLast_depth = []

        for i, row in self.pileup.iterrows():
            if ((self.gff_LastE['start'] <= row['pos']) & (self.gff_FirstE['end'] >= row['pos'])).any(): 
                Exon1_depth.append(row['depth'])
            if ((self.gff_LastE['start'] <= row['pos']) & (self.gff_FirstE['end'] >= row['pos'])).any():
                ExonLast_depth.append(row['depth'])
        
        return Exon1_depth, ExonLast_depth

# Plotting the coverage depth at the first and last exon for every pileup file
    def plotting_coverage(self, Exon1_depth, ExonLast_depth):
        plt.hist(Exon1_depth, bins=50, alpha=0.5, label='First Exon', color='violet')
        plt.hist(ExonLast_depth, bins=50, alpha=0.5, label='Last Exon', color='green')
        plt.axvline(np.mean(Exon1_depth), color='red', linestyle='dashed', linewidth=1, label = 'Mean Coverage')
        plt.axvline(np.mean(ExonLast_depth), color='red', linestyle='dashed', linewidth=1)

        # Perform Mann-Whitney U test
        stat, p = mannwhitneyu(Exon1_depth, ExonLast_depth)
        plt.text(0.1, 0.9, f'Mann-Whitney p-value: {p:.3f}', transform=plt.gca().transAxes)

        plt.xlabel('Coverage Depth')
        plt.ylabel('Frequency')
        plt.title(f'Coverage Depth at First and Last Exon for {self.pileup}')
        plt.legend()

        # Save the plot
        plt.savefig(f'coverageDepthHistogram_{self.pileup}.png', dpi = 900)

        plt.show()

pileup_gloobList = glob.glob('*.mpileup')
gff_FirstE_gloobList = glob.glob('*.gff3')
gff_LastE_gloobList = glob.glob("*.last.exon.gff3")

for pileup_file in pileup_gloobList:
    base_name = os.path.splitext(os.path.basename(pileup_file))[0]
    scaffold_number_pileup = "_".join(base_name.split("_")[1:3])

    for gff_FirstE_file in gff_FirstE_gloobList:
        base_name = os.path.splitext(os.path.basename(gff_FirstE_file))[0]
        scaffold_number_FirstE = base_name.split(".")[0]  

        for gff_LastE_file in gff_LastE_gloobList:
            base_name = os.path.splitext(os.path.basename(gff_LastE_file))[0]
            scaffold_number_LastE = base_name.split(".")[0]  

            if scaffold_number_pileup == scaffold_number_FirstE == scaffold_number_LastE:
                run = Coverage(pileup_file, gff_FirstE_file, gff_LastE_file)
