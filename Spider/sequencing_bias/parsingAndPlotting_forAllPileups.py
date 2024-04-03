#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
import glob
import os
from tqdm import tqdm
import psutil
import sys

#### ADD SCAFFOLD NUMBER TO SALL BED FILES, PLOTS, ETC THAT ARE BEING SAVED

class Coverage:
    
    def __init__(self, pileup_glooblist, gff_FirstE_gloobList, gff_LastE_gloobList):

        self.pileup_column_names = ['chrom', 'pos', 'NT', 'depth', 'irr_1', 'irr_2']
        self.gff_column_names = ['chrom', 'irr_1', 'exon', 'start', 'end', 'irr_2', 'irr_3', 'irr_4', 'ID']

        self.pileup = pd.read_csv(pileup_glooblist, sep='\t', header = None, index_col = False, names = self.pileup_column_names)
        self.gff_FirstE = pd.read_csv(gff_FirstE_gloobList, sep='\s+', index_col = False, header = None, names = self.gff_column_names)
        self.gff_LastE = pd.read_csv(gff_LastE_gloobList, sep='\s+', index_col = False, header = None, names = self.gff_column_names)

    def check_memory(self):
        """Check the current memory usage, and raise an exception if it's too high."""
        threshold = 0.8  # If we're using more than 80% of memory, we consider it an error
        usage = psutil.virtual_memory().percent / 100
        if usage > threshold:
            raise MemoryError(f"Memory usage is over {threshold * 100}%")

    # Extract coverage depth at the first and last exon for every gene
    def extract_coverage(self):
        Exon1_df = pd.DataFrame(columns=['chrom', 'start', 'end', 'ID', 'pos', 'depth'])
        ExonLast_df = pd.DataFrame(columns=['chrom', 'start', 'end', 'ID', 'pos', 'depth'])

        #Extract coverage depth at the first and last exon for every gene
        for i, row in tqdm(self.pileup.iterrows(), total=self.pileup.shape[0], desc="Iterating through .mpileup file", position=0, leave=True): #tqdm is a progress bar that appears in the terminal
            
            #Check if position is within the start and end of the first exon
            exon1_mask = (self.gff_FirstE['start'] <= row['pos']) & (self.gff_FirstE['end'] >= row['pos'])
            if exon1_mask.any(): #if there is any True value in the mask (a boolean)
                #Get the exon data, indexing based on the mask
                exon1_data = self.gff_FirstE[exon1_mask]
                #Append the data to the bed dataframe
                Exon1_df = pd.concat([Exon1_df, pd.DataFrame({
                    'chrom': [exon1_data['chrom'].values[0]], #get the value of the first element in the array (the only value)
                    'start': [exon1_data['start'].values[0]],
                    'end': [exon1_data['end'].values[0]],
                    'ID': [exon1_data['ID'].values[0]],
                    'pos': [row['pos']], #getting position from the pileup file
                    'depth': [row['depth']]
                })], ignore_index=True)

            #Check if position is within the start and end of the last exon
            exonLast_mask = (self.gff_LastE['start'] <= row['pos']) & (self.gff_LastE['end'] >= row['pos'])
            if exonLast_mask.any():
                #Get the exon data, indexing based on the mask
                exonLast_data = self.gff_LastE[exonLast_mask]
                #Append the data to the bed dataframe
                ExonLast_df = pd.concat([ExonLast_df, pd.DataFrame({
                    'chrom': [exonLast_data['chrom'].values[0]],
                    'start':[exonLast_data['start'].values[0]],
                    'end': [exonLast_data['end'].values[0]],
                    'ID': [exonLast_data['ID'].values[0]],
                    'pos': [row['pos']],
                    'depth': [row['depth']]
                })], ignore_index=True)
        
        return Exon1_df, ExonLast_df
    
    def concatenate_coverage(self, Exon1_df, ExonLast_df):
        
        self.check_memory()

        #Average coverage depth at the first and last exon for every gene
        Exon1_bed = Exon1_df.groupby('ID').agg({'chrom':'first', 'start':'first', 'end':'first', 'pos':'first', 'depth':'mean'}).reset_index() #group by ID, aggregating the data by taking the first value of each column, but averaging depth
        Exon1_bed.rename(columns={'depth': 'average_depth'}, inplace=True) #renaming the 'depth' column to 'average_depth'
        Exon1_bed = Exon1_bed.drop(columns=['pos']) #removing posiiton (as it is irrelevant after averaging)

        ExonLast_bed = ExonLast_df.groupby('ID').agg({'chrom':'first', 'start':'first', 'end':'first', 'pos':'first', 'depth':'mean'}).reset_index()
        ExonLast_bed.rename(columns={'depth': 'average_depth'}, inplace=True)
        ExonLast_bed = ExonLast_bed.drop(columns=['pos'])

        Exon1_bed['average_depth'] = pd.to_numeric(Exon1_bed['average_depth'], errors='coerce') #coerce will turn any non-numeric values into NaN
        ExonLast_bed['average_depth'] = pd.to_numeric(ExonLast_bed['average_depth'], errors='coerce')

        Exon1_bed.to_csv(f'/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/spider/sequencing_bias/output_data/{"_".join(self.pileup.split("_")[0:3])}_Exon1.bed', sep='\t', index=False)
        ExonLast_bed.to_csv(f'/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/spider/sequencing_bias/output_data/{"_".join(self.pileup.split("_")[0:3])}_ExonLast.bed', sep='\t', index=False)

        return Exon1_bed, ExonLast_bed  

# Plotting the coverage depth at the first and last exon for every pileup file
    def plotting_coverage(self, Exon1_bed, ExonLast_bed):
        
        plt.hist(Exon1_bed['average_depth'], bins=50, alpha=0.5, label='First Exon', color='pink')
        plt.hist(ExonLast_bed['average_depth'], bins=50, alpha=0.5, label='Last Exon', color='cyan')
        plt.axvline(np.mean(Exon1_bed['average_depth']), color='red', linestyle='dashed', linewidth=1, label = 'Mean Coverage')
        plt.axvline(np.mean(ExonLast_bed['average_depth']), color='blue', linestyle='dashed', linewidth=1, label = 'Mean Coverage')

        #Perform Mann-Whitney U test (non-parametic))
        stat, p = mannwhitneyu(Exon1_bed['average_depth'], ExonLast_bed['average_depth'])
        plt.text(0.1, 0.9, f'Mann-Whitney p-value: {p:.3e}', transform=plt.gca().transAxes) #the e specifies scientific notation

        plt.xlabel('Average Coverage Depth for Each Gene')
        plt.ylabel('Frequency')
        plt.yscale('log')
        plt.title(f'Coverage Depth at First and Last Exon for for Every Gene in {"_".join(self.pileup.split("_")[0:3])}')
        plt.legend()

        # Save the plot
        plt.savefig(f'/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/Spider/sequencing_bias/plots/coverageDepthHistogram_{"_".join(self.pileup.split("_")[0:3])}.png', dpi = 900)

        plt.show()

# pileup_gloobList = glob.glob('/media/will/mimic/k2024/pileups/*.mpileup')
# gff_FirstE_gloobList = glob.glob('/media/will/mimic/k2024/exons/*.gff3')
# gff_LastE_gloobList = glob.glob('/media/will/mimic/k2024/exons/*.last.exon.gff3')

pileup_gloobList = glob.glob('/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/Spider/sequencing_bias/pileups/*.mpileup')
gff_FirstE_gloobList = glob.glob('/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/Spider/sequencing_bias/exons_gff_files/*.gff3')
gff_LastE_gloobList = glob.glob('/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/Spider/sequencing_bias/exons_gff_files/*.last.exon.gff3')


if pileup_gloobList and gff_FirstE_gloobList and gff_LastE_gloobList:
        for pileup_file in tqdm(pileup_gloobList, desc="Processing pileup files"):
        base_name = os.path.splitext(os.path.basename(pileup_file))[0]
        scaffold_number_pileup = "_".join(base_name.split("_")[1:3]) # ex: WholeFemale_scaffold_5.mpileup

        for gff_FirstE_file in gff_FirstE_gloobList:
            base_name = os.path.splitext(os.path.basename(gff_FirstE_file))[0]
            scaffold_number_FirstE = base_name.split(".")[0]  # ex: scaffold_10.last.exon.gff3

            for gff_LastE_file in gff_LastE_gloobList:
                base_name = os.path.splitext(os.path.basename(gff_LastE_file))[0]
                scaffold_number_LastE = base_name.split(".")[0]  

                if scaffold_number_pileup == scaffold_number_FirstE == scaffold_number_LastE:
                    print(f"Matching files found for {os.path.basename(pileup_file)}, analysis is running!")
                    run = Coverage(pileup_file, gff_FirstE_file, gff_LastE_file)
                    Exon1_df, ExonLast_df = run.extract_coverage()
                    Exon1_bed, ExonLast_bed = run.concatenate_coverage(Exon1_df, ExonLast_df)
                    run.plotting_coverage(Exon1_bed, ExonLast_bed)
else:
    print("No matching files found")