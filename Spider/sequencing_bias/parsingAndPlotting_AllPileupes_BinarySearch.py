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
        self.pileup.set_index('pos', inplace=True)
        self.pileup.sort_index(inplace=True)

        self.gff_FirstE = pd.read_csv(gff_FirstE_gloobList, sep='\s+', index_col = False, header = None, names = self.gff_column_names)
        self.gff_LastE = pd.read_csv(gff_LastE_gloobList, sep='\s+', index_col = False, header = None, names = self.gff_column_names)

        self.pileup_name_0 = os.path.basename(pileup_glooblist)
        self.pileup_name_temp = os.path.splitext(self.pileup_name_0)[0] # example looks like: WholeFemale_scaffold_5_SORTED
        self.pileup_name = "_".join(self.pileup_name_temp.split("_")[0:3]) # example looks like: WholeFemale_scaffold_5
        #print(self.pileup_name)

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

        rows_list_E1 = []
        rows_list_EL = []
        for i, row in tqdm(self.gff_FirstE.iterrows(), total=self.gff_FirstE.shape[0], desc="Iterating through first Exon gff file", position=0, leave=True):
            exon_data = self.pileup.loc[row['start']:row['end']]
            for pos, data in exon_data.iterrows():
                dict_E1 = {
                    'chrom': row['chrom'],
                    'start': row['start'],
                    'end': row['end'],
                    'ID': row['ID'],
                    'pos': pos,
                    'depth': data['depth']
                }
                rows_list_E1.append(dict_E1)

        for i, row in tqdm(self.gff_LastE.iterrows(), total=self.gff_LastE.shape[0], desc="Iterating through last Exon gff file", position=0, leave=True):
            exon_data = self.pileup.loc[row['start']:row['end']]
            for pos, data in exon_data.iterrows():
                dict_EL = {
                    'chrom': row['chrom'],
                    'start': row['start'],
                    'end': row['end'],
                    'ID': row['ID'],
                    'pos': pos,
                    'depth': data['depth']
                }
                rows_list_EL.append(dict_EL)
        
        
        Exon1_df = pd.DataFrame(rows_list_E1)
        ExonLast_df = pd.DataFrame(rows_list_EL)
           
        return Exon1_df, ExonLast_df
    

    def concatenate_coverage(self, Exon1_df, ExonLast_df):

        self.check_memory()

        #Average coverage depth at the first and last exon for every gene
        Exon1_bed = Exon1_df.groupby('ID').agg({'chrom':'first', 'start':'first', 'end':'first', 'pos':'first', 'depth':'median'}).reset_index() #group by ID, aggregating the data by taking the first value of each column, but averaging depth
        Exon1_bed.rename(columns={'depth': 'median_depth'}, inplace=True) #renaming the 'depth' column to 'median_depth'
        Exon1_bed = Exon1_bed.drop(columns=['pos']) #removing posiiton (as it is irrelevant after averaging)

        ExonLast_bed = ExonLast_df.groupby('ID').agg({'chrom':'first', 'start':'first', 'end':'first', 'pos':'first', 'depth':'mean'}).reset_index()
        ExonLast_bed.rename(columns={'depth': 'median_depth'}, inplace=True)
        ExonLast_bed = ExonLast_bed.drop(columns=['pos'])

        Exon1_bed['median_depth'] = pd.to_numeric(Exon1_bed['median_depth'], errors='coerce') #coerce will turn any non-numeric values into NaN
        ExonLast_bed['median_depth'] = pd.to_numeric(ExonLast_bed['median_depth'], errors='coerce')

        # Exon1_bed.to_csv(f'/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/spider/sequencing_bias/output_data/{self.pileup_name}_Exon1_BINARYTEST.bed', sep='\t', index=False)
        # ExonLast_bed.to_csv(f'/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/spider/sequencing_bias/output_data/{self.pileup_name}_ExonLast_BINARYTEST.bed', sep='\t', index=False)

        Exon1_bed.to_csv(f'/media/will/mimic/k2024/lance/output_data/{self.pileup_name}_Exon1.bed', sep='\t', index=False)
        ExonLast_bed.to_csv(f'/media/will/mimic/k2024/lance/output_data/{self.pileup_name}_ExonLast.bed', sep='\t', index=False)


        return Exon1_bed, ExonLast_bed  

# Plotting the coverage depth at the first and last exon for every pileup file
    def plotting_coverage(self, Exon1_bed, ExonLast_bed):
        
        plt.style.use("ggplot") #to have plot look more like ggplot2 in R
        plt.hist(Exon1_bed['median_depth'], bins=50, alpha=0.5, label='First Exon (FE)', color='red')
        plt.hist(ExonLast_bed['median_depth'], bins=50, alpha=0.5, label='Last Exon (LE)', color='blue')
        plt.axvline(np.mean(Exon1_bed['median_depth']), color='orange', linestyle='dashed', linewidth=1, label = f'FE Mean Coverage: {Exon1_bed["median_depth"].mean():.2f}')
        plt.axvline(np.mean(ExonLast_bed['median_depth']), color='cyan', linestyle='dashed', linewidth=1, label = f'LE Mean Coverage: {ExonLast_bed["median_depth"].mean():.2f}')

        #Perform Mann-Whitney U test (non-parametic))
        stat, p = mannwhitneyu(Exon1_bed['median_depth'], ExonLast_bed['median_depth'])
        plt.text(0.1, 0.9, f'Mann-Whitney p-value: {p:.3e}', transform=plt.gca().transAxes) #the e specifies scientific notation

        plt.xlabel('Median Coverage Depth for Each Gene')
        plt.ylabel('Frequency')
        plt.yscale('log')
        plt.title(f'{self.pileup_name}')
        plt.legend(facecolor='white', loc = 'center right')

        # Save the plot
        #plt.savefig(f'/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/Spider/sequencing_bias/plots/coverageDepthHistogram_{self.pileup_name}_BINARYTEST.png', dpi = 900)
        plt.savefig(f'/media/will/mimic/k2024/lance/plots/coverageDepthHistogram_{self.pileup_name}.png', dpi = 900)

        plt.close()

        #plt.show()

pileup_gloobList = glob.glob('/media/will/mimic/k2024/pileups/*.mpileup')
gff_FirstE_gloobList = glob.glob('/media/will/mimic/k2024/exons/*.exon1.gff3')
gff_LastE_gloobList = glob.glob('/media/will/mimic/k2024/exons/*.last.exon.gff3')

# pileup_gloobList = glob.glob('/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/Spider/sequencing_bias/pileups/SORTED_mpileup/*.mpileup')
# gff_FirstE_gloobList = glob.glob('/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/Spider/sequencing_bias/exons_gff_files/*exon1.gff3')
# gff_LastE_gloobList = glob.glob('/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/Spider/sequencing_bias/exons_gff_files/*.last.exon.gff3')


if pileup_gloobList and gff_FirstE_gloobList and gff_LastE_gloobList:
        for pileup_file in tqdm(pileup_gloobList, desc="Processing pileup files"):
            base_name = os.path.splitext(os.path.basename(pileup_file))[0]
            scaffold_number_pileup = "_".join(base_name.split("_")[1:3]) # ex: WholeFemale_scaffold_5.mpileup would be the original

            for gff_FirstE_file in gff_FirstE_gloobList:
                base_name = os.path.splitext(os.path.basename(gff_FirstE_file))[0]
                scaffold_number_FirstE = base_name.split(".")[0]  # ex: scaffold_10.last.exon.gff3 would be the original

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
    print("No matching files found, check that all files have be processed")