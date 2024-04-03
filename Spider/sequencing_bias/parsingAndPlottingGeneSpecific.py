#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
from tqdm import tqdm
import multiprocessing as mp


#Create empty list  to store the dataframes of coverage depth at the first and last exon for every gene
# Exon1_dfs = []
# ExonLast_dfs = []

# #Extract coverage depth at the first and last exon for every gene
# for i, row in tqdm(pileup.iterrows(), total=pileup.shape[0]): #tqdm is a progress bar that appears in the terminal
    
#     #Check if position is within the start and end of the first exon
#     exon1_mask = (first_gff['start'] <= row['pos']) & (first_gff['end'] >= row['pos'])
#     if exon1_mask.any(): #if there is any True value in the mask (a boolean)
#         #Get the exon data, indexing based on the mask
#         exon1_data = first_gff[exon1_mask]
#         #Append the data to the bed dataframe
#         Exon1_dfs.append(pd.DataFrame({
#             'chrom': [exon1_data['chrom'].values[0]], #get the value of the first element in the array (the only value)
#             'start': [exon1_data['start'].values[0]],
#             'end': [exon1_data['end'].values[0]],
#             'ID': [exon1_data['ID'].values[0]],
#             'pos': [row['pos']], #getting position from the pileup file
#             'depth': [np.float32(row['depth'])]
#         }), ignore_index=True)

#     #Check if position is within the start and end of the last exon
#     exonLast_mask = (last_gff['start'] <= row['pos']) & (last_gff['end'] >= row['pos'])
#     if exonLast_mask.any():
#         #Get the exon data, indexing based on the mask
#         exonLast_data = last_gff[exonLast_mask]
#         #Append the data to the bed dataframe
#         ExonLast_dfs.append(pd.DataFrame({
#             'chrom': [exonLast_data['chrom'].values[0]],
#             'start':[exonLast_data['start'].values[0]],
#             'end': [exonLast_data['end'].values[0]],
#             'ID': [exonLast_data['ID'].values[0]],
#             'pos': [row['pos']],
#             'depth': [np.float32(row['depth'])]
#         }), ignore_index=True)


# def extract_coverege(row):
#     Exon1_df = pd.DataFrame()
#     ExonLast_df = pd.DataFrame()
    
#     #Check if position is within the start and end of the first exon
#     exon1_mask = (first_gff['start'] <= row['pos']) & (first_gff['end'] >= row['pos'])
#     if exon1_mask.any(): 
#         exon1_data = first_gff[exon1_mask]
#         Exon1_df = pd.DataFrame({
#             'chrom': [exon1_data['chrom'].values[0]], 
#             'start': [exon1_data['start'].values[0]],
#             'end': [exon1_data['end'].values[0]],
#             'ID': [exon1_data['ID'].values[0]],
#             'pos': [row['pos']], 
#             'depth': [np.float32(row['depth'])]
#         })

#     #Check if position is within the start and end of the last exon
#     exonLast_mask = (last_gff['start'] <= row['pos']) & (last_gff['end'] >= row['pos'])
#     if exonLast_mask.any():
#         exonLast_data = last_gff[exonLast_mask]
#         ExonLast_df = pd.DataFrame({
#             'chrom': [exonLast_data['chrom'].values[0]],
#             'start':[exonLast_data['start'].values[0]],
#             'end': [exonLast_data['end'].values[0]],
#             'ID': [exonLast_data['ID'].values[0]],
#             'pos': [row['pos']],
#             'depth': [np.float32(row['depth'])]
#         })

#     return Exon1_df, ExonLast_df

# results = list([extract_coverege(row) for _, row in pileup.head(100).iterrows()])
# # print(results)

# Exon1_dfs = pd.concat([item[0] for item in results])
# ExonLast_dfs = pd.concat([item[1] for item in results])
# print(Exon1_dfs)

# #Using a multiprocessing Pool to process the DataFrame in parallel
# SIZE = pileup.shape[0] // (10*mp.cpu_count())  #not sure what to set chunksize to....
# results = []
# if __name__ == '__main__':
#     with mp.Pool(mp.cpu_count()) as pool: #using all CPUs, in my case 8
#        # results = list(tqdm(pool.imap(extract_coverege, pileup.head(100).iterrows(), chunksize=1), total=pileup.head(100).shape[0])) #total is the number of rows in the pileup file
#         results = results.append((tqdm(pool.imap(extract_coverege, pileup.head(100).iterrows(), chunksize=1), total=pileup.head(100).shape[0]))) #total is the number of rows in the pileup file



# # Concatenate the results
# Exon1_dfs = pd.concat([item[0] for item in results if not item[0].empty])
# ExonLast_dfs = pd.concat([item[1] for item in results if not item[1].empty])





#Load in files
##columns I am naming "irr" are columns that I am not interested in (irrelevant)
pileup_column_names = ['chrom', 'pos', 'NT', 'depth', 'irr_1', 'irr_2']
pileup = pd.read_csv("/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/Spider/sequencing_bias/pileups/WholeFemale.scaffold_10.mpileup", sep='\t', header = None, index_col = False, names = pileup_column_names)

gff_column_names = ['chrom', 'irr_1', 'exon', 'start', 'end', 'irr_2', 'irr_3', 'irr_4', 'ID']
first_gff = pd.read_csv("/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/Spider/sequencing_bias/exons_gff_files/scaffold_10.exon1.gff3", sep='\s+', index_col = False, header = None, names = gff_column_names)
last_gff = pd.read_csv("/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/Spider/sequencing_bias/exons_gff_files/scaffold_10.last.exon.gff3", sep='\s+', index_col = False, header = None, names = gff_column_names)


# #Create empty list  to store the dataframes of coverage depth at the first and last exon for every gene
# Exon1_dfs = []
# ExonLast_dfs = []

# #Extract coverage depth at the first and last exon for every gene
# for i, row in tqdm(pileup.iterrows(), total=pileup.shape[0]): #tqdm is a progress bar that appears in the terminal
    
#     #Check if position is within the start and end of the first exon
#     exon1_mask = (first_gff['start'] <= row['pos']) & (first_gff['end'] >= row['pos'])
#     if exon1_mask.any(): #if there is any True value in the mask (a boolean)
#         #Get the exon data, indexing based on the mask
#         exon1_data = first_gff[exon1_mask]
#         #Append the data to the bed dataframe
#         Exon1_dfs.append(pd.DataFrame({
#             'chrom': [exon1_data['chrom'].values[0]], #get the value of the first element in the array (the only value)
#             'start': [exon1_data['start'].values[0]],
#             'end': [exon1_data['end'].values[0]],
#             'ID': [exon1_data['ID'].values[0]],
#             'pos': [row['pos']], #getting position from the pileup file
#             'depth': [np.float32(row['depth'])]
#         }))

#     #Check if position is within the start and end of the last exon
#     exonLast_mask = (last_gff['start'] <= row['pos']) & (last_gff['end'] >= row['pos'])
#     if exonLast_mask.any():
#         #Get the exon data, indexing based on the mask
#         exonLast_data = last_gff[exonLast_mask]
#         #Append the data to the bed dataframe
#         ExonLast_dfs.append(pd.DataFrame({
#             'chrom': [exonLast_data['chrom'].values[0]],
#             'start':[exonLast_data['start'].values[0]],
#             'end': [exonLast_data['end'].values[0]],
#             'ID': [exonLast_data['ID'].values[0]],
#             'pos': [row['pos']],
#             'depth': [np.float32(row['depth'])]
#         }))

# #concatenate the results
# Exon1_concatenated = pd.concat(Exon1_dfs, ignore_index = True)
# ExonLast_concatenated = pd.concat(ExonLast_dfs, ignore_index = True)

def process_row(row):
    exon1_mask = (first_gff['start'] <= row['pos']) & (first_gff['end'] >= row['pos'])
    if exon1_mask.any():
        exon1_data = first_gff[exon1_mask]
        return pd.DataFrame({
            'chrom': [exon1_data['chrom'].values[0]],
            'start': [exon1_data['start'].values[0]],
            'end': [exon1_data['end'].values[0]],
            'ID': [exon1_data['ID'].values[0]],
            'pos': [row['pos']],
            'depth': [np.float32(row['depth'])]
        })
    else:
        return pd.DataFrame()

Exon1_dfs = pileup.apply(process_row, axis=1).tolist()

def process_row(row):
    exonLast_mask = (last_gff['start'] <= row['pos']) & (last_gff['end'] >= row['pos'])
    if exonLast_mask.any():
        exonLast_data = last_gff[exonLast_mask]
        return pd.DataFrame({
            'chrom': [exonLast_data['chrom'].values[0]],
            'start':[exonLast_data['start'].values[0]],
            'end': [exonLast_data['end'].values[0]],
            'ID': [exonLast_data['ID'].values[0]],
            'pos': [row['pos']],
            'depth': [np.float32(row['depth'])]
        })
    else:
        return pd.DataFrame()


tqdm.pandas(desc="Processing rows")
ExonLast_dfs = pileup.progress_apply(process_row, axis=1)
ExonLast_dfs = [item for sublist in ExonLast_dfs for item in sublist]


#Average coverage depth at the first and last exon for every gene
        ##should I be taking the median coverage depth instead of the mean???? 
Exon1_bed = Exon1_dfs.groupby('ID').agg({'chrom':'first', 'start':'first', 'end':'first', 'pos':'first', 'depth':'mean'}).reset_index() #group by ID, aggregating the data by taking the first value of each column, but averaging depth
Exon1_bed.rename(columns={'depth': 'average_depth'}, inplace=True) #renaming the 'depth' column to 'average_depth'
Exon1_bed = Exon1_bed.drop(columns=['pos']) #removing posiiton (as it is irrelevant after averaging)

ExonLast_bed = ExonLast_dfs.groupby('ID').agg({'chrom':'first', 'start':'first', 'end':'first', 'pos':'first', 'depth':'mean'}).reset_index()
ExonLast_bed.rename(columns={'depth': 'average_depth'}, inplace=True)
ExonLast_bed = ExonLast_bed.drop(columns=['pos'])

Exon1_bed['average_depth'] = pd.to_numeric(Exon1_bed['average_depth'], errors='coerce', downcast = 'float') #coerce will turn any non-numeric values into NaN
ExonLast_bed['average_depth'] = pd.to_numeric(ExonLast_bed['average_depth'], errors='coerce', downcast = 'float')

#Plot the coverage depth at the first and last exon for every gene
plt.hist(Exon1_bed['average_depth'], bins=50, alpha=0.5, label='First Exon', color='pink')
plt.hist(ExonLast_bed['average_depth'], bins=50, alpha=0.5, label='Last Exon', color='cyan')
mean_exon1coverage = np.mean(Exon1_bed['average_depth'])
plt.axvline(mean_exon1coverage, color='red', linestyle='dashed', linewidth=1, label = f'Mean Coverage ({mean_exon1coverage}))')
mean_lastExoncoverage = np.mean(ExonLast_bed['average_depth'])
plt.axvline(mean_lastExoncoverage, color='blue', linestyle='dashed', linewidth=1, label = f'Mean Coverage ({mean_lastExoncoverage})')

#Perform Mann-Whitney U test (non-parametic))
stat, p = mannwhitneyu(Exon1_bed['average_depth'], ExonLast_bed['average_depth'])
plt.text(0.1, 0.9, f'Mann-Whitney p-value: {p:.3e}', transform=plt.gca().transAxes) #the e specifies scientific notation

plt.xlabel('Average Coverage Depth for Each Gene')
plt.ylabel('Frequency')
plt.yscale('log')
plt.title('Coverage Depth at First and Last Exon for Every Gene')
plt.legend()

#Save the plot
#plt.savefig('/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/spider/sequencing_bias/plots/coverageDepthHistogramGENElevel.png', dpi = 900)
plt.savefig('/Users/cmdb/Desktop/GORDUS_rotation/Gordus_scripts/spider/sequencing_bias/plots/coverageDepthHistogramGENElevelTEST.png', dpi = 900)


plt.show()
