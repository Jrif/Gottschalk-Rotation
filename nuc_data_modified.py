# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 10:35:13 2020

@author: rober
"""

import pandas as pd


#read in STAT1 nucleotide sequences

nuc_data_df = pd.read_csv('STAT1_nuc_data.txt', delimiter='\t', header=None, sep=' ')



#formating nuc_data_df

#dropping NAN
nuc_data_df.dropna(inplace = True)

#column spliting all spaces

nuc_data_df = nuc_data_df.iloc[:,0].str.split(' ', n=-1, expand=True)

#saving df to modified csv

nuc_data_df.to_csv('nuc_data_modified1.csv')


nuc_data_df = nuc_data_df.iloc[:,0].str.split('=|:|-', n=-1, expand=True) #column splitting at =, :, and -

#merging dataframes to complete data frame

nuc_data_df_modified = pd.read_csv('nuc_data_modified1.csv')
nuc_data_df = nuc_data_df.drop(columns=[0])
nuc_data_df = pd.merge(nuc_data_df, nuc_data_df_modified, right_index=True, left_index=True)
nuc_data_df = nuc_data_df.drop(columns=['Unnamed: 0', '0'])
nuc_data_df.columns = ['Chromosome Number', 'Starting Base Loc', 'Ending Base Loc', 'Strand', 'Sequence']
nuc_data_df['Strand'] = nuc_data_df['Strand'].str.replace('strand=', '') #change the strand from 'strand=+ or -' to just + or -


nuc_data_df.to_csv('nuc_data_modified.csv') #saving modified nuc_data_df

print (nuc_data_df.head(5))

    


