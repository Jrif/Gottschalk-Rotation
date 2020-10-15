# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 10:38:29 2020

@author: rober
"""

import numpy as np
import pandas as pd
import timeit

start = timeit.default_timer()

nuc_data_df = pd.read_csv('nuc_data_modified.csv') # read in nuc_data
nuc_data = np.asarray(nuc_data_df) #convert nuc_data to array

gene_loc_df = pd.read_csv('mm10_annotation.csv') #read in gene location data

# read in gene clusters - gene names only
clus1_genes_df = pd.read_csv('il6il10_clus_1.txt', header=None) 
clus2_genes_df = pd.read_csv('il6il10_clus_2.txt', header=None)
clus3_genes_df = pd.read_csv('il6il10_clus_3.txt', header=None)
clus4_genes_df = pd.read_csv('il6il10_clus_4.txt', header=None)
clus5_genes_df = pd.read_csv('il6il10_clus_5.txt', header=None)
clus6_genes_df = pd.read_csv('il6il10_clus_6.txt', header=None)
clus7_genes_df = pd.read_csv('il6il10_clus_7.txt', header=None)
clus8_genes_df = pd.read_csv('il6il10_clus_8.txt', header=None)
clus9_genes_df = pd.read_csv('il6il10_clus_9.txt', header=None)

# write a function to find the STAT1 motifs near genes in each cluster
# this function can be used for every cluster
# inputs: cluster of interest, nuc_data_df, gene_loc_df
# outputs: array of information noted in doc. 
#(gene name, location of gene, location of STAT1 motif and nucleotide sequence)

# find the location of genes in cluster from gene_data
    # check to see if it is near any stat1 motifs by looking at nuc_data
    # if it is near a stat1 motif, save the appropriate information to an array
    # the above steps will need to be repeated for all the gene in the cluster

def find_motifs(cluster,nuc_data):
    
    
    clusteredgenes = cluster[0].values.tolist() #turn clusters into a list
    clusteredgenes_df = pd.DataFrame() #create empty dataframe 
    clusteredgenes_df = gene_loc_df[pd.DataFrame(gene_loc_df.SYMBOL.tolist()).isin(clusteredgenes).any(1)]
    # clusteredgenes_df = clusteredgenes_df.drop(['annotation.GeneID', 'annotation.Strand', 
                                                # 'annotation.Length'], axis = 1)
    clusteredgenes_df = gene_loc_df[pd.DataFrame(gene_loc_df.SYMBOL.tolist()).isin(clusteredgenes).any(1)] # find clustered genes within gene location

    annotationstartlist = clusteredgenes_df['annotation.Start'].str.split(';') # split the annotation start column into a list
    clusteredgenes_df = clusteredgenes_df.drop('annotation.Start', axis = 1) # drop the original annotation start column 
    clusteredgenes_df['annotationstart'] = annotationstartlist # add the splitted annotation start column list
    clusteredgenes_df = clusteredgenes_df.explode('annotationstart') # covert the splitted annotation start list into rows


    annotationendlist = clusteredgenes_df['annotation.End'].str.split(';') 
    clusteredgenes_df = clusteredgenes_df.drop('annotation.End', axis = 1)
    clusteredgenes_df['annotationend'] = annotationendlist 
    clusteredgenes_df = clusteredgenes_df.explode('annotationend') 
 

    chromosomesplitlist = clusteredgenes_df['annotation.Chr'].str.split(';') 
    clusteredgenes_df = clusteredgenes_df.drop('annotation.Chr', axis = 1) 
    clusteredgenes_df['chromosome annotations'] = chromosomesplitlist 
    clusteredgenes_df = clusteredgenes_df.explode('chromosome annotations') 


    strandlist = clusteredgenes_df['annotation.Strand'].str.split(';') 
    clusteredgenes_df = clusteredgenes_df.drop('annotation.Strand', axis = 1) 
    clusteredgenes_df['annotationstrand'] = strandlist 
    clusteredgenes_df = clusteredgenes_df.explode('annotationstrand') 


    clusteredgenes_df = clusteredgenes_df.drop_duplicates()
    
    clusteredgenesstart_df = clusteredgenes_df.drop_duplicates(subset='annotationstart', keep='first')
    a = clusteredgenesstart_df['annotationstart']
    clusteredgenesend_df = clusteredgenes_df.drop_duplicates(subset='annotationend', keep='first')
    b = clusteredgenesend_df['annotationend']
    
    clusteredgenes_df = clusteredgenes_df.drop_duplicates(subset='annotationstart', keep='first')
    
    clusteredgenes_df = clusteredgenes_df.drop('annotationstart', axis = 1) 
    clusteredgenes_df = clusteredgenes_df.drop('annotationend', axis = 1)
    
    frames = [clusteredgenes_df, a, b]
    
    clusteredgenes_df = pd.concat(frames, axis=1)
    
    
    
    
    genes = np.asarray(clusteredgenes_df) # convert clusteredgenes dataframe to array

    genesstartingannotation = genes[:,5].astype(float) #get a list of starting annotation
    genesendingannotation = genes[:,6].astype(float) #get a list of ending annotation
    
    
    genestrand = genes[:,4].astype(str) #get a list of gene strand 
    geneschr = genes[:,3].astype(str) #get a list of gene chr
    genename = genes[:,0].astype(str) #get a list of genenames
    nucchr = nuc_data [:,1].astype(str) #get a list of chr for nuc_data to compare
    nucstrand = nuc_data[:,4] # get the list of nuclear strand to compare
    
    
    #determine the range of base pairs
    maxdistance = nuc_data[:,3] + 20
    mindistance = nuc_data[:,2] - 20
    maxd = maxdistance.astype(float)
    mind = mindistance.astype(float)
    

    #create empty arrays to put motifs in range
    x = [] #for the nuc data
    y = [] #for the gene cluster data

    #iterate through both gene and nuc data and find == strand, chromosome, range +- 20
    for j in range(len(nuc_data)):
        for i in range(len(genes)):
            
            # if  nucstrand[j] == genestrand[i] and geneschr[i] == nucchr[j] and ((maxd[j] >= genesstartingannotation[i] and mind[j] <= genesstartingannotation[i]) or (maxd[j] >= genesendingannotation[i] and mind[j] <= genesendingannotation[i])):
            if  nucstrand[j] == genestrand[i] and geneschr[i] == nucchr[j] and np.any(maxd[j] >= np.arange(genesstartingannotation[i], genesendingannotation[i])) and np.any(mind[j] <= np.arange(genesstartingannotation[i], genesendingannotation[i])):
                x.append(nuc_data[j])
                y.append(genename[i])
            else:
                x = x
    n = len(y)      
    y2 = np.reshape(y, (n,1)) #turning 1d array into 2d array
    print(x)
    print(y2)
    if n > 0:
        output_array = np.concatenate((x,y2), axis=1) #combine the arrays!
    else:
        output_array = ['There are no motifs :(']
    
    return output_array

clus1_output = find_motifs(clus1_genes_df,nuc_data)
np.savetxt('clus1output.txt', clus1_output, fmt='%s')

clus2_output = find_motifs(clus2_genes_df,nuc_data)
np.savetxt('clus2output.txt', clus2_output, fmt='%s')

clus3_output = find_motifs(clus3_genes_df,nuc_data)
np.savetxt('clus3output.txt', clus3_output, fmt='%s')

clus4_output = find_motifs(clus4_genes_df,nuc_data)
np.savetxt('clus4output.txt', clus4_output, fmt='%s')

clus5_output = find_motifs(clus5_genes_df,nuc_data)
np.savetxt('clus5output.txt', clus5_output, fmt='%s')

clus6_output = find_motifs(clus6_genes_df,nuc_data)
np.savetxt('clus6output.txt', clus6_output, fmt='%s')

clus7_output = find_motifs(clus7_genes_df,nuc_data)
np.savetxt('clus7output.txt', clus7_output, fmt='%s')

clus8_output = find_motifs(clus8_genes_df,nuc_data)
np.savetxt('clus8output.txt', clus8_output, fmt='%s')

clus9_output = find_motifs(clus9_genes_df,nuc_data)
np.savetxt('clus9output.txt', clus9_output, fmt='%s')



stop = timeit.default_timer()
execution_time = stop - start
print("Program Executed in "+str(execution_time))
print('Good Job Rob :)')