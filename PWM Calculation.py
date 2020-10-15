# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 21:45:04 2020

@author: rober
"""

import numpy as np
import pandas as pd
import timeit
import logomaker
import PIL
import matplotlib.pyplot as plt

start = timeit.default_timer()

#read in clusters
clus1 = pd.read_csv('clus1output.txt', delimiter=' ', header=None) 
clus2 = pd.read_csv('clus2output.txt', delimiter=' ', header=None)
clus3 = pd.read_csv('clus3output.txt', delimiter=' ', header=None)
clus4 = pd.read_csv('clus4output.txt', delimiter=' ', header=None)
clus5 = pd.read_csv('clus5output.txt', delimiter=' ', header=None)
clus6 = pd.read_csv('clus6output.txt', delimiter=' ', header=None)
clus7 = pd.read_csv('clus7output.txt', delimiter=' ', header=None)
clus8 = pd.read_csv('clus8output.txt', delimiter=' ', header=None)
clus9 = pd.read_csv('clus9output.txt', delimiter=' ', header=None)


# #read in aligned clusters
# clus4aligned = pd.read_csv('clus4_aligned.txt', delimiter=' ', header=None)
# clus4motifs = clus4aligned[6].str.split("\t13", n = 1, expand = True)


# clus6aligned = pd.read_csv('clus6_aligned.txt', delimiter=' ', header=None)
# clus6motifs = clus6aligned[6].str.split("\t13", n = 1, expand = True)


# clus7aligned = pd.read_csv('clus7_aligned.txt', delimiter=' ', header=None)
# clus7motifs1 = clus7aligned[6].str.split("\t13", n = 1, expand = True)
# clus7motifs2 = clus7aligned[7].str.split("\t13", n = 1, expand = True)
# clus7motifs = clus7motifs1.fillna('') + clus7motifs2.fillna('')


# clus8aligned = pd.read_csv('clus8_aligned.txt', delimiter=' ', header=None)
# clus8motifs = clus8aligned[6].str.split("\t13", n = 1, expand = True)



# clus9aligned = pd.read_csv('clus9_aligned.txt', delimiter=' ', header=None)
# clus9motifs = clus9aligned[6].str.split("\t13", n = 1, expand = True)
               


# checking motif length 
def checkequal(geneclusters):
    allstrlen = geneclusters[5].str.len()
    mean = allstrlen.mean()
    return mean 


# clus1strlen = checkequal(clus1)
# clus2strlen = checkequal(clus2)
# clus3strlen = checkequal(clus3)
# clus4strlen = checkequal(clus4)
# clus5strlen = checkequal(clus5)
# clus6strlen = checkequal(clus6)
# clus7strlen = checkequal(clus7)
# clus8strlen = checkequal(clus8)
# clus9strlen = checkequal(clus9)


# print(clus1strlen)
# print(clus2strlen)
# print(clus3strlen)
# print(clus4strlen)
# print(clus5strlen)
# print(clus6strlen)
# print(clus7strlen)
# print(clus8strlen)
# print(clus9strlen)


def split(word): 
    return [char for char in word]

def PWM(geneclusters):
    
    #create PWM of zeros
    
    string = (geneclusters.iloc[0,5])
    pwm = np.zeros((4, len(string)), dtype=int)

    
    # count the occurences of A,G,T,C

    for j in range(len(string)):
        for i in range(len(geneclusters)):
            
            listofgenes = np.array(split(geneclusters.iloc[:,5])) #split the sequences depending on their column
            nucleotides = split(listofgenes[i])
            
            if 'A' == nucleotides[j]:
                pwm[0,j] += 1
            if 'G' == nucleotides[j]:
                pwm[1,j] +=1
            if 'T' == nucleotides[j]:
                pwm[2,j] +=1
            if 'C' == nucleotides[j]:
                pwm[3,j] +=1
            
    pwm = pwm/(len(geneclusters))
    pwm_df = pd.DataFrame(pwm, index=('A', 'G', 'T', 'C'))
    pwm = pwm_df.transpose()
    return pwm


# saving pwms

clus1pwm = PWM(clus1)
print(clus1pwm)
clus1pwm.to_csv('clus1pwm.csv')
clus1logo = logomaker.Logo(clus1pwm, font_name = 'Arial Rounded MT Bold')


                
clus2pwm = PWM(clus2)
clus2pwm.to_csv('clus2pwm.csv')
clus2logo = logomaker.Logo(clus2pwm, font_name = 'Arial Rounded MT Bold')
                
clus3pwm = PWM(clus3)
clus3pwm.to_csv('clus3pwm.csv')
clus3logo = logomaker.Logo(clus3pwm, font_name = 'Arial Rounded MT Bold')
                
clus4pwm = PWM(clus4)
clus4pwm.to_csv('clus4pwm.csv')
clus4logo = logomaker.Logo(clus4pwm, font_name = 'Arial Rounded MT Bold')

clus5pwm = PWM(clus5)
clus5pwm.to_csv('clus5pwm.csv')
clus5logo = logomaker.Logo(clus5pwm, font_name = 'Arial Rounded MT Bold')

clus6pwm = PWM(clus6)
clus6pwm.to_csv('clus6pwm.csv')
clus6logo = logomaker.Logo(clus6pwm, font_name = 'Arial Rounded MT Bold')
                
clus7pwm = PWM(clus7)
clus7pwm.to_csv('clus7pwm.csv')
clus7logo = logomaker.Logo(clus7pwm, font_name = 'Arial Rounded MT Bold')

clus8pwm = PWM(clus8)
clus8pwm.to_csv('clus8pwm.csv')
clus8logo = logomaker.Logo(clus8pwm, font_name = 'Arial Rounded MT Bold')

clus9pwm = PWM(clus9)
clus9pwm.to_csv('clus9pwm.csv')
clus9logo = logomaker.Logo(clus9pwm, font_name = 'Arial Rounded MT Bold')





# saving pwms for aligned clusters

# clus4pwmaligned = PWM(clus4motifs)
# clus4pwmaligned.to_csv('clus4pwmaligned.csv')

# clus7pwmaligned = PWM(clus7motifs)
# clus7pwmaligned.to_csv('clus7pwmaligned.csv')

# clus6pwmaligned = PWM(clus6motifs)
# clus6pwmaligned.to_csv('clus6pwmaligned.csv')

# clus8pwmaligned = PWM(clus8motifs)
# clus8pwmaligned.to_csv('clus8pwmaligned.csv')

# clus9pwmaligned = PWM(clus9motifs)
# clus9pwmaligned.to_csv('clus9pwmaligned.csv')



stop = timeit.default_timer()
execution_time = stop - start

print("Program Executed in "+str(execution_time))