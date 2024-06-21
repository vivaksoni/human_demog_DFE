import argparse
import sys
import pandas as pd
import math
import os
import numpy as np
from scipy import stats
import allel

parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-inPath', dest = 'inPath', action='store', nargs = 1, type = str, help = 'path to input files (suffixes will be added)')
parser.add_argument('-pop', dest = 'pop', action='store', nargs = 1, type = str, help = 'population')
parser.add_argument('-model', dest = 'model', action='store', nargs = 1, type = int, help = 'model')
parser.add_argument('-rep', dest = 'rep', action='store', nargs = 1, type = int, help = 'replicate')
parser.add_argument('-win_size', dest = 'win_size', action='store', nargs = 1, type = int, help = 'one tailed window size (bp)$')
parser.add_argument('-outPath', dest = 'outPath', action='store', nargs = 1, type = str, help = 'path to output files (suffixes will be added)')

args = parser.parse_args()
inPath = args.inPath[0]
pop = args.pop[0]
model = args.model[0]
rep = args.rep[0]
win_size = args.win_size[0]
outPath = args.outPath[0]


#Function to parse .ms file data into nested list of positions and genotypes
#Function takes as input open reading file object and no. of samples
def get_nested_data_list(f_ms, pop):
    pops = {'AFR':'p1', 'EUR':'p3', 'EAS':'p5', 'SAS':'p4'}
    samples = {"AFR": 198, "EUR": 1004, "EAS": 208, "SAS": 978}
    l_Pos = [] #list of positions of SNPs
    l_Genos = [] #list of alleles
    d_tmp = {} #dict to store individual allele info for each individual (values) at each site (keys)

    #positions on line 2
    pos_lines = [2]
    for position, line in enumerate(f_ms):
        if position in pos_lines:
            #store positions in list
            pos_list  = line.split()
    #Set file pointer to start of file
    f_ms.seek(0)

    i = 0
    #Loop through positions, storing in list
    for pos in pos_list[1:]:
        #Append position to l_Pos (after converting to float)
        l_Pos.append(float(pos))
        #Add dictionary key for each position, with empty value
        d_tmp[str(i)] = ""
        i += 1 
        
    
    #genotypes on line 3 onwards (use samples argument to determine length of file)
    g_lines = [x for x in range(3, samples[pop] + 4)]
    #Loop through lines (ie individuals)
    for position, line in enumerate(f_ms):
        if position in g_lines:
            #Remove newline character
            line1 = line.strip('\n')
            i = 0
            #For each individual, loop through each site, appending allele information for that individual to 
            #the site number in dict
            while i < len(line1):
                d_tmp[str(i)] = d_tmp[str(i)] + line1[i]
                i = i + 1

    f_ms.seek(0)

    #Create nested list of positions and genotypes
    l_data = [[j, d_tmp[str(i)]] for i,j in enumerate(l_Pos)]
    return(l_data)


regionLen = 198345
samples = {"AFR": 198, "EUR": 1004, "EAS": 208, "SAS": 978}
#inPath = r"/ospool/ap21/data/vsoni11/hdd_res/sweep_detection/DFE/sims/"
#outPath = r"/ospool/ap21/data/vsoni11/hdd_res/sweep_detection/DFE/h12/"

#Read in .ms file
f_ms = open(inPath + pop + "_rep" + str(rep) + "_model" + str(model) + ".ms", 'r')
l_data = get_nested_data_list(f_ms, pop)

snps = [x[1] for x in l_data]
snps = [[float(x) for x in y] for y in snps]
snps_pos = [int(np.round(x[0]*regionLen))-1 for x in l_data]

hapArr = np.zeros(shape=[regionLen,samples[pop]])

for i,j in enumerate(snps_pos):
    hapArr[j] = snps[i]

hapArr = allel.HaplotypeArray(hapArr, dtype='i1')

#Create df of snp postions and windows around snps
df = pd.DataFrame([[x, x-win_size, x+win_size] for x in snps_pos], columns=['snp_position', 'win_start', 'win_end'])
df['win_start'] = np.where(df.win_start < 0, 0, df.win_start)
df['win_end'] = np.where(df.win_end > regionLen, regionLen, df.win_end)

lst = []
#Loop through df, subsetting haploytype array and estimating H12 on subsetted array
for x in range(0, len(df)):
    subhap = hapArr[df.win_start[x]:df.win_end[x]]
    lst.append(allel.garud_h(subhap)[1])

df['H12'] = lst

df.to_csv(outPath + pop + "_rep" + str(rep) + "_model" + str(model) + "_" + str(int(win_size/500)) + "kb.h12", sep='\t', header=True, index=False)
