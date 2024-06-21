import sys
import pandas as pd
import math
import os
import numpy as np
import msprime
import tskit
import libsequence
import allel
import argparse

parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-inPath', dest = 'inPath', action='store', nargs = 1, type = str, help = 'path to input files (suffixes will be added)')
parser.add_argument('-region', dest = 'region', action='store', nargs = 1, type = int, help = 'Region')
parser.add_argument('-seq_len', dest = 'seq_len', action='store', nargs = 1, type = float, help = 'Sequence length (bp)')
parser.add_argument('-unmasked_len', dest = 'unmasked_len', action='store', nargs = 1, type = int, help = 'Unmasked sequence length (bp)$')
parser.add_argument('-coords', dest = 'coords', action='store', nargs = 1, type = str, help = 'path to masking coords file')
parser.add_argument('-outPath', dest = 'outPath', action='store', nargs = 1, type = str, help = 'path to output files (suffixes will be added)')

args = parser.parse_args()
inPath = args.inPath[0]
region = args.region[0]
seq_len = args.seq_len[0]
unmasked_len = args.unmasked_len[0]
coords = args.coords[0] 
outPath = args.outPath[0]

#Function to parse .ms file data into nested list of positions and genotypes
#Function takes as input open reading file object and no. of samples
def get_nested_data_list(inPath, region, rep, chr_len):
    p1 = ['AFR', 'AFR', 'AFR', 'EUR', 'EUR', 'EAS']
    p2 = ['EUR', 'EAS', 'SAS', 'EAS', 'SAS', 'SAS']
    d_pops = {"AFR": 198, "EUR": 1004, "EAS": 208, "SAS": 978}
    d_res = {}
    d_fix = {}
    for pop in d_pops:
        f_ms = open(inPath + "/" + pop + "_region" + str(region) + "_rep" + str(rep) + ".ms", 'r')
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
        g_lines = [x for x in range(3, d_pops[pop] + 4)]
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
        df = pd.DataFrame(l_data)
        df = df[(df[1].str.contains("0")) & (df[1].str.contains("1"))]
        df.columns = ['pos','genotype']
        df['pos'] = (df.pos * chr_len)
        d_res[pop] = df

    for i,pop in enumerate(p1):
        df = pd.merge(d_res[p1[i]], d_res[p2[i]], on='pos', how='outer')
        df['genotype_x'] = np.where(df.genotype_x.isna(), '0'*d_pops[p1[i]], df.genotype_x)
        df['genotype_y'] = np.where(df.genotype_y.isna(), '0'*d_pops[p2[i]], df.genotype_y)
        df['genotype'] = df.genotype_x + df.genotype_y
        df = df[['pos', 'genotype']]
        d_res[p1[i] + "_" + p2[i]] = df
    return(d_res)

#Function to mask sites as per empirical masks (phastcons and accessibility)
def mask_sites(d_data, coords, length):
    coords = pd.read_csv(coords, names=['start','end'], sep='\t')
    for pop in d_data.keys():
        df = d_data[pop]
        #Convert masking df into array
        arr = coords.to_numpy()
        #:Use np.ma.mask_where to mask regions then add masking status as column in df
        x = np.array(df.pos)
        for l, u in arr:
            x = np.ma.masked_where((x > l) & (x < u), x)
        df['mask'] = x
        #Remove masked sites
        df = df.dropna()
        #Convert positions to relative qty
        l = df['pos'] / length
        df = df.drop(columns=['pos'])
        df['pos'] = l
        df = df[['pos','genotype']]
        d_data[pop] = np.array(df)
    return(d_data)

#Function to get per site summary statistics from output of get_libseq_matrix
def get_S(d_data, region, rep, unmasked_length):
    d_res = {'region' : region, 'replicate' : rep}
    #Get per site summary statistics
    for pop in ['AFR', 'EUR', 'EAS', 'SAS']:
        S = len(d_data[pop]) / unmasked_length
        d_res['S_' + pop] = S
    return(d_res)

def get_perSite_Fst(d_data, d_res):
    #Create for population names, and one for sample sizes
    d_samples1 = [99, 99, 99, 502, 502, 104]
    d_samples2 = [502, 104, 489, 104, 489, 489]
    #Loop through populations
    for p, pop in enumerate(['AFR_EUR', 'AFR_EAS', 'AFR_SAS', 'EUR_EAS', 'EUR_SAS', 'EAS_SAS']):
        gn = [x[1] for x in d_data[pop]]
        for i,j in enumerate(gn):
            gn[i] = list(j)
            gn[i] = [int(x) for x in gn[i]]
        
        h = allel.HaplotypeArray(gn, dtype='i1')
        h.n_haplotypes
        
        res = []
        #Loop through variants
        for hap in h:
            lst = []
            #Loop through samples combining genotypes for each individual into single list (thereby making a 3d lst which can be converted to gt array)
            for i in range(0, len(hap)-1, 2):
                lst.append([hap[i], hap[i+1]])
            res.append(lst)
        #Convert to gt array and calculate weir and cockerham Fst
        g = allel.GenotypeArray(res)
        a, b, c = allel.weir_cockerham_fst(g, [[x for x in range(0, d_samples1[p])], [x for x in range(d_samples1[p], d_samples1[p]+d_samples2[p])]])
        Fst = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))
        d_res['Fst_' + pop] = Fst
    return(d_res)

#Function to create dictionary of polySIM summary statistics. Returns dictionary of summary stats
def get_polySIM_stats(sd):
    #Create polysim object
    ps = libsequence.PolySIM(sd)
    #Create list of methods (ie polySIM summaryStats)
    a = [method for method in dir(ps) if callable(getattr(ps, method)) if not method.startswith('_')]
    #Loop through methods, storing names as keys, and estimates as values in dict
    ss_dict = {}
    for method in a:
        ss_dict[method] = getattr(ps, method)()
        
    return(ss_dict)

#Function to create dictionary of LD stats. Returns dictionary.
def get_LD_stats(sd):
    ld = libsequence.ld(sd)
    df = pd.DataFrame(ld)
    ss_dict = {}
    ss_dict['meanrsq'] = sum(df['rsq'])/len(df['rsq'])
    return(ss_dict)

#Function to get all summary stats except for divergence and Fst
def get_summary_stats(l_data, region, rep, length):
    df = pd.DataFrame()
    for pop in ['AFR', 'EUR', 'EAS', 'SAS']:
        sd = libsequence.SimData(l_data[pop])
        #ss = libsequence.PolySIM(sd)
        d = {'population': pop, 'region': region, 'rep':rep}
        d = {**d, **get_polySIM_stats(sd)}
        if len(sd.pos()) >= 5: #LD stats are pairwise. If only 1 site exists, it'll show an error.
            d = {**d, **get_LD_stats(sd)}
        else:
            d['meanrsq'] = 'NA'
        df2 = pd.DataFrame.from_dict(d, orient='index').T
        df2 = df2[['population', 'region', 'rep', 'numpoly','numsingletons', 'thetapi', 'tajimasd', 'meanrsq']]
        for stat in ['numpoly','numsingletons', 'thetapi']:
            df2[stat] = df2[stat] / length
        df = pd.concat([df, df2])
    #invert df columns    
    rdf = df[df.population=='AFR']
    rdf = rdf.drop(columns='population')
    rdf.columns = ['region', 'rep', 'S_AFR', 'singletons_AFR', 'pi_AFR', 'tajimasd_AFR', 'meanrsq_AFR']

    for pop in ['EUR', 'EAS', 'SAS']:
        tdf = df[df.population==pop]
        tdf = tdf.drop(columns='population')
        tdf.columns = ['region', 'rep', 'S_'+pop, 'singletons_'+pop, 'pi_'+pop, 'tajimasd_'+pop, 'meanrsq_'+pop]
        rdf = pd.merge(rdf, tdf, on=['region', 'rep'], how='inner')

    return(rdf)

def get_divergence(inPath, region, rep, length, coords):
    coords = pd.read_csv(coords, names=['start','end'], sep='\t')
    pops = {'AFR':'p1', 'EUR':'p3', 'EAS':'p5', 'SAS':'p4'}
    res = []
    for pop in pops:
        f = pd.read_csv(inPath + "/" + pop + "_region" + str(region) + "_rep" + str(rep) + ".fixed", sep=' ', 
                names=['OUT', 'tick', 'cycle', 'T', 'subpopID', 'mutID', 'mutType', 'pos', 's', 'h', 'originPop', 'originTick', 'AC'])

        #Subset for only relevant population and mutations that occur after burnin
        f = f[(f.subpopID == pops[pop]) & (f.originTick>=94000)]
        #Keep only position column
        f = f[['pos']]
        #Convert masking df into array
        arr = coords.to_numpy()
        #:Use np.ma.mask_where to mask regions then add masking status as column in df
        x = np.array(f.pos)
        for l, u in arr:
            x = np.ma.masked_where((x > l) & (x < u), x)
        f['mask'] = x
        #Remove masked sites
        f = f.dropna()
        #Calculate divergence
        div = len(f)/(length)
        res.append([pop,region,rep,div])
    df = pd.DataFrame(res)
    df.columns=['population', 'region', 'rep', 'divergence']
    
    #Invert df
    rdf = df[df.population=='AFR']
    rdf = rdf.drop(columns='population')
    rdf.columns = ['region', 'rep', 'divergence_AFR']
    
    for pop in ['EUR', 'EAS', 'SAS']:
        tdf = df[df.population==pop]
        tdf = tdf.drop(columns='population')
        tdf.columns = ['region', 'rep', 'divergence_'+pop]
        rdf = pd.merge(rdf, tdf, on=['region', 'rep'], how='inner')
    return(rdf)


def master_function(inPath, coords, region, rep, length, unmasked_length):
    d_data = get_nested_data_list(inPath, region, rep, length)
    d_data = mask_sites(d_data, coords, length)
    d_res = get_summary_stats(d_data, region, rep, unmasked_length)
    d_res = get_perSite_Fst(d_data, d_res)
    div = get_divergence(inPath, region, rep, unmasked_length, coords)
    df = pd.merge(d_res, div, on=['region','rep'], how='inner')
    return(df)


df = pd.DataFrame()
for rep in range(1, 101):
    df1 = master_function(inPath, coords, region, rep, seq_len, unmasked_len)
    df = pd.concat([df, df1])

df.to_csv(outPath + ".stats", sep='\t', header=True, index=False)
