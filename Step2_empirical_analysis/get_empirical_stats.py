import libsequence
import sys
import pandas as pd
import math
import argparse
import vcf
import numpy as np

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


#Read in region information
df = pd.read_csv(r"/home/vivak/human_demog_DFE/GRCh37_functionalMasked_15kb.bed", sep='\t', header=0)
regions = list(df.region)
lengths = list(df.length)
starts = list(df.start)
ends = list(df.end)
u_lengths = list(df.unmasked_length)


#Create empty dfs to store summary stats (wdf for windowed stats)
sdf = pd.DataFrame()
wdf = pd.DataFrame()
#Loop through populations
pops = {'LWK':'masked/LWK/', 'EUR':'masked/EUR/', 'EAS':'masked/EAS/', 'SAS':'masked/SAS/'}
for pop in pops:
    #Loop through regions
    for i,j in enumerate(regions):
        #Read in vcf data and store positions and genotypes in pylibseq friendly format
        vcf_reader = vcf.Reader(filename=r'/home/vivak/human_demog_DFE/vcf/regions_only/' + pops[pop] + 'region' + str(j) + '.vcf.gz')
        l_data = []
        for record in vcf_reader:
            l_data.append([(record.POS-starts[i])/lengths[i], ''.join([x['GT'] for x in record.samples]).replace('|','')])
            l.append(''.join([x['GT'] for x in record.samples]).replace('|',''))
        #Run libsequence to calculate summary stats
        sd = libsequence.SimData(l_data)
        ss = get_polySIM_stats(sd)
        d = {'population': pop, 'region': j}
        d = {**d, **get_polySIM_stats(sd)}
        df2 = pd.DataFrame.from_dict(d, orient='index').T
        df2 = df2[['population', 'region', 'numpoly','numsingletons', 'thetapi']]
        for stat in ['numpoly','numsingletons', 'thetapi']:
            df2[stat] = df2[stat] / u_lengths[i]
        sdf = pd.concat([sdf, df2])

        #Repeat but for 10kb windows for tajima's D and rsq
        #define sliding windows
        wins = libsequence.Windows(sd,window_size=10000/lengths[i],step_len=5000/lengths[i],starting_pos=0.0,ending_pos=1.0)
        num_wins = len(wins)
        d = {'population': pop, 'region': j}
        
        for k, win in enumerate(wins):
            d['window'] = k
            ps = libsequence.PolySIM(win)
            d['tajimasd'] = ps.tajimasd()
            if len(win.pos()) >= 5: #LD stats are pairwise. If only 1 site exists, it'll show an error.
                d = {**d, **get_LD_stats(win)}
            else:
                d['meanrsq'] = 'NA'

            df3 = pd.DataFrame.from_dict(d, orient='index').T
            wdf = pd.concat([wdf, df3])

wdf['meanrsq'] = wdf['meanrsq'].mask(wdf['meanrsq'] == "NA")
rdf = pd.merge(sdf, wdf.groupby(['population', 'region']).mean().reset_index().drop(columns=['window']), on=['population','region'], how='inner')

#Read in Fst stas
fdf = pd.read_csv(r'/home/vivak/human_demog_DFE/empirical_summary_stats/empirical_Fst.txt', sep='\t', header=0)

#Rearrange dfs to order by regions
edf = rdf[rdf.population=='LWK']
edf = edf.drop(columns='population')
edf.columns = ['region', 'S_AFR', 'singletons_AFR', 'pi_AFR', 'tajimasd_AFR', 'meanrsq_AFR', ]

for pop in ['EUR', 'EAS', 'SAS']:
    tdf = rdf[rdf.population==pop]
    tdf = tdf.drop(columns='population')
    tdf.columns = ['region', 'S_'+pop, 'singletons_'+pop, 'pi_'+pop, 'tajimasd_'+pop, 'meanrsq_'+pop]
    edf = pd.merge(edf, tdf, on='region', how='inner')

edf3 = fdf[fdf['populations']=='LWK_EUR']
edf3 = edf3.drop(columns=['populations'])
edf3.columns = ['region', 'Fst_AFR_EUR']

#Repeat for Fst data
pops = {'LWK_EAS':'AFR_EAS', 'LWK_SAS':'AFR_SAS', 'EAS_EUR': 'EUR_EAS', 'SAS_EUR': 'EUR_SAS', 'EAS_SAS': 'EAS_SAS'}
for pop in pops:
    tdf = fdf[fdf['populations']==pop]
    tdf = tdf.drop(columns=['populations'])
    tdf.columns = ['region', 'Fst_' + pops[pop]]
    edf3 = pd.merge(edf3, tdf, on='region', how='inner')

#Merge into single df and output to file
edf = pd.merge(edf, edf3, on='region', how='inner')
edf.to_csv(r'/home/vivak/human_demog_DFE/empirical_summary_stats/empirical_summary_stats_final.txt', header=True, index=False, sep='\t')
