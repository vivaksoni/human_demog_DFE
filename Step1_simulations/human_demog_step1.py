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


#Example input
#python human_demog_msprime.py -seq_len 50000 -rec_rate 1 -mut_rate 2e-8 -num_replicates 100 \
#-T_AFR_EURASI 45  -T_EURASI_EUR 17.2 -T_EURASI_EAS 17.2 -T_EURASI_SAS 17.2 \
#-T_BANTU 5 -r_AFR 0.1 -r_EURASI 0.1 -r_EUR 0.1 -r_EAS 0.1 -r_SAS 0.1 \
#-N_AFR_ANC 4400 -B_EURASI 0.01 -B_EUR 0.01 -B_EAS 0.01 -B_SAS 0.01 \
#-m_AFR_EURAS 1 -m_AFR_EUR 1 -m_AFR_EAS 1 -m_AFR_SAS 1 -m_EURASI_EUR 1 -m_EURASI_EAS 1 -m_EURASI_SAS 1 -m_EUR_EAS 1 -m_EUR_SAS 1 -m_EAS_SAS 1 \
#-outPath "/home/vivak/human_demog_DFE/"

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-region', dest = 'region', action='store', nargs = 1, type = int, help = 'Region')
parser.add_argument('-seq_len', dest = 'seq_len', action='store', nargs = 1, type = float, help = 'Sequence length (bp)')
parser.add_argument('-rec_rate', dest = 'rec_rate', action='store', nargs = 1, type = float, help = 'Recombination rate (cM/Mb)')
parser.add_argument('-mut_rate', dest = 'mut_rate', action='store', nargs = 1, type = float, help = 'Mutation rate (per bp per generation)')
parser.add_argument('-num_replicates', dest = 'num_replicates', action='store', nargs = 1, type = int, help = 'Number of replicates')
parser.add_argument('-unmasked_len', dest = 'unmasked_len', action='store', nargs = 1, type = int, help = 'Unmasked sequence length (bp)$')
parser.add_argument('-T_AFR_EURASI', dest = 'T_AFR_EURASI', action='store', nargs = 1, type = float, help = 'Time of African-Eurasian split')
parser.add_argument('-T_EURASI_EUR', dest = 'T_EURASI_EUR', action='store', nargs = 1, type = float, help = 'Time of Eurasian-European split')
parser.add_argument('-T_EURASI_EAS', dest = 'T_EURASI_EAS', action='store', nargs = 1, type = float, help = 'Time of Eurasian-East Asian split')
parser.add_argument('-T_EURASI_SAS', dest = 'T_EURASI_SAS', action='store', nargs = 1, type = float, help = 'Time of Eurasian-South Asian split')
parser.add_argument('-T_BANTU', dest = 'T_BANTU', action='store', nargs = 1, type = float, help = 'Time of start of Bantu expansion')

parser.add_argument('-r_AFR', dest = 'r_AFR', action='store', nargs = 1, type = float, help = 'African growth rate')
parser.add_argument('-r_EURASI', dest = 'r_EURASI', action='store', nargs = 1, type = float, help = 'Eurasian growth rate')
parser.add_argument('-r_EUR', dest = 'r_EUR', action='store', nargs = 1, type = float, help = 'European growth rate')
parser.add_argument('-r_EAS', dest = 'r_EAS', action='store', nargs = 1, type = float, help = 'East Asian growth rate')
parser.add_argument('-r_SAS', dest = 'r_SAS', action='store', nargs = 1, type = float, help = 'South Asian growth rate')

parser.add_argument('-N_AFR_ANC', dest = 'N_AFR_ANC', action='store', nargs = 1, type = int, help = 'Initial population size')
parser.add_argument('-B_EURASI', dest = 'B_EURASI', action='store', nargs = 1, type = float, help = 'Severity of Eurasian bottleneck')
parser.add_argument('-B_EUR', dest = 'B_EUR', action='store', nargs = 1, type = float, help = 'Severity of European bottleneck')
parser.add_argument('-B_EAS', dest = 'B_EAS', action='store', nargs = 1, type = float, help = 'Severity of East Asian bottleneck')
parser.add_argument('-B_SAS', dest = 'B_SAS', action='store', nargs = 1, type = float, help = 'Severity of South Asian bottleneck')

parser.add_argument('-m_AFR_EURASI', dest = 'm_AFR_EURASI', action='store', nargs = 1, type = float, help = 'AFR-EURASI migration rate')
parser.add_argument('-m_AFR_EUR', dest = 'm_AFR_EUR', action='store', nargs = 1, type = float, help = 'AFR-EUR migration rate')
parser.add_argument('-m_AFR_EAS', dest = 'm_AFR_EAS', action='store', nargs = 1, type = float, help = 'AFR-EAS migration rate')
parser.add_argument('-m_AFR_SAS', dest = 'm_AFR_SAS', action='store', nargs = 1, type = float, help = 'AFR-SAS migration rate')
parser.add_argument('-m_EURASI_EUR', dest = 'm_EURASI_EUR', action='store', nargs = 1, type = float, help = 'EURASI-EUR migration rate')
parser.add_argument('-m_EURASI_EAS', dest = 'm_EURASI_EAS', action='store', nargs = 1, type = float, help = 'EURASI-EAS migration rate')
parser.add_argument('-m_EURASI_SAS', dest = 'm_EURASI_SAS', action='store', nargs = 1, type = float, help = 'EURASI-SAS migration rate')
parser.add_argument('-m_EUR_EAS', dest = 'm_EUR_EAS', action='store', nargs = 1, type = float, help = 'EUR-EAS migration rate')
parser.add_argument('-m_EUR_SAS', dest = 'm_EUR_SAS', action='store', nargs = 1, type = float, help = 'EUR-SAS migration rate')
parser.add_argument('-m_EAS_SAS', dest = 'm_EAS_SAS', action='store', nargs = 1, type = float, help = 'EAS-SAS migration rate')

parser.add_argument('-coords', dest = 'coords', action='store', nargs = 1, type = str, help = 'path to masking coords file')
parser.add_argument('-outPath', dest = 'outPath', action='store', nargs = 1, type = str, help = 'path to output files (suffixes will be added)')

args = parser.parse_args()
region = args.region[0]
seq_len = args.seq_len[0]
unmasked_len = args.unmasked_len[0]
rec_rate = args.rec_rate[0]
mut_rate = args.mut_rate[0]
num_replicates = args.num_replicates[0]
T_AFR_EURASI = args.T_AFR_EURASI[0]
T_EURASI_EUR = args.T_EURASI_EUR[0]
T_EURASI_EAS = args.T_EURASI_EAS[0]
T_EURASI_SAS = args.T_EURASI_SAS[0]
T_BANTU = args.T_BANTU[0] 
r_AFR = args.r_AFR[0]
r_EURASI = args.r_EURASI[0]
r_EUR = args.r_EUR[0]
r_EAS = args.r_EAS[0]
r_SAS = args.r_SAS[0]
N_AFR_ANC = args.N_AFR_ANC[0]
B_EURASI = args.B_EURASI[0]
B_EUR = args.B_EUR[0]
B_EAS = args.B_EAS[0]
B_SAS = args.B_SAS[0] 
m_AFR_EURASI = args.m_AFR_EURASI[0]
m_AFR_EUR = args.m_AFR_EUR[0]
m_AFR_EAS = args.m_AFR_EAS[0]
m_AFR_SAS = args.m_AFR_SAS[0]
m_EURASI_EUR = args.m_EURASI_EUR[0]
m_EURASI_EAS = args.m_EURASI_EAS[0]
m_EURASI_SAS = args.m_EURASI_SAS[0]
m_EUR_EAS = args.m_EUR_EAS[0]
m_EUR_SAS = args.m_EUR_SAS[0]
m_EAS_SAS = args.m_EAS_SAS[0]
coords = args.coords[0] 
outPath = args.outPath[0]

def human_demog_replicates(seq_len, rec_rate, mut_rate, num_replicates,
                           T_AFR_EURASI, T_EURASI_EUR, T_EURASI_EAS, T_EURASI_SAS, T_BANTU, 
                           r_AFR, r_EURASI, r_EUR, r_EAS, r_SAS, 
                           N_AFR_ANC, B_EURASI, B_EUR, B_EAS, B_SAS, 
                           m_AFR_EURASI, m_AFR_EUR, m_AFR_EAS, m_AFR_SAS, m_EURASI_EUR, m_EURASI_EAS, m_EURASI_SAS, m_EUR_EAS, m_EUR_SAS, m_EAS_SAS):
                           

    #Adjust rec rate to per site
    rec_rate = rec_rate * 1e-8
    # Times are provided in years, so we convert into generations.
    generation_time = 26.9
    T_AFR_EURASI = np.round((T_AFR_EURASI * 1e3) / generation_time,0)
    T_EURASI_EUR = np.round((T_EURASI_EUR * 1e3) / generation_time,0)
    T_EURASI_EAS = np.round((T_EURASI_EAS * 1e3) / generation_time,0)
    T_EURASI_SAS = np.round((T_EURASI_SAS * 1e3) / generation_time,0)
    T_BANTU = np.round((T_BANTU * 1e3) / generation_time,0)
    
    #Growth rates converted from %s
    r_AFR = r_AFR / 100
    r_EURASI = r_EURASI / 100
    r_EUR = r_EUR / 100
    r_EAS = r_EAS / 100
    r_SAS = r_SAS / 100
    
    #Calculate size of EURASI population at time when daughter populations split
    N_EURASI_EUR = (B_EURASI * N_AFR_ANC) / math.exp(-r_EURASI * (T_AFR_EURASI-T_EURASI_EUR))
    N_EURASI_EAS = (B_EURASI * N_AFR_ANC) / math.exp(-r_EURASI * (T_AFR_EURASI-T_EURASI_EAS))
    N_EURASI_SAS = (B_EURASI * N_AFR_ANC) / math.exp(-r_EURASI * (T_AFR_EURASI-T_EURASI_SAS))
    
    #If population at time of split is <1 set to 1.
    if(N_EURASI_EUR < 1):
        N_EURASI_EUR = 1
    if(N_EURASI_EAS < 1):
        N_EURASI_EAS = 1
    if(N_EURASI_SAS < 1):
        N_EURASI_SAS = 1
    
    #Initial population size (ie current day)
    N_AFR = N_AFR_ANC / math.exp(-r_AFR * T_BANTU)
    N_EURASI = (B_EURASI * N_AFR_ANC) / math.exp(-r_EURASI * T_AFR_EURASI)
    
    #If bottlenecks results in sizes less than 1, set to 1
    if ((B_EUR * N_EURASI_EUR) < 1):
        N_EUR = 1 / math.exp(-r_EUR * T_EURASI_EUR)
    else:
        N_EUR = (B_EUR * N_EURASI_EUR) / math.exp(-r_EUR * T_EURASI_EUR) 
    
    if ((B_EAS * N_EURASI_EAS) < 1):
        N_EAS = 1 / math.exp(-r_EAS * T_EURASI_EAS)
    else:
        N_EAS = (B_EAS * N_EURASI_EAS) / math.exp(-r_EAS * T_EURASI_EAS)
    
    if ((B_SAS * N_EURASI_SAS) < 1):
        N_SAS = 1 / math.exp(-r_SAS * T_EURASI_SAS)
    else:
        N_SAS = (B_SAS * N_EURASI_SAS) / math.exp(-r_SAS * T_EURASI_SAS)
    
    demography = msprime.Demography()
    demography.add_population(
        name="AFR",
        description="African population. Specifically LWK, from Kenya in 1kG",
        initial_size=N_AFR,
        growth_rate=r_AFR,
        default_sampling_time=0, 
        initially_active=True,
    )
    demography.add_population(
        name="EURASI",
        description=(
            "Eurasian population. This has no empirical data to compare to in 1kG"
        ),
        initial_size=N_EURASI,
        growth_rate=r_EURASI,
        default_sampling_time=0, 
        initially_active=True,
    )
    demography.add_population(
        name="EUR",
        description="European population. Comprised of CEU, TSI, GBR, FIN and IBS populations in 1kG",
        initial_size=N_EUR,
        growth_rate=r_EUR,
    )
    
    demography.add_population(
        name="EAS",
        description="East asian population. Comprised of CHB, JPT, CHS, CDX, KHV, and CHD populations in 1kG",
        initial_size=N_EAS,
        growth_rate=r_EAS,
    )
    
    demography.add_population(
        name="SAS",
        description="South Asian population. Comprised of GIH, PJL, BEB, STU, and ITU populations in 1kG",
        initial_size=N_SAS,
        growth_rate=r_SAS,
    )
    
    # Set the migration rates between extant populations
    demography.set_symmetric_migration_rate(["AFR", "EURASI"], m_AFR_EURASI / (4*N_AFR_ANC))
    demography.set_symmetric_migration_rate(["AFR", "EUR"], m_AFR_EUR / (4*N_AFR_ANC))
    demography.set_symmetric_migration_rate(["AFR", "EAS"], m_AFR_EAS / (4*N_AFR_ANC))
    demography.set_symmetric_migration_rate(["AFR", "SAS"], m_AFR_SAS / (4*N_AFR_ANC))
    demography.set_symmetric_migration_rate(["EURASI", "EUR"], m_EURASI_EUR / (4*N_AFR_ANC))
    demography.set_symmetric_migration_rate(["EURASI", "EAS"], m_EURASI_EAS / (4*N_AFR_ANC))
    demography.set_symmetric_migration_rate(["EURASI", "SAS"], m_EURASI_SAS / (4*N_AFR_ANC))
    demography.set_symmetric_migration_rate(["EUR", "EAS"], m_EUR_EAS / (4*N_AFR_ANC))
    demography.set_symmetric_migration_rate(["EUR", "SAS"], m_EUR_SAS / (4*N_AFR_ANC))
    demography.set_symmetric_migration_rate(["EAS", "SAS"], m_EAS_SAS / (4*N_AFR_ANC))
    
    #Add events
    demography.add_population_parameters_change(time=T_BANTU, population="AFR", growth_rate=0, initial_size=N_AFR_ANC)
    demography.add_population_split(time=T_EURASI_EUR, derived=["EUR"], ancestral="EURASI")
    demography.add_population_split(time=T_EURASI_EAS, derived=["EAS"], ancestral="EURASI")
    demography.add_population_split(time=T_EURASI_SAS, derived=["SAS"], ancestral="EURASI")
    demography.add_population_split(time=T_AFR_EURASI, derived=["EURASI"], ancestral="AFR")
    demography.sort_events()


    ancestry_reps = msprime.sim_ancestry(
        {"AFR": 99, "EUR": 502, "EAS": 104, "SAS": 489}, 
        demography=demography, 
        sequence_length = seq_len,
        recombination_rate = rec_rate,
        num_replicates=num_replicates)

    for ts in ancestry_reps:
            mutated_ts = msprime.sim_mutations(ts, rate=mut_rate)
            yield mutated_ts

#Function gets libseqeunce type matrices and divergence estimates (based on fixations) for each population
def get_libseq_matrix(ts):
    p1 = [0, 0, 0, 2, 2, 3]
    p2 = [2, 3, 4, 3, 4, 4]
    n = ['AFR_EUR', 'AFR_EAS', 'AFR_SAS', 'EUR_EAS', 'EUR_SAS', 'EAS_SAS']
    dnames = {'AFR':0, 'EUR':2, 'EAS':3, 'SAS':4}
    d_data = {}
    for key, value in dnames.items():
        lst = []
        for var in ts.variants(samples=ts.samples(population=[value])):
            lst.append([var.site.position, ''.join([str(x) for x in var.genotypes])])
        
        df = pd.DataFrame(lst)
        df = df[(df[1].str.contains("0")) & (df[1].str.contains("1"))]
        df.columns = ['pos','genotype']
        d_data[key] = df
    
    #Repeat for population pairs for Fst
    for i,j in enumerate(p1):
        lst = []
        for var in ts.variants(samples=np.concatenate((ts.samples(population=[j]), (ts.samples(population=[p2[i]]))))):
            lst.append([var.site.position, ''.join([str(x) for x in var.genotypes])])
        df = pd.DataFrame(lst)
        df = df[(df[1].str.contains("0")) & (df[1].str.contains("1"))]
        df.columns = ['pos','genotype']
        d_data[n[i]] = df
    return(d_data)

#Function to mask sites as per empirical masks (phastcons and accessibility)
def mask_sites(d_data, coords, length):
    coords = pd.read_csv(coords, names=['start','end'], sep='\t')
    for pop in ['AFR', 'EUR', 'EAS', 'SAS', 'AFR_EUR', 'AFR_EAS', 'AFR_SAS', 'EUR_EAS', 'EUR_SAS', 'EAS_SAS']:
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

#Function to get sliding window stats
def get_windowed_stats(d_data, length, region, rep):

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
    
    dnames = {'AFR':0, 'EUR':2, 'EAS':3, 'SAS':4}
    wdf = pd.DataFrame()
    for pop in dnames:
        sd = libsequence.SimData(d_data[pop])
        #define sliding windows
        wins = libsequence.Windows(sd,window_size=10000/length,step_len=5000/length,starting_pos=0.0,ending_pos=1.0)
        d = {'population': pop, 'region': region, 'replicate': rep}
        
        for k, win in enumerate(wins):
            d['window'] = k
            d = {**d, **get_polySIM_stats(win)}
            if len(win.pos()) >= 5: #LD stats are pairwise. If only 1 site exists, it'll show an error.
                d = {**d, **get_LD_stats(win)}
            else:
                d['meanrsq'] = 'NA'

            wdf = pd.concat([wdf, pd.DataFrame.from_dict(d, orient='index').T])
    wdf = wdf[['population', 'region', 'replicate', 'window', 'tajimasd', 'meanrsq']].reset_index(drop=True)
    return(wdf)

def master_function(ts, coords, region, rep, length, unmasked_length):
    d_data = get_libseq_matrix(ts)
    d_data = mask_sites(d_data, coords, length)
    d_res = get_S(d_data, region, rep, unmasked_length)
    d_res = get_perSite_Fst(d_data, d_res)    
    df = pd.DataFrame.from_dict(d_res, orient='index').T
    wdf = get_windowed_stats(d_data, length, region, rep)
    return(df, wdf)



ts_list = []
for replicate_index, ts in enumerate(human_demog_replicates(seq_len, rec_rate, mut_rate, num_replicates,
                           T_AFR_EURASI, T_EURASI_EUR, T_EURASI_EAS, T_EURASI_SAS, T_BANTU, 
                           r_AFR, r_EURASI, r_EUR, r_EAS, r_SAS, 
                           N_AFR_ANC, B_EURASI, B_EUR, B_EAS, B_SAS, 
                           m_AFR_EURASI, m_AFR_EUR, m_AFR_EAS, m_AFR_SAS, m_EURASI_EUR, m_EURASI_EAS, m_EURASI_SAS, m_EUR_EAS, m_EUR_SAS, m_EAS_SAS)):
    ts_list.append(ts)


df = pd.DataFrame()
wdf = pd.DataFrame()
for rep in range(0, num_replicates):
    df1, df2 = master_function(ts_list[rep], coords, region, rep, seq_len, unmasked_len)
    df = pd.concat([df, df1])
    wdf = pd.concat([wdf, df2])


df.to_csv(outPath + "perSite.txt", sep='\t', header=True, index=False)
wdf.to_csv(outPath + "windows.txt", sep='\t', header=True, index=False)
