#script takes .ms file and calculates summary statistics across sliding windows using libsequence

#for d in "${demog[@]}"; do for m in "${model[@]}"; do for s in $(seq 1 6); do for n in "${N[@]}"; do for i in $(seq 1 200); do \
#python3 stats_sliding_window.py -msFile ../results/"$d"/"$m"/demog_DFE/100/hSap_DFE"$s"_rep"$i"_"$n".ms \
#-fixedFile ../results/"$d"/"$m"/demog_DFE/100/hSap_DFE"$s"_rep"$i"_"$n".fixed \
#-outFile ../summary_stats/"$d"/"$m"/demog_DFE/100/hSap_DFE \
#"$s"_rep"$i"_"$n".stats -winSize 2000 -stepSize 1000 -regionLen 85005 -samples 100 -N 10000 -burnIn 10; done; done; done; done; done

from __future__ import print_function
import libsequence
import sys
import pandas as pd
import math
import argparse
import vcf

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



df = pd.read_csv(r"/home/vivak/human_demog_DFE/GRCh37_functionalMasked_15kb.bed", sep='\t', header=0)
regions = list(df.region)
lengths = list(df.length)
starts = list(df.start)
ends = list(df.end)


for i,j in enumerate(regions):
    vcf_reader = vcf.Reader(filename=r'/home/vivak/human_demog_DFE/vcf/regions_only/4pops/region' + str(j) + '.vcf.gz')

    l_data = []
    for record in vcf_reader:
        try:
            l_data.append([(record.POS-starts[i])/lengths[i], ''.join([x['GT'] for x in record.samples]).replace('|','')])
        except Exception:
            break


    sd = libsequence.SimData(l_data)
    d = {'region': j}
    f = libsequence.Fst(sd,[99,502,104,489])
    d['Fst_hsm'] = f.hsm()
    d['Fst_slatkin'] = f.slatkin()
    d['Fst_hbk'] = f.hbk()
    df = pd.DataFrame.from_dict(d, orient='index').T
    if((i==0) & (pop=='all')):
        df.to_csv(r'/home/vivak/human_demog_DFE/empirical_Fst.txt', sep='\t', index=False, header=True, mode='a')
    else:
        df.to_csv(r'/home/vivak/human_demog_DFE/empirical_Fst.txt', sep='\t', index=False, header=False, mode='a')

    print ("Stats for region " + str(j) + " output to file")
