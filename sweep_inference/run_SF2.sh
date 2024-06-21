#!/bin/bash
# slim_DFE.sh
set -x
# untar software
tar -xzf SF2.tar.gz
tar -xzf SF_files.tar.gz

# modify environment variables 
export PATH=$_CONDOR_SCRATCH_DIR/SF2:$PATH
#rmap=rr_mu_project2/rr_maps/10kb_mean/$6.txt

SweepFinder2 \
-lru \
SF_files/$1/$2_rep$4_model$3.grid \
SF_files/$1/$2_rep$4_model$3.aff \
SF_files/DFE2_$2.sfs \
SF_files/$1/$2_rep$4_model$3.rec \
./$1_$2_rep$4_model$3.clr

#sf_expansion/rr_variable_mu_fixed/benProp_0.0005_2Nes_100/1Mb_rep$1.aff \
#sf_expansion/rr_variable_mu_fixed/benProp_0.0005_2Nes_100/1Mb_rep$1.sfs \
#sf_expansion/rr_variable_mu_fixed/benProp_0.0005_2Nes_100/1Mb_rep$1.rec \
#./expansion_rr_variable_mu_fixed_benProp_0.0005_2Nes_100_1Mb_rep$1.out
cp *.clr .
