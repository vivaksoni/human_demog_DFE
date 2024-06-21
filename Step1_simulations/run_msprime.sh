#!/bin/bash

# have job exit if any command returns with non-zero exit status (aka failure)
set -e

# replace env-name on the right hand side of this line with the name of your conda environment
ENVNAME=env5
# if you need the environment directory to be named something other than the environment name, change this line
ENVDIR=$ENVNAME

# these lines handle setting up the environment; you shouldn't have to modify them
export PATH
mkdir $ENVDIR
tar -xzf $ENVNAME.tar.gz -C $ENVDIR
tar -xzf simulation_masks.tar.gz
. $ENVDIR/bin/activate

# modify environment variables
export PATH=$_CONDOR_SCRATCH_DIR/build:$PATH

# untar software
#tar -xzf dadi_mff.tar.gz
#cd dadi_sfv


# modify this line to run your desired Python script and any other work you need to do
python3 human_demog_step1.py \
-region $1 \
-seq_len $2 \
-unmasked_len ${30} \
-rec_rate $3 \
-mut_rate $4 \
-num_replicates 100 \
-T_AFR_EURASI $5  \
-T_EURASI_EUR $6 \
-T_EURASI_EAS $7 \
-T_EURASI_SAS $8 \
-T_BANTU $9 \
-r_AFR ${10} \
-r_EURASI ${11} \
-r_EUR ${12} \
-r_EAS ${13} \
-r_SAS ${14} \
-N_AFR_ANC ${15} \
-B_EURASI ${16} \
-B_EUR ${17} \
-B_EAS ${18} \
-B_SAS ${19} \
-m_AFR_EURAS ${20} \
-m_AFR_EUR ${21} \
-m_AFR_EAS ${22} \
-m_AFR_SAS ${23} \
-m_EURASI_EUR ${24} \
-m_EURASI_EAS ${25} \
-m_EURASI_SAS ${26} \
-m_EUR_EAS ${27} \
-m_EUR_SAS ${28} \
-m_EAS_SAS ${29} \
-coords "simulation_masks/region$1.coords" \
-outPath "./run${31}_region$1_"

cp *.txt .
