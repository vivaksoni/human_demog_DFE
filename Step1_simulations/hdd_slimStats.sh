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
tar -zxf demog_sims_div.tar.gz
. $ENVDIR/bin/activate

# modify environment variables
export PATH=$_CONDOR_SCRATCH_DIR/build:$PATH

# untar software
#tar -xzf dadi_mff.tar.gz
#cd dadi_sfv


# modify this line to run your desired Python script and any other work you need to do
python3 hdd_slimStats.py \
-inPath "demog_sims_div" \
-region $1 \
-seq_len $2 \
-unmasked_len $3 \
-coords "simulation_masks/region$1.coords" \
-outPath "./region$1"

cp *.txt .
