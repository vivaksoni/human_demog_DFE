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
tar -zxf demog_DFE_div_sims.tar.gz
. $ENVDIR/bin/activate

# modify environment variables
export PATH=$_CONDOR_SCRATCH_DIR/build:$PATH

# untar software
#tar -xzf dadi_mff.tar.gz
#cd dadi_sfv


# modify this line to run your desired Python script and any other work you need to do
python3 hdd_DFE_stats.py \
-inPath "demog_DFE_div_sims" \
-region $1 \
-exon_len $2 \
-outPath "./region$1"

cp *.stats .
