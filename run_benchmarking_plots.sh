#!/bin/bash
set -e

# Extract arguments
maxProcessors=$1
dirName=$2

# Create directory to store benchmarking data
mkdir $dirName

# Load modules for plotting in python
module load Python/2.7.10-gimkl-2.11.5
module load matplotlib/3.2.1-intel-2020a-Python-3.8.2

# Call the plotting script
python benchmarkPlots.py $maxProcessors $dirName

# Unload all modules and re-load those used by the tmpi code only
module purge
source UoN_HPC_SetUp.sh

# Move .txt files to storage dir and remove .err files
for ((i=1; i<=$maxProcessors; i++)); do
  mv $i-data.txt $dirName/$i-data.txt
  rm $i-data.err
done
