#!/bin/bash
set -e

# Extract arguments
maxProcessors=$1
name=$2

# Load modules for plotting in python
module load Python/2.7.10-gimkl-2.11.5
module load matplotlib/3.2.1-intel-2020a-Python-3.8.2

# Call the plotting script
python benchmarkPlots.py $maxProcessors $name

# Unload all modules and re-load those usedby the tmpi code only
module purge
source UoN_HPC_SetUp.sh
