#!/bin/bash
set -e

# Extract arguments
maxProcessors=$1
name=$2

# Create directory to store benchmarking data
mkdir $name

# Loop over no. of procs and run job for each N_proc
for ((i=1; i<=$maxProcessors; i++)); do
  sed -ie 's/[0-50] /'$i' /' myJob.sh # Update myJob.sh
  sbatch --partition=defq --nodes=1 --ntasks-per-node=$i --mem=40g --time=1:00:00 --output=$i-data.txt myJob.sh
  mv $i-data.txt $name/$i-data.txt
done
