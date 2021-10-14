#!/bin/bash
set -e

# Extract arguments
maxProcessors=$1

# Loop over no. of procs and run job for each N_proc
for ((i=1; i<=$maxProcessors; i++)); do
  sed -ie 's/[0-9]\+ /'$i' /' myJob.sh; # Update myJob.sh
  #sbatch --partition=defq --nodes=1 --ntasks-per-node=$i --mem=40g --time=1:00:00 --output=$i-data.txt myJob.sh
done
