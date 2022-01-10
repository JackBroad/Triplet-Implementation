#!/bin/bash
set -e

# Extract arguments
maxProcessors=$1
maxNodes=$2

# Set constants
iMax=$((maxProcessors/maxNodes))

# Loop over no. of procs and run job for each N_proc
for ((i=1; i<=$iMax; i++)); do
  num=$((i*maxNodes)) # Wrap any arithmetic expressions in parentheses
  maxMem=$((i*8))
  if [[ i -eq 1 ]]; then
    sed -ie 's/[0-9]\+ /'$i' /' myJob.sh;
    sbatch --partition=mmemq --nodes=1 --ntasks-per-node=$i --mem=8g --time=1:00:00 --output=$i-data.txt --error=$i-data.err myJob.sh
    if [[ $maxNodes -gt 1 ]]; then
      sed -ie 's/[0-9]\+ /'$num' /' myJob.sh;
      sbatch --partition=mmemq --nodes=$maxNodes --ntasks-per-node=$i --mem=$((maxMem))g --time=1:00:00 --output=$num-data.txt --error=$num-data.err myJob.sh
    fi
  else
    sed -ie 's/[0-9]\+ /'$num' /' myJob.sh;
    sbatch --partition=mmemq --nodes=$maxNodes --ntasks-per-node=$i --mem=$((maxMem))g --time=1:00:00 --output=$num-data.txt --error=$num-data.err myJob.sh
  fi
done
