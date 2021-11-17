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
  if [[ i -eq 1 ]]; then
    sed -ie 's/[0-9]\+ /'$i' /' myJob.sh;
    sbatch --partition=defq --nodes=1 --ntasks-per-node=$i --mem=4g --time=1:00:00 --output=$i-data.txt --error=$i-data.err myJob.sh
    if [[ $maxNodes -gt 1 ]]; then
      sed -ie 's/[0-9]\+ /'$num' /' myJob.sh;
      sbatch --partition=defq --nodes=$maxNodes --ntasks-per-node=$i --mem=8g --time=1:00:00 --output=$num-data.txt --error=$num-data.err myJob.sh
    fi
  elif [[ i -gt 1 ]]&& [[ i -le 20/$maxNodes ]]; then
    sed -ie 's/[0-9]\+ /'$num' /' myJob.sh;
    sbatch --partition=defq --nodes=$maxNodes --ntasks-per-node=$i --mem=80g --time=1:00:00 --output=$num-data.txt --error=$num-data.err myJob.sh
  else
    sed -ie 's/[0-9]\+ /'$num' /' myJob.sh;
    sbatch --partition=defq --nodes=$maxNodes --ntasks-per-node=$i --mem=160g --time=1:00:00 --output=$num-data.txt --error=$num-data.err myJob.sh
  fi
done
