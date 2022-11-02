To benchmark the speed-up that accompanies increasing the number of processes for a given tmpi or toy function, use the following steps:

1) Set up the main.F90 file to load in the desired file of atomic positions and call the function(s) being tested.

2) Run './build_and_run.sh' to generate the 'triplet.out' file. (If there are lots of atoms, use more processes for this step! The number of processes 
   can be editied in 'build_and_run.sh' by changing the number in the file.)

3) Run './run_benchmarking_jobs.sh N_proc N_nodes', where N_proc is the maximum desired number of processes and N_nodes is the number of different nodes
   to spread them over (both int); this generates (N_proc/N_nodes +1) .txt files (containing times) and the same number of .err files (containing any 
   slurm error messages).

4) When all jobs from (3) are complete, run './run_benchmarking_plots.sh N_proc N_nodes storageDir', where N_proc is unchanged from (3), as is N_nodes,
   and storageDir is the name of the new directory in which all data and plots will be stored (string); this script also removes all of the .err files
   produced during step (3).

Any additional scripts to be called from 'main.F90' should be written in the 'Triplet-Implementation/src' directory and added to the makefiles in 'src'
and 'Triplet-Implementation'.
