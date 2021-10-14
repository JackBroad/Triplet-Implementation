To benchmark the speed increase with the number of processes of a given tmpi or toy function, use the following steps:

1) Set up the main.F90 file to load in the desired file of atomic positions and call the function(s) being tested.

2) Run './build_and_run.sh' to generate the 'triplet.out' file (if there are lots of atoms, use more processes for this step!)

3) Run './run_benchmarking_jobs.sh N_proc', where N_proc is the maximum desired number of processes (int), to run the jobs.

4) When all jobs from the previous step are complete, run './run_benchmarking_plots.sh N_proc storageDir', where N_proc is unchanged from (3) and
   storgaeDir is the name of the new directory where all data and plots will be stored (string).