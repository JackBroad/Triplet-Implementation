To benchmark the speed increase with the number of processes of a given tmpi or toy function, use the following steps:

1) Set up the main.F90 file to load in the desired file of atomic positions and call the function(s) being tested.

2) Run './build_and_run.sh' to generate the 'triplet.out' file (if there are lots of atoms, use more processes for this step!)

3) Run './run_benchmarking_jobs.sh N_proc storageDir', where N_proc is the maximum desired number of processes (int) and storageDir is the name of the 
   directory in which data will be stored, to run the jobs.

4) When all from the previous step are complete, run './run_benchmarking_plots.sh N_proc storageDir', where N_proc and storageDir are exactly as they 
   were in (3), to generate the plots and stor them in storageDir.
