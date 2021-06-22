make
#mpirun --mca shmem posix --oversubscribe -np  2 ./triplet.out
mpirun -np  2 ./triplet.out
