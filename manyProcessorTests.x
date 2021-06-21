#!/bin/bash -f
make -f testsMakefile
mpirun --mca shmem posix --oversubscribe -np 1 tests/test_regression
mpirun --mca shmem posix --oversubscribe -np 2 tests/test_regression
mpirun --mca shmem posix --oversubscribe -np 3 tests/test_regression
mpirun --mca shmem posix --oversubscribe -np 4 tests/test_regression
mpirun --mca shmem posix --oversubscribe -np 6 tests/test_regression
mpirun --mca shmem posix --oversubscribe -np 8 tests/test_regression
mpirun --mca shmem posix --oversubscribe -np 12 tests/test_regression
#mpirun --mca shmem posix --oversubscribe -np 16 tests/test_regression

