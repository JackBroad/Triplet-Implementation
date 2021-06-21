#!/bin/bash -f

rm tests/test_regression

make -f testsMakefile
mpirun --mca shmem posix --oversubscribe -np 2 tests/test_regression
