#!/bin/bash -f
make -f testsMakefile
mpirun -np 1 tests/test_regression
mpirun -np 2 tests/test_regression
mpirun -np 3 tests/test_regression
mpirun -np 4 tests/test_regression
