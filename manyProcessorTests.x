#!/bin/bash -f
make -f testsMakefile
mpirun -np 1 tests/test_regression
mpirun -np 2 tests/test_regression
mpirun -np 3 tests/test_regression
mpirun -np 4 tests/test_regression
mpirun -np 5 tests/test_regression
mpirun -np 6 tests/test_regression
mpirun -np 7 tests/test_regression
mpirun -np 8 tests/test_regression
