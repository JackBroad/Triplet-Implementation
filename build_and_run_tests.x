#!/bin/bash -f

rm tests/test_regression

make -f testsMakefile
mpirun -np 1 tests/test_regression
