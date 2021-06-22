#!/bin/bash -f

rm tests/test_regression

make -f testsMakefile
mpirun -np 2 tests/test_regression
