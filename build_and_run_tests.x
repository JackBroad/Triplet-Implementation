#!/bin/bash -f

rm tests/test_regression

make -f testsMakefile
mpirun -np 3 tests/test_regression
