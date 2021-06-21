#!/bin/bash -f

rm tests/test_regression

make -f testsMakefile

./tests/test_regression
