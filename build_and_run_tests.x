#!/bin/bash -f

#rm testDir/tests

make -f testsMakefile

./testDir/tests
