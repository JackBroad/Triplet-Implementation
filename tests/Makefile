SHELL := /bin/bash

FC := mpif90


LATEST_PFUNIT_DIR=$(PFUNIT_DIR)
#LATEST_PFUNIT_DIR := /Users/pmzrsg/pFUnit/build/installed/PFUNIT-4.1
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk

all: test_regression

FFLAGS =  -c -cpp #-Wall -DDEBUG -Wextra -Wconversion  -ffpe-trap=invalid -ffpe-trap=zero,overflow,underflow -fbacktrace -fdump-core -fcheck=bounds  -Wno-tabs #-fmax-errors=5

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I../src

test_regression_TESTS := test_regression.pf 
test_regression_OTHER_LIBRARIES := -L../src -lsut
$(eval $(call make_pfunit_test,test_regression))



error_handling_driver.o: throw_with_pfunit.o

clean:
	echo $(LATEST_PFUNIT_DIR)
	$(RM) *.o *.mod *.a  *.inc
	$(RM) test_regression.F90 

