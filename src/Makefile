FC :=  mpif90
SRCS := MPI_Variables.F90 GP_Variables.F90 triplet_mod.F90 triplet_mpi.F90
OBJS := $(SRCS:%.F90=%.o)

MODS=$(wildcard mod*.f90)
MOD_OBJS=$(patsubst %.f90,%.o,$(MODS))

FCFLAGS =  -c -cpp -Wno-argument-mismatch #-c -cpp -DDEBUG -DASSERTS -Wall -DDEBUG -Wextra -Wconversion  -ffpe-trap=invalid -ffpe-trap=zero,overflow,underflow -fbacktrace -fdump-core -fcheck=bounds  -Wno-tabs  #-fmax-errors=5

FLFLAGS = # -ffpe-trap=invalid -ffpe-trap=zero,overflow,underflow -fbacktrace -fdump-core -fcheck=bounds #-fmax-errors=5

all: libsut.a

libsut.a: $(OBJS)
	$(AR) $(FLFLAGS) -r $@ $?

%.o : %.F90
	$(FC) $(FCFLAGS)  -c $(FFLAGS) $<

clean:
	$(RM) *.o *.mod *.a

