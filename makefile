VPATH = ./ src/
# Find all source files, create a list of corresponding object files
SRCS=  assert.F90 DataStructures.F90 MPI_Variables.F90 GP_Variables.F90 triplet_mod.F90 tmpi_calcFullSimBoxEnergy.F90 tmpi_calcAtomMoveEnergy.F90  main.F90
OBJS=$(patsubst %.F90,%.o,$(SRCS))

# Ditto for mods (They will be in both lists)
MODS=$(wildcard mod*.F90)
MOD_OBJS=$(patsubst %.F90,%.o,$(MODS))

# Compiler/Linker settings
FC = mpif90
FCFLAGS =  -c -cpp  -Wall -Wextra -Wno-do-subscript  -DASSERTS -Wno-argument-mismatch #-fallow-argument-mismatch -Wconversion -Wno-unused-parameter -ffpe-trap=invalid -ffpe-trap=zero,overflow,underflow -fbacktrace -fdump-core -fcheck=bounds -Wno-tabs  #-fmax-errors=5
FLFLAGS = # -g -Wall -DDEBUG -Wextra -Wconversion  -ffpe-trap=invalid -ffpe-trap=zero,overflow,underflow -fbacktrace -fdump-core -fcheck=bounds   #-fmax-errors=5
PROGRAM = triplet.out
PRG_OBJ = $(PROGRAM).o

# make without parameters will make first target found.
default : $(PROGRAM)

# Compiler steps for all objects
$(OBJS) : %.o : %.F90
	$(FC) $(FCFLAGS) -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $^

# If something doesn't work right, have a 'make debug' to 
# show what each variable contains.
debug:
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJS)"
	@echo "MODS = $(MODS)"
	@echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "PRG_OBJ = $(PRG_OBJ)"

clean:
	rm -rf $(OBJS) $(PROGRAM) $(patsubst %.o,%.mod,$(MOD_OBJS)) *.mod
	$(RM) src/*.o src/*.mod src/*.a  src/*.inc
	$(RM) tests/*.o tests/*.mod tests/*.a  tests/*.inc
	$(RM) tests/test_regression.F90 


.PHONY: debug default clean

# Dependencies
#main.o:  testing_2CO2_Ar.o

# Main program depends on all modules
$(PRG_OBJ) : $(MOD_OBJS)

# Blocks and allocations depends on shared
mod_blocks.o mod_allocations.o : mod_shared.o
