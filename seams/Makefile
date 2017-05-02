# ======================================================================
# Makefile for seams module
# Author: Luis Garreta
# Date  : Jun 02/2013
# ======================================================================

# The compiler
FC = gfortran
FCFLAGS = -g -fbounds-check -fbacktrace
#FCFLAGS = -O2
#FC = g95 -O2
# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -O2
# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

SRC = .
OBJ = obj
BIN = bin

# List of executables to be built within the package
PROGRAMS = xpdb2seams xhb2cij

# "make" builds all
all: dir $(PROGRAMS)

dir:
	mkdir -p obj

xhb2cij: hb2cij.f90 contact_obj.o
	$(FC) $(FCFLAGS) -I$(OBJ) -J$(OBJ) -o $@ $(LDFLAGS) hb2cij.f90 $(OBJ)/contact_obj.o

contact_obj.o: $(SRC)/contact_obj.f90

xpdb2seams: pdb2seams.f90 seams_main.o seams_graph.o seams_sequences.o seams_sequence_generic.o seams_pdbtools.o seams_utils.o seams_debug.o
	$(FC) $(FCFLAGS) -I$(OBJ) -J$(OBJ) -o $@ $(LDFLAGS) pdb2seams.f90 $(OBJ)/seams_main.o $(OBJ)/seams_graph.o \
                               $(OBJ)/seams_sequences.o $(OBJ)/seams_sequence_generic.o $(OBJ)/seams_pdbtools.o \
							   $(OBJ)/seams_utils.o $(OBJ)/seams_debug.o

seams_main.o: $(SRC)/seams_main.f90 seams_sequences.o seams_pdbtools.o seams_graph.o seams_utils.o seams_debug.o
seams_sequences.o: $(SRC)/seams_sequences.f90 seams_sequence_generic.o
seams_graph.o: $(SRC)/seams_graph.f90 seams_sequences.o
seams_sequence_generic.o: $(SRC)/seams_sequence_generic.f90
seams_pdbtools.o: $(SRC)/seams_pdbtools.f90 seams_sequences.o seams_utils.o
seams_utils.o: $(SRC)/seams_utils.f90
seams_debug.o: $(SRC)/seams_debug.f90

# ======================================================================
# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -I$(OBJ) -J$(OBJ) -o $(OBJ)$@ $^ $(LDFLAGS) $<

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: $(SRC)/%.f90
	$(FC) $(FCFLAGS) -I$(OBJ) -J$(OBJ) -o $(OBJ)/$@ -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -fv *.o *.mod *.MOD $(SRC)/*.mod $(SRC)/*.o $(PROGRAMS) $(OBJ)/*.mod $(OBJ)/*.o

veryclean: clean
	rm -fv *~ $(PROGRAMS)
