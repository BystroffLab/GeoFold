## MASKER MAKE FILE
##
## modified: Tue Apr 15 20:04:49 EDT 2008
## C. Bystroff  (bystrc@rpi.edu)
##
## This is the makefile for compiling MASKER
## a progam for numerical estimation of the 
## solvent excluded surface  area.
## Currently, it needs fortran90 compiler 'pgf90'
## (available from Portland Group. www.pgroup.com )
## or the gnu fortran compiler g95 (www.g95.org)
##
## Use the INSTALL script to run this Makefile and to
## run the programs necessary to generate the binary masks.
##
## 'make all' compiles the masker module and programs
## 'make prof' compiles everything with the profiler turned on (much slower). Output: xpdbmask_p
## 'make setNNN' modifies the source of masker.f90 to generate masks of NNN bit where NNN={128, 256, 512, 1024, 2048 or 4096}
## 'make clean' deletes objects and executables
## 'make package' creates a tar package like the one you used to install this.
##
## CB Thu Apr 12 13:49:20 EDT 2007
## =============== NOTES ===================
## An insoluble bug in triangles.f90 force me to shut off  rendering
## routines for raster3D. These are routines starting with "render" in
## masker.f90.   Look for "BUG in " lines in masker.f90
##    Sun Jul 25 17:31:18 EDT 2010
##

F90prof = pgf90 -tpp6 -Mextend  -g -Mprof=lines
##F90 = pgf90 -Mextend   -tpp6 -fast -Bstatic 
# F90 = g95 -ftrace=full
F90 = gfortran
## F90 = /home/bystrc/src/g95-install/bin/i686-suse-linux-gnu-g95 -O2

MODS =  vectormath.o masker.o 
OBJS =  ${MODS} pdbmask.o 
UBJS =  ${MODS} contactmask.o 

all : xcontactmask xsas2ps xplottheta xviewmask xmakemask xpdbmask xbinarymask \
  xsortpdbmask xcountbits xmakecollars xr3drotate xvoidmask
 
set128 : masker.f90 INSTALL
	cp masker.f90 masker.f90.bck; \
	sed -e "s/MAXATOM=..../MAXATOM=128 /" masker.f90.bck > masker.f90; \
	sed -e "s/setenv MASKSIZE 4096/setenv MASKSIZE  128/" INSTALL > install.csh

set256 : masker.f90 INSTALL
	cp masker.f90 masker.f90.bck; \
	sed -e "s/MAXATOM=..../MAXATOM=256 /" masker.f90.bck > masker.f90; \
	sed -e "s/setenv MASKSIZE 4096/setenv MASKSIZE  256/" INSTALL > install.csh

set512 : masker.f90 INSTALL
	cp masker.f90 masker.f90.bck; \
	sed -e "s/MAXATOM=..../MAXATOM=512 /" masker.f90.bck > masker.f90; \
	sed -e "s/setenv MASKSIZE 4096/setenv MASKSIZE  512/" INSTALL > install.csh

set1024 : masker.f90 INSTALL
	cp masker.f90 masker.f90.bck; \
	sed -e "s/MAXATOM=..../MAXATOM=1024/" masker.f90.bck > masker.f90; \
	sed -e "s/setenv MASKSIZE 4096/setenv MASKSIZE  1024/" INSTALL > install.csh

set2048 : masker.f90 INSTALL
	cp masker.f90 masker.f90.bck; \
	sed -e "s/MAXATOM=..../MAXATOM=2048/" masker.f90.bck > masker.f90; \
	sed -e "s/setenv MASKSIZE 4096/setenv MASKSIZE  2048/" INSTALL > install.csh

set4096 : masker.f90 INSTALL
	cp masker.f90 masker.f90.bck; \
	sed -e "s/MAXATOM=..../MAXATOM=4096/" masker.f90.bck > masker.f90; \
	sed -e "s/setenv MASKSIZE 4096/setenv MASKSIZE  4096/" INSTALL > install.csh

masker.o : masker.f90 vectormath.o
	$(F90) -c masker.f90

vectormath.o : vectormath.f90
	$(F90) -c vectormath.f90

pdbmask.o : pdbmask.f90
	$(F90) -c pdbmask.f90

mcmask2.o : mcmask2.f90 
	$(F90) -c mcmask2.f90 

contactmask.o : contactmask.f90
	$(F90) -c contactmask.f90

voidmask.o : voidmask.f90
	$(F90) -c voidmask.f90

xvoidmask : voidmask.o ${MODS}
	$(F90) -o xvoidmask voidmask.o $(MODS)

viewmask.o : viewmask.f90 
	$(F90) -c viewmask.f90 

xviewmask : viewmask.o ${MODS}
	$(F90) -o xviewmask viewmask.o $(MODS)

xcontactmask : ${UBJS}
	$(F90) -o xcontactmask $(UBJS)

xpdbmask : ${OBJS}
	$(F90) -o xpdbmask $(OBJS)

xmcmask : mcmask2.o masker.o vectormath.o
	$(F90) -o xmcmask masker.o vectormath.o mcmask2.o 

xsas2ps : sas2ps.f90
	$(F90) -o xsas2ps sas2ps.f90

xmakemask : makemask.f90 $(MODS)
	$(F90) -o xmakemask makemask.f90 $(MODS)

xbinarymask : binarymask.f90 $(MODS)
	$(F90) -o xbinarymask binarymask.f90 $(MODS)

xsortpdbmask : sortpdbmask.f90 $(MODS)
	$(F90) -o xsortpdbmask sortpdbmask.f90 $(MODS)

xcollartest : collartest.f90 $(MODS)
	$(F90) -o xcollartest collartest.f90 $(MODS)

xplottheta : plottheta.f90 $(MODS)
	$(F90) -o xplottheta plottheta.f90 $(MODS)

xtriangles : triangles.f90 $(MODS)
	$(F90) -o xtriangles triangles.f90 $(MODS)

xmakecollars : makecollars.f90 $(MODS)
	$(F90) -o xmakecollars makecollars.f90 $(MODS)

xshowmask : showmask.f90 $(MODS)
	$(F90) -o xshowmask showmask.f90 $(MODS)

xcountbits : countbits.f90 
	$(F90) -o xcountbits countbits.f90

prof : ${OBJS}
	$(F90prof) -o xpdbmask_p $(OBJS)

xr3drotate : r3drotate.f90 vectormath.o
	$(F90) -o xr3drotate r3drotate.f90 vectormath.o

clean:
	-rm -f *.mod *.o x* install.csh slopes*

package: 
	pushd ..; \
	tar -zcvf masker.tgz masker/*.f90 masker/Makefile masker/INSTALL \
                  masker/atoms.lib masker/setup.csh \
                  masker/128.pdb masker/256.pdb masker/512.pdb masker/1024.pdb \
                  masker/2048.pdb masker/4096.pdb masker/slopes* \
                  masker/cbit* masker/README masker/test.pdb \
                  masker/gromacs_ff.prm masker/test.results masker/voidmask \
                  masker/pdbmask masker/contactmask

install.csh:
	echo "Please run : make setNNN, where NNN = 128, 256, 512, 1024, 2048 or 4096" 

install: install.csh
	csh install.csh
