#| ======================== GEOFOLD =================================
#| GEOFOLD is a suite of programs for simulating protein
#| unfolding pathways and kinetics. This installation package
#| is designed for Linux systems only. Other systems
#| are surrently unsupported.
#|-------------------------------------------------------------------
#| To install GeoFold programs, type:
#|
#| make clean
#| make all
#|
#| To test, use
#| make test
#|
#| To install MASKER programs only, use
#| make masker
#|
#| To run:
#|  ./RUNGEOFOLD.csh PDBFILE chaincode
#|
#| FORTRAN compiler:
#|  set FF to the fortran compiler on your system.
#|  For best results use the g95 compiler (www.g95.org)
#|
#| C.Bystroff  14 NOV 2008
#| www.bioinfo.rpi.edu/bystrc/
#| plase cite: Ramakrishnan V, Salem SM, Zaki MJ , Matthews SJ,
#| Srinivasan S, Colon W, & Bystroff C. (2009) Developing
#| a detailed mechanistic model for protein unfolding.
#| Ramakrishnan V, Srinivasan S, Salaem SM, Zaki MJ , Matthews SJ,
#| Colon W, & Bystroff C. (2012) GeoFold: Topology-based protein
#| unfolding pathways capture the effects of engineered disulfides
#| on kinetic stability. Proteins 80(3):920-934.
#|-------------------------------------------------------------------
#|    MODIFY  these settings to fit your system         --------------
#| Modification:
#|     Jun/05/2013: Changed the Fortran compiler to "gfortran" (before g95)
#|     Jun/05/2013: Added seams rule and geofold_seams module
#|
##---- ImageMagick convert
# CONVERT  = /ext2/www/html/applications/hybrid/bin/convert
CONVERT  = /usr/bin/convert

# The compiler
FF = gfortran -g -fno-range-check
#FF = gfortran -O3
FFdebug = gfortran -g -fbacktrace
#FCFLAGS = -g -fbacktrace
##---- FORTRAN 90 compiler (preferably g95)
# FF = g95 -O2
# FF = /home/bystrc/src/g95-install/bin/i686-suse-linux-gnu-g95 -O2
##---- C++ compiler
CPP = g++
##|---- don't change these ----------
TESTPDB = 1cis
TESTCHN = _
SCRIPT = RUNGEOFOLD.csh

help :
	more INSTALL

bin32 :
	make all

bin64 :
	make maxTraffic

all  : xgeofold xvoidmask xunfoldsim xcontactmask xgetchain xrenumber_one xpdb2cij \
       xageplot xpathway2ps x3to1 xpdb2hb xfit_poly maxTraffic xsplitseams seams/xpdb2seams \
       hxpathway2ps

debug : geofold_pivots.f90 geofold_global.f90 geofold.f90 vectormath.f90 geofold_masker.f90 geofold_hbonds.f90
	$(FFdebug)  -c vectormath.f90
	$(FFdebug)  -c geofold_global.f90
	$(FFdebug)  -c geofold_hbonds.f90
	$(FFdebug)  -c geofold_masker.f90
	$(FFdebug)  -c geofold_pivots.f90
	$(FFdebug)  -o xgeofold  geofold.f90 geofold_pivots.o geofold_global.o geofold_hbonds.o geofold_masker.o vectormath.o

xgeofold_split2 : geofold_pivots.o
	sed -e "s/MAXSPLIT=.*/MAXSPLIT=2/" geofold_global.f90 > geofold_global.f90.2
	mv geofold_global.f90.2 geofold_global.f90
	make xgeofold
	cp xgeofold xgeofold_split2

xgeofold_split4 : geofold_pivots.o
	sed -e "s/MAXSPLIT=.*/MAXSPLIT=4/" geofold_global.f90 > geofold_global.f90.2
	mv geofold_global.f90.2 geofold_global.f90
	make xgeofold
	cp xgeofold xgeofold_split4

xgeofold_split8 : geofold_pivots.o
	sed -e "s/MAXSPLIT=.*/MAXSPLIT=8/" geofold_global.f90 > geofold_global.f90.2
	mv geofold_global.f90.2 geofold_global.f90
	make xgeofold
	cp xgeofold xgeofold_split8

xgeofold_split16 : geofold_pivots.o
	sed -e "s/MAXSPLIT=.*/MAXSPLIT=16/" geofold_global.f90 > geofold_global.f90.2
	mv geofold_global.f90.2 geofold_global.f90
	make xgeofold
	cp xgeofold xgeofold_split16

xgeofold_split32 : geofold_pivots.o
	sed -e "s/MAXSPLIT=.*/MAXSPLIT=32/" geofold_global.f90 > geofold_global.f90.2
	mv geofold_global.f90.2 geofold_global.f90
	make xgeofold
	cp xgeofold xgeofold_split32

seams/xpdb2seams :  seams/pdb2seams.f90
	cd seams; make all -f Makefile.seams; cd ../

masker : masker.o

masker.o : masker/masker.o
	cp masker/masker.o .
	cp masker/masker.mod .

masker/masker.o : vectormath.o
	cd masker; make clean; cd ../
	cd masker; make set512; cd ../
	cd masker; make all; cd ../
	cd masker; make install; cd ../

voidmask.o : voidmask.f90 masker.o vectormath.o
	$(FF) -c voidmask.f90

xageplot : ageplot2ps.f90
	$(FF) -o xageplot ageplot2ps.f90

xpathway2ps : pathway2ps.f90
	$(FF) -o xpathway2ps pathway2ps.f90

hxpathway2ps : hxpathway2ps.f90
	$(FF) -o hxpathway2ps hxpathway2ps.f90

x3to1 : 3to1.f90
	$(FF) -o x3to1 3to1.f90

xgetchain : getchain.f90
	$(FF) -o xgetchain getchain.f90

xrenumber_one : renumber_one.f90
	$(FF) -o xrenumber_one renumber_one.f90

xpdb2cij : pdb2cij.f90
	$(FF) -o xpdb2cij pdb2cij.f90

xvoidmask : masker/xvoidmask
	cp masker/xvoidmask .

xcontactmask : masker/xcontactmask
	cp masker/xcontactmask .

xfit_poly : fit_poly.f90 lsq.o
	$(FF) -o xfit_poly fit_poly.f90 lsq.o

lsq.o : lsq.f90
	$(FF) -c lsq.f90

masker/xvoidmask : masker/masker.o
	cd masker; make xvoidmask; cd ../

masker/xcontactmask : masker/masker.o
	cd masker; make xcontactmask; cd ../

unfoldsim.o : unfoldsim.f90
	$(FF)  -fno-range-check -c unfoldsim.f90

pdb2hb.o : pdb2hb.f90
	$(FF) -c pdb2hb.f90

xunfoldsim : unfoldsim.o geofold_global.o geofold_pivots.o geofold_seams.o vectormath.o
	$(FF)  -fno-range-check -o xunfoldsim unfoldsim.o geofold_global.o geofold_pivots.o geofold_seams.o vectormath.o

xpdb2hb : pdb2hb.o
	$(FF) -o xpdb2hb pdb2hb.o

vectormath.o : vectormath.f90
	$(FF) -c vectormath.f90

geofold_global.o : geofold_global.f90 vectormath.o
	$(FF) -c geofold_global.f90

geofold_hbonds.o : geofold_hbonds.f90 vectormath.o geofold_pivots.o
	$(FF) -c geofold_hbonds.f90

geofold_masker.o : geofold_masker.f90 geofold_seams.o vectormath.o
	$(FF) -c geofold_masker.f90

geofold_seams.o: geofold_seams.f90 geofold_global.o
	$(FF) -c geofold_seams.f90

geofold_pivots.o : geofold_pivots.f90 geofold_global.o vectormath.o geofold_seams.o
	$(FF) -c geofold_pivots.f90

geofold_flory.o : geofold_flory.f90 geofold_global.o
	$(FF) -c geofold_flory.f90

xgeofold : geofold.f90 geofold_global.o geofold_flory.o geofold_pivots.o geofold_hbonds.o geofold_seams.o geofold_masker.o vectormath.o vb.incl
	$(FF) -o xgeofold geofold.f90 geofold_global.o geofold_flory.o geofold_pivots.o geofold_hbonds.o geofold_masker.o vectormath.o \
	      geofold_seams.o

xsplitseams: splitseams.o geofold_masker.o geofold_hbonds.o geofold_global.o geofold_pivots.o geofold_seams.o vectormath.o
	$(FF) -o xsplitseams splitseams.o geofold_seams.o geofold_masker.o geofold_global.o geofold_hbonds.o geofold_pivots.o vectormath.o

splitseams.o: splitseams.f90
	$(FF)  -c splitseams.f90

xmakemask : masker.o
	cp masker/xmakemask .

maxTraffic : maxTraffic.cpp graph.h token.h
	$(CPP) -o maxTraffic maxTraffic.cpp

## TEST function. Run GeoFold programs starting from just a PDB file.
## Last step (convert) requires ImageMagick http://www.imagemagick.org/

test : xgeofold xvoidmask xunfoldsim x3to1 xcontactmask xageplot
	cp parameters  $(TESTPDB)$(TESTCHN).par
	./xgetchain $(TESTCHN) < $(TESTPDB).pdb > $(TESTPDB)$(TESTCHN).pdb
	./xrenumber_one $(TESTPDB)$(TESTCHN).pdb junk.pdb ; mv junk.pdb $(TESTPDB)$(TESTCHN).pdb
	./x3to1 "?" < $(TESTPDB)$(TESTCHN).pdb > $(TESTPDB)$(TESTCHN).seq
	./xpdb2cij $(TESTPDB)$(TESTCHN).pdb 8. > $(TESTPDB)$(TESTCHN).cij
	./xpdb2hb $(TESTPDB)$(TESTCHN).par $(TESTPDB)$(TESTCHN).pdb  > $(TESTPDB)$(TESTCHN).hb
	echo "HBONDS $(TESTPDB)$(TESTCHN).hb" >> $(TESTPDB)$(TESTCHN).par
	source masker/masker_setup.sh; masker/xcontactmask $(TESTPDB)$(TESTCHN).pdb $(TESTPDB)$(TESTCHN).sas 1.4
	echo "CONTACTS $(TESTPDB)$(TESTCHN).sas" >> $(TESTPDB)$(TESTCHN).par
	source masker/masker_setup.sh; masker/xvoidmask $(TESTPDB)$(TESTCHN).pdb $(TESTPDB)$(TESTCHN).void 1.4 1.2 1.4
	./xgeofold $(TESTPDB)$(TESTCHN).void $(TESTPDB)$(TESTCHN).dag $(TESTPDB)$(TESTCHN).par
	./xunfoldsim $(TESTPDB)$(TESTCHN).dag $(TESTPDB)$(TESTCHN).par
	./xpathway2ps $(TESTPDB)$(TESTCHN).seq $(TESTPDB)$(TESTCHN).dag.path $(TESTPDB)$(TESTCHN).cij $(TESTPDB)$(TESTCHN).ps 4
	$(CONVERT) -geometry 600x600 $(TESTPDB)$(TESTCHN).ps $(TESTPDB)$(TESTCHN).jpg

package : masker.tgz
	cd ..; \
	tar -zcvf geofold/geofold.tgz geofold/geofold.f90 geofold/geofold_pivots.f90 \
                 geofold/vectormath.f90 geofold/2ptl_.sas  \
                 geofold/Makefile geofold/1dwk+.pdb geofold/1dwk+.sas \
                 geofold/setupgeofold.csh geofold/geofold_global.f90 \
                 geofold/geofold.f90.bck \
                 geofold/geofold_global.f90.bck geofold/1dwkEG.pdb \
                 geofold/1dwkEG.sas geofold/1cis.pdb  \
                 geofold/vb.incl \
                 geofold/vectorball.csh  geofold/masker.tgz \
                 geofold/unfoldsim.f90 geofold/pdb2cij.f90 geofold/pdb2hb.f90 \
                 geofold/ageplot2ps.f90 geofold/renumber_one.f90 \
		 geofold/example geofold/getchain.f90 geofold/$(TESTPDB).pdb \
                 geofold/3to1.f90 geofold/pathway2ps.f90 geofold/melting.csh \
		 geofold/unfold.csh geofold/rungeofold.csh geofold/parameters \
		 geofold/rungeofold.csh geofold/ku.csh geofold/fit_poly.f90 \
		 geofold/lsq.f90 geofold/maxTraffic.cpp geofold/graph.h geofold/token.h \
		 geofold/RUNGEOFOLD.csh geofold/creategnuplot.csh \
		 geofold/header.html geofold/header_refresh.html \
		 geofold/unfoldsim2.f90 geofold/splitseams.f90   \
		 geofold/seams/*.f90

clean :
	rm -vrf *.o xgeofold x3to1 xgeofold xvoidmask xunfoldsim xcontactmask \
	xgetchain xrenumber_one xpdb2cij xageplot xpathway2ps x3to1 xpdb2hb \
	xfit_poly maxTraffic xsplitseams  *.mod *.o *.dSYM

cleaner :
	rm -vf x* *.o *.mod masker.tgz masker/*.o masker/*.mod masker/x*

## UTILITY tests. Modify as needed.

test2 : xgeofold
	./xgeofold 1dwkEG.void 1dwkEG.sas 1dwkEG.dag 0.4 0.05 0.05

test3 : xgeofold
	./xgeofold 1dwk+.void 1dwk+.sas 1dwk+.dag

test4 : xgeofold
	./xgeofold 1ten.void 1ten.sas 1ten.dag

test5 : xunfoldsim
	./xunfoldsim 1ten.dag

test6 : xunfoldsim
	./xunfoldsim 1csp_.dag

setup :
	./setupgeofold.csh

checkout :
	cp geofold.f90.bck geofold.f90
	cp geofold_global.f90.bck geofold_global.f90
	cp geofold_pivots.f90.bck geofold_pivots.f90
	cp unfoldsim.f90.bck unfoldsim.f90

checkin :
	cp geofold.f90 geofold.f90.bck
	cp geofold_pivots.f90 geofold_pivots.f90.bck
	cp geofold_global.f90 geofold_global.f90.bck
	cp unfoldsim.f90 unfoldsim.f90.bck
