!*=========================================================================
!** Copyright 1997-2002 Chris Bystroff
!**
!** Permission is hereby granted to use, execute, copy, distribute, modify,
!** and distribute in modified form this software for non-profit purposes,
!** provided this notice is retained in any and all copies.
!**
!** For for-profit usage, please contact bystrc@rpi.edu
!*=========================================================================
module ARCHIVERO
  !*=========================================================================
  !** These are the initial coordinates for the 20 amino acids.
  !** All PHI, PSI, and OMG angles are set to 180.00 degrees.
  !** (even for PRO)
  !** J is the number of amino acid 
  !**           in alphabetical order of one-letter codes. 
  !** jchi is the number of sidechain dihedrals
  !** jside is the number of atoms in the side chain, numbered
  !**           from 5. [ Atoms 1,2,3,4 are N,CA,C,O ]
  !** The numbers of the atoms (J1,J2) that define
  !**           the chi rotation are 2,5 for chi1, 5,6 for chi2
  !**           if it exists, and 6,7 for chi3 if it exists.
  !**           atoms from J2+1 to the end of
  !**           the residue can move relative to atoms after that and
  !**           before J1-1 [ unless J1 = 2 (C-alpha)
  !**           in which case atoms before 5, except 2, are movable.]
  !** names(14) are the atom names, not counting the last three
  !**           'dummy' atoms
  !** base(20)  are the residue names
  !** dang(8)   are the values for phi,psi,omg,and chi(1-5)
  !**           as they are in the database coordinates.
  !** xyz(3,17) are the coordinates (converted to real*8) of the residue
  !**           atoms including three 'dummy' atoms to serve as
  !**           a template for the next addition.
  !** charge() is the atomic charge
  !** htype() is 1=no h-bonding, 2=donor only, 3=receptor only, 4=either
  !** iorder = how many atoms covalently bonded (1,2,or 3)
  !** iat    = 0=hydrogen,1=carbon,2=nitrogen,3=oxygen,4=sulfur
  !** CFFTtors(MAXFFT,I)= 1st MAXFFT (i.e.32) terms of the fourier transform of 
  !**          a 128-bin torsion angle frequency profile. 
  !** tposit(20) = position (I) of AA=1,20 in array CFFTtors(), 2nd index.
  !**
  !CB Changed 11-DEC-94 to add N-hydrogen. It will be numbered after
  !CB the sidechain (nside(I)+4). The charges are set for the no-hydrogen
  !CB option. If hydrogens are used (USEH=.TRUE.) then the charges
  !CB will be reset for N and H. These are still dummy values until I get
  !CB better info! --CB
  !CB
  !CB Modified to mark un-read template atoms with x=y=z='999.0' 21-apr-95
  !CB
  !CB Added 'bfby' data lines. 26-apr-95
  !*****
        use prot_globals
        private
        !** atomic weight of atom(iat)
        data mw     /1.,12.,13.,14.,15.,15.,17.,16.,17.,32./
        data vdw1   /VDWSIZ*0.0/
        data vdw0   /VDWSIZ*0.0/
        data hydr1  /VDWSIZ*0.0/
        data hydr0  /VDWSIZ*0.0/
        data h_bond /MAXRES*0/
        data OHnrg  /MAXRES*0.0/
        data xvalue / 1.27, 1.32, 1.21, 0.92, 2.35, 1.82, 2.51, &
                      1.48, 0.67, 2.95, 0.73, 3.00, 2.29, 2.98, 2.94/
        data yvalue /15*5.0/
        data zvalue /15*0.0/
        data wvalue /15*0.0/
        data gvalue /15*0.0/
        data xvold  / 1.27, 1.32, 1.21, 0.92, 2.35, 1.82, 2.51, &
                      1.48, 0.67, 2.95, 0.73, 3.00, 2.29, 2.98, 2.94/
        data yvold  /15*5.0/
        data zvold  /15*0.0/
        data wvold  /15*0.0/
        data gvold  /15*0.0/
        data hrval  /15*1000.0/
        data torrstrwt /MAXANG*0.0/
!** VDW well depths (but not radii) from Levitt et al., Comp Phys. Comm. (1994)
        data vdwwell   /0.03763,0.07382,0.07382,0.07382,0.41315,0.41315, &
                        0.18479,0.18479,0.07382/
!** subsets of 'CFFTtors':
!** These arrays contain the fourier coefficients for calculating
!** the torsion angle energy and force. MAXFFT coeffs / angle
!** these are the positions in the CFFTtors array (note: ARG has only 4 chi)
!**          nchi    /0,1,2, 3, 2, 0, 2, 2, 4, 2, 3, 2, 0, 3, 4, 1, 1, 1, 2, 2/
!**            AA     A C D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
!** # of dihedrals:   2 3 4  5  4  2  4  4  6  4  5  4  1  5  6  3  3  3  4  4
!**           tposit /1,3,6,10,15,19,21,25,29,35,39,44,48,49,54,60,63,66,70,74/
!** complex*8 CFFTtors(32,78) = 20kB of memory
        data names(1:siz( 1), 1)  /'N   ','CA  ','C   ','O   ','CB  '/                      ! A N,CA,CB,C,O                       5
        data names(1:siz( 2), 2)  /'N   ','CA  ','C   ','O   ','CB  ','SG  '/               ! C N,CA,C,O,CB,SG                    6
        data names(1:siz( 3), 3)  /'N   ','CA  ','C   ','O   ','CB  ','CG  ','OD1 ',  &     ! D N,CA,C,O,CB,CG,OD1,OD2            8
                                   'OD2 '/
        data names(1:siz( 4), 4)  /'N   ','CA  ','C   ','O   ','CB  ','CG  ','CD  ',  &     ! E N,CA,C,O,CB,CG,CD,OE1,OE2         9
                                   'OE1 ','OE2 '/
        data names(1:siz( 5), 5)  /'N   ','CA  ','C   ','O   ','CB  ','CG  ','CD1 ',  &     ! F N,CA,C,O,CB,CG,CD1,CD2,CE1,CE2,CZ 11
                                   'CD2 ','CE1 ','CE2 ','CZ  '/
        data names(1:siz( 6), 6)  /'N   ','CA  ','C   ','O   '/                             ! G N,CA,C,O                          4
        data names(1:siz( 7), 7)  /'N   ','CA  ','C   ','O   ','CB  ','CG  ','ND1 ',  &     ! H N,CA,C,O,CB,CG,ND1,CD2,CE1,NE2    10
                                   'CD2 ','CE1 ','NE2 '/
        data names(1:siz( 8), 8)  /'N   ','CA  ','C   ','O   ','CB  ','CG1 ','CG2 ','CD1 '/ ! I N,CA,C,O,CB,CG1,CG2,CD1           8
        data names(1:siz( 9), 9)  /'N   ','CA  ','C   ','O   ','CB  ','CG  ','CD  ',  &     ! K N,CA,C,O,CB,CG,CD,CE,NZ           9
                                   'CE  ','NZ  '/
        data names(1:siz(10),10)  /'N   ','CA  ','C   ','O   ','CB  ','CG  ','CD1 ','CD2 '/ ! L N,CA,C,O,CB,CG,CD1,CD2            8
        data names(1:siz(11),11)  /'N   ','CA  ','C   ','O   ','CB  ','CG  ','SD  ','CE  '/ ! M N,CA,C,O,CB,CG,SD,CE              8
        data names(1:siz(12),12)  /'N   ','CA  ','C   ','O   ','CB  ','CG  ','OD1 ','ND2 '/ ! N N,CA,C,O,CB,CG,OD1,ND2            8
        data names(1:siz(13),13)  /'N   ','CA  ','C   ','O   ','CB  ','CG  ','CD  '/        ! P N,CA,C,O,CB,CG,CD                 7
        data names(1:siz(14),14)  /'N   ','CA  ','C   ','O   ','CB  ','CG  ','CD  ',  &     ! Q N,CA,C,O,CB,CG,CD,OE1,NE2         9
                                   'OE1 ','NE2 '/    
        data names(1:siz(15),15)  /'N   ','CA  ','C   ','O   ','CB  ','CG  ','CD  ',  &     ! R N,CA,C,O,CB,CG,CD,NE,CZ,NH1,NH2   11
                                   'NE  ','CZ  ','NH1 ','NH2 '/
        data names(1:siz(16),16)  /'N   ','CA  ','C   ','O   ','CB  ','OG  '/               ! S N,CA,C,O,CB,OG                    6
        data names(1:siz(17),17)  /'N   ','CA  ','C   ','O   ','CB  ','OG1 ','CG2 '/        ! T N,CA,C,O,CB,OG1,CG2               7
        data names(1:siz(18),18)  /'N   ','CA  ','C   ','O   ','CB  ','CG1 ','CG2 '/        ! V N,CA,C,O,CB,CG1,CG2               7
        data names(1:siz(19),19)  /'N   ','CA  ','C   ','O   ','CB  ','CG  ','CD1 ',  &     ! W N,CA,C,O,CB,CG,CD1,CD2,NE1,       14
                                   'CD2 ','NE1 ','CE2 ','CE3 ','CZ2 ','CZ3 ','CH2 '/        !            CE2,CE3,CZ2,CZ3,CH2
        data names(1:siz(20),20)  /'N   ','CA  ','C   ','O   ','CB  ','CG  ','CD1 ',  &     ! Y N,CA,C,O,CB,CG,CD1,CD2,CE1,       12
                                   'CD2 ','CE1 ','CE2 ','CZ  ','OH  '/                      !            CE2,CZ,OH
        data names(1:1 ,21)  /'OXT '/
!cb These arrays were incremented by one 11-DEC-94 to include the N-hydrogen.
!cb Except PRO & OXT (Pro has no H on N, except when it's residue 1, but
!cb in that case the H is not used in this program.)
!******  radii for H C N O S => radii(0:4)
!**      data radii / 1.20, 2.00, 1.50, 1.40, 1.85 /
!**      data radii / 0.55, 1.55, 1.50, 1.35, 1.80 /
!**        include 'prot_wzsvdw.inc'
!** 17 WZS group types COODE,CNK,CNR,CNH,COST,CONQN,CP,Calpha,
!**                    COY,CSM,CON,CNW,CY,Caliph,CW,CF,CSC
!**                    in rough order of hydrophilicity.
!** Inluded are VDW radi for SAS calculation. 
        data wzsvdwrad /  &
         3.45, 3.35, 3.40, 3.25, 3.00, 3.45, 2.55, 3.05, 3.30, &
         3.35, 3.55, 3.40, 3.30, 3.15, 3.55, 2.95, 3.05, 3.15, &
         3.40, 3.40, 3.35, 3.10, 3.10, 3.00, 2.75, 2.90, 2.65, &
         3.25, 3.30, 3.10, 2.95, 3.00, 3.50, 2.50, 2.85, 3.15, &
         3.00, 3.15, 3.10, 3.00, 3.05, 3.25, 2.40, 2.60, 3.10, &
         3.45, 3.55, 3.00, 3.50, 3.25, 5.00, 2.65, 3.25, 4.00, &
         2.55, 2.95, 2.75, 2.50, 2.40, 2.65, 2.60, 2.40, 2.85, &
         3.05, 3.05, 2.90, 2.85, 2.60, 3.25, 2.40, 2.70, 3.60, &
         3.30, 3.15, 2.65, 3.15, 3.10, 4.00, 2.85, 3.60, 1.65/
        data wzssolrad /1.80,1.85,1.925,2.00,1.75,1.85,1.60,1.65,2.00/
        data wzssig /-0.010981,-0.012290,-0.007259,-0.003320,-0.000936,  &
                      0.001399, 0.004994, 0.002688, 0.003315, 0.009280,  &
                      0.003466, 0.006145, 0.006439, 0.009961, 0.012050,  &
                      0.017880, 0.028740/
        data wzsgrp  /'COO(DE) ','CN(K)   ','CN(R)   ','CN(H)   ','CO(ST)  ',  &
                      'CON(QN) ','C(P)    ','Calpha  ','CO(Y)   ','CS(M)   ',  &
                      'CON     ','CN(W)   ','C(Y)    ','Caliph  ','C(W)    ',  &
                      'C(F)    ','CS(C)   '/
!** Atom type numbers based on WZS
!** 9  WZS atom types:     C,   CH,  CH2, CH3, NH,  NH3,  O,   OH,  S
!** iat=0 is for hydrogens.
        data iat(1:siz( 1)+1, 1)   /5,2,1,7,4,0/                     ! A
        data iat(1:siz( 2)+1, 2)   /5,2,1,7,3,9,0/                   ! C
        data iat(1:siz( 3)+1, 3)   /5,2,1,7,3,1,7,7,0/               ! D
        data iat(1:siz( 4)+1, 4)   /5,2,1,7,3,3,1,7,7,0/             ! E
        data iat(1:siz( 5)+1, 5)   /5,2,1,7,3,1,2,2,2,2,2,0/         ! F
        data iat(1:siz( 6)+1, 6)   /5,2,1,7,0/                       ! G
        data iat(1:siz( 7)+1, 7)   /5,2,1,7,3,1,5,2,2,5,0/           ! H 
        data iat(1:siz( 8)+1, 8)   /5,2,1,7,2,4,3,4,0/               ! I 
        data iat(1:siz( 9)+1, 9)   /5,2,1,7,3,3,3,3,6,0/             ! K
        data iat(1:siz(10)+1,10)   /5,2,1,7,3,2,4,4,0/               ! L
        data iat(1:siz(11)+1,11)   /5,2,1,7,3,3,9,4,0/               ! M
        data iat(1:siz(12)+1,12)   /5,2,1,7,3,1,7,5,0/               ! N
        data iat(1:siz(13)  ,13)   /5,2,1,7,3,3,3/                   ! P
        data iat(1:siz(14)+1,14)   /5,2,1,7,3,3,1,7,5,0/             ! Q
        data iat(1:siz(15)+1,15)   /5,2,1,7,3,3,3,5,1,5,5,0/         ! R
        data iat(1:siz(16)+1,16)   /5,2,1,7,3,8,0/                   ! S
        data iat(1:siz(17)+1,17)   /5,2,1,7,2,8,4,0/                 ! T
        data iat(1:siz(18)+1,18)   /5,2,1,7,2,4,4,0/                 ! V
        data iat(1:siz(19)+1,19)   /5,2,1,7,3,1,2,1,5,1,2,2,2,2,0/   ! W
        data iat(1:siz(20)+1,20)   /5,2,1,7,3,1,2,2,2,2,1,8,0/       ! Y
        data iat(1, 21)   /7/                                        ! X
!** WZS group type by atom
        data wzs(1:siz( 1)+1, 1)   /11,8,11,11,14,0/
        data wzs(1:siz( 2)+1, 2)   /11,8,11,11,17,17,0/
        data wzs(1:siz( 3)+1, 3)   /11,8,11,11,14,1,1,1,0/
        data wzs(1:siz( 4)+1, 4)   /11,8,11,11,14,14,1,1,1,0/
        data wzs(1:siz( 5)+1, 5)   /11,8,11,11,14,16,16,16,16,16,16,0/
        data wzs(1:siz( 6)+1, 6)   /11,8,11,11,0/
        data wzs(1:siz( 7)+1, 7)   /11,8,11,11,14,4,4,4,4,4,0/
        data wzs(1:siz( 8)+1, 8)   /11,8,11,11,14,14,14,14,0/
        data wzs(1:siz( 9)+1, 9)   /11,8,11,11,14,14,14,2,2,0/
        data wzs(1:siz(10)+1,10)   /11,8,11,11,14,14,14,14,0/
        data wzs(1:siz(11)+1,11)   /11,8,11,11,14,10,10,10,0/
        data wzs(1:siz(12)+1,12)   /11,8,11,11,14,6,6,6,0/
        data wzs(1:siz(13)  ,13)   /11,8,11,11,7,7,7/
        data wzs(1:siz(14)+1,14)   /11,8,11,11,14,14,6,6,6,0/
        data wzs(1:siz(15)+1,15)   /11,8,11,11,14,14,14,3,3,3,3,0/
        data wzs(1:siz(16)+1,16)   /11,8,11,11,5,5,0/
        data wzs(1:siz(17)+1,17)   /11,8,11,11,5,5,14,0/
        data wzs(1:siz(18)+1,18)   /11,8,11,11,14,14,14,0/
        data wzs(1:siz(19)+1,19)   /11,8,11,11,14,12,12,15,12,15,15,15,15,15,0/
        data wzs(1:siz(20)+1,20)   /11,8,11,11,14,13,13,13,13,13,9,9,0/
        data wzs(1          ,21)   /1/
!** number of bonded neighbors (subtract 1 for N-term)
        data iorder(1:siz( 1), 1)  /2,3,3,1,1/
        data iorder(1:siz( 2), 2)  /2,3,3,1,2,1/
        data iorder(1:siz( 3), 3)  /2,3,3,1,2,3,1,1/
        data iorder(1:siz( 4), 4)  /2,3,3,1,2,2,3,1,1/
        data iorder(1:siz( 5), 5)  /2,3,3,1,2,3,2,2,2,2,2/
        data iorder(1:siz( 6), 6)  /2,2,3,1/
        data iorder(1:siz( 7), 7)  /2,3,3,1,2,3,2,2,2,2/
        data iorder(1:siz( 8), 8)  /2,3,3,1,3,2,1,1/
        data iorder(1:siz( 9), 9)  /2,3,3,1,2,2,2,2,1/
        data iorder(1:siz(10),10)  /2,3,3,1,2,3,1,1/
        data iorder(1:siz(11),11)  /2,3,3,1,2,2,2,1/
        data iorder(1:siz(12),12)  /2,3,3,1,2,3,1,1/
        data iorder(1:siz(13),13)  /3,3,3,1,2,2,2/
        data iorder(1:siz(14),14)  /2,3,3,1,2,2,3,1,1/
        data iorder(1:siz(15),15)  /2,3,3,1,2,2,2,2,3,1,1/
        data iorder(1:siz(16),16)  /2,3,3,1,2,1/
        data iorder(1:siz(17),17)  /2,3,3,1,3,1,1/
        data iorder(1:siz(18),18)  /2,3,3,1,3,1,1/
        data iorder(1:siz(19),19)  /2,3,3,1,2,3,2,3,2,3,2,2,2,2/
        data iorder(1:siz(20),20)  /2,3,3,1,2,3,2,2,2,2,3,1/
        data iorder(1,21)          /1/                          ! C-terminal OXT
!** number of 1-3 nieghbors (add one to C followed by PRO, subtract 2 for N-term
!** subtract 1 for C of C-term)
        data iorder2(1:siz( 1), 1)  /4,3,3,2,2/
        data iorder2(1:siz( 2), 2)  /4,4,3,2,2,1/
        data iorder2(1:siz( 3), 3)  /4,4,3,2,4,1,2,2/
        data iorder2(1:siz( 4), 4)  /4,4,3,2,3,3,1,2,2/
        data iorder2(1:siz( 5), 5)  /4,4,3,2,4,3,3,3,2,2,2/
        data iorder2(1:siz( 6), 6)  /3,3,2,2/
        data iorder2(1:siz( 7), 7)  /4,4,3,2,4,3,3,3,2,2/
        data iorder2(1:siz( 8), 8)  /4,5,3,2,3,2,2,1/
        data iorder2(1:siz( 9), 9)  /4,4,3,2,3,2,2,1,1/
        data iorder2(1:siz(10),10)  /4,4,3,2,4,1,2,2/
        data iorder2(1:siz(11),11)  /4,4,3,2,3,2,1,1/
        data iorder2(1:siz(12),12)  /4,4,3,2,4,1,2,2/
        data iorder2(1:siz(13),13)  /5,5,3,2,3,2,3/
        data iorder2(1:siz(14),14)  /4,4,3,2,3,3,1,2,2/
        data iorder2(1:siz(15),15)  /4,4,3,2,3,2,2,3,1,2,2/
        data iorder2(1:siz(16),16)  /4,4,3,2,2,1/
        data iorder2(1:siz(17),17)  /4,5,3,2,2,2,2/
        data iorder2(1:siz(18),18)  /4,5,3,2,2,2,2/
        data iorder2(1:siz(19),19)  /4,4,3,2,4,4,3,5,3,4,3,3,2,2/
        data iorder2(1:siz(20),20)  /4,4,3,2,4,3,3,3,3,3,2,2/
        data iorder2(1        ,21)  /2/
!** SASAtypes by atom. These are based on 'order' and 'iorder2'
        data sastyp(1:siz( 1), 1)  /7,12,12,2,2/
        data sastyp(1:siz( 2), 2)  /7,13,12,2,4,1/
        data sastyp(1:siz( 3), 3)  /7,13,12,2,7,8,2,2/
        data sastyp(1:siz( 4), 4)  /7,13,12,2,6,6,8,2,2/
        data sastyp(1:siz( 5), 5)  /7,13,12,2,7,11,6,6,5,5,5/
        data sastyp(1:siz( 6), 6)  /6,6,10,2/
        data sastyp(1:siz( 7), 7)  /7,13,12,2,7,13,6,6,4,4/
        data sastyp(1:siz( 8), 8)  /7,15,12,2,12,4,2,1/
        data sastyp(1:siz( 9), 9)  /7,13,12,2,6,5,5,3,1/
        data sastyp(1:siz(10),10)  /7,13,12,2,7,8,2,2/
        data sastyp(1:siz(11),11)  /7,13,12,2,6,5,3,1/
        data sastyp(1:siz(12),12)  /7,13,12,2,7,8,2,2/
        data sastyp(1:siz(13),13)  /15,15,12,2,6,5,6/
        data sastyp(1:siz(14),14)  /7,13,12,2,6,6,8,2,2/
        data sastyp(1:siz(15),15)  /7,13,12,2,6,5,5,6,8,2,2/
        data sastyp(1:siz(16),16)  /7,13,12,2,4,1/
        data sastyp(1:siz(17),17)  /7,15,12,2,9,2,2/
        data sastyp(1:siz(18),18)  /7,15,12,2,9,2,2/
        data sastyp(1:siz(19),19)  /7,13,12,2,7,13,6,15,6,15,6,6,5,5/
        data sastyp(1:siz(20),20)  /7,13,12,2,7,11,6,6,6,6,10,2/
        data sastyp(1, 1)          /2/
!****** residue 3-letter names
        data base /'ALA ','CYS ','ASP ','GLU ','PHE ','GLY ', &
                   'HIS ','ILE ','LYS ','LEU ','MET ','ASN ','PRO ', &
                   'GLN ','ARG ','SER ','THR ','VAL ','TRP ','TYR '/
!****** atom numbers defining chi-angle axes, in pairs.
        data chiposit /2,5,5,6,6,7,7,8,8,9/
        data angle /MAXANG*0.0/
        data tau   /MAXANG*0.0/
        data accel /MAXANG*0.0/
        data oldaccel /MAXANG*0.0/
        data veloc /MAXANG*0.0/
        data sahsa /MAXAT*0.0/
!** 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
!** A C D E F G H I K  L  M  N  P  Q  R  S  T  V  W  Y OXT
        data nside /1,2,4,5,7,0,6,4,5,4,4,4,3,5,7,2,3,3,10,8,-3,0/
        data nchi  /0,1,2,3,2,0,2,2,4,2,3,2,0,3,4,1,1,1,2, 2/
!****** sidechain rotamer angles. 1st 3. 
!** These are the defaults when angles when not read in or measured.
        data mer / 0., 0., 0.,     &  !  A
           -65.2,    0.,    0.,    &  !  C
           -68.3, -25.7,   0.0,    &  !  D
           -69.6,-177.2, -11.4,    &  !  E
           -66.3,  94.3,   0.0,    &  !  F
             0.0,   0.0,   0.0,    &  !  G
           -62.8, -74.3,   0.0,    &  !  H
           -60.9, 168.7,   0.0,    &  !  I
           -68.9,-178.4, 180.0,    &  !  K
           -64.9, 176.0,   0.0,    &  !  L
           -64.5, -68.5, -75.6,    &  !  M
           -68.3, -36.8,   0.0,    &  !  N
             0.0,   0.0,   0.0,    &  !  P
           -58.7, -63.8, -46.3,    &  !  Q
           -67.6, 176.9, 180.0,    &  !  R
            64.7,   0.0,   0.0,    &  !  S
            62.7,   0.0,   0.0,    &  !  T
           173.5,   0.0,   0.0,    &  !  V
           -70.4, 100.5,   0.0,    &  !  W
           -66.5,  96.6,   0.0/       !  Y
!** ANGS are the current dihedrals in 
        data ANGS(1:7,1)  /180.0,180.0,180.0 ,   0.0,   0.0,  0.0,  0.0/
        data ANGS(1:7,2)  /180.0,180.0,180.0 , 180.0,   0.0,  0.0,  0.0/
        data ANGS(1:7,3)  /180.0,180.0,180.0 ,-175.0,  33.0,  0.0,  0.0/
        data ANGS(1:7,4)  /180.0,180.0,180.0 ,  60.3,-174.5,  0.0,  0.0/
        data ANGS(1:7,5)  /180.0,180.0,180.0 , -57.0,  87.3,  0.0,  0.0/
        data ANGS(1:7,6)  /180.0,180.0,180.0 ,   0.0,   0.0,  0.0,  0.0/
        data ANGS(1:7,7)  /180.0,180.0,180.0 ,-120.0,   0.1,  0.0,  0.0/
        data ANGS(1:7,8)  /180.0,180.0,180.0 , 178.0, 173.0,  0.0,  0.0/
        data ANGS(1:7,9)  /180.0,180.0,180.0 , 180.0, 180.0,180.0, 180.0/
        data ANGS(1:7,10) /180.0,180.0,180.0 , -56.6,-174.3,  0.0,  0.0/
        data ANGS(1:7,11) /180.0,180.0,180.0 , -67.0,-165.0,-166.,  0.0/
        data ANGS(1:7,12) /180.0,180.0,180.0 , -54.1,  99.3,  0.0,  0.0/
        data ANGS(1:7,13) /-74.7,180.0,180.0 ,   0.0,   0.0,  0.0,  0.0/
        data ANGS(1:7,14) /180.0,180.0,180.0 , 180.0, 180.0,180.0,  0.0/
        data ANGS(1:7,15) /180.0,180.0,180.0 , 180.0, 180.0,180.0,  0.0/
        data ANGS(1:7,16) /180.0,180.0,180.0 ,  70.0,   0.0,  0.0,  0.0/
        data ANGS(1:7,17) /180.0,180.0,180.0 , -60.0,   0.0,  0.0,  0.0/
        data ANGS(1:7,18) /180.0,180.0,180.0 , -69.2,   0.0,  0.0,  0.0/
        data ANGS(1:7,19) /180.0,180.0,180.0 , -65.3,  98.7,  0.0,  0.0/
        data ANGS(1:7,20) /180.0,180.0,180.0 , -59.0, 101.0,  0.0,  0.0/
 !** Partial charges were derived from MOE for a peptide of sequence GGACDEFGHIKLMNPQRSTVWYGG using the MMFF force field
        data charge(1:siz( 1),  1)  /  -0.360,  0.361,  0.000,  0.569, -0.570  /                  !  A
        data charge(1:siz( 2),  2)  /  -0.360,  0.361,  0.569, -0.570,  0.230, -0.410  /          !  C
        data charge(1:siz( 3),  3)  /  -0.360,  0.361,  0.569, -0.570, -0.106,  0.906, -0.900,  &
                                       -0.900  /                                                  !  D
        data charge(1:siz( 4),  4)  /  -0.360,  0.361,  0.569, -0.570,  0.000, -0.106,  0.906,  &
                                       -0.900, -0.900  /                                          !  E
        data charge(1:siz( 5),  5)  /  -0.360,  0.361,  0.569, -0.570,  0.143, -0.143, -0.150,  &
                                       -0.150, -0.150, -0.150, -0.150  /                          !  F
        data charge(1:siz( 6),  6)  /  -0.360,  0.361,  0.569, -0.570  /                          !  G
        data charge(1:siz( 7),  7)  /  -0.360,  0.361,  0.569, -0.570,  0.181,  0.046, -0.565,  &
                                       -0.302, 0.037,  0.033  /                                   !  H
        data charge(1:siz( 8),  8)  /  -0.360,  0.361,  0.569, -0.570,  0.000,  0.000,  0.000,  &
                                        0.000  /                                                  !  I
        data charge(1:siz( 9),  9)  /  -0.360,  0.361,  0.569, -0.570,  0.000,  0.000,  0.000,  &
                                        0.503, -0.853  /                                          !  K
        data charge(1:siz(10), 10)  /  -0.360,  0.361,  0.569, -0.570,  0.000,  0.000,  0.000,  &
                                        0.000  /                                                  !  L
        data charge(1:siz(11), 11)  /  -0.360,  0.361,  0.569, -0.570,  0.000,  0.230, -0.460,  &
                                        0.230  /                                                  !  M
        data charge(1:siz(12), 12)  /  -0.360,  0.361,  0.569, -0.570,  0.061,  0.569, -0.570,  &
                                       -0.800  /                                                  !  N
        data charge(1:siz(13), 13)  /  -0.660,  0.361,  0.569, -0.570,  0.000,  0.000,  0.300  /  !  P
        data charge(1:siz(14), 14)  /  -0.360,  0.361,  0.569, -0.570,  0.000,  &
                                        0.061,  0.569, -0.570, -0.800  /                          !  Q
        data charge(1:siz(15), 15)  /  -0.360,  0.361,  0.569, -0.570,  0.000,  0.000,  &
                                        0.328, -0.844,  1.200, -0.967, -0.967  /                  !  R
        data charge(1:siz(16), 16)  /  -0.360,  0.361,  0.569, -0.570,  0.280, -0.680  /          !  S
        data charge(1:siz(17), 17)  /  -0.360,  0.361,  0.569, -0.570,  0.280, -0.680,  0.000  /  !  T
        data charge(1:siz(18), 18)  /  -0.360,  0.361,  0.569, -0.570,  0.000,  0.000,  0.000  /  !  V
        data charge(1:siz(19), 19)  /  -0.360,  0.361,  0.569, -0.570,  0.181, -0.181, -0.302,  &
                                        0.000,  0.033, -0.152, -0.150, -0.150, -0.150, -0.150  /  !  W
        data charge(1:siz(20), 20)  /  -0.360,  0.361,  0.569, -0.570,  0.143, -0.143, -0.150,  &
                                       -0.150, -0.150, -0.150,  0.083, -0.532   /                 !  Y
        data charge(1, 21)  /  -0.360 /
!** Hbond type by atom
        data htype(1:siz( 1)+1, 1)  /2,1,1,3,1,1/
        data htype(1:siz( 2)+1, 2)  /2,1,1,3,1,4,1/
        data htype(1:siz( 3)+1, 3)  /2,1,1,3,1,1,3,3,1/
        data htype(1:siz( 4)+1, 4)  /2,1,1,3,1,1,1,3,3,1/
        data htype(1:siz( 5)+1, 5)  /2,1,1,3,1,1,1,1,1,1,1,1/
        data htype(1:siz( 6)+1, 6)  /2,1,1,3,1/
        data htype(1:siz( 7)+1, 7)  /2,1,1,3,1,1,2,1,1,2,1/
        data htype(1:siz( 8)+1, 8)  /2,1,1,3,1,1,1,1,1/
        data htype(1:siz( 9)+1, 9)  /2,1,1,3,1,1,1,1,2,1/
        data htype(1:siz(10)+1,10)  /2,1,1,3,1,1,1,1,1/
        data htype(1:siz(11)+1,11)  /2,1,1,3,1,1,3,1,1/
        data htype(1:siz(12)+1,12)  /2,1,1,3,1,1,3,2,1/
        data htype(1:siz(13)  ,13)  /1,1,1,3,1,1,1/
        data htype(1:siz(14)+1,14)  /2,1,1,3,1,1,1,3,2,1/
        data htype(1:siz(15)+1,15)  /2,1,1,3,1,1,1,2,1,2,2,1/
        data htype(1:siz(16)+1,16)  /2,1,1,3,1,4,1/
        data htype(1:siz(17)+1,17)  /2,1,1,3,1,4,1,1/
        data htype(1:siz(18)+1,18)  /2,1,1,3,1,1,1,1/
        data htype(1:siz(19)+1,19)  /2,1,1,3,1,1,1,1,2,1,1,1,1,1,1/
        data htype(1:siz(20)+1,20)  /2,1,1,3,1,1,1,1,1,1,1,4,1/
        data htype(1,21) /3/
!c  this is the Engh and Huber amino acid, all trans.
      data Axyz/   0.682, -0.382, -1.068,   & !n
                  -0.370, -0.760, -0.128,   & !ca
                  -1.060,  0.466,  0.445,   & !c
                  -0.718,  1.598,  0.104,   & !o
                   0.197, -1.614,  0.980,   & !cb
                  -2.039,  0.237,  1.318,   & !n+1
                  -2.780,  1.327,  1.940,   & !ca+1 
                  -3.839,  0.799,  2.893,   & !c+1
                   0.852,  0.657, -1.256/     ! 1H
!** Template amino acids coordinates with i+1 anchor atoms
        data AALIB(1:3,1:siz( 1)+4, 1)   / 0.000,  0.000,  0.000,  & ! A
                                           1.400,  0.000,  0.000,  & ! CA
                                           2.048,  1.346,  0.000,  & ! C
                                           1.305,  2.427,  0.000,  & ! O
                                           1.976, -0.706, -1.223,  & ! CB
                                           3.522,  1.460,  0.000,  & ! N+1
                                           4.193,  2.689,  0.000,  & ! CA+1
                                           5.685,  2.613,  0.000,  & ! C+1
                                           -0.588,  0.906,  0.000/  & ! 1H
        data AALIB(1:3,1:siz( 2)+4, 2)   /  0.000,  0.000,  0.000,  ! C
                                           1.470,  0.000,  0.000,  ! CA
                                           1.993,  1.438,  0.000,  ! C
                                           1.213,  2.401,  0.000,  ! O
                                           2.018, -0.714, -1.237,  ! CB
                                           3.848, -0.698, -1.209,  ! SG
                                           3.452,  1.682,  0.000,  ! N+1
                                           4.038,  3.030,  0.000,  ! CA+1
                                           5.565,  2.936,  0.000,  ! C+1
                                           -0.588,  0.906,  0.000/  ! 1H
        data AALIB(1:3,1:siz( 3)+4, 3)   /  0.000,  0.000,  0.000,  ! D
                                           1.430,  0.000,  0.000,  ! CA
                                           1.903,  1.455,  0.000,  ! C
                                           1.021,  2.369,  0.000,  ! O
                                           1.964, -0.696, -1.205,  ! CB
                                           3.478, -0.829, -1.191,  ! CG
                                           4.182, -0.812, -2.100,  ! OD1
                                           3.898, -0.969,  0.020,  ! OD2
                                           3.344,  1.788,  0.000,  ! N+1
                                           3.833,  3.132,  0.000,  ! CA+1
                                           5.362,  3.078,  0.000,  ! C
                                           -0.588,  0.906,  0.000/  ! 1H
        data AALIB(1:3,1:siz( 4)+4, 4)   /   0.000,  0.000,  0.000,  ! E
                                           1.494,  0.000,  0.000,  ! CA
                                           1.945,  1.465,  0.000,  ! C
                                           1.116,  2.367,  0.000,  ! O
                                           2.046, -0.719, -1.246,  ! CB
                                           1.626, -0.087, -2.551,  ! CG
                                           2.067, -0.882, -3.759,  ! CD
                                           1.833, -0.552, -4.882,  ! OE1
                                           2.726, -1.961, -3.419,  ! OE2
                                           3.387,  1.792,  0.000,  ! N+1
                                           3.737,  3.244,  0.000,  ! CA+1
                                           5.267,  3.339,  0.000,  ! C+1
                                           -0.588,  0.906,  0.000/  ! 1H
        data AALIB(1:3,1:siz( 5)+4, 5)  /   0.000,  0.000,  0.000,  ! F
                                           1.436,  0.000,  0.000,  ! CA
                                           1.942,  1.422,  0.000,  ! C
                                           1.092,  2.356,  0.000,  ! O
                                           2.008, -0.722, -1.251,  ! CB
                                           1.543, -2.112, -1.397,  ! CG
                                           0.380, -2.450, -2.048,  ! CD1
                                           2.281, -3.150, -0.868,  ! CD2
                                           -0.042, -3.745, -2.177,  ! CE1
                                           1.905, -4.470, -0.968,  ! CE2
                                           0.730, -4.760, -1.631,  ! CZ
                                           3.392,  1.713,  0.000,  ! N+1
                                           3.922,  3.047,  0.000,  ! CA+1
                                           5.430,  2.993,  0.000,  ! C+1
                                           -0.588,  0.906,  0.000/  ! 1H
         data AALIB(1:3,1:siz( 6)+4, 6)  / -0.023, -0.051,  0.000,  ! G
                                           1.425,  0.030,  0.000,  ! Ca
                                           1.996,  1.460,  0.000,  ! C
                                           1.095,  2.355,  0.000,  ! O
                                           3.444,  1.760,  0.000,  ! N+1
                                           3.641,  3.197,  0.000,  ! CA+1
                                           5.112,  3.653,  0.000,  ! C+1
                                           -0.440,  0.945,  0.000/  ! 1H
        data AALIB(1:3,1:siz( 7)+4, 7)   /  0.000,  0.000,  0.000,  ! H
                                           1.470,  0.000,  0.000,  ! CA
                                           1.993,  1.438,  0.000,  ! C
                                           1.213,  2.401,  0.000,  ! O
                                           1.917, -0.732, -1.267,  ! CB
                                           2.690, -1.861, -1.025,  ! CG
                                           3.062, -2.344,  0.225,  ! ND1
                                           3.225, -2.707, -1.989,  ! CD2
                                           3.828, -3.488,  0.032,  ! CE1
                                           3.929, -3.712, -1.336,  ! NE2
                                           3.452,  1.682,  0.000,  ! N+1
                                           4.038,  3.030,  0.000,  ! CA+1
                                           5.565,  2.936,  0.000,  ! C+1
                                           -0.588,  0.906,  0.000/  ! 1H
        data AALIB(1:3,1:siz( 8)+4, 8)   /  0.003,  0.006,  0.000,  ! I
                                           1.463, -0.007,  0.000,  ! CA
                                           1.931,  1.439,  0.000,  ! C
                                           1.100,  2.360,  0.000,  ! O
                                           2.043, -0.719, -1.224,  ! CB
                                           1.529, -2.148, -1.285,  ! CG1
                                           3.560, -0.743, -1.142,  ! CG2
                                           1.903, -2.916, -2.555,  ! CD1
                                           3.377,  1.749,  0.000,  ! N+1
                                           3.898,  3.113,  0.000,  ! CA+1
                                           5.415,  3.020,  0.000,  ! C+1
                                           -0.577,  0.917,  0.000/  ! 1H
        data AALIB(1:3,1:siz( 9)+4, 9)  /   0.000,  0.000,  0.000,  ! K
                                           1.470,  0.000,  0.000,  ! CA              
                                           1.993,  1.438,  0.000,  ! C              
                                           1.213,  2.401,  0.000,  ! O              
                                           2.018, -0.714, -1.237,  ! CB              
                                           3.408, -0.696, -1.206,  ! CG              
                                           3.993, -1.403, -2.430,  ! CD              
                                           5.523, -1.383, -2.395,  ! CE              
                                           6.089, -2.066, -3.579,  ! NZ
                                           3.452,  1.682,  0.000,  ! N+1              
                                           4.038,  3.030,  0.000,  ! CA+1              
                                           5.565,  2.936,  0.000,  ! C+1              
                                           -0.588,  0.906,  0.000/  ! 1H
        data AALIB(1:3,1:siz(10)+4,10)  /  -0.006, -0.028,  0.000,  ! L
                                           1.433,  0.012,  0.000,  ! CA
                                           1.971,  1.455,  0.000,  ! C
                                           1.197,  2.424,  0.000,  ! O
                                           2.022, -0.677, -1.221,  ! CB
                                           1.596, -2.135, -1.409,  ! CG
                                           2.288, -2.864, -2.549,  ! CD1
                                           1.800, -2.841, -0.056,  ! CD2
                                           3.431,  1.688,  0.000,  ! N+1
                                           3.789,  3.083,  0.000,  ! CA+1
                                           5.315,  3.287,  0.000,  ! C+1
                                           -0.440,  0.901,  0.000/  ! 1H
        data AALIB(1:3,1:siz(11)+4,11)  /   0.000,  0.000,  0.000,  ! M
                                           1.430,  0.000,  0.000,  ! CA              
                                           1.877,  1.463,  0.000,  ! C              
                                           1.048,  2.385,  0.000,  ! O              
                                           1.982, -0.719, -1.245,  ! CB              
                                           1.688, -2.261, -1.275,  ! CG              
                                           2.678, -3.062, -2.460,  ! SD
                                           1.912, -4.542, -2.568,  ! CE              
                                           3.321,  1.783,  0.000,  ! N+1
                                           3.822,  3.123,  0.000,  ! CA+1
                                           5.349,  3.029,  0.000,  ! C+1
                                           -0.588,  0.906,  0.000/  ! 1H
        data AALIB(1:3,1:siz(12)+4,12)  /   0.000,  0.000,  0.000,  ! N              
                                           1.450,  0.000,  0.000,  ! CA              
                                           2.007,  1.436,  0.000,  ! C              
                                           1.250,  2.405,  0.000,  ! O              
                                           1.983, -0.717, -1.242,  ! CB              
                                           1.389, -2.112, -1.354,  ! CG              
                                           0.440, -2.340, -2.102,  ! OD1              
                                           1.961, -3.033, -0.601,  ! ND2
                                           3.470,  1.652,  0.000,  ! N+1
                                           4.074,  2.970,  0.000,  ! CA+1
                                           5.611,  2.879,  0.000,  ! C+1
                                           -0.588,  0.906,  0.000/  ! 1H
        data AALIB(1:3,1:siz(13)+3,13)  /   0.112,  0.684,  0.124,  ! P
                                           1.622,  0.530, -0.039,  ! CA              
                                           2.255,  1.528, -1.005,  ! C              
                                           1.502,  2.364, -1.565,  ! O              
                                           1.850, -0.869, -0.593,  ! CB              
                                           0.598, -1.039, -1.470,  ! CG              
                                           -0.523, -0.429, -0.622,  ! CD              
                                           3.710,  1.512, -1.274,  ! N+1
                                           4.038,  2.559, -2.222,  ! CA+1
                                           5.529,  2.656, -2.594/  ! C+1
        data AALIB(1:3,1:siz(14)+4,14)  /   0.000,  0.000,  0.000,  ! Q
                                           1.449,  0.000,  0.000,  ! CA              
                                           1.970,  1.430,  0.000,  ! C              
                                           1.186,  2.377,  0.000,  ! O              
                                           1.997, -0.714, -1.237,  ! CB              
                                           3.527, -0.714, -1.237,  ! CG              
                                           4.076, -1.428, -2.474,  ! CD              
                                           5.273, -1.545, -2.676,  ! OE1
                                           3.135, -1.899, -3.288,  ! NE2
                                           3.427,  1.682,  0.000,  ! N+1
                                           3.922,  3.044,  0.000,  ! CA+1
                                           5.444,  3.044,  0.000,  ! C+1
                                           -0.477,  0.833,  0.000/  ! 1H
        data AALIB(1:3,1:siz(15)+4,15)  /   0.000,  0.000,  0.000,  ! R
                                           1.449,  0.000,  0.000,  ! CA              
                                           1.970,  1.430,  0.000,  ! C              
                                           1.186,  2.377,  0.000,  ! O              
                                           1.997, -0.714, -1.237,  ! CB              
                                           3.387, -0.696, -1.206,  ! CG              
                                           3.972, -1.403, -2.430,  ! CD              
                                           5.291, -1.432, -2.480,  ! NE
                                           5.998, -2.006, -3.475,  ! CZ
                                           7.317, -1.977, -3.423,  ! NH1
                                           5.100, -2.594, -4.494,  ! NH2
                                           3.427,  1.682,  0.000,  ! N+1
                                           3.922,  3.044,  0.000,  ! CA+1
                                           5.444,  3.044,  0.000,  ! C+1
                                           -0.477,  0.833,  0.000/  ! 1H
        data AALIB(1:3,1:siz(16)+4,16)  /   0.000,  0.000,  0.000,  ! S
                                           1.448,  0.000,  0.000,  ! CA              
                                           2.010,  1.404,  0.000,  ! C              
                                           1.310,  2.417,  0.000,  ! O              
                                           1.980, -0.715, -1.238,  ! CB              
                                           1.778, -0.035, -2.437,  ! OG
                                           3.479,  1.569,  0.000,  ! N+1
                                           3.807,  2.979,  0.000,  ! CA+1
                                           5.302,  3.208,  0.000,  ! C+1
                                           -0.361,  1.018,  0.000/  ! 1H
        data AALIB(1:3,1:siz(17)+4,17)  /   0.000,  0.000,  0.000,  ! T
                                           1.452,  0.000,  0.000,  ! CA              
                                           2.043,  1.421,  0.000,  ! C              
                                           1.321,  2.398,  0.000,  ! O              
                                           1.977, -0.713, -1.235,  ! CB              
                                           1.454, -2.043, -1.142,  ! OG1
                                           3.489, -0.729, -1.262,  ! CG2
                                           3.511,  1.600,  0.000,  ! N+1
                                           4.149,  2.904,  0.000,  ! CA+1
                                           5.685,  2.812,  0.000,  ! C+1
                                           -0.588,  0.906,  0.000/  ! 1H
        data AALIB(1:3,1:siz(18)+4,18)  /   0.000,  0.000,  0.000,  ! V
                                           1.461,  0.000,  0.000,  ! CA              
                                           1.937,  1.438,  0.000,  ! C              
                                           1.133,  2.379,  0.000,  ! O              
                                           2.004, -0.726, -1.257,  ! CB              
                                           1.735, -2.220, -1.212,  ! CG1
                                           3.483, -0.442, -1.389,  ! CG2
                                           3.388,  1.723,  0.000,  ! N+1
                                           3.933,  3.079,  0.000,  ! CA+1
                                           5.445,  2.984,  0.000,  ! C+1
                                           -0.588,  0.906,  0.000/  ! 1H
        data AALIB(1:3,1:siz(19)+4,19)  /   0.310,  8.246,  0.187,  ! W
                                           1.758,  8.246,  0.187,  ! CA              
                                           2.276,  9.677,  0.156,  ! C              
                                           1.508, 10.638,  0.136,  ! O              
                                           2.269,  7.494, -1.054,  ! CB              
                                           1.910,  6.041, -1.022,  ! CG              
                                           0.876,  5.432, -1.653,  ! CD1
                                           2.618,  5.018, -0.301,  ! CD2
                                           0.884,  4.073, -1.373,  ! NE1
                                           1.942,  3.803, -0.548,  ! CE2
                                           3.746,  5.034,  0.517,  ! CE3
                                           2.409,  2.597,  0.035,  ! CZ2
                                           4.191,  3.854,  1.078,  ! CZ3
                                           3.516,  2.657,  0.828,  ! CH2
                                           3.735,  9.917,  0.151,  ! N+1
                                           4.316, 11.243,  0.123,  ! CA+1
                                           5.835, 11.144,  0.125,  ! C+1
                                           -0.235,  9.085,  0.169/  ! 1H
        data AALIB(1:3,1:siz(20)+4,20)  /   0.000,  0.000,  0.000,  ! Y              
                                           1.452,  0.000,  0.000,  ! CA              
                                           2.068,  1.403,  0.000,  ! C              
                                           1.305,  2.372,  0.000,  ! O              
                                           1.994, -0.721, -1.248,  ! CB              
                                           1.559, -2.159, -1.394,  ! CG              
                                           0.531, -2.506, -2.265,  ! CD1
                                           2.158, -3.174, -0.675,  ! CD2
                                           0.115, -3.813, -2.416,  ! CE1
                                           1.753, -4.495, -0.814,  ! CE2
                                           0.731, -4.815, -1.685,  ! CZ
                                           0.346, -6.135, -1.805,  ! OH
                                           3.535,  1.589,  0.000,  ! N+1
                                           4.166,  2.897,  0.000,  ! CA+1
                                           5.697,  2.842,  0.000,  ! C+1
                                           -0.588,  0.906,  0.000/  ! 1H
!** Terminal oxygen. restype=21
        data AALIB(1:3,1,21)  /   0.000,  0.000,  0.000 /             
!** Complile with -ccd (MPW) if you DON'T want to see the H's
!** MONITOR needs to be modified if this is changed!!
!D       data nseg/7,8,11,12,14,5,13,11,11,11,10,11,10,12,14,8,10,10,19,16/
!D       data (Aconn(I,1),I=1,7 )  /1,4,-1,2,3,-2,6/
!D       data (Aconn(I,2),I=1,8 )  /1,4,5,-1,2,3,-2,7/
!D       data (Aconn(I,3),I=1,11)  /1,4,5,6,-5,7,-1,2,3,-2,9/
!D       data (Aconn(I,4),I=1,12)  /1,4,5,6,7,-6,8,-1,2,3,-2,10/
!D       data (Aconn(I,5),I=1,14)  /1,4,5,6,8,10,9,7,5,-1,2,3,-2,12/
!D       data (Aconn(I,6),I=1,5 )  /1,2,3,-2,5/
!D       data (Aconn(I,7),I=1,13)  /1,4,5,6,8,9,7,5,-1,2,3,-2,11/
!D       data (Aconn(I,8),I=1,11 ) /1,4,5,7,-4,6,-1,2,3,-2,9/
!D       data (Aconn(I,9),I=1,11 ) /1,4,5,6,7,8,-1,2,3,-2,10/
!D       data (Aconn(I,10),I=1,11) /1,4,5,6,-5,7,-1,2,3,-2,9/
!D       data (Aconn(I,11),I=1,10) /1,4,5,6,7,-1,2,3,-2,9/
!D       data (Aconn(I,12),I=1,11) /1,4,5,6,-5,7,-1,2,3,-2,9/
!D       data (Aconn(I,13),I=1,10) /1,4,5,6,0,-1,2,3,-2,7/
!D       data (Aconn(I,14),I=1,12) /1,4,5,6,7,-6,8,-1,2,3,-2,10/
!D       data (Aconn(I,15),I=1,14) /1,4,5,6,7,8,9,-8,10,-1,2,3,-2,12/
!D       data (Aconn(I,16),I=1,8 ) /1,4,5,-1,2,3,-2,7/
!D       data (Aconn(I,17),I=1,10) /1,4,5,-4,6,-1,2,3,-2,8/
!D       data (Aconn(I,18),I=1,10) /1,4,5,-4,6,-1,2,3,-2,8/
!D       data (Aconn(I,19),I=1,19) /1,4,5,6,8,9,11,13,12,10,7,9,-7,5,-1,2,3,-2,15/
!D       data (Aconn(I,20),I=1,16) /1,4,5,6,8,10,11,-10,9,7,5,-1,2,3,-2,13/
!** Complile with -ccx (MPW) if you DO want to see the H's
!CX      data nseg/9,10,13,14,16,7,15,13,13,13,12,13,12,14,16,10,12,12,21,18/
!CX      data (Aconn(I,1),I=1,9 )  /-5,0,1,4,-1,2,3,-2,6/
!CX      data (Aconn(I,2),I=1,10 )  /-6,0,1,4,5,-1,2,3,-2,7/
!CX      data (Aconn(I,3),I=1,13)  /-8,0,1,4,5,6,-5,7,-1,2,3,-2,9/
!CX      data (Aconn(I,4),I=1,14)  /-9,0,1,4,5,6,7,-6,8,-1,2,3,-2,10/
!CX      data (Aconn(I,5),I=1,16)  /-11,0,1,4,5,6,8,10,9,7,5,-1,2,3,-2,12/
!CX      data (Aconn(I,6),I=1,7 )  /-4,0,1,2,3,-2,5/
!CX      data (Aconn(I,7),I=1,15)  /-10,0,1,4,5,6,8,9,7,5,-1,2,3,-2,11/
!CX      data (Aconn(I,8),I=1,13 ) /-8,0,1,4,5,7,-4,6,-1,2,3,-2,9/
!CX      data (Aconn(I,9),I=1,13 ) /-9,0,1,4,5,6,7,8,-1,2,3,-2,10/
!CX      data (Aconn(I,10),I=1,13) /-8,0,1,4,5,6,-5,7,-1,2,3,-2,9/
!CX      data (Aconn(I,11),I=1,12) /-8,0,1,4,5,6,7,-1,2,3,-2,9/
!CX      data (Aconn(I,12),I=1,13) /-8,0,1,4,5,6,-5,7,-1,2,3,-2,9/
!CX      data (Aconn(I,13),I=1,12) /-6,0,1,4,5,6,0,-1,2,3,-2,7/
!CX      data (Aconn(I,14),I=1,14) /-9,0,1,4,5,6,7,-6,8,-1,2,3,-2,10/
!CX      data (Aconn(I,15),I=1,16) /-11,0,1,4,5,6,7,8,9,-8,10,-1,2,3,-2,12/
!CX      data (Aconn(I,16),I=1,10 ) /-6,0,1,4,5,-1,2,3,-2,7/
!CX      data (Aconn(I,17),I=1,12) /-7,0,1,4,5,-4,6,-1,2,3,-2,8/
!CX      data (Aconn(I,18),I=1,12) /-7,0,1,4,5,-4,6,-1,2,3,-2,8/
!CX      data (Aconn(I,19),I=1,21) /-14,0,1,4,5,6,8,9,11,13,12,10,7,9,-7,5,-1,2,3,-2,15/
!CX      data (Aconn(I,20),I=1,18) /-12,0,1,4,5,6,8,10,11,-10,9,7,5,-1,2,3,-2,13/
!** bfby replaces the DIHBYTE routine. Each integer in binary
!** represents the dihedral angles before (bits 0:7) and after (bits 8:15)
!** a given atom, bposit(seq(I))+IO, where I is the residue number,
!** seq(I) is the residue type and IO is
!** the number of the atom in that residue. If the bit J or J+8 is set, the
!** angle phiposit(I)+J rotates that atom and either points away (J) or
!** toward (J+8) atom nposit(I)+IO.
!cb New bfby 31-may-95
        data bposit /1,7,14,23,33,45,50,61,70,80,   &
                    89,98,107,114,124,136,143,151,159,174/
        data bfby /  512,    0,    1,    3,  513,  768,             &
         1536,    0,  513,  517, 1025, 1027, 1792, 3584, 1024, 1537,&
         1545, 2049, 2051, 2055, 2055, 3840, 7680, 3072, 3585, 3601,&
         6145, 4099, 4103, 4111, 4111, 7936, 3584, 1024, 1537, 1545,&
         2049, 2051, 2055, 2055, 2055, 2055, 2051, 3840,  512,    0,&
            1,    3,  768, 3584, 1024, 1537, 1545, 2049, 2051, 2055,&
         2055, 2055, 2055, 3840, 3584, 1024, 1537, 1545, 2049, 2051,&
         3075, 2055, 3840,15872, 7168, 7681, 7713,14337,12291, 8199,&
         8207, 8223,16128, 3584, 1024, 1537, 1545, 2049, 2051, 2055,&
         2055, 3840, 7680, 3072, 3585, 3601, 6145, 4099, 4103, 4111,&
         7936, 3584, 1024, 1537, 1545, 2049, 2051, 2055, 2055, 3840,&
          256,    0,    0,    1,  256,  256,  256, 7680, 3072, 3585,&
         3601, 6145, 4099, 4103, 4111, 4111, 7936,15872, 7168, 7681,&
         7713,14337,12291, 8199, 8207, 8223, 8223, 8223,16128, 1536,&
            0,  513,  517, 1025, 1027, 1792, 1536,    0,  513,  517,&
         1025, 1027, 1027, 1792, 1536,    0,  513,  517, 1025, 1027,&
         1027, 1792, 3584, 1024, 1537, 1545, 2049, 2051, 2055, 2055,&
         2055, 2055, 2055, 2055, 2055, 2055, 3840, 3584, 1024, 1537,&
         1545, 2049, 2051, 2055, 2055, 2055, 2055, 2051, 2051, 3840/
        data tposit /1,3,6,10,15,19,21,25,29,35,39,44,48,49,54,60,63,66,69,73/
!** tortyp(J) = dihedral angle type for J = tposit(iseq) + (iang - phiposit(ires))
        data tortyp /1,2,1,3,2,1,4,5,2,1,4,6,5,2,1,4,7,2,8,9,1,4,10,2,1,   &
        11,12,2,1,4,6,6,13,2,1,4,14,2,1,4,15,16,2,1,4,17,2,2,1,            &
        4,6,17,2,1,4,6,13,18,2,1,19,2,1,20,2,1,11,2,1,4,21,2,1,4,7,2/
!cb        include 'prot_fft.chem.inc'  ! this is now runtime
CONTAINS
!************************************************************************
  subroutine database(J,dang,xyz)
    implicit none
    !****** local variables
    integer(kind=2) ::  J,j1
    real(kind=8) :: xyz(3,18), dang(8)
    !******
    do j1=1,8+nside(J)
      xyz(1,j1) = dble(AALIB(1,j1,J))
      xyz(2,j1) = dble(AALIB(2,j1,J))
      xyz(3,j1) = dble(AALIB(3,j1,J))
    enddo
  end subroutine database
end module ARCHIVERO
