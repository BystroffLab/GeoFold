!!========================================================================================
!! UNFOLDSIM (unfoldsim.f90)
!! This program carries out a finite difference calculation on
!! a directed acyclic graph of bifurcated edges (DAG), generated as
!! output by the program GEOFOLD (geofold.f90)
!!
!! The graph is composed of nodes (ISEGM, intermediates of unfolding) and
!! bifurcating edges (TSTATE, transition states). TSTATEs connect node "f" with two nodes
!! "u1" and "u2". Each node consists of a set of residues in the protein sequence.
!! The residues in f are bisected to give two non-overlapping sets of residues
!! u1 and u2. ISEGMs f, u1, and u2 have associated energies. TSTATE has a type (h,b,p,m)
!! depending on whether the unfolding motion that separates u1 from u2 is a hinge (h),
!! a pivot (p), a break (b), or a melting (m) step that cleaves residues off one at a time.
!!
!! At each time step, for each node, a new concentration is calculated.
!! The change in concentration is calculated based on the amount lost due to
!! unfolding or refolding, and the amount gained by the unfolding or refolding
!! of other nodes. The simulation ends when the concentrations no longer change.
!!
!! Usage: xunfoldsim dagfile paramsfile > output'
!! Output files:
!!   <stdout> contains TIMECOURSE lines for plotting.
!!   <dagfile>.out contains  a new DAG with traffic values for each TSTATE
!!     and other additions
!!   <dagfile>.age contains a contact map with age values for each contact.
!!     The age of a contact is calculated in writeageplot as the time at which
!!     the contact is less than 1/2 populated.
!!   <dagfile>.path contains a contact map with path values for each contact.
!!     The path value is the concentration of a contact at a given point in time.
!!     Both this file and the .age file can be plotted as postscript using
!!     the program pathway2ps.f90.
!!   <dagfile>.nrg contains the ground state (ISTATE) and transition state (TSTATE)
!!     energies for the pathway tree of highest traffic, not counting melting steps.
!!     Each ground state energy is the sum of all ground state energies for
!!     nodes on the maxTraffic tree with the same number of steps from Folded.
!!     Each transition state energy is assigned the highest value of all
!!     barrier heights over all transitions in the maxTraffic tree with the
!!     same number of steps from Folded.
!! 
!! |===== INPUT PARAMETERS for UNFOLDSIM ====='
!! |  1. dagfile = output of GeoFold. DAG format'
!! |  2. paramsfile = keyworded parameters. One per line.'
!! |  Recognized keywords
!! |  BREAKPOINTENTROPY ",breakpointentropy
!! |    Entropy in kJ/mol/deg to be applied to a BREAK move"
!! |  HINGEPOINTENTROPY ",hingepointentropy
!! |    Entropy in kJ/mol/deg to be applied to a HINGE move"
!! |  PIVOTPOINTENTROPY ",pivotpointentropy
!! |    Entropy in kJ/mol/deg to be applied to a PIVOT move"
!! |  TEMPERATURE ",Temperature
!! |    Temperature in Kelvins"
!! |  OMEGA ",omega
!! |    Desolvation energy in kJ/mol/A^2"
!! |  INTERMEDIATE ",intm
!! |    Optional: the number of an intermediate to monitor."
!! |  CONCENTRATION ",conc
!! |    Total protein concentration, in mol/L"
!! |  FOLDING ",fing
!! |    1=folding, 0=unfolding"
!! |  VOIDENTROPY ",spervoid
!! |    Entropy assigned to the loss of a void, in kJ/mol/deg"
!! |  SOLIDITY ",flexibility
!! |    Parameter for scaling configurational entropy to zero as SAS << SOLIDITY"
!! |  HBONDENERGY ",eperhbond
!! |    Energy per hydrogen bond broken, in kJ/mol"
!! |  HAMMONDSCALE ",tmscale
!! |    Inreasing this scales down Hammond shift, "
!! |    which is the position of T-state relative to higher energy ground state."
!! |  SIDECHAINENTROPY ",scfac
!! !===== Keywords recognized by other GEOFOLD programs
!! | BREAKCUT, HINGECUT, PIVOTCUT, ORANGE, REDUCING
!!
!!=============================================================================================
!! MODIFICATIONS:
!! ---------------- Tue Jan  7 11:53:26 EST 2014
!! Increased DCCUTOFF by a factor of 10 and increased the maximum concentration change
!! criterea by a factor of ten. Hopefully this will speed convergence without
!! seriously affecting the results.
!! ---------------- Sat Jul  6 15:48:35 EDT 2013
!! Adding GeoFold seam move. Done Fri Aug 23 17:34:25 EDT 2013
!! ---------------- Mon Jan 18 15:08:12 EST 2010
!! Energy profiles are too flat.  Following Zwanzig, roughly, entropy is dealt
!! out on a per residue basis. A residue is unfolded if it is within pivottail of
!! a pivot, or if it is an unfolded segment.
!! ---------------- Fri Jan  1 17:47:06 EST 2010
!! Pre-equilibration of melting steps installed. Melting steps are ingored when unflding.
!! Melting steps are pre-equilibrated when folding. Thus all concentration is in
!! initiation sites (I-sites) before time starts.
!! ---------------- Thu Dec 31 07:07:20 EST 2009
!! getthetam() has been returned to its former self, the function described in
!! the paper, which uses similar triangles.
!! ---------------- Wed Dec 30 13:46:22 EST 2009 C.Bystroff
!! Added energy output (writeenergyprofile) for plotting using gnuplot.
!! Added entropic barrier heights (getSbarrier, keywords PIVOTBARRIER, HINGEBARRIER)
!! This changes the hingebarrier, making it dependent on geofold entropy, and adds
!! a similar barrier for pivots. 
!! ---------------- Tue Sep 29 06:15:57 EDT 2009  C.Bystroff
!! Tracing version: find out what intermediates are dominant
!! INTERMEDIATE -1 sets the program to output the node number of the highest
!! concentration at each output step. (see intm)
!! ---------------15-Jul-2009 C.Bystroff
!! Added MELTINGENERGY, the energy and energy barrier for a melting step
!! Unmelting is given no barrier
!! ---------------17-NOV-2008 C.Bystroff
!! Input of parameters changed so that they come from a file insted
!! of from the command line. This should make it easier
!! to optimize the program.
!! ---------------21-May-2009 C.Bystroff
!! pivotpointentropy  = ( hingepointentropy + breakpointentropy) / 2.0 is now
!! enforced. This is because entropy is a state function and the 
!! total entropy added cannot be pathway dependent.
!! ---------------3-JUL-2008 C.Bystroff
!! t1/2 changed to mean time at which protein is half-unfolded (unfolding)
!! or half-folded (folding).  This will give the t1/2 of an intermdiate state
!! if the unfolding is not two-state.
!! ---------------Thu Jun 26 21:19:42 EDT 2008 C.Bystroff
!! Solidity = entropy (from geofold)
!! ---------------2007 C.Bystroff, V. Ramakrishnan
!! ---- The following are fixed parameters ------
!! nuk = 10^6 from Fersht book
!! timestep = variable, but should be << nuk
!! NOTE: changing the timestep should not effect the rate.
!!================================================================================
program unfoldsim
  use geofold_global, only : barrels_array, seam_type, button_type , barrel_type,&
      intermediate, MAXBARREL
  use geofold_seams, only : geofold_seams_read
  use geofold_pivots, only : geofold_pivots_queryinseam
  implicit none
  real(8),parameter    :: nuk=1000000
  real(8),parameter    :: Rvalue = 8.314  !! J/mol/K
  REAL(8),parameter    :: DCCUTOFF = 0.00000001, TINYCONC=0.00000001
  REAL(8),parameter    :: SASBOND = 30.0
  REAL(8),parameter    :: FCUTOFF=0.60, UCUTOFF=3000.,speedfac=2.0, speedlimit=10.0
  REAL(8),parameter    :: MELTINGENERGY=1.00
  REAL(8),parameter    :: FLUIDITY=0.00
  !integer,parameter :: logfreq=10000, agefreq=5000, lifefreq=1000, speedfreq=10000
  integer,parameter :: logfreq=10000, agefreq=500, lifefreq=50, speedfreq=10000
  integer,parameter :: MAXLEN=2000, SUMALLFLAG=99, convfreq=logfreq*3
  real(8),parameter    :: TIMESTEPSTART=1.0*10.**(-11)
  !!----- variable declarations ------------------
  real(8)              :: timestep=TIMESTEPSTART, totaltime=0., maxtime=1000.0, mintime=0.001 
  real(8)              :: Temperature = 300   
  real(8)              :: breakpointentropy = 100.   
  real(8)              :: pivotpointentropy = 0.  
  real(8)              :: hingepointentropy = 10
  real(8)              :: seampointentropy = 0.
  REAL(8)              :: SPERRESID=10.00
  real(8)              :: omega=30., decrement, spervoid=4.0, eperhbond=1.0
  real(8)              :: RT,flexibility, sumf,sumu,uconc=0.,iconc=0.,lastf=0.
  real(8)              :: fconc=100., conc=100.00, maxdc, x, rate
  real(8)              :: hingebarrier=0., pivotbarrier=0., breakbarrier=0., seambarrier=0.
  real(8)              :: sumi, sumtot
  character(len=200)   :: inputfile,aline,paramfile
  real(8)              :: dsas,lnk, dnrg, tmscale, scfac=4.183, lnkf, lnku, trf(4), dhb, dlam, dvoid
  real(8),dimension(0:100) :: halflife
  logical              :: halflifeonly=.true., verbose=.true., trafficking=.true.
  !!----- data structures ------------------------
  type tstype
    integer    :: f,u1,u2,pf   !! index to inter
    character  :: cuttype
    real(8)    :: entropy,thetam,ku,kf,TdS,dsas,solid,traffic,dlambda,dvoid,dhb
    real(8)    :: dnrg,dGddu,dGddf,dSdd
  end type
  type (tstype),dimension(:),allocatable,target :: tstate
  type (tstype),pointer                         :: tsptr
  type intype
    character(len=200)    :: residuelist
    real(8)               :: solv,sas,conc,lambda
    integer               :: folded,nres,void,hb,sym,barrel(MAXBARREL)
    character(len=MAXLEN) :: flags
    type (intermediate), pointer   :: f
  end type
  type (intype),dimension(:),allocatable,target :: inter
  type (intype),pointer                         :: u1ptr,u2ptr,fptr
  !!----- miscelaneous variables -----------------
  integer                                :: iargc, jarg, its, fstate, iflag,global_nres,fing, ihalf, ix,j
  integer                                :: mtstate, minter, intm=0, ios, inputstatus, icycle 
  integer                                :: i, ii, pwayout=0, f, u1, u2
  real(8),dimension(:),allocatable       :: dc
  real(8),dimension(:,:),allocatable     :: age
  real(4),dimension(3,MAXLEN)            :: calpha
  logical :: hasbarrels=.false.
  real(4) :: cpustart, cpuend, xin, xout, tmonitor, tmonitor2, y
  !! ---------- initialization ---------------------------------
  RT = Temperature*Rvalue
  flexibility = 1000.  !! 100. A^2 drops one log unit from solid to liquid.
  breakpointentropy = 100.
  Temperature = 270.
  omega = 5.
  conc = 1.00
  fconc = conc
  uconc = 0.
  iconc = 0.
  intm = 0
  fing = 0
  tmscale = 10000.0   !! possibly obsolete
  spervoid = 0.0
  eperhbond = 0.0
  scfac = 0.0
  hingebarrier = 0.
  !! ---------- command line ---------------------------------
  jarg = iargc()
  if (jarg < 2) then
    write(*,*) ""
    write(*,*) 'Usage: xunfoldsim dagfile paramsfile > output'
    write(*,*) ""
    write(*,*) '|===== input parameters for UnfoldSim ====='
    write(*,*) '|  1. dagfile = output of GeoFold. Revised for GeoFold2 .dag format'
    write(*,*) '|  2. paramsfile = keyworded parameters. One per line.'
    write(*,*) '|  ------ Default Parameters ----- '
    write(*,*) "|  BREAKPOINTENTROPY ",breakpointentropy
    write(*,*) "|    Entropy in kJ/mol/deg to be applied to a BREAK move"
    write(*,*) "|  HINGEPOINTENTROPY ",hingepointentropy
    write(*,*) "|    Entropy in kJ/mol/deg to be applied to a HINGE move"
    write(*,*) "|  PIVOTPOINTENTROPY ",pivotpointentropy
    write(*,*) "|    Entropy in kJ/mol/deg to be applied to a PIVOT move "
    write(*,*) "|    (Ignored! PIVOT set to average of BREAK and HINGE) "
    write(*,*) "|  SEAMENTROPY ",seampointentropy
    write(*,*) "|    Entropy in kJ/mol/deg to be applied to a SEAM move "
    write(*,*) "|  TEMPERATURE ",Temperature
    write(*,*) "|    Temperature in Kelvins"
    write(*,*) "|  OMEGA ",omega
    write(*,*) "|    Desolvation energy in kJ/mol/A^2"
    write(*,*) "|  INTERMEDIATE ",intm
    write(*,*) "|    Optional: if /=0, monitor the intermediate of largest conc."
    write(*,*) "|  CONCENTRATION ",conc
    write(*,*) "|    Total protein concentration, in mol/L"
    write(*,*) "|  FOLDING ",fing
    write(*,*) "|    1=folding, 0=unfolding"
    write(*,*) "|  VOIDENTROPY ",spervoid
    write(*,*) "|    Entropy assigned to the loss of a void, in kJ/mol/deg"
    write(*,*) "|  SOLIDITY ",flexibility
    write(*,*) "|    Parameter for scaling configurational entropy to zero as SAS << SOLIDITY"
    write(*,*) "|  HBONDENERGY ",eperhbond
    write(*,*) "|    Energy per hydrogen bond broken, in kJ/mol"
    write(*,*) "|  HAMMONDSCALE ",tmscale
    write(*,*) "|    Increasing this scales down Hammond shift, "
    write(*,*) "|    which is the position of T-state relative to higher energy ground state."
    write(*,*) "|  SIDECHAINENTROPY ",scfac
    write(*,*) "|    Scale factor for sidechain entropy change upon unfolding."
    write(*,*) "|  HINGEBARRIER  ",hingebarrier
    write(*,*) "|    Energy to be added to the barrier height for hinge moves."
    write(*,*) "|  PIVOTBARRIER  ",pivotbarrier
    write(*,*) "|    Energy to be added to the barrier height for pivot moves."
    write(*,*) "|  BREAKBARRIER  ",breakbarrier
    write(*,*) "|    Energy to be added to the barrier height for break moves."
    write(*,*) "|  SEAMBARRIER  ",seambarrier
    write(*,*) "|    Energy to be added to the barrier height for seam moves."
    write(*,*) "|  SPERRESID  ",SPERRESID
    write(*,*) "|    An additional per-residue entropy, added to unfolded states by residue."
    write(*,*) '|===== unfoldsim.f90 v.  Wed Jun 18 16:10:49 EDT 2014'
    stop
  endif
  call getarg(1,inputfile)
  call getarg(2,paramfile)
  write(*,'(a,a," ",a)') "xunfoldsim ",trim(inputfile),trim(paramfile)
  !!--------- read parameters file ---------------------------------------------
  call cpu_time(cpustart)
  call readparams(paramfile)
  !!---------- initialize concentrations ---------------------------------------
  RT = Temperature*Rvalue
  if (fing==1) then
    fconc = 0.
    uconc = conc
  else
    uconc = 0.
    fconc = conc
  endif
  !!-------- read DAG file -----------------------------------------------------
  call readDAGfile(inputfile)
  call geofold_seams_read(inputfile)
  call setfolded()
  !! ------------------------------  allocations --------------------------------
  global_nres = inter(1)%nres
  allocate(dc(minter),stat=ios) ; if (ios/=0) stop 'Error allocating dc'
  allocate(age(global_nres,global_nres),stat=ios) ; if (ios/=0) stop 'Error allocating age'
  age = HUGE(age)
  call cpu_time(cpuend)
  write(*,'("CPU_TIME spent reading input files   =",f9.4," seconds.")') cpuend - cpustart
  cpustart = cpuend
  !!------------------ Create LOG file for inspection of tstates.  --------------
  if (verbose) then
    open(99,file=trim(inputfile)//".LOG",status="replace",form='formatted',iostat=ios)
    if (ios/=0) stop 'unfoldsim.f90 :: ERROR cant creat file. Permissions?'
  endif
  icycle = 0
  !! -------------- Initialize all transition state energies and rates -----------
  do its=1,mtstate
    u1 = tstate(its)%u1
    u2 = tstate(its)%u2
    f = tstate(its)%f
    !!------------ Initialize traffic
    tstate(its)%traffic = 0
    !!------------ Calculate the changes in various values at the transition state.
    if (tstate(its)%cuttype=="s") then !! If move is a seam, then there is only u1, not u2.
      dsas = inter(u1)%solv
      dlam = inter(u1)%lambda 
      dvoid = inter(u1)%void 
      dhb = inter(u1)%hb  
      tstate(its)%u2=0
      u2 = 0
    else !! If move is not a seam, then sum the values for u1 and u2.
      dsas = inter(u1)%solv + inter(u2)%solv
      dlam = inter(u1)%lambda + inter(u2)%lambda
      dvoid = inter(u1)%void + inter(u2)%void
      dhb = inter(u1)%hb   + inter(u2)%hb  
    endif
    !! ------------------- calculate delta-SAS 
    tstate(its)%dsas    = inter(f)%solv - dsas
    !! ------------------- calculate change in SC entropy
    tstate(its)%dlambda = inter(f)%lambda - dlam
    !! ------------------- calculate change in Void entropy
    tstate(its)%dvoid   = (inter(f)%void - dvoid)*spervoid
    !! ------------------- calculate change in H-bond enthalpy
    tstate(its)%dhb     = (inter(f)%hb   - dhb)*eperhbond
    !! ------------------- calculate solidity, the fraction of TdS expressed before the t-state. Max is 0.5.
    tstate(its)%solid   = Solidity(tstate(its)%dsas,tstate(its)%entropy,flexibility)
    !! ------------------- calculate change in entropic term, TdS
    tstate(its)%TdS     = Temperature*dSfunction(tstate(its),inter(1)%hb)
    !! ------------------- calculate equilibrium delta-Gu = wa dA + wh dHB - L wl T  - dV wv T - T dS 
    dnrg =    tstate(its)%dsas*omega                 &  !! desolvation free energy, enthalpic part
            + tstate(its)%dhb                        &  !! change in H-bonds, an enthalpic term.
            - Temperature*(tstate(its)%dlambda*scfac &  !! side chain entropy gain
                           + tstate(its)%dvoid)         !! void filling entropy gain
    if (tstate(its)%cuttype/="m".and.dnrg==0.00) then
      write(*,*) 'unfoldsim.f90:: WARNING: zero energy difference. TSTATE=',its,' cuttype ',tstate(its)%cuttype
      dnrg = 0.0001
      ! stop 'unfoldsim.f90:: main: ERROR. zero energy difference.'
    endif
    tstate(its)%dnrg = dnrg - tstate(its)%TdS           !! move-dependent entropy gain.
    !! ------------------- calculate thetam  using similar triangles (see paper)
    tstate(its)%thetam  = getthetam(dnrg,tstate(its)%TdS)
    !! ------------------- calculate UNFOLDING barrier and rate: dGddu, ku, -ln(Ku) = dGddu/RT 
    !! NOTE: lnku is -ln(ku)
    lnku = ( ( tstate(its)%dsas*omega +                    &   !!  desolvation free energy (>0)
               tstate(its)%dhb -                           &   !!  H-bonds lost (>0)
               Temperature*( tstate(its)%dlambda*scfac  +  &   !!  side chain entropy gained (>0)
                             tstate(its)%dvoid )           &   !!  entropy of voids filled  (>0)
             )*tstate(its)%thetam -                        &   !!  position of the transition state (0..1)
             tstate(its)%TdS*tstate(its)%solid             &   !!  size,move-dependent entropy gain (>0)
           ) /RT
    !if (dnrg>0.) then
    !  lnku = dnrg
    !else
    !  lnku = 0
    !endif
    !! ------------------- additional barrier height, for exploratory purposes.
    ! if (tstate(its)%cuttype == "h") then
    !   lnku = lnku + hingebarrier
    ! elseif (tstate(its)%cuttype == "m") then
    !   lnku = MELTINGENERGY
    ! endif
    !! --------- add entropic barrier to unfolding. 
    if (tstate(its)%cuttype == "m") then
      !! --------- Melting is a special case. Energy is ignored.
      lnku = MELTINGENERGY - FLUIDITY
      tstate(its)%dSdd = 0.
    else
      !! --------- additional transition state entropy is defined by getSbarrier()
      tstate(its)%dSdd = getSbarrier(tstate(its))
      lnku = lnku + tstate(its)%dSdd/Rvalue
    endif
    tstate(its)%dGddu = lnku*RT
    !! ------------------- calculate FOLDING barrier and rate: dGddf, kf, -ln(Kf) = dGddf/RT 
    !! NOTE: lnkf is -ln(kf)
    !! dGddu - dGddf = dnrg  (dd is double-dagger, u is unfolding, f is folding)
    !! therefore, RT*lnku - RT*lnkf = dnrg !! sas*tm + sas*(1-tm) = sas, check!  
    !! -dlT*tm - dlT(1-tm) = -dlT, check!  !! -dvT*tm - dvT(1-tm) = -dvT, check!
    !! hb*tm + hb*(1-tm) = hb,  check!  !! -TdS*sol -  TdS(1-sol) = -TdS, check!
    lnkf = ( (-tstate(its)%dsas*omega -                    &   !! solvation free energy (<0)
               tstate(its)%dhb +                           &   !! H-bonds made energy (<0)
               Temperature*( tstate(its)%dlambda*scfac  +  &   !! -side chain entropy gained (>0)
                             tstate(its)%dvoid )           &   !! -entropy of void filled (>0)
             )*(1. - tstate(its)%thetam) +                 &   !! position of t-state of folding (0..1)
             tstate(its)%TdS*(1-tstate(its)%solid)         &   !! size,move-dependent entropy gain (>0)
           ) /RT
    !if (dnrg>0.) then
    !  lnkf = 0.
    !else
    !  lnkf = -dnrg
    !endif
    !! ------------------- additional barrier height, for exploratory purposes. both sides.
    ! if (tstate(its)%cuttype == "h") then
    !   lnkf = lnkf + hingebarrier
    ! elseif (tstate(its)%cuttype == "m") then
    !   lnkf = 0.
    ! endif
    !! --------- add entropic barrier to folding. 
    if (tstate(its)%cuttype == "m") then
      !! --------- Melting is a special case. Energy is ignored.
      lnkf = -FLUIDITY
      tstate(its)%dSdd = 0.
    else
      lnkf = lnkf + tstate(its)%dSdd/Rvalue
      !! ----- cases where folding/unfolding is diffusion controlled ------- !!
      if (lnku < 0.) then 
        lnkf = lnkf - lnku
        lnku = 0.
        icycle = icycle + 1
      elseif (lnkf < 0.) then
        lnku = lnku - lnkf
        lnkf = 0.
        icycle = icycle + 1
      endif
    endif
    tstate(its)%dGddf = lnkf*RT
    !! ---- calculate folding/unfolding rates, accounting for multiple syms  ----- !!
    tstate(its)%ku      = nuk * exp(-lnku) * (inter(tstate(its)%f)%sym + 1)
    if (tstate(its)%cuttype=="s") then  !! seam move
      tstate(its)%kf      = nuk * exp(-lnkf) * (inter(tstate(its)%u1)%sym + 1)
    else
      tstate(its)%kf      = nuk * exp(-lnkf) * (inter(tstate(its)%u1)%sym + 1) * &
      (inter(tstate(its)%u2)%sym + 1) 
    endif
    !! diagnostic------------------------ !! write out TSTATE data ------------------
    if (verbose) then
      if (mod(its-1,10)==0) then
        write(99,'(100("-"))') 
        write(99,'(a,a)') &
        "      n    f   u1   u2 C   entropy    thetam        ku        kf       TdS      ", &
        "dSAS  solidity   dlambda     dvoid     dHbnd   Nf  Nu1  Nu2       dnrg  res_list"
        write(99,'(100("-"))') 
      endif
      ! integer    :: f,u1,u2,pf   !! index to inter
      ! character  :: cuttype
      ! real(8)    :: entropy, thetam, ku, kf, TdS,dsas,solid,traffic,dlambda,dvoid,dhb
      write(99,'(i7,3i5,1x,a1,10(1pe10.2e2),3i5,1pE11.3e2,a)') &
        its,tstate(its)%f,  &
        tstate(its)%u1,  &
        tstate(its)%u2,  &
        tstate(its)%cuttype,  &
        tstate(its)%entropy,tstate(its)%thetam,tstate(its)%ku,tstate(its)%kf,tstate(its)%TdS,  &
        tstate(its)%dsas,tstate(its)%solid,tstate(its)%dlambda,tstate(its)%dvoid,  &
        tstate(its)%dhb, &
        inter(tstate(its)%f)%nres, &
        inter(tstate(its)%u1)%nres, &
        inter(tstate(its)%u2)%nres, &
        dnrg, &  
        trim(inter(tstate(its)%f)%residuelist)
    endif
  enddo
  if (verbose) close(99)
  write(*,*) icycle," of ",mtstate," transitions are diffusion controlled."
  call cpu_time(cpuend)
  write(*,'("CPU_TIME spent initializing params   =",f9.4," seconds.")') cpuend - cpustart
  cpustart = cpuend
  tmonitor = 0.000
  tmonitor2 = 0.000
  !! diagnostic----------------------------- !! write out INTERMEDIATE data ------------------------
  !do i=1,minter
  !  write(*,*) "INTERMEDIATE",i,inter(i)%folded, inter(i)%conc, inter(i)%nres, inter(i)%sym
  !enddo
  !! ------------------------------- pre-equilibrate melting I-sites if FOLDING --------------------
  !! This block of code pushes unfolded state concentration to I-sites,
  !! which are the smallest intermediates that result from "b", "p" or "h" moves.
  !! The idea was to establish a more realistic model of the unfolded state as an
  !! ensemble of short fragments, not individaul residues.
  !! The code fails to push all concentration to I-sites, presumably because
  !! each tstate requires both u1 and u2, and if one is depleted to zero the other
  !! can't move. As a result, some melted states go to zero, but not all.
  !! In a quick experiment, folding runs were done with and without preequilibrating the
  !! melted states. It made little difference in the half-life of folding. 
  !! However, if the melted states were made inaccessible, the protein could not fold.
  !! Also, if the MELTINGENERGY was set to a very high number (>10.) then folding
  !! was significanntly slowed. If MELTINGENERGY was set to zero, on the other hand,
  !! folding was slow. A low value for MELTINGENERGY (2.0) was optimal for
  !! giving the shortest halflife of folding in the testcase. Why?
  !! I think this is because unfolded states need to fluidly interchange in response
  !! to depletion by folding. If they cannot melt and reform as a different I-site,
  !! then folding stalls because u1 or u2 is depleted. 
  !! Perhaps increasing the FLUIDITY will help? (currently set to zero)
  !! -----------------------------------------------------------------------------------------
  !if (fing==3) then   !! fing is never 3, This block is permanently off !
  !  timestep = 0.0001
  !  icycle = 0
  !  PREEQ: do
  !    icycle = icycle + 1
  !    !! -------- calculate changes in concentration
  !    dc = 0 
  !    do its=1,mtstate
  !      if (tstate(its)%cuttype == "m") then  
  !        if (inter(tstate(its)%u1)%conc>TINYCONC.AND.inter(tstate(its)%u2)%conc>TINYCONC) then 
  !          rate = inter(tstate(its)%u1)%conc*inter(tstate(its)%u2)%conc*tstate(its)%kf
  !          decrement = rate
  !          dc(tstate(its)%f) = dc(tstate(its)%f) + decrement 
  !          dc(tstate(its)%u1) = dc(tstate(its)%u1) - decrement
  !          dc(tstate(its)%u2) = dc(tstate(its)%u2) - decrement
  !        endif
  !        if (inter(tstate(its)%f)%conc > TINYCONC) then  
  !          rate = inter(tstate(its)%f)%conc * tstate(its)%ku 
  !          decrement = rate 
  !          dc(tstate(its)%f) = dc(tstate(its)%f) - decrement 
  !          dc(tstate(its)%u1) = dc(tstate(its)%u1) + decrement
  !          dc(tstate(its)%u2) = dc(tstate(its)%u2) + decrement
  !        endif
  !      endif
  !    enddo
  !    !! ----- would concentration be negative? If so, reduce timestep
  !    x = timestep
  !    do ii=1,minter
  !      if (-dc(ii)*x > inter(ii)%conc ) then
  !        x = inter(ii)%conc/(-dc(ii)) 
  !      endif
  !    enddo
  !    if (x>0.00) timestep = x
  !    !! ----- update concentrations. This is where the timestep comes in.
  !    maxdc = 0.
  !    do ii=1,minter
  !      inter(ii)%conc = inter(ii)%conc + dc(ii)*timestep
  !      x = abs(dc(ii))
  !      if (x > maxdc) then; maxdc = x ; endif
  !    enddo 
  !    !! -------- check for convergence, are concentrations changing?
  !    if ((maxdc*100000.)<1) then
  !       write(0,*) "PRE-EQUILIBRATED before FOLDING. ",&
  !       "Maximum change =",maxdc*1000000000.," ppb. Cycles=",icycle
  !       exit PREEQ
  !    elseif (icycle > 1000000) then
  !       write(0,*) "WARNING: failed to pre-equilibrate I-sites. ",&
  !       "Going on anyway. Cycles=",icycle
  !       exit PREEQ
  !    endif
  !    if (mod(icycle,logfreq)==0.and.intm > 0) then
  !      write(*,'(i9,i5,f8.3,i12)') icycle, intm, &
  !      inter(intm)%conc*100./conc,nint(maxdc*1000000000)
  !    endif
  !  enddo PREEQ
  !  !!---- print out molten intermediates. They should all have zero conc.
  !  if (verbose) then
  !    do its=1,mtstate
  !      if (tstate(its)%cuttype == "m") then  
  !        write(*,'("TSTATE ",4i5,3f9.6," fsize=",i4)') &
  !            ii,tstate(its)%f,tstate(its)%u1,tstate(its)%u2,&
  !            inter(tstate(its)%f)%conc,   &
  !            inter(tstate(its)%u1)%conc,   &
  !            inter(tstate(its)%u2)%conc,   &
  !            inter(tstate(its)%f)%nres   
  !      endif
  !    enddo
  !  endif
  !endif
  !! ------------------------------- main loop : SIMULATE UNFOLDING --------------------------------
  timestep=TIMESTEPSTART
  icycle = 0
  totaltime = 0.
  iflag = 0
  halflife = 0
  ihalf = 0
  LOOP1: do
    icycle = icycle + 1
    !! ------- increase timestep if possible
    if (icycle >= 10000 .and. mod(icycle,speedfreq)==0) then
      if (timestep < speedlimit) then
        timestep = speedfac*timestep
      endif
      iflag = 0
    endif
    totaltime = totaltime + timestep
    !! -------- calculate changes in concentration
    dc = 0 
    call cpu_time(xin)
    do its=1,mtstate
      f = tstate(its)%f
      u1 = tstate(its)%u1
      u2 = tstate(its)%u2
      fptr => inter(f)
      u1ptr => inter(u1)
      u2ptr => inter(u2)
      tsptr => tstate(its)
      if (inter(f)%conc > TINYCONC) then  
        rate = fptr%conc * tsptr%ku 
        decrement = rate 
        dc(f) = dc(f) - decrement 
        dc(u1) = dc(u1) + decrement
        if (u2/=0) dc(u2) = dc(u2) + decrement
        tsptr%traffic = tsptr%traffic + decrement
      endif 
      if (u1ptr%conc>TINYCONC) then
        rate = 0.0
        if (tsptr%cuttype=="s") then     !! seam is a unimolecular reaction
          rate = u1ptr%conc*tsptr%kf
        else  !! All other cuttypes m, b, h, p are (virtually) bimolecular.
          if (u2ptr%conc>TINYCONC) then 
            rate = u1ptr%conc*u2ptr%conc*tsptr%kf
          endif
        endif
        decrement = rate
        dc(f) = dc(f) + decrement 
        dc(u1) = dc(u1) - decrement
        if (u2/=0) dc(u2) = dc(u2) - decrement
        tsptr%traffic = tsptr%traffic + decrement 
      endif
    enddo
    iflag = 0
    !! ----- would concentration be negative? If so, reduce timestep
    x = timestep
    do ii=1,minter
      !! find minimum timestep
      if (-dc(ii)*x > inter(ii)%conc ) then
        x = inter(ii)%conc/(-dc(ii)) 
      endif
    enddo
    if (x>0.0.and.x<timestep) timestep = x
    !! -------- apply changes in concentration, track maximum relative change
    maxdc = 0.
    ILOOP2: do ii=1,minter
      u1ptr => inter(ii)
      u1ptr%conc = u1ptr%conc + dc(ii)*timestep
      if (u1ptr%conc < 0.00000) then
        !! Warning: Negative conc is possible source of loss-of-mass bug.
        !! diagnostic: look for significant loss of mass.
        if (u1ptr%conc < -TINYCONC) then
          write(0,'("WARNING: negative concentration ",1pe10.3e2," for interme'//&
          'diate ",i4," SETTING IT TO ZERO.")') u1ptr%conc, ii
        endif
        u1ptr%conc = 0.
      endif
      if (dc(ii) < 0.) then
        if ( u1ptr%conc > 0.00) x = -timestep*dc(ii)/u1ptr%conc
        if (x > maxdc) then; maxdc = x; ix = ii; endif
      endif
    enddo ILOOP2
    !! -------- check for convergence, are concentrations changing?
    if (totaltime>mintime .and. (maxdc*100000.)<1) then
       write(0,*) "CONVERGED. Maximum concentration change =",maxdc*1000000000.," ppb"
       exit    !! 1 ppm change. converged
    elseif (totaltime>mintime .and. sum(abs(dc)) < DCCUTOFF) then
       write(0,*) "CONVERGED. sum(abs(dc))=",sum(abs(dc))
       exit
    endif
    call cpu_time(xout)
    tmonitor = tmonitor + (xout - xin)
    xin = xout
    !! ----------- Report the half-life of inter-residue contacts
    if (mod(icycle,agefreq)==0) call ageplot(totaltime,fing)
    if (mod(icycle,10*agefreq)==0) call writeageplot(inputfile)
    !! ----------- Report time and percent folded, unfolded, intemediate.
    if (mod(icycle,logfreq)==0) then
      ! if (intm/=0)  write(*,'(a,999e10.3e2)') "ALL ISEGMT",inter(:)%conc
      ! sumf = sum(inter(:)%conc, dim=1, mask=(inter(:)%folded==1))
      sumf = getsumf()*100./conc
      sumu = getsumu()*100./conc
      sumi = getsumc(0)*100./conc
      sumtot = getsumc(SUMALLFLAG)*100./conc
      write(*,'(a,E12.5e1,3(1x,f9.5),$)') "TIMECOURSE ", totaltime, sumf, sumu , sumi
      !write(*,'(a,E12.5e1,4(1x,f9.5),i12,$)') "TIMECOURSE ", totaltime, sumf, sumu, &
      !                                         sumi, sumtot, nint(maxdc*100000000.)
      !! ----------- report running halflife
      if (fing==1) then
        if (nint(sumf/2.)> 100) stop 'unfoldsim.f90:: BUG. sumf > 200%'
        write(*,'(E12.5e1,$)') halflife(nint(sumf/2.))
      else
        if (nint(sumu/2.)> 100) stop 'unfoldsim.f90:: BUG. sumu > 200%'
        write(*,'(E12.5e1,$)') halflife(nint(sumu/2.))
      endif
      if (trafficking) then
        call currenttraffic(trf)  !! report the relative traffic through b p h .
        write(*,'(4f6.2,$)') trf
      endif
      !! ----------- Report contact concentrations at 50% folded.
      if (pwayout==0) then
        if (fing==1.and.sumf>50.) then
          call writepathway(inputfile)
          pwayout = 1
        elseif (fing==0.and.sumu>50.) then
          call writepathway(inputfile)
          pwayout = 1
        endif
      endif
      !! ----------- monitor residue intm
      !! report largest intermediate concentration (excluding unfolded)
      if (intm /= 0) then
        if (intm<0 .or. intm>minter) then
          j = maxloc(inter(:)%conc, mask=(inter(:)%folded>=0),dim=1)
        else
          j = intm
        endif
        write(*,'(i5,f8.3,$)') j, inter(j)%conc*100./conc
      endif
      write(*,*)
      !! ----------- check for completion
      if (totaltime > maxtime) then
        write(0,*) "MAXIMUM simulation time reached. maxtime = ",maxtime, " seconds."
        exit
      elseif (fing==1.and.sumf>99.) then
        write(0,*) "FOLDING simulation complete. sumf= ", sumf
        exit
      elseif (fing==0.and.sumu>99.) then
        write(0,*) "UNFOLDING simulation complete. sumu= ", sumu
        exit
      elseif (mod(icycle,convfreq)==0) then
        if (totaltime>mintime.and.abs(lastf-sumf)<0.00000001)  then
          write(0,*) "CONVERGENCE reached. change in F= ", abs(lastf-sumf)
          exit  
        endif
        lastf = sumf
      endif
      if (fing==1.and.halflifeonly.and.sumf>50.0) exit
      if (fing==0.and.halflifeonly.and.sumu>50.0) exit
    endif
    !! ----------- Update estimate of half-life. OBSELETE?
    if (mod(icycle,lifefreq)==0) then
      sumf = getsumf()*100./conc
      sumu = getsumu()*100./conc
      if (fing==1) then
        do while (ihalf < nint(sumf))
          ihalf = ihalf + 1
          if (ihalf>100) exit
          halflife(ihalf) = totaltime
        enddo
        ! if (halflifeonly.and.sumf>50.0) exit
      else
        do while (ihalf < nint(sumu))
          ihalf = ihalf + 1
          if (ihalf>100) exit
          halflife(ihalf) = totaltime
        enddo
        !! -------- Useful for quitting early. Just get half-life
        ! if (halflifeonly.and.sumu>50.0) exit
      endif
    endif
    call cpu_time(xout)
    tmonitor2 = tmonitor2 + (xout - xin)
  enddo LOOP1
  call cpu_time(cpuend)
  y = cpuend - cpustart
  write(*,'("CPU_TIME spent simulating            =",f9.4," seconds.",'//&
          'f9.2,"% calculating dc",f9.2,"% bookkeeping")') &
           y, 100*(tmonitor/y),100*(tmonitor2/y) 
  cpustart = cpuend
  !! --------------------------- all done. ----------------------------
  sumf = getsumc(1)*100./conc
  sumu = getsumc(-1)*100./conc
  sumi = getsumc(0)*100./conc
  sumtot = getsumc(SUMALLFLAG)*100./conc
  write(*,'(a,E12.5e1,3(1x,f9.5),$)') "TIMECOURSE ", totaltime, sumf, sumu , sumi
  !write(*,'(a,E12.5e1,4(1x,f9.5),i12,$)') "TIMECOURSE ", totaltime, sumf, sumu, &
  !                                         sumi, sumtot, nint(maxdc*100000000.)
  !! ----------- final halflife
  if (fing==1) then
    write(*,'(E12.5e1,$)') halflife(nint(sumf/2.))
  else
    write(*,'(E12.5e1,$)') halflife(nint(sumu/2.))
  endif
  !! ----------- final residue intm
  if (intm /= 0) then
    if (intm<0 .or. intm>minter) then
      j = maxloc(inter(:)%conc, mask=(inter(:)%folded>=0),dim=1)
    else
      j = intm
    endif
    write(*,'(i5,f8.3,$)') j, inter(j)%conc*100./conc
  endif
  write(*,*)
  !!---- print out molten intermediates. They should all have zero conc.
  if (verbose.and.fing==1) then
    write(*,'("====== MOLTEN INTERMEDIATES (diagnostic, should be zero) ============")')
    do its=1,mtstate
      if (tstate(its)%cuttype == "m") then  
        write(*,'("TSTATE ",4i5,3f9.6," fsize=",i4,2(1x,1pe10.4e2))') &
            ii,tstate(its)%f,tstate(its)%u1,tstate(its)%u2,&
            inter(tstate(its)%f)%conc,   &
            inter(tstate(its)%u1)%conc,   &
            inter(tstate(its)%u2)%conc,   &
            inter(tstate(its)%f)%nres,   &
            tstate(its)%ku, tstate(its)%kf
        endif
    enddo
    write(*,'("====== END OF MOLTEN INTERMEDIATES ============")')
  endif
  !! ----------- read/write new DAG file with traffic values
  call readDAGfile(inputfile,traffic_out=1)
  !! ----------- write age plot data
  call writeageplot(inputfile)
  !! ----------- write age plot data
  call writeenergyprofile(inputfile,mtstate)
  !! ----------- clean up
  if (allocated(age)) deallocate(age)
  if (allocated(dc)) deallocate(dc)
  call cpu_time(cpuend)
  write(*,'("CPU_TIME spent reporting results     =",f9.4," seconds.")') cpuend - cpustart
  cpustart = cpuend
CONTAINS
  !------------------------------------------------------------------------------!
  subroutine readDAGfile(inputfile, traffic_out)
    !! read and write DAG file
    implicit none
    character(len=*),intent(in) :: inputfile
    integer :: ios,i,j,ii,jj, ounit=99, k, nu,iunit
    integer, optional :: traffic_out
    character(len=MAXLEN) :: aline, bline
    character             :: ach
    real(8) :: sumt
    real    :: discard
    integer :: nres, iseg, ires, iseam
    !!!!!
    iunit = pickunit(10)
    open(iunit,file=inputfile,status='old',action='read',iostat=ios)
    if (ios/=0) stop 'unfoldsim.f90:: ERROR file not found.'
    ounit = pickunit(11)
    if (present(traffic_out)) then
      open(ounit,file=trim(inputfile)//".out",status='replace',action='write',iostat=ios)
      if (ios/=0) stop 'Error opening inputfile'
      !! normalize traffic values
      !! ------------------------------ changed 12-MAY-2009
      !! traffic is "relative unfolding traffic", which is defined
      !! as the fraction of each f that goes through tstate i.
      !! So the output value is often 1.0, when i is the only tstate
      !! and sometimes < 1, when there are multiple tstates
      !! from the same f.
      !!------------------------------------------------------
      do i=1,mtstate
        sumt = 0.
        do j=1,mtstate
          if (tstate(i)%f==tstate(j)%f) sumt = sumt + tstate(j)%traffic
        enddo
        if (sumt==1.00) cycle
        if (sumt<=0.00) cycle
        do j=1,mtstate
          if (tstate(i)%f==tstate(j)%f) tstate(j)%traffic = tstate(j)%traffic/sumt
        enddo
      enddo
      !do i=1,mtstate
      !  tstate(i)%traffic = tstate(i)%traffic*mtstate/sumt
      !enddo
    endif
    !! read file, get mtstate, minter
    mtstate = 0
    !!!!!
    minter = 0
    !!!!!
    Do
      read (iunit, '(a)', iostat = ios) aline
      if(ios/= 0) exit
      if (aline(1:7)=="TSTATE ") mtstate = mtstate + 1
      !!!!!
      if (aline(1:7)=="ISEGMT") then
        minter = minter + 1
        read(iunit, '(a)', iostat = ios) aline
        if (ios/=0) stop 'Premature end. Expecting flags line.'
      endif
      !!!!!
    enddo
    !! allocate memory
    if (.not.allocated(tstate)) then
      allocate(tstate(mtstate),stat=ios)
      if (ios/=0) stop 'Error 1'
    endif
    write(*,*) mtstate,' t-states allocated.'
    !!!!!
    !rewind(iunit)
    !minter = 0
    !Do
      !read (iunit, '(a)', iostat = ios) aline
      !if(ios/= 0) exit
      !if (aline(1:7)=="ISEGMT") then 
      !minter = minter + 1
        !read (iunit, '(a)', iostat = ios) aline
        !if (ios/=0) stop 'Premature end. Expecting flags line.'
      !endif 
    !enddo
    !!!!!
    if (.not.allocated(inter)) then
      allocate(inter(minter),stat=ios)
      if (ios/=0) stop 'Error allocating intermediates' !!!!!Error 2
    endif
    write(*,*) minter,' intermediates allocated.'
    rewind(iunit)
    !! read again, save data
    i = 0
    ires = 0
    nu = 0
    iseg = 0
    aline = " "
    Do
      read(iunit,'(a)',iostat=ios) aline
      if(ios/= 0) exit
      if (aline(1:5)=="ATOM ") then
        ires = ires + 1
        read(aline(31:54),'(3f8.3)',iostat=ios) calpha(1:3,ires)
        if (ios/=0) stop 'BUG in unfoldsim.f90 ATOM lines.'
        if (present(traffic_out)) then
          write(aline(55:60),'(f6.2)') getsumc(SUMALLFLAG,ires=ires)
          write(aline(61:66),'(f6.2)') getsumc(-1,ires=ires)  !! conc unfolded in B-factor column
          write(aline(67:72),'(f6.2)') getsumc(1,ires=ires)  !! conc folded after B-factor column
          write(ounit, '(a)') trim(aline) 
        endif
      elseif (aline(1:7)=="TSTATE ") then
        !! The TSTATE lines are written by geofold.f90, and re-written by this program.
        !! They are read by this program , isegment.f90, and maxTraffic.cpp
        i = i + 1
        iseam = 0
        read(aline(8:),*,iostat=inputstatus) j, &  
          tstate(i)%f,tstate(i)%u1,tstate(i)%u2,tstate(i)%entropy,tstate(i)%cuttype,iseam
        if (inputstatus/=0) then
           write(0,'(a)') trim(aline)
           stop 'Parsing error 1'
        endif
        if (tstate(i)%cuttype=="s") then
          if (tstate(i)%u2/=0) then
            write(0,'("WARNING: cuttype=seam and u2=",i9)') tstate(i)%u2
            tstate(i)%u2 = 0
          endif
          !! Read seam to the end of the line. From geofold.f90. This number is also read by isegment.f90
          !read(aline(8:),*,iostat=inputstatus) j, &  
          !  tstate(i)%f,tstate(i)%u1,tstate(i)%u2,tstate(i)%entropy,tstate(i)%cuttype,iseam
          if (verbose)  write(*,'("SEAM move encountered ",i5,"-->",i5, "seam=",i3)') &
            tstate(i)%f,tstate(i)%u1, iseam
          !! Add this seam to the list
        endif
      if (present(traffic_out)) then 
          !! Add traffic to the end of the line. WARNING: this line is read by maxTraffic.cpp
          j = len_trim(aline)+2
          write(aline(j:),'(f8.4)') tstate(i)%traffic
          write(aline,'(a7,i7,i7,i7,i7,f11.2,3x,a1,i5,f13.8)') "TSTATE ",i,tstate(i)%f,tstate(i)%u1, &
             tstate(i)%u2,tstate(i)%entropy,tstate(i)%cuttype,iseam,tstate(i)%traffic
        write(ounit, '(a)') trim(aline) 
          !TSTATE    2636      1   2774   2763        0.22   p    0
      endif 
      elseif (aline(1:7)=="ISEGMT ") then
        !! The ISEGMT lines are written by geofold.f90, and re-written by this program.
        !! They are read by this program , isegment.f90, and maxTraffic.cpp
        iseg = iseg + 1
        if (.not.associated(inter(iseg)%f)) then
          allocate(inter(iseg)%f)
          nullify(inter(iseg)%f%next)
          inter(iseg)%f%state = 0
          inter(iseg)%f%state = inter(iseg)%sym 
          inter(iseg)%f%axis = 0
          inter(iseg)%f%barrel = 0
        endif
        read(aline(8:),*,iostat=inputstatus) j, k, inter(iseg)%sym, &
          inter(iseg)%solv,inter(iseg)%sas,inter(iseg)%lambda ,inter(iseg)%void, inter(iseg)%hb, &
          discard, inter(iseg)%f%barrel(:)
          !! discard is a placeholder for concentration.
        if (inputstatus/=0) then
          write(*,*) 'iseg=',iseg
          write(*,*) 'size of barrels=',size(inter(iseg)%f%barrel)
          write(*,'(a)') trim(aline)
          stop 'unfoldsim.f90:: Parsing error 2: ISEGMT line'
        endif
      if (present(traffic_out)) then 
          !! This is the line from isegment.f90 that reads the ISEGMT lines. Make sure these sync!
          !isegment.f90:284     read(line1(7:),*) nseg,i,nsym,sas,ntrp,nvoid,nhb,ftype,conc
          aline = " "
          write(aline(1:),'("ISEGMT ",2i5,i3,f12.3,f8.3,i5,i5)') j, k, inter(iseg)%sym, &
            inter(iseg)%sas,inter(iseg)%lambda,inter(iseg)%void, inter(iseg)%hb
          j = len_trim(aline)+2
          write(aline(j:),'(f12.8)') inter(iseg)%conc
          j = len_trim(aline)+2
          write(aline(j:),'(10i4)') inter(iseg)%f%barrel(:)
          j = len_trim(aline)+2
          select case (inter(iseg)%folded)
          case (1)
            aline(j:j) = "F"
          case (0)
            aline(j:j) = "I"
          case (-1)
            aline(j:j) = "U"
          end select 
        write(ounit, '(a)') trim(aline) 
      endif 
        !! The line after ISEGM must be the GeoFold flags. Aline must be len=MAXLEN
        read(iunit,'(a)',iostat=ios) aline
        if (ios/=0) stop 'Premature end. Expecting flags line.'
        inter(iseg)%f%iflag = trim(aline)
        k = 1
        nres = 0
        do while (aline(k:k)/=' ')
          if (aline(k:k)/='.') nres = nres + 1
          k = k + 1
          if (k>MAXLEN) stop 'unfoldsim:: ERROR. Flag line too long.'
        enddo
        inter(iseg)%nres = nres
        if (nres==1) then
          inter(iseg)%conc = uconc
          nu = nu + 1
        else
          inter(iseg)%conc = iconc
        endif
        ach = '.'
        inter(iseg)%residuelist = " "
        inter(iseg)%flags = aline
        k = 1
        do while (aline(k:k)/=' ')
          if (aline(k:k)/='.'.and.ach=='.') then
            write(inter(iseg)%residuelist,'(a,i4,a)') trim(inter(iseg)%residuelist), k, '-'
          elseif (aline(k:k)=='.'.and.ach/='.') then
            write(inter(iseg)%residuelist,'(a,i4,a)') trim(inter(iseg)%residuelist), k-1, ','
          endif
          ach = aline(k:k)
          k = k + 1
        enddo
        if (ach/='.') then
          write(inter(iseg)%residuelist,'(a,i4)') trim(inter(iseg)%residuelist), k-1
        endif
        if (present(traffic_out)) &
          write(ounit, '(a)') trim(aline) // " " // trim(inter(iseg)%residuelist)
      else  !! HEADER, REMARK lines
        if (present(traffic_out)) write(ounit, '(a)') trim(aline) 
      endif
    enddo
    inter(1)%conc = fconc
    !! close
    if (present(traffic_out)) then
      write(*,*) "Traffic and final concentration values sent to output DAG file ",&
        trim(inputfile)//".out"
      close(ounit)
    else
      write(*,'(i9," native state node initialized to conc =",f8.3,"M")') 1, fconc
      write(*,'(i9," leaf nodes initialized to conc =",f8.3,"M")') nu, uconc
    endif
    close(iunit)
!
  end subroutine readDAGfile
  !-----------------------------------------------------------------------------!
  ! real(8) function dSfunction(ctype,nhb,totalhb)
  real(8) function dSfunction(tstate,totalhb)
    type (tstype),intent(in) :: tstate
    !character,intent(in) :: ctype
    integer,intent(in) :: totalhb
    character :: ctype
    integer :: nhb, nu
    real(8) :: x
    !! (p=pivot,h=hinge b=break m=melt)
    dSfunction = 0
    ctype = tstate%cuttype
    nhb = tstate%dhb
    nu = 0
    !! get change in unfolded residues
    if (inter(tstate%f)%folded/=0) then
      if (inter(tstate%u1)%folded==0) nu = nu + inter(tstate%u1)%nres
      if (tstate%u2/=0) then; if (inter(tstate%u2)%folded==0) nu = nu + inter(tstate%u2)%nres; endif
    endif
    select case (ctype)
      case ("p") 
        ! BUG found 21-MAY-2009. entropy used in two places, here and Solidity().
        ! dSfunction = tstate%entropy*pivotpointentropy
        dSfunction = pivotpointentropy + nu*SPERRESID/2.
      case ("h") 
        ! dSfunction = tstate%entropy*hingepointentropy
        dSfunction = hingepointentropy
      case ("b") 
        ! dSfunction = tstate%entropy*breakpointentropy
        dSfunction = breakpointentropy + nu*SPERRESID
      case ("m")
        !! Melting rates are fixed using MELTINGENERGY
        dSfunction = 0.0     
      case ("s")
        dSfunction = seampointentropy
      case default
        write(0,*) 'unfoldsim:: unknown cuttype: ',tstate%cuttype
      dSfunction = 0.0     
    end select
    !! experimental -- try locking TdS to the number of Hbonds broken.
    !if (totalhb==0 .or. nhb<0) then
    !  x = 0
    !else
    !  x= sqrt(real(nhb)/real(totalhb))
    !endif
    !dSfunction = dSfunction*x
  end function dSfunction
  !------------------------------------------------------------------------------!
  real(8) function Solidity(sas,dS,flex)
    ! type (tstype),intent(in) :: tstate
    real(8),intent(in) :: sas, dS
    real(8),intent(in) :: flex
    Real(8) :: dd, meltrange=500
    real(8), parameter :: ovrflo=5.0
    !! -----------------------------------------------------------------------
    !! Solidity should be 1.00 when the surface area is small (approaching zero)
    !! and small when the surface area is large. Any function that has these limits
    !! is fair game. For example Solidity = exp(-dSAS)
    ! if (tstate%dsas < 0.) stop 'BUG!! dSAS < 0. in Solidity()'
    ! Solidity = exp(-tstate%dsas/flex)  !! dsas must be a positive number
        ! dd = (tstate%dsas/flex)
    !! -----------------------------------------------------------------------
    !! Solidity is the fraction of configurational entropy expressed before the
    !! transtion state. It makes sense that it would be less for a tight
    !! pivot, which would have to more or less completely rotate out to be
    !! free. A free pivot, on the other hand, has already gained a good
    !! bit of its entropy before being completely dissociated.
    !! Does that meake sense? If so, then Solidity = entropy (from geofold)
    !! C.B.  Thu Jun 26 21:19:42 EDT 2008
    !! -----------------------------------------------------------------------
    !    Solidity = dS*0.5
    !    return
    !! -----------------------------------------------------------------------
    !    dd = (sas*dS/flex)
    !    if (dd==0.) then
    !        Solidity = 0 
    !    elseif(dd.gt.ovrflo) then
    !        Solidity = 0.5
    !    elseif(dd.lt.0) then
    !        Solidity = 0
    !    else
    !        Solidity = tanh(dd)/2.0
    !    endif
    ! print *, "Solidity is ", Solidity, dd, tstate%dsas, flex
    !! -----------------------------------------------------------------------
    !! Solidity, the fraction of configurational entropy expressed before the 
    !! transition state of unfolding, should be zero, since the degrees of
    !! freedom gained are manifest only when 100% of the interaction is
    !! broken. Here we set Solidity to zero. In order to capture the 
    !! entropic component of the unfolding barrier, we use a new function
    !! dSdd.
    !! -----------------------------------------------------------------------
    Solidity = 0.25
    return
  end function Solidity
  !------------------------------------------------------------------------------!
  subroutine sed(str)
    character(len=*),intent(inout) :: str
    integer :: i
    do i=1,len(str)
      if (str(i:i)=="-") str(i:i) = " "
    enddo
  end subroutine sed
  !------------------------------------------------------------------------------!
  real(8) function getsumu()
    integer :: i,j,k
    real(8) :: x, y
    y = 0
    do k=1,global_nres
      x = 0
      do i=1,minter
        if (inter(i)%flags(k:k)/='.'.and.inter(i)%folded==-1) then
          x = x + inter(i)%conc
        endif
      enddo
      y = y + x
    enddo
    getsumu = y/real(global_nres)
  end function getsumu
  !------------------------------------------------------------------------------!
  !   sumf = sum(inter(:)%conc, dim=1, mask=(inter(:)%folded==1))
  real(8) function getsumf()
    integer :: i,j,k
    real(8) :: x, y
    y = 0
    do k=1,global_nres
      x = 0
      do i=1,minter
        if (inter(i)%flags(k:k)/='.'.and.inter(i)%folded==1) then
          x = x + inter(i)%conc
        endif
      enddo
      y = y + x
    enddo
    y = y/real(global_nres)
    getsumf = y
  end function getsumf
  !------------------------------------------------------------------------------!
  real(8) function getsumc(fstate,ires)
    implicit none
    integer,intent(in) :: fstate
    integer,intent(in),optional :: ires
    integer :: i,j,k
    real(8) :: x, y
    y = 0.00
    if (present(ires)) then
      !! return the amount of ires in state fstate
      k = ires
      if (k<1 .or. k>global_nres) stop 'BUG in calling getsumc. ires'
      x = 0
      do i=1,minter
        if (fstate/=inter(i)%folded .and. fstate/=SUMALLFLAG) cycle
        if (inter(i)%flags(k:k)=='.') cycle
        x = x + inter(i)%conc
      enddo
      getsumc = x
      return
    endif
    do k=1,global_nres
      x = 0
      do i=1,minter
        if (fstate/=inter(i)%folded .and. fstate/=SUMALLFLAG) cycle
        if (inter(i)%flags(k:k)=='.') cycle
        x = x + inter(i)%conc
      enddo
      y = y + x
    enddo
    y = y/real(global_nres)
    getsumc = y
  end function getsumc
  !------------------------------------------------------------------------------!
  integer function pickunit(iunit)
    !!-- select a unit number that is not in use
    implicit none
    integer,intent(in) :: iunit
    logical :: alreadyused
    integer :: ounit
    ounit = iunit
    inquire(unit=ounit,opened=alreadyused)
    do while (alreadyused)
      ounit = ounit + 1
      inquire(unit=ounit,opened=alreadyused)
    enddo
    pickunit = ounit
  end function pickunit
  !!------------------------------------------------------------------
  real(8) function getthetam(dnrg,TdS)
    !!----------------------------
    !! thetam is the fraction of SAS exposed at the transition
    !! state of the elemental unfolding reaction.
    !! If the energy difference between f and u1 + u2 is
    !! positive, then thetam (using the Hammond postulate) is
    !! closer to u1+u2, that is, thetam > 0.5
    !! If the energy is negative, then thetam is closer to
    !! f, that is, thetam < 0.5
    !!----------------------------
    !! thetam = (2*dnrg - TdS) / (2*dnrg), where dnrg is the
    !! interaction energy not including configurational entropy.
    !!----------------------------
    implicit none
    ! real(8),intent(in) :: dnrg,tmscale
    real(8),intent(in) :: dnrg,TdS
    real(8), parameter :: ovrflo=2.0
    real(8) :: dd, gt 
    real(8),parameter :: TMSHIFT=0.5  !!  minimum theta-m
    !! getthetam = tmscale
    !! return
    !dd = dnrg/tmscale
    !if(dd > ovrflo) then
    !   gt = (tanh(ovrflo)+1)/2.
    !elseif(dd < -ovrflo) then
    !   gt = (tanh(-ovrflo)+1)/2.
    !else
    !   gt = (tanh(dd)+1)/2.
    !!endif
    !getthetam = TMSHIFT + (1-TMSHIFT)*gt
    !--------------- modified 31-dec-2009
    gt = (2*dnrg - TdS)/(2*dnrg)
    if (gt > 0.8) gt = 0.8
    if (gt < 0.2) gt = 0.2
    getthetam = gt
    return
  end function getthetam
  !!------------------------------------------------------------------
  subroutine ageplot(t,folding)
    integer, save :: debug = 31
    real(8),intent(in) :: t
    integer,intent(in) :: folding
    type (intype), pointer :: iptr
    real(8) :: x
    integer :: i,j,im,ib,ios
    logical, save :: is_open = .false.
    !!
    if(.not. is_open) then
      debug = pickunit(debug)
      open(debug,file="unfold_debug.log",form='formatted',status='replace',iostat=ios)
      if(ios/= 0) then
        debug = 0
      endif
      is_open = .true.
    else
      open(debug,file="unfold_debug.log",form='formatted',access='append',iostat=ios)
      if(ios /=0) debug = 0
    endif
    do i=1,global_nres
      do j=i+3,global_nres
        ib = 999
        if (age(i,j)<t) cycle   !! age is initialized to HUGE()
        x = 0.
        do im=1,minter
          iptr => inter(im)
          if (iptr%flags(i:i)=='.') cycle
          if (iptr%flags(j:j)=='.') cycle
!          if (any(iptr%f%barrel==ib)) cycle
          if (ignorethiscontact(iptr%f,i,j,ib)) cycle
          x = x + iptr%conc
        enddo
        x = x/conc
        if (folding==1) then
          if (x>0.5) age(i,j)=t
        else
          if (x<0.5) then
            age(i,j)=t
          endif
        endif
      enddo
    enddo
    close(debug)
  end subroutine ageplot
  !!------------------------------------------------------------------
  logical function ignorethiscontact(f,ires,jres, ib) result (inseam)
    type(intermediate),pointer :: f
    integer,intent(in) :: ires, jres
    integer,intent(inout) :: ib
    real :: xin, xout
    inseam = .false.
    inseam = geofold_pivots_queryinseam(f, ires, jres, ib)
  end function ignorethiscontact
  !!------------------------------------------------------------------
  !!------------------------------------------------------------------
  subroutine writepathway(infile)
  !! Modifications:
  !!   29-NOV-2008
  !!     Sums only over imtermediate states. Not folded==1, or folded=-1
  !!     C.Bystroff
  !!--------------------------------------
  !! Record the concentration of each contact Cij at the current
  !! stage of the simulation. The concentration of a contact
  !! is the sum of the concentrations of the folding intermediates
  !! that have Cij . The output file from this routine
  !! can be converted to a postscript image using
  !! pathway2ps.f90 
  !!--------------------------------------
  character(len=*),intent(in) :: infile
  integer :: iunit, ios=0,i,j,im, ib
  real(8) :: sumi, x, sumf
  real(4) :: y, vec(3), dd
  real,parameter :: MAXDD=12.
  type(intype),pointer :: iptr
  !!
  iunit = pickunit(12)
  open(iunit,file=trim(infile)//".path",form='formatted',status='replace',iostat=ios)
  if (ios/=0) stop 'unfoldsim.f90 :: ERROR writing to pathway file. '
  sumi = getsumc(0)
  sumf = getsumc(1)
  do i=1,global_nres
    do j=i+3,global_nres
      vec = calpha(1:3,i) - calpha(1:3,j)
      dd = sqrt(sum(vec*vec))
      if (dd > MAXDD) cycle
      ib = 999
      x = 0.
      do im=1,minter
        iptr => inter(im)
        if (iptr%flags(i:i)=='.') cycle
        if (iptr%flags(j:j)=='.') cycle
        if (any(iptr%f%barrel==ib)) cycle
        if (ignorethiscontact(iptr%f,i,j,ib)) cycle
        x = x + iptr%conc
      enddo
      ! x = x/sumi + TINY(y)
      if (x>0.0) write(iunit,'(2i6,e12.4e2)') i,j,x
    enddo
  enddo
  close(iunit)
  end subroutine writepathway
  !!------------------------------------------------------------------
  subroutine writeageplot(infile)
  character(len=*),intent(in) :: infile
  character(len=200) :: agefile
  integer :: iunit, ios=0,i,j
  !!
  iunit = pickunit(11)
  open(iunit,file=trim(infile)//".age",form='formatted',status='replace',iostat=ios)
  do i=1,global_nres
    do j=i+3,global_nres
      if (age(i,j) < 999999.0) then
        write(iunit,'(2i6,e12.3e2)') i,j,age(i,j)
      endif
    enddo
  enddo
  close(iunit)
  end subroutine writeageplot
  !!------------------------------------------------------------------
  subroutine writeenergyprofile(infile,nt)
    !!------------------------------------
    !! Output the energy profile of the tree with highest traffic.
    !! Each step in the profile (the reaction coordinate) is a recursion
    !! depth in the pathway tree. If more than one branch has the same
    !! recursion depth, the energies are summed. 
    !! The transition state energy and position are plotted for
    !! just the tstate with the highest barrier. Thus, information is lost,
    !! but the plot should provide be qualitatively correct.
    !! Output is suitable for Gnuplot. (energyprofile.csh)
    !! C.Bystroff 26-DEC-2009
    !!------------------------------------
    character(len=*),intent(in) :: infile
    integer,intent(in) :: nt
    integer :: i, intr, idepth, iunit, ios
    real(8),dimension(nt) :: nrg, sas
    integer,dimension(nt) ::  bigf,bigts
    real(8) :: sumnrg,dsas
    idepth = 0
    nrg = 0
    sas = 0
    bigf = 0
    bigts = 0
    intr = 1
    call sumpathtree(idepth,intr,nrg,sas,bigf,bigts,nt)
    iunit = pickunit(11)
    open(iunit,file=trim(infile)//".nrg",form='formatted',status='replace',iostat=ios)
    if (ios/=0) stop 'unfoldsim:: writeenergyprofile : failed to open output file.'
    sumnrg = 0.
    do i=1,nt
      if (bigf(i)==0) exit
      write(iunit,'("ISTATE ",f8.3,2(1pe12.3e2),2i9)') &
        real(i),0.001*sas(i),0.001*sumnrg,bigf(i),inter(bigf(i))%nres
      if (i<nt) then
        dsas = sas(i+1) - sas(i)
        write(iunit,'("TSTATE ",f8.3,2(1pe12.3e2),i9,1x,a1)') &
          real(i)+tstate(bigts(i))%thetam, &
          0.001*(sas(i)+dsas*tstate(bigts(i))%thetam), &
          0.001*(sumnrg+tstate(bigts(i))%dGddu), &
          bigts(i),tstate(bigts(i))%cuttype
      endif
      sumnrg = sumnrg + nrg(i)
    enddo
    close(iunit)
  end subroutine writeenergyprofile
  !!------------------------------------------------------------------
  recursive subroutine sumpathtree(idepth,fint,nrg,sas,bigf,bigts,nt)
    !!---------  Find the next step in the greedy sub-tree of the DAG.
    !!---------  Ignore melting steps
    integer,intent(inout) :: idepth
    integer,intent(in) :: fint, nt
    real(8),intent(inout),dimension(nt) :: nrg,sas
    integer,intent(inout),dimension(nt) :: bigf,bigts
    integer :: u1int, u2int, maxts, maxnres
    real(8) :: maxtr
    !--------- Return if idepth is equal to nt, can't go further
    if (idepth == nt) then
      idepth = idepth - 1
      return
    endif
    !---------
    idepth = idepth + 1
    !! ----- find max traffic tstates of fint
    maxts = 0
    maxtr = 0.
    do i=1,mtstate
      if (tstate(i)%f==fint) then
        if (tstate(i)%cuttype /= "m") then
          if (tstate(i)%traffic > maxtr) then
            maxtr = tstate(i)%traffic
            maxts = i
          endif
        endif
      endif
    enddo
    if (maxts == 0) then
      idepth = idepth - 1
      return
    endif
    !! ----- sum energies for all ts of the same tree depth
    nrg(idepth) = nrg(idepth) + tstate(maxts)%dnrg
    !! ----- sum SAS
    sas(idepth) = sas(idepth) + inter(fint)%sas
    !! ----- keep highest barrier
    if (bigts(idepth)==0 .and. tstate(maxts)%cuttype/="m") bigts(idepth) = maxts
    if (tstate(maxts)%dGddu > tstate(bigts(idepth))%dGddu) then
      if (tstate(maxts)%cuttype/="m") bigts(idepth) = maxts
    endif
    !! ----- keep largest fragment
    if (bigf(idepth)==0) bigf(idepth) = fint
    if (inter(fint)%nres >  inter(bigf(idepth))%nres) bigf(idepth) = fint
    !! ----- go deeper in DAG
    u1int = tstate(maxts)%u1
    if (u1int/=0) call sumpathtree(idepth,u1int,nrg,sas,bigf,bigts,nt)
    u2int = tstate(maxts)%u2
    if (u2int/=0) call sumpathtree(idepth,u2int,nrg,sas,bigf,bigts,nt)
    !! ----- hit bottom, trace back
    idepth = idepth - 1
    return
  end subroutine sumpathtree
  !!------------------------------------------------------------------
  subroutine currenttraffic(trf)
    real(8),intent(out),dimension(4) :: trf
    real(8),dimension(4) :: ttrf
    real(8),save,dimension(4) :: sumtraf=(/0.0,0.0,0.0,0.0/)
    integer :: i
    ttrf = 0.
    do i=1,mtstate
      if (tstate(i)%cuttype=="b") ttrf(1) = ttrf(1) + tstate(i)%traffic
      if (tstate(i)%cuttype=="p") ttrf(2) = ttrf(2) + tstate(i)%traffic
      if (tstate(i)%cuttype=="h") ttrf(3) = ttrf(3) + tstate(i)%traffic
      if (tstate(i)%cuttype=="s") ttrf(4) = ttrf(4) + tstate(i)%traffic
    enddo
    trf = ttrf - sumtraf
    trf = trf/sum(trf)
    sumtraf = ttrf
  end subroutine currenttraffic
  !!------------------------------------------------------------------
  subroutine readparams(paramfile)
    !!-----------------------------------
    !! Paramaters file is a keyworded file with each line containing
    !! KEYWORD value
    !! If the keyword is not one of those listed below, it will be ignored. 
    !! So comments can be freely added in the parameters file.
    !!-----------------------------------
    character(len=*),intent(in) :: paramfile
    integer :: iunit,ios, i
    !------
    iunit = pickunit(12)
    open(iunit,file=paramfile,status='old',form='formatted',iostat=ios)
    if (ios/=0) stop 'unfoldsim.f90:: ERROR param file not found.'
    do
      read(iunit,'(a)',iostat=ios) aline
      if (ios/=0) exit
      if (aline(1:1)=='!') cycle
      if (aline(1:1)=='#') cycle
      if (verbose) write(*,'(a)') trim(aline)
      select case (aline(1:index(aline,' ')-1))
      case ("REMARK")
      case ("TEMPERATURE")
        read(aline(index(aline,' '):),*,iostat=ios) Temperature
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for TEMPERATURE'
        if (verbose) write(*,'(a,f8.3)') "REMARK system temperature set to ",Temperature
      case ("BREAKPOINTENTROPY")
        read(aline(index(aline,' '):),*,iostat=ios) breakpointentropy
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for BREAKPOINTENTROPY'
        if (verbose) write(*,'(a,f8.3)') "REMARK break entropy set to ",breakpointentropy
      case ("HINGEPOINTENTROPY")
        read(aline(index(aline,' '):),*,iostat=ios) hingepointentropy
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for HINGEPOINTENTROPY'
        if (verbose) write(*,'(a,f8.3)') "REMARK hinge entropy set to ",hingepointentropy
      case ("OMEGA")
        read(aline(index(aline,' '):),*,iostat=ios) omega
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for OMEGA'
        if (verbose) write(*,'(a,f10.3)') "REMARK omega set to ",omega
      case ("INTERMEDIATE")
        read(aline(index(aline,' '):),*,iostat=ios) intm
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for INTERMEDIATE'
        if (verbose) write(*,'(a,i5)') "REMARK following intermediate number ",intm
      case ("CONCENTRATION")
        read(aline(index(aline,' '):),*,iostat=ios) conc
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for CONCENTRATION'
        if (verbose) write(*,'(a,f10.5)') "REMARK total protein concentration ",conc
      case ("FOLDING")
        read(aline(index(aline,' '):),*,iostat=ios) fing
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for FOLDING'
        if (verbose) then; if (fing==0) then;write(*,'(a)') "REMARK we are UNFOLDNG";
        else; write(*,'(a)') "REMARK we are FOLDING";endif; endif
      case ("VOIDENTROPY")
        read(aline(index(aline,' '):),*,iostat=ios) spervoid
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for VOIDENTROPY'
        if (verbose) write(*,'(a,f10.5,a)') "REMARK each void counts for ",spervoid," kj/mol/deg"
      case ("SOLIDITY")
        read(aline(index(aline,' '):),*,iostat=ios) flexibility
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for SOLIDITY'
        if (verbose) write(*,'(a,f10.5)') "REMARK solidity parameter= ",flexibility
      case ("HBONDENERGY")
        read(aline(index(aline,' '):),*,iostat=ios) eperhbond
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for HBONDENERGY'
        if (verbose) write(*,'(a,f10.5,a)') "REMARK each Hbond counts for ",eperhbond," kj/mol"
      case ("HAMMONDSCALE")
        read(aline(index(aline,' '):),*,iostat=ios) tmscale
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for HAMMONDSCALE'
        if (verbose) write(*,'(a,f10.5)') "REMARK transition state placement parameter= ",tmscale
      case ("SIDECHAINENTROPY")
        read(aline(index(aline,' '):),*,iostat=ios) scfac
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for SIDECHAINENTROPY'
        if (verbose) write(*,'(a,f10.5)') "REMARK side chain entropy scale factor = ",scfac
      case ("HINGEBARRIER")
        read(aline(index(aline,' '):),*,iostat=ios) hingebarrier
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for HINGEBARRIER'
        if (verbose) write(*,'(a,f10.5)') "REMARK transition state barrier for hinges= ",hingebarrier
      case ("PIVOTBARRIER")
        read(aline(index(aline,' '):),*,iostat=ios) pivotbarrier
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for PIVOTBARRIER'
        if (verbose) write(*,'(a,f10.5)') "REMARK transition state barrier for pivots= ",pivotbarrier
      case ("BREAKBARRIER")
        read(aline(index(aline,' '):),*,iostat=ios) breakbarrier
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for BREAKBARRIER'
        if (verbose) write(*,'(a,f10.5)') "REMARK transition state barrier for breaks= ",breakbarrier
      case ("SEAMBARRIER")
        read(aline(index(aline,' '):),*,iostat=ios) seambarrier
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for SEAMBARRIER'
        if (verbose) write(*,'(a,f10.5)') "REMARK transition state barrier for seams= ",seambarrier
      case ("SPERRESID")
        read(aline(index(aline,' '):),*,iostat=ios) SPERRESID
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for SPERRESID'
      case ("HALFLIFE")
        read(aline(index(aline,' '):),*,iostat=ios) i
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for HALFLIFE'
        if (i==0) then
          halflifeonly = .false.
          if (verbose) write(*,'(a)') "REMARK simulation will go to CONVERGENCE  "
        else
          halflifeonly = .true.
          if (verbose) write(*,'(a)') "REMARK simulation will go to HALFLIFE  "
        endif
      case ("VERBOSE")
        read(aline(index(aline,' '):),*,iostat=ios) i
        if (ios/=0) write(0,*) 'unfoldsim.f90:: WARNING: bad value for VERBOSE'
        verbose = (i>0)
        if (verbose) write(*,*) "VERBOSE ON"
      case ("MAXTIME")
        read(aline(index(aline,' '):),*,iostat=ios) maxtime
        if (ios/=0) STOP 'unfoldsim.f90:: bad value for MAXTIME'
        if (verbose) write(*,'(a,f10.5)') "REMARK maximum simulation time  ",maxtime
      case default
      end select 
    enddo
    close(iunit)
    !! pivotpointentropy must be the average of breakpointentropy and hingepointentropy
    !! since (1) the total entropy is a state function, 
    !! and (2) for every hinge there is one break, and 
    !! (3) a hinge plus a break is equivalent to two pivots
    pivotpointentropy = (breakpointentropy + hingepointentropy)/2.0
  end subroutine readparams
  !------------------------------------------------------------------

 logical function isnan(x)
!====NAN.F90 illustrates what works and what doesn't when
!    detecting a NaN
! Platforms: Windows 9x/Me/NT/2000, AIX 4.3.
! Compilers:I)Compaq Visual Fortran 6.6a with default
!             fpe settings (/fpe:3 /nocheck) and /NOPTIMIZE.
!             (ISNAN is an Elemental Intrinsic Function)
!             Options /fpe:0 /traceback will cause this
!             program to stop with error message,
!          II) AIX XLF90 without optimization.
!              (ISNAN is part of a BOS Runtime C library;
!               thus ISNAN must be declared LOGICAL.)
!
! Author: hdkLESS at SPAM psu dot edu
! Date: March, 2002
!
! Output:
! Y = Plus Infinity
! i=           1  Y= Infinity
!
! Y = Minus Infinity
! i=           2  Y= -Infinity
!
! Y = Minus Zero
! i=           3  Y=  0.0000000E+00
!
! 1) Y is a NaN
! 2) Y is a NaN
! 3) Y is a NaN
! i=           4  Y= NaN
!
! 1) Y is a NaN
! 2) Y is a NaN
! 3) Y is a NaN
! i=           5  Y= NaN
!
! References:
! http://www.psc.edu/general/software/packages/ieee/ieee.html
! http://homepages.borland.com/efg2lab/Mathematics/NaN.htm
!
     implicit none
       ! logical :: ISNAN
       integer :: i
       real, intent(in) :: x
       real, Dimension (6) :: y
       real ::  PInf, MInf, MZero, DivNan
       data PInf/B'01111111100000000000000000000000'/    ! +Infinity
       data MInf/B'11111111100000000000000000000000'/    ! -Infinity
       data MZero/B'10000000000000000000000000000000'/   ! -0
       ! data y(1)/B'01111111100000000000000000000000'/       ! +Infinity
       ! data y(2)/B'11111111100000000000000000000000'/       ! -Infinity
       ! data y(3)/B'10000000000000000000000000000000'/       ! -0
       ! data y(4)/B'01111111100000100000000000000000'/       ! NaN
       ! data y(5)/B'11111111100100010001001010101010'/       ! NaN
       ! DivNan=0
       ! y(6)=DivNan/DivNan
       isnan = .true.
       if (x.eq.PInf) return
       if (x.eq.MInf) return
       isnan = .false.
 end function isnan

!!------------------------------------------------------------------
  subroutine setfolded()
  implicit none
  integer i,j,k,iseg,its
  !!
  do iseg=1,minter
    if (inter(iseg)%sas > FCUTOFF*inter(1)%sas) then
       inter(iseg)%folded = 1
    elseif (inter(iseg)%sas < UCUTOFF) then
       inter(iseg)%folded = -1
    else
      inter(iseg)%folded = 0
    endif
  enddo
  !! Also set as unfolded any intermediate that is to be melted.
  !do its=1,mtstate
  !  if (tstate(its)%cuttype == "m") then
  !    iseg = tstate(its)%f
  !    if (inter(iseg)%folded == 0) inter(iseg)%folded = -1
  !    iseg = tstate(its)%u1
  !    if (inter(iseg)%folded == 0) inter(iseg)%folded = -1
  !    iseg = tstate(its)%u2
  !    if (inter(iseg)%folded == 0) inter(iseg)%folded = -1
  !  endif
  !enddo
  !! Commented out because it will label SS-trapped states as unfolded.
  !! They are not and should not be labeled that way.
  !! Changes in GeoFOLD eliminate large melting segments, so this is mute.
  end subroutine setfolded
!!------------------------------------------------------------------
  real(8) function getSbarrier(tst)
  type (tstype),intent(in) :: tst
  !!-------------------------------
  !! The entropic barrier to folding, unfolding is destinct from
  !! the configurational entropy. It is a function of the 
  !! width of the unfolding/folding pathway. If there are many
  !! ways to break the interactions, then getSbarrier() is less.
  !! If there are few ways to break the interaction, getSbarrier()
  !! is higher. GeoFOLD returns a value (tst%entropy) that is the fraction
  !! of all possible rotations/hinges/break vectors that lead to
  !! loss of interaction with no collisions. This is a number between
  !! 0 and 1. As tst%entropy approaches 1, getSbarrier() approaches zero.
  !! As tst%entropy approaches 0, getSbarrier() approaches infinity.
  !! A function that fits is log(1/tst%entropy)
  !! Also, there is a qualitative difference between different cuttypes.
  !! Breaks represent diffusion in Cartesion space, 3 degrees of freedom.
  !! Pivots represent diffusion in angle space, with the number of angles
  !! being unknown, a function of chain stiffness and orientation.
  !! Hinges are doubly dependent on stiffness and orientation and involve
  !! roughly twice as many angles. 
  !! Seams add no explicit degrees of freedome.
  !! All four moves are different. Therefore, we set a user-dependent
  !! entropic barrier hieght for each cuttype in the limit as tst%entropy
  !! = 1, and add this to the log term. getSbarrier = dSdd(cuttype) + log(1/dS)
  !!-------------------------------
  select case (tst%cuttype)
    case ("h")
      getSbarrier = hingebarrier - log(tst%entropy)
    case ("p")
      getSbarrier = pivotbarrier - log(tst%entropy)
    case ("b")
      getSbarrier = breakbarrier - log(tst%entropy)
    case ("m")
      getSbarrier = 0.
    case ("s")
      getSBarrier = seambarrier 
    case default
      stop 'unfoldsim:: getSbarrier: unknown tstate type.'
  end select
  end function getSbarrier
end program unfoldsim
