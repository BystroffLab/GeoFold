module geofold_pivots
  use vectormath
  use geofold_global
  use geofold_seams   
  private
  !! ------  vectorball data ------
  !! The following is a set of N uniformly-spaced vectors of length 10. These are
  !! used to generate a set of rotation matrices (getmat). the rotation matrices are used to
  !! test for pivots. makemask.f90 (xmakemask) is available within the MASKER package on the Bystroff
  !! lab website.
  !!
  !! Generate vectorball as follows:
  !! xmakemask N vb.pdb
  !! echo 'integer,parameter :: NVB=200 ' > vb.incl   ## 200="N" in comment above
  !! echo 'real,dimension(3,NVB) :: vectorball=reshape((/ &' >> vb.incl
  !! awk '{printf "%7.3f,%7.3f,%7.3f,",$6,$7,$8}' vb.pdb | fold | sed -e "s/\(.*\)/   \1 \&/"  >> vb.incl
  !! echo '     /),(/3,NVB/))' >> vb.incl
  !! Edit the file vb.incl to remove the last comma from the list.
  !!
  ! include "vb1000.incl"
  include "vb.incl"
  ! integer,parameter :: NVB=200  !! "N" in comment above
  ! real,dimension(3,NVB) :: vectorball=reshape((/ &
  ! 0.348, -4.277, -9.032,  1.969,  0.950,  9.758,  2.575,  9.508,  1.721,  2.384,  4.126, -8.792,  &
  ! 1.449,  8.966,  4.185,  7.759, -6.120,  1.530, -2.434, -9.608, -1.330,  2.090,  8.145, -5.412,  &
  ! 0.976,  9.940, -0.491,  7.743, -4.981, -3.904,  1.334, -1.914, -9.724,  3.576,  5.225,  7.740, &
  !-9.028, -3.420,  2.608,  3.026, -5.502,  7.782, -0.155,  5.174, -8.556,  9.093, -4.139,  0.437,  &
  ! 3.176, -9.401,  1.238,  5.728,  8.106, -1.222, -7.911, -5.816, -1.896, -9.646, -2.636, -0.070, &
  !-4.551,  8.884,  0.606,  2.643,  7.454,  6.120, -9.030,  0.537, -4.262,  8.406, -4.312,  3.278,  &
  ! 9.609, -0.690, -2.683, -0.351,  8.718, -4.886, -0.091,  7.821,  6.231, -6.181, -3.229,  7.167,  &
  ! 6.978,  7.097,  0.974,  5.334, -2.720, -8.009,  5.691, -6.353, -5.220,  5.528,  0.130, -8.332, &
  !-2.328,  0.650, -9.703, -2.261, -4.200, -8.789, -6.556,  1.348, -7.430, -7.209,  6.375, -2.718, &
  !-9.668, -1.039, -2.336, -7.036,  1.535,  6.938,  5.404, -8.256,  1.624,  1.709,  9.371, -3.043, &
  !-4.233, -7.134,  5.584, -4.832,  7.995,  3.569, -3.164,  4.030,  8.587, -9.993,  0.291, -0.217,  &
  ! 9.197, -3.290, -2.142,  9.910,  1.202, -0.589,  0.509, -8.341,  5.492,  6.753, -4.072, -6.150,  &
  ! 3.944, -0.499,  9.176,  2.532,  1.122, -9.609, -3.974, -3.019,  8.666, -4.962,  1.660,  8.522, &
  !-7.397,  5.832,  3.356,  4.110, -8.235,  3.910, -7.031, -3.191, -6.355, -1.305, -2.514,  9.591, &
  !-2.724,  7.855,  5.557,  2.971, -3.790, -8.764, -4.607,  0.757, -8.843, -6.079, -7.132, -3.488, &
  !-7.695, -5.369,  3.457, -1.927, -6.806,  7.069,  6.463,  5.265, -5.524,  0.427,  2.778, -9.597, &
  !!-1.246, -1.984, -9.722, -4.782,  5.975, -6.437,  1.971, -6.067, -7.701,  8.842, -1.699,  4.351,  &
  ! 3.636, -1.140, -9.246,  6.891,  4.915,  5.325, -7.928, -3.286,  5.134,  5.879, -1.689,  7.911,  &
  ! 6.940, -5.001,  5.179,  7.429, -2.771,  6.094,  7.249,  1.031, -6.811,  4.204, -8.314, -3.633,  &
  ! 2.398, -9.626, -1.261,  7.672, -0.197,  6.411,  0.568, -6.537,  7.546,  4.391,  5.083, -7.409, &
  !-8.520,  3.813, -3.586,  8.597,  3.259,  3.934, -0.961,  9.013,  4.224, -5.753, -7.363,  3.562, &
  !-2.008,  3.167, -9.270, -8.686,  4.780, -1.305, -3.276, -0.774,  9.417,  8.646,  4.984,  0.645, &
  !-2.769,  1.839,  9.431, -0.778, -9.460,  3.147, -4.807, -3.865, -7.871, -5.802, -5.479, -6.027, &
  !-3.415,  9.230, -1.776, -2.817,  9.325,  2.261,  1.831, -9.081, -3.765,  9.513, -2.279,  2.078, &
  !-9.556,  1.953, -2.206,  3.340, -7.374, -5.870, -7.118, -6.856,  1.526, -2.110,  6.657,  7.158,  &
  ! 6.378, -6.830,  3.561,  1.820, -9.221,  3.414, -8.579,  1.469,  4.923,  6.248,  6.941, -3.576,  &
  ! 7.676,  6.244, -1.444, -1.749,  9.846,  0.042, -8.642, -2.311, -4.470,  8.685,  0.074, -4.956,  &
  ! 3.491, -3.163,  8.821,  9.799,  0.268,  1.975, -1.943, -8.376,  5.105,  6.697, -7.417, -0.361, &
  !-3.689, -1.791, -9.120, -5.653, -0.863,  8.203,  3.952,  8.324,  3.885, -6.258, -5.424,  5.605, &
  !-3.267, -6.006, -7.298,  0.016,  9.802,  1.981,  5.366, -4.188,  7.326,  0.519, -9.940,  0.965, &
  !-2.678,  7.378, -6.197, -6.612, -7.466, -0.736, -3.245, -8.725, -3.653,  6.385, -7.129, -2.900,  &
  ! 7.249, -1.554, -6.711, -7.391,  3.969,  5.442, -8.775,  3.513,  3.263,  4.387,  2.443, -8.648, &
  ! 7.249, -1.554, -6.711, -7.391,  3.969,  5.442, -8.775,  3.513,  3.263,  4.387,  2.443, -8.648, &
  !-0.678, -0.914,  9.935,  0.845, -7.974, -5.975,  4.091,  8.463, -3.412, -0.721, -9.180, -3.899, &
  !-1.769, -5.299,  8.294, -2.310, -9.670,  1.075, -8.542,  5.060,  1.194,  8.093,  4.751, -3.453,  &
  ! 0.958, -4.270,  8.992,  9.292,  1.922, -3.156, -5.130,  7.460, -4.246,  2.133,  6.401, -7.381, &
  !-5.628,  8.087, -1.710,  4.941,  8.586,  1.368, -6.184,  7.585,  2.055,  7.879,  5.409,  2.945, &
  !-1.766, -7.867, -5.915, -4.237,  5.866,  6.902, -0.273,  2.627,  9.645,  8.152, -5.617, -1.413, &
  !-9.030, -3.646, -2.275, -7.589, -4.894, -4.295, -8.636, -5.017,  0.498, -9.753, -1.098,  1.916, &
  !-0.285,  7.169, -6.966,  4.948, -6.386,  5.894, -0.632, -6.332, -7.714,  4.909,  6.430,  5.879,  &
  ! 4.695, -8.759, -1.112,  5.922,  1.006,  7.995, -9.578,  2.859,  0.288, -6.883,  5.311, -4.941, &
  !-4.717, -8.673, -1.592, -4.888, -8.659,  1.063, -4.300, -7.299, -5.314,  4.033,  2.249,  8.870,  &
  ! 1.280,  6.066,  7.846,  5.661,  3.760,  7.336, -7.263,  6.872, -0.123, -0.037, -9.884, -1.516, &
  !-3.556, -8.764,  3.246, -5.703,  3.706,  7.331, -7.951,  2.443, -5.551, -5.841,  6.306,  5.111,  &
  ! 2.701, -7.522,  6.010,  1.447, -1.785,  9.732,  9.213,  3.653, -1.335, -4.326,  3.302, -8.390, &
  !-2.861,  8.661, -4.098,  9.479,  2.772,  1.569, -7.575, -0.970,  6.456, -2.708,  5.420, -7.956, &
  !-8.977, -1.013,  4.288,  1.830,  3.711,  9.104, -4.240, -5.258,  7.374,  0.111,  0.309, -9.995, &
  !-7.808, -0.540, -6.225,  8.523, -2.544, -4.570, -0.830,  9.653, -2.474,  6.198,  3.267, -7.135,  &
  ! 6.068,  7.116,  3.543,  8.040,  2.915, -5.183,  4.393,  6.997, -5.634, -9.608,  1.339,  2.426,  &
  ! 7.467,  2.417,  6.197, -0.452,  4.393,  8.972,  9.891, -1.446, -0.277,  8.982,  0.835,  4.316, &
  !-5.962, -1.381, -7.909, -6.157,  3.857, -6.872,  3.466,  9.339, -0.880,  4.515, -5.220, -7.236/),(/3,NVB/))
  !! ----- end of vectorball data ------
  !! ----- chain IDs may be read from another module eventually 
  real,private,parameter :: pi=3.1415927
  real,parameter :: BALLLEN=10., MINROT=5*pi/180., SCUT=0.05
  real,dimension(3,3,NVB),save :: ballmat=0.
  real,parameter :: DISTCUT=5.0     !! If atoms are closer than this, they are colliding.
  real,parameter :: DMATCUTOFF=12.0     !! ignore atom pairs this far or farther apart
  real,parameter :: KAPPASTEP=5.0       !! angular pivot step in degrees
  real,parameter :: HINGESTEP=5.0       !! 5 degree rotations
  real,dimension(:,:),allocatable :: distancemat
  logical,dimension(:,:),allocatable :: chaincontact
  integer,parameter :: pbuffer=0    !! ignore pbuffer residues around hinge (obselete. keep at 0)
  integer :: pivottail=4  !! leave at least 5 residues
  integer,parameter :: hingeloop=5  !! at least 5 residues in between hingepoints
  integer,parameter :: cutmod=3   !! allow only pivots and hinges at mod(i,cutmod)==0
  !! pallow = ignore collisions for residues within pallow seq separation.
  !!  pallow = 3 allows a hing to be found in a helix in 1sbp. pallow=3 allows
  !! collisions between i -> i+2. What happens if we set pallow=2?
  integer,parameter :: pallow=2   !! ignore collisions with |i-j| < pallow
  integer,parameter ::  maxsym=10
  integer, dimension(maxsym, maxres) :: geofold_symop  !! see geofold_initsymmetry
  integer :: geofold_nsym
  integer,parameter :: MAXSEAM=50
  integer :: seamlist(MAXSEAM), seamcontact(2,50,MAXSEAM)
  !!
  public :: geofold_getnextpivot, geofold_initpivot, geofold_getnextbreak, geofold_getnexthinge, &
            pivottail, geofold_nsym, geofold_symop, getallchains, geofold_cleanup, NVB, &
            geofold_getballmat, geofold_gethingemat,geofold_getballvec,  &
            geofold_pivots_queryinseam
CONTAINS  !! public routines start with geofold_ not geofold_pivot_
  !!------------------------------------------------------------------------------
  logical function geofold_pivots_disulfidepresent(u1,u2,nres)
    !! Ask whether a disulfide connects u1 and u2.
    implicit none
    integer,intent(in) :: nres
    character,dimension(nres),intent(in) :: u1,u2
    integer :: i,j
    geofold_pivots_disulfidepresent = .false.
    do i=1,nres
      if (u1(i)==".") cycle
      if (geofold_ss(i)==0) cycle
      do j=1,nres
        if (u2(j)==".") cycle
        if (geofold_ss(j)==i) then
          geofold_pivots_disulfidepresent = .true.
          return
        endif
      enddo
    enddo
    return
  end function geofold_pivots_disulfidepresent
  !!========================================================================
  logical function geofold_pivots_ssinseg(ires,jres)
    implicit none
    integer,intent(in) :: ires,jres
    integer :: i
    geofold_pivots_ssinseg = .false.
    do i=ires,jres
      if (geofold_ss(i)/=0) then
        geofold_pivots_ssinseg = .true.
        return
      endif
    enddo
    return
  end function geofold_pivots_ssinseg
 !!========================================================================
  subroutine getdistancemat(calpha,nres)
    implicit none
    integer,intent(in) :: nres
    real,dimension(3,nres),intent(in) :: calpha
    integer :: ires,jres
    if (allocated(distancemat)) &
       stop 'getpivots.f90 :: BUG. Call getdistancemat only once.'
    allocate(distancemat(nres,nres))
    distancemat = 0.
    do ires=1,nres
      do jres=1,ires-1
        distancemat(ires,jres) = sqrt(sum((calpha(1:3,ires)-calpha(1:3,jres))**2))
        distancemat(jres,ires) = distancemat(ires,jres)
      enddo
    enddo
    if (verbose) then
      write(*,*) "Done with getdistancemat"
    endif
  end subroutine getdistancemat
  !!------------------------------------------------------------------------------
  subroutine geofold_initpivot(calpha,nres,chainid)
    implicit none
    integer,intent(in) :: nres
    real,dimension(3,nres),intent(in) :: calpha
    character,dimension(nres) :: chainid
    character,dimension(MAXCHAIN) :: uniqchains
    integer :: nchain,clen, ios=0
    !!
    call getallmatrices
    call getdistancemat(calpha,nres)
    !! optionally rescale vectorball here.
    call geofold_initsymmetry(calpha,nres,chainid)
    call getallchains(chainid,uniqchains,nchain,nres)  !! returns chain chars in uniqchains
    ! clen = count(chainid==uniqchains(1))
    ! if (allocated(chaincontact)) deallocate(chaincontact); if (ios/=0) stop 'geofold_initpivot:: error deallocating chaincontact'
    ! allocate(chaincontact(nchain,nchain),stat=ios); if (ios/=0) stop 'geofold_initpivot:: error allocating chaincontact'
    ! call getchaincontact(chainid,uniqchains,nres,clen,nchain)
  end subroutine geofold_initpivot
  !!------------------------------------------------------------------------------
  subroutine getchaincontact(chainid,uniqchains,nres,clen,nchain)
    implicit none
    integer,intent(in) :: nres, clen, nchain
    !real,dimension(3,nres),intent(in) :: calpha
    character,dimension(nres),intent(in) :: chainid
    character,dimension(MAXCHAIN),intent(in) :: uniqchains
    integer :: ires,i,j,ich, jres, jch
    chaincontact(:,:) = .false.
    if (.not.allocated(distancemat)) &
      stop 'geofold.f90:: BUG. call getdistancemat before getchaincontact'
    do ich=1,nchain-1
      do jch=ich+1,nchain
        iresloop: do ires=1,geofold_nres
          if (chainid(ires)/=uniqchains(ich)) cycle
          do jres=1,geofold_nres
            if (chainid(jres)/=uniqchains(jch)) cycle
            if (distancemat(ires,jres)<12.0) then
              chaincontact(ich,jch) = .true.
              chaincontact(jch,jch) = .true.
              exit iresloop
            endif
          enddo
        enddo iresloop
      enddo
    enddo
  end subroutine getchaincontact
  !!------------------------------------------------------------------------------
  subroutine getallmatrices
    implicit none
    real,save :: phi,psi,kappa=KAPPASTEP*pi/180.,x
    real,dimension(3,3) :: mat
    real,dimension(3) :: vec
    integer :: i
    do i=1,NVB
      psi = acos(vectorball(3,i)/BALLLEN)
      x = sqrt(vectorball(1,i)**2 + vectorball(2,i)**2)
      phi = atan2((vectorball(1,i)/x),(vectorball(2,i)/x))
      call getmat_S(phi,psi,kappa,mat)  !! angles in radians
      ballmat(:,:,i) = mat
    enddo
    !! done populating ballmat
  end subroutine getallmatrices
  !!------------------------------------------------------------------------------
  subroutine geofold_getnextpivot(f,calpha,chainid,u1,u2,nres,pivotpoint,entropy,bvec)
    implicit none
    !!----------------------------------------------------
    !! calpha = coordinates of all alpha carbons.  From native state.
    !! chainID = chain identifiers for current intermediate (static)  !! actual PDB chain ID stored separately
    !! u1,u2 = subset of chainID (conforming). Not-in-set set to '.'
    !  NOTE: Associate non-optional args for readability on calling line. i.e.
    !! call geofold_getnextpivot(calpha=nativef,chainid=fchainid,u1=u1chainid,u2=u2chainid,nres=nres,pivotpoint=i,entropy=s)
    !!----------------------------------------------------
    type(intermediate), POINTER :: f
    integer,intent(in) :: nres
    integer,intent(inout) :: pivotpoint
    real,dimension(3,nres),intent(in) :: calpha
    character,dimension(nres),intent(out) :: u1,u2
    character,dimension(nres),intent(in) :: chainid
    character,dimension(nres) :: u3
    integer,intent(out) :: bvec
    logical :: foundit
    character,dimension(MAXCHAIN) :: uniqchains
    !! integer,parameter :: pivottail=5  !! leave at least 5 residues
    !! integer,parameter :: pbuffer=0  !! ignore 1 residues around hinge
    real,intent(out) :: entropy
    logical,dimension(:),allocatable :: assigned
    !!----------------------------------------------------
    ! real,parameter :: scutoff=0.1
    integer :: nchain,ires,jres,ich,mres,ios=0,nn,n1,n2,jch,kch,kres
    logical,dimension(MAXCHAIN) :: goeswithu1=.false.,goeswithu2=.false.
    !!----------------------------------------------------
    entropy = -1    
    bvec = 0
    !! 
    !! To signal that no pivot was found, send back entropy < 0
    !! and pivotpoint > nres.
    !! If F is too short to pivot, then melt it.
    !! If melting, set pivotpoint = pivotpoint + 1
    !! 
    if (.not.allocated(distancemat)) &
      stop 'geofold.f90:: BUG. call getdistancemat before getnextpivot'
    mres = count(chainid/='.')
    if (mres==0) then      !! BUG?
      pivotpoint = nres+1
      entropy = -1
      u1 = '.'
      u2 = '.'
      return !! nothing left to pivot
    endif
    if (mres==1) then !! finished
      pivotpoint = nres+1
      entropy = -1
      u1 = '.'
      u2 = '.'
      return !! nothing left to pivot
    endif
    call getallchains(chainid,uniqchains,nchain,nres)  !! returns chain chars in uniqchains
    if (allocated(assigned)) deallocate(assigned)
    allocate(assigned(nchain),stat=ios) ; if (ios/=0) stop 'geofold_getnextpivot:: error allocating assigned'
    entropy = -1    
    !if (nchain > 1) then
    !  if (allocated(chaincontact)) deallocate(chaincontact); if (ios/=0) stop 'geofold_getnextpivot:: error deallocating chaincontact'
    !  allocate(chaincontact(nchain,nchain),stat=ios); if (ios/=0) stop 'geofold_getnextpivot:: error allocating chaincontact'
    !  call getchaincontact(chainid,uniqchains,nres,clen,nchain)
    !endif
    resloop: do ires=pivotpoint+pivottail,nres-pivottail
      if (.not.mod(ires,cutmod)==0) cycle resloop
      !! ask if this position is in F
      if (chainid(ires)=='.') cycle resloop
      !! ask if this point is at least pivottail from the N-terminius of the current chain
      jres =  ires-pivottail
      !! ...or has a disulfide link between it and the N-terminus
      foundit = geofold_pivots_ssinseg(jres,ires-1)
      if (chainid(ires)/=chainid(jres).and..not.foundit) cycle resloop
      !! ask if this point is at least pivottail from the C-terminius of the current chain
      jres =  ires+pivottail
      !! ...or has a disulfide link between it and the C-terminus
      foundit = geofold_pivots_ssinseg(ires,jres)
      if (chainid(ires)/=chainid(jres).and..not.foundit) cycle resloop
      !!
      !! If here, the chain is long enough to pivot. First find out if ires is a pivot 
      !! in this chain, ignoring all other chains.
      !!
      u1 = '.'
      where (chainid(1:ires-pbuffer)==chainID(ires)) u1(1:ires-pbuffer) = chainID(ires)
      u2 = '.'
      where (chainid(ires+pbuffer+1:nres)==chainID(ires)) u2(ires+pbuffer+1:nres) = chainID(ires)
      !moved to checkpivot
!      if((count(u1(1:nres)/='.')==0).or.(count(u2(1:nres)/='.')==0))then !children are empty
!        pivotpoint = nres + 1
!        entropy = -1
!        u1 = '.'
!        u2 = '.'
!        return
!      endif
      call checkpivot (f,calpha,u1,u2,calpha(1:3,ires),entropy,nres,bvec)
      if (entropy<pcutoff) cycle resloop
      if (nchain == 1 ) exit resloop
      !!
      !! If there are multiple chains, try all combinations, assigning entire chain to either u1 or u2
      !!
      !! REVISED: Assign chain with smallest distance to u1, u2 in greedy fashion., using the contact distances.
      !! Assign chain to whichever u1/u2 has more contacts.
      !! 
      assigned = .false.
      do ich = 1,nchain
        if (chainid(ires)==uniqchains(ich)) exit
      enddo
      assigned(ich) = .true.
      nn = 0
      do while (.not.all(assigned))
        nn = nn + 1
        if (nn > 2*nchain) then
          !  stop 'BUG in getnextpivots. 1'
          !! Remaining chains are disjoint. Disallow the creation
          !! of a disjoint intermediate.
          cycle resloop
        endif
        jchloop: do jch=1,nchain
          if (assigned(jch)) cycle
          !! Check for disulfide linkages
          u3 = '.'
          where (chainid==uniqchains(jch)) u3 = uniqchains(jch)
          if (geofold_pivots_disulfidepresent(u1,u3,nres)) then
            where (chainid==uniqchains(jch)) u1 = uniqchains(jch)
            assigned(jch) = .true.
          elseif (geofold_pivots_disulfidepresent(u2,u3,nres)) then
            where (chainid==uniqchains(jch)) u2 = uniqchains(jch)
            assigned(jch) = .true.
          else
            !! check u1 u2 distances
            n1 = 0
            n2 = 0
            do kres=1,geofold_nres
              if (chainid(kres)/=uniqchains(jch)) cycle
              do jres=1,geofold_nres
                if (distancemat(jres,kres)<DISTCUT) then
                  if (u1(jres)/='.') n1 = n1 + 1
                  if (u2(jres)/='.') n2 = n2 + 1
                endif
              enddo
            enddo
            if (n1==0.and.n2==0) cycle jchloop  !! disjoint chain
            if (n1 > n2) then
              where (chainid==uniqchains(jch)) u1 = uniqchains(jch)
              assigned(jch) = .true.
            else
              where (chainid==uniqchains(jch)) u2 = uniqchains(jch)
              assigned(jch) = .true.
            endif
          endif
        enddo jchloop
      enddo
      !! at this point all chains have been assigned
      call checkpivot (f, calpha,u1,u2,calpha(1:3,ires),entropy,nres,bvec)
      if (entropy>pcutoff) exit resloop   
      !!!!!!!!!!!!!!!!!
      !goeswithu1 = .false.
      !goeswithu2 = .false.
      !!! assign chains to either u1 or u2 or both
      !do ich=1,nchain
      !  if (chainid(ires)==uniqchains(ich)) then
      !    goeswithu1(ich) = .true.
      !    cycle  !! dont check the pivoting chain
      !  endif
      !  !!
      !  !! check whether ich pivots with u1, u2, both or neither
      !  !! first check whether ich goes with u1
      !  !!
      !  u1 = '.'   !! assign u1 to be ich and the N-trminal part of the pivoting chain
      !  where (chainid(1:ires-pbuffer)==chainID(ires)) u1(1:ires-pbuffer) = chainID(ires)
      !  where (chainid==uniqchains(ich)) u1 = chainid
      !  u2 = '.'   !! assign u2 to be the C-terminal part of the pivoting chain
      !  where (chainid(ires+pbuffer+1:nres)==chainID(ires)) u2(ires+pbuffer+1:nres) = chainid(ires)
      !  !! where (chainid/=uniqchains(ich).and.chainid/=chainid(ires)) u2 = chainid
      !  call checkpivot(calpha,u1,u2,calpha(1:3,ires),entropy,nres)
      !  goeswithu1(ich) = (entropy>scutoff) !! chain ich goes with u1
      !  !! now check whether ich goes with u2
      !  u1 = '.'   !! assign u1 to be the N-trminal part of the pivoting chain
      !  where (chainid(1:ires-pbuffer)==chainID(ires)) u1(1:ires-pbuffer) = chainID(ires)
      !  u2 = '.'   !! assign u2 to be the C-terminal part of the pivoting chain and ich
      !  where (chainid(ires+pbuffer+1:nres)==chainID(ires)) u2(ires+pbuffer+1:nres) = chainid(ires)
      !  where (chainid==uniqchains(ich)) u2 = chainid
      !  call checkpivot(calpha,u1,u2,calpha(1:3,ires),entropy,nres)
      !  goeswithu2(ich) = (entropy>scutoff) !! chain ich goes with u2
      !enddo
      !!!
      !!! if any chain cannot go with either u1 or u2, then the pivot is not possible.
      !!! (ires is probably a hinge)
      !!!
      !if (.not.all(goeswithu1(1:nchain).or.goeswithu2(1:nchain))) cycle
      !!!
      !!! now find the first combination of assignments of chains to u1 and u2 that allows a pivot
      !!! This is not necessarilly the best pivot, but oh well.
      !!!
      !u1 = '.'   !! assign u1 to be the N-terminal part of the pivoting chain
      !where (chainid(1:ires-pbuffer)==chainID(ires)) u1(1:ires-pbuffer) = chainID(ires)
      !u2 = '.'   !! assign u2 to be the C-terminal part of the pivoting chain and ich
      !where (chainid(ires+pbuffer+1:nres)==chainID(ires)) u2(ires+pbuffer+1:nres) = chainid(ires)
      !do ich=1,nchain
      !  if (chainid(ires)==uniqchains(ich)) cycle
      !  if (goeswithu1(ich).and..not.goeswithu2(ich)) then
      !    !! assign this chain to u1
      !    where (chainid==uniqchains(ich)) u1 = chainid
      !  elseif (goeswithu2(ich).and..not.goeswithu1(ich)) then
      !    !! assign this chain to u2
      !    where (chainid==uniqchains(ich)) u2 = chainid
      !  elseif (goeswithu2(ich).and.goeswithu1(ich)) then
      !    !! dont assign this chain yet
      !  else
      !      stop 'getpivots.f90:: BUG in chain assignment '  !! see: (.not.all(goeswithu1.or.goeswithu2))
      !  endif
      !enddo
      !!!
      !!! Now all obligatory chain assignments have been made. Check again to see
      !!! If pivoting is allowed and the entropy is still good.
      !!!
      !call checkpivot(calpha,u1,u2,calpha(1:3,ires),entropy,nres)
      !if (entropy<scutoff) cycle  !! if there is no pivot when the unambiguous chains are added, give up.
      !!!
      !!! If there are any chains left unassigned, try to place them.
      !!!
      !do ich=1,nchain
      !  if (chainid(ires)==uniqchains(ich)) cycle
      !  !! Consider only unassigned chains
      !  if (goeswithu2(ich).and.goeswithu1(ich)) then !! assign to u1 first
      !    where (chainid==uniqchains(ich)) u1 = chainid
      !    call checkpivot(calpha,u1,u2,calpha(1:3,ires),entropy,nres)
      !    if (entropy<scutoff) then                   !! If it wont go with u1, try assigning it to u2
      !      where (u1==uniqchains(ich)) u1 = '.'
      !      where (chainid==uniqchains(ich)) u2 = chainid
      !      call checkpivot(calpha,u1,u2,calpha(1:3,ires),entropy,nres)
      !      if (entropy<scutoff) cycle resloop        !! If we cant assign to either u1 or u2, give up.
      !    endif
      !  endif
      !enddo
      !! IF HERE: we found a pivot, skip out and remember where we were. Next time around, we will
      !! increment the pivotpoint.
    enddo resloop
    pivotpoint = ires    
    if (entropy<pcutoff) entropy = -1
    deallocate(assigned)
    return
  end subroutine geofold_getnextpivot
  !------------------------------------------------------------------------------
  subroutine getallchains(chainid,uniqchains,nchain,nres)
    integer,intent(in) :: nres
    character,dimension(nres),intent(in) :: chainid
    character,dimension(MAXCHAIN),intent(out) :: uniqchains
    integer,intent(out) :: nchain
    character :: lastchain
    integer :: i
    !! 
    lastchain = "?"
    nchain = 0
    do i=1,nres
      if (chainid(i)=='.') cycle
      if (chainid(i)==lastchain) cycle
      nchain = nchain + 1
      uniqchains(nchain) = chainid(i)
      lastchain = chainid(i)
    enddo
  end subroutine getallchains
  !------------------------------------------------------------------------------
  subroutine checkbreak (f, calpha,u1,u2,entropy,nres,bvec)
    !! u1 and u2 are set, now count how many ways u1 can be translated from u2
    !! for a minimum distance of 10A without a collision. (1 step)
    implicit none
    type(intermediate), POINTER :: f
    integer,intent(in) :: nres
    character,dimension(nres),intent(in) :: u1,u2
    real,dimension(3,nres),intent(in) :: calpha
    real,intent(out) :: entropy
    integer,intent(out) :: bvec
    real,dimension(3) :: vec, v1
    real :: dd
    integer :: ires,jres,ivec, ivb
    !! Note: all vectors in vectorball are length BALLLEN=10A. 
    !! to avoid unnecessary computations, we assume 10A in the following code
    !! To change the distance of translation, pre-process vectorball to the
    !! length (BALLLEN) desired, within geofold_init().
    if (.not.allocated(distancemat)) &
      stop 'geofold.f90:: BUG. call getdistancemat before checkbreak'
    entropy = 0.
    if (geofold_pivots_disulfidepresent(u1,u2,nres)) return
    bvec = 0
    vec = 0.
    ivb = 0
    vecloop: do ivec=1,NVB
      vec(1:3) = vectorball(1:3,ivec)
      resloop: do ires=1,nres
        if (u1(ires)=='.') cycle resloop
        v1 = calpha(1:3,ires) + vec
        do jres=1,nres
          if (u2(jres)=='.') cycle 
          if (ires==jres) stop 'BUG: overlapping u1, u2 in checkbreak.'
          if (distancemat(ires,jres)>DMATCUTOFF) cycle !! too far apart to consider
		  !! barrel:f%barrel, nBarrel
          if (geofold_pivots_queryinseam (f, ires,jres)) cycle
          if (geofold_pivots_inclosedbarrel(f,ires,jres)) then
            entropy = -1
            return
          endif
          dd  = sqrt(sum((v1(1:3)-calpha(1:3,jres))**2))
          if (dd < DISTCUT.and.dd < distancemat(ires,jres) ) cycle vecloop !! collision
        enddo
      enddo resloop  !! if here, then no collisions were found.
      ivb = ivb + 1
      bvec = ivec
    enddo vecloop
    entropy = real(ivb)/real(NVB)
    !if (verbose) then
    !  write(*,*) "bvec = ",bvec
    !endif
  end subroutine checkbreak
  !------------------------------------------------------------------------------
  !     call checkpivot(calpha,u1,u2,calpha(1:3,ires),entropy,nres)
  subroutine checkpivot (f,calpha,u1,u2,center,entropy,nres,bvec)
    implicit none
    type(intermediate), POINTER :: f
    integer,intent(in) :: nres
    character,dimension(nres),intent(in) :: u1,u2
    real,dimension(3,nres),intent(in) :: calpha
    real,dimension(3),intent(in) :: center
    real,intent(out) :: entropy
    integer,intent(out) :: bvec
    integer :: ires,jres,ivec
    real :: v1(3),v2(3),vec(3),ang,mat(3,3), dd
    real,parameter :: MAXPIVOTANG=30.0 ! atoms must pivot at least this much

    ! real,parameter :: DISTCUT=5.0 ! geofold_global
    !!
    if((count(u1(1:nres)/='.')==0).or.(count(u2(1:nres)/='.')==0))then
      entropy = -1
      return
    endif
    if (.not.allocated(distancemat)) &
      stop 'geofold.f90:: BUG. call getdistancemat before checkpivot'
    entropy = 0.
    !added entropy = -1 line here, attempt to fix disulfide problem
    if (geofold_pivots_disulfidepresent(u1,u2,nres)) then
      entropy = -1
      write (0,*) "disulfide present"
      return
    endif
    vec = 0.
    bvec = 0
    vecloop: do ivec=1,NVB
      mat(1:3,1:3) = ballmat(1:3,1:3,ivec)
      call rotate_S(mat,center,vec)
      vec = center - vec
      resloop: do ires=1,nres
        if (u1(ires)=='.') cycle resloop
        v1 = calpha(1:3,ires)
        ang = 0.
        do while (ang < MAXPIVOTANG)
          call move_S(v1,mat,vec)
          ang = ang + KAPPASTEP
          do jres=1,nres
            if (u2(jres)=='.') cycle
            !write(0,*)"u2(jres)!='.'"
            if (abs(ires-jres)<pallow) cycle 
            !write(0,*)"abs(ires-jres)>=pallow"
            !! omit adjacents (**possible bug if ires,jres in different chains**)
            if (distancemat(ires,jres)>DMATCUTOFF) cycle !! too far apart to consider
            !write(0,*)"distancemat(ires,jres)<=DMATCUTOFF"
            v2 = calpha(1:3,jres)
            dd = sqrt(sum((v1-v2)**2))
            if (dd > DISTCUT.or.dd > distancemat(ires,jres) ) cycle 
            !write(0,*)"dd<=DISTCUT or dd <= distancemat(ires,jres)"
            if (geofold_pivots_queryinseam (f, ires, jres)) cycle
            if (geofold_pivots_inclosedbarrel(f,ires,jres)) then
              entropy = -1
              return
            endif
            !write(0,*)"cycle vecloop"
            cycle vecloop !! count as collision
          enddo
        enddo
      enddo resloop
      entropy = entropy + 1   !! counts of no-collision
      bvec = ivec
    enddo vecloop  !! ivec

    !! return entropy = the fraction of all rotations that are possible.
    entropy = entropy / real(NVB)
  end subroutine checkpivot
!------------------------------------------------------------------------------
!!  checkhinge(calpha,u1,u2,calpha(1:3,ires),calpha(1:3,kres),entropy,nres)
  subroutine checkhinge (f, calpha,u1,u2,hp1,hp2,entropy,nres)
    implicit none
    type(intermediate), POINTER :: f
    integer,intent(in) :: nres
    character,dimension(nres),intent(in) :: u1,u2
    real,dimension(3,nres),intent(in) :: calpha
    real,dimension(3),intent(in) :: hp1,hp2
    real,intent(out) :: entropy
    integer :: ires,jres,ivec
    real :: v1(3),v2(3),vec(3),ang,mat(3,3),chi,mang,dd,mdd
    ! real,parameter :: MAXHINGEANG=30.0 ! atoms must hinge at least this much (geofold_global)
    real,parameter :: HINGEDISPL=5.0 ! minimum c-alpha displacement at max angle
    ! real,parameter :: HINGESTEP=5.0 ! 5 degree rotations (in module)
    ! real,parameter :: DISTCUT=5.0 ! geofold_global
    !!
    if (.not.allocated(distancemat)) &
      stop 'geofold.f90:: BUG. call getdistancemat before checkhinge'
    entropy = 0.
    if (geofold_pivots_disulfidepresent(u1,u2,nres)) then
      entropy = -1
      write(0,*) "disulfide found, hinge"
      return
    endif
    vec = 0.
    !! try hinging both ways.
    chi = HINGESTEP*pi/180.   !! HINGESTEP degree incr
    call getrotS(hp1,hp2,chi,mat,vec)
    mang = MAXHINGEANG
    mdd = 0.
    resloop: do ires=1,nres
      if (u1(ires)=='.') cycle
      v1 = calpha(1:3,ires)
      ang = 0.
      angloop: do while (ang < MAXHINGEANG)
        call move_S(v1,mat,vec)
        ang = ang + HINGESTEP
        do jres=1,nres
          if (u2(jres)=='.') cycle
          if (abs(ires-jres)<2) cycle
          if (distancemat(ires,jres)>DMATCUTOFF) cycle !! too far apart to consider
          v2 = calpha(1:3,jres)
          dd = sqrt(sum((v1-v2)**2))
          if (dd < DISTCUT.and.dd < distancemat(ires,jres) ) exit angloop
        enddo
      enddo angloop
      if (ang < mang) mang = ang
      v2 = calpha(1:3,ires)
      dd = sqrt(sum((v1-v2)**2))
      if (dd > mdd) mdd = dd
    enddo resloop
    entropy = mang
    call getrotS(hp2,hp1,chi,mat,vec)
    mang = MAXHINGEANG
    resloop2: do ires=1,nres
      if (u1(ires)=='.') cycle resloop2
      v1 = calpha(1:3,ires)
      ang = 0.
      angloop2: do while (ang < MAXHINGEANG)
        call move_S(v1,mat,vec)
        ang = ang + HINGESTEP
        do jres=1,nres
          if (u2(jres)=='.') cycle
          if (abs(ires-jres)<pallow) cycle
          if (distancemat(ires,jres)>DMATCUTOFF) cycle !! too far apart to consider
          v2 = calpha(1:3,jres)
          dd = sqrt(sum((v1-v2)**2))
          if (dd > DISTCUT.or.dd > distancemat(ires,jres) ) cycle
          if (geofold_pivots_queryinseam (f, ires, jres)) cycle
          if (geofold_pivots_inclosedbarrel(f,ires,jres)) then
            entropy = -1
            return
          endif
          exit angloop2
        enddo
      enddo angloop2
      if (ang < mang) mang = ang
      v2 = calpha(1:3,ires)
      dd = sqrt(sum((v1-v2)**2))
      if (dd > mdd) mdd = dd
    enddo resloop2
    entropy = entropy + mang
    entropy = entropy/(2*MAXHINGEANG)  !! ranges from 0 - 1.0
    if (mdd < HINGEDISPL) entropy = 0.
  end subroutine checkhinge
  !------------------------------------------------------------------------------
  subroutine geofold_getnexthinge(f,calpha,chainid,u1,u2,nres,hinge1,hinge2,entropy)
    implicit none
    !!----------------------------------------------------
    !! calpha = coordinates of all alpha carbons.  From native state.
    !! chainID = chain identifiers for current intermediate (static)  !! actual PDB chain ID stored separately
    !! u1,u2 = subset of chainID (conforming). Not-in-set set to '.'
    !  NOTE: Associate non-optional args for readability on calling line. i.e.
    !! call geofold_getnexthinge(calpha=nativef,chainid=fchainid,u1=u1chainid,u2=u2chainid,nres=nres,hinge1=22,hinge2=33,entropy=s)
    !!----------------------------------------------------
    type(intermediate), POINTER :: f
    integer,intent(in) :: nres
    integer,intent(inout) :: hinge1,hinge2
    real,dimension(3,nres),intent(in) :: calpha
    character,dimension(nres),intent(out) :: u1,u2
    character,dimension(nres) :: u3
    character,dimension(nres),intent(in) :: chainid
    character,dimension(nres) :: u1save,u2save
    character,dimension(MAXCHAIN) :: uniqchains
    !! integer,parameter :: pivottail=5  !! leave at least 5 residues
    !! integer,parameter :: pbuffer=0  !! ignore 1 residues around hinge, geofold_global
    real,intent(out) :: entropy
    logical :: foundit
    !!----------------------------------------------------
    ! real,parameter :: hcutoff=0.5  !! radians of rotation required
    integer :: nchain,ires,jres,ich,kres,i,j,kk,lres,mres,jch,kch,nn,n1,n2,ios
    logical,dimension(MAXCHAIN) :: goeswithu1=.false.,goeswithu2=.false.
    character :: newchain
    logical,dimension(:),allocatable :: assigned
    !!----------------------------------------------------
    if (.not.allocated(distancemat)) &
      stop 'geofold.f90:: BUG. call getdistancemat before getnexthinge'
    entropy = -1
    if (all(chainid=='.')) return !! nothing left 
    if (count(chainid/='.')<=2*pivottail+hingeloop+1) return !! not enough left to hinge
    call getallchains(chainid,uniqchains,nchain,nres)
    if (allocated(assigned)) deallocate(assigned)
    allocate(assigned(nchain),stat=ios) ; if (ios/=0) stop 'geofold_getnextpivot:: error allocating assigned'
    if (nchain>1) then
      !!
      !! If there is more than one chain, check to see that at least one is long enough to hinge
      !! If it is not, return. Chains should pivot, break or melt instead.
      !!
      i = 0
      do ich=1,nchain
        if (count(chainid==uniqchains(ich))>=2*pivottail+hingeloop+1)  i = i + 1
      enddo
      if (i==0) return  !! All chains are too short to hinge (h). Pivot (p) or Melt (m) instead.
    endif
    ILOOP: do ires=hinge1+pivottail,nres-pivottail
      ! if (.not.mod(ires,cutmod)==0) cycle ILOOP
      !! ask if this position is in F
      if (chainid(ires)=='.') cycle
      !! ask if this point is at least pivottail from the N-terminius of the current chain
      jres =  ires-pivottail
      foundit = geofold_pivots_ssinseg(jres,ires-1)
      if (chainid(ires)/=chainid(jres).and..not.foundit) cycle
      !! ask if this point is at least pivottail from the C-terminius of the current chain
      jres =  ires+pivottail
      !! ...or has a disulfide link between it and the C-terminus
      foundit = geofold_pivots_ssinseg(ires,jres)
      if (chainid(ires)/=chainid(jres).and..not.foundit) cycle
      !! 
      !! Check for a hinge at (ires, kres)
      !! kres must be at least hingeloop residues away and must be
      !! at least pivotttail residues from the end of the chain
      !!
      kk = ires+hingeloop
      if ((hinge2>hinge1).and.(ires==hinge1+pivottail)) kk = hinge2+hingeloop
      KLOOP: do kres=kk,  nres-pivottail
        ! if (.not.mod(kres,cutmod)==0) cycle KLOOP
        !! ask if this position is in F
        if (chainid(kres)=='.') cycle
        !! ask if this point is at least pivottail from the N-terminius of the current chain
        jres =  kres-pivottail
        !! ...or has a disulfide link between it and the N-terminus
        foundit = geofold_pivots_ssinseg(jres,kres-1)
        if (chainid(kres)/=chainid(jres).and..not.foundit) cycle
        !! ask if this point is at least pivottail from the C-terminius of the current chain
        jres =  kres+pivottail
        !! ...or has a disulfide link between it and the C-terminus
        foundit = geofold_pivots_ssinseg(kres,jres)
        if (chainid(kres)/=chainid(jres).and..not.foundit) cycle
        !! 
        !! Are the two hinge points are in the same chain ?
        !! then u1 = inner segment, and
        !! u2 = outer segments split into two chains
        !!
        if (chainid(kres)/=chainid(ires)) cycle KLOOP
        !----else
        !! Hinge axis goes between two chains, chainid(ires) and chainid(kres)
        !! This type of hinge is possible. However, there will probably be
        !! a break or pivot that does the same job and with higher entropy.
        !! For now, disallow between-chain hinges.
        !----endif
        !! u1 is the inner segment
        !! label the inner segment with the current chain ID
        u1 = '.'
        u1(ires+pbuffer:kres-pbuffer) = chainID(ires)
        !! u2 is the outer segment
        !! label the outer segment, N-terminus with the current chain ID
        u2 = '.'
        where (chainid==chainID(ires)) u2 = chainID(ires)
        u2(ires-pbuffer:kres+pbuffer) = '.'
        !! Find a new chain letter for the outer segment, C-terminal segment
        chloop: do i=1,MAXCHAIN
          newchain = chainletters(i)
          !! if (.not.any(newchain==uniqchains(1:nchain))) exit
          do j=1,nchain
            if (newchain==uniqchains(j)) cycle chloop
          enddo
          exit chloop
        enddo chloop
        if (i>=MAXCHAIN) stop 'geofold_getnexthinge:: cant assign an unused chain letter.'
        where (chainid(kres+pbuffer+1:nres)==chainID(ires)) u2(kres+pbuffer+1:nres) = newchain
        !! flags are set. Check whether this u1, u2 is a hinge, ignoring other chains.
        call checkhinge (f, calpha,u1,u2,calpha(1:3,ires),calpha(1:3,kres),entropy,nres)
        if (entropy<hcutoff) then
          cycle KLOOP
        endif
        if (nchain == 1 ) exit KLOOP  !! found hinge in single chain
        !!-------------------------------------------------------------------------
        !!-------------------------------------------------------------------------
        !!
        !! If there are multiple chains, try all combinations, assigning entire chain to either u1 or u2
        !!
        !! REVISED: Assign chain with smallest distance to u1, u2 in greedy fashion., using the contact distances.
        !! Only allow a multichain hinge if there are no chains crossing from u1 to u2.
        assigned = .false.
        do ich = 1,nchain
          if (chainid(ires)==uniqchains(ich)) exit
        enddo
        assigned(ich) = .true.
        nn = 0
        do while (.not.all(assigned))
          nn = nn + 1
          if (nn > 3*nchain) then
            !! stop 'BUG in getnexthinge. 1'
            !! One of the chains is disjoint. Disallow it.
            cycle KLOOP
          endif
          jchloop: do jch=1,nchain
            if (assigned(jch)) cycle
            !! check for disulfide-linked chains
            u3 = '.'
            where (chainid==uniqchains(jch)) u3 = uniqchains(jch)
            if (geofold_pivots_disulfidepresent(u1,u3,nres)) then
              where (chainid==uniqchains(jch)) u1 = uniqchains(jch)
              assigned(jch) = .true.
            elseif (geofold_pivots_disulfidepresent(u2,u3,nres)) then
              where (chainid==uniqchains(jch)) u2 = uniqchains(jch)
              assigned(jch) = .true.
            else
              !! check u1 u2 distances
              n1 = 0
              n2 = 0
              do lres=1,geofold_nres
                if (chainid(lres)/=uniqchains(jch)) cycle
                do mres=1,geofold_nres
                  if (distancemat(mres,lres)<DISTCUT) then
                    if (u1(mres)/='.') n1 = n1 + 1
                    if (u2(mres)/='.') n2 = n2 + 1
                  endif
                enddo
              enddo
              if (n1==0.and.n2==0) cycle jchloop  !! disjoint chain
              if (n1 > n2) then
                where (chainid==uniqchains(jch)) u1 = uniqchains(jch)
                assigned(jch) = .true.
              else
                where (chainid==uniqchains(jch)) u2 = uniqchains(jch)
                assigned(jch) = .true.
              endif
            endif
          enddo jchloop
        enddo
        !! at this point all chains have been assigned
        call checkhinge (f, calpha,u1,u2,calpha(1:3,ires),calpha(1:3,kres),entropy,nres)
        if (entropy>=hcutoff) exit KLOOP
        !!-------------------------------------------------------------------------
        !!-------------------------------------------------------------------------
        !u1save= u1
        !u2save= u2
        !!! 
        !!!! If there are multiple chains, try all combinations
        !!! 
        !goeswithu1 = .false.
        !goeswithu2 = .false.
        !!! assign chains to either u1 or u2 or both
        !do ich=1,nchain
        !  if (chainid(ires)==uniqchains(ich)) then
        !    goeswithu1(ich) = .true.
        !    cycle  !! dont check the pivoting chain
        !  endif
        !  !! check whether ich hinges with u1, u2, both or neither
        !  !! first check whether ich goes with u1
        !  u1 = u1save
        !  where (chainid==uniqchains(ich)) u1 = chainid
        !  u2 = u2save
        !  !! flags are set. Check whether this u1, u2 definintion is a hinge
        !  call checkhinge(calpha,u1,u2,calpha(1:3,ires),calpha(1:3,kres),entropy,nres)
        !  goeswithu1(ich) = (entropy>hcutoff) !! chain ich goes with u1
        !  !! now check whether ich goes with u2
        !  u1 = u1save
        !  u2 = u2save
        !  where (chainid==uniqchains(ich)) u2 = chainid
        !  !! flags are set. Check whether this u1, u2 definintion is a hinge
        !  call checkhinge(calpha,u1,u2,calpha(1:3,ires),calpha(1:3,kres),entropy,nres)
        !  goeswithu2(ich) = (entropy>hcutoff) !! chain ich goes with u2
        !enddo
        !!! if any chain cannot go with either u1 or u2, then the hinge is not possible.
        !if (.not.all(goeswithu1.or.goeswithu2)) cycle
        !!! now find the first combination of assignments of chains to u1 and u2 that allows a hinge
        !!! This is not necessarilly the best hinge
        !u1 = u1save
        !u2 = u1save
        !!! If any chain MUST go with u1 or u2 (not both), then assign them as such.
        !do ich=1,nchain
        !  if (goeswithu1(ich).and..not.goeswithu2(ich)) then
        !    where (chainid==uniqchains(ich)) u1 = chainid
        !  elseif (goeswithu2(ich).and..not.goeswithu1(ich)) then
        !    where (chainid==uniqchains(ich)) u2 = chainid
        !  elseif (goeswithu2(ich).and.goeswithu1(ich)) then
        !    !! do these later
        !  else
        !    !! This should not happen, since we checked.  if (.not.all(goeswithu1.or.goeswithu2))
        !    stop 'getpivots.f90:: BUG in chain assignment -- checkhinge'
        !  endif
        !enddo
        !call checkhinge(calpha,u1,u2,calpha(1:3,ires),calpha(1:3,kres),entropy,nres)
        !if (entropy<hcutoff) cycle  !! if there is no pivot when the unambiguous chains are added, give up.
        !u1save= u1
        !u2save= u2
        !!! Try to place the remaining chains
        !do ich=1,nchain
        !  if (goeswithu2(ich).and.goeswithu1(ich)) then ! assign to u1 first
        !    where (chainid==uniqchains(ich)) u1 = chainid
        !    call checkhinge(calpha,u1,u2,calpha(1:3,ires),calpha(1:3,kres),entropy,nres)
        !    if (entropy<hcutoff) then  ! No good. Try assigning it to u2 instead
        !      where (u1==uniqchains(ich)) u1 = '.'
        !      where (chainid==uniqchains(ich)) u2 = chainid
        !      call checkhinge(calpha,u1,u2,calpha(1:3,ires),calpha(1:3,kres),entropy,nres)
        !      if (entropy<hcutoff) cycle KLOOP  ! cant assign this chain to either u1 or u2
        !    endif
        !  endif
        !enddo
      enddo KLOOP
      !! IF HERE: we found a hinge, skip out 
      if (entropy>=hcutoff) exit ILOOP
    enddo ILOOP
    if (entropy>=hcutoff) then
      hinge1 = ires
      hinge2 = kres
    else
      hinge1 = nres
      hinge2 = nres
      entropy = -1
    endif
    deallocate(assigned)
  end subroutine geofold_getnexthinge
  !!========================================================================
  ! call geofold_getnextbreak(calpha=allcoords,chainid=f%iflag,u1=u1%iflag,u2=u2%iflag, &
  !                           nres=ires,breakpoint=lastpos,entropy=entropy)
  subroutine geofold_getnextbreak(f,calpha,chainid,u1,u2,nres,breakpoint,entropy,bvec)
    !! Starting from breakpoint, check each remaining chain
    type(intermediate), POINTER :: f
    integer,intent(in) :: nres
    integer,intent(inout) :: breakpoint
    integer,intent(out) :: bvec
    real,dimension(3,nres),intent(in) :: calpha
    character,dimension(nres),intent(in) :: chainid
    character,dimension(nres),intent(out) :: u1,u2
    integer :: ires,ibits,i,j,k,kk
    character,dimension(MAXCHAIN) :: uniqchains
    ! real,parameter :: bcutoff=0.2
    integer,save :: nchain
    !!
    entropy = -1
    bvec = 0
    call getallchains(chainid,uniqchains,nchain,nres)  !! returns chain chars in uniqchains
    if (breakpoint==0) then
      if (nchain > 32) then
        write(0,*) 'geofold_getnextbreak:: WARNING! there are more than 32 chains in chainid'
        nchain = 32
      endif
    endif
    if (nchain<=1) then
      breakpoint = 1
      entropy = -1
      return !! no breaks possible if only one chain
    endif
    kk = 2**(nchain-2)
    if (breakpoint>=kk) then
      entropy = -1
      return
    endif
    !! check all 2^(n-2) possible breaks
    do i=breakpoint+1,kk
      !! set u1 according to the bits set in i, the counter
      u1 = '.'
      u2 = '.'
      do ibits=0,nchain-1
        if (btest(i,ibits)) then
          where (chainid==uniqchains(ibits+1)) u1 = chainid
        else
          where (chainid==uniqchains(ibits+1)) u2 = chainid
        endif
      enddo
      call checkbreak (f, calpha,u1,u2,entropy,nres,bvec)
      if (entropy > bcutoff) exit
    enddo
    if (entropy < bcutoff) then
      entropy = -1
    else
      ! write(*,*)"Found a break with entropy = ",entropy
    endif
    breakpoint = i
    return
  end subroutine geofold_getnextbreak
!!====================================================================================
  subroutine geofold_initsymmetry(coords,nres,chainid)
    implicit none
    integer,intent(in) :: nres
    real,dimension(3,nres),intent(in) :: coords
    character,dimension(nres),intent(in) :: chainid
    character,dimension(MAXCHAIN)        :: uniqchains
    integer :: nchain, clen, i,j,k,ires,jres, pair(MAXCHAIN), nsym, ii,jj,ij,ik, isym, bestj
    real :: err, rmsd, mat(3,3),vec(3), mat2(3,3), vec2(3), dd, bestdd
    real,parameter :: symdistcut=2.0  !! accuracy of symmetry
    !!
    nsym = 0
    geofold_nsym = nsym
    call getallchains(chainid,uniqchains,nchain,nres)  !! returns chain chars in uniqchains
    clen = count(chainid==uniqchains(1))
    ! write(*,*) nchain," chains. Length=", clen
    if (nchain<=1) return
    !! assume all chains are identical.
    do i=1,nchain
      k = count(chainid==uniqchains(i))
      if (verbose) then
        write(*,*) "Initsymmetry:  Chain ",i,"  length=",k
      endif
      if (i>1) then
        if (k/=clen) then
          ! write(*,*) "Chains are not equal length. Symmetry not allowed."
          write(*,*) "No symmetry."
          return
        endif
      endif
      clen = k
    enddo
    nsym = 0
    do i=2,nchain
      !! pair chain 1 with chain i
      pair = 0
      j = (i-1)*clen+1
      k = i*clen
      call ovrlap(coords(1:3,1:clen),coords(1:3,j:k),err,clen,mat,vec)     !! vectormath.f90
      ! write(*,'("Superimposing chain 1 on chain ",i2," RMSD=",f8.3)') i,err
      if (err > 3.0) cycle
      pair(1) = i
      do ii=2,nchain
        if (pair(ii)/=0) cycle
        j = (ii-1)*clen+1
        k = ii*clen
        bestdd = 999.
        bestj = 0
        do jj=1,nchain
          if (ii==jj) cycle
          ij = (jj-1)*clen+1
          ik = jj*clen
          call ovrlap(coords(1:3,j:k),coords(1:3,ij:ik),err,clen,mat2,vec2)     !! vectormath.f90
          if (err > 3.0) write(0,*) "Warning! High RMSD aligning ",ii,jj,err
          dd = matrixdiff(mat,vec,mat2,vec2)
          if (dd<symdistcut) then
             bestdd = dd
             pair(ii) = jj
             bestj = jj
             exit
          else
            if (dd<bestdd) then
               bestdd = dd
               bestj = jj
            endif
          endif
        enddo
        if (pair(ii)==0)  then
           write(*,'("Symmetry operator rotating chain 1 to chain: ",i5, '//&
               '" doesnt apply to chain ",i3, " best dev=",f9.3," for chain",i3)') i, ii, bestdd, bestj
        endif
      enddo
      if (all(pair(1:nchain)/=0)) then
         nsym = nsym + 1
         write(*,'("Symmetry operator found: ",i5)') nsym
         do ii=1,nchain
            write(*,'("SYM ",i5," chain ",i3," ",a1," chain ",i3," ",a1)') nsym, ii, uniqchains(ii), pair(ii), uniqchains(pair(ii))
            do ires=1,clen
               geofold_symop(nsym,(ii-1)*clen+ires) = (pair(ii)-1)*clen+ires
            enddo
         enddo
      else
         write(*,'("Symmetry operator NOT found rotating chain 1 to chain: ",i5,f9.3)') i, dd
      endif
    enddo
    geofold_nsym = nsym
  end subroutine geofold_initsymmetry
!!========================================================================
   real function matrixdiff(m1,v1,m2,v2)
     !! what is the absolute distance after applying m1, v1 and the inverse of
     !! m2, v2 
     real,intent(in) :: m1(3,3),v1(3),m2(3,3),v2(3)
     real :: t1(3),t2(3), xx,yy,zz
     t1 = (/10.,0.,0./)
     t2 = t1
     call move_S(t1,m1,v1)
     call move_S(t2,m2,v2)
     xx = sqrt(sum( (t1-t2)**2 ))
     t1 = (/0.,10.,0./)
     t2 = t1
     call move_S(t1,m1,v1)
     call move_S(t2,m2,v2)
     yy = sqrt(sum( (t1-t2)**2 ))
     t1 = (/0.,0.,10./)
     t2 = t1
     call move_S(t1,m1,v1)
     call move_S(t2,m2,v2)
     zz = sqrt(sum( (t1-t2)**2 ))
     matrixdiff = max(xx,yy,zz)
  end function matrixdiff
!!========================================================================
  subroutine geofold_getballmat(bvec,mat)
    integer,intent(in) :: bvec
    real,dimension(3,3),intent(out) :: mat
    mat = ballmat(:,:,bvec)
  end subroutine geofold_getballmat
  !----------------------------------------------------------------------------------
  subroutine geofold_gethingemat(ij,mat,vec)
    integer,intent(in) :: ij
    real,dimension(3) :: hp1,hp2
    real,dimension(3,3),intent(out) :: mat
    real,dimension(3),intent(out) :: vec
    real :: chi
    integer :: ires,jres
    chi = HINGESTEP*pi/180.   !! HINGESTEP degree incr
    ires = ij/10000       !! tricky business
    jres = mod(ij,10000)  !! tricky business
    hp1 = allcoords(1:3,ires)
    hp2 = allcoords(1:3,jres)
    call getrotS(hp1,hp2,chi,mat,vec)
  end subroutine geofold_gethingemat
  !----------------------------------------------------------------------------------
  subroutine geofold_getballvec(bvec,vec)
    integer,intent(in) :: bvec !! index to vectorball()
    real,dimension(3),intent(out) :: vec
    vec = vectorball(:,bvec)
  end subroutine geofold_getballvec
  !----------------------------------------------------------------------------------
  logical function geofold_pivots_queryinseam(f,ires,jres,ib,seamchar) result (inseam)
    type(intermediate),pointer :: f
    integer,intent(in) :: ires, jres
    integer,intent(out),optional :: ib
    character,intent(in), optional :: seamchar
	  type (seam_type),pointer    :: aseam
    character :: seamch="A"
    integer  :: nb, iseam, fb, i,j
	
	  nb = size(barrels_array)
    if (present(seamchar)) seamch=seamchar
    inseam= .false.
	  do i=1, nb
      fb = f%barrel(i) 
	    if (fb == 0) cycle !! if fb==0, seam fb is still closed
      aseam => barrels_array(i)%seams(fb)
!      if ( ((aseam%u1flag(ires)==seamch).and. &
!            (aseam%u2flag(jres)==seamch) ).or. &
!           ((aseam%u2flag(ires)==seamch).and. &
!            (aseam%u1flag(jres)==seamch) )) then
!        inseam = .true.
!      endif
!      if ( ((aseam%u1flag(ires) /= '.') .and. &
!          (aseam%u1flag(jres) == '.')) .or. &
!          ((aseam%u1flag(ires) == '.') .and. &
!          (aseam%u1flag(jres) /= '.')) .or. &
!          ((aseam%u2flag(ires) /= '.') .and. &
!          (aseam%u2flag(jres) == '.')) .or. &
!          ((aseam%u2flag(ires) == '.') .and. &
!          (aseam%u2flag(jres) /= '.')) ) then
!          inseam = .true.
!       endif 
      if( ((aseam%u1flag(ires) /= ".") .and. &
           (aseam%u2flag(jres) /= ".")) .or. &
          ((aseam%u2flag(ires) /= ".") .and. &
           (aseam%u1flag(jres) /= "."))) then
        inseam = .true.
      endif
    enddo
    if (present(ib)) ib = fb
    return
  endfunction geofold_pivots_queryinseam
  !----------------------------------------------------------------------------------
  logical function geofold_pivots_inclosedbarrel(f,ires,jres) result (inbarrel)
    type(intermediate),pointer :: f
    integer, intent(in) :: ires, jres
    type (seam_type), pointer :: aseam
    integer :: nb, iseam, fb, i, j
    
    nb = size(barrels_array)
    inbarrel = .false.
    loop1: do i = 1, nb
      if(f%barrel(i) /= 0) cycle
      do fb = 1, barrels_array(i)%nSeams
        aseam => barrels_array(i)%seams(fb)
        if(aseam%segments(1)<ires .and. ires<aseam%segments(2)) then
          inbarrel = .true.
        elseif(aseam%segments(1)<jres .and. jres<aseam%segments(2)) then
          inbarrel = .true.
        elseif(aseam%segments(3)<ires .and. ires<aseam%segments(4)) then
          inbarrel = .true.
        elseif(aseam%segments(3)<jres .and. jres<aseam%segments(4)) then
          inbarrel = .true.
        endif
        if(inbarrel) exit loop1
      enddo
    enddo loop1
    return
  endfunction geofold_pivots_inclosedbarrel  
  !----------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------
!----------------------------------------------------------------------------
! Return true/false if the residues (ires,jres) are part of a seam
! Uses the global structure of seams: "barrels_array"
!----------------------------------------------------------------------------
function geofold_pivots_inseam (f, ires, jres ) result (inseam)
  !!!!! SUPERCEDED by geofold_pivots_queryinseam()
	implicit none
	  type(intermediate), POINTER :: f
	  integer, intent (in)        :: ires, jres
	  logical                     :: inseam
	  integer                     :: i, j, nBarrels, nSeams, x1, x2, y1, y2
	  type (barrel_type)          :: barrel
	  type (seam_type)            :: seam
	  integer                     :: nbuttons, allButtons(200)  ! maximum 200 buttons
	
	inseam = .false.
	!call getAllButtons (barrels_array, allButtons)


	nBarrels = size (barrels_array)
	do i=1, nBarrels
		if (f%barrel(i) == 0 ) cycle !! barrel closed, skip.
        !! diagnostic
        !  write(*,*) "Open barrel found i=",i," iseam=",f%barrel(i)  !! 60!!
	    call getAllButtons (barrels_array(i), f%barrel(i), allButtons, nbuttons)
		j = f%barrel(i) !! seam j is open in barrel i, if j < 0, side1, else side2
		seam = barrels_array(i)%seams(abs(f%barrel(i)))
		y1=seam%segments(1)  ! ini Beta 1
		y2=seam%segments(2)  ! end Beta 1
		x1=seam%segments(3)  ! ini Beta 2
		x2=seam%segments(4)  ! end Beta 2
		if (y1<=ires.and.ires<=y2) then ! ires in Beta1
			if (x1<=jres.and.jres<= x2) then  !jres in B2
				inseam = .true.
				return
			endif
			if (any (allButtons(1:nbuttons) == jres)) then ! jres is a button
				if (j < 0) then  ! broken on beta1 side
					inseam = .true.
					return
				endif
			endif
		elseif (x1<=ires.and.ires<=x2) then ! ires in Beta2
			if (y1<=jres.and.jres<= y2) then  !jres in B1
				inseam = .true.
				return
			endif
			if (any (allButtons(1:nbuttons) == jres)) then ! jres is a button
				if (j > 0) then  ! broken on beta2 side
					inseam = .true.
					return
				endif
			endif
		elseif (any (allButtons(1:nbuttons) == ires)) then ! ires is a button
			if (y1<=jres.and.jres<= y2) then  !jres in B1
				if (j < 0) then  ! broken on beta1 side
				  inseam = .true.
				  return
				endif
			endif
			if (x1<=jres.and.jres<= x2) then  !jres in B2
				if (j > 0) then  ! broken on beta2 side
				  inseam = .true.
				  return
				endif
			endif
		endif
	enddo
endfunction

!!========================================================================
  subroutine geofold_cleanup()
    if (allocated(distancemat)) &
        deallocate(distancemat)
    if (allocated(barrels_array)) &
        deallocate(barrels_array)
  end subroutine geofold_cleanup
!!========================================================================
end module geofold_pivots
