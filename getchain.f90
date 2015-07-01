!!-------------------------------------------------------------------
program getchain
!! A simple program to extract just one chain from a PDB file.
!! CB Tue Feb 15 23:20:42 EST 2005
!! Modified. If chain=_ use all chains.
!! Modified. altloc bug fixed.  Tue Jul 13 17:39:11 EDT 2010
!! Modified. Takes multiple chain IDs
!! CB Mon Jul  7 16:27:23 EDT 2008
!! merged using sdiff Fri Jul 30 14:07:37 EDT 2010
!!-------------------------------------------------------------------
character(len=255) :: aline
character :: chain,altloc=" ",chainid
character(len=20) :: chainids
integer :: jarg,iargc,ios,natm,istart,iend,nhet,npro,lastres
logical :: onlyp
!!-------------------------------------------------------------------
jarg = iargc()
istart = -9999
iend = 9999
onlyp = .false.
if (jarg < 1) then
  write(*,*) "Usage: xgetchain chainIDs [start end protein_only] < pdbfile > pdbfile"
  write(*,*) "Selects ATOM lines with character 22 = chain"
  write(*,*) "Use chain = <letter> to use only the chain with that letter."
  write(*,*) "Use chain = . to accept the first chain, regardless of chainID"
  write(*,*) "Use chain = _ to accept the chain with a <blank> character"
  write(*,*) "Use chain = + to accept all protein chains "
  write(*,*) "Use chain = @ to accept all protein chains and all waters (HOH)"
  write(*,*) "Use chain = % to accept all protein chains and all HETATM records"
  write(0,*) 'getchain.f90 Tue Jul 13 17:39:37 EDT 2010'
  stop 
endif
call getarg(1,chainids)
if (chainids(1:1)=="+") onlyp = .true.
chain = chainids(1:1)
if (jarg > 1) then
  call getarg(2,aline)
  read(aline,*,iostat=ios) istart
  if (ios/=0) stop 'Error parsing arg 2'
endif
if (jarg > 2) then
  call getarg(3,aline)
  read(aline,*,iostat=ios) iend
  if (ios/=0) stop 'Error parsing arg 3'
endif
if (jarg > 3) then  !! second way to specify only protein
  call getarg(4,aline)
  if (aline(1:1)=="y") then
    onlyp = .true.
  endif
endif
natm = 0
nhet = 0
npro = 0
chainid = '.'
lastres = -999
do
  read(*,'(a)',iostat=ios) aline
  if (ios/=0) exit
  ! if ((aline(1:4)=="TER ").and.(natm > 0)) exit  !! ignore TER
  if ((aline(1:7)=="ENDMDL ").and.(natm > 0)) exit  !! prevents continuation through a NMR ensemble.
  if (.not.(aline(1:5)=="ATOM ".or.aline(1:6)=="HETATM")) cycle
  if (aline(17:17)/=" ".and.aline(17:17)/=altloc.and.altloc/=" ") cycle  !! ignore alternative locations
  if (aline(14:14)=="H") cycle  !! ignore hydrogens
  read(aline(23:26),*) ires
  if (ires < istart .or. ires > iend) cycle
  !! special case for seleno-methionine MSE
  !! Sun Mar 30 12:39:43 EDT 2008
  if (aline(1:6)=="HETATM".and.aline(18:20)=="MSE") then
    aline(1:6) = "ATOM  "
    aline(18:20) = "MET"
    if (aline(14:14)=="S") aline(13:16)=" SD "
  endif
  select case (chain)
  case ('.')   !!-------accept first chain encountered
    if (chainid=='.') chainid = aline(22:22)
    if (chainid/=aline(22:22)) cycle
  case ('_')   !!-------accept chain with blank character
    if (aline(22:22)/=' ') cycle
  case ('+')   !!-------accept all protein chains
  case ('@')   !!-------accept all protein and water
    !if (aline(1:6)=="HETATM".and.aline(18:20)/="HOH") cycle
  case ('%')   !!-------accept all protein and water and ligands
  case default   !!-------accept only specified chain, ATOM and HETATM
    if (index(trim(chainids),aline(22:22))==0) cycle
  end select
  if (onlyp.and.aline(1:6)=="HETATM") cycle
  altloc = aline(17:17)
  aline(17:17) = " "
  write(*,'(a)') trim(aline)
  if (aline(1:6)=="HETATM") nhet = nhet + 1
  if (aline(1:5)=="ATOM ") npro = npro + 1
  natm = natm + 1
enddo
write(0,'(i9," atoms extracted from PDB file,",i9," protein atoms,",i9," other atoms.")') natm,npro,nhet
end program getchain
