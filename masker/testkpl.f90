program testkpl
  use masker
  use karpluslazaridis
  implicit none

  character (len=809) :: pdb = "ATOM      1  N   ASN A   1      -8.901   4.127  -0.555  1.00  0.00           N&
ATOM      2  CA  ASN A   1      -8.608   3.135  -1.618  1.00  0.00           C&
ATOM      3  C   ASN A   1      -7.117   2.964  -1.897  1.00  0.00           C&
ATOM      4  O   ASN A   1      -6.634   1.849  -1.758  1.00  0.00           O&
ATOM      5  CB  ASN A   1      -9.437   3.396  -2.889  1.00  0.00           C&
ATOM      6  CG  ASN A   1     -10.915   3.130  -2.611  1.00  0.00           C&
ATOM      7  OD1 ASN A   1     -11.269   2.700  -1.524  1.00  0.00           O&
ATOM      8  ND2 ASN A   1     -11.806   3.406  -3.543  1.00  0.00           N&
ATOM      9  H1  ASN A   1      -8.330   3.957   0.261  1.00  0.00           H&
ATOM     10  H2  ASN A   1      -8.740   5.068  -0.889  1.00  0.00           H  "
  character (len=80) :: aline
  real :: gsolv
  integer :: i,count
  real,dimension(3,10) :: xyz = 0
  write(*,*) "pdb: ",pdb

  xyz(:,1) = (/-8.901,4.127,-0.555/)
  xyz(:,2) = (/-8.608,3.135,-1.618/)
  xyz(:,3) = (/-7.117,2.964,-1.897/)
  xyz(:,4) = (/-6.634,1.849,-1.758/)
  xyz(:,5) = (/-9.437,3.396,-2.889/)
  xyz(:,6) = (/-10.915,3.130,-2.611/)
  xyz(:,7) = (/-11.269,2.700,-1.524/)
  xyz(:,8) = (/-11.806,3.406,-3.543/)
  xyz(:,9) = (/-8.330,3.957,0.261/)
  xyz(:,10) = (/-8.740,5.068,-0.889/)


  if(allocated(atype)) deallocate(atype)
  allocate(atype(10))
  atype = 0

  count = 1
  do i = 1,10
    aline = pdb(count:count+77)
    !read(pdb(count:count+77),*)aline
    write(*,*)  "aline: ",aline
    atype(i) = getKplType(aline)
    write(*,'("The type of atom",i2," is ",i2)') i,atype(i)
    count = count + 78
  enddo
  write(*,*) atype
  gsolv = kpl_gsolv(xyz,10,1)
  write(*,*) "GSOLV is ",gsolv

end program testkpl
