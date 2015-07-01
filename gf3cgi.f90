program gf3cgi
!! GeoFOld CGI server program 3
!! Handles expertsettings.php
!! C.Bystroff Thu Nov 20 21:27:50 EST 2008
!!----------------------------------------
  use gfcgi
  implicit none
  character(len=2000) :: aline
  character(len=200) :: email, jobfile,urlname,command,pdbdir,codeline
  character(len=4) :: code
  character(len=24) :: chains, pdbchains
  character(len=50) :: outdir, jobdir
  integer :: pid, ios, i,j
  !!---------
  pdbchains = " "
  codeline = " "
  urlname = " "
  email = " "
  pdbdir = "/ext1/bystrc/server/data/pdb/"
  !! call getenv("PDB",pdbdir)
  outdir="output/"
  jobdir="jobs/"
  !! uses method="POST" 
  !! Write HTML output
  write(*,'(a,//)') 'Content-type: text/html'
  write(*,'(a)') "<html>"
  write(*,'(a)') "<head>"
  ! write(*,'(a,a,a)') '<meta http-equiv="refresh" content="3;url=',trim(urlname),'">'
  write(*,'(a)') "</head>"
  write(*,'(a)') "<body>"
  do
    read(*,*,iostat=ios) aline
    if (ios/=0) exit
    call gfcgi_replace(aline,"params=","")
    call gfcgi_replace(aline,"+"," ")
    call gfcgi_replace(aline,"%0A","")
    ! write(*,'(a)') trim(aline)
    call gfcgi_split(aline,"%0D",6,tag="<br>")
  enddo
  write(*,'(a)') "</body>"
  write(*,'(a)') "</html>"
  !! Write the job file

end program gf3cgi
