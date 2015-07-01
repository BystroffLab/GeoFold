module gfcgi
!! GeoFOld CGI module for gf1cgi and gf1cgi
!! C.Bystroff Fri Nov 21 16:38:37 EST 2008
!!--------------------------------------------
  implicit none
  private
  character(len=200),parameter :: tmpdir="/bach1/home/bystrc/server/geofold/tmp"
  character(len=200),parameter :: curlpdb = "/bach1/home/bystrc/bin/curlpdb "
  public :: gfcgi_parseit, gfcgi_errormsg, gfcgi_findpdbfile, &
            gfcgi_replace, gfcgi_split, gfcgi_getpdbchains, gfcgi_stdinsave, &
            gfcgi_replacefield, gfcgi_untaint
  !!---------
CONTAINS
  !---------------------------------------------
  subroutine gfcgi_stdinsave(iunit,pid)
    implicit none
    integer,intent(out) :: iunit
    integer,optional :: pid
    logical :: isopen=.false.
    integer :: ios, ipid
    character(len=1000) :: aline, tmpfile
    !!
    iunit = 21
    inquire(unit=iunit,opened=isopen)
    do while (isopen)
      iunit = iunit + 1
      inquire(unit=iunit,opened=isopen)
    enddo
    if (present(pid)) then
      ipid = pid
    else
     ipid = getpid()
    endif
    write(tmpfile,'(a,i8.8,a)') trim(tmpdir)//"/",pid,".stdin"
    open(iunit,file=trim(tmpfile),status='replace',form='formatted',iostat=ios)
    if (ios/=0) stop 'gfcgi:: stdinsave: error opening output file. Permissions?'
    do
      read(*,'(a)',iostat=ios) aline
      if (ios/=0) exit
      write(iunit,'(a)') trim(aline)
    enddo
    rewind(iunit)  !! ready to read... just close it when you're done.
  end subroutine gfcgi_stdinsave
  !---------------------------------------------
  subroutine gfcgi_parseit(aline,atag,result)
    !! Find the string between atag and "=" and "&" in aline.
    !! Return it in result.
    !! If there are more than one instance of atag, concatenate
    !! the strings in result.
    !!----
    character(len=*),intent(in) :: aline, atag
    character(len=*),intent(out) :: result
    integer :: i,j,ich, k, n
    ich = 1
    result = " "
    do
      i = index(aline(ich:),atag)
      !! --- if tag not found, return empty result
      if (i==0) exit
      j = i + len(atag)+ich
      k = index(aline(j:),'&')
      result = trim(result)//aline(j:j+k-2)
      ich = j + k
    enddo
  end subroutine gfcgi_parseit
  !---------------------------------------------
  subroutine gfcgi_errormsg(atag,entry)
    character(len=*),intent(in) :: atag,entry
    write(*,'(a,//)') 'Content-type: text/html'
    write(*,'(a)') '<meta http-equiv="refresh" content="30;url=server.php">'
    write(*,'(a)') "<html>"
    write(*,'(a)') "<body>"
    write(*,'(a,a)') '<font color="#FF0000">ERROR! Please enter a valid </font>',trim(atag)
    write(*,'(a,a,a)') 'You entered "',trim(entry),'" </font>'
    write(*,'(a)') '<h2><A HREF="server.php">try again</A></h2>'
    write(*,'(a)') "</body>"
    write(*,'(a)') "</html>"
    stop
  end subroutine gfcgi_errormsg
  !---------------------------------------------
  subroutine gfcgi_findpdbfile(code,pdbdir,pdbchains,biolunit)
    character(len=*),intent(in) :: pdbdir, code, biolunit
    character(len=*),intent(out) :: pdbchains
    character(len=200) :: pdbfile, aline, command
    character :: cnchar
    logical :: isthere
    integer :: ios, ich
    ! curlpdb = "/home/bystrc/bin/curlpdb "
    ich = 0
    !!
    ! write(*,'(a)') "Looking locally...<br>"
    pdbfile = trim(pdbdir)//trim(code)//".pdb"
    inquire(file=pdbfile,exist=isthere)
    if (isthere) then
      ! write(*,'(a)') "Found. <br>"
    else
      ! write(*,'(a)') "Not found. <br>Looking offsite (www.pdb.org)...<br>"
      if (biolunit(1:8)=="biolunit") then
        command = trim(curlpdb)//" "//trim(code)//" 1 "
        ! write(*,'(a)') trim(command)
        call system(command)
      elseif (biolunit(1:6)=="select") then
        !! NOTE. curlpdb puts the file in $PDB. The setting of $PDB
        !! in this program must match the setting in curlpdb
        command = trim(curlpdb)//" "//trim(code)
        call system(command)
      else    !! if chain is any other string, dont send a job to curlpdb
      endif
    endif
    !! see if PDB file landed
    inquire(file=pdbfile,exist=isthere)
    if (isthere) then
      call gfcgi_getpdbchains(pdbfile,pdbchains)
    else
      write(*,'(a)') "<p>gfcgi.f90:: gfcgi_findpdbfile: PDB file not found. Check your input. <br>"
      pdbchains = " "
    endif
  end subroutine gfcgi_findpdbfile
  !---------------------------------------------
    subroutine gfcgi_getpdbchains(pdbfile,pdbchains)
      character(len=*),intent(in) :: pdbfile
      character(len=*),intent(out) :: pdbchains
      integer :: ios, ich
      character(len=200) :: aline
      character :: cnchar
      !!
      open(11,file=pdbfile,status='old',form='formatted',iostat=ios)
      if (ios/=0) then
        pdbchains = " "
        return
        ! stop "gfcgi::getpdbchains: A file was found but could not be opened! BUG?"
      else
        ! write(*,'(a)') "... found. <br>"
        ich = 0
        do 
          read(11,'(a)',iostat=ios) aline
          if (ios/=0) exit
          !if (aline(1:6)=="HEADER") then
          !  write(*,'(a)') trim(aline)//"<br>"
          !endif
          if (aline(1:5)/="ATOM ") cycle
          cnchar = aline(22:22)
          if (cnchar==" ") cnchar = "_"
          if (index(pdbchains,cnchar)==0) then
            ich = ich + 1
            pdbchains(ich:ich) = cnchar
          endif
        enddo
        ! write(*,'(i9,a)') ich," chains found.<br>"
        if (ich==0) pdbchains = " "
      endif
      close(11)
    end subroutine gfcgi_getpdbchains
  !---------------------------------------------
  subroutine gfcgi_replace(str,srch,rplc)
    !! Replace all occurrences of srch in str with rplc
    !! change nothing else
    character(len=*),intent(inout) :: str
    character(len=*),intent(in) :: srch,rplc
    integer :: i,j,k
    i = 1
    k = len(srch)
    do
      i = index(str(1:),srch)
      if (i == 0) exit
      str = str(1:i-1)//rplc//str(i+k:)
    enddo
  end subroutine gfcgi_replace
  !---------------------------------------------
  subroutine gfcgi_untaint(str)
    !! Remove all occurrences of unallowed characters
    !! Replace spaces with underscores
    character(len=*),intent(inout) :: str
    character(len=200) :: rplc
    integer :: i,j,k
    rplc = " "
    j = 0
    do i=1,len_trim(str)
      select case (str(i:i))
        case (' ','+')
          j = j + 1
          rplc(j:j) = "_"
        case ('"',"'","/","\",",","?","&","%","*","@","!","#","$",";",":",'~','`','(',')','|')
        case default
          j = j + 1
          rplc(j:j) = str(i:i)
      end select
    enddo
    str = rplc
  end subroutine gfcgi_untaint
  !---------------------------------------------
  subroutine gfcgi_replacefield(str,srch,tag,rplc)
    !! Remove the field bracketed by srch that starts with tag, from str.
    character(len=*),intent(inout) :: str
    character(len=*),intent(in) :: srch
    character(len=*),intent(in) :: tag, rplc
    integer :: i,j,k,m,n
    i = 1
    k = len(srch)
    m = len(tag)
    j = 0
    !! diagnostic
    !  write(0,'(a)') "Replacing field >>"//srch//"<< with tag>>"//tag//"<< with >>"//rplc//"<<"
    do
      i = j + index(str(j+1:),srch) 
      if (i == j) exit                    !! no more srch
      !! diagnostic
      if (i > len(str)) exit              !! no more srch
      !write(0,'(a)') str(i:)
      if (str(i+k:i+k+m-1)==tag(1:m)) then
        !! make i,n bracket field
        i = i + k  !! 1st char after tag
        n = i + index(str(i:),srch) - 1   !! last char before tag
        if (n<=i) then                    !! or end of str
          str = str(1:i-1)//rplc
        else
          str = str(1:i-1)//rplc//str(n:)
        endif
        exit
      endif
      j = i
    enddo
  end subroutine gfcgi_replacefield
  !---------------------------------------------
  subroutine gfcgi_split(str,srch,iunit,tag)
    integer,intent(in) :: iunit
    character(len=*),intent(inout) :: str
    character(len=*),intent(in) :: srch
    character(len=*),optional,intent(in) :: tag
    integer :: i,j,k
    i = 1
    k = len(srch)
    do
      i = index(str(1:),srch)
      if (i == 0) exit
      if (present(tag)) then
        write(iunit,'(a)') str(1:i-1)//tag
      else
        write(iunit,'(a)') str(1:i-1)
      endif
      str = str(i+k:)
    enddo
  end subroutine gfcgi_split
  !---------------------------------------------

end module gfcgi
