program gf1cgi
!!----------------------------------------
!! GeoFOld CGI server program 1
!! This program parses a file uploaded by server2.php
!! using the "POST" method.
!! C.Bystroff Thu Nov 20 21:27:50 EST 2008
!!----------------------------------------
!====== FORM content from server2.html ======
!<INPUT TYPE="text" NAME="email_address" VALUE="" SIZE=40>
!<INPUT TYPE="text" NAME="keyword" VALUE="" SIZE=80>
!<p> Enter a 4-letter <A HREF="http://www.pdb.org">PDB</A> identifier
!<INPUT TYPE="text" NAME="pdbid" VALUE="" SIZE=4 >
!...or upload a file in PDB format
!<INPUT TYPE="file" NAME="pdbfile" VALUE="" MAXLENGTH=10000 SIZE=20 >
!<p><INPUT TYPE="radio" NAME="chains" VALUE="biolunit" >Use biological unit
!<INPUT TYPE="radio" CHECKED NAME="chains" VALUE="select" >Select chains
!<INPUT TYPE="hidden" NAME="params" VALUE="/ext1/bystrc/server/geofold/bin/parameters">
!<INPUT TYPE="hidden" NAME="plus" VALUE="+">
!<INPUT TYPE="hidden" NAME="slash" VALUE="/">
!<INPUT TYPE="hidden" NAME="backslash" VALUE="\">
!<INPUT TYPE="hidden" NAME="colon" VALUE=":">
!<INPUT TYPE="hidden" NAME="space" VALUE=" ">
!!=============================================
! EXAMPLE OF FILE UPLOADED FROM A BROWSER BY server2.php
!-----------------------------17957399112655092392090245054
!Content-Disposition: form-data; name="email_address"
!
!
!-----------------------------17957399112655092392090245054
!Content-Disposition: form-data; name="pdbid"
!
!
!-----------------------------17957399112655092392090245054
!Content-Disposition: form-data; name="pdbfile"; filename="1crk.pdb"
!Content-Type: text/plain
!
!
!!! ...followed by PDB records...
!!-===========================================
  use gfcgi
  implicit none
  character(len=2000) :: aline,paramline,bline
  character(len=200) :: email, jobfile,urlname,command,pdbdir,codeline,pdbupload,pdbunit
  character(len=200) :: paramfile, lname, separator, settings, emailaddress,pdbfile
  character(len=200) :: dbgfile, tmpdir, keyword
  character(len=4) :: code,pluschar,slashchar,spacechar,backslash,colon,atchar
  character(len=1) :: ach
  character(len=24) :: chains, pdbchains, ctype, biolunit
  character(len=50) :: outdir, jobdir, basedir, urldir
  integer :: pid, ios, i,j, punit=23, bunit, ich, ounit=33, dunit=43, iunit=21
  integer :: ilen, jlen
  logical,parameter :: debug=.false.
  !!--------- initialize character constants
  biolunit = "biolunit"
  pdbchains = " "
  codeline = " "
  urlname = " "
  email = " "
  lname = " "
  code = " "
  separator = " "
  pluschar = " "
  slashchar = " "
  backslash = " "
  colon = " "
  spacechar = " "
  pdbupload = " "
  pdbfile = " "
  keyword = " "
  basedir = "/bach1/home/bystrc"
  urldir = "/bach1/home/bystrc/public_html/geofold"
  tmpdir = trim(basedir)//"/server/geofold/tmp/"
  dbgfile = trim(urldir)//"/output/debug.out"
  !! --------- location of default GeoFOLD parameters file
  paramfile = trim(basedir)//"/server/geofold/bin/parameters"
  !! --------- location of local PDB repository
  pdbdir = trim(basedir)//"/server/data/pdb/"
  pdbunit = trim(basedir)//"/server/data/pdb1/"
  settings = "settings.html"
  pid = getpid()
  call gfcgi_stdinsave(iunit=iunit,pid=pid)
  !! call getenv("PDB",pdbdir)
  outdir="output/"
  jobdir="jobs/"
  if (debug) then
    open(dunit,file=trim(dbgfile),status='replace',form='formatted',iostat=ios)
    if (ios/=0) then
      write(0,'("gf1cgi.f90:: error opening ",a,i9)') trim(dbgfile), ios
      stop "gf1cgi.f90:: error opening dbgfile"
    endif
  endif
  aline = " "
  !! --------------------------- read remaining stdin data and parse
  !! NOTE: error messages go to /var/log/httpd/error_log
  !!--------- BEGIN READING STDIN from clientside
  !! First line is a unique separator, 
  !! e.g. -----------------------------1075260298824938981595028635
  read(iunit,'(a)',iostat=ios) separator
  if (ios/=0) stop 'gf1cgi2:: empty input file!!!'
  if (debug) write(dunit,'(a)') trim(separator)
  MULTILOOP: do 
    !! ------- The following reads the first line after separator
    !! ------- This must ALWAYS be the Content-Disposition line, e.g.
    !!------- !Content-Disposition: form-data; name="pdbid"
    !--if (aline(1:19)/="Content-Disposition") then
    read(iunit,'(a)',iostat=ios) aline
    if (ios/=0) stop 'gf1cgi2.f90:: Unexpected EOF after separator line.'
    if (debug) write(dunit,'(a)') trim(aline)
    !--endif
    if (aline(1:19)/="Content-Disposition") then
      stop 'gf1cgi2.f90:: ERROR. Content-Disposition line should be after every separator line.'
      !! This error message will show up in /var/log/httpd/error_log
    endif
    !! ---------- get the content type (ctype) from the name= tag
    ctype = " "
    ich = index(aline,"name=") + 6
    ilen = index(aline(ich:),'"') - 2
    ctype = aline(ich:ich+ilen) 
    !! ---------- read the blank line that always follows the ctype line
    read(iunit,'(a)',iostat=ios) bline
    if (ios/=0) stop 'gf1cgi2.f90:: Unexpected EOF after Content-Disposition line.'
    if (debug) write(dunit,'(a)') trim(bline)
    !! ---------------- read lines from stdin depending on content type
    select case (trim(ctype))
      case ("chains")
        !! ------ read text indicating chainIDs or a word indicating how chains will be chosen: select, biolunit
        read(iunit,'(a)',iostat=ios) aline
        if (ios/=0) stop 'gf1cgi2.f90:: Error reading chains line'
        if (debug) write(dunit,'(a)') trim(aline)
        if (trim(aline)==trim(separator)) cycle MULTILOOP
        biolunit = trim(aline)
      case ("lname")
        !! ---- keyword lname: read unique name for files -------!!!
        !! ------- optional 
        read(iunit,'(a)',iostat=ios) aline
        if (ios/=0) stop 'gf1cgi2.f90:: Error reading lname line'
        if (debug) write(dunit,'(a)') trim(aline)
        if (trim(aline)==trim(separator)) cycle MULTILOOP
        lname = trim(aline)
      case ("params")
        !! -------- read the filename for the GEOFOLD parameters file
        read(iunit,'(a)',iostat=ios) aline
        if (ios/=0) stop 'gf1cgi2.f90:: Error reading paramline line'
        if (debug) write(dunit,'(a)') trim(aline)
        if (trim(aline)==trim(separator)) cycle MULTILOOP
        paramline = trim(aline)
      case ("plus")
        ! ------- read the character that the clientside browser sends
        ! ------- as the plus character "+". This may be a special
        ! ------- character that must be translated back to "+".
        read(iunit,'(a1)',iostat=ios) pluschar
        if (ios/=0) stop 'gf1cgi2.f90:: Error reading pluschar line'
        if (debug) write(dunit,'(a)') trim(pluschar)
        if (trim(aline)==trim(separator)) cycle MULTILOOP
        pluschar = trim(aline)
      case ("slash")
        ! ------- same, for "/" character.
        read(iunit,'(a1)',iostat=ios) aline
        if (ios/=0) stop 'gf1cgi2.f90:: Error reading slashchar line'
        if (debug) write(dunit,'(a)') trim(aline)
        if (trim(aline)==trim(separator)) cycle MULTILOOP
        slashchar = trim(aline)
      case ("backslash")
        ! ------- same, for "\" character.
        read(iunit,'(a1)',iostat=ios) aline
        if (ios/=0) stop 'gf1cgi2.f90:: Error reading backslash line'
        if (debug) write(dunit,'(a)') trim(aline)
        if (trim(aline)==trim(separator)) cycle MULTILOOP
        backslash = trim(aline)
      case ("colon")
        ! ------- same, for ":" character.
        read(iunit,'(a1)',iostat=ios) aline
        if (ios/=0) stop 'gf1cgi2.f90:: Error reading colon line'
        if (debug) write(dunit,'(a)') trim(aline)
        if (trim(aline)==trim(separator)) cycle MULTILOOP
        colon = trim(aline)
      case ("space")
        ! ------- same, for " " character.
        read(iunit,'(a1)',iostat=ios) aline
        if (ios/=0) stop 'gf1cgi2.f90:: Error reading spacechar line'
        if (debug) write(dunit,'(a)') trim(aline)
        if (trim(aline)==trim(separator)) cycle MULTILOOP
        spacechar = trim(aline)
      case ("at")
        ! ------- same, for "@" character.
        read(iunit,'(a1)',iostat=ios) aline
        if (ios/=0) stop 'gf1cgi2.f90:: Error reading spacechar line'
        if (debug) write(dunit,'(a)') trim(aline)
        if (trim(aline)==trim(separator)) cycle MULTILOOP
        atchar = trim(aline)
      case ("email_address")
        ! ------- Read emailaddress form entry
        read(iunit,'(a)',iostat=ios) aline
        if (ios/=0) stop 'gf1cgi2.f90:: Error reading email line'
        if (debug) write(dunit,'(a)') trim(aline)
        if (trim(aline)==trim(separator)) cycle MULTILOOP
        emailaddress = trim(aline)
      case ("keyword")
        ! ------- Read keywords.
        read(iunit,'(a)',iostat=ios) aline
        if (ios/=0) stop 'gf1cgi.f90:: Error reading keyword line'
        if (debug) write(dunit,'(a)') trim(aline)
        if (trim(aline)==trim(separator)) cycle MULTILOOP
        ! keyword = 'gf_'//trim(aline)
        keyword = trim(aline)
        !! Replace any spaces with underscores, untaint.
        call gfcgi_untaint(keyword)
      case ("pdbid")
        ! ------- parse out PDB code
        read(iunit,'(a)',iostat=ios) aline
        if (ios/=0) stop 'gf1cgi2.f90:: Error reading pdbid line'
        if (debug) write(dunit,'(a)') trim(aline)
        if (trim(aline)==trim(separator)) cycle MULTILOOP
        code = aline(1:4)  !! this may be an empty line
        codeline = trim(aline)
      case ("pdbfile")
        !! ---- get file name for uploaded file, then upload it.
        ich = index(aline,"filename=") + 10
        !! --- if 1st character is close quote, there is no uploaded file.
        if (aline(ich:ich)=='"') then
          ilen = 0
          pdbupload = " "
          pdbfile = " "
        else
          !! --- use first 4 letters of uploaded file
          do ilen=1,3
            if (aline(ich+ilen:ich+ilen)=='.') exit
          enddo
          ilen = ilen - 1
          biolunit = "upload"
          pdbupload = aline(ich:ich+ilen) 
          !! ---- assume  lname precedes pdbfile in stdin --- !!!
          if (lname==" ") then
            codeline = " "
            write(codeline,'(i5.5)') pid 
          else
            codeline = trim(lname)   
          endif
          code = " " !! --- This signals NOT to get from PDB. --- !!
          pdbfile = trim(tmpdir)//trim(codeline)//".pdb"
          open(ounit,file=pdbfile,form='formatted',status='replace',iostat=ios)
          if (ios/=0) stop 'gf1cgi2.f90:: Error creating uploaded PDB file.'
        endif
        !! --------------- READ PDB FILE  from stdin
        do
          read(iunit,'(a)',iostat=ios) aline
          if (ios/=0) exit 
          if (debug) write(dunit,'(a)') trim(aline)
          if (trim(aline)==trim(separator)) then
            if (pdbfile/=" ") close(ounit)
            cycle MULTILOOP
          endif
          if (trim(aline)=="") cycle
          if (pdbfile/=" ") write(ounit,'(a)') trim(aline)
        enddo
        if (pdbfile/=" ") then
          close(ounit)
          command = "chmod 0777 "//trim(pdbfile)
          call system(command)
        endif
        if (ios/=0) exit MULTILOOP
      case (".submit")
        exit MULTILOOP
      case default
        !!---- httpd sends error messages to error_log
        write(0,'(a)') "geofold server:: gf1cgi.f90:: ERROR, unrecognized content type"
        write(0,'(a)') trim(aline)
    end select 
    do while (trim(aline)/=trim(separator))
      read(iunit,'(a)',iostat=ios) aline
      if (ios/=0) exit MULTILOOP
    enddo
  enddo MULTILOOP
  !! ------- parse out username
  i = index(emailaddress,trim(atchar))
  if (i==0) then
    i = index(emailaddress,"@")
    if (i==0) then
      i = len_trim(emailaddress)+1
      !! no @ sign in the email address. Add one.
      emailaddress = trim(emailaddress)//'@null.null'
    endif
  endif
  email = emailaddress(1:i-1)
  ! if (lname==" ") then
  !   write(lname,'(a,i5.5)') trim(email),pid 
  ! endif
  !! ----------------------------- retreive the PDB files, and get chain IDs ----
  if (code == " ") then
    !!  If code is an empty string, use the uploaded file to find the chain IDs.
    call gfcgi_getpdbchains(pdbfile,pdbchains)
  else
    !!  code is not an empty string, find the file locally or wget it from RCSB.
    select case (trim(biolunit))
    case ("select")
      !! gfcgi_findpdbfile returns empty string for pdbchains if the file was not found
      !! otherwise chain ID (_=" ")
      call gfcgi_findpdbfile(code,pdbdir,pdbchains,biolunit)
    case ("biolunit")
      call gfcgi_findpdbfile(code,pdbunit,pdbchains,biolunit)
    case ("upload")
      !! if code is not an empty string, we are not uploading. So this signals a bug.
      write(0,'(a)') "gf1cgi.f90:: WARNING, bug 1 upload, continuing"
    case default
      !! if chains are already set, then biolunit contains the chain IDs.
      !! This is the case if the program is running a "do over"
      pdbchains = trim(biolunit)
    end select
  endif
  chains = trim(pdbchains)
  !! ------------------------------- write HTML FORM output to stdout ---------------
  ! write(urlname,'(a,i5.5,".html")') trim(outdir)//trim(email),pid
  write(*,'(a,//)') 'Content-type: text/html'
  write(*,'(a)') "<html>"
  write(*,'(a)') "<head>"
  write(*,'(a)') "<title>GeoFOLD: "//trim(lname)//"</title>"
  write(*,'(a)') "</head>"
  write(*,'(a)') "<h3>GeoFOLD: "//trim(lname)//"</h3>"
  write(*,'(a)') "<body>"
  if (pdbchains==" ") then
    !! This is possible if the PDB file contains no coordinates, or if
    !! No chains were selected. Further informative error messages would be good here...
    write(*,'(a)') "<p>PDB code "//trim(code)//" not found. Please check your input. <br>"
  else
    write(*,'(a)') '<br>Select chains for PDB-code=<strong>'//trim(codeline)//'</strong>'
    !!  --- NOTE: sym-link geofold2.cgi to xgf2cgi, executable of gf2cgi.f90 ---
    write(*,'(a)') '<FORM METHOD="POST" ACTION="./geofold2.cgi" >'  
    write(*,'(a)') '<INPUT TYPE="hidden" NAME="pdbid" VALUE="'//trim(codeline)//'">'
    write(*,'(a)') '<INPUT TYPE="hidden" NAME="unit" VALUE="0">'
    write(*,'(a)') '<INPUT TYPE="hidden" NAME="lname" VALUE="'//trim(lname)//'">'
    write(*,'(a)') '<INPUT TYPE="hidden" NAME="at" VALUE="@">'
    write(*,'(a)') '<INPUT TYPE="hidden" NAME="email_address" VALUE="'//trim(emailaddress)//'">'
    write(*,'(a)') '<INPUT TYPE="hidden" NAME="keyword" VALUE="'//trim(keyword)//'">'
    !! --- NOTE: new FORM inherits pid from this job. ---
    write(*,'(a,i5.5,a)') '<INPUT TYPE="hidden" NAME="pid" VALUE="',pid,'">'
    do ich=1,len_trim(pdbchains)
      write(*,'(a)') '<INPUT TYPE="checkbox" CHECKED NAME="chains" VALUE="'&
                     //pdbchains(ich:ich)//'">'//pdbchains(ich:ich)//'<br>'
    enddo
    write(*,'(a)') "<p>Expert settings "
    write(*,'(a)') '(<A HREF="'//trim(settings)//'">Click for more info.</A>):<br>'
    write(*,'(a)') '<textarea NAME="params" rows=20 cols=60 WRAP="off">'
    ! if (paramfile==" ") paramfile = "/ext1/bystrc/server/geofold/bin/parameters"
    !! ----------- write the parameter file to the textarea
    open(punit,file=trim(paramfile),status='old',form='formatted',iostat=ios)
    if (ios/=0) then
      write(0,'(a,a)') "!----- ERROR parameter file not found",trim(paramfile)
      write(*,'(a,a)') "!----- ERROR parameter file not found",trim(paramfile)
      write(*,'(a)') "!----- Using default parameters."
    else
      j = 0
      do
        read(punit,'(a)',iostat=ios) aline
        if (ios/=0) exit
        if (aline(1:7)=="PDBCODE") then
          aline(9:) = trim(codeline)
        elseif (aline(1:5)=="LNAME") then
          aline(7:) = lname
        elseif (aline(1:7)=="KEYWORD") then
          aline(9:) = trim(keyword)
        elseif (aline(1:5)=="CHAIN") then
          aline(7:) = pdbchains
        elseif (aline(1:5)=="EMAIL") then
          aline(7:) = trim(emailaddress)
          j = 1
        endif
        write(*,'(a)') trim(aline)
      enddo
      !! If for some reason this keyword is missing...
      if (j==0) write(*,'(a)') "EMAIL "//trim(emailaddress)
      close(punit)
    endif
    write(*,'(a)') '</textarea>'
    write(*,'(a)') '<p><INPUT TYPE="submit" NAME=".submit" VALUE="Submit GeoFOLD job">'
  endif
  write(*,'(a)') "</body>"
  write(*,'(a)') "</html>"
  if (debug) close(dunit)

end program gf1cgi
