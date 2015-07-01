program isegment
!! Capture the label from a map file
!! and display it
!! GeoFOld CGI server
!! C.Bystroff Fri Jul 25 08:06:35 EDT 2008
!! Last modified: Tue Aug 27 12:36:46 EDT 2013
  implicit none
  character(len=200) :: aline, basename
  character(len=200) :: label, dagfile
  character(len=4) :: code
  integer :: pid, ios, i,j
  character(len=7),parameter :: outputurl = 'output/'
  !!---------
  !! uses method="POST" 
  !!---------------------
  !!========== UNCOMMENT TO DEBUG =============
  !write(*,'(a,//)') 'Content-type: text/html'
  !  write(*,'(a)') "<html>"
  !  write(*,'(a)') "<body>"
  !  write(*,'(a)') "<pre>"
  !call system("set")
  !  write(*,'(a)') "</pre>"
  !  write(*,'(a)') "</body>"
  !  write(*,'(a)') "</html>"
  !stop
  !!---------------------
  ! read(*,*) aline
  call getenv("QUERY_STRING",aline)
  call parseit(aline,"iseg",label)
  call parseit(aline,"dag",dagfile)
  ! call singlelineoutput(label)
  call readdagfile(label,dagfile)
  i = index(dagfile,'/',.true.)+1
  !! NOTE: naming convention in RUNGEOFOLD.csh is <basename>_<N>.dag.out
  j = index(dagfile,'.dag',.true.)-1
  ! j = index(dagfile(i:j),'_',.true.)+i-2
  basename = dagfile(i:j-2)
  write(*,'("<h4><A HREF=",a,">Back to main</A></h4>")') &
          '"'//trim(outputurl)//trim(basename)//'.html"'
  write(*,'(a)') "</body>"
  write(*,'(a)') "</html>"
CONTAINS
  !---------------------------------------------
  subroutine parseit(aline,atag,result)
    character(len=*),intent(in) :: aline, atag
    character(len=*),intent(out) :: result
    integer :: i,j,ich, k
    i = index(aline,atag)
    ! write(*,*) aline(1:i)
    j = i + len(atag)+1
    k = index(aline(j:),'&')
    result = aline(j:j+k-2)
    ! write(*,*) trim(result)
  end subroutine parseit
  !---------------------------------------------
  subroutine errormsg(atag,entry)
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
  end subroutine errormsg
  !---------------------------------------------
  subroutine istatelineoutput(line1,line2,dagfile)
    !! Output HTML for ISEGMT, plus links to TSTATEs
    !! Example ISEGMT line:
    !!ISEGMT     104    5    0      624.24       11.28       0     0       U  0.01686023
    character(len=*),intent(in) :: line1, line2, dagfile
    integer :: i,j,k,m,nseg,nvoid,nhb,dunit=22,ios=0,nts,f,u1,u2,nsym,nres,nch,iblock
    real :: conc, sas, ntrp
    integer,dimension(8) :: seams
    character(len=1) :: ftype, lstch
    character(len=1000) :: aline, filename
    integer,parameter :: blsize=80  !! blsize must be a multiple of 10, and <= 100
    character(len=100) :: &
     tics="....,....|....,....|....,....|....,....|....,....|....,....|....,....|....,....|....,....|....,....|"
    character(len=101) :: &
     nums="1--------10--------20--------30--------40--------50--------60--------70--------80--------90-------100"
    read(line1(7:),*) nseg,i,nsym,sas,ntrp,nvoid,nhb,ftype,conc
    write(*,'(a,//)') 'Content-type: text/html'
    write(*,'(a)') "<html>"
    write(*,'(a)') "<body>"
    write(*,'(a,a,a)') '<font color="#11AA55"><pre>',trim(line1),'</pre></font>'
    j = index(line2," ")-1
    write(*,'(a,a,a)') '<font color="#112211"><pre>',nums(1:blsize+1),'</pre></font>'
    write(*,'(a,a,a)') '<font color="#112211"><pre>',tics(1:blsize),'</pre></font>'
    do iblock=1,int(j/blsize +1)
      k = (iblock-1)*blsize+1
      m = iblock*blsize
      if (m > j) m = j
      write(*,'(a,a," ",i0,a)') &
      '<font color="#AA1155"><pre>',line2(k:m),iblock*blsize,'</pre></font>'
    enddo
    write(*,'("<table border=1>")') 
    write(*,'(a,i9,a)') '<tr><td bgcolor="#BBAAFF">Segment node</td><td bgcolor="#FE8877"> ',nseg,'</td></tr>'
    select case (ftype)
      case ("F")
        write(*,'(a)') '<tr><td bgcolor="#CC99FF">Folded state</td><td bgcolor="#CAFABA">FOLDED</td></tr>'
      case ("I")
        write(*,'(a)') '<tr><td bgcolor="#CC99FF">Folded state</td><td bgcolor="#BACAFA">INTERMEDIATE</td></tr>'
      case ("U")
        write(*,'(a)') '<tr><td bgcolor="#CC99FF">Folded state</td><td bgcolor="#FACADA">UNFOLDED</td></tr>'
    end select
    nres = 0
    nch = 0
    lstch = "."
    do i=1,index(line2," ")-1
      if (line2(i:i)/='.') then
        nres = nres + 1
        if (lstch/=line2(i:i)) nch = nch + 1
      endif
      lstch = line2(i:i)
    enddo
    write(*,'(a,i0,a)')   &
         '<tr><td bgcolor="#BBAAFF">Size     </td><td bgcolor="#DACAFA"> ',nres,' </td></tr>'
    write(*,'(a,i0,a)')   &
         '<tr><td bgcolor="#AABBFF">Number of chains </td><td bgcolor="#CADAFA"> ',nch,' </td></tr>'
    write(*,'(a,i0,a)')   &
         '<tr><td bgcolor="#BBAAFF">Symmetry-generated copies    </td><td bgcolor="#DACAFA"> ',nsym,' </td></tr>'
    write(*,'(a,f9.2,a)') &
         '<tr><td bgcolor="#AABBFF">Total buried surface area    </td><td bgcolor="#CADAFA"> ',sas,' </td></tr>'
    write(*,'(a,f9.2,a)') &
         '<tr><td bgcolor="#BBAAFF">Total unexpressed sidechain entropy    </td><td bgcolor="#DACAFA"> ',ntrp,' </td></tr>'
    write(*,'(a,i0,a)')   &
         '<tr><td bgcolor="#AABBFF">Total remaining buried voids </td><td bgcolor="#CADAFA"> ',nvoid,' </td></tr>'
    write(*,'(a,i0,a)')   &
         '<tr><td bgcolor="#BBAAFF">Total remaining H-bonds </td><td bgcolor="#DACAFA"> ',nhb,' </td></tr>'
    write(*,'(a,E8.2e2,a)') &
         '<tr><td bgcolor="#AABBFF">Concentration at end of simulation</td><td bgcolor="#CADAFA"> ',conc,' </td></tr>'
    write(*,'(a)') '</table>' 
    filename = trim(outputurl)//trim(dagfile)
    open(dunit,file=filename,form='formatted',status='old',iostat=ios)
    if (ios/=0) then
      call singlelineoutput("Error opening dagfile line  139 "//trim(filename))
    endif
    write(*,'("<h5>Transition states (folding)</h5>")') 
    write(*,'("<pre>")') 
    write(*,'(7x,4a7)') "tstate","f","u1","u2"
    do
      read(dunit,'(a)',iostat=ios) aline
      if (ios/=0) exit
      if (aline(1:6)/="TSTATE") cycle
      read(aline(7:),*) nts,f,u1,u2
      if (u1==nseg.or.u2==nseg) then
         write(*,'(a7,4(a,i0,a,i7,a),a)') aline(1:7),&
           '<A HREF="isegment.cgi?iseg=tu',nts,'&dag='//trim(dagfile)//'&">',nts,'</A>', &
           '<A HREF="isegment.cgi?iseg=n',f,'&dag='//trim(dagfile)//'&">',f,'</A>', &
           '<A HREF="isegment.cgi?iseg=n',u1,'&dag='//trim(dagfile)//'&">',u1,'</A>', &
           '<A HREF="isegment.cgi?iseg=n',u2,'&dag='//trim(dagfile)//'&">',u2,'</A>', &
           trim(aline(36:))
      endif
    enddo
    write(*,'("</pre>")') 
    rewind(dunit)
    write(*,'("<h5>Transition states (unfolding)</h5>")')
    write(*,'("<pre>")') 
    write(*,'(7x,4a7)') "tstate","f","u1","u2"
    do
      read(dunit,'(a)',iostat=ios) aline
      if (ios/=0) exit
      if (aline(1:6)/="TSTATE") cycle
      read(aline(7:),*) nts,f,u1,u2
      if (f==nseg) then
         write(*,'(a7,4(a,i0,a,i7,a),a)') aline(1:7),&
           '<A HREF="isegment.cgi?iseg=tu',nts,'&dag='//trim(dagfile)//'&">',nts,'</A>', &
           '<A HREF="isegment.cgi?iseg=n',f,'&dag='//trim(dagfile)//'&">',f,'</A>', &
           '<A HREF="isegment.cgi?iseg=n',u1,'&dag='//trim(dagfile)//'&">',u1,'</A>', &
           '<A HREF="isegment.cgi?iseg=n',u2,'&dag='//trim(dagfile)//'&">',u2,'</A>', &
           trim(aline(36:))
      endif
    enddo
    write(*,'("</pre>")') 
    close(dunit)
    !write(*,'(a)') "</body>"
    !write(*,'(a)') "</html>"
  end subroutine istatelineoutput
  !---------------------------------------------
  subroutine tstatelineoutput(entry,dagfile)
    character(len=*),intent(in) :: entry
    character(len=*),intent(in) :: dagfile
    character(len=20) :: cuttype
    integer :: nts,f,u1,u2,iseam,iunit=98,ios=0,i,nb,x1,x2,y1,y2
    real :: ntrp, traf, nrg
    character :: ctype
    character(len=200) :: aline, filename
    !!----
    filename = trim(outputurl)//trim(dagfile)
    read(entry(7:),*,iostat=ios) nts,f,u1,u2,ntrp,ctype,iseam,traf
    if (ios/=0)  then
      call singlelineoutput("ERROR parsing TSTATE line : "//trim(entry))
    endif
    select case (ctype)
      case ("h")
        cuttype = "HINGE"
      case ("p")
        cuttype = "PIVOT"
      case ("b")
        cuttype = "BREAK"
      case ("m")
        cuttype = "MELTING"
      case ("s")
        cuttype = "SEAM"
      case default
        cuttype ="UNKNOWN CUTTYPE"
    end select
    write(*,'(a,//)') 'Content-type: text/html'
    write(*,'(a)') "<html>"
    write(*,'(a)') "<body>"
    write(*,'(a,a,a)') '<font color="#000055"><pre>',trim(entry),'</pre></font>'
    write(*,'("<br><b>Transition state ",i9,"</b>&nbsp;",a)') nts,trim(cuttype)
    write(*,'("<br>folded state node number     =<A HREF=",a,i0,a,">",i0,"</A>")') &
           '"'//'isegment.cgi?iseg=n',f,'&dag='//trim(dagfile)//'&"',f
    write(*,'("<br>unfolded state 1 node number =<A HREF=",a,i0,a,">",i0,"</A>")') &
           '"'//"isegment.cgi?iseg=n",u1,"&dag="//trim(dagfile)//'&"',u1
    if (cuttype/="SEAM") then
      write(*,'("<br>unfolded state 2 node number =<A HREF=",a,i0,a,">",i0,"</A>")') &
           '"'//"isegment.cgi?iseg=n",u2,"&dag="//trim(dagfile)//'&"',u2
    endif
    if (scan(ctype,"hpb")/=0) then
        write(*,'("<br>",a,": entropy = ",f8.3)') trim(cuttype),ntrp
    endif
    write(*,'("<br>Relative traffic =",f8.3)') traf
    if (cuttype=="SEAM") then
      open(iunit,file=filename,form='formatted',status='old',iostat=ios)
      if (ios/=0)  then
         write(*,'(a,a)') '<font color="#FF0000">BUG! DAG file open statement in tstatelineoutput</font>',trim(filename)
         write(*,'(a)') "</body>"
         write(*,'(a)') "</html>"
         stop
      endif
      write(*,'("<br>Seam ",i5," limits=")') iseam
      iseam = abs(iseam)
      do
        read(iunit,'(a)',iostat=ios) aline
        if (ios/=0) exit
        if (aline(1:5)/="SEAM ") cycle
        !! SEAM line format must conform with pdb2seams and geofold.f90
        read(aline(6:),*,iostat=ios) i,nb,nrg,x1,x2,y1,y2
        if (ios/=0)  then
           write(*,'(a,a)') '<font color="#AA6633">BUG! ERROR parsing SEAM line</font>',trim(aline)
           write(*,'(a)') "</body>"
           write(*,'(a)') "</html>"
           stop
        endif
        if (i/=iseam) cycle
        write(*,'(4i5)') x1,x2,y1,y2
        exit
      enddo
      close(iunit)
    endif
    !write(*,'(a)') "</body>"
    !write(*,'(a)') "</html>"
  end subroutine tstatelineoutput
  !---------------------------------------------
  subroutine singlelineoutput(entry)
    character(len=*),intent(in) :: entry
    write(*,'(a,//)') 'Content-type: text/html'
    write(*,'(a)') "<html>"
    write(*,'(a)') "<body>"
    write(*,'(a,a,a)') '<font color="#000055"><pre>',trim(entry),'</pre></font>'
    write(*,'(a)') "</body>"
    write(*,'(a)') "</html>"
    stop
  end subroutine singlelineoutput
  !---------------------------------------------
  subroutine readdagfile(label,dagfile)
    character(len=*),intent(in) :: label, dagfile
    character(len=2000) :: aline,bline
    integer :: ios, inode, i, dunit=22
    logical :: iststate
    character(len=200) :: filename
    !! 
    filename = trim(outputurl)//trim(dagfile)
    open(dunit,file=filename,form='formatted',status='old',iostat=ios)
    if (ios/=0 ) call singlelineoutput(entry="DAG file is missing: "//trim(dagfile))
    !! The label structure for intermediates is n<inode> (a1,i) and for transition states is t<cuttype><inode> (a1,a1,i)
    !! These lines must conform with maxTraffic.cpp
    if (label(1:1)=="t") then
      read(label(3:),*,iostat=ios) inode
      if (ios/=0 ) call singlelineoutput(entry="ERROR parsing tstate node label: "//trim(label))
      iststate = .true.
    elseif (label(1:1)=="n") then
      read(label(2:),*,iostat=ios) inode
      if (ios/=0 ) &
      call singlelineoutput(entry="ERROR parsing intermediate node label: "//trim(label))
      iststate = .false.
    else
      call singlelineoutput(entry="ERROR bad node label: "//trim(label))
    endif
    !! read DAG file
    do
      read(dunit,'(a)',iostat=ios) aline
      if (ios/=0) exit
      if (iststate) then
        if (aline(1:6)/="TSTATE") cycle
        read(aline(7:14),*,iostat=ios) i
        if (ios/=0 ) call singlelineoutput("ERROR parsing dagfile : "//trim(aline))
        if (i==inode) then
          close(dunit)  !! must close, because file will be opened in the subroutine
          call tstatelineoutput(aline,dagfile)
          return    !! must return, because file is closed.
        endif
      else
        if (aline(1:6)/="ISEGMT") cycle
        read(aline(7:14),*,iostat=ios) i
        read(dunit,'(a)',iostat=ios) bline
        if (ios/=0 ) call singlelineoutput("ERROR parsing dagfile : "//trim(aline))
        if (i==inode) then
          close(dunit)  !! must close, because file will be opened in the subroutine
          call istatelineoutput(aline,bline,dagfile)
          return    !! must return, because file is closed.
        endif
      endif
    enddo
    call singlelineoutput("ERROR parsing dagfile : "//trim(dagfile)//"  Node line not found.")
    close(dunit)  !! If here, there is a bug.
  end subroutine readdagfile
end program isegment
