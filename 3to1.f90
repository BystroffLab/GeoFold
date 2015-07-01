        program three2one
!!* Modified to stop when a ENDMDL line is encountered.
!!*  CB. Wed Sep  1 10:02:45 EDT 2004
!!* new version. output is fasta  4-MAR-99
        implicit none
!!* convert 3-letter sequence in SEQRES lines to one-letter code.
!!* 27-AUG-96 --CB
!!8
        character(len=80) :: aline
        character(len=5000) :: bline, cline
        character(len=3) :: three(20)
        character(len=4) :: code
        character(len=1) :: res1(21),chain,altloc,allchain
        integer :: I,J,K,L,lcount,bch,N,iline,nres , jarg
        character(len=6) :: last
!!* function
        integer :: iargc, ios
        three= (/'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE', &
        'LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL', &
        'TRP','TYR'/)
        res1=(/'A','C','D','E','F','G','H','I','K','L','M','N', &
        'P','Q','R','S','T','V','W','Y','X'/)
!!**
        bline = ' '
        code = " "
        bch = 0
        nres = 0
        N = 0
        iline = 0
        chain = '?'
        altloc = ' '
        allchain = 'n'
        jarg = iargc()
        if (jarg.ge.1) then
          call getarg(1,chain)
          if (chain.eq.'_') chain = ' '
          if (chain.eq.'+') allchain = 'y'
          if (chain.eq.'.') allchain = 'y'
        endif
!!* find and read SEQRES lines
      do
        read(*,'(a80)',iostat=ios) aline
        if (ios/=0) exit
        if (aline(1:6).eq.'HEADER') then
          code = aline(63:66)
        endif
        if (aline(1:6).eq.'ATOM  ') exit  ! no SEQRES lines, use ATOM lines instead
        if (aline(1:6).ne.'SEQRES') cycle
        if (chain.eq.'?') chain = aline(12:12)
        if (allchain.eq.'y') chain = aline(12:12)
        if (aline(12:12).ne.chain) cycle
        if (iline.eq.0) then
          read(aline(13:17),*) nres
          if (chain.eq.' ') chain = '_'
          if (code==" ") code = aline(73:76)
          write(*,'("> Code: ",a4," Chain: ",a1,i9," residues. From SEQRES lines.")')  &
          code,chain,nres
          if (chain.eq.'_') chain = ' '
        endif
        iline = iline + 1
        lcount = 1
        I = 20
        do while ((aline(I:I+3).ne.'    ').and.(lcount.le.13))
          K = 21
          do J=1,20
            if (aline(I:I+2).eq.three(J)) K=J
          enddo
          bch = bch + 1
          N = N + 1
          bline(bch:bch) = res1(K)
          if (bch.ge.50) then
            write(*,'(a50)') bline(1:50)
            bline = ' '
            bch = 0
          endif
          lcount = lcount + 1
          I = I + 4
        enddo
      enddo
      if (bch.gt.0) then
        write(*,'(a50)') bline(1:50)
        bch = 0
      endif
      if (N.ne.nres.and.nres.ne.0) then
        write(0,'("# ERROR: number read not equal to number advertized in SEQRES lines.",' &
        //'2i9)') N,nres
      endif
      last = ' '
      if (iline.eq.0) then
        bch = 0
        if (chain.eq.'_') chain = ' '
        nres = 0
        do 
          if (iline.gt.0) read(*,'(a80)',iostat=ios) aline
          if (ios/=0) exit
          iline = iline + 1
          if (aline(1:6).eq.'ENDMDL') exit  !! use first model of NMR structures
          if (aline(1:6).ne.'ATOM  '.and.aline(1:6).ne.'HETATM') cycle
          if (chain.eq.'?') chain = aline(22:22) ! pick first chain
          if (allchain=='y') chain = aline(22:22) ! use all chains
          if (aline(22:22).ne.chain) cycle
          if (aline(13:16).ne.' CA ') cycle
          if (aline(22:26).eq.last) cycle
          last = aline(22:26)
          nres = nres + 1
          K = 21
          do J=1,20
            if (aline(18:20).eq.three(J)) K=J
          enddo
          bch = bch + 1
          N = N + 1
          bline(bch:bch) = res1(K)
          cline(bch:bch) = chain
        enddo
        if (chain.eq.' ') chain = '_'
        if (code==" ") code = aline(73:76)
        if (allchain=='y') then
          write(*,'("# Code: ",a4," Chain: all",i9,'//&
                  '" residues. From ATOM lines.")')  &
          code,nres
        else
          write(*,'("# Code: ",a4," Chain: ",a1,i9,'//&
                  '" residues. From ATOM lines.")')  &
          code,chain,nres
        endif
        N = 0
        chain = cline(1:1)
        do k=1,nres
          if (cline(k:k)/=chain) then
            write(*,"(/)")
            N = 0
            chain = cline(k:k)
          endif
          write(*,'(a,$)') bline(k:k)
          N = N + 1
          if (N==50) then
            write(*,*)
            N = 0
          endif
        enddo
      endif
      write(*,*)
      !! if (bch.gt.0) write(*,'(a50)') bline
    end program three2one
