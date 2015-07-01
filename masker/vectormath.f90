!****** PROTEAN_MATH routines for PROTEAN
!* Modified for F90. *x converted to real kind
!******************************************************************
module vectormath
          !!use stochastic
          implicit none
          integer,parameter :: kind_8=selected_real_kind(P=8)
          real,private,parameter :: pi=3.1415927
CONTAINS
!******************************************************************
        subroutine align(a1,a2,a3,b1,mat,vec)
!* rotate around vector a1->a2 to put a3
!* in the plane defined by a1,a2,b1
!* return the matrix and vector for the operation.
        implicit none
        real(kind=kind_8) a1(3),a2(3),a3(3),b1(3),vec(3)
        real(kind=kind_8) mat(3,3)
        real(kind=kind_8) a(3),b(3),c(3),ab(3),ac(3),aba(3),dab,dac,daba, &
               ab_ac,aba_ac,chi
        integer :: i
!* get difference vectors
        do I=1,3
        a(I)=a1(I)-a2(I)
        b(I)=a3(I)-a2(I)
        c(I)=b1(I)-a2(I)
        enddo
        call cross(b,a,ab,dab)
        call cross(c,a,ac,dac)
        call cross(ab,a,aba,daba)
!* get vector lengths
!* get ab.dot.ac
        ab_ac = ab(1)*ac(1) +ab(2)*ac(2) +ab(3)*ac(3) 
!*     aba.dot.ac
        aba_ac = aba(1)*ac(1) +aba(2)*ac(2) +aba(3)*ac(3) 
!* cosine ab_ac
        ab_ac=ab_ac/dab
!* cosine aba_ac = sine ab_ac
        aba_ac=aba_ac/daba
!* get chi
        chi=datan2(aba_ac,ab_ac)
!* get matrix for chi rotation
        call getrot(a1(1),a2(1),chi,mat,vec)
        return
        end subroutine align
!***********************************************************************
        subroutine lineup(a1,a2,b1,b2,mat,vec)
!* find the matrix (mat) and vector (vec) that
!* line up a1->a2 with b1->b2, with a1 superimposed 
!* on b1.
!        implicit none
        real(kind=kind_8) :: a1(3),a2(3),b1(3),b2(3),vec(3)
        real(kind=kind_8) :: mat(3,3)
        real(kind=kind_8) :: a(3),b(3),c(3),d(3),da,db,dab,phi,psi,chi
        integer :: i
!*
        chi=0.0
        phi=0.0
        psi=0.0
!* get difference vectors
        do I=1,3
        a(I)=a2(I)-a1(I)
        b(I)=b2(I)-b1(I)
        enddo
!* get vector lengths
!* get a.dot.b
        da=a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
        db=b(1)*b(1) + b(2)*b(2) + b(3)*b(3)
        dab=a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
        da=dsqrt(da)
        db=dsqrt(db)
!* get chi
!cb The following is to protect against underflows on the Mac:
        if (abs(dab-(da*db)).lt.0.000001) then
          chi = 0.0
        else
          chi=dacos(dab/(da*db))
        endif
        if (chi.ge.1.0E-4) then
!* cross product A X B = C
          call cross(a,b,c,dab)
!* get phi, psi
          psi=dacos(c(3)/dab)
          if (c(1).ne.0.0.or.c(2).ne.0.0) phi=datan2(-c(1),c(2))
        endif
!* get matrix
        call getmat(phi,psi,chi,mat)
!* get vector, = -Mat*a1 + b1
        call rotate(mat,a1,c)
        do I=1,3
        vec(I)=-c(I)+b1(I)
        enddo
        return
        end subroutine lineup
!***********************************************************************
        subroutine mover(c,mat,vec)
        implicit none
        real(kind=kind_8) c(3),t(3),mat(3,3),vec(3)
        call rotate(mat,c,t)
        c(1) = t(1) + vec(1)
        c(2) = t(2) + vec(2)
        c(3) = t(3) + vec(3)
        return
        end subroutine mover
!***********************************************************************
        subroutine move_S(c,mat,vec)
        implicit none
        real c(3),t(3),mat(3,3),vec(3)
        call rotate_S(mat,c,t)
        c(1) = t(1) + vec(1)
        c(2) = t(2) + vec(2)
        c(3) = t(3) + vec(3)
        return
        end subroutine move_S
!***********************************************************************
        subroutine rotate(mat,v1,v2)
        implicit none
        real(kind=kind_8) mat(3,3),v1(3),v2(3)
        integer :: i
        do I=1,3
        v2(I)=mat(1,I)*v1(1)+mat(2,I)*v1(2)+mat(3,I)*v1(3)
        enddo
        return
        end subroutine rotate
!***********************************************************************
        subroutine rotate_S(mat,v1,v2)
        implicit none
        real mat(3,3),v1(3),v2(3)
        integer :: i
        do I=1,3
        v2(I)=mat(1,I)*v1(1)+mat(2,I)*v1(2)+mat(3,I)*v1(3)
        enddo
        return
        endsubroutine rotate_S
!***********************************************************************
        subroutine getmat(phi,psi,kappa,aa)
        implicit none
!* get matrix for Tanaka convention polar angles
!* CB  14-JAN-1991, subroutine 3-JUN-1991
!* Double precision  15-AUG-1991
!* >>> phi,psi,kappa in radians <<<
        real(kind=kind_8),intent(inout) :: phi,psi,kappa
        real(kind=kind_8),intent(out) :: aa(3,3)
        real(kind=kind_8) sf,cf,ss,cs,sk,ck,ss2,sf2,cs2,cf2
        integer :: i,j
!***
!cb Try to avoid underflows...
        if (dabs(phi).lt.1.0E-4) phi = 0.d0
        if (dabs(psi).lt.1.0E-4) psi = 0.d0
        if (dabs(kappa).lt.1.0E-4) then
          aa=0.d0
          do I=1,3
            aa(I,I)=1.d0
          enddo
          return
        endif
        sf = dsin(phi)
        cf = dcos(phi)
        ss = dsin(psi)
        cs = dcos(psi)
        sk = dsin(-kappa)
        ck = dcos(-kappa)
        ss2 = ss*ss
        sf2 = sf*sf
        cf2 = cf*cf
        cs2 = cs*cs
!* now calculate the product matrix of the five rotation matrices
!! revised Wed Dec 31 13:12:58 EST 2003
        aa(1,1)=(cf2+sf2*cs2)*ck+sf2*ss2
        aa(1,2)=ss2*sf*cf*(ck-1.0)-cs*sk
        aa(1,3)=ss*(sf*cs*(ck-1.0)+cf*sk)
        aa(2,1)=ss2*cf*sf*(ck-1.0)+cs*sk
        aa(2,2)=(sf2+cf2*cs2)*ck+cf2*ss2
        aa(2,3)=ss*(cf*cs*(1.0-ck)+sf*sk)
        aa(3,1)=ss*(cs*sf*(ck-1.0)-cf*sk)
        aa(3,2)=ss*(cf*cs*(1.0-ck)-sf*sk)
        aa(3,3)=ss2*ck+cs2
!* done.
!cb      write(ounit,'(/,3(3f12.8,/))') ((aa(I,J),I=1,3),J=1,3)
        return
        end subroutine getmat
!***********************************************************************
        subroutine getrot(a1,a2,chi,mat,vec)
        implicit none
!* find the matrix (mat) and vector (vec) that
!* rotate chi rad about an axis defined by a1-->a2
!* rotation is right-handed. That is, it is a
!* clockwise rotation when looking from a1 to a2
!* === double precision ===
!* chi in radians
        real(kind=kind_8) a1(3),a2(3),vec(3)
        real(kind=kind_8) mat(3,3)
        real(kind=kind_8) a(3),c(3),da,phi,psi,chi,rad
        integer :: i,j,k
!*
        if (dabs(chi).lt.1.0E-4) goto 100
        phi=0.0
        psi=0.0
!* get difference vector: a
        call Dsubvec(a2,a1,a)
!* get length of vector A
        da=a(1)*a(1) +a(2)*a(2) +a(3)*a(3) 
        da=dsqrt(da)
!* get phi, psi
        psi=dacos(a(3)/da)
        if (a(1).ne.0.0.or.a(2).ne.0.0) phi=datan2(-a(1),a(2))
!* get matrix
        call getmat(phi,psi,chi,mat)
!* get vector, = -Mat*a1 + a1
        call rotate(mat,a1,c)
        call Dsubvec(a1,c,vec)
        return
 100    do I=1,3
          do K=1,3
            mat(K,I)=0.d0
          enddo
          mat(I,I)=1.d0
          vec(I)=0.d0
        enddo
        return
        end subroutine getrot
!***********************************************************************
        subroutine cross(a,b,c,d)
        implicit none
!* a X b = c, length of c  = d
        real(kind=kind_8) :: a(3),b(3),c(3),d
!* cross product A X B
        c(1)=a(2)*b(3)-a(3)*b(2)
        c(2)=a(3)*b(1)-a(1)*b(3)
        c(3)=a(1)*b(2)-a(2)*b(1)
        d=sqrt(c(1)*c(1)+c(2)*c(2)+c(3)*c(3))
        return
        end subroutine cross
!***********************************************************************
        subroutine cros(v1,v2,v3)
        implicit none
!* v3 = v1 X v2
        real v1(3), v2(3), v3(3)
        v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
        v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
        v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
        return
        end subroutine cros
!***********************************************************************
        subroutine Dsubvec(v1,v2,v3)
        implicit none
        real(kind=kind_8) v1(3),v2(3),v3(3)
!* v3 = v1 - v2
        v3(1) = v1(1) - v2(1)
        v3(2) = v1(2) - v2(2)
        v3(3) = v1(3) - v2(3)
        return
        end subroutine Dsubvec
!***********************************************************************
        subroutine subvec(v1,v2,v3)
        implicit none
        real v1(3),v2(3),v3(3)
!* v3 = v1 - v2
        v3(1) = v1(1) - v2(1)
        v3(2) = v1(2) - v2(2)
        v3(3) = v1(3) - v2(3)
        return
        end subroutine subvec
!***********************************************************************
        subroutine divvec(v1,x,v2)
        implicit none
        real v1(3),x,v2(3)
        v2(1) = v1(1)/x
        v2(2) = v1(2)/x
        v2(3) = v1(3)/x
        return
        end subroutine divvec
!***********************************************************************
        subroutine mulvec(v1,x,v2)
        implicit none
        real v1(3),x,v2(3)
        v2(1) = v1(1)*x
        v2(2) = v1(2)*x
        v2(3) = v1(3)*x
        return
        end subroutine mulvec
!***********************************************************************
        subroutine addvec(v1,v2,v3)
        implicit none
        real v1(3),v2(3),v3(3)
        v3(1) = v1(1) + v2(1)
        v3(2) = v1(2) + v2(2)
        v3(3) = v1(3) + v2(3)
        return
        end subroutine addvec
!***********************************************************************
        subroutine unitvec(v1,v2)
        implicit none
        real v1(3),v2(3) !! ,cc
        !! cc = sqrt(dotprod(v1,v1))
        !! call divvec(v1,cc,v2)
        v2 = v1/sqrt(dotprod(v1,v1))
        return
        end subroutine unitvec
!***********************************************************************
        subroutine unitvecD(v1,v2)
        implicit none
        real(kind=kind_8) v1(3),v2(3),cc
        !! cc = dsqrt(Ddotprod(v1,v1))  !! new: sqrt in f90 takes type of args.
        cc = sqrt(Ddotprod(v1,v1))
        v2 = v1/cc
        return
        end subroutine unitvecD
!***********************************************************************
        real function dotprod(v1,v2)
        implicit none
        real v1(3),v2(3)
        dotprod = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)  
        return
        end function dotprod
!***********************************************************************
        real(kind=kind_8) function Ddotprod(v1,v2)
        implicit none
        real(kind=kind_8) v1(3),v2(3)
        Ddotprod = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)  
        return
        end function Ddotprod
!***********************************************************************
        subroutine MM_D(mat1,mat2,mat3)
        implicit none
!* Multiply two double precision rotation matrixes
!* mat3 = mat2*mat1
!* NOTE:the order of rotations is therefore mat1 -> mat2
!* 5-DEC-94
!* CB
        real(kind=kind_8) mat1(3,3),mat2(3,3),mat3(3,3)
        integer I,J,K
!
        do I=1,3
          do J=1,3
            mat3(J,I) = 0.
            do K=1,3
              mat3(J,I) = mat3(J,I) + mat2(K,I)*mat1(J,K)
            enddo
          enddo
        enddo
        return
        end subroutine MM_D
!***********************************************************************
        subroutine MM_S(mat1,mat2,mat3)
        implicit none
!* Multiply two single precision rotation matrixes
!* mat3 = mat2*mat1
!* NOTE:the order of rotations is therefore mat1 -> mat2
!* 6-AUG-01
!* CB
        real mat1(3,3),mat2(3,3),mat3(3,3)
        integer I,J,K
!
        do I=1,3
          do J=1,3
            mat3(J,I) = 0.
            do K=1,3
              mat3(J,I) = mat3(J,I) + mat2(K,I)*mat1(J,K)
            enddo
          enddo
        enddo
        return
        end subroutine MM_S
!***********************************************************************
!        subroutine copymat(mat1,mat2)
!        implicit none
!* 5-DEC-94 double precision
!* CB
!! Obselete.
!        real(kind=kind_8) mat1(3,3),mat2(3,3)
!        integer I,J
!
!        do I=1,3
!          do J=1,3
!            mat2(J,I) = mat1(J,I) 
!          enddo
!        enddo
!        return
!        end
!***********************************************************************
        subroutine MATMULT_S(mat1,vec1,mat2,vec2)
        implicit none
!* This routine applies the operation contained in
!* matrix mat2 and vector vec2 and applies it to mat1, vec1.
!* The resulting operation is returned in mat1 and vec1.
!* mat2 ( mat1 * a + vec1) + vec2 = b
!* = mat2 * mat1 * a + (mat2 * vec1  +  vec2) = b
!* Thus the new matrix is mat2 * mat1 and the new vector is
!* (mat2 * vec1  +  vec2).
!* 19-OCT-94 single-precision
!* CB
        real mat1(3,3),vec1(3),mat2(3,3),vec2(3)
        real tmp(3,3),t(3)
        integer I,J,K
!
        do I=1,3
          t(I) = 0.
          do J=1,3
            tmp(J,I) = 0.
            do K=1,3
              tmp(J,I) = tmp(J,I) + mat2(K,I)*mat1(J,K)
            enddo
            t(I) = t(I) + mat2(J,I)*vec1(J)
          enddo
          vec1(I) = t(I) + vec2(I)
        enddo
        do I=1,3
          do J=1,3
            mat1(J,I) = tmp(J,I)
          enddo
        enddo
        return
        end subroutine MATMULT_S
!***********************************************************************
        subroutine LUDCMP(A,N,INDX,D)
        implicit none
!* This is the first part of the LU decomposition for solving
!* a set of linear equations. The matrix is input in A, and
!* the output matrices (upper and lower triangle) are output in the
!* same array, the diagonal containing elements of the upper (the diagonal
!* element of the lower are all ones.)
!* N is 'nang', the dimension of the square matrix. INDX is the output
!* permutation vector for row interconversions. D is 1 or -1 on output
!* and is only used for calculating the determinant later (which this
!* program does not do, so we throw it away).
!*
!* From Numerical Recipes in C, 2nd Ed. pp.46-47
!*
        real A(N,N), D, TINY
        parameter(TINY = 1.0E-20)
        integer N, INDX(N)
        integer i,imax,j,k
        real big,dum,sum,temp,vv(2000)
!*
        D = 1.0
        do i=1,N
            big = 0.0
            do j=1,N
                temp = abs(A(j,i))
                 if (temp.gt.big) big = temp
            enddo
            if (big.eq.0.0) stop ' Singular matrix in LUDCMP'
            vv(i) = 1.0/big
        enddo
        do j=1,N
            do i=1,j-1
                sum = A(j,i)
                 do k=1,i-1
                     sum = sum - A(k,i)*A(j,k)
                 enddo
                 A(j,i) = sum
            enddo
            big = 0.0
            do i=j,N
                sum = A(j,i)
                 do k=1,j-1
                     sum = sum -A(k,i)*A(j,k)
                 enddo
                 A(j,i) = sum
                 dum = vv(i)*abs(sum)
                 if (dum.ge.big) then
                     big = dum
                     imax = i
                 endif
            enddo
            if (j.ne.imax) then
                do k=1,N
                    dum = A(k,imax)
                    A(k,imax) = A(k,j)
                    A(k,j) = dum
                enddo
                d = -d
                vv(imax) = vv(j)
            endif
            INDX(j) = imax
            if (A(j,j).eq.0.0) A(j,j) = TINY
            if (j.ne.N) then
                dum = 1.0/A(j,j)
                do i=j+1,n
                    A(j,i) = A(j,i)*dum
                enddo
            endif
        enddo
        return
        end subroutine LUDCMP
!***********************************************************************
        subroutine LUBKSB(A,N,INDX,B)
        implicit none
!* This is the second part of the PU decomposition,
!* in which the solution vector is calculated by
!* back-substitution. A is the matrix output by LUDCMP,
!* N is 'nang', INDX is the permution vector output by
!* LUDCMP, B is input as the right-hand side vector 'tau'
!* and output as the solution vector.
        real A(N,N), B(N)
        integer N,INDX(N)
        integer i,ii,ip,j
        real sum
!*
        ii = 0
        do i=1,N
            ip = INDX(i)
            sum = B(ip)
            B(ip) = B(i)
            if (ii.ne.0) then
                do j=ii,i-1
                     sum = sum - A(j,i)*B(j)
                 enddo
            else
                if (sum.ne.0.0) ii=i
            endif
            B(i) = sum
        enddo
        do i=N,1,-1
            sum = B(i)
            do j=i+1,N
                sum = sum - A(j,i)*B(j)
            enddo
            B(i) = sum/A(i,i)
        enddo
        return
        end subroutine LUBKSB
!***********************************************************************
        real function MEASURE(a1,a2,a3,a4)
!!* This function returns the dihedral angle (degrees) between
!!* planes a1,a2,a3 and a2,a3,a4, using the right-handed-is-positive
!!* convention. A positive angle rotates a3,a4 to align with a2,a1.
        implicit none
!***** local variables
        real a1(3),a2(3),a3(3),a4(3)
        real v23(3),v21(3),v34(3),vy(3),vx(3),uy(3),ux(3),cx,cy,cc
        real,parameter :: deg=(180./pi)
!*****
        !! call subvec(a3,a2,v23)
        v23 = a3 - a2
        !! call subvec(a1,a2,v21)
        v21 = a1 - a2
        call cros(v23,v21,vy)
        call unitvec(vy,uy)
        call cros(uy,v23,vx)
        call unitvec(vx,ux)
        !! call subvec(a4,a3,v34)
        v34 = a4 - a3
        cx = dotprod(v34,ux)
        cy = dotprod(v34,uy)
        !! cc = atan2d(cy,cx)
        cc = deg*atan2(cy,cx)
        MEASURE = cc
        return
        end function MEASURE
!***********************************************************************
        real(kind=kind_8) function MEASUREd(a1,a2,a3,a4)
!!* This function returns the dihedral angle (degrees) between
!!* planes a1,a2,a3 and a2,a3,a4, using the right-handed-is-positive
!!* convention. A positive angle rotates a3,a4 to align with a2,a1.
!!** double precision **
        implicit none
!***** local variables
        real(kind=kind_8),intent(in) :: a1(3),a2(3),a3(3),a4(3)
        real(kind=kind_8) :: dd
        real(kind=kind_8) :: v23(3),v21(3),v34(3),vy(3),vx(3),uy(3),ux(3),cx,cy,cc
        real(kind=kind_8),parameter :: deg=(180./pi)
!*****
        call Dsubvec(a3,a2,v23)
        call Dsubvec(a1,a2,v21)
        call cross(v23,v21,vy,dd)
        call unitvecD(vy,uy)
        call cross(uy,v23,vx,dd)
        call unitvecD(vx,ux)
        call Dsubvec(a4,a3,v34)
        cx = Ddotprod(v34,ux)
        cy = Ddotprod(v34,uy)
        !! cc = datan2d(cy,cx)
        cc = deg*atan2(cy,cx)
        MEASUREd = cc
        return
        end function MEASUREd
!***********************************************************************
        subroutine CHOLDC(A,N,P,izero)
        implicit none
!* The first part of the Cholesky Decomposition  (used in protean)
!* Numerical Recipes, p.97
!* A = matrix
!* N = size
!* P = diagonal
        integer N
        real A(N,N), P(N)
        integer i,j,k,izero
        real sum,VerySmall
        parameter(VerySmall=-0.1E-10)
!*
        izero = 0
        do i=1,N
          do j=i,N
            sum = A(j,i)
            do k=i-1,1,-1
              sum = sum - A(k,i)*A(k,j)
            enddo
            if (i.eq.j) then
              if (sum.le.0.0000) then
                if (sum.lt.0.0000) then
                  write(*,*) '  CHOLDC failed on row= ',i
                  write(*,*) '  Sum= ',sum
                  write(*,*) '  A(i,i) before summing= ',A(j,i)
                  !! if (sum.lt.VerySmall) stop 'CHOLDC failed!! '
                endif
                izero = i
                return 
              endif
              P(i) = sqrt(sum)
            else
              A(i,j) = sum/P(i)
            endif
          enddo
        enddo
        return
        end subroutine CHOLDC
!*****************************************************
        subroutine CHOLSL(A,N,P,B,X)
        implicit none
!* The second part of the Cholesky Decomposition
!* See Numerical Recipes, p.97
!* A is from CHOLDC, N is size, P is from CHOLDC
!* B is known, X is the output vector 
        integer N
        real A(N,N), P(N), B(N), X(N)
        real sum
        integer i,k
!*
        do i=1,N
          sum = B(i)
          do k=i-1,1,-1
            sum = sum - A(k,i)*X(k)
          enddo
          X(i) = sum/P(i)
        enddo
        do i=N,1,-1
          sum = X(i)
          do k=i+1,N
            sum = sum - A(i,k)*X(k)
          enddo
          X(i) = sum/P(i)
        enddo
        return
        end subroutine CHOLSL
!***********************************************************************
        subroutine getrotS(a1,a2,chi,mat,vec)
        implicit none
!* find the matrix (mat) and vector (vec) that
!* rotate chi rad about an axis defined by a1-->a2
!* rotation is right-handed. That is, it is a
!* clockwise rotation when looking from a1 to a2
!* chi in radians
        real a1(3),a2(3),vec(3)
        real mat(3,3)
        real a(3),c(3),da,phi,psi,chi,rad
        integer :: i,j,k
!*
        if (abs(chi).lt.1.0E-4) goto 100
        phi=0.0
        psi=0.0
!* get difference vector: a
        !!  call subvec(a2,a1,a)
        a = a2 - a1
!* get length of vector A
        da=a(1)*a(1) +a(2)*a(2) +a(3)*a(3) 
        da= sqrt(da)
!* get phi, psi
        psi= acos(a(3)/da)
        if (a(1).ne.0.0.or.a(2).ne.0.0) phi= atan2(-a(1),a(2))
!* get matrix
        call getmat_S(phi,psi,chi,mat)
!* get vector, = -Mat*a1 + a1
        call rotate_S(mat,a1,c)
        call subvec(a1,c,vec)
        return
 100    do I=1,3
          do K=1,3
            mat(K,I)=0.0
          enddo
          mat(I,I)=1.0
          vec(I)=0.0
        enddo
        return
        end subroutine getrotS
!***********************************************************************
        subroutine getmat_S(phi,psi,kappa,aa)
        implicit none
!* get matrix for Tanaka convention polar angles
!* CB  14-JAN-1991, subroutine 3-JUN-1991
!* >>> phi,psi,kappa in radians <<<
        real phi,psi,kappa,sf,cf,ss,cs,sk,ck,aa(3,3)
        integer :: i,j,k
!***
!cb Try to avoid underflows...
        if (abs(phi).lt.1.0E-4) phi = 0.0
        if (abs(psi).lt.1.0E-4) psi = 0.0
        if (abs(kappa).lt.1.0E-4) then
          aa = 0.
          do I=1,3
            aa(I,I)=1.0
          enddo
          return
        endif
        sf = sin(phi)
        cf = cos(phi)
        ss = sin(psi)
        cs = cos(psi)
        sk = sin(-kappa)
        ck = cos(-kappa)
!* now calculate the product matrix of the five rotation matrices
!! not yet revised...
        aa(1,1)=(cf*cf+sf*sf*cs*cs)*ck+sf*sf*ss*ss
        aa(1,2)=ss*ss*sf*cf*(ck-1.0)-cs*sk
        aa(1,3)=sf*cs*ss*(ck-1.0)+cf*sk*ss
        aa(2,1)=ss*ss*cf*sf*(ck-1.0)+cs*sk
        aa(2,2)=(sf*sf+cf*cf*cs*cs)*ck+cf*cf*ss*ss
        aa(2,3)=cf*ss*cs*(1.0-ck)+sf*ss*sk
        aa(3,1)=cs*ss*sf*(ck-1.0)-cf*ss*sk
        aa(3,2)=cf*ss*cs*(1.0-ck)-sf*ss*sk
        aa(3,3)=ss*ss*ck+cs*cs
!* done.
!cb      write(*,'(/,3(3f12.8,/))') ((aa(I,J),I=1,3),J=1,3)
        return
        end subroutine getmat_S
!***********************************************************************
        real function mytanh(dd)
        implicit none
        real dd,y,ovrflo
        parameter(ovrflo=5.)
        if (dd.gt.ovrflo) then
          y = 1.
        elseif (dd.lt.-ovrflo) then
          y = -1.
        else
          y = tanh(dd)
        endif
        mytanh = y
        return
        end function mytanh
!***********************************************************************
        subroutine switchvec(v1,v2)
        implicit none
        real v1(3),v2(3),v3(3)
!* 
        v3(1) = v1(1) 
        v3(2) = v1(2)
        v3(3) = v1(3)
        v1(1) = v2(1)
        v1(2) = v2(2)
        v1(3) = v2(3)
        v2(1) = v3(1)
        v2(2) = v3(2)
        v2(3) = v3(3)
        return
        end subroutine switchvec
!***********************************************************************
        real function chirality(v0,v1,v2,v3)
        implicit none
        real v0(3),v1(3),v2(3),v3(3)
        real v01(3),v02(3),v03(3),v12(3)
        call subvec(v1,v0,v01)
        call subvec(v2,v0,v02)
        call cros(v01,v02,v12)
        call subvec(v3,v0,v03)
        chirality = dotprod(v12,v03)
        return
        end function chirality
!***********************************************************************
        real function dist(v0,v1)
        implicit none
        real v0(3),v1(3),vv(3)
        call subvec(v1,v0,vv)
        dist = sqrt(dotprod(vv,vv))
        return
        end function dist
!***********************************************************************
        real(kind=kind_8) function Ddsqr(v0,v1)
        implicit none
        real(kind=kind_8) :: v0(3),v1(3),vv(3)
        vv = v1 - v0
        Ddsqr = sum(vv*vv)
        end function Ddsqr
!***********************************************************************
        real(kind=kind_8) function Ddist(v0,v1)
        implicit none
        real(kind=kind_8) v0(3),v1(3),vv(3)
        !call Dsubvec(v1,v0,vv)
        vv = v1 - v0
        Ddist = dsqrt(sum(vv*vv))
        return
        end function Ddist
!***********************************************************************
        real function diffdist(v0,v1,v2,v3)
        implicit none
        real v0(3),v1(3),v2(3),v3(3),x,y
        real v01(3),v23(3)
        call subvec(v1,v0,v01)
        call subvec(v3,v2,v23)
        x = sqrt(dotprod(v01,v01))
        y = sqrt(dotprod(v23,v23))
        diffdist = x - y
        return
        end function diffdist
!***********************************************************************
        real(kind=kind_8) function dDdX(ri,rh,rk,rj)
        implicit none
!* given four atom positions, a1 a2 a3 a4 ,
!* this routine finds the slope of the change
!* in the distance between ri and rj with rotation
!* about an axis defined by rh and rk. 'dDdX' is the
!* slope in angstrom/radian.
!* A "positive" rotation is a right-handed rotation.
        real(kind=kind_8),intent(in) :: ri(3),rh(3),rk(3),rj(3)
        real(kind=kind_8) :: v1(3),vij(3),uij(3),vhk(3),uhk(3),vkj(3)
        real(kind=kind_8) :: dd
!**************
        ! call Dsubvec(rj,rk,vkj)
        vkj = rj - rk
        ! call Dsubvec(rj,ri,vij)
        vij = rj - ri
        ! call Dsubvec(rk,rh,vhk)
        vhk = rk - rh
        ! call unitvecD(vij,uij)
        uij = vij/sqrt(vij(1)**2 + vij(2)**2 + vij(3)**2)
        ! call unitvecD(vhk,uhk)
        uhk = vhk/sqrt(vhk(1)**2 + vhk(2)**2 + vhk(3)**2)
        call cross(uhk,vkj,v1,dd)
        ! dd = Ddotprod(uij,v1)
        dd = uij(1)*v1(1)+uij(2)*v1(2)+uij(3)*v1(3)
        dDdX = dd
        end function dDdX
!***********************************************************************
        subroutine getframe(vi,vj,vk,mat,vec)
        implicit none
!!
!!      The matrix and vector passed back by this routine
!!      will rotate and translate a point into the global reference frame
!!      from a local reference frame defined by the three atoms vi, vj, vk. 
!!      In the local frame, the X-axis is defined as the vi->vj and the
!!      XY plane is the plane of the three points with vk on the
!!      positive side (i.e. a R-handed frame).
!!      CB Mon Aug 13 17:13:47 EST 2001
!!
        real vi(3),vj(3),vk(3),mat(3,3),vec(3)
        real vij(3),vik(3),vvk(3),vv(3),m1(3,3),m2(3,3),m3(3,3)
        real jxk(3),vvj(3)
        real aa,dj,phi,psi,chi
        integer i,j,k
        
        ! write(*,*) 'IN new GETFRAME'
        !-------------------- new getframe -----------------------!
        vij = vj - vi
        vik = vk - vi
        call cros(vij,vik,jxk)
        phi = atan2(jxk(2),jxk(1))  !! rotate phi around Z
        ! write(*,*) "y, x, phi", jxk(2),jxk(1), phi
        dj = sqrt(dotprod(jxk,jxk))
        psi = acos(jxk(3)/dj)       !! rotate psi around Y
        ! write(*,*) "z, psi, d", jxk(3), psi, dj
        m1 = reshape((/cos(phi), sin(phi),  0.,  &
                      -sin(phi), cos(phi),  0.,  &
                       0.,        0.,       1./),(/3,3/))        
         ! write(*,'(3f10.5)') m1(1,1),m1(2,1),m1(3,1)
         ! write(*,'(3f10.5)') m1(1,2),m1(2,2),m1(3,2)
         ! write(*,'(3f10.5)') m1(1,3),m1(2,3),m1(3,3)
        ! call rotate_S(m1,jxk,vvk)
        ! write(*,'("Should be in XZ plane:",3f10.5)') vvk(1:3)
        m2 = reshape((/cos(psi),  0., -sin(psi),  &
                        0.,       1.,       0.,  &
                       sin(psi), 0.,  cos(psi)/),(/3,3/))        
        call MM_S(m1,m2,m3)
        ! call rotate_S(m2,vvk,vvj)
        ! ! write(*,*) "d rotated m1m2=", sqrt(dotprod(vvj,vvj))
        ! write(*,'("Should be Z-axis     :",3f10.5)') vvj(1:3)
        ! call rotate_S(m3,jxk,vvk)
        ! ! write(*,'("This is the cros-prod:",3f10.5)') jxk(1:3)
        ! ! write(*,'(3f10.5)') m3(1,1),m3(2,1),m3(3,1)
        ! ! write(*,'(3f10.5)') m3(1,2),m3(2,2),m3(3,2)
        ! ! write(*,'(3f10.5)') m3(1,3),m3(2,3),m3(3,3)
        ! ! write(*,*) "d rotated m3=", sqrt(dotprod(vvk,vvk))
        ! write(*,'("This should be Z-axis:",3f10.5)') vvk(1:3)
        call rotate_S(m3,vij,vvj)
        chi = atan2(vvj(2),vvj(1))  !! rotate chi around Z
        m1 = reshape((/cos(chi), sin(chi),  0.,  &
                      -sin(chi), cos(chi),  0.,  &
                       0.,        0.,       1./),(/3,3/))        
        call MM_S(m3,m1,m2)
        !call rotate_S(m2,vij,vvk)
        !write(*,'("This should be x-axis:",3f10.5)') vvk(1:3)
        !call rotate_S(m2,vik,vvk)
        !write(*,'("This in XY-plane     :",3f10.5)') vvk(1:3)
        mat = transpose(m2)        !! remove this line to get the inverse operation
        vec = vi
        end subroutine getframe
!***********************************************************************
        subroutine cosinerule(a,b,c,alph,flag)
        implicit none
        !! returns angle opposite c (between a and b)
        !! NOTE: distances a,b,c cannot be zero or negative
        real,parameter :: pi=3.1415927
        real,intent(out) :: alph
        real,intent(in) :: a,b,c
        real :: x
        integer :: flag
        x = ( b*b + a*a - c*c)/(2*a*b)
        if (x >= 1.0) then
          alph = 0.00    !! bad triangle, c too short
          flag = 1
        elseif (x < (-1.0)) then
          alph = pi      !! bad triangle, c too long
          flag = 1
        else
          alph = acos(x)
          flag = 0
        endif
        end subroutine cosinerule
!***********************************************************************
        real function sphericalexcess(a,b,c)
        implicit none
        !! Find the "spherical excess", whic is the sum of angles
        !! of a spherical trienagle minus pi (the sum of angles for a triangle)
        !! Spherical excess times radius^2 is the area of a spherical
        !! triangle.
        !! Uses l'Huilier's Theorem. From mathworld.wolfram.com/SphericalExcess.html
        !! Input: spherical triangle sides: a,b,c (radians)
        !! Output: spherical excess = A + B + C - pi, where A,B,C are the angles 
        !!    of a spherical triangle, not to be confused with a,b,c, the sides of
        !!    the spherical triangle in radians. A,B,C are the angles made by tangents
        !!    to the triangle verteces.  a,b,c are the angles made by the verteces
        !!    with respect to the center of the sphere.
        !!  Sat Dec 27 11:02:32 EST 2003  C.Bystroff
        real,intent(in) :: a,b,c !! radians
        real :: x, sp  !! semiperimeter
        sp = ( a + b + c )/2
        x = tan(sp/2)
        x = x * tan((sp-a)/2)
        x = x * tan((sp-b)/2)
        x = x * tan((sp-c)/2)
        x = sqrt(x)
        sphericalexcess = 4*atan(x)
        end function sphericalexcess
!***********************************************************************
        subroutine heronsrule(a,b,c,ar,alph)
        implicit none
        !! returns area of a triangle and the angle between the first two sides
        !! NOTE: distances a,b,c cannot be zero or negative
        real a,b,c,ar,alph,s,xh
        real,parameter :: pi=3.1415927
        s = (a+b+c)/2.
        ar = sqrt(s*(s-a)*(s-b)*(s-c))
        xh = 2*ar/a
        alph = (pi/2.) - acos(xh/b)    !! BUG. This angle is always acute
        end subroutine heronsrule
!***********************************************************************
        subroutine getdme(c0,c1,dme,natm)
        implicit none
        integer,intent(in) ::  natm
        real,intent(out) :: dme
        real,dimension(3,natm),intent(in):: c0,c1
        real :: x,y
        integer :: I,J,nlocal=3,nn
        !!
        x = 0
        nn = 0
        do I=1,natm-nlocal
          do J=I+nlocal,natm
            x = dist(c0(1:3,I),c0(1:3,J))
            y = dist(c1(1:3,I),c1(1:3,J))
            nn = nn + 1
            dme = dme + (x-y)**2
          enddo
        enddo
        dme = sqrt(dme/real(nn))
        end subroutine getdme
!***********************************************************************
        subroutine ovrlap(c0,c1,err,natm,mtx,vec)
        !!  22-JUN-03
        !** coord0 is template.
        !** coord1 is target, returned
        !** err is rmsd of natm atoms
        !** get the rms deviation between two arbitrarilly oriented
        !** sets of 'natm' atoms.
        implicit none
        integer,intent(in) ::  natm
        real,intent(out) :: err,mtx(3,3),vec(3)
        real,dimension(3,natm),intent(in):: c0,c1
        real,dimension(:,:),allocatable :: x0,x1
        real(kind=kind_8) :: aa(3,3),rot(3,3),t(3),sig,sg,gam,bet,cc,bb
        real ::  mat(3,3),xc0(3),xc1(3)
        integer :: ix,iy,iz,i,j,k,iflag,ict,nnn,ios
        real :: y
        real,parameter :: tol=0.00001
!*
        !! if (natm>1000) stop 'Increase memory for ovrlap'
        if (allocated(x0)) deallocate(x0)
        allocate(x0(3,natm),stat=ios)
        if (ios/=0) stop 'vectormath:: ovrlap: allocation failed for x0'
        if (allocated(x1)) deallocate(x1)
        allocate(x1(3,natm),stat=ios)
        if (ios/=0) stop 'vectormath:: ovrlap: allocation failed for x1'
        xc0 = 0.
        xc1 = 0.
        do I=1,natm
          xc0 = xc0 + c0(1:3,I)
          xc1 = xc1 + c1(1:3,I)
        enddo
        if (natm < 5) then
          write(*,'("VECTORMATH::ovrlap: Number of atoms   =",i9,"  <===too small.")') natm
          err = -1
          return
        endif
        xc0 = xc0/real(natm)   !! center of mass, target
        xc1 = xc1/real(natm)   !! center of mass, template
        do I=1,natm
          x0(:,I) = c0(:,I) - xc0
          x1(:,I) = c1(:,I) - xc1
        enddo
        !* x0 and x1 are centered on the origin
        aa = 0.
        do i=1,3
          do j=1,3
            aa(i,j) = dot_product(x1(i,:),x0(j,:))
          enddo
        enddo
        rot = 0.
        do i=1,3
          rot(i,i) = 1.
        enddo
        !! new f90 version
        ict=0
        iflag=0
        ix=1
        do while (ict < 1000)
          ict=ict+1
          do ix=1,3
            iy=mod(ix,3)+1
            iz=6-ix-iy
            sig=aa(iz,iy)-aa(iy,iz)
            gam=aa(iy,iy)+aa(iz,iz)
            sg=sqrt(sig*sig+gam*gam)
            if(sg == 0.) cycle
            sg=1./sg
            if(abs(sig) <= tol*abs(gam)) cycle
            do k=1,3
              bb=gam*aa(iy,k)+sig*aa(iz,k)
              cc=gam*aa(iz,k)-sig*aa(iy,k)
              aa(iy,k)=bb*sg
              aa(iz,k)=cc*sg
              bb=gam*rot(iy,k)+sig*rot(iz,k)
              cc=gam*rot(iz,k)-sig*rot(iy,k)
              rot(iy,k)=bb*sg
              rot(iz,k)=cc*sg
            enddo
            iflag=1
          enddo
          if (iflag==0) exit
        enddo
        err=0.
        mtx = sngl(rot)
        do j=1,3
          t(j) = dot_product(rot(j,:),xc1(:))
        enddo
        vec = sngl(t)
        vec = xc0 - vec
        do I=1,natm
          do j=1,3
            t(j)=dot_product(rot(j,:),x1(:,I))
          enddo
          x1(:,I)= t + xc0  !! x1 is now superimposed on x0
          t = c0(:,I) - x1(:,I)
          err = err + dot_product(t,t)
        enddo
        err=sqrt(err/natm)
        ! do I=1,natm
        !   write(*,'("Superimposed ATOMS ",i5,6f8.3)') &
        !   I,c0(1:3,I),x1(1:3,I)
        ! enddo
        !write(*,*) 'RMSD===============',err
        deallocate(x0,x1)
        end subroutine ovrlap
!***********************************************************************
        subroutine superimp(coord0,coord1,err,natm)
        !!  06-SEP-02
        !** coord0 is template.
        !** coord1 is target, returned
        !** err is rmsd of natm atoms
        !** get the rms deviation between two arbitrarilly oriented
        !** sets of 'natm' atoms.
        implicit none
        integer,intent(in) ::  natm
        real,intent(out) :: err
        real,dimension(3,natm),intent(in):: coord0
        real,dimension(3,natm),intent(inout):: coord1
        real,dimension(:,:),allocatable :: x0,x1
        real ::  mat(3,3), vec(3),t(3)
        integer :: I,nnn,j,k
!*
        allocate(x0(3,natm),x1(3,natm),stat=i)
        !! write(*,'("Number of atoms to allocate=",i9)') natm
        if (i/=0) then
          stop 'Failed to allocate space in SUPERIMP'
        endif
! note: if any atom is missing, omit it from the calculation
!       if less than 5 atoms left, return err=-1 as a flag
        nnn = 0
        do I=1,natm
          if (coord0(1,I) > 998.0) cycle
          if (coord1(1,I) > 998.0) cycle
          nnn = nnn + 1
          x0(1:3,nnn) = coord0(1:3,I)
          x1(1:3,nnn) = coord1(1:3,I)
          !! write(*,'(6f8.3)') x0(1:3,nnn),x1(1:3,nnn)
        enddo
        if (nnn < 5) then
          !! write(*,'("VECTORMATH::superimp:Number of atoms   =",i9,"  <===too small.")') nnn
          err = -1
          deallocate(x0,x1)
          return
        endif
        call ovrlap(x0,x1,err,nnn,mat,vec)
          !!  write(*,'("ERR=",f9.5)') err
          ! write(*,'("nnn=",i9)') nnn
        if (err < 0.) write(0,'("VECTORMATH::superimp: ovrlap failed.")')
        do I=1,natm
          if (coord1(1,I) > 998.0) cycle
          do j=1,3
            t(j) = 0.
            do k=1,3
              t(j) = t(j) + mat(j,k)*coord1(k,I)
            enddo
          enddo
          do k=1,3
            coord1(k,I) = t(k) + vec(k)
          enddo
          ! write(*,'("Superimposed ATOMS ",i5,6f8.3)') &
          ! I,coord0(1:3,I),coord1(1:3,I)
        enddo
        deallocate(x0,x1)
        end subroutine superimp
!!***********************************************************************
!* find the limits of up to N fields of non-blank
!* characters in the string 'aline' . The start of
!* field I is A(I), the end is B(I)
        integer function getfields(N,aline,A,B)
        implicit none
        integer,intent(in) ::  N
        integer,intent(inout) ::  A(N),B(N)
        integer ::  I,FF,LL
        character(len=*),intent(in) :: aline
!**
        I = 1
        FF = 0
        LL = len(aline)
        do while (FF.lt.N.and.I.le.LL)
          do while (aline(I:I).eq.' '.and.I.le.LL)
            I = I + 1
          enddo
          if (I.le.LL) then
            FF = FF+1
            A(FF) = I
            do while (aline(I:I).ne.' '.and.I.le.LL)
              I = I + 1
            enddo
            B(FF) = I-1
          endif
        enddo
        getfields = FF
        end function getfields
!************************************************************************
        real function TmplFrc(dd,tt)
          implicit none
          real dd,tt,x,y
          x = tt-dd
          if (x.gt.8.) then
            y = 0
          elseif (x.lt.-8.) then
            y = 0
          else
            y = tanh(x)*exp(-(x/4.)**2)
          endif
          TmplFrc = y
          !!write(*,'("TmplFrc=",f9.4)') y
        end function TmplFrc
!************************************************************************
!!======================== Reentrant surface derivatives =============================
!!
!! The gradient for reentrant surfaces depends on just the three atoms forming the
!! base of the tetrahedron with a water. If there is no intersecting reentrant surface,
!! the derivitives are just the pairwise distance-based derivatives of the 
!! spherical triangle (strianglederiv). If part of the triangle is blocked by
!! intersecting surfaces, then the triangle is truncated by the base of the
!! tetrahedron. The surface lost is the spherical cap below the ijk plane if no part
!! of the circle made by the ijk plane and the sphere is outside of the spherical
!! triangle. Then the derivative is the sum of the derivatives of the triangle
!! and the cap (dCapdDij). If the circle extends beyond the boundaries of the spherical
!! triangle, the amount of the cap that is within the triangle is estimated
!! numerically by subtracting the masked triangle (masked by intersecting
!! probes) from the triangle. The difference is the part of the triangle
!! lost upon intersection. The derivative in this case is still the sum of the
!! derivatives of the triangle and the truncated spherical cap. The derivative
!! of the truncated cap is the sum of the derivative of the untruncated cap (dCapdDij)
!! minus the derivatives of all spherical zones or frustums that fall outside
!! the triangle. This is hard to calculate, and might lead to bugs due to
!! the need to check for frustums. Instead, we can assume that the derivative 
!! is proportional to the portion
!! of the cap's rim that falls inside the triangle. This is the arc fraction that can
!! be calculated with one additional masking step. It is not likely to make the program
!! much slower since it is only calculated in cases where an intersecting reentrant surface exists.
!! In short, we calculate the change in the size of the triangle analytically,
!! then check for intersecting, then subtract the change in the intersecting surface,
!! which is the change in area of a spherical cap times the fraction of the cap rim exposed.
!!
!!==================================================================================

subroutine dCapdDij(abur,rw,dij,dik,djk,diw,djw,dkw,dAdij,dAdjk,dAdik)
  implicit none
  !!
  !! This is dA/d(dij)  where A is the area of a spherical cap of radius rw
  !! where the base is defined by the distances dij,dik,djk,diw,djw,dkw for a tetrahedron
  !! with w centered in the sphere. There is no guarantee that the base (ijk) of the 
  !! tetrahedron intersects the sphere. If it does not, then this number is NaN.
  !! NOTE: the equation for the height of the cap is the same as the
  !! equation for the z-coordinate (d) of the water in the reference frame
  !! of the tree-atom base, which was calculated for the routine tetrahedron() in
  !! masker.f90   Here is an excerpt of that routine:
  !!----
  !! x = diw*cos(ajiw)
  !! yp = diw*cos(akiw)
  !! y = yp/sin(ajik) - x/tan(ajik)
  !! d = sqrt(diw*diw - x*x - y*y)
  !!----
  !! It is the equation for d that we differenciate with respect to the ij edge, below.
  !! Sorry for the long equation. (to be simplified...)
  !! 30-DEC-03  C.Bystroff
  !!-----------------------------------------
  !! The height d is also:
  !! dotprod(viw, cross(vij,vik)) / abs(cross(vij,vik))
  !!----------------------------------------
  real,intent(inout) :: abur
  real,intent(in) :: rw,dij,dik,djk,diw,djw,dkw
  real,intent(out) :: dAdij,dAdik,dAdjk
  real :: ddd(3),dd2(3),dww(3),dw2(3)
  real :: perim, dddikmj, dddijmk, dddjkmi, dAdd(3), acap, dddabmc(3)
  integer :: ijk,i,j,k

  !! rw,pi are global
  !!
  !! abur is A^2 of buried surface in the cap region
  !! divide this by the total area of the cap to get the fraction that is variable
  !! with respect to the distances..
  !! Multiply this fraction by the full cap deriv and return.
  !! ------ Area of cap : from Mathematica, using eq for d in tetrahedron()
  !! 2*Pi*rw*(-Sqrt[diw^2 - (-dij^2 - diw^2 + djw^2)^2/(4*dij^2) - 
  !!                diw^2*(-(((-dij^2 - dik^2 + djk^2)*(-dij^2 - diw^2 + djw^2))/
  !!          (4*dij^2*dik*diw*Sqrt[1 - (-dij^2 - dik^2 + djk^2)^2/(4*dij^2*
  !!               dik^2)])) - (-dik^2 - diw^2 + dkw^2)/(2*dik*diw*
  !!          Sqrt[1 - (-dij^2 - dik^2 + djk^2)^2/(4*dij^2*dik^2)]))^2] + rw)
  !!  ---------  Simplified form:
  !! 2*Pi*rw*
  !! (-Sqrt[(dik^4*djw^2 + dij^4*dkw^2 + dij^2*((dik - diw)*(dik + diw)*
  !!       (djk - djw)*(djk + djw) - (dik^2 + diw^2 + djk^2 + djw^2)*dkw^2 + 
  !!      dkw^4) + djk^2*(diw^4 + djw^2*dkw^2 + diw^2*(djk^2 - djw^2 - dkw^2)) - 
  !!   dik^2*(diw^2*(djk^2 + djw^2 - dkw^2) + djw^2*(djk^2 - djw^2 + dkw^2)))/
  !!   ((dij - dik - djk)*(dij + dik - djk)*(dij - dik + djk)*
  !!    (dij + dik + djk))] + rw)
  !! 
  ddd = (/dij,djk,dik/)
  dd2 = ddd*ddd  !! =(/dij*dij,djk*djk,dik*dik/)
  dww = (/diw,djw,dkw/)
  dw2 = dww*dww  !! =(/diw*diw,djw*djw,dkw*dkw/)
  perim = dij + dik + djk
  ijk = 1
  i = ijk; j=mod(ijk,3)+1; k=mod(ijk+1,3)+1
  dddijmk = ddd(i) + ddd(j) - ddd(k)
  dddjkmi = ddd(j) + ddd(k) - ddd(i)
  dddikmj = ddd(k) + ddd(i) - ddd(j)
  dddabmc = (/dddijmk,dddjkmi,dddikmj/)
  !! acap is the total area of the spherical cap
  !acap = 2*pi*rw*                                                                               &  
  !(-sqrt(((dd2(3)**2)*dw2(2) + (dd2(1)**2)*dw2(3) + dd2(1)*((dik - diw)*(dik + diw)*                  &
  !        (djk - djw)*(djk + djw) - (dd2(3) + dw2(1) + dd2(2) + dw2(2))*dw2(3) +                &
  !         dw2(3)**2) + dd2(2)*(dw2(1)**2 + dw2(2)*dw2(3) + dw2(1)*(dd2(2) - dw2(2) - dw2(3))) -    &
  !         dd2(3)*(dw2(1)*(dd2(2) + dw2(2) - dw2(3)) + dw2(2)*(dd2(2) - dw2(2) + dw2(3))))/     &
  !        ((dij - dik - djk)*(dddikmj)*(dddijmk)*(perim))) + rw)
  !!
  !!  ijk=1:  1,2,3,  ijk=2: 2,3,1,  ijk=3:  3,1,2
  do ijk=1,3
    i = ijk; j=mod(ijk,3)+1; k=mod(ijk+1,3)+1
    dddijmk = dddabmc(i)
    dddjkmi = dddabmc(j)
    dddikmj = dddabmc(k)
    !!================== SIMPLIFY THIS ======================!!
    dAdd(ijk) =                                                         &
    ( 2*ddd(i)*pi*rw*                                                   &
     ( dd2(j)*(dw2(i)-dw2(k))-                                          &
       dd2(k)*dd2(k) +                                                  &
       dd2(i)*(dd2(k)-dw2(i)+dw2(k))+                                   &
       dd2(k)*(dw2(i)+dd2(j)-2*dw2(j)+dw2(k)) )*                        &
     ( dd2(k)*(dd2(j)+dw2(j)-dw2(k)) +                                  &
       dd2(i)*(dd2(j)-dw2(j)+dw2(k)) +                                  &
       dd2(j)*(dw2(j)+dw2(k)-2*dw2(i)-dd2(j)) )                         &
    )/                                                                  &
    ( dddikmj*dddikmj*dddijmk*dddijmk*dddjkmi*dddjkmi*perim*perim*      &
      sqrt( ( dd2(k)*dd2(k)*dw2(j) + dd2(i)*dd2(i)*dw2(k) +             &
              dd2(i)*( (dd2(k)-dw2(i))*(dd2(j)-dw2(j)) -                &
                       (dd2(k)+dw2(i)+dd2(j)+dw2(j)+dw2(k))*dw2(k) ) +  &
              dd2(j)*( dw2(i)*dw2(i) + dw2(j)*dw2(k) +                  &
                       dw2(i)*(dd2(j)-dw2(j)-dw2(k)) ) -                &
              dd2(k)*( dw2(i)*(dd2(j)+dw2(j)-dw2(k)) +                  &
                       dw2(j)*(dd2(j)-dw2(j)+dw2(k)) )                  &
            )/( (ddd(i)-ddd(k)-ddd(j))*dddikmj*dddijmk*perim )          &
          )                                                             &
    )
  enddo
  dAdd = abur*dAdd
  dAdij = dAdd(1)
  dAdjk = dAdd(2)
  dAdik = dAdd(3)
  !! write(0,'("buried fraction=",f10.4)') (abur/acap)
  !! abur = abur/acap
end subroutine dcapddij

!! ========= The following is d(Acap)/dDij pasted directly from Mathematica. In terms of the six distances.
!! ========= This has been factored and slightly rearranged to get the code above.
!! (2*dij*(-dik^4 + djk^2*(diw - dkw)*(diw + dkw) + 
!!     dij^2*(dik^2 - diw^2 + dkw^2) + dik^2*(diw^2 + djk^2 - 2*djw^2 + dkw^2))*
!!    (dik^2*(djk^2 + djw^2 - dkw^2) + dij^2*(djk^2 - djw^2 + dkw^2) + 
!!     djk^2*(-2*diw^2 - djk^2 + djw^2 + dkw^2))*Pi*rw)/
!!   ((dij + dik - djk)^2*(dij - dik + djk)^2*(-dij + dik + djk)^2*
!!    (dij + dik + djk)^2*Sqrt[(dik^4*djw^2 + dij^4*dkw^2 + 
!!       dij^2*((dik - diw)*(dik + diw)*(djk - djw)*(djk + djw) - 
!!         (dik^2 + diw^2 + djk^2 + djw^2)*dkw^2 + dkw^4) + 
!!       djk^2*(diw^4 + djw^2*dkw^2 + diw^2*(djk^2 - djw^2 - dkw^2)) - 
!!       dik^2*(diw^2*(djk^2 + djw^2 - dkw^2) + djw^2*(djk^2 - djw^2 + dkw^2)))/
!!      ((dij - dik - djk)*(dij + dik - djk)*(dij - dik + djk)*(dij + dik + djk))])
!!====================================================================================

!! ==================== dstriangle, from masker_striangle.nb on Panama.
!! This is the change in the area of a spherical triangle with repect to the
!! distance between two atoms, i and j.  To get the derivative with respect to
!! atom pairs ik and jk just substitute.
!! 
real function dAngdDij(dij,diw,djw)
  !! this is the derivative of the cosine rule w/respect to the opposite side, Dij.
  real,intent(in) :: dij,diw,djw
  real :: iw2,jw2,ij2,x
  iw2 = diw*diw
  jw2 = djw*djw
  ij2 = dij*dij
  x = 1 - ((ij2 - iw2 - jw2)**2)/(4*iw2*jw2)
  if (x<=0.) stop 'vectormath:dAngdDij: BUG x<=0.'
  x = diw*djw*sqrt(x)
  if (x<=0.) stop 'vectormath:dAngdDij: BUG sqrt(x)==0.'
  dAngdDij = dij/x
end function dAngdDij

!!
!!-----------------------------
!subroutine strianglederiv(rw,dij,dik,djk,diw,djw,dkw,ddij,ddjk,ddik)
!  implicit none
!  !! This routine calculates the pairwise distance-based partial
!  !! derivatives of the spherical triangle defined by three atoms in 
!  !! contact with a probe sphere.
!  !!  First six arguments are input distances in angstroms
!  !! global distance rw is the radius of the probe water.
!  !! The last three arguments are the output derivatives in
!  !! square angstroms surface area per angstrom change in distance.
!  !! In other words: ddij = change in spherical triangle with change in ij atom distance.
!  !! ...and so on.
!  !!
!  !! Sat Dec 27 14:49:03 EST 2003  C.Bystroff
!  real,intent(in) :: rw,dij,dik,djk,diw,djw,dkw
!  real,intent(out) :: ddij,ddjk,ddik
!  real :: ijang, ikang, jkang
!  integer :: flag
!  real :: semip , spij , spik , spjk , daangddij , dbangddik , dcangddjk , tsp2 
!  real :: tspij , tspik , tspjk , tsp2ij , tsp2ik , tsp2jk , tspikjk , tspkjij , tspijik 
!  real :: tsp2ijk , rwrw , denom , scsemip2 , scspij2 , scspik2 , scspjk2 
!  !! real :: sec
!
!  call cosinerule(diw,djw,dij,ijang,flag)  !! ijang = ArcCos[-((dij^2 - djw^2 - diw^2)/(2*diw*djw))
!  call cosinerule(diw,dkw,dik,ikang,flag)  !! ikang = ArcCos[-((dik^2 - diw^2 - dkw^2)/(2*diw*dkw))
!  call cosinerule(djw,dkw,dij,jkang,flag)  !! jkang = ArcCos[-((djk^2 - djw^2 - dkw^2)/(2*djw*dkw))
!
!  semip = (ijang + ikang + jkang ) /2.     !! radians
!  spij = (semip-ijang)/2.     !! radians
!  spik = (semip-ikang)/2.     !! radians
!  spjk = (semip-jkang)/2.     !! radians
!  daangddij = dAngdDij(dij,diw,djw)  !!  dij/(diw*djw*sqrt(1 - ((ij2 - iw2 - jw2)**2)/(4*iw2*jw2)))
!  dbangddik = dAngdDij(dik,dkw,diw)  !!   units are 1/angstroms
!  dcangddjk = dAngdDij(djk,djw,dkw)  !!  
!  tsp2 = tan(semip/2)
!  tspij = tan(spij)
!  tspik = tan(spik)
!  tspjk = tan(spjk)
!  tsp2ij = tsp2*tspij
!  tsp2ik = tsp2*tspik
!  tsp2jk = tsp2*tspjk
!  tspikjk = tspik*tspjk
!  tspkjij = tspij*tspjk
!  tspijik = tspij*tspik
!  tsp2ijk = tsp2ij*tspikjk
!  rwrw = 0.5*rw*rw       !! angstroms^2
!  denom = (sqrt(tsp2ijk)*(1+tsp2ijk))
!  write(0,'("in strianglederiv: denom: ",5f8.4)') denom,tsp2,tsp2ijk
!  scsemip2 = sec(semip/2)**2
!  scspij2 = sec(spij)**2
!  scspik2 = sec(spik)**2
!  scspjk2 = sec(spjk)**2
!  
!  !! derivative of spherical triangle area w/respect to Dij
!  !! Calculated using Mathematica.
!  !! units are angstroms^2/angstroms
!  ddij = rwrw*daangddij*( tsp2ij*(  scspjk2*tspik + scspik2*tspjk    ) -    &
!           tspikjk*( scspij2*tsp2  - scsemip2*tspij ) )/denom
!  !!
!  !! to get dDjk from dDij. i-->j, j-->k, k-->i
!  ddjk = rwrw*dcangddjk*( tsp2jk*(  scspik2*tspij + scspij2*tspij    ) -    &
!           tspijik*( scspjk2*tsp2  - scsemip2*tspjk ) )/denom
!  !!
!  !! to get dDik from dDij, i-->k, j-->i, k-->j
!  ddik = rwrw*dbangddik*( tsp2ik*(  scspij2*tspjk + scspjk2*tspij    ) -    &
!           tspkjij*( scspik2*tsp2  - scsemip2*tspik ) )/denom
!
!
!end subroutine strianglederiv
!!!---- secant function. Comment this out is it is intrinsic.
real function sec(x)
  real,intent(in) :: x
  sec = 1./cos(x)
end function sec
!!===========new strianglederiv ==========

! subroutine strianglederiv(rw,dij,dik,djk,diw,djw,dkw,ddij,ddjk,ddik)
!   implicit none
  !! Second attempt:
  !  Fri Jan  2 12:24:28 EST 2004
  !! This routine calculates the pairwise distance-based partial
  !! derivatives of the spherical triangle defined by three atoms in 
  !! contact with a probe sphere.
  !!  First six arguments are input distances in angstroms
  !! global distance rw is the radius of the probe water.
  !! The last three arguments are the output derivatives in
  !! square angstroms surface area per angstrom change in distance.
  !! In other words: ddij = change in spherical triangle area with change in ij atom distance.
  !! ...and so on.
  !!
  !! Sat Dec 27 14:49:03 EST 2003  C.Bystroff
!   real,intent(in) :: rw,dij,dik,djk,diw,djw,dkw
!   real,intent(out) :: ddij,ddjk,ddik
!   real,dimension(3) :: ddd,dww,dd2,dw2,dadd
!   real::acjkjwkw,acijiwjw,acikiwkw,sqijiwjw,acacac,ohacacac,toq3ac,tnacij
!   real :: tnacik,tnacjk,tntntn,toqtntn,dijwsq,diodij
!   integer :: ijk,i,j,k

  !! dstriangle, from masker_striangle.nb on Panama.
  !! This is the change in the area of a spherical triangle with repect to the
  !! distance between two atoms, i and j.  To get the derivative with respect to
  !! atom pairs ik and jk just substitute.
!  ddd = (/dij,djk,dik/)
!  dww = (/diw,djw,dkw/)
!  dd2 = ddd*ddd
!  dw2 = dww*dww
!  do ijk=1,3
!    i = ijk; j=mod(ijk,3)+1;  k=mod(ijk+1,3)+1
!    acjkjwkw = acos(-((dd2(j) - dw2(j) - dw2(k))/(2*dww(j)*dww(k))))
!    acijiwjw = acos(-((dd2(i) - dw2(i) - dw2(j))/(2*dww(i)*dww(j))))  
!    acikiwkw = acos(-((dd2(k) - dw2(i) - dw2(k))/(2*dww(i)*dww(k))))  
!    sqijiwjw = Sqrt(1 - ((dd2(i) - dw2(i) - dw2(j))**2)/(4*dw2(i)*dw2(j)))
!    acacac = acijiwjw + acikiwkw + acjkjwkw
!    ohacacac = 0.5*(acacac)
!    toq3ac = Tan(0.25*(acacac))
!    tnacij = Tan(0.5*(-acijiwjw + ohacacac))
!    tnacik = Tan(0.5*(-acikiwkw + ohacacac))
!    tnacjk = Tan(0.5*(-acjkjwkw + ohacacac))
!    tntntn = tnacij* tnacik* tnacjk
!    toqtntn = toq3ac* tnacij* tnacik
!    dijwsq = 4*dww(i)*dww(j)*sqijiwjw
!    diodij = 2*rw*rw*ddd(i)/dijwsq
!  
!    dadd(ijk)= diodij*(((Sec(0.5*(ohacacac-acjkjwkw))**2)*toq3ac*tnacij*tnacik)+  &
!                       ((Sec(0.5*(ohacacac-acikiwkw))**2)*toq3ac*tnacij*tnacjk) -  &
!                       ((Sec(0.5*(ohacacac-acijiwjw))**2)*toq3ac*tnacik*tnacjk) +  &
!                       (((Sec(0.25*(acacac))**2)*tntntn)))/  &
!                      (Sqrt(toq3ac*tntntn)*(1+toq3ac*tntntn))
!  enddo
!  ddij = -dadd(1)
!  ddjk = -dadd(2)
!  ddik = -dadd(3)
!! The above calculation has a bug, maybe due to the improper use of lHuiler's theorem?
!! Instead, we use the spherical excess formula. The derivative of the spjerical triangle
!! area is rw**2 times the derivative of the dihedral angle opposite the distance
!! in question. The derivative of the dihedral angle is the derivative of the 
!!  rule, but the three distance a,b,c, must be for the projection of the ijk triangle
!! on a plane perpendicular to iw.
subroutine strianglederiv(rw,vi,vj,vk,vw,ddij,ddjk,ddik)
  implicit none
  real,intent(in) :: rw
  real,intent(out) :: ddij,ddjk,ddik
  real,dimension(3),intent(in) :: vi,vj,vk,vw
  real,dimension(3) :: vij,vjk,vki,viw,vjw,vkw,cij,cjk
  real,dimension(3,3) :: vxx,vxw,uxx,uxw
  real,dimension(3) :: dxx,ciw,ckw,vv,ua,ub,der
  real :: adist,bdist,cdist,z
  integer :: ijk,i,k,j
  !!
  vij = vj - vi
  vjk = vk - vj
  vki = vi - vk
  vxx(1:3,1) = vij; vxx(1:3,2)=vjk; vxx(1:3,3)=vki
  vxw(1:3,1) = vw - vi
  vxw(1:3,2) = vw - vj
  vxw(1:3,3) = vw - vk
  do ijk=1,3
    dxx(ijk) = sqrt(dotprod(vxx(1:3,ijk),vxx(1:3,ijk)))
    uxx(1:3,ijk) = vxx(1:3,ijk)/dxx(ijk)
    uxw(1:3,ijk) = vxw(1:3,ijk)/sqrt(dotprod(vxw(1:3,ijk),vxw(1:3,ijk)))
  enddo

  do ijk=1,3
     i = ijk; j=mod(ijk,3)+1;  k=mod(ijk+1,3)+1
     !! get triangle distances
     call cros(uxx(1:3,i),uxw(1:3,i),ciw)
     ciw = ciw/sqrt(dotprod(ciw,ciw))
     call cros(uxw(1:3,i),ciw,ua)
     ua = ua/sqrt(dotprod(ua,ua))
     adist = dotprod(vxx(1:3,i),ua)
     call cros(uxx(1:3,k),uxw(1:3,i),ckw)
     ckw = ckw/sqrt(dotprod(ckw,ckw))
     call cros(ckw,uxw(1:3,i),ub)
     ub = ub/sqrt(dotprod(ub,ub))
     bdist = dotprod(-vxx(1:3,k),ub)
     vv = dxx(i)*ua - dxx(k)*ub
     cdist = sqrt(dotprod(vv,vv))
     !! calculate the cosine rule
     z = dotprod(ciw,ckw)  !! cos(kij angle)
     ! z = (cdist**2 - bdist**2 - adist**2)/(-2*adist*bdist)
     z = 1 - z*z
     if (z == 0.) then
       write(0,'(i4,11f12.4)') ijk,z,ciw,uxx(1:3,i),vxx(1:3,i),dxx(i)
       stop 'BUG in strianglederiv().  z == 0'
     endif
     if (z <= 0) then
       write(0,'(i4,11f12.4)') ijk,z,ciw,uxx(1:3,i),vxx(1:3,i),dxx(i)
        stop 'BUG in strianglederiv().  z <= 0'
     endif
     der(i) = (cdist/(adist*bdist))*(1/sqrt(z))
  enddo
  der = rw*rw*der
  ddjk = der(1)
  ddik = der(2)
  ddij = der(3)
end subroutine strianglederiv
end module vectormath
