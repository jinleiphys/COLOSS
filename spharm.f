      module spharm
      contains


      subroutine SPYLM(L,M,cth,Y)
       implicit none
       real*8 :: theta
       integer :: L,M
       real*8 :: cth,Y
       real*8,dimension(0:l,0:l) :: PL
       call PLM(cth,l,l,l+1,PL)
       Y=YLMC(L,M)*PL(l,abs(M)) !  *exp(-imag*abs(ma)*phi)

!       write(*,*) "YLMC=",YLMC(L,M)
      end subroutine



c-----------------------calculate the norm(l,m)-------------------------
!
!  CALCULATE THE COEFFICIENT OF P(L,M)*E(I*M*PHI) IN Y(L,M)
!  AMoro: Uses |M|
!
      FUNCTION YLMC(L,M)
      use cleb_coefficient
c	  use factorials
!      use constants,only:pi
      use precision
      IMPLICIT REAL*8(A-H,O-Z)
      integer,intent(IN):: L,M
      PHASE(I) = (-1)**I
c      pi=acos(-1d0)
      LF1 = L + 1
      MA = ABS(M)
      R =  FLOG(LF1-MA)-FLOG(LF1+MA)
      R = SQRT((2*L+1)/(4*PI)*EXP(R))* PHASE(M)
      IF(M.LT.0) R = R * PHASE(MA)
      YLMC = R
      RETURN
      END function



c-------------------------end calculate norm(l,m)-----------------------

c-----------------------calculate the p_{lm}(\theta)--------------------
c from fresco
c  x: cos(\theta) -1<=x<=1
c  N:   lmax
c  M:   m_max -> lmax
c  NA:  lmax+1
      SUBROUTINE PLM(X,N,M,NA,PL)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 PL(NA,M+1),X
      if(X>1) X=1.0d0
      if(X<-1) X=-1.0d0
      if(n==0) then
        PL=1.0d0
        return
      end if
      N1 = N+1
      M1 = M+1
      DO 10 J=1,M1
      DO 10 I=1,N1
10    PL(I,J)=0.
      PL(1,1) = 1.
      PL(2,1) = X
      SX = SQRT(1.-X*X)
      PL(2,2) = SX
      FACT=1.
      PMM=1.
      DO 15 J=2,min(M1,N1)
        mm = J-1
        PMM = PMM*FACT*SX
        FACT=FACT+2.
        PL(J,J) = PMM
        if(J+1.le.N1) PL(J+1,J) = X*(2*mm+1.) * PL(J,J)
15      CONTINUE
      DO 20 J=1,M1
       mm = J-1
      DO 20 I=J+2,N1
       ll = I-1
      PL(I,J)=((2.*ll-1.)*X*PL(I-1,J) - (ll+mm-1.)*PL(I-2,J))/(ll-mm)
20    CONTINUE
      RETURN
      END subroutine


c  AMoro
c  Factor to convert from P(l,m,x) = plmaux*P(l,-m,x)                    !!! without
c      function plmaux(l,m)
c      use factorials
c      IMPLICIT REAL*8(A-H,O-Z)
c      integer l,m
c      real*8 plmaux,r
c
c      if (m.gt.l) then
c         write(*,*)' plmaux:m>l!'; stop
c      endif
c      if (m.ge.0) then
c        plmaux=1
c      else
c        R =  DLFAC(L+M)-DLFAC(L-M)
c        plmaux=(-1)**m*exp(r)
c      endif
c      return
c      end function



c---------------------end calculate the p_{lm}(\theta)------------------


      end module
