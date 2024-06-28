      module spharm
      contains

c  Subroutine from fresco to calculate the p_{lm}(\theta)
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









      end module
