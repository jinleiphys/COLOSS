      module cleb_coefficient
      use precision
        real*8:: dlfac(0:10000)
        real*8:: dl2fac(0:10000)
        real*8:: FACT(0:10000),FFAKINV(0:10000),WFAK(0:10000)
        real*8 :: WFAKINV(0:10000)
      contains
c-----------------------------------------------------------------------
c   call subroutine factorialgen before using this mod
c-----------------------------------------------------------------------
      subroutine factorialgen(n)           !!!!call factorialgen(2*lmax)
! n>0
! ln(i!) and ln((2i+1)!)
      implicit real*8 (a-h,o-z)
      dlfac(0)=0.
      fact(0)=1.
      FFAKINV(0)=1.0d0
      WFAK(0)=1.0d0
      WFAKINV(0)=1.0d0


      do 1 i=1,n
      a=i
1      dlfac(i)=dlfac(i-1)+log(a)
C      fact(i)=I*fact(I-1)
C      FFAKINV(I)=1.0/fact(I)
C      WFAK(I)=SQRT(fact(I))
C1     WFAKINV(I)=1.0/WFAK(I)

      dl2fac(0)=0.d0
      do 2 i=1,n
      a=2.d0*i+1.d0
2     dl2fac(i)=dl2fac(i-1)+log(a)
      continue
      end subroutine
c---------------------------------------------
c******
      function flog(i)
        implicit none
        integer::i
        real*8:: flog
!        real*8 :: logfac
!        flog=logfac(i-1)
        if(i==0) then
           flog=0.0d0
        else
           flog=dlfac(i-1)
        end if
      end function flog

c *** ---------------------------------------------
c Factorial LOG
c *** --------------------------------------------
      real*8 function logfac(n) ! FL(N)
      implicit real*8(a-h,o-z),integer*4(i-n)
       fl=0
       if(n>1) then
       FN = 1.
       DO 10 I = 2,N
       FN = FN + 1.
   10  FL = FL +  LOG(FN)
      endif
      logfac=fl
      END FUNCTION





cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Clebsch-Gordan coefficient  <l',m',l,m|JM>
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function cleb(ia,id,ib,ie,ic,if)
      implicit real*8(a-h,o-z)
C      real*8:: ria,rid,rib,rie,ric,rif
!      common/clebma/faclog(500)
!      COMMON /PRAHA/ FLOG(100), GM(100), DG(25)
C      ia=2.0d0*(ria+.0001d0)
C      ib=2.0d0*(rib+.0001d0)
C      ic=2.0d0*(ric+.0001d0)
C      id=int(sign(1.0_dpreal,rid)*2.0_dpreal*(abs(rid)+.0001d0))
C      ie=int(sign(1.0_dpreal,rie)*2.0_dpreal*(abs(rie)+.0001d0))
C      if=int(sign(1.0_dpreal,rif)*2.0_dpreal*(abs(rif)+.0001d0))
      wwww=-1.0d0
      cleb=0.0d0
      if(id+ie-if) 7000,105,7000
  105 k1=ia+ib+ic
      if((-1)**k1) 7000,107,107
  107 if(.not.((id.eq.0).and.(ie.eq.0))) go to 110
      k1=k1/2
      if((-1)**k1) 7000,110,110
  110 k1=ia+ib-ic
      k2=ic-iabs(ia-ib)
      k3=min0(k1,k2)
      if(k3) 7000,130,130
  130 if((-1)**(ib+ie)) 7000,7000,140
  140 if((-1)**(ic+if)) 7000,7000,150
  150 if(ia-iabs (id)) 7000,152,152
  152 if(ib-iabs (ie)) 7000,154,154
  154 if(ic-iabs (if)) 7000,160,160
  160 if(ia) 7000,175,165
  165 if(ib) 7000,175,170
  170 if(ic) 7000,180,250
  175 cleb=1.0d0
      go to 7000
  180 fb=float(ib+1)
      cleb=((wwww)**((ia-id)/2))/sqrt(fb)
      go to 7000
  250 fc2=ic+1
      iabcp=(ia+ib+ic)/2+1
      iabc=iabcp-ic
      icab=iabcp-ib
      ibca=iabcp-ia
      iapd=(ia+id)/2+1
      iamd=iapd-id
      ibpe=(ib+ie)/2+1
      ibme=ibpe-ie
      icpf=(ic+if)/2+1
      icmf=icpf-if
      vvv=0.5d0
      sqfclg=vvv*(log(fc2)-flog(iabcp+1)
     1      +flog(iabc)+flog(icab)+flog(ibca)
     2      +flog(iapd)+flog(iamd)+flog(ibpe)
     3      +flog(ibme)+flog(icpf)+flog(icmf))
      nzmic2=(ib-ic-id)/2
      nzmic3=(ia-ic+ie)/2
      nzmi= max0(0,nzmic2,nzmic3)+1
      nzmx= min0(iabc,iamd,ibpe)
      if(nzmx.lt.nzmi) go to 7000
      s1=(wwww)**(nzmi-1)
      do 400 nz=nzmi,nzmx
      nzm1=nz-1
      nzt1=iabc-nzm1
      nzt2=iamd-nzm1
      nzt3=ibpe-nzm1
      nzt4=nz-nzmic2
      nzt5=nz-nzmic3
      termlg=sqfclg-flog(nz)-flog(nzt1)-flog(nzt2)
     1           -flog(nzt3)-flog(nzt4)-flog(nzt5)
      ssterm=s1*exp (termlg)
      cleb=cleb+ssterm
  400 s1=-s1
 7000 return
      end function
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Clebsch-Gordan coefficient from FRESCO <l',m',l,m|JM>
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION CLEB6(A,AL,B,BE,C,M)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8,intent(IN):: A,AL,B,BE,C,M
!      real*8 :: flog
      LOGICAL FAIL3,FRAC
      PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-INT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      GA = -M
C
      IF(AL+BE+GA) 11,10,11
C11    WIG3J = 0.0
11    CLEB6 = 0.0
      RETURN
10    IF(FAIL3(C,A,B)) GO TO 11
      IF(A-ABS(AL)) 11,14,14
14    IF(B-ABS(BE)) 11,15,15
15    IF(C-ABS(GA)) 11,13,13
13    IA = C-B+AL
      IB = C-A-BE
      IF(IA) 20,21,21
21    IF(IB) 20,23,23
23    MIN = 0
      GO TO 24
20    MIN = -IA
      IF(MIN+IB) 25,24,24
25    MIN = -IB
24    IC = A-AL
      ID = B+BE
      IE = A+B-C
      NIN = MIN
      T=PHASE(MIN)
      S = T
30    MIN = MIN+1
      IZ = IC-MIN+1
      IF(IZ) 29,29,26
26    IY = ID-MIN+1
      IF(IY) 29,29,27
27    IX = IE-MIN+1
      IF(IX) 29,29,28
28    TA = real(IX)*real(IY)*real(IZ)
      TB = MIN*real(IA+MIN)*real(IB+MIN)
      T = -T*TA/TB
      S = S+T
      GO TO 30
29    I = B-A+GA
      IF(S.EQ.0.0) GO TO 11
      IH = A+AL
      II = B-BE
      IJ = C+GA
      IK = C-GA
      IL = A+C-B
      IM = B+C-A
      IN = A+B+C+1.0
      XDDD = 0.5*(flog(IH+1)+flog(IC+1)+flog(ID+1)
     $           +flog(II+1)+flog(IJ+1)+flog(IK+1)
     $           +flog(IE+1)+flog(IM+1)+flog(IL+1)-flog(IN+1))
     $      - ( flog(IC-NIN+1)+flog(IA+NIN+1)+flog(ID-NIN+1)+
     $          flog(IB+NIN+1)+flog(NIN+1)+flog(IE-NIN+1))
C     WIG3J = (-1.0)**I *  EXP(XDDD) * S
      CLEB6 = SQRT(2*C+1.)*EXP(XDDD) * S
      RETURN
      END  FUNCTION
c------
CG coefficients from Andreas Nogga
c------
      function CG (I,L,J,M,K,N)
      use precision
      implicit real(dpreal) (a-h,o-z)
C*******************************************************************************
C                    <I/2 L/2 J/2 M/2 K/2 N/2 >
C  CLEBSCH-GORDON
C                     <2l',2m',2l,2m|2J2 M>
C*******************************************************************************
      INTEGER T
      IF (L+M .NE. N .OR. ABS(I-J) .GT. K .OR. K .GT. I+J .OR.
     .    ABS(L) .GT. I .OR. ABS(M) .GT. J .OR. ABS(N) .GT. K) THEN
       CG=0._dpreal
      RETURN
      ENDIF
      I1=(I+L)/2
      I2=(I-L)/2
      I3=(J+M)/2
      I4=(J-M)/2
      I5=(K-N)/2
      I6=(K+N)/2
      I7=(I+J-K)/2
      I8=(J+K-I)/2
      I9=(K+I-J)/2
      I10=(I+J+K+2)/2

      XX=WFAK(I7)*WFAK(I8)*WFAK(I9)*WFAK(I1)*WFAK(I2)*
     X     WFAK(I3)*WFAK(I4)*WFAK(I5)*
     .   WFAK(I6)*WFAKINV(I10)
      J1=(K-J+L)/2
      J2=(K-I-M)/2
      NT=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9)+1
      IT=0
      SUM=0._dpreal
      T=-1
   10   T=T+1
        L1=J1+T
        L2=J2+T
        L3=I7-T
        L4=I2-T
        L5=I3-T
        IF (MIN(L1,L2,L3,L4,L5) .LT. 0) GOTO 10
       SUM=SUM+RM1H(T)*FFAKINV(T)*FFAKINV(L1)*FFAKINV(L2)
     X       *FFAKINV(L3)*FFAKINV(L4)*FFAKINV(L5)
       IT=IT+1
       IF (IT .LT. NT) GOTO 10
      CG=XX*SUM*SQRT(K+1._dpreal)
      END function



      FUNCTION RM1H (I)
       use precision
      implicit real(dpreal) (a-h,o-z)
C*******************************************************************************
C
C  RM1H = (-1) ** I
C
C*******************************************************************************
      RM1H=(-1.0_dpreal)**I    ! 1-2*MOD(ABS(I),2)
      END FUNCTION
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     spherical harmonics Y_lm(\theta,\psi)
c     Y_{lm}(\theta,\psi)=norm(l,m)*p_{lm}(\theta)*exp(i*m*psi)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      subroutine YLM(l,m,theta,psi)
c
c      implicit none
c      integer,intent(in) :: l,m
c      real*8,intent(in) :: theta,psi       ! in Radians
c      real*8 :: norm
c
c      norm=YLMC(l,m)
c
c      end subroutine









cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ------------------------------------------
c     Six-J coefficient {{j1,j2,j3},{j4,j5,j6}}
c     ------------------------------------------
      real*8 function sixj(r1,r2,r3,r4,r5,r6)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-h,o-z)
        real*8:: phase
        phase=(-1)**(r1+r2+r5+r4)
        sixj=phase*rac(r1,r2,r5,r4,r3,r6)
      end function
c-----------------------------------------------------------------------


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function rac(ria,rib,ric,rid,rie,rif)
c     -------------------------------------------------------------------
c     subroutine calculates the racah coefficient w(abcd;ef) defined
c     according to the convention of brink and satchler.
c     the arguments are real and are the actual values of the angular
c     momenta, ( i.e. they can take half integer values )
c     -------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-h,o-z)
!      COMMON /PRAHA/ FLOG(100), GM(100), DG(25)

!      common/clebma/faclog(500)
      dimension lt(6)

!      write(*,*)'In rac:'
!      write(*,*) "2",flog(2),exp(flog(2))
!      write(*,*)"10",flog(10),exp(flog(10))
! 	write(*,*)flog(10),exp(flog(10))

      rac=0.0d0
      ia=2.d0*(ria+.0001d0)
      ib=2.d0*(rib+.0001d0)
      ic=2.d0*(ric+.0001d0)
      id=2.d0*(rid+.0001d0)
      ie=2.d0*(rie+.0001d0)
      if=2.d0*(rif+.0001d0)
      k1=ia+ib-ie
      k2=ie-iabs (ia-ib)
      k3=ic+id-ie
      k4=ie-iabs (ic-id)
      k5=ia+ic-if
      k6=if-iabs (ia-ic)
      k7=ib+id-if
      k8=if-iabs(ib-id)
      k9= min0 (k1,k2,k3,k4,k5,k6,k7,k8)
      if(k9) 7000,20,20
   20 k2=k1-2*(k1/2)
      k4=k3-2*(k3/2)
      k6=k5-2*(k5/2)
      k8=k7-2*(k7/2)
      if(max0(k2,k4,k6,k8)) 7000,25,7000
   25 ltmin=min0(ia,ib,ic,id,ie,if)
      if(ltmin) 7000,30,150
   30 lt(1)=ia
      lt(2)=ib
      lt(3)=ic
      lt(4)=id
      lt(5)=ie
      lt(6)=if
      ltmin=lt(1)
      kmin=1
      do 40 n=2,6
      if(lt(n)-ltmin) 35,40,40
   35 ltmin=lt(n)
      kmin=n
   40 continue
      s1=1.0d0
      f1=ie
      f2=if
      go to (55,55,55,55,45,50),kmin
   45 f1=ia
      f2=ic
      s1=(-1.d0)**(k5/2)
      go to 55
   50 f1=ia
      f2=ib
      s1=(-1.d0)**(k1/2)
   55 rac=s1/dsqrt((f1+1.d0)*(f2+1.d0))
      go to 7000
  150 iabep=(ia+ib+ie)/2+1
      icdep=(ic+id+ie)/2+1
      iacfp=(ia+ic+if)/2+1
      ibdfp=(ib+id+if)/2+1
      iabe=iabep-ie
      ieab=iabep-ib
      ibea=iabep-ia
      icde=icdep-ie
      iecd=icdep-id
      idec=icdep-ic
      iacf=iacfp-if
      ifac=iacfp-ic
      icfa=iacfp-ia
      ibdf=ibdfp-if
      ifbd=ibdfp-id
      idfb=ibdfp-ib
      iabcd1=(ia+ib+ic+id+4)/2
      iefmad=(ie+if-ia-id)/2
      iefmbc=(ie+if-ib-ic)/2
      nzmax=min0(iabe,icde,iacf,ibdf)
      nzmi1=-iefmad
      nzmi2=-iefmbc
      nzmin=max0(0,nzmi1,nzmi2)+1
      if(nzmax.lt.nzmin) go to 7000
      sqlog=flog(iabe)+flog(ieab)+flog(ibea)+flog(icde)+
     &      flog(iecd)+flog(idec)+flog(iacf)+flog(ifac)+
     &      flog(icfa)+flog(ibdf)+flog(ifbd)+flog(idfb)-
     &      flog(iabep+1)-flog(icdep+1)-flog(iacfp+1)-flog(ibdfp+1)
      sqlog=0.5d0*sqlog
      do 200 nz=nzmin,nzmax
      nzm1=nz-1
      k1=iabcd1-nzm1
      k2=iabe-nzm1
      k3=icde-nzm1
      k4=iacf-nzm1
      k5=ibdf-nzm1
      k6=nz
      k7=iefmad+nz
      k8=iefmbc+nz
      sslog=sqlog+flog(k1)-flog(k2)-flog(k3)-flog(k4)
     &           -flog(k5)-flog(k6)-flog(k7)-flog(k8)
      ssterm=((-1.d0)**nzm1)*exp(sslog)
      rac=rac+ssterm
!      write(96,'(8i3,20g12.6)') k1,k2,k3,k4,k5,k6,k7,k8,
!     &                      flog(k1),flog(k2),flog(k3),flog(k4),
!     &                      flog(k5),flog(k6),flog(k7),flog(k8)
  200 continue
 7000 return
      end function
c-----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function u9(ra,rb,rc,rd,re,rf,rg,rh,ri)
c     nine-j symbol. definition as in brink and satchler.
c      | a b c |
c      | d e f |
c      | g h i |
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-h,o-z)
      u9=0.0d0
      k1=idint(2.d0*abs(ra-ri)+0.01d0)
      k2=idint(2.d0*abs(rb-rf)+0.01d0)
      k3=idint(2.d0*abs(rd-rh)+0.01d0)
      minrda=max0(k1,k2,k3)
      k1=idint(2.d0*(ra+ri)+0.01d0)
      k2=idint(2.d0*(rb+rf)+0.01d0)
      k3=idint(2.d0*(rd+rh)+0.01d0)
      maxrda=min0(k1,k2,k3)
      if(minrda-maxrda) 30,30,20
   30 do 50 n1=minrda,maxrda,2
      r1=float(n1)/2.d0
      ramda2=n1
      y9=(ramda2+1.d0)*rac(ra,ri,rd,rh,r1,rg)*rac(rb,rf,rh,rd,r1,re)
      u9=u9+y9*rac(ra,ri,rb,rf,r1,rc)
   50 continue
   20 return
      end function
c-----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc









      SUBROUTINE WHIT(HETA,R,XK,E,LL,F,FD,IE)
C
C     CALCULATES  WHITTAKER  FUNCTION  WL(K,R)  WITH
C     ASYMPTOTIC  FORM  EXP(-(KR + ETA(LOG(2KR)))
C     E  IS  NEGATIVE
C     If IE = 0, allowed to return result e**IE larger than Whittaker,
C                for the IE value returned.
C     If IE > 0, must scale results by that amount.
C
!	use drier !  AMoro
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(LL+1),FD(LL+1) ,T(12),S(7)
!! AMoro: to replace drier module
      REAL*8 FPMAX
!	acc8 = epsilon(acc8);
      fpmax = huge(acc8)**0.8d0
!! ------------------------------

      L = LL+1
C              NOW L = NO. OF VALUES TO FIND
      EE=-1.0
      AK=XK
      ETA=HETA
      LP1=L+1
      RHO=AK*R
	S(:) = 0
      IF(L-50)1,1,2
    1 LM=60
      GO TO 3
    2 LM=L+10
    3 LMP1=LM+1
      IS=7
      PJE=30.0*RHO+1.0
      H=max(INT(PJE),4)
      H=RHO/H
      RHOA=10.0*(ETA+1.0)
      IF(RHOA-RHO)13,13,14
   13 IFEQL=1
      RHOA=RHO
      GO TO 15
   14 IFEQL=0
   15 PJE=RHOA/H+0.5
      RHOA=H*INT(PJE)
      IF(IFEQL)16,16,18
   16 IF(RHOA-RHO-1.5*H)17,18,18
   17 RHOA=RHO+2.0*H
   18 IF(EE)55,55,19
   19 STOP 'WHIT'
   27 A=2.0-10.0/12.0*H*H*EE
      B=1.0/6.0*H*ETA
      C=1.0+1.0/12.0*H*H*EE
      M1=INT(RHOA/H-0.5)
      M2=INT(RHO/H-1.5)
      T(2)=B/FLOAT(M1+1)
      T(3)=B/FLOAT(M1)
      JS=M1
      DO 29 IS=M2,M1
      DO 28 I=1,6
      S(I)=S(I+1)
   28 CONTINUE
      T(1)=T(2)
      T(2)=T(3)
      T(3)=B/FLOAT(JS-1)
      S(7)=((A+10.0*T(2))*S(6)-(C-T(1))*S(5))/(C-T(3))
      JS=JS-1
      IF(ABS(S(7)).LE.FPMAX) GO TO 29
       DO 285 I=2,7
  285   S(I) = S(I) / FPMAX
   29 CONTINUE
      T(1)=S(4)
      T(2)=(1.0/60.0*(S(1)-S(7))+0.15*(S(6)-S(2))+0.75*(S(3)-S(5)))/H
      GO TO 60
   55 C=1.0/RHOA
      A=1.0
      B=1.0-C*ETA
      F(1)=A
      FD(1)=B
      DO 56 M=1,26
      D=0.5*(ETA+FLOAT(M-1))*(ETA+FLOAT(M))*C/FLOAT(M)
      A=-A*D
      B=-B*D-A*C
      F(1)=F(1)+A
      FD(1)=FD(1)+B
   56 CONTINUE
      A=-ETA*LOG(2.0*RHOA)-RHOA
      FPMINL = -LOG(FPMAX)
      if(IE.eq.0.and.A.LT.FPMINL) IE = INT(FPMINL-A)
      A=EXP(A+IE)
      F(1)=A*F(1)
c      FD(1)=A*FD(1)
      FD(1)=A*FD(1) * (-1d0 - 2*ETA/(RHOA))
      IF(IFEQL)57,57,61
   57 S(IS)=F(1)
      IF(IS-7)27,58,27
   58 IS=6
      RHOA=RHOA+H
      GO TO 55
   60 F(1)=T(1)
      FD(1)=T(2)
   61 C=1.0/RHO
      DO 63 M=1,L-1
      A=ETA/FLOAT(M)
      B=A+C*FLOAT(M)
      F(M+1)=(B*F(M)-FD(M))/(A+1.0)
      FD(M+1)=(A-1.0)*F(M)-B*F(M+1)
   63 CONTINUE
      DO 65 M=1,L
      FD(M)=AK*FD(M)
   65 CONTINUE
      RETURN
      END  SUBROUTINE


      end module
