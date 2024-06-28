        module npcc

        use constants
        use precision
        use mesh
        use matrix_element
        use slove_eigen

        implicit NONE



            contains 
            subroutine slove_np_cc(ecm)

                real*8 :: ecm

                integer :: i,j

                real*8 :: npmu

                complex*16, dimension(1:nr*2,1:nr*2) :: Hcc
                complex*16, dimension(1:nr*2,1:nr*2) :: Ecc
                complex*16, dimension(1:nr*2,1:nr*2) :: Acc
                complex*16, dimension(1:2*nr,1:2) :: bvec,xvec

                complex*16, dimension(1:nr,1:nr) :: T1,T2
                complex*16, dimension(1:nr,1:nr) :: V11,V22,V12

                complex*16, dimension(1:nr,1:2,1:2) :: cpotr, cpoto


                complex*16 :: rr,xx,zz
                complex*16 :: vcen,vtens,vls
                complex*16 :: vthis 

                complex*16,dimension(1:2,1:2) :: f1,f2,ff,ss,phs

                j=1!total angular momentum


                ! generate the potential array diagonal and nondiagonal
                do i=1,nr
                    ! Reid neutron-proton potential (T=1, soft core)
                    rr=mesh_rr(i)*eitheta!dont forget to rotate
                    xx=rr*0.7d0
                    zz=exp(-xx)
                    vcen=(-10.463d0*zz+105.468d0*zz**2-3187.8d0*zz**4+9924.3d0*zz**6)/xx
                    vtens=-10.463d0*((1+3/xx+3/xx**2)*zz-(12/xx+3/xx**2)*zz**4)/xx+
     &                      351.77d0*zz**4/xx-1673.5d0*zz**6/xx
                    vls=708.91d0*zz**4/xx-2713.1d0*zz**6/xx
                    cpotr(i,1,1)=vcen-2*(j-1)*vtens/(2*j+1)+(j-1)*vls
                    cpotr(i,1,2)=6*vtens*sqrt(j*(j+1.0d0))/(2*j+1)
                    cpotr(i,2,1)=cpotr(i,1,2)
                    cpotr(i,2,2)=vcen-2*(j+2)*vtens/(2*j+1)-(j+2)*vls

                    rr=mesh_rr(i)!no rotation
                    xx=rr*0.7d0
                    zz=exp(-xx)
                    vcen=(-10.463d0*zz+105.468d0*zz**2-3187.8d0*zz**4+9924.3d0*zz**6)/xx
                    vtens=-10.463d0*((1+3/xx+3/xx**2)*zz-(12/xx+3/xx**2)*zz**4)/xx+
     &                      351.77d0*zz**4/xx-1673.5d0*zz**6/xx
                    vls=708.91d0*zz**4/xx-2713.1d0*zz**6/xx
                    cpoto(i,1,1)=vcen-2*(j-1)*vtens/(2*j+1)+(j-1)*vls
                    cpoto(i,1,2)=6*vtens*sqrt(j*(j+1.0d0))/(2*j+1)
                    cpoto(i,2,1)=cpoto(i,1,2)
                    cpoto(i,2,2)=vcen-2*(j+2)*vtens/(2*j+1)-(j+2)*vls
                end do

                npmu = 0.5d0*amu
                call generate_T0(T1,0,npmu)
                call generate_T0(T2,2,npmu)
                ! for benchmark with pierre
                T1=T1/hbarc/hbarc*2d0*mu*20.736d0*2d0
                T2=T2/hbarc/hbarc*2d0*mu*20.736d0*2d0

                do i = 1,nr
                    V11(i,i) = cpotr(i,1,1)
                    V22(i,i) = cpotr(i,2,2)
                    V12(i,i) = cpotr(i,1,2)
                enddo

                Hcc(1:nr,1:nr) = T1 + V11
                Hcc(nr+1:2*nr,nr+1:2*nr) = T2 + V22
                Hcc(1:nr,nr+1:2*nr)= V12
                Hcc(nr+1:2*nr,1:nr)= V12

                call cal_N0()

                Ecc(1:nr,1:nr) = ecm*Nmat
                Ecc(nr+1:2*nr,nr+1:2*nr) = ecm*Nmat


                do i=1,nr
                    vthis = cpotr(i,1,1)
                    bvec(i,1) = mesh_rw(i)**(1d0/2d0)*vthis *fc_rotated(0,i)
                end do

                do i=1,nr
                    vthis = cpotr(i,1,2)
                    bvec(i+nr,1) = mesh_rw(i)**(1d0/2d0)*vthis*fc_rotated(0,i)
                end do

                do i=1,nr
                    vthis = cpotr(i,2,1)
                    bvec(i,2) = mesh_rw(i)**(1d0/2d0)*vthis *fc_rotated(2,i)
                end do

                do i=1,nr
                    vthis = cpotr(i,2,2)
                    bvec(i+nr,2) = mesh_rw(i)**(1d0/2d0)*vthis*fc_rotated(2,i)
                end do


                bvec = bvec* sqrt(eitheta)

                xvec = bvec

                Acc = Ecc - Hcc
                call complex_linear_eq(2*nr,Acc,xvec(:,1))
                Acc = Ecc - Hcc
                call complex_linear_eq(2*nr,Acc,xvec(:,2))



                f1=0d0
                do i=1,2*nr
                    f1(1,1) = f1(1,1) + xvec(i,1)*bvec(i,1)
                    f1(1,2) = f1(1,2) + xvec(i,1)*bvec(i,2)
                    f1(2,1) = f1(2,1) + xvec(i,2)*bvec(i,1)
                    f1(2,2) = f1(2,2) + xvec(i,2)*bvec(i,2)
                end do
                f1 = -f1/ecm

                f2=0d0
                do i=1,nr
                    vthis = cpoto(i,1,1)
                    f2(1,1) = f2(1,1) + mesh_rw(i) * vthis * fc(0,i)**2
                    vthis = cpoto(i,1,2)
                    f2(1,2) = f2(1,2) + mesh_rw(i) * vthis * fc(0,i)*fc(2,i)
                    vthis = cpoto(i,2,1)
                    f2(2,1) = f2(2,1) + mesh_rw(i) * vthis * fc(0,i)*fc(2,i)
                    vthis = cpoto(i,2,2)
                    f2(2,2) = f2(2,2) + mesh_rw(i) * vthis * fc(2,i)**2
                end do
                f2 = -f2/ecm

                ff=f1+f2

                ss=1.d0 + 2.d0*iu*k*ff

                ss(1,2) = ss(1,2) - 1d0
                ss(2,1) = ss(2,1) - 1d0

                phs = atan2( aimag(ss),real(ss) )/2


                !write(*,*) "fborn is", f2
                !write(*,*) "fsc is", f1
                write(*,*) "Smat is", ss
                !write(*,*) "Abs of Smat is", abs(ss)
100             FORMAT('phase shift (rad.)=',2es12.4,' eta_12=',es12.4)    
                WRITE(*,100) real(phs(1,1)),real(phs(2,2)),abs(ss(1,2))



                stop

            end subroutine

        end module