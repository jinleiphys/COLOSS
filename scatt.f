        module scatt
            use mesh
            use system
            use precision
            use constants
            use matrix_element
            use csteed
            use bound
            use rot_pot
            use pot_class
            use coulfunc

            implicit none

            complex*16, dimension(:), allocatable :: scatt_amp_nuc_l


            contains

            subroutine solve_scatt_green(l)
                
                integer :: l
                integer :: ir 
            
                complex*16, dimension(1:nr) :: dn
                complex*16 :: vmod

                complex*16 :: f_sc,f_born,ftot!different amplitudes
                integer :: i_eigen, i_cor!iteration index
                complex*16 :: smat
                
                dn = 0.d0
                do i_eigen = 1,nr
                    do i_cor = 1, nr
                        vmod = V_nuc(i_cor) + V_coul(i_cor) -e2*zp*zt/mesh_rr(i_cor)/eitheta
c                        write(*,*) "vmod is",vmod
                        dn(i_eigen) = dn(i_eigen) + 
     &                  vmod * fc_rotated(l,i_cor) * eigen_vec(i_cor,i_eigen)*sqrt(mesh_rw(i_cor))
                    end do

c                    write(*,*) i_eigen,"is",dn(i_eigen)
                end do

c                write(*,*)"coul is", sin(k*100.d0*eitheta)

                !calculate the f_sc
                f_sc = 0.d0
                do i_eigen = 1, nr
                    f_sc = f_sc + dn(i_eigen)**2/(ecm - eigen_val(i_eigen))
c                    write(*,*) i_eigen,"is",dn(i_eigen)**2/(ecm - eigen_val(i_eigen))
                end do
                f_sc = -eitheta/ecm*f_sc
                
                !calculate the f_born
                f_born = 0.d0
                do ir = 1, nr
                    vmod = V_nuc_origin(ir) + Vcoul_origin(ir) - e2*zp*zt/mesh_rr(ir)
                    f_born = f_born + mesh_rw(ir) * vmod * fc(l,ir)**2
                end do
                f_born = -f_born/ecm

                ftot = f_sc + f_born
                scatt_amp_nuc_l(l) = ftot
                smat = 1.d0 + 2.d0*iu*k*ftot

                write(*, 300) l, real(smat), aimag(smat)
300             FORMAT(I3,' | (',F10.6,', ',F10.6,')  | ')
C               write(*,100) l, real(smat), aimag(smat)
C100             FORMAT(' l=',I4,' S-matrix is (',f10.6,' , ',f10.6, ' )')
                write(60,101) l,real(smat),aimag(smat)
101             FORMAT('l=',I3,2f10.6) 
                write(61,102) l,real(ftot),aimag(ftot)
102             FORMAT('l=',I3,2f14.8)


                
            end subroutine

            subroutine xsec(eta,k,fl,jmax,thetah,thetanmax)
                use coulfunc
                use spharm

                real*8 :: eta
                real*8 :: k
                integer :: jmax
                complex*16, dimension(0:jmax) :: fl
                real*8 :: thetah
                integer :: thetanmax

                real*8, dimension(0:jmax) :: cph
                real*8 :: cph0

                complex*16, dimension(1:thetanmax) :: fc_theta, fn_theta
                complex*16, dimension(1:thetanmax) :: ftot_theta,ftot_rel
                integer :: itheta
                real*8 :: ang_rad

                real*8, dimension(0:jmax,0:jmax) :: legendre_poly
                real*8 :: costheta
                integer :: ll

                if(eta>0) then

                    call coulph(eta,cph,jmax)
                    cph0 = cph(0)

                    !calculate the coulomb scattering amp for different angles
                    fc_theta = 0d0
                    do itheta = 1, thetanmax
                        ang_rad = itheta*thetah/180d0*pi
                        fc_theta(itheta) = -eta/k/2d0/sin(0.5d0*ang_rad)**2d0
     &                                      *exp( 2d0*iu*(cph0 - eta*log(sin(0.5d0*ang_rad))) )
                    end do

                endif

                !calculate the short range scattering amp for different angles
                fn_theta = 0d0
                do itheta = 1, thetanmax
                    ang_rad = itheta*thetah/180d0*pi
                    costheta = cos(ang_rad)
                    call PLM( costheta, jmax, jmax, jmax+1, legendre_poly )
                    do ll = 0, jmax
                        fn_theta(itheta) = fn_theta(itheta) + 
     &                          (2d0*ll+1)*exp(2d0*iu*cph(ll))*fl(ll)*legendre_poly(ll,0)
                    end do
                end do

                if(eta>0) then
                    ftot_theta = fn_theta + fc_theta!add up two amplitudes
                    ftot_rel = ftot_theta/fc_theta

                    do itheta = 1, thetanmax
                        write(67,110) itheta*thetah, abs(ftot_rel(itheta))**2
                    end do

                else
                    
                do itheta = 1, thetanmax
                    write(67,110) itheta*thetah, abs(fn_theta(itheta))**2
                end do

                endif
110             FORMAT(f6.2, es12.4)

            end subroutine


            subroutine solve_scatt(l,para)
                use coulfunc
                use slove_eigen
                type(pot_para) :: para
                integer :: l

                complex*16,dimension(1:nr,1:nr) :: Hmat
                complex*16,dimension(1:nr) :: B_vec,X_vec
                integer :: ir, i_cor
                complex*16 :: vmod
                real*8 :: reac_xsec

                complex*16, dimension(1:nr,1:nr) :: A_mat

                complex*16 :: f_born, f_sc, ftot, smat

                call cal_Hmat(l,para,Hmat)
                call cal_N0()

                B_vec = 0d0
                call cal_b(l,B_vec)

                X_vec = B_vec! copy the inhomo term, and zgesv will return the result into it
                
                A_mat = ecm*Nmat - Hmat

                call z_lineq(nr,A_mat,X_vec)

                !calculate the f_born
                f_born = 0.d0
                do ir = 1, nr
                    vmod = V_nuc_origin(ir) + Vcoul_origin(ir) - e2*zp*zt/mesh_rr(ir)
                    f_born = f_born + mesh_rw(ir) * vmod * fc(l,ir)**2
                end do
                f_born = -f_born/ecm

                !calculate the f_sc
                f_sc = 0d0
                do ir = 1, nr
                    f_sc = f_sc + X_vec(ir)*B_vec(ir)
                end do
                f_sc = -f_sc/ecm/exp(2d0*iu*cph(l))

                ftot = f_born + f_sc
                scatt_amp_nuc_l(l) = ftot
                
                smat = 1.d0 + 2.d0*iu*k*ftot
                reac_xsec = pi/k/k*(2d0*l+1d0)*(1d0 - abs(smat)**2)

                write(*, 300) l, real(smat), aimag(smat), reac_xsec
300             FORMAT(I3,' | (',F10.6,', ',F10.6,')  | ',F10.4)
C               write(*,99) l
99              FORMAT(' l = ',I3,':')
C               write(*,100) real(smat), aimag(smat)
100             FORMAT(' S-matrix is (',f10.6,' , ',f10.6, ' )')
C               write(*,111) reac_xsec
111             FORMAT(' Partial wave reaction cross section is',f10.4,' mb')
                write(60,101) l,real(smat),aimag(smat)
101             FORMAT('l=',I3,2f10.6) 
                write(61,102) l,real(ftot),aimag(ftot)
102             FORMAT('l=',I3,2f14.8)


            end subroutine



        end module