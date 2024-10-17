        module matrix_element
            use system
            use mesh
            use precision
            use constants
            use pot_class
            use generate_laguerre
            use rot_pot
            use channels
             
            implicit none

            complex*16, allocatable,dimension(:,:) :: Nmat

            contains
c-----------------------------------------------------------------
            subroutine generate_T0(Tm,l,mass)
c           This subroutine is used to generate the 
c           kinetic energy and centrifugal barrier matrix element
c           for partial wave <l>
c           INPUT:
c               <l>: partail wave of the centrifugal barrier
c           OUTPUT:
c               <Tm>: kinetic energy and centrifugal barrier matrix
c-----------------------------------------------------------------

                complex*16,dimension(1:nr,1:nr) :: Tm
                integer :: l
                real*8 :: mass

                complex*16 :: aux1, aux2
                integer :: i, j
                real*8 :: ri, rj
                complex*16 :: coeff

                aux1 = 0.d0; aux2 = 0.d0

                Tm = 0.d0

                do i = 1, nr

                    ri = laguerre_rr(i)
                    aux1 = ri**2 - 2.d0*( 2.d0*nr + alpha + 1.d0 )*ri + alpha**2 - 4.d0
                    aux2 = 12.d0*ri**2
                    Tm(i,i) = -aux1/aux2 - 1.d0/4.d0/ri + l*(l + 1.d0)/ri/ri ! dont forget centrifugal barrier

                    do j = i+1, nr

                        rj = laguerre_rr(j)
                        aux1 =  ri + rj
                        aux2 = sqrt(ri*rj)*( ri - rj )**2
                        Tm(i,j) = (-1.d0)**(i-j)*( aux1/aux2 -1.d0/4.d0/sqrt(ri*rj))

                        Tm(j,i) = Tm(i,j)

                    end do
                end do

                !mutiply with the coeffecient
                coeff = hbarc**2/2.d0/mass/eitheta**2/scaling_factor**2
                Tm = Tm * coeff

            end subroutine

c--------------------------------------------------------------
            subroutine generate_V0(V,para,ich)
c           This subroutine is used to generate the 
c           potential matrix element with parameters <para>
c           INPUT:
c               <para>: the parameter for optical potential
c           OUTPUT:
c               <V>: the potential matrix
c--------------------------------------------------------------
                complex*16, dimension(1:nr,1:nr) :: V
                integer :: ich
                complex*16, dimension(1:nr ,1:nr) :: coul_matrix,nuc_matrix
                type(pot_para) :: para

                integer :: ir,jr
                real*8 :: a13
                complex*16 :: rr
                real*8 :: rrc,z12

                a13 = masst**(1./3.)
                rrc = a13*para%rc
                z12 = zt*zp
                
                !calculate Coulomb potentials and matrix elements
                call rot_V_coul(z12,rrc)!calculate the array of coulomb potential on laggurre mesh
                if(matgauss) then
                    call coul_mat_gauss(z12,rrc,coul_matrix)!use gauss quad. to evaluate the coulomb matrix elements
                else
                    do ir = 1, nr
                        coul_matrix(ir,ir) = V_coul(ir)!use the approx. for lagrange func. to evaluate the coulomb marix elements
                    end do
                endif !end if for matgauss

                !calculate nuclear potentials and matrix elements
                if (nonlocal) then
                    call rot_V_nl(para,ich)
                    if(matgauss) then
                        call nuc_mat_nl(para,ich,nuc_matrix)
                    else
                        do ir = 1, nr
                        do jr = ir, nr
                                nuc_matrix(ir,jr) = eitheta*sqrt(mesh_rw(ir) * mesh_rw(jr))*V_nl(ir,jr)
                                nuc_matrix(jr,ir) = nuc_matrix(ir,jr)
                        end do !end for jr
                        end do ! end for ir
                    end if!end if for matgauss


                else
                    call rot_V_nuc(para,ich)
                    if(matgauss) then
                        call nuc_mat_gauss(para,ich,nuc_matrix)
                    else
                        do ir = 1, nr
                            nuc_matrix(ir,ir) = V_nuc(ir)
                        end do
                    end if!end if for matgauss
                end if!end if for nonlocal
 
                !add up the matrices
                V = 0.d0
                V = nuc_matrix + coul_matrix

            end subroutine
c--------------------------------------------------------------
            subroutine cal_Hmat(ich,p,Hmat)
c           This subroutine is used to generate the 
c           Hamiltonian matrix element for partial wave<l>
c           with potential parameters <p>
c           INPUT:
c               <l>: partial wave for the Hamiltonian
c               <p>: the parameter for optical potential
c           OUTPUT:
c               <H>: the Hamiltonian matrix
c--------------------------------------------------------------
                integer :: ich 
                type(pot_para) :: p
                complex*16, dimension(1:nr,1:nr) :: Hmat


                integer :: l
                complex*16,dimension(1:nr,1:nr) :: Tmat
                complex*16,dimension(1:nr,1:nr) :: Vmat

                Tmat = 0.d0; Vmat = 0.d0
                
                l = channel_index%L(ich)
                call generate_T0(Tmat,l,mu)

                call generate_V0(Vmat,p,ich)

                Hmat = Tmat + Vmat

            end subroutine
c--------------------------------------------------------------
            subroutine cal_b(ich,B_vec)
c           This subroutine is used to generate the 
c           inhomogenous term of the linear equation
c           INPUT:
c               <l>: partial wave number
c           OUTPUT:
c               <B_vec>: the inhomogenous term array 
c--------------------------------------------------------------
                integer :: ich
                complex*16, dimension(1:nr) :: B_vec

                complex*16 :: vmod

                integer :: ir, i_cor
                integer :: l

                l = channel_index%L(ich)
                
                if(bgauss) then!use gauss quad to evaluate the inhomo term

                    if(backrot) then

                        do ir = 1,nr 
                            do i_cor = 1, numgauss
                                vmod = Vnuc_gauss(i_cor) + Vcoul_gauss(i_cor) -e2*zp*zt/gauss_rr(i_cor)!no eitheta
                                B_vec(ir) = B_vec(ir) + lag_func_br(i_cor,ir)*
     &                                      gauss_rw(i_cor) *vmod *fc_gauss(l,i_cor)
                            end do
                        end do
                        B_vec = B_vec * exp(iu*cph(l)) / sqrt(eitheta)

                    else

                        do ir = 1,nr 
                            do i_cor = 1, numgauss
                                vmod = Vnuc_gauss(i_cor) + Vcoul_gauss(i_cor) -e2*zp*zt/gauss_rr(i_cor)/eitheta
                                B_vec(ir) = B_vec(ir) + lag_func(i_cor,ir)*
     &                                      gauss_rw(i_cor) *vmod *fc_rot_gauss(l,i_cor)
                            end do
                        end do
                        B_vec = B_vec * exp(iu*cph(l)) * sqrt(eitheta)

                    endif

                else !use lagrange condition to evaluate the inhomo term

                    do ir = 1,nr 
                        vmod = V_nuc(ir) + V_coul(ir) -e2*zp*zt/mesh_rr(ir)/eitheta
                        B_vec(ir) = mesh_rw(ir)**(1d0/2d0) *vmod *fc_rotated(l,ir)
                    end do
                    B_vec = B_vec * exp(iu*cph(l)) * sqrt(eitheta)

                endif

            end subroutine

            subroutine cal_b_nl(ich,B_vec)
                integer :: ich
                complex*16, dimension(1:nr) :: B_vec

                complex*16, dimension(1:nr) :: B_n,B_CS
                complex*16 :: Btemp
                integer :: ii
                integer :: ir,jr
                integer :: l

                complex*16 :: VCS,rr

                l = channel_index%L(ich)

                !==========================================================================================
                if(bgauss) then
                    if(backrot) then
                        do ii = 1, nr
                            do ir = 1,numgauss
                                VCS = Vcoul_gauss(ir) -e2*zp*zt/gauss_rr(ir)
                                B_CS(ii) = B_CS(ii) + gauss_rw(ir)*lag_func_br(ir,ii)*VCS*fc_gauss(l,ir)/eitheta
                            end do
                        end do
                    else
                        do ii = 1, nr
                            do ir = 1,numgauss
                                VCS = Vcoul_gauss(ir) -e2*zp*zt/gauss_rr(ir)/eitheta
                                B_CS(ii) = B_CS(ii) + gauss_rw(ir)*lag_func(ir,ii)*VCS*fc_rot_gauss(l,ir)
                            end do
                        end do
                    end if
                else 
                    write(*,*) "Non local potentials mush use gauss quadrature!"
                    stop
                endif!end if for bgauss
                !==========================================================================================


                !==========================================================================================
                do ii = 1, nr 
                    if(bgauss) then
                        if(backrot) then
                            Btemp = 0.d0
                            do ir = 1,numgauss
                            do jr = 1,numgauss
                                Btemp = Btemp + gauss_rw(ir)*gauss_rw(jr)*lag_func_br(jr,ii)*Vnl_gauss(jr,ir)*fc_gauss(l,ir)/eitheta
                            end do
                            end do
                            B_n(ii) = Btemp
                        else
                            Btemp = 0.d0
                            do ir = 1,numgauss
                            do jr = 1,numgauss
                                Btemp = Btemp + gauss_rw(ir)*gauss_rw(jr)*lag_func(ir,ii)*Vnl_gauss(ir,jr)*fc_rot_gauss(l,jr)* eitheta
                            end do
                            end do
                            B_n(ii) = Btemp
                        end if
                    else
                        write(*,*) "Non local potentials mush use gauss quadrature!"
                        stop
                    end if 
                end do!end do for different ii
                !==========================================================================================

                B_vec = B_n + B_CS
                B_vec = B_vec * exp(iu*cph(l)) * sqrt(eitheta)

            end subroutine

c--------------------------------------------------------------
            subroutine cal_N0()
c           This subroutines generates the inner product matrix
c           for x regularized lagrange-laguerre functions,
c           which is used in the generalized eigenvalue problem 
c--------------------------------------------------------------
                 
                integer :: i,j
                real*8 :: ri,rj
                
                if(allocated(Nmat)) deallocate(Nmat)
                allocate(Nmat(1:nr, 1:nr))

                Nmat = 0.d0

                do i = 1, nr

                    ri = laguerre_rr(i)
                    Nmat(i,i) = 1.d0 + 1.d0/ri

                    do j = i + 1, nr

                        rj = laguerre_rr(j)
                        Nmat(i,j) = (-1.d0)**(i-j)/sqrt(ri*rj)

                        Nmat(j,i) = Nmat(i,j)

                    end do
                end do
            end subroutine

        end module matrix_element