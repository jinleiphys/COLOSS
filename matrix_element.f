        module matrix_element
            use system
            use mesh
            use precision
            use constants
            use pot_class
            use generate_laguerre
            use rot_pot
             
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
            subroutine generate_V0(V,para)
c           This subroutine is used to generate the 
c           potential matrix element with parameters <para>
c           INPUT:
c               <para>: the parameter for optical potential
c           OUTPUT:
c               <V>: the potential matrix
c--------------------------------------------------------------
                complex*16, dimension(1:nr,1:nr) :: V
                complex*16, dimension(1:nr ,1:nr) :: coul_matrix,nuc_matrix
                type(pot_para) :: para

                integer :: ir
                real*8 :: a13
                complex*16 :: rr
                real*8 :: rrc,z12

                a13 = masst**(1./3.)
                rrc = a13*para%rc
                z12 = zt*zp
                
                !calculate the array of nuclear potential on laggurre mesh
                call rot_V_nuc(para)
                !calculate the array of coulomb potential on laggurre mesh
                call rot_V_coul(z12,rrc)

                if(matgauss) then
                    ! use gauss quadtrature to evaluate the matrix element with more grids
                    call coul_mat_gauss(z12,rrc,coul_matrix)
                    call nuc_mat_gauss(para,nuc_matrix)
                else
                    ! use the approximation for lagrange functions to evaluate the marix element
                    do ir = 1, nr
                        coul_matrix(ir,ir) = V_coul(ir)
                        nuc_matrix(ir,ir) = V_nuc(ir)
                    end do
                endif
 
                !add up the matrices
                V = 0.d0
                V = nuc_matrix + coul_matrix

            end subroutine
c--------------------------------------------------------------
            subroutine cal_Hmat(l,p,Hmat)
c           This subroutine is used to generate the 
c           Hamiltonian matrix element for partial wave<l>
c           with potential parameters <p>
c           INPUT:
c               <l>: partial wave for the Hamiltonian
c               <p>: the parameter for optical potential
c           OUTPUT:
c               <H>: the Hamiltonian matrix
c--------------------------------------------------------------
                integer :: l 
                type(pot_para) :: p
                complex*16, dimension(1:nr,1:nr) :: Hmat


                complex*16,dimension(1:nr,1:nr) :: Tmat
                complex*16,dimension(1:nr,1:nr) :: Vmat

                Tmat = 0.d0
                Vmat = 0.d0

                call generate_T0(Tmat,l,mu)

                call generate_V0(Vmat,p)

                Hmat = Tmat + Vmat

            end subroutine
c--------------------------------------------------------------
            subroutine cal_b(l,B_vec)
c           This subroutine is used to generate the 
c           inhomogenous term of the linear equation
c           INPUT:
c               <l>: partial wave number
c           OUTPUT:
c               <B_vec>: the inhomogenous term array 
c--------------------------------------------------------------
                integer :: l
                complex*16, dimension(1:nr) :: B_vec

                complex*16 :: vmod

                integer :: ir, i_cor
                

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