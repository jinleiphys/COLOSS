        module bound

            use constants
            use precision
            use system
            use mesh
            use matrix_element
            use slove_eigen
            use pot_class   
            use rot_pot
            implicit none

            complex*16, dimension(:), allocatable :: eigen_val 
            complex*16, dimension(:,:), allocatable :: eigen_vec
            
            contains
c------------------------------------------------------------------------
            subroutine solve_bound(l,para)
c           This subroutine generate the H matrix with 
c           potential parameters <para> for partial wave <l>,
c           diagonalize the Hamiltonian matrix,
c           and get the eigenvalue and eigenvectors.
c------------------------------------------------------------------------

                integer :: l
                type(pot_para) :: para

                complex*16, dimension(1:nr,1:nr) :: Hmat
                integer :: i
                integer :: cor1, cor2
                complex*16 :: r1,r2
                complex*16 :: norm
                real*8 :: a13
                real*8 :: tan_angle
                integer :: i_eigen
                complex*16 :: wf
                
                if(allocated(eigen_val)) deallocate(eigen_val)
                allocate(eigen_val(1:nr))
                if(allocated(eigen_vec)) deallocate(eigen_vec)
                allocate(eigen_vec(1:nr, 1:nr))

                call cal_Hmat(l,para,Hmat)

                call cal_N0()

                call zg_eigen(nr,Hmat,Nmat,eigen_val,eigen_vec)

                !filter bound and resonance state
                do i=1,nr
                    if( real(eigen_val(i)) < 0.d0 ) then
c                       write(*,*) "bound is", real(eigen_val(i)),aimag(eigen_val(i))
                        cycle
                    end if
                    tan_angle = -aimag(eigen_val(i)) / real(eigen_val(i)) 
                    if( tan_angle < 0.9d0*tan(2.d0*theta_rad) ) then
c                       write(*,*) "resonance is", real(eigen_val(i)), aimag(eigen_val(i))
                    end if
                end do

                ! normalize the wave funtion/eigenvector
                do i = 1, nr 
                    norm = 0.d0
                    do cor1 = 1,nr
                        do cor2 = 1, nr
                            r1 = laguerre_rr(cor1)
                            r2 = laguerre_rr(cor2)
                            norm = norm + eigen_vec(cor1,i)*eigen_vec(cor2,i)*
     &                     ( kronecker_delta(cor1, cor2) + (-1.d0)**(cor1-cor2) /sqrt( r1*r2 ) )
                            !note the the regularized laguerre are not orthogonal
                        end do
                    end do
c                    write(*,*) i," th eigenvector's norm is ", norm
                    eigen_vec(:,i) = eigen_vec(:,i)/sqrt(norm)
                end do


            end subroutine


            function kronecker_delta(i, j) result(delta)
                implicit none
                integer, intent(in) :: i, j
                integer :: delta
            
                if (i == j) then
                    delta = 1
                else
                    delta = 0
                end if
            end function kronecker_delta
            


        end module
