        module slove_eigen

            implicit none

            contains
c---------------------------------------------------------------------------------------------------
            subroutine zg_eigen(Num,Matrix,Nmatrix, eigenval, eigenvec)
c           This subroutine solves the generalized eigenvalue problem using lapack
c           The ouput eigenvalues and eigenvectors are sorted according to real value of eigenvalues
c           The generalized eigenvalue problem takes the form:
c           AX = lambda*BX
c           INPUT:
c               <Num>: size of the Matrix
c               <Matrix>: A matrix
c               <Nmatrix>: B matrix
c           OUTPUT:
c               <eigenval>: the sorted eigenvalue
c               <eigenvec>: the sorted eigenvector
c---------------------------------------------------------------------------------------------------
                use, intrinsic :: iso_c_binding
                implicit none
                integer, intent(in) :: Num  ! The size of the matrix
                complex(kind=c_double_complex), dimension(Num,Num), intent(inout) :: Matrix, Nmatrix  ! The matrices to be diagonalized
                complex(kind=c_double_complex), dimension(Num), intent(out) :: eigenval  ! The eigenvalues
                complex(kind=c_double_complex), dimension(Num,Num), intent(out) :: eigenvec  ! The eigenvectors
                integer :: info, i
                integer :: lwork
                complex(kind=c_double_complex), dimension(Num) :: alpha,beta  ! Intermediate variables for LAPACK routine
                complex(kind=c_double_complex), dimension(Num,Num) ::vl,vr  ! Left and right eigenvectors
                real(kind=c_double) :: rcond
                integer, dimension(Num) :: ipiv
                character(len=1) :: jobvl, jobvr  ! Options for LAPACK routine
                
                integer, dimension(Num) :: iwork
                complex*16,dimension(1:2*Num) :: work
                real*8,dimension(1:4*Num) :: rwork ! Real workspace for LAPACK routine

                integer :: ie
            
                jobvl = 'N'  ! Do not compute left eigenvectors
                jobvr = 'V'  ! Compute right eigenvectors

                lwork = 2*num
            

                call zggev(jobvl,jobvr,Num,Matrix,Num,Nmatrix,Num,
     &                      alpha,beta,vl,Num,vr,Num,work,lwork,rwork,info)
                if (info /= 0) stop 'Error in zggev'  ! Check for errors in zggev
            
                eigenval = alpha / beta  ! Compute eigenvalues
                eigenvec = vr  ! Store right eigenvectors
            
c                call zlacpy('Full', Num, Num, vr, Num, eigenvec, Num)  ! Copy right eigenvectors to output array
                call sort_eigenpairs(Num,eigenval,eigenvec)



            end subroutine zg_eigen
            
c----------------------------------------------------------------------------------
            subroutine sort_eigenpairs(n, eigenvalues, eigenvectors)
c           This subroutien sorts the eigenval and eigenvectors according to the 
c           real value of the eigenvalues
c----------------------------------------------------------------------------------
                implicit none
                integer, intent(in) :: n
                complex*16, dimension(n), intent(inout) :: eigenvalues
                complex*16, dimension(n,n), intent(inout) :: eigenvectors
                integer :: i, j
                complex*16 :: temp_val
                complex*16, dimension(n) :: temp_vec
            
                do i = 1, n
                    do j = 1, n-i
                        if (real(eigenvalues(j)) > real(eigenvalues(j+1))) then
                            ! Swap eigenvalues
                            temp_val = eigenvalues(j)
                            eigenvalues(j) = eigenvalues(j+1)
                            eigenvalues(j+1) = temp_val
            
                            ! Swap corresponding eigenvectors
                            temp_vec = eigenvectors(:,j)
                            eigenvectors(:,j) = eigenvectors(:,j+1)
                            eigenvectors(:,j+1) = temp_vec
                            ! The second dimension of the array represents which eigenvector it is. 
                        end if
                    end do
                end do
            end subroutine sort_eigenpairs

            subroutine z_lineq(N, A, B)
                ! solve AX = B linear equation
                implicit none
                integer, intent(in) :: N
                complex*16, dimension(1:N,1:N), intent(inout) :: A
                complex*16, dimension(1:N), intent(inout) :: B
                integer, dimension(1:N) :: ipiv
                integer :: info
            
                ! Call LAPACK's zgesv routine
                call zgesv(N, 1, A, N, ipiv, B, N, info)
            
                if (info /= 0) then
                    print*, 'Error: The ', info, '-th parameter had an illegal value.'
                end if
            end subroutine z_lineq
            
            

        end module