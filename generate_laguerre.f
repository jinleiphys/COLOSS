        module generate_laguerre
        
        use mesh 
        use precision
        use system
        use constants
        use gauss_mesh
        use gauss
        use pot_class
        implicit none

        complex*16, dimension(:,:), allocatable :: lag_func!lagrange function on Gauss-Legendre mesh
        complex*16, dimension(:,:), allocatable :: lag_func_br!back rotated lagrange function on Gauss-Legendre mesh 

        real*8, dimension(:), allocatable :: cph!coulomb phase shift

        complex*16,dimension(:,:),allocatable :: fc_rotated
        complex*16,dimension(:,:),allocatable :: fc_rot_gauss
        complex*16,dimension(:,:), allocatable :: fc
        complex*16,dimension(:,:), allocatable :: fc_gauss
        real*8,dimension(:),allocatable :: nfc,ngc!regular and irregular of Coulumb function F and G
        real*8,dimension(:),allocatable :: nfcp,ngcp!derivatives of F and G

        contains 
c---------------------------------------------------------------
        subroutine init_laguerre_mesh()
c       This subroutine is used to generate the 
c       (scaled) laguerre mesh arrays and weight arrays
c---------------------------------------------------------------

        
                real*8 :: beta = 0.d0

                integer :: ir

                if(allocated(mesh_rr)) deallocate(mesh_rr)
                if(allocated(mesh_rw)) deallocate(mesh_rw)
                if(allocated(laguerre_rr)) deallocate(laguerre_rr)
                if(allocated(laguerre_rw)) deallocate(laguerre_rw)

                allocate(mesh_rr(1:nr))
                allocate(mesh_rw(1:nr))
                allocate(laguerre_rr(1:nr))
                allocate(laguerre_rw(1:nr))

                !first generate the original mesh points and weights
                call cdgqf (nr, 5, alpha, beta, laguerre_rr, laguerre_rw)

                !now transform the orginal mesh and weights into the actual mesh and weights
                write(*,*) '-------------- Lagrange-Laguerre Mesh --------------'
                write(*,99) nr
99              format(' Laguerre Polynomial Order:',I4)
                write(*,100) laguerre_rr(nr), Rmax
100             Format(' Laguerre Mesh Max Value: Scaled from ',f10.4,' to ',f10.4)
                scaling_factor = Rmax/laguerre_rr(nr)
                write(*,101) scaling_factor
101             FORMAT(' Scaling Factor:', f10.4, ' (fm)')
                write(*,*) ''
                

                mesh_rr=laguerre_rr * scaling_factor

                write(10,*) '& The scaled Lagrane-Laguerre mesh and weights are'
                do ir=1, nr
                        !the weight function is: (x-a)^alpha*exp(-(x-a)) 
                        mesh_rw(ir) =  laguerre_rw(ir) * scaling_factor / laguerre_rr(ir)**alpha /exp( -laguerre_rr(ir))  
                        write(10,*)  mesh_rr(ir), mesh_rw(ir) 
                end do 
                
        end subroutine

                

c---------------------------------------------------------------
        subroutine basis_func()
c       This subroutine is used to generate the 
c       lagrange-laguerre basis on gauss mesh
c---------------------------------------------------------------
                      
                integer :: i,j
                complex*16 :: xi,xj
                integer :: i_cor, i_basis
                complex*16 :: rr
                complex*16 :: prod
                real*8 :: rrc
                        
                if(allocated(lag_func)) deallocate(lag_func)
                allocate( lag_func(1:numgauss, 1:nr) ) 
        
                if(allocated(gauss_rr)) deallocate(gauss_rr)
                if(allocated(gauss_rw)) deallocate(gauss_rw)
                allocate( gauss_rr(1:numgauss), gauss_rw(1:numgauss) )
                        
                rrc = input_pot%rc*masst**(1d0/3d0)

                call TRNS(numgauss/2, numgauss/2, numgauss,rrc,20d0,rmaxgauss,gauss_rr, gauss_rw)
c                call gauleg(numgauss, 0.d0, rmaxgauss, gauss_rr, gauss_rw)
        
                do i_basis = 1,nr
                        do i_cor = 1, numgauss
                                xi = laguerre_rr(i_basis)
                                rr = gauss_rr(i_cor)/scaling_factor
                                lag_func(i_cor,i_basis) = mesh_rw(i_basis)**(-1.d0/2.d0)
     &                                  *(rr/xi)**(alpha/2.d0+1.d0) *exp( -(rr-xi)/2.d0 )
                                prod = 1.d0
                                do j = 1, nr
                                        if(j .eq. i_basis) cycle
                                        xj = laguerre_rr(j)
                                        prod = prod * (rr-xj)/(xi-xj)
                                end do
                                lag_func(i_cor,i_basis)=lag_func(i_cor,i_basis)*prod
                        end do
                end do

                if(backrot) then

                if(allocated(lag_func_br)) deallocate(lag_func_br)
                allocate( lag_func_br(1:numgauss, 1:nr) )

                do i_basis = 1,nr
                        do i_cor = 1, numgauss
                                xi = laguerre_rr(i_basis)
                                rr = gauss_rr(i_cor)/scaling_factor/eitheta
                                lag_func_br(i_cor,i_basis) = mesh_rw(i_basis)**(-1.d0/2.d0)
     &                                  *(rr/xi)**(alpha/2.d0+1.d0) *exp( -(rr-xi)/2.d0 )
                                prod = 1.d0
                                do j = 1, nr
                                        if(j .eq. i_basis) cycle
                                        xj = laguerre_rr(j)
                                        prod = prod * (rr-xj)/(xi-xj)
                                end do
                                lag_func_br(i_cor,i_basis)=lag_func_br(i_cor,i_basis)*prod
                        end do
                end do

                endif
                      
        end subroutine

c---------------------------------------------------------------
        subroutine init_coul()
c       This subroutine is used to generate the 
c       complex coulomb wave function on rotated/original
c       coordinate on laguerre/legendre mesh
c---------------------------------------------------------------
                use mesh
                use cwf_cpp
                use coulfunc

                
                integer :: i
                complex*16 :: rho

                !some variable for COULCC
                integer :: NL
                complex*16 :: zlmin
                integer :: mode
                integer :: ifail
                
                complex*16, dimension(0:lmax) :: gc,fcp,gcp,sig
                integer :: il

                complex*16 :: cmplxzero
                real*8 :: t1,t2

                if(allocated(cph)) deallocate(cph)
                allocate(cph(0:lmax))
                cph = 0d0
                call coulph(real(eta),cph,lmax)

                if(allocated(fc_rotated)) deallocate(fc_rotated)
                allocate(fc_rotated(0:lmax,1:nr))
                fc_rotated=0d0

                if(allocated(fc)) deallocate(fc)
                allocate(fc(0:lmax,1:nr))
                fc = 0d0
            
                if(bgauss) then
                    if(allocated(fc_rot_gauss)) deallocate(fc_rot_gauss)
                    allocate(fc_rot_gauss(0:lmax,1:numgauss))
                    fc_rot_gauss=0d0
                end if

                if(backrot) then
                    if(allocated(fc_gauss)) deallocate(fc_gauss)
                    allocate(fc_gauss(0:lmax,1:numgauss))
                    fc_gauss = 0d0
                endif

                mode =4
                zlmin = 0.d0

                write(*,*)"------------------------------------------------"
                write(*,*)'Initializing Complex Coulomb Function'
                write(*,12) real(eta)
12              Format(' Sommerfeld Parameter:', f10.4) 
                


                NL = 1+lmax


                write(*,*)'Generating Rotated Coulomb Wave Function on Laguerre Mesh'
                call cpu_time(t1)
                do i = 1,nr

                    rho = k*mesh_rr(i)*eitheta! rotate the coordinate

                    select case(cwftype)

                    case(1)
                        call COULCC(rho,eta,zlmin,NL,fc_rotated(:,i),gc,fcp,gcp,sig,mode,0,ifail)
                        if (ifail/=0) then
                           write(*,*) 'coul90: ifail=',ifail
                        endif

                    case(2)
                        call compute_coulomb_wave_functions('true', eta,rho,cmplxzero,NL,fc_rotated(:,i),fcp,gc,gcp)
                        
                    end select
                    

                    rho = k*mesh_rr(i)! no rotation always use COULCC
                    call COULCC(rho,eta,zlmin,NL,fc(:,i),gc,fcp,gcp,sig,mode,0,ifail) 
                    if (ifail/=0) then
                        write(*,*) 'coul90: ifail=',ifail
                    endif


                end do
                call cpu_time(t2)
                write(*,11) t2-t1
11              format(20x,'(CPU  time =',F12.8,2x,'seconds)')
                write(*,*) ''

                if(bgauss) then


                call cpu_time(t1)
                write(*,*)"Generating the rotated coulomb wave function on gauss mesh"
                select case(cwftype)

                case(1)

                    do i = 1, numgauss
                        rho = k*gauss_rr(i)*eitheta
                        call COULCC(rho,eta,zlmin,NL,fc_rot_gauss(:,i),gc,fcp,gcp,sig,mode,0,ifail)
                        if (ifail/=0) then
                            write(*,*) 'coul90: ifail=',ifail
                        endif
                    end do

                case(2)
                    do i = 1, numgauss
                        rho = k*gauss_rr(i)*eitheta
                        call compute_coulomb_wave_functions('true', eta,rho,cmplxzero,NL,fc_rot_gauss(:,i),fcp,gc,gcp)
                    end do

                end select

                call cpu_time(t2)
                
                endif

                if(backrot) then
                    write(*,*) 'Generating the back-rotated coulomb wave function on gauss mesh'
                    call cpu_time(t1)
                    do i = 1, numgauss
                        rho = k*gauss_rr(i)
                        call COULCC(rho,eta,zlmin,NL,fc_gauss(:,i),gc,fcp,gcp,sig,mode,0,ifail)
                        if (ifail/=0) then
                            write(*,*) 'coul90: ifail=',ifail
                        endif
                    end do
                    call cpu_time(t2)
                    write(*,11) t2-t1
                endif


            end subroutine


        


        end module