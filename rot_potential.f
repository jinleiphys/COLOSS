        module rot_pot

            use constants
            use precision
            use system
            use pot_class
            use mesh
            use generate_laguerre


            implicit none


            complex*16,dimension(:), allocatable :: V_nuc!rotated nuclear potential on laguerre mesh
            complex*16,dimension(:), allocatable :: V_coul!rotated coulomb potential on laguerre mesh
            complex*16,dimension(:), allocatable :: Vcoul_origin!original nuclear potential on laguerre mesh
            complex*16,dimension(:), allocatable :: V_nuc_origin!original nuclear potential laguerre mesh
            complex*16,dimension(:), allocatable :: Vcoul_gauss!rotated coulomb potential on gauss mesh
            complex*16,dimension(:), allocatable :: Vnuc_gauss!rotated nuclear potential on gauss mesh

            contains

            subroutine get_pot_para(E,para)
                real*8 :: E
                type(pot_para) :: para

                real*8 :: mass
                real*8 :: charge

                para%a1 = massp
                para%a2 = masst
                para%z1 = zp
                para%z2 = zt
                para%vv=0d0; para%rvv=0d0; para%avv=0d0
                
                para%wv=0d0; para%rw=0d0; para%aw=0d0
                
                para%vs=0d0; para%rvs=0d0; para%avs=0d0
                
                para%ws=0d0; para%rws=0d0; para%aws=0d0
                
                para%rc=0d0




                    para%vv=vv
                    para%rvv=rv
                    para%avv=av

                    para%wv=wv
                    para%rw=rw
                    para%aw=aw

                    para%vs=vs
                    para%rvs=rvs
                    para%avs=avs

                    para%ws=ws
                    para%rws=rws
                    para%aws=aws

                    para%rc=rc







                write(*,*) '------------- Optical Model Potential -------------'
                write(*,20) masst, zt, elab, para%rc
                write(*,111) para%vv,para%rvv,para%avv,para%wv,para%rw,para%aw
                write(*,112) para%vs,para%rvs,para%avs,para%ws,para%rws,para%aws

20              format(' Optical Potential for A = ',f5.1,' Z = ',f4.1,' at ',f7.3,' MeV ','rc =',f7.3)
111             format('  Vv     rvv    avv    Wv     rw     aw', /,6f7.3,/)
112             format('  Vs     rvs    avs    Ws     rws    aws',/,6f7.3,/)
            end subroutine




            subroutine rot_V_nuc(para)

                type(pot_para) :: para

                integer :: ir 
                complex*16 :: r
                real*8 :: a13

                complex*16 :: real_volume, img_volume, real_surf, img_surf

                a13 = masst**(1./3.)

                if(allocated(V_nuc)) deallocate(V_nuc)
                allocate(V_nuc(1:nr))
                V_nuc = 0.d0

                if(allocated(V_nuc_origin)) deallocate(V_nuc_origin)
                allocate(V_nuc_origin(1:nr))
                V_nuc_origin = 0.d0
                

                do ir = 1, nr
                    r  = mesh_rr(ir)*eitheta
                    real_volume = -w_s(r,para%vv,para%rvv*a13,para%avv)
                    img_volume = -w_s(r,para%wv,para%rw*a13,para%aw)
                    real_surf = 4.0d0*para%avs*dws(r,para%vs,para%rvs*a13,para%avs) !usually thers is no real_surf
                    img_surf = 4.0d0*para%aws*dws(r,para%ws, para%rws*a13, para%aws)
                    V_nuc(ir) = real_volume + real_surf + iu*(img_surf+img_volume)                               
                    !write(901,*) mesh_rr(ir), real(V_nuc(ir)), aimag(V_nuc(ir))
                end do

                do ir = 1, nr
                    r  = mesh_rr(ir)
                    real_volume = -w_s(r,para%vv,para%rvv*a13,para%avv)
                    img_volume = -w_s(r,para%wv,para%rw*a13,para%aw)
                    real_surf = 4.0d0*para%avs*dws(r,para%vs,para%rvs*a13,para%avs) !usually thers is no real_surf
                    img_surf = 4.0d0*para%aws*dws(r,para%ws, para%rws*a13, para%aws)
                    V_nuc_origin(ir) = real_volume + real_surf + iu * (img_surf+img_volume)                              
                    !write(902,*) mesh_rr(ir), real(V_nuc_origin(ir)), aimag(V_nuc_origin(ir))
                end do

            end subroutine


            subroutine rot_V_coul(z12,rc)

                real*8 :: z12
                real*8 :: rc
                integer :: ir
                complex*16 :: rr


                if(allocated(V_coul)) deallocate(V_coul)
                allocate(V_coul(1:nr))
                V_coul = 0.d0

                if(allocated(Vcoul_origin)) deallocate(Vcoul_origin)
                allocate(Vcoul_origin(1:nr))
                Vcoul_origin = 0.d0

                do ir = 1,nr 
                    rr = mesh_rr(ir)*eitheta
                    V_coul(ir) = VCOUL(rr,z12,rc)
                    rr = mesh_rr(ir)
                    Vcoul_origin(ir) = VCOUL(rr,z12,rc) 
                    !write(801,*) mesh_rr(ir),real(V_coul(ir)),aimag(V_coul(ir))
                    !write(802,*) mesh_rr(ir),real(Vcoul_origin(ir)),aimag(Vcoul_origin(ir))
                end do



            end subroutine

            subroutine coul_mat_gauss(z12,rc,coulmat)
c           This subroutine calculate the rotated coulomb potential
c           on the gauss mesh, and then evalutate the matrix element
c           by doing the numerical integral with Gauss quadtrature method.

                real*8 :: rc
                real*8 :: z12
                complex*16,dimension(1:nr,1:nr) :: coulmat

                integer :: ir
                complex*16 :: rr
                complex*16 :: coulij
                integer :: i,j,i_cor


                if(allocated(Vcoul_gauss)) deallocate(Vcoul_gauss)
                allocate(Vcoul_gauss(1:numgauss))

                if(backrot) then

                    do ir = 1, numgauss
                        rr = gauss_rr(ir)
                        Vcoul_gauss(ir) = VCOUL(rr,z12,rc)
                    end do

                else

                    do ir = 1, numgauss
                        rr = gauss_rr(ir)* eitheta
                        Vcoul_gauss(ir) = VCOUL(rr,z12,rc)
                    end do

                endif



                  

                !calculate the coulomb matrix
                coulmat = 0.d0
                if(backrot) then
                
                do i = 1, nr
                    do j = 1, nr
                        coulij = 0.d0
                        do i_cor = 1, numgauss!iteration for gauss quadrature
                            rr = gauss_rr(i_cor)
                            coulij = coulij + 
     &                            gauss_rw(i_cor)*Vcoul_gauss(i_cor)
     &                            *lag_func_br(i_cor,i)*lag_func_br(i_cor,j)
                        end do!end of iteration for gauss quadrature
                        coulmat(i,j) = coulij!/eitheta
                    end do
                end do
                coulmat = coulmat/eitheta

                else

                do i = 1, nr
                    do j = 1, nr
                        coulij = 0.d0
                        do i_cor = 1, numgauss!iteration for gauss quadrature
                            rr = gauss_rr(i_cor)
                            coulij = coulij + 
     &                            gauss_rw(i_cor)*Vcoul_gauss(i_cor)
     &                              *lag_func(i_cor,i)*lag_func(i_cor,j)
                        end do!end of iteration for gauss quadrature
                        coulmat(i,j) = coulij!/eitheta
                    end do
                end do

                endif

                


            end subroutine


            subroutine nuc_mat_gauss(para,nucmat)
c           This subroutine calculate the rotated optical potential
c           on the gauss mesh, and then evalutate the matrix element
c           by doing the numerical integral with Gauss quadtrature method.

                type(pot_para) :: para
                complex*16, dimension(1:nr,1:nr) :: nucmat

                integer :: ir
                complex*16 :: r
                real*8 :: a13
                complex*16 :: real_volume, img_volume, real_surf, img_surf

                complex*16 :: nucij
                integer :: i,j,i_cor
                complex*16 :: rr

                a13 = masst**(1d0/3d0)

                if(allocated(Vnuc_gauss)) deallocate(Vnuc_gauss)
                allocate(Vnuc_gauss(1:numgauss))

                if(backrot) then
    
                !calculate the rotated optical potential on the gauss mesh
                do ir = 1, numgauss
                    r  = gauss_rr(ir)
                    real_volume = -w_s(r,para%vv,para%rvv*a13,para%avv)
                    img_volume = -w_s(r,para%wv,para%rw*a13,para%aw)
                    real_surf = 4.0d0*para%avs*dws(r,para%vs,para%rvs*a13,para%avs) !usually thers is no real_surf
                    img_surf = 4.0d0*para%aws*dws(r,para%ws, para%rws*a13, para%aws)
                    Vnuc_gauss(ir) = real_volume + real_surf + iu*(img_surf+img_volume)                               
                end do

                else

                do ir = 1, numgauss
                    r  = gauss_rr(ir)*eitheta
                    real_volume = -w_s(r,para%vv,para%rvv*a13,para%avv)
                    img_volume = -w_s(r,para%wv,para%rw*a13,para%aw)
                    real_surf = 4.0d0*para%avs*dws(r,para%vs,para%rvs*a13,para%avs) !usually thers is no real_surf
                    img_surf = 4.0d0*para%aws*dws(r,para%ws, para%rws*a13, para%aws)
                    Vnuc_gauss(ir) = real_volume + real_surf + iu*(img_surf+img_volume)                               
                end do

                endif

                ! calculate the nuclear potential matrix
                nucmat = 0.d0

                if(backrot) then

                do i = 1, nr
                    do j = 1, nr
                        nucij = 0.d0
                        do i_cor = 1, numgauss!iteration for gauss quadrature
                            rr = gauss_rr(i_cor)
                            nucij = nucij + 
     &                            gauss_rw(i_cor)*Vnuc_gauss(i_cor)
     &                              *lag_func_br(i_cor,i)*lag_func_br(i_cor,j)
                        end do!end of iteration for gauss quadrature
                        nucmat(i,j) = nucij
                    end do
                end do
                nucmat = nucmat/eitheta

                else 

                do i = 1, nr
                    do j = 1, nr
                        nucij = 0.d0
                        do i_cor = 1, numgauss!iteration for gauss quadrature
                            rr = gauss_rr(i_cor)
                            nucij = nucij + 
     &                            gauss_rw(i_cor)*Vnuc_gauss(i_cor)
     &                              *lag_func(i_cor,i)*lag_func(i_cor,j)
                        end do!end of iteration for gauss quadrature
                        nucmat(i,j) = nucij
                    end do
                end do

                endif

            end subroutine

            


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c              functions for different types of potential
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c *** Woods-Saxon (volume term)
        function w_s(r,v0,r0,a)
            implicit none
            real*8 v0,r0,a
            complex*16 :: r
            complex*16 :: w_s,aux
            w_s=0d0
            if (abs(v0) < 1e-6) return
            if (a>1e-6) then
                aux=exp(-(r-r0)/a)
                w_s=v0/(1d0+1d0/aux)
            else
                write(0,*)'WS: a too small!'
                stop
            end if
        return
        end function
c-----------------------------------------------------------------------
c *** Gaussian
        function gausspot(r,v0,r0,a)
            implicit none
            real*8 v0,r0,a
            complex*16 :: r,gausspot
            gausspot = 0.d0
            if(abs(v0) < 1e-6) return
            if (a.gt.1e-6) then
                gausspot=V0*exp(-(r-r0)**2/a**2)
            else
               write(*,*)'a too small in gausspot!'
               stop
            endif
            return
        end function

c-----------------------------------------------------------------------
c Coulomb potential
        FUNCTION VCOUL(R,z12,Rc)
            use constants
            implicit none
            real*8 :: rc,rc2,z12
            complex*16 :: r,vcoul,aux

            RC2=RC*2d0
            aux=e2*Z12
            vcoul=0.d0
            if (z12.lt.1e-4) return
            if (rc.lt.1e-6) rc=1e-6

            IF(abs(R).GT.RC)GO TO 1
            VCOUL=AUX*(3.-(R/RC)**2)/RC2
            RETURN
1           VCOUL=AUX/R
            RETURN
        END function
c-----------------------------------------------------------------------
c *** WS derivative
        function dws(r,v0,r0,a)
            implicit none
            real*8 v0,r0,a
            complex*16 :: r,dws,aux
            dws = 0.d0
            if (abs(r)<1e-6) r=1e-6
            if(abs(v0)< 1e-6) return
            if (a>1e-6) then
                aux = exp(-(r-r0)/a)
                dws = -v0*aux/(1d0+aux)**2/a
            else
                write(0,*)'derivative WS: a too small!';stop
            endif
            return
        end function

        end module rot_pot