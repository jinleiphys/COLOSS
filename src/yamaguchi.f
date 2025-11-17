        module yamaguchi

        use precision
        use constants
        use mesh
        use matrix_element
        use slove_eigen

        implicit none

        contains
        subroutine nonlocal_test()

        complex*16,dimension(1:nr,1:nr) :: Hyama,Tyama,Vyama,Cyama
        complex*16,dimension(1:nr) :: Byama, Xyama

        complex*16 :: ri,rj

        complex*16 :: fsc,fborn,smat
        real*8 :: phase_angle


        integer :: i,j
        write(*,*) "const is", hbarc**2/amu
        call generate_T0(Tyama,0,mu)

        do i=1,nr
        do j=1,nr
                Vyama(i,j) = Vyama(i,j) + sqrt(mesh_rw(i)*mesh_rw(j))*yama_pot(mesh_rr(i)*eitheta, mesh_rr(j)*eitheta)
        end do
        end do
        Vyama = Vyama * eitheta

        Hyama = Vyama + Tyama

        call cal_N0()

        do i=1,nr
                Byama(i) = 0.d0
                do j=1,nr
                        Byama(i) = Byama(i) + mesh_rw(j)
     &                         *sqrt(mesh_rw(i)) *yama_pot(mesh_rr(i)*eitheta,mesh_rr(j)*eitheta)
     &                         *fc_rotated(0,j)
                end do
        end do

        Byama = Byama * sqrt(eitheta)
        Byama = Byama * eitheta

        Xyama = Byama
        Cyama = ecm*Nmat - Hyama

        call z_lineq(nr,Cyama,Xyama)

        fborn = 0d0
        do i=1,nr
        do j=1,nr
                ri = mesh_rr(i)
                rj = mesh_rr(j)
                fborn = fborn + mesh_rw(i)*mesh_rw(j)*fc(0,i)*fc(0,j)*yama_pot(ri,rj)
        end do
        end do
        fborn = -fborn/ecm
        write(*,*) "fborn is", fborn

        fsc = 0d0
        do i=1,nr
                fsc = fsc + Xyama(i)*Byama(i)
        end do
        fsc = -fsc/ecm
        write(*,*) "fsc is", fsc

        smat = 1.d0 + 2.d0*iu*k*(fborn + fsc)

        phase_angle = atan2(aimag(smat), real(smat))/2d0
        phase_angle = phase_angle*180d0/pi

        write(*,*) "S-matrix is", smat
        write(*,*) "delta is",phase_angle
        write(*,*) "abs is", abs(smat) 

        stop

        end subroutine

        function yama_pot(r1,r2)
                complex*16 :: yama_pot
                complex*16 :: r1,r2

                real*8 :: al=0.2316053d0
                real*8 :: be=1.3918324d0
                real*8 :: hbarcd2mu = 41.472d0

                yama_pot = -hbarcd2mu*2d0*be*(al+be)**2*exp(-be*(r1+r2))

        end function
        
        end module