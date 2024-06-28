        program COLOSS
        
        use mesh
        use system
        use constants
        use precision
        use generate_laguerre
        use input
        use bound
        use scatt
        use coulfunc
        use npcc

        implicit none

        integer :: ll
        real*8 :: t1,t2
        
        call read_input()! read the input from the standard IO

        call init_laguerre_mesh()! generate the laguerre mesh and the scaled mesh

        call get_pot_para(elab,input_pot)!initialize the parameter of optical potential

        call basis_func()!generate the lagrange function basis

        call init_coul()!generate the coulomb function to caculate scatt

        allocate(scatt_amp_nuc_l(0:jmax))

        !call slove_np_cc(ecm)

        select case(method)
        
        case(1)
                call cpu_time(t1)
                write(*,*) "Using Linear Equation Method to Solve"
                write(*,200)
                write(*,201)
                do ll = jmin, jmax
                        call solve_scatt(ll,input_pot)
                end do
                call cpu_time(t2)
        case(2)
                call cpu_time(t1)
                write(*,*) "Using Green's Function Method to Solve"
                write(*, 202)
                write(*,201)
                do ll = jmin, jmax
                        call solve_bound(ll,input_pot)! generate the matrix value and solve the eigenvalue problem
                        call solve_scatt_green(ll)! caclulate the scatt amplitude
                end do
                call cpu_time(t2)
        end select
200   FORMAT(' l  |   S-matrix (real, imag)   | Partial Wave Reaction Cross Section (mb)')
201   FORMAT(' -----------------------------------------------------------------------')
202   FORMAT(' l  | S-matrix (real, imag)     | ')
        write(*,*) '-----------------------------------------------------'
        write(*,11) t2-t1
11              format(20x,'(CPU  time =',F12.6,2x,'seconds)')

        call xsec(real(eta),k,scatt_amp_nuc_l,jmax,thetah,theta_n_max)

        call outinfo()

        call printconstants
        call write_verison()

        end program COLOSS