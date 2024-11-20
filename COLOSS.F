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
        use yamaguchi
        use channels

        implicit none

        integer :: ich
        real*8 :: t1,t2,t0,te

        call cpu_time(t0)
        
        call read_input()! read the input from the standard IO

        call generate_channels()!generate the channel index

        call init_laguerre_mesh()! generate the laguerre mesh and the scaled mesh

        call get_pot_para(elab,input_pot)!initialize the parameter of optical potential

        call basis_func()!generate the lagrange function basis

        call init_coul()!generate the coulomb function to caculate scatt

        allocate(scatt_amp_nuc_channel(1:channel_index%ch_numbers))

        !call slove_np_cc(ecm)
        !call nonlocal_test()

        write(60,*) "& S-matrix for differnt angular momentum channels {LSJ}"
        write(61,*) "& Scattering amplitude f for differnt angular momentum channels {LSJ}"

        write(60,300) "&  Real(S)", "  Aimag(S)"
        write(61,300) "&  Real(f)", "  Aimag(f)"
300     FORMAT(a10,2x,a10,"  angular momentum channel") 

        select case(method)
        
        case(1)
                call cpu_time(t1)
                write(*,*) "Using Linear Equation Method to Solve"
                write(*,200) "L", "S", "J", "Real(Smat)", "Imag(Smat)", "Partial Wave Reac_Xsec (mb)"
                write(*,201)
                do ich = 1, channel_index%ch_numbers
                        call solve_scatt(ich,input_pot)
                end do
                call cpu_time(t2)
        case(2)
                if(nonlocal) then 
                        write(*,*) "For nonlocal potential, Please use linear equation method"
                        stop
                endif
                call cpu_time(t1)
                write(*,*) "Using Green's Function Method to Solve"
                write(*,200) "L", "S", "J", "Real(Smat)", "Imag(Smat)", "Partial Wave Reac_Xsec (mb)"
                write(*,201)
                do ich = 1, channel_index%ch_numbers
                        call solve_bound(ich,input_pot)! generate the matrix value and solve the eigenvalue problem
                        call solve_scatt_green(ich)! caclulate the scatt amplitude
                end do
                call cpu_time(t2)
        end select
        write(*,201)

200   FORMAT(A3,3X,A3,2X,A5,' |  (',A10,', ',A10,')  | ',A28)
201   FORMAT(' ----------------------------------------------------------------------------')
202   FORMAT(' l  | S-matrix (real, imag)     | ')

        write(*,11) t2-t1
11    FORMAT(20x,'(CPU  time =',F12.6,2x,'seconds)')
12    FORMAT(10x,'Total CPU  time =',F12.6,2x,'seconds')

        call xsec(real(eta),k,scatt_amp_nuc_channel,lmax,thetah,theta_n_max)

        call outinfo()
        call cpu_time(te)
        write(*,12) te-t0

        call printconstants
        call write_verison()

        end program COLOSS