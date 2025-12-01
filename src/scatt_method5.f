        module scatt_method5
c-----------------------------------------------------------------------
c       Module for Method 5: Scattering using Baye's exact D and T
c       matrices with CORRECT usage on expansion coefficients c_j
c
c       Key insight: Baye's D and T matrices act on expansion
c       coefficients c_j, NOT on function values phi(x_j)!
c
c       For x-regularized basis: f_j(x) = alpha_j * x * L_j(x)
c       Wave function: phi(x) = sum_j c_j * f_j(x)
c       At mesh points: phi(x_j) = c_j * alpha_j * x_j
c       So: c_j = phi(x_j) / (alpha_j * x_j)
c
c       Baye's matrix equation:
c         sum_j [T_ij + (k^2 - l(l+1)/x_i^2 - U(x_i))*delta_ij] c_j = b_i
c
c       Reference: D. Baye, Physics Reports 565 (2015) 1-107
c                  Section 3.4.5, Equations (3.122)-(3.126)
c-----------------------------------------------------------------------
            use mesh
            use system
            use precision
            use constants
            use matrix_element
            use pot_class
            use channels
            use generate_laguerre
            use rot_pot
            use scatt, only: scatt_amp_nuc_channel
            use coulfunc
            use gauss

            implicit none

            ! Lagrange-Legendre mesh arrays
            real*8, dimension(:), allocatable :: leg5_x      ! Mesh points in (0,1)
            real*8, dimension(:), allocatable :: leg5_r      ! Mesh points in (0,R)
            real*8, dimension(:), allocatable :: leg5_w      ! Gauss weights on (0,1)
            real*8, dimension(:), allocatable :: leg5_alpha  ! Normalization alpha_j

            ! Baye's exact D and T matrices (act on coefficients c_j)
            real*8, dimension(:,:), allocatable :: D1_baye   ! Baye's D matrix
            real*8, dimension(:,:), allocatable :: T_baye    ! Baye's T matrix

            contains

c-----------------------------------------------------------------------
            subroutine init_legendre_mesh5(N, R)
c           Initialize shifted Legendre mesh on (0, R)
c           and compute normalization factors alpha_j
c-----------------------------------------------------------------------
            implicit none
            integer, intent(in) :: N
            real*8, intent(in) :: R

            real*8, dimension(N) :: t_gauss, w_gauss
            integer :: j

            if (allocated(leg5_r)) deallocate(leg5_r)
            if (allocated(leg5_w)) deallocate(leg5_w)
            if (allocated(leg5_x)) deallocate(leg5_x)
            if (allocated(leg5_alpha)) deallocate(leg5_alpha)
            allocate(leg5_r(N), leg5_w(N), leg5_x(N), leg5_alpha(N))

            ! Get Gauss-Legendre points on (-1, 1)
            call gauleg(N, -1.d0, 1.d0, t_gauss, w_gauss)

            ! Map from (-1, 1) to (0, 1): x_j = (t_j + 1)/2
            do j = 1, N
                leg5_x(j) = (t_gauss(j) + 1.d0) / 2.d0
                ! Gauss weight on (0,1): lambda_j = w_j^GL / 2
                leg5_w(j) = w_gauss(j) / 2.d0
            end do

            ! Scale to (0, R): r_j = R * x_j
            do j = 1, N
                leg5_r(j) = R * leg5_x(j)
            end do

            ! Compute f_j(x_j) for x-regularized basis
            ! From Baye (2015) Eq. (3.122):
            ! f_j(x) = (-1)^{N-j} * sqrt((1-x_j)/x_j) * x * P_N(2x-1) / (x - x_j)
            !
            ! At mesh point x_j (using L'Hopital's rule):
            !   f_j(x_j) = 1 / sqrt(lambda_j)
            !
            ! The wave function phi = sum_j c_j * f_j satisfies:
            !   phi(x_j) = c_j * f_j(x_j) = c_j / sqrt(lambda_j)
            !
            ! Therefore: c_j = phi(x_j) * sqrt(lambda_j)
            !
            ! We store f_j(x_j) = 1/sqrt(lambda_j) in leg5_alpha
            do j = 1, N
                leg5_alpha(j) = 1.d0 / sqrt(leg5_w(j))
            end do

            end subroutine

c-----------------------------------------------------------------------
            subroutine generate_baye_matrices(N, R)
c           Generate Baye's exact D and T matrices for x-regularized basis
c
c           From Baye (2015) Section 3.4.5:
c
c           D_{i!=j} = (-1)^{i-j} * sqrt(x_i(1-x_j)/(x_j(1-x_i))) / (x_i-x_j)
c           D_{ii} = 1 / (2*x_i*(1-x_i))
c
c           T_{i!=j} = (-1)^{i-j} * (x_i+x_j-2*x_i^2) / (x_j*(x_j-x_i)^2)
c                      * sqrt(x_j*(1-x_j) / (x_i*(1-x_i)^3))
c           T_{ii} = [N(N+1)*x_i*(1-x_i) - 3*x_i + 1] / [3*x_i^2*(1-x_i)^2]
c
c           These are on the UNIT interval (0,1).
c           For physical interval (0,R):
c             D -> D/R  (first derivative scales as 1/R)
c             T -> T/R^2 (second derivative scales as 1/R^2)
c
c           INPUT:
c               N: Number of mesh points
c               R: Outer boundary
c-----------------------------------------------------------------------
            implicit none
            integer, intent(in) :: N
            real*8, intent(in) :: R

            integer :: i, j
            real*8 :: xi, xj, sign_ij, sqrt_fac, numer, denom

            if (allocated(D1_baye)) deallocate(D1_baye)
            if (allocated(T_baye)) deallocate(T_baye)
            allocate(D1_baye(N, N), T_baye(N, N))

            ! Build Baye D matrix (d/dx on unit interval)
            D1_baye = 0.d0
            do i = 1, N
                xi = leg5_x(i)
                do j = 1, N
                    xj = leg5_x(j)
                    if (i /= j) then
                        sign_ij = (-1.d0)**(i - j)
                        sqrt_fac = sqrt(xi*(1.d0-xj) / (xj*(1.d0-xi)))
                        D1_baye(i,j) = sign_ij * sqrt_fac / (xi - xj)
                    else
                        D1_baye(i,i) = 1.d0 / (2.d0*xi*(1.d0-xi))
                    endif
                end do
            end do
            ! Scale to physical interval: d/dr = (1/R) d/dx
            D1_baye = D1_baye / R

            ! Build Baye T matrix (-d^2/dx^2 on unit interval)
            ! Note: T represents the kinetic energy operator -d^2/dx^2
            T_baye = 0.d0
            do i = 1, N
                xi = leg5_x(i)
                do j = 1, N
                    xj = leg5_x(j)
                    if (i == j) then
                        numer = dble(N*(N+1))*xi*(1.d0-xi) - 3.d0*xi + 1.d0
                        denom = 3.d0 * xi*xi * (1.d0 - xi)**2
                        T_baye(i,i) = numer / denom
                    else
                        sign_ij = (-1.d0)**(i - j)
                        numer = xi + xj - 2.d0*xi*xi
                        denom = xj * (xj - xi)**2
                        sqrt_fac = sqrt(xj*(1.d0-xj) / (xi*(1.d0-xi)**3))
                        T_baye(i,j) = sign_ij * numer / denom * sqrt_fac
                    endif
                end do
            end do
            ! Scale: T is -d^2/dx^2, we need d^2/dr^2 = (1/R^2)*d^2/dx^2
            ! So d^2/dr^2 matrix = -T/R^2
            T_baye = T_baye / (R * R)

            end subroutine

c-----------------------------------------------------------------------
            subroutine compute_gamma_s5(k_wave, rN, l, gamma_s, eta_som)
c           Compute the logarithmic derivative of outgoing Hankel function
c           gamma_s = H_l^+'(k*r_N) / H_l^+(k*r_N)
c-----------------------------------------------------------------------
            implicit none
            real*8, intent(in) :: k_wave, rN, eta_som
            integer, intent(in) :: l
            complex*16, intent(out) :: gamma_s

            real*8 :: x
            real*8 :: xlmin
            integer :: ifail, kfn
            real*8, dimension(0:l) :: fc_loc, gc_loc, fcp_loc, gcp_loc
            complex*16 :: hhat, hhatp

            x = k_wave * rN
            xlmin = 0.d0
            kfn = 0

            call COUL90(x, eta_som, xlmin, l, fc_loc, gc_loc,
     &                  fcp_loc, gcp_loc, kfn, ifail)
            if (ifail /= 0) then
                write(*,*) 'compute_gamma_s5: COUL90 ifail=', ifail
            endif

            ! u^{C(+)} = (G + iF) * e^{-i*sigma_l}
            ! But for logarithmic derivative, the phase cancels:
            ! gamma = u'/u = (G'+iF')/(G+iF) (phase factor cancels)
            hhat = gc_loc(l) + iu * fc_loc(l)
            hhatp = k_wave * (gcp_loc(l) + iu * fcp_loc(l))

            if (abs(hhat) < 1.d-20) then
                write(*,*) "Warning: H^+ very small in compute_gamma_s5"
                gamma_s = cmplx(0.d0, k_wave, kind=8)
            else
                gamma_s = hhatp / hhat
            endif

            end subroutine

c-----------------------------------------------------------------------
            subroutine solve_scatt_method5(ich, para)
c           Solve scattering problem using Baye's exact D and T matrices
c           with CORRECT action on expansion coefficients c_j
c
c           Key relations:
c             phi(x) = sum_j c_j * f_j(x)
c             phi(x_j) = c_j * alpha_j * x_j
c             c_j = phi(x_j) / (alpha_j * x_j)
c
c           Scattering equation in coefficient space:
c             sum_j [-T_ij + (k^2 - l(l+1)/r_i^2 - U_i)*delta_ij] c_j = b_i
c
c           where b_i = U_i * F_l(k*r_i) / (alpha_i * r_i)
c           (source term converted to coefficient space)
c
c           Boundary condition at r = R:
c             phi'(R) = gamma_s * phi(R)
c
c           Using phi(R) = sum_j c_j * f_j(R) and similar for phi'(R)
c-----------------------------------------------------------------------
            use slove_eigen

            implicit none
            integer, intent(in) :: ich
            type(pot_para), intent(in) :: para

            integer :: l
            real*8 :: S, J
            integer :: ir, jr
            integer :: N_leg

            complex*16, dimension(:,:), allocatable :: M_matrix
            complex*16, dimension(:), allocatable :: b_vec, c_vec

            complex*16 :: vmod, vmod_short
            complex*16 :: f_l, smat
            real*8 :: reac_xsec

            real*8 :: a13, rrc, z12
            real*8 :: R_outer

            complex*16 :: gamma_s
            real*8 :: eta_real

            real*8 :: x_mesh, xlmin_mesh
            real*8, dimension(:), allocatable :: fc_loc, gc_loc
            real*8, dimension(:), allocatable :: fcp_loc, gcp_loc
            integer :: kfn_mesh, ifail_mesh

            complex*16 :: phi_R, phi_R_deriv, hhat_R
            real*8 :: coeff_kin

            ! For scattering amplitude extraction
            complex*16 :: f_l_direct, f_l_integral
            complex*16 :: f_born, f_sc_int
            complex*16 :: psi_sc_ir

            ! For boundary condition
            real*8 :: f_j_at_R, f_j_deriv_at_R
            real*8 :: L_j_at_R, L_j_deriv_at_R
            integer :: i, kk
            real*8 :: prod

            ! Get channel quantum numbers
            l = channel_index%L(ich)
            S = channel_index%S(ich)
            J = channel_index%J(ich)

            R_outer = Rmax
            N_leg = nr

            allocate(fc_loc(0:lmax), gc_loc(0:lmax))
            allocate(fcp_loc(0:lmax), gcp_loc(0:lmax))

            ! Initialize mesh and compute alpha_j
            call init_legendre_mesh5(N_leg, R_outer)

            ! Generate Baye's exact D and T matrices
            call generate_baye_matrices(N_leg, R_outer)

            a13 = para%a2**(1./3.) + para%a1**(1./3.)
            rrc = a13 * para%rc
            z12 = zt * zp

            eta_real = real(eta)
            call compute_gamma_s5(k, R_outer, l, gamma_s, eta_real)


            allocate(M_matrix(N_leg, N_leg))
            allocate(b_vec(N_leg), c_vec(N_leg))

            M_matrix = cmplx(0.d0, 0.d0, kind=8)
            b_vec = cmplx(0.d0, 0.d0, kind=8)

            coeff_kin = hbarc**2 / (2.d0 * mu)

            ! ============================================================
            ! BUILD MATRIX AND SOURCE VECTOR IN COEFFICIENT SPACE
            ! ============================================================
            ! Equation: [-d^2/dr^2 + l(l+1)/r^2 + U - k^2] phi = -U * F_l
            !
            ! In Baye's formulation with T = -d^2/dr^2:
            !   [T + l(l+1)/r^2 + U - k^2] c = b
            !
            ! Matrix: M_ij = T_ij + [l(l+1)/r_i^2 + U_i - k^2] * delta_ij
            ! Source: b_i = -U_i * F_l(k*r_i) / (alpha_i * r_i)
            !         (converted from phi space to c space)

            do ir = 1, N_leg - 1
                x_mesh = k * leg5_r(ir)
                xlmin_mesh = 0.d0
                kfn_mesh = 0

                call COUL90(x_mesh, eta_real, xlmin_mesh, lmax,
     &                      fc_loc, gc_loc, fcp_loc, gcp_loc,
     &                      kfn_mesh, ifail_mesh)

                ! Get V_short = V_nuc + V_coul_finite - V_coul_point
                call compute_potential_at_r5(leg5_r(ir), para, ich, vmod_short)

                ! For matrix: use full potential V_nuc + V_coul_finite
                ! vmod = V_short + V_coul_point = V_nuc + V_coul_finite
                vmod = vmod_short
                if (abs(z12) > 1.d-10) then
                    vmod = vmod + e2 * z12 / leg5_r(ir)
                endif
                vmod = vmod / coeff_kin  ! U_full = 2*mu*V_full/hbar^2

                do jr = 1, N_leg
                    ! Baye's T matrix (note: T represents -d^2/dr^2)
                    ! We need d^2/dr^2 in Schrodinger eq, so use -T
                    M_matrix(ir, jr) = -T_baye(ir, jr)

                    if (ir == jr) then
                        M_matrix(ir, ir) = M_matrix(ir, ir) + k*k
                        M_matrix(ir, ir) = M_matrix(ir, ir) -
     &                      dble(l*(l+1)) / leg5_r(ir)**2
                        M_matrix(ir, ir) = M_matrix(ir, ir) - vmod
                    endif
                end do

                ! Source term: use SHORT-RANGE potential Ṽ_N = V_N + V_C^S
                ! From Eq. (13): [E - H] ψ^sc = e^{iσ_l} Ṽ_N F_l
                ! Try WITHOUT exp(i*sigma_l) first to match no-Coulomb structure
                ! b = U_short * F_l * sqrt(lambda)
                vmod_short = vmod_short / coeff_kin  ! U_short = 2*mu*V_short/hbar^2
                b_vec(ir) = vmod_short * fc_loc(l) * sqrt(leg5_w(ir))
            end do

            ! ============================================================
            ! BOUNDARY CONDITION AT r = R (Row N)
            ! ============================================================
            ! phi'(R) - gamma_s * phi(R) = 0
            !
            ! For Baye's x-regularized basis f_j(x):
            !   f_j(x) = (-1)^{N-j} * sqrt((1-x_j)/x_j) * x*P_N(2x-1)/(x-x_j)
            !
            ! This can be written as: f_j(r) = (r/r_j) * L_j(r) / sqrt(lambda_j)
            ! where L_j(r) is standard Lagrange polynomial
            !
            ! At r = R:
            !   f_j(R) = (R/r_j) * L_j(R) / sqrt(lambda_j)
            !   f'_j(R) = [L_j(R)/r_j + (R/r_j)*L'_j(R)] / sqrt(lambda_j)

            do jr = 1, N_leg
                ! Compute L_j(R) on physical r coordinates
                L_j_at_R = 1.d0
                do i = 1, N_leg
                    if (i /= jr) then
                        L_j_at_R = L_j_at_R * (R_outer - leg5_r(i)) /
     &                                        (leg5_r(jr) - leg5_r(i))
                    endif
                end do

                ! Compute L'_j(R)
                L_j_deriv_at_R = 0.d0
                do i = 1, N_leg
                    if (i /= jr) then
                        prod = 1.d0
                        do kk = 1, N_leg
                            if (kk /= jr .and. kk /= i) then
                                prod = prod * (R_outer - leg5_r(kk)) /
     &                                        (leg5_r(jr) - leg5_r(kk))
                            endif
                        end do
                        L_j_deriv_at_R = L_j_deriv_at_R +
     &                      prod / (leg5_r(jr) - leg5_r(i))
                    endif
                end do

                ! f_j(R) = (R/r_j) * L_j(R) / sqrt(lambda_j)
                f_j_at_R = (R_outer / leg5_r(jr)) * L_j_at_R /
     &                     sqrt(leg5_w(jr))

                ! f'_j(R) = [L_j(R)/r_j + (R/r_j)*L'_j(R)] / sqrt(lambda_j)
                f_j_deriv_at_R = (L_j_at_R / leg5_r(jr) +
     &              (R_outer / leg5_r(jr)) * L_j_deriv_at_R) /
     &              sqrt(leg5_w(jr))

                ! BC: f'_j(R) - gamma_s * f_j(R) = 0
                M_matrix(N_leg, jr) = f_j_deriv_at_R - gamma_s * f_j_at_R
            end do
            b_vec(N_leg) = cmplx(0.d0, 0.d0, kind=8)

            ! Debug: print b_vec before solve
            if (ich == 1) then
                write(*,*) "=== DEBUG Method 5 (L=", l, ") ==="
                write(*,*) "k =", k, "eta =", eta_real, "cph(l) =", cph(l)
                write(*,*) "R_outer =", R_outer, "N_leg =", N_leg
                write(*,*) "gamma_s =", gamma_s
                write(*,*) "Sum |b_vec| =", sum(abs(b_vec))
            endif

            ! ============================================================
            ! SOLVE LINEAR SYSTEM FOR COEFFICIENTS c_j
            ! ============================================================
            c_vec = b_vec
            call z_lineq(N_leg, M_matrix, c_vec)

            ! Debug: print c_vec after solve
            if (ich == 1) then
                write(*,*) "Sum |c_vec| =", sum(abs(c_vec))
                write(*,*) "c_vec(1:3) =", c_vec(1), c_vec(2), c_vec(3)
                write(*,*) "c_vec(N-2:N) =", c_vec(N_leg-2), c_vec(N_leg-1),
     &                     c_vec(N_leg)
            endif

            ! ============================================================
            ! EXTRACT SCATTERING AMPLITUDE
            ! ============================================================
            ! phi(R) = sum_j c_j * f_j(R)
            ! f_l = phi(R) / (k * H^+(kR))

            phi_R = cmplx(0.d0, 0.d0, kind=8)
            do jr = 1, N_leg
                ! Compute L_j(R)
                L_j_at_R = 1.d0
                do i = 1, N_leg
                    if (i /= jr) then
                        L_j_at_R = L_j_at_R * (R_outer - leg5_r(i)) /
     &                                        (leg5_r(jr) - leg5_r(i))
                    endif
                end do
                ! f_j(R) = (R/r_j) * L_j(R) / sqrt(lambda_j)
                f_j_at_R = (R_outer / leg5_r(jr)) * L_j_at_R /
     &                     sqrt(leg5_w(jr))
                phi_R = phi_R + c_vec(jr) * f_j_at_R
            end do

            ! Compute H^+ at boundary
            x_mesh = k * R_outer
            call COUL90(x_mesh, eta_real, 0.d0, l,
     &                  fc_loc, gc_loc, fcp_loc, gcp_loc, 0, ifail_mesh)

            ! H^+ = G + iF
            hhat_R = gc_loc(l) + iu * fc_loc(l)

            ! ============================================================
            ! METHOD A: Direct matching at boundary
            ! psi_sc -> k * f_l * O_l^{(+)} (without source phase factor)
            ! f_l^{direct} = phi(R) / (k * H^+(kR))
            ! ============================================================
            f_l_direct = phi_R / (k * hhat_R)

            ! ============================================================
            ! METHOD B: Integral formula
            ! f_born = -1/ecm * integral[V_short * F_l^2] dr
            ! f_sc   = -1/ecm * integral[V_short * F_l * psi_sc] dr
            ! ============================================================

            ! f_born: integral of V_short * F_l^2
            f_born = cmplx(0.d0, 0.d0, kind=8)
            do ir = 1, N_leg - 1
                x_mesh = k * leg5_r(ir)
                call COUL90(x_mesh, eta_real, 0.d0, l,
     &                      fc_loc, gc_loc, fcp_loc, gcp_loc, 0, ifail_mesh)

                ! Get V_short
                call compute_potential_at_r5(leg5_r(ir), para, ich, vmod)

                ! f_born: V_short * F_l^2 * dr
                f_born = f_born + leg5_w(ir) * R_outer *
     &                        vmod * fc_loc(l) * fc_loc(l)
            end do
            f_born = -f_born / ecm

            ! f_sc: integral of V_short * F_l * psi_sc
            ! From Eq. (15): f = -2mu/hbar^2/k^2 * exp(-i*sigma_l) * int[F_l*V*psi]
            ! Since source term has exp(i*sigma_l), psi_sc = exp(i*sigma_l) * psi_sc_tilde
            ! So: f_sc = -1/ecm * exp(-i*sigma_l) * int[V*F_l*psi_sc]
            !         = -1/ecm * exp(-i*sigma_l) * int[V*F_l*exp(i*sigma_l)*psi_sc_tilde]
            !         = -1/ecm * int[V*F_l*psi_sc_tilde]
            ! No additional phase factor needed!
            f_sc_int = cmplx(0.d0, 0.d0, kind=8)
            do ir = 1, N_leg - 1
                x_mesh = k * leg5_r(ir)
                call COUL90(x_mesh, eta_real, 0.d0, l,
     &                      fc_loc, gc_loc, fcp_loc, gcp_loc, 0, ifail_mesh)

                ! Get V_short
                call compute_potential_at_r5(leg5_r(ir), para, ich, vmod)

                ! psi_sc at mesh point: phi(r_i) = c_i / sqrt(lambda_i)
                psi_sc_ir = c_vec(ir) / sqrt(leg5_w(ir))

                ! integral: V_short * F_l * psi_sc * dr
                f_sc_int = f_sc_int + leg5_w(ir) * R_outer *
     &                     vmod * fc_loc(l) * psi_sc_ir
            end do
            f_sc_int = -f_sc_int / ecm

            f_l_integral = f_born + f_sc_int

            ! Output comparison for L=0 channel
            if (ich == 1) then
                write(*,*) "======== Rmax =", R_outer, " L =", l,
     &                     " ========"
                write(*,*) "Direct match: f_l =", f_l_direct
                write(*,*) "Integral:     f_l =", f_l_integral
                write(*,*) "  f_born =", f_born
                write(*,*) "  f_sc   =", f_sc_int
                write(*,*) "Difference: ", abs(f_l_direct - f_l_integral)
            endif

            ! Use integral method result
            f_l = f_l_integral

            scatt_amp_nuc_channel(ich) = f_l

            ! S = 1 + 2ik*f_l (same formula as Method 1)
            smat = 1.d0 + 2.d0 * iu * k * f_l

            reac_xsec = pi/k/k/(2.d0*S+1.d0)*(2.d0*J+1.d0)
     &                 *(1.d0 - abs(smat)**2) * 10.d0

            write(*, 300) l, S, J, real(smat), aimag(smat), reac_xsec
300         FORMAT(I3,3x,F3.1,2x,F5.1,' |  (',F10.6,', ',F10.6,')  | ',
     &             F14.4)

            write(60, 101) real(smat), aimag(smat), l, S, J
            write(61, 101) real(f_l), aimag(f_l), l, S, J
101         FORMAT(F10.6,2x,F10.6,"  (L S J):",I3,3x,F3.1,2x,F5.1)

            deallocate(M_matrix, b_vec, c_vec)
            deallocate(fc_loc, gc_loc, fcp_loc, gcp_loc)

            end subroutine

c-----------------------------------------------------------------------
            subroutine compute_potential_at_r5(r, para, ich, V_short)
c           Compute the short-range potential at radius r
c-----------------------------------------------------------------------
            implicit none
            real*8, intent(in) :: r
            type(pot_para), intent(in) :: para
            integer, intent(in) :: ich
            complex*16, intent(out) :: V_short

            real*8 :: a13, rrc, z12
            complex*16 :: V_nuc_r, V_coul_r
            real*8 :: V_coul_point

            a13 = para%a2**(1./3.) + para%a1**(1./3.)
            rrc = a13 * para%rc
            z12 = zt * zp

            call eval_nuclear_potential5(r, para, ich, V_nuc_r)
            call eval_coulomb_potential5(r, z12, rrc, V_coul_r)

            if (abs(z12) > 1.d-10) then
                V_coul_point = e2 * z12 / r
            else
                V_coul_point = 0.d0
            endif

            V_short = V_nuc_r + V_coul_r - V_coul_point

            end subroutine

c-----------------------------------------------------------------------
            subroutine eval_nuclear_potential5(r, para, ich, V_nuc)
c           Evaluate nuclear optical potential at radius r
c-----------------------------------------------------------------------
            implicit none
            real*8, intent(in) :: r
            type(pot_para), intent(in) :: para
            integer, intent(in) :: ich
            complex*16, intent(out) :: V_nuc

            real*8 :: a13
            real*8 :: Vv_p, rv_p, av_p, Wv_p, rwv_p, awv_p
            real*8 :: Vs_p, rvs_p, avs_p, Ws_p, rws_p, aws_p
            real*8 :: Vso_p, rso_p, aso_p, Wso_p, rwso_p, awso_p
            real*8 :: R_v, R_wv, R_s, R_ws, R_so, R_wso
            real*8 :: f_v, f_wv, f_s, f_ws
            real*8 :: df_s, df_ws, df_so, df_wso
            complex*16 :: V_central, V_surface, V_so_term
            integer :: l
            real*8 :: S_ch, J_ch
            real*8 :: ls_factor

            l = channel_index%L(ich)
            S_ch = channel_index%S(ich)
            J_ch = channel_index%J(ich)

            a13 = para%a2**(1./3.) + para%a1**(1./3.)

            Vv_p = para%vv; rv_p = para%rvv; av_p = para%avv
            Wv_p = para%wv; rwv_p = para%rw; awv_p = para%aw
            Vs_p = para%vs; rvs_p = para%rvs; avs_p = para%avs
            Ws_p = para%ws; rws_p = para%rws; aws_p = para%aws
            Vso_p = para%vsov; rso_p = para%rsov; aso_p = para%asov
            Wso_p = para%vsow; rwso_p = para%rsow; awso_p = para%asow

            R_v = rv_p * a13
            R_wv = rwv_p * a13
            R_s = rvs_p * a13
            R_ws = rws_p * a13
            R_so = rso_p * a13
            R_wso = rwso_p * a13

            if (av_p < 1.d-6) av_p = 0.65d0
            if (awv_p < 1.d-6) awv_p = 0.65d0
            if (avs_p < 1.d-6) avs_p = 0.65d0
            if (aws_p < 1.d-6) aws_p = 0.65d0
            if (aso_p < 1.d-6) aso_p = 0.65d0
            if (awso_p < 1.d-6) awso_p = 0.65d0

            f_v = 1.d0 / (1.d0 + exp((r - R_v)/av_p))
            f_wv = 1.d0 / (1.d0 + exp((r - R_wv)/awv_p))

            f_s = 1.d0 / (1.d0 + exp((r - R_s)/avs_p))
            df_s = -exp((r - R_s)/avs_p) / avs_p /
     &             (1.d0 + exp((r - R_s)/avs_p))**2

            f_ws = 1.d0 / (1.d0 + exp((r - R_ws)/aws_p))
            df_ws = -exp((r - R_ws)/aws_p) / aws_p /
     &              (1.d0 + exp((r - R_ws)/aws_p))**2

            if (r > 1.d-6) then
                df_so = -exp((r - R_so)/aso_p) / aso_p /
     &                  (1.d0 + exp((r - R_so)/aso_p))**2 / r

                df_wso = -exp((r - R_wso)/awso_p) / awso_p /
     &                   (1.d0 + exp((r - R_wso)/awso_p))**2 / r
            else
                df_so = 0.d0
                df_wso = 0.d0
            endif

            V_central = -Vv_p * f_v - iu * Wv_p * f_wv

            V_surface = 4.d0 * Vs_p * avs_p * df_s +
     &                  4.d0 * iu * Ws_p * aws_p * df_ws

            ls_factor = (J_ch*(J_ch+1.d0) - dble(l*(l+1)) -
     &                   S_ch*(S_ch+1.d0)) / 2.d0

            V_so_term = 2.d0 * Vso_p * df_so * ls_factor +
     &                  2.d0 * iu * Wso_p * df_wso * ls_factor

            V_nuc = V_central + V_surface + V_so_term

            end subroutine

c-----------------------------------------------------------------------
            subroutine eval_coulomb_potential5(r, z12, rc, V_coul)
c           Evaluate finite-size Coulomb potential at radius r
c-----------------------------------------------------------------------
            implicit none
            real*8, intent(in) :: r, z12, rc
            complex*16, intent(out) :: V_coul

            if (abs(z12) < 1.d-10) then
                V_coul = cmplx(0.d0, 0.d0, kind=8)
                return
            endif

            if (rc < 1.d-6) then
                V_coul = e2 * z12 / r
            else if (r >= rc) then
                V_coul = e2 * z12 / r
            else
                V_coul = e2 * z12 / (2.d0 * rc) *
     &                   (3.d0 - (r/rc)**2)
            endif

            end subroutine

        end module scatt_method5
