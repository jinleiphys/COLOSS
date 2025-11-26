        module scatt_method3
c-----------------------------------------------------------------------
c       Module for Method 3: Scattering using Lagrange-Legendre mesh
c
c       This module solves the scattering problem WITHOUT complex scaling
c       using Lagrange-Legendre basis functions on a finite interval [0, R].
c
c       KEY INSIGHT: Lagrange-Laguerre basis has built-in exp(-r/2h) decay
c       which makes it IMPOSSIBLE to represent oscillating scattering waves.
c       Lagrange-Legendre basis is PURE POLYNOMIAL with no intrinsic
c       asymptotic behavior, allowing proper boundary condition matching.
c
c       Physical decomposition: psi_l = F_l(kr) + phi_l
c       where F_l = regular Coulomb/spherical Bessel function
c       phi_l satisfies outgoing wave boundary condition at r = R
c
c       Boundary conditions:
c         - Origin (r=0): phi_l(0) = 0
c         - Outer (r=R): phi'_l(R) = gamma_s * phi_l(R)
c           where gamma_s = H'^+_l(kR) / H^+_l(kR)
c
c       We use standard Lagrange interpolation with Gauss-Legendre points.
c       The differential equation is discretized using collocation method.
c
c       Reference: D. Baye, Physics Reports 565 (2015) 1-107
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

            ! Lagrange-Legendre mesh arrays (shifted Legendre on (0,1) scaled to (0,R))
            real*8, dimension(:), allocatable :: leg_x      ! Mesh points in (0,1)
            real*8, dimension(:), allocatable :: leg_r      ! Mesh points in (0,R)
            real*8, dimension(:), allocatable :: leg_w      ! Integration weights

            ! Derivative matrices from Baye (2015) exact formulas
            real*8, dimension(:,:), allocatable :: D1_mat   ! First derivative d/dr
            real*8, dimension(:,:), allocatable :: D2_mat   ! Kinetic -d²/dr²

            contains

c-----------------------------------------------------------------------
            subroutine init_legendre_mesh_shifted(N, R)
c           Initialize shifted Legendre mesh on (0, R)
c           Following Baye (2015) Section 3.4.5
c
c           Mesh points x_i satisfy P_N(2*x_i - 1) = 0
c           where P_N is the Legendre polynomial
c
c           These are Gauss-Legendre points mapped from (-1,1) to (0,1)
c           then scaled to (0, R)
c
c           INPUT:
c               N: Number of mesh points
c               R: Outer boundary radius
c-----------------------------------------------------------------------
            implicit none
            integer, intent(in) :: N
            real*8, intent(in) :: R

            real*8, dimension(N) :: t_gauss, w_gauss
            integer :: i

            if (allocated(leg_r)) deallocate(leg_r)
            if (allocated(leg_w)) deallocate(leg_w)
            if (allocated(leg_x)) deallocate(leg_x)
            allocate(leg_r(N), leg_w(N), leg_x(N))

            ! Get Gauss-Legendre points on (-1, 1)
            call gauleg(N, -1.d0, 1.d0, t_gauss, w_gauss)

            ! Map from (-1, 1) to (0, 1): x_i = (t_i + 1)/2
            do i = 1, N
                leg_x(i) = (t_gauss(i) + 1.d0) / 2.d0
                ! Gauss weight on (0,1): lambda_i/2
                leg_w(i) = w_gauss(i) / 2.d0
            end do

            ! Scale to (0, R): r_i = R * x_i
            do i = 1, N
                leg_r(i) = R * leg_x(i)
            end do

c            write(*,*) '-------------- Lagrange-Legendre Mesh --------------'
c            write(*,'(A,I4)') ' Number of mesh points: ', N
c            write(*,'(A,F10.4,A)') ' Outer boundary R = ', R, ' fm'
c            write(*,'(A,F10.6)') ' First mesh point r_1 = ', leg_r(1)
c            write(*,'(A,F10.6)') ' Last mesh point r_N = ', leg_r(N)
c            write(*,*) ''

            end subroutine

c-----------------------------------------------------------------------
            subroutine generate_derivative_matrices_baye(N, R)
c           Generate first and second derivative matrices for
c           standard Lagrange-Legendre collocation method
c
c           Using barycentric formula for first derivative:
c           D1_{i≠j} = w_j / (w_i * (r_i - r_j))
c           D1_{ii} = -sum_{k≠i} D1_{ik}  (row sum = 0)
c
c           Second derivative: D2 = D1 * D1
c
c           This gives d²/dr² matrix for the standard wave function,
c           which we use with regularized u = r*ψ formulation.
c
c           INPUT:
c               N: Number of mesh points
c               R: Outer boundary
c-----------------------------------------------------------------------
            implicit none
            integer, intent(in) :: N
            real*8, intent(in) :: R

            integer :: i, j, kk
            real*8 :: bary_w(N)
            real*8 :: prod

            if (allocated(D1_mat)) deallocate(D1_mat)
            if (allocated(D2_mat)) deallocate(D2_mat)
            allocate(D1_mat(N, N), D2_mat(N, N))

            ! Compute barycentric weights: w_i = 1 / prod_{k≠i} (r_i - r_k)
            do i = 1, N
                prod = 1.d0
                do kk = 1, N
                    if (kk /= i) then
                        prod = prod * (leg_r(i) - leg_r(kk))
                    endif
                end do
                bary_w(i) = 1.d0 / prod
            end do

            ! First derivative matrix using barycentric formula
            D1_mat = 0.d0
            do i = 1, N
                do j = 1, N
                    if (i /= j) then
                        D1_mat(i,j) = bary_w(j) / (bary_w(i) *
     &                               (leg_r(i) - leg_r(j)))
                    endif
                end do
                ! Diagonal: sum of row = 0 for derivative of constant
                D1_mat(i,i) = 0.d0
                do kk = 1, N
                    if (kk /= i) then
                        D1_mat(i,i) = D1_mat(i,i) - D1_mat(i,kk)
                    endif
                end do
            end do

            ! Second derivative: D2 = D1 * D1
            D2_mat = matmul(D1_mat, D1_mat)

            end subroutine

c-----------------------------------------------------------------------
            subroutine compute_gamma_s(k_wave, rN, l, gamma_s, eta_som)
c           Compute the logarithmic derivative of outgoing Hankel function
c           gamma_s = H_l^+'(k*r_N) / H_l^+(k*r_N)
c
c           Uses COUL90 to compute F and G, then H^+ = G + i*F
c
c           INPUT:
c               k_wave: wave number
c               rN: outer boundary radius
c               l: angular momentum
c               eta_som: Sommerfeld parameter
c           OUTPUT:
c               gamma_s: complex logarithmic derivative
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

            ! Compute Coulomb wave functions at r_N using COUL90
            x = k_wave * rN
            xlmin = 0.d0
            kfn = 0  ! 0 = Coulomb functions

            call COUL90(x, eta_som, xlmin, l, fc_loc, gc_loc,
     &                  fcp_loc, gcp_loc, kfn, ifail)
            if (ifail /= 0) then
                write(*,*) 'compute_gamma_s: COUL90 ifail=', ifail
            endif

            ! H^+ = G + i*F (outgoing Hankel)
            hhat = gc_loc(l) + iu * fc_loc(l)
            ! dH^+/d(kr) is returned by COUL90
            ! dH^+/dr = k * dH^+/d(kr)
            hhatp = k_wave * (gcp_loc(l) + iu * fcp_loc(l))

            if (abs(hhat) < 1.d-20) then
                write(*,*) "Warning: H^+ very small in compute_gamma_s"
                gamma_s = cmplx(0.d0, k_wave, kind=8)
            else
                gamma_s = hhatp / hhat
            endif

            end subroutine

c-----------------------------------------------------------------------
            subroutine solve_scatt_method3(ich, para)
c           Solve scattering problem using Lagrange-Legendre collocation
c
c           Equation: [d²/dr² - l(l+1)/r² - 2μV/ℏ² + k²] φ = 2μV/ℏ² F_l
c
c           Discretized on Gauss-Legendre points r_1, ..., r_N with:
c           - Rows 1 to N: collocation equations
c           - Boundary conditions incorporated via penalty or projection
c
c           Actually we solve: φ satisfies outgoing BC at R
c           and φ → 0 as r → 0 (automatically satisfied for φ = V*F*Green)
c
c           INPUT:
c               ich: channel index
c               para: potential parameters
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

            complex*16 :: vmod
            complex*16 :: f_l, smat
            real*8 :: reac_xsec

            real*8 :: a13, rrc, z12
            real*8 :: R_outer

            ! Variables for boundary conditions
            complex*16 :: gamma_s
            real*8 :: eta_real

            ! Variables for Coulomb functions
            real*8 :: x_mesh, xlmin_mesh
            real*8, dimension(:), allocatable :: fc_loc, gc_loc
            real*8, dimension(:), allocatable :: fcp_loc, gcp_loc
            integer :: kfn_mesh, ifail_mesh

            ! Variables for amplitude extraction
            complex*16 :: phi_R, hhat_R
            real*8 :: coeff_kin

            ! Variables for boundary condition implementation
            real*8 :: D1_at_R(nr)  ! D1 evaluated at r=R
            real*8 :: phi_at_0     ! φ(0) from interpolation
            real*8 :: L_j_at_R, L_j_at_0
            integer :: i
            real*8 :: prod

            ! Get channel quantum numbers
            l = channel_index%L(ich)
            S = channel_index%S(ich)
            J = channel_index%J(ich)

            ! Use Rmax from mesh module as outer boundary
            R_outer = Rmax
            N_leg = nr

            ! Allocate local Coulomb function arrays
            allocate(fc_loc(0:lmax), gc_loc(0:lmax))
            allocate(fcp_loc(0:lmax), gcp_loc(0:lmax))

            ! Initialize shifted Legendre mesh on (0, R) following Baye (2015)
            call init_legendre_mesh_shifted(N_leg, R_outer)

            ! Generate derivative matrices using exact Baye formulas
            call generate_derivative_matrices_baye(N_leg, R_outer)

            ! Calculate potential parameters
            a13 = para%a2**(1./3.) + para%a1**(1./3.)
            rrc = a13 * para%rc
            z12 = zt * zp

            ! Compute gamma_s (outgoing wave log-derivative at R)
            eta_real = real(eta)
            call compute_gamma_s(k, R_outer, l, gamma_s, eta_real)

c            write(*,*) "DEBUG Method 3 (Legendre):"
c            write(*,*) "  l =", l, " k =", k
c            write(*,*) "  R_outer =", R_outer
c            write(*,*) "  gamma_s =", gamma_s
c            write(*,*) "  eta =", eta_real

            ! Allocate matrices
            allocate(M_matrix(N_leg, N_leg))
            allocate(b_vec(N_leg), c_vec(N_leg))

            M_matrix = cmplx(0.d0, 0.d0, kind=8)
            b_vec = cmplx(0.d0, 0.d0, kind=8)

            ! Kinetic energy coefficient: hbar²/(2*mu) in MeV.fm²
            coeff_kin = hbarc**2 / (2.d0 * mu)

            ! ============================================================
            ! COMPUTE DERIVATIVE AT BOUNDARY r=R
            ! ============================================================
            ! We need dφ/dr at r=R using φ values at mesh points
            ! dφ/dr|_{r=R} = sum_j D1_at_R(j) * φ(r_j)
            !
            ! Compute D1_at_R: derivative of L_j(r) evaluated at r=R
            ! D1_at_R(j) = L'_j(R) = sum over formula

            do jr = 1, N_leg
                ! L'_j(R) = d/dr [prod_{k≠j} (r-r_k)/(r_j-r_k)] at r=R
                D1_at_R(jr) = 0.d0
                do i = 1, N_leg
                    if (i /= jr) then
                        prod = 1.d0
                        do ir = 1, N_leg
                            if (ir /= jr .and. ir /= i) then
                                prod = prod * (R_outer - leg_r(ir)) /
     &                                        (leg_r(jr) - leg_r(ir))
                            endif
                        end do
                        D1_at_R(jr) = D1_at_R(jr) +
     &                    prod / (leg_r(jr) - leg_r(i))
                    endif
                end do
            end do

            ! ============================================================
            ! BUILD MATRIX AND SOURCE VECTOR
            ! ============================================================
            ! Equation for scattered wave φ (NOT regularized):
            ! [d²/dr² - l(l+1)/r² - U(r) + k²] φ = U(r) * F_l(kr)
            ! where U(r) = 2μV_short/ℏ² is the reduced potential
            !       F_l = regular Coulomb/Bessel function
            !
            ! D2_mat = D1*D1 represents d²/dr²
            !
            ! Matrix form: M * c = b
            ! M_ij = D2_ij + [k² - l(l+1)/r_i² - U(r_i)] * δ_ij
            ! b_i = U(r_i) * F_l(k*r_i)

            do ir = 1, N_leg - 1  ! Rows 1 to N-1: differential equation
                ! Compute Coulomb functions at this mesh point
                x_mesh = k * leg_r(ir)
                xlmin_mesh = 0.d0
                kfn_mesh = 0

                call COUL90(x_mesh, eta_real, xlmin_mesh, lmax,
     &                      fc_loc, gc_loc, fcp_loc, gcp_loc,
     &                      kfn_mesh, ifail_mesh)

                ! Compute short-range potential at this point
                call compute_potential_at_r(leg_r(ir), para, ich, vmod)

                ! Scale potential: U = 2μV/ℏ² (units of 1/fm²)
                vmod = vmod / coeff_kin

                ! Fill matrix row: M = D2 + (k² - l(l+1)/r² - U)*δ
                do jr = 1, N_leg
                    ! Second derivative: d²/dr²
                    M_matrix(ir, jr) = D2_mat(ir, jr)

                    ! Diagonal terms
                    if (ir == jr) then
                        ! Energy term: +k²
                        M_matrix(ir, ir) = M_matrix(ir, ir) + k*k

                        ! Centrifugal: -l(l+1)/r²
                        M_matrix(ir, ir) = M_matrix(ir, ir) -
     &                      dble(l*(l+1)) / leg_r(ir)**2

                        ! Short-range potential: -U
                        M_matrix(ir, ir) = M_matrix(ir, ir) - vmod
                    endif
                end do

                ! Source term: U * F_l(k*r)
                b_vec(ir) = vmod * fc_loc(l)
            end do

            ! Row N: Outgoing wave boundary condition at r = R
            ! φ'(R) - γ_s * φ(R) = 0
            ! where γ_s = H'^+_l(kR) / H^+_l(kR)
            !
            ! Using: φ'(R) = sum_j D1_at_R(j) * φ(r_j)
            !        φ(R) = sum_j L_j(R) * φ(r_j)
            !
            ! L_j(R) = prod_{k≠j} (R - r_k) / (r_j - r_k)

            do jr = 1, N_leg
                ! Compute L_j(R)
                L_j_at_R = 1.d0
                do i = 1, N_leg
                    if (i /= jr) then
                        L_j_at_R = L_j_at_R * (R_outer - leg_r(i)) /
     &                                        (leg_r(jr) - leg_r(i))
                    endif
                end do

                ! BC: φ' - γ_s*φ = 0
                M_matrix(N_leg, jr) = D1_at_R(jr) - gamma_s * L_j_at_R
            end do
            b_vec(N_leg) = cmplx(0.d0, 0.d0, kind=8)

c            ! ============================================================
c            ! DEBUG: Check matrix condition
c            ! ============================================================
c            write(*,*) "DEBUG: M_matrix diagonal sample:"
c            write(*,*) "  M(1,1) =", M_matrix(1,1)
c            write(*,*) "  M(N/2,N/2) =", M_matrix(N_leg/2, N_leg/2)
c            write(*,*) "  M(N,N) =", M_matrix(N_leg, N_leg)
c            write(*,*) "DEBUG: b_vec sample:"
c            write(*,*) "  b(1) =", b_vec(1)
c            write(*,*) "  b(N/2) =", b_vec(N_leg/2)
c            write(*,*) "  b(N) =", b_vec(N_leg)

            ! ============================================================
            ! SOLVE LINEAR SYSTEM
            ! ============================================================
            c_vec = b_vec
            call z_lineq(N_leg, M_matrix, c_vec)

            ! ============================================================
            ! EXTRACT SCATTERING AMPLITUDE
            ! ============================================================
            ! φ(R) = sum_j L_j(R) * c_j
            ! where c_j are the expansion coefficients (function values at mesh points)
            !
            ! Asymptotically: φ ~ f_l * H^+(kr)
            ! So: f_l = φ(R) / H^+(kR)

            phi_R = cmplx(0.d0, 0.d0, kind=8)
            do jr = 1, N_leg
                L_j_at_R = 1.d0
                do i = 1, N_leg
                    if (i /= jr) then
                        L_j_at_R = L_j_at_R * (R_outer - leg_r(i)) /
     &                                        (leg_r(jr) - leg_r(i))
                    endif
                end do
                phi_R = phi_R + L_j_at_R * c_vec(jr)
            end do

            ! Compute H^+ at boundary
            x_mesh = k * R_outer
            call COUL90(x_mesh, eta_real, 0.d0, l,
     &                  fc_loc, gc_loc, fcp_loc, gcp_loc, 0, ifail_mesh)
            hhat_R = gc_loc(l) + iu * fc_loc(l)

            ! f_l = φ(R) / H^+(kR)
            if (abs(hhat_R) < 1.d-20) then
                write(*,*) "Warning: H^+ very small at boundary"
                f_l = cmplx(0.d0, 0.d0, kind=8)
            else
                f_l = phi_R / hhat_R
            endif

c            write(*,*) "DEBUG: phi(R) =", phi_R
c            write(*,*) "DEBUG: H^+(kR) =", hhat_R
c            write(*,*) "DEBUG: f_l =", f_l

            ! Store scattering amplitude
            scatt_amp_nuc_channel(ich) = f_l

            ! Compute S-matrix: S = 1 + 2*i*f_l
            smat = 1.d0 + 2.d0 * iu * f_l

            ! Compute reaction cross section
            reac_xsec = pi/k/k/(2.d0*S+1.d0)*(2.d0*J+1.d0)
     &                 *(1.d0 - abs(smat)**2) * 10.d0  ! fm² to mb

            ! Output results
            write(*, 300) l, S, J, real(smat), aimag(smat), reac_xsec
300         FORMAT(I3,3x,F3.1,2x,F5.1,' |  (',F10.6,', ',F10.6,')  | ',
     &             F14.4)

            write(60, 101) real(smat), aimag(smat), l, S, J
            write(61, 101) real(f_l), aimag(f_l), l, S, J
101         FORMAT(F10.6,2x,F10.6,"  (L S J):",I3,3x,F3.1,2x,F5.1)

            ! Cleanup
            deallocate(M_matrix, b_vec, c_vec)
            deallocate(fc_loc, gc_loc, fcp_loc, gcp_loc)

            end subroutine

c-----------------------------------------------------------------------
            subroutine compute_potential_at_r(r, para, ich, V_short)
c           Compute the short-range potential at radius r
c           V_short = V_nuclear + V_coulomb_finite - Z1*Z2*e²/r
c
c           INPUT:
c               r: radius (fm)
c               para: potential parameters
c               ich: channel index
c           OUTPUT:
c               V_short: short-range potential (MeV)
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

            ! Nuclear potential (Woods-Saxon form)
            call eval_nuclear_potential(r, para, ich, V_nuc_r)

            ! Finite-size Coulomb potential
            call eval_coulomb_potential(r, z12, rrc, V_coul_r)

            ! Point Coulomb (to be subtracted)
            if (abs(z12) > 1.d-10) then
                V_coul_point = e2 * z12 / r
            else
                V_coul_point = 0.d0
            endif

            ! Short-range potential
            V_short = V_nuc_r + V_coul_r - V_coul_point

            end subroutine

c-----------------------------------------------------------------------
            subroutine eval_nuclear_potential(r, para, ich, V_nuc)
c           Evaluate nuclear optical potential at radius r
c           Includes central, surface, and spin-orbit terms
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

            ! Get potential parameters (using correct pot_para field names)
            Vv_p = para%vv; rv_p = para%rvv; av_p = para%avv
            Wv_p = para%wv; rwv_p = para%rw; awv_p = para%aw
            Vs_p = para%vs; rvs_p = para%rvs; avs_p = para%avs
            Ws_p = para%ws; rws_p = para%rws; aws_p = para%aws
            Vso_p = para%vsov; rso_p = para%rsov; aso_p = para%asov
            Wso_p = para%vsow; rwso_p = para%rsow; awso_p = para%asow

            ! Compute radii
            R_v = rv_p * a13
            R_wv = rwv_p * a13
            R_s = rvs_p * a13
            R_ws = rws_p * a13
            R_so = rso_p * a13
            R_wso = rwso_p * a13

            ! Avoid division by zero for zero diffuseness
            if (av_p < 1.d-6) av_p = 0.65d0
            if (awv_p < 1.d-6) awv_p = 0.65d0
            if (avs_p < 1.d-6) avs_p = 0.65d0
            if (aws_p < 1.d-6) aws_p = 0.65d0
            if (aso_p < 1.d-6) aso_p = 0.65d0
            if (awso_p < 1.d-6) awso_p = 0.65d0

            ! Woods-Saxon form factors (volume)
            f_v = 1.d0 / (1.d0 + exp((r - R_v)/av_p))
            f_wv = 1.d0 / (1.d0 + exp((r - R_wv)/awv_p))

            ! Surface (derivative) form factors
            f_s = 1.d0 / (1.d0 + exp((r - R_s)/avs_p))
            df_s = -exp((r - R_s)/avs_p) / avs_p /
     &             (1.d0 + exp((r - R_s)/avs_p))**2

            f_ws = 1.d0 / (1.d0 + exp((r - R_ws)/aws_p))
            df_ws = -exp((r - R_ws)/aws_p) / aws_p /
     &              (1.d0 + exp((r - R_ws)/aws_p))**2

            ! Spin-orbit form factors
            if (r > 1.d-6) then
                df_so = -exp((r - R_so)/aso_p) / aso_p /
     &                  (1.d0 + exp((r - R_so)/aso_p))**2 / r

                df_wso = -exp((r - R_wso)/awso_p) / awso_p /
     &                   (1.d0 + exp((r - R_wso)/awso_p))**2 / r
            else
                df_so = 0.d0
                df_wso = 0.d0
            endif

            ! Central potential (volume): -V*f
            V_central = -Vv_p * f_v - iu * Wv_p * f_wv

            ! Surface potential: +4*a*V*df/dr (note sign convention)
            V_surface = 4.d0 * Vs_p * avs_p * df_s +
     &                  4.d0 * iu * Ws_p * aws_p * df_ws

            ! Spin-orbit potential
            ! <l.s> = (J(J+1) - l(l+1) - S(S+1))/2
            ls_factor = (J_ch*(J_ch+1.d0) - dble(l*(l+1)) -
     &                   S_ch*(S_ch+1.d0)) / 2.d0

            V_so_term = 2.d0 * Vso_p * df_so * ls_factor +
     &                  2.d0 * iu * Wso_p * df_wso * ls_factor

            ! Total nuclear potential
            V_nuc = V_central + V_surface + V_so_term

            end subroutine

c-----------------------------------------------------------------------
            subroutine eval_coulomb_potential(r, z12, rc, V_coul)
c           Evaluate finite-size Coulomb potential at radius r
c           Uses uniform sphere model
c-----------------------------------------------------------------------
            implicit none
            real*8, intent(in) :: r, z12, rc
            complex*16, intent(out) :: V_coul

            if (abs(z12) < 1.d-10) then
                V_coul = cmplx(0.d0, 0.d0, kind=8)
                return
            endif

            if (rc < 1.d-6) then
                ! Point Coulomb if no finite size
                V_coul = e2 * z12 / r
            else if (r >= rc) then
                V_coul = e2 * z12 / r
            else
                ! Uniform sphere: V = (e²Z1Z2/2Rc) * (3 - (r/Rc)²)
                V_coul = e2 * z12 / (2.d0 * rc) *
     &                   (3.d0 - (r/rc)**2)
            endif

            end subroutine

        end module scatt_method3
