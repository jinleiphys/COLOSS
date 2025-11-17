        module mesh
c-------Some general parameters-------------------
                real*8 :: scaling_factor! its value should be less than 1!
                real*8 :: rmax! value of the biggest zero
                real*8 :: ctheta!rotation angle in degree
                real*8 :: theta_rad!rotation angle in rad
                complex*16 :: eitheta!value for exp(i \theta)
                real*8 :: alpha!parameter for the laguerre polynomial
                integer :: nr!number of grid points in gauss quadrature method
                real*8, dimension(:),allocatable :: laguerre_rr, laguerre_rw!original weight and mesh points for gauss quadrature
                real*8, dimension(:),allocatable :: mesh_rr, mesh_rw!transformed weight and mesh points for gauss quadrature
                
                logical :: backrot

                integer :: cwftype!type for complex coulomb functions

                ! parameters for guassian quadrature when evaluating matrix elements
                integer :: numgauss 
                real*8 :: rmaxgauss
                real*8, dimension(:),allocatable :: gauss_rr, gauss_rw
                logical :: matgauss
                logical :: bgauss
                
                ! parameters for cross-section angular distribution
                real*8 :: thetah, thetamax
                integer :: theta_n_max

                logical :: readinpot
                integer :: method

                logical :: nonlocal
                real*8 :: nlbeta
                
        end module