# COLOSS
**COLOSS** is a program designed to address the scattering problem using a bound-state technique known as complex scaling. In this method, the oscillatory boundary conditions of the wave function are transformed into exponentially decaying ones, accommodating the long-range Coulomb interaction. This program implements the Woods-Saxon form of a realistic optical potential, with all potential parameters included in a well-designed input format for ease of use. This design offers users straightforward access to compute \(S\)-matrices and cross-sections of the scattering process. Here we use the Lagrange-Laguerre functions as the basis functions, and provide two method to perform the numerical integration. One is based on the approximation of Lagrange functions, and the other is the direct Gauss-Legendre quadrature. Additionally, COLOSS incorporates two distinct rotation methods, making it adaptable to potentials without analytical expressions. 

## Input description
COLOSS utilizes the FORTRAN namelist to construct a user-friendly input format. 

### General namelist: nr, alpha, Rmax, ctheta, matguass, bguass, numgauss, rmaxgauss, method,backrot, cwftype, thetah, thetamax
- nr(integer*4): number of the Lagrange-Laguerre basis/ order of the generalized Laguerre polynomial in the calculation
- alpha(real*8): parameter of the generalized Laguerre polynomial
- Rmax(real*8): maximum value of the points in the scaled Lagrange-Laguerre mesh
- ctheta(real*8): rotation angle for complex scaling in degrees
- matgauss(logical): determines whether to use the Gauss-Legendre quadrature to evaluate the matrix elements of the Hamiltonian
- bgauss(logical): determines whether to use the Gauss-Legendre quadrature to evaluate the inhomogeneous terms in the linear equation
- numgauss(integer*4): number of Gauss-Legendre mesh points
- rmaxgauss(real*8): maximum radius of the Gauss-Legendre mesh points
- method(integer*4): option for 2 different method to calculate the scattering amplitude, set it as 1 for linear equation method, and 2 for the Green's function method
- backrot(logical): determines whether to rotate the basis backward or to rotate the potential directly in the calculation
- cwftype(integer*4): options for 2 different subroutines to calculate Coulomb wave functions. Set it as 1 to call COULCC, and 2 to call cwfcomplex
- thetah(real*8): step size for the angle in the output differential cross section
- thetamax(real*8): maximum angle in the output differential cross section

### System namelist: zp, zt, massp, masst, namep, namet, lmin, lmax, elab
- zp, zt(real*8): charge number of the projectile and the target
- massp, masst(real*8): mass number of the projectile and the traget
- namep, namet(character): name of the projectile and the target
- lmin, lmax(integer*4): minimum/maximum total angular momentum of the reaction system considered in the calculation
- sp(real*8): spin of the projectile
- elab(real*8): incident kinetic energy of the projectile in the laboratory frame

### Pot namelist: vv, rv, av, wv, rw, aw, vs, rvs, avs, ws, rws, aws, vsov, rsov, asov, vsow, rsow, asow, rc
- vv, rv, av(real*8): depth, radius, and width parameters of the real volume term in OMP
- wv, rw, aw(real*8): depth, radius, and width parameters of the imaginary volume term in OMP
- vs, rvs, avs(real*8): depth, radius, and width parameters of the real surface term in OMP
- ws, rws, aws(real*8): depth, radius, and width parameters of the imaginary surface term in OMP
- vsov, rsov, asov(real*8): depth, radius, and width parameters of the real spin-orbit coupling term in OMP
- vsow, rsow, asow(real*8): depth, radius, and width parameters of the imaginary spin-orbit coupling term in OMP
- rc(real*8): charge radius for Coulomb interaction in OMP

### nonlocalpot namelist: nonlocal, nlbeta
- nonlocal(logical): determines whether to use a nonlocal form of OMP proposed by F. PEREY and B. BUCK(Nuclear Physics,32,353-380)
- nlbeta(real*8): beta parameter in F. PEREY and B. BUCK's nonlocal form of OMP

## Output description
- fort.1: The local copy of the input file
- fort.2: The table of all the angular momentum channles consided in the calculation.
- fort.10: The scaled Lagrange-Laguerre mesh points and weights
- fort.60: The S-matrices for different angular momentum channels
- fort.61: The nuclear scattering amplitudes for different angular momentum channels
- fort.67: The angular distribution of the differential cross section

## Getting started
In **COLOSS**, we use GNU Compiler Collection(GCC). Make sure you have installed GCC. For linear algebra subroutines, we use LAPACK package.
To complie the program:
- Edit the Makefile and specify your path of the LAPACK package on your machine.
- Use the provided Makefile to compile all the source code.
`prompt> make`
- Transfer the executable program, **COLOSS**, to the test directory, and initiate program execution through standard input:
`prompt> ./COLOSS < inpufile`

## Examples for 40Ca(n,n)40Ca at 20MeV
### Input
```
&general  
        nr=60  alpha=0 Rmax=40 ctheta=7 
        matgauss=f bgauss=f method=1
        thetah=1.0 thetamax=180 /

&system 
        zp=0    massp=1   namep='n'
        zt=20   masst=40  namet='40Ca'
        jmin=0 jmax=10  elab=20 sp=0.5/  

&pot
        vv=46.553   rv=1.185    av=0.672
        wv=1.777    rw=1.185    aw=0.672
        vs=0        rvs=0       avs=0
        ws=7.182    rws=1.288   aws=0.538
        vsov=5.343  rsov=0.996  asov=0.590
        vsow=-0.110 rsow=0.996  asow=0.590
        rc=1.698 /

&nonlocalpot
        nonlocal=f nlbeta=0.0 /

```

### Output
```
 ************************************************************
 *                  COLOSS: Complex-scaled                  *
 *          Optical and couLOmb Scattering Solver           *
 ************************************************************
 
 Rotation Angle for Complex Scaling:  7.00 degrees
 Rotation Operation: Applied to Potential
 
 -------------------Reaction system-------------------
   Projectile:  n     (A:     1.00 Z:    0.00)
       Target:  40Ca  (A:    40.00 Z:   20.00)
 Lab  Frame Energy:   20.00 MeV
 C.M. Frame Energy:   19.51 MeV
 
 Generating the channel index
 There are   20 different channels 
 Total J Range:  0.0 <= J <=  10.0

 -------------- Lagrange-Laguerre Mesh --------------
 Laguerre Polynomial Order:  60
 Laguerre Mesh Max Value: Scaled from   219.3181 to    40.0000
 Scaling Factor:    0.1824 (fm)
 
 ------------- Optical Model Potential -------------
 Optical Potential for A =  40.0 Z = 20.0 at  20.000 MeV rc =  1.698
  Vv     rvv    avv    Wv     rw     aw
 46.553  1.185  0.672  1.777  1.185  0.672

  Vs     rvs    avs    Ws     rws    aws
  0.000  0.000  0.000  7.182  1.288  0.538

  Vsov   rsov   asov   vsow   rsow   asow
  5.343  0.996  0.590 -0.110  0.996  0.590

 ------------------------------------------------
 Initializing Complex Coulomb Function
 Sommerfeld Parameter:    0.0000
 Generating Rotated Coulomb Wave Function on Laguerre Mesh
                    (CPU  time =  0.00025700  seconds)
 
 Using Linear Equation Method to Solve
  L     S      J |  (Real(Smat), Imag(Smat))  |  Partial Wave Reac_Xsec (mb)
 ----------------------------------------------------------------------------
  0   0.5    0.5 |  (  0.349804,   0.223235)  |        28.5534
  1   0.5    0.5 |  (  0.508537,   0.041797)  |        25.5125
  1   0.5    1.5 |  (  0.479616,   0.195699)  |        50.4749
  2   0.5    1.5 |  (  0.321689,  -0.152658)  |        60.2393
  2   0.5    2.5 |  (  0.392012,   0.086108)  |        86.8096
  3   0.5    2.5 |  ( -0.032075,  -0.496351)  |        77.8789
  3   0.5    3.5 |  (  0.317075,  -0.310480)  |       110.8004
  4   0.5    3.5 |  ( -0.134338,  -0.189716)  |       130.5159
  4   0.5    4.5 |  (  0.171986,  -0.508775)  |       122.7206
  5   0.5    4.5 |  (  0.448912,   0.217764)  |       129.5308
  5   0.5    5.5 |  (  0.219593,   0.066785)  |       196.0549
  6   0.5    5.5 |  (  0.878907,   0.126747)  |        43.7629
  6   0.5    6.5 |  (  0.839528,   0.167094)  |        64.5332
  7   0.5    6.5 |  (  0.977509,   0.037757)  |        10.3944
  7   0.5    7.5 |  (  0.975107,   0.047042)  |        12.9563
  8   0.5    7.5 |  (  0.995748,   0.009952)  |         2.3144
  8   0.5    8.5 |  (  0.995614,   0.011765)  |         2.6742
  9   0.5    8.5 |  (  0.999183,   0.002546)  |         0.5053
  9   0.5    9.5 |  (  0.999178,   0.002913)  |         0.5636
 10   0.5    9.5 |  (  0.999841,   0.000644)  |         0.1093
 ----------------------------------------------------------------------------
                    (CPU  time =    0.003733  seconds)
 
----------------- Files Created -----------------
  1: local copy of input.                         
  2: list of channel index.                       
 10: scaled Laguerre mesh points and weights.     
 60: S matrix LSJ distribution.                   
 61: scat amplitude LSJ distribution.             
 67: cross section angular distribution.          
          Total CPU  time =    0.029168  seconds
  
 ***************** CONSTANTS **********************
 * hbarc=197.32697 MeV.fm     e^2= 1.43997 MeV.fm *
 * so, alpha= 1/137.03547       amu=931.4943 MeV  *
 **************************************************
 
```
