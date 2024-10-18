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

### System namelist: zp, zt, massp, masst, namep, namet, jmin, jmax, elab
- zp, zt(real*8): charge number of the projectile and the target
- massp, masst(real*8): mass number of the projectile and the traget
- namep, namet(character): name of the projectile and the target
- jmin, jmax(integer*4): minimum/maximum total angular momentum of the reaction system considered in the calculation
- elab(real*8): incident kinetic energy of the projectile in the laboratory frame

### Pot namelist: vv, rv, av, wv, rw, aw, vs, rvs, avs, ws, rws, aws, vsov, rsov, asov, vsow, rsow, asow, rc
- vv, rv, av(real*8): depth, radius, and width parameters of the real volume term in OMP
- wv, rw, aw(real*8): depth, radius, and width parameters of the imaginary volume term in OMP
- vs, rvs, avs(real*8): depth, radius, and width parameters of the real surface term in OMP
- ws, rws, aws(real*8): depth, radius, and width parameters of the imaginary surface term in OMP
- vsov, rsov, asov(real*8): depth, radius, and width parameters of the real spin-orbit coupling term in OMP
- vsow, rsow, asow(real*8): depth, radius, and width parameters of the imaginary spin-orbit coupling term in OMP
- rc(real*8): charge radius for Coulomb interaction in OMP

### nonlocal_pot namelist: nonlocal, nonlocal_beta
- nonlocal(logical): determines whether to use a nonlocal form of OMP given in F. PEREY and B. BUCK(Nuclear Physics,32,353-380)
- nonlocal_beta(real*8): beta parameter in F. PEREY and B. BUCK's nonlocal form of OMP

## Output description
- fort.1: the local copy of the input file
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

## Examples for 93Nb(d,d)93Nb at 20MeV
### Input
```
&general  
    nr=60  alpha=0 Rmax=40 ctheta=6 
    matgauss=f bgauss=f method=1
    thetah=1.0 thetamax=180 /

&system 
    zp=1    massp=2   namep='2H'
    zt=41   masst=93  namet='93Nb'
    jmin=0 jmax=20  elab=20   /  

&pot 
    vv=84.323 rv=1.174 av=0.809
    wv=0.351 rw=1.563 aw=0.904
    vs=0 rvs=0 avs=0
    ws=14.247 rws=1.328 aws=0.669 
    rc=1.698 /

```

### Output
```
 ************************************************************
 *                  COLOSS: Complex-scaled                  *
 *          Optical and couLOmb Scattering Solver           *
 ************************************************************
 
 Rotation Angle for Complex Scaling:  6.00 degrees
 Rotation Operation: Applied to Potential
 
 -------------------Reaction system-------------------
   Projectile:  2H    (A:     2.00 Z:    1.00)
       Target:  93Nb  (A:    93.00 Z:   41.00)
 Lab  Frame Energy:   20.00 MeV
 C.M. Frame Energy:   19.58 MeV
 Total J Range:  0 <= J <=  20
 
 -------------- Lagrange-Laguerre Mesh --------------
 Laguerre Polynomial Order:  60
 Laguerre Mesh Max Value: Scaled from   219.3181 to    40.0000
 Scaling Factor:    0.1824 (fm)
 
 ------------- Optical Model Potential -------------
 Optical Potential for A =  93.0 Z = 41.0 at  20.000 MeV rc =  1.698
  Vv     rvv    avv    Wv     rw     aw
 84.323  1.174  0.809  0.351  1.563  0.904

  Vs     rvs    avs    Ws     rws    aws
  0.000  0.000  0.000 14.247  1.328  0.669

 ------------------------------------------------
 Initializing Complex Coulomb Function
 Sommerfeld Parameter:    2.0419
 Generating Rotated Coulomb Wave Function on Laguerre Mesh
                    (CPU  time =  0.00051300  seconds)
 
 Using Linear Equation Method to Solve
 l  |   S-matrix (real, imag)   | Partial Wave Reaction Cross Section (mb)
 -----------------------------------------------------------------------
  0 | ( -0.105887,  -0.084723)  |     1.6814
  1 | (  0.060705,   0.092842)  |     5.0755
  2 | (  0.118453,  -0.060172)  |     8.4134
  3 | ( -0.084936,  -0.028826)  |    11.8939
  4 | ( -0.120555,   0.095725)  |    15.0509
  5 | (  0.013458,  -0.016254)  |    18.8336
  6 | (  0.041922,  -0.194055)  |    21.3901
  7 | (  0.117786,  -0.102841)  |    25.0654
  8 | (  0.385315,  -0.007029)  |    24.7947
  9 | (  0.569733,   0.074756)  |    21.7993
 10 | (  0.800504,   0.118034)  |    12.4194
 11 | (  0.926720,   0.072771)  |     5.3538
 12 | (  0.975106,   0.034918)  |     2.0533
 13 | (  0.991320,   0.015328)  |     0.7885
 14 | (  0.996870,   0.006580)  |     0.3083
 15 | (  0.998845,   0.002806)  |     0.1222
 16 | (  0.999566,   0.001194)  |     0.0490
 17 | (  0.999835,   0.000507)  |     0.0198
 18 | (  0.999936,   0.000215)  |     0.0081
 19 | (  0.999975,   0.000091)  |     0.0033
 20 | (  0.999990,   0.000039)  |     0.0014
 -----------------------------------------------------
                    (CPU  time =    0.004426  seconds)
 
----------------- Files Created -----------------
  1: local copy of input.                         
 10: scaled Laguerre mesh points and weights.     
 60: cross section j distribution.                
 61: scat amplitude j distribution.               
 67: cross section angular distribution.          
  
 ***************** CONSTANTS **********************
 * hbarc=197.32697 MeV.fm     e^2= 1.43997 MeV.fm *
 * so, alpha= 1/137.03547       amu=931.4943 MeV  *
 **************************************************
```
