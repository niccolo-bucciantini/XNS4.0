# Files and outputs

Here is a list of all the files and subroutines included int the XNS package, together with a brief
description of what they do and how they operates. We start with the Fortran 90 files of the code (with
extension *.f90*), then we describe the output files produced by a run (with extension *.dat*), and we conclude
with the Python 3 files needed for visualisation (with extension *.py* or *.ipynb*, for use in Jupyter). The code
must run in double precision for convergence. The *makefile* is provided for the gfortran (GNU) compiler. For more
info on the compilation of XNS, see section "Compiling XNS".

## Files

- **XNS.f90** - main program. Makes some consistency checks, and invokes XNSMAIN. Depending on
the pre-compiling option it simply calls XNSMAIN (if make serial is used to compile the code), or it
performs a Newton-Raphson search for an equilibrium model with a given value for a desired quantity
of interest, i.e. a certain value of the central density or gravitational mass (if make nwtrps is used). Moreover, it can compute many models in parallel with different initial conditions (if make parspace is used).
<br><br>
- **XNSMAIN.f90**
    - subroutine **XNSMAIN** - the main kernel of the code: it defines the grid, builds a 2D initial guess based on the 1D TOV output of TOVINIMOD.f90, performs the convergence loop calling all the various metric solvers and procedures in the appropriate order. When the loop is over, it writes all the outputs.
    - subroutine **CONFORMAL** - solves for the scalar Poisson-like equation for $\psi$.
    - subroutine **LAPSE** - solves for the scalar Poisson-like equation for $\alpha \psi$.
    - subroutine **SHIFTPHI** - solves the $\phi$ component of the two vector Poisson equations for $W^i$ and $\beta ^i$, given the corresponding source terms.
    - subroutine **CURV1** - computes the curvature source term in the routines for $\psi$ and $\alpha \psi$.
    - subroutine **CURV2** - computes the curvature source term in the routine for $\beta ^\phi$.
    - subroutine **SOURCECHI** - computes the metric terms and derivatives used in the source terms of the equation for $\chi$.
    - subroutine **SOLVECHI** - solves for the scalar field equation for $\chi$, given the source term. Note that since $\chi$ also appears in the source term, a single call of SOLVECHI won't yield a true solution.
    - subroutine **CHISOL** - wrapper that calls SOURCECHI and CHISOL in a relaxation loop, until it converges to the true solution.
    - subroutine **CHIDERIVS** - computes the derivatives of the scalar field ($Q^\mu$) and the scalar coupling function $\mathcal{A}(\chi)$.
<br><br>
- **SYSTEMXNS.f90**
    - module **SYSTEMXNS** - contains various parameters of the run, to be specified by the user, and definitions of common arrays (see section "User parameters").
<br><br>
- **HYDROEQ.f90**
    - subroutine **HYDROEQ** - given the CFC metric and a value of $\rho _\mathrm{c}$ it computes the equilibrium
    configuration for the corresponding Bernoulli integral. It finally calls
    hydrovar_<x> depending on the physical parameter set in SYSTEMXNS.f90 to compute local
    equilibrium quantities.
    - subroutine **HYDROVAR**, **HYDROVAR_TOR**, **HYDROVAR_POL** - they compute local equilibrium
    quantities such as $\rho$, $p$, $v^\phi$, $B^i$ and $E_i$ depending on the specific choice for the magnetization (respectively unmagnetized case, purely toroidal magnetic field and poloidal magnetic field).
    - subroutine **COVTERM** - computes the local terms of the metric tensor
    - subroutine **CONS_TO_PRIM** - computes the inversion from conserved to primitive variables.
    - subroutine **CONS_TO_PRIM_POL** - computes the inversion from conserved to primitive variables for the specific case of poloidal field.
    - subroutine **QUANTITIES** - computes several quantities (e.g. mass, energy, angular momentum, scalar charge, magnetic deformation) at the end of the convergence loop, according to standard definitions. See also SDB20 for some definitions in the case of STTs.
    - subroutine **SOURCEPOT** - compute source terms (currents and metric) for the Grad-Shafranov
    Equation or Maxwell equations depending if the rotational rate OMG is set to zero or not.
    - subroutine **VECPOTPHI** - called by the subroutine hydrovar_pol when OMG.EQ.0, it solves
    the Grad-Shafranov Equation.
    - subroutine **MXWLSOL** - called by subroutine hydrovar_pol when OMG.NE.0, it solves iteratively the Maxwell-Ampere and the Maxwell-Gauss equation. It finally corrects the solution for the electric potential $\Phi$ in order to guarantee that the MHD condition $\Phi = -\Omega \Psi + C$ is valid
    inside the star. Indeed, as explained in Pili et al. (2017), the solution for $\Phi$ obtained by solving
    the non-homogeneus Maxwell equations, does not satisfy the perfect conducting relation inside
    the star, but differs from the MHD solution solution $\Phi _\mathrm{MHD} = -\Omega \Psi + C$ by an harmonic function
    $\Phi _\mathrm{a}$ so that $\Phi = \Phi _\mathrm{MHD} + \Phi _\mathrm{a}$ with $\Delta \Phi _\mathrm{a} = 0$. The harmonic function is obtained evoking the laplace subroutine.
    - subroutine **SOLVEAPHI** and subroutine **SOLVEATIME** - solve respectively for the Maxwell-Ampere and Maxwell-Gauss equations.
    - subroutine **LAPLACE** - solves the equations $\Phi _\mathrm{a} \big |_{S _\mathrm{NS}} = \Sigma _l Y(\theta) (a_l r^l) |_{S _\mathrm{NS}}$ (inside the star) and $\Phi _\mathrm{a} \big |_{S _\mathrm{NS}} = \Sigma _l Y(\theta) (b_l r^{-(l+1)}) |_{S _\mathrm{NS}}$ (outside the star), where $S _\mathrm{NS}$ is stellar surface and $\Phi _\mathrm{a} \big |_{S _\mathrm{NS}} = (\Phi + \Omega \Psi + C) |_{S _\mathrm{NS}}$. Each system of equations is solved with a LU decomposition and a subsequent backward substitution adopting the routines provided in the Numerical Recipes (ludcmp and
    lubksb). Notice that, in order to avoid spurious effects, the surface terms are evaluated on top
    of the super-ellipsoid that best fit the numerical surface.
<br><br>
- **ROTATION.f90**
    - subroutine **CHECKROTDIFF** - performs several checks for consistency in the choice of the options and parameters for differential rotation. 
    - subroutine **OMEGAVALUE** - derives the function $\Omega = \Omega (r,\theta)$ for the differential rotation.
    - subroutine **OMEGA3LVALUE** - derives the function $\Omega = \Omega (r,\theta)$ and A3L ($\int jd\Omega$ or $\int \Omega dj$, depeding if omega or J-space, the parameter the enters the Bernoulli Equation for the differential rotation) for the differential rotation.
    - subroutine **FODFO_OS** - evaluate FO and DFO required to compute the local rotation rate, in the Omega-space, for a differential rotator.
    - subroutine **A3L_OS** - evaluate A3L $=\int jd\Omega$  in the Omega-space for a differential rotator. This is the parameter the enters the Bernoulli Equation for the differential rotation.
    - subroutine **PARS_VALUE_JS** evaluate JMAX, XVAL and YVAL (only needed for Uryu Laws, see Franceschetti et al. 2022 A8, A9, A18).
    - subroutine **FODFO_JS** - evaluate FO and DFO in the J-space for a differential rotator (only needed for Uryu Laws, see Franceschetti et al. 2022 A2, A3).
    - subroutine **A3L_JS** evaluate A3L $=\int \Omega dj$ in the J-space for a differential rotator. This is the parameter the enters the Bernoulli Equation for the differential rotation.
<br><br>
- **TOVINIMOD.f90**
    - subroutine **TOVINIMOD** - solves the 1D TOV (either in GR or in STTs) equations in isotropic coordinates to provide the initial guess. It uses a relaxation method to achieve convergence.
    - subroutine **EXPANSION** - a Taylor expansion of the TOV equations at small initial radii (they are
    singular for $r \rightarrow 0$).
    - subroutine **TOVEQS** - provides the derivatives needed to integrate the TOV equations via the
    RK4 method.
    - subroutine **RK4** - the 4th order RK integrator (modified from the Numerical Recipes).
    - subroutine **MASSFIND** - computes the ADM mass of the scalarised TOV solution at a given radius (either at the middle of the grid or at its outer edge) by knowing the scalar charge and the value of $\mu$ at that point. It is used in order to achieve convergence, as explained in SBD20 Appendix B.
<br><br>
- **PHYSICS.f90**
    - subroutine **GRIDBUILD** - computes the radial grid (either uniform or stretched) and derivative terms used in the DGTSV subroutine.
    - subroutine **FUNCD_STRETCH** - used by the root-finding subroutine to derive the stretching factor for the grid, if it is stretched.
    - subroutine **EOSTABLEREAD** - reads the EoS file specified by the FILEEOS parameter. Note that the first line in the file must be the number of points present in the file, and this must be equal to the parameter NPTRHO. The subsequent lines must contain the minimum and maximum density and their indexes (second line), the minimum and maximum pressure and their indexes (third line), the minimum and maximum internal energy and their indexes (fourth line), the minimum and maximum enthalpy and their indexes (fifth line). Then, the table is read. Please, refer to the routine code to see the specific structure that the EoS file must have.
    - subroutine **RHO2EOS** - given $\rho$, it computes the pressure $p$, the internal energy $\varepsilon$ and the enthalpy $h$ according to the tabulated EoS.
    - subroutine **PRS2EOS** - given $p$, it computes the $\rho$ according to the tabulated EoS.
    - subroutine **ENT2EOS** - given $h$, it computes the $\rho$ according to the tabulated EoS.
    - subroutine **EOS** - computes the density and the internal energy given the pressure, both in case the EoS is tabulated or an analytical polytropic.
    - subroutine **FUNCD_EOS** - used by the root-finding subroutine to derive the central pressure given
    the central density.
<br><br>
- **FUNCTIONS.f90**
    - subroutine **DGTSV** - solves the linear system $AX = B$, where $A$ is a tridiagonal matrix, by
    Gaussian elimination with partial pivoting (taken from the LAPACK routines).
    - subroutine **LUSOLVER** - solves the linear system $AX = B$, where $A$ is an $N \times N$ matrix and $B$ is a vector of length $N$, with an LU decomposition.
    - subroutine **MYSWAP** - exchanges two arrays, adapted from LAPACK 3.2
    - function **MYISAMAX** - returns the index of maximum value of an array, adapted from LAPACK 3.2
    - subroutine **SGER** - performs the rank 1 operation A := $x \times y^T$ + A, where  x is an nn element vector, y is an n element vector and A is an nn by n matrix. Adapted from LAPACK 3.2 
    - subroutine **SLASWP** - exchange elements of an array,  adapted from LAPACK 3.2
    - subroutine **STRSM** solve A*X = B where A is a matrix in LU form, adapted from LAPACK 3.2
    - subroutine **LAGRANGEINT** - interpolate YP(XP) in XINT using (N-1)th-Polinomial according to Lagrange Interpolation with Neville algorithm
    - subroutine **LAGRANGEINT2D** - perform 2D interpolation using (N-1)th-Polinomial according to Lagrange Interpolation with Neville algorithm, on structured grid
    - subroutine **LEGZO** - computes the zeros of Legendre polynomials and the corresponding
    weights for Gaussian quadrature integration. From: SPECIAL_FUNCTIONS by Shanjie Zhang, Jianming Jin.
    - subroutine **LPN** - computes the Legendre polynomials and their derivatives. From: SPECIAL_FUNCTIONS by Shanjie Zhang, Jianming Jin.
    - subroutine **POLINT** - a polynomial 2nd order interpolation routine (e.g. Numerical Recipes).
    - subroutine **XNS2ECHO_OUT** - output file for X-ECHO
## Outputs

- **Grid.dat** - contains the grid mesh points.
<br><br>
- **TOVINIMOD_PROFILES.dat** - contains the 1D TOV solution $(r,\mu,\rho,\nu,p,\rho\varepsilon,\chi)$.
<br><br>
- **Source.dat** - contains 2D source term for the metric solver $(\rho,p,\rho\varepsilon)$.
<br><br>
- **XShiftphi.dat** - contains the $W^\phi$ component and the related source term of its vector Poisson
equation.
<br><br>
- **Conformal.dat** - contains $\psi$ and the two (matter and curvature) source terms of its scalar Poisson
equation.
<br><br>
- **Primitive.dat** - contains the primitive variables $(\rho,p,\rho\varepsilon,v^\phi,B^\phi)$ recovered self-consistently from
the metric and the conserved variables.
<br><br>
- **Primitive_mag.dat** - contains the magnetic primitive variables $(B^\phi,B^r,B^\theta)$ recovered self-consistently
from the metric and the conserved variables.
<br><br>
- **Lapse.dat** - contains $\alpha$ and the two (matter and curvature) source terms of the related scalar Poisson
equation.
<br><br>
- **Shiftphi.dat** - $\beta ^\phi$ vector and the related source term of the vector Poisson equation.
<br><br>
- **Chi.dat** - $\chi$ and the related source term of the scalar field equation.
<br><br>
- **Hydroeq.dat** - contains the new equilibrium configuration $(\rho,p,\psi,v^\phi,\alpha,\beta^\phi,\chi,Q^r,Q^t)$.
<br><br>
- **Hydroeq_mag.dat** - contains the new equilibrium configuration for magnetic field $(B^\phi,B^r,B^\theta,\tilde A ^\phi,E_\phi,E_r,E_\theta,A^t,J^\phi,J^r,J^\theta)$.
<br><br>
- **Mxwll_test.dat** - cointains data related to the source term of both Maxwell-Ampere and Maxwell-Gauss equation $(\rho _\mathrm{e},J^\phi,\Phi _\mathrm{int},\Phi _\mathrm{ext},\Phi _\mathrm{a},\omega,\Gamma)$.
<br><br>
- **Apconv.dat** - maximum value of $\Psi$ at each step of the XNSMAIN subroutine.
<br><br>
- **Atconv.dat** - maximum value of $\Phi$ at each step of the XNSMAIN subroutine.
<br><br>
- **Chiconv.dat** - maximum value of $\chi$ at each step of the XNSMAIN subroutine.
<br><br>
- **Rhovec.dat** - central density and other quantities at each step of the XNSMAIN subroutine. It is used to check convergence.
<br><br>
- **Surf.dat** - contains the radius of the NS surface at all angles $\theta$.
<br><br>
- **LogFile.dat** - summary of the run (input and output quantities).

## Visualisation

- **starplot.py** - plots fluid, magnetic or metric quantities in a section of the star in the $x-z$ plane (i.e. contours of either the poloidal magnetic field, with or without the field lines, toroidal magnetic field, the scalar field, the density, etc...).
- **starcut.py** - plots the radial profiles, both at the pole and at the equator, of several quantities: $\rho$, $p$, $\psi$, $\alpha$, $\chi$, $B_\mathrm{pol}$, etc....
- **startov.py** - plots the radial profiles, of several quantities: $\rho$, $p$, $\psi$, $\alpha$, $\chi$, for the TOV initial solution.

## Reformatting the Equation of State

- **resample.py** - PYTHON script that takes an EoS (typically in COMPSTAR format, but can be any other), and rebuilds a table in the XNS 4.0 format.

