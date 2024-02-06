# User parameters

The input parameters that the user might want to change are all set in the module systemxns inside the file *SYSTEMXNS.f90*. There are other parameters in other parts of the code that deal with specific routines (root finding, convergence etc ...), but those should not need to be changed. Here is a list of the parameters of the model as they appear in the module system:

## Convergence

### Parameters for convergence of the Newton-Raphson scheme of XNS

- **NVALUE** - the maximum number of loops employable by the Newton-Raphson scheme in the search
for a equilibrium solution, having a target value for a desired quantity (central density, total mass,
etc..) by the program XNS. Usually convergence is reached within about 10 steps, unless the NS is
strongly distorted (fast rotation, and/or strong magnetic field). The default value is set to 100.
<br><br>
- **CONVF** - a convergence parameter for the Newton-Raphson scheme in XNS. It is given in relative terms
(beware that the code accuracy is $\sim 10^{-3}$).
<br><br>

### Parameters for convergence of XNSMAIN over a model

- **MAXLOOP** - the maximum number of loops employable in the search for a converged equilibrium
solution by the subroutine XNSMAIN. Usually convergence is reached within the first 100 steps,
unless the NS is strongly distorted (fast rotation, and/or strong magnetic field). The default value is
set to 100.
<br><br>
- **CONVRHO** - sets if convergence by XNSMAIN should be checked on central density or on central lapse function. CONVRHO should be set to .FALSE. if the central density is held fixed, as in the case of Tabulated EoS 
<br><br>
- **OSCCONV** - desired precision in case the solution is reached with oscillatory convergence in XNSMAIN. In this case the check on the convergence is done between a step and the previous few, in such a way to bound the oscillation range.
<br><br>
- **MONCONV** - desired precision in case the solution is reached with monotonous convergence in XNSMAIN. In this case the check on the convergence is done between a step and the previous. 
<br><br>

### Parameters for convergence of the initial TOV solution

- **RELIT** - the maximum number of loops in TOVINIMOD for the convergence of the full TOV solution, including the metric-matter distribution and, in case of STT, the scalar field.
<br><br>
- **CONVT** - desired precision for the convergence of the full TOV solution and also of the scalar field in TOVINIMOD, in STT.
<br><br>
- **CONV** - desired precision for the convergence of the GR TOV part (or in STT the matter-metric part at fixed scalar field) in TOVINIMOD.
<br><br>
- **MAXSTEPTV** - the maximum number of loops in TOVINIMOD for the convergence of the TOV system to a solution with a fixed scalar field. Note that the full TOV solution is obtained when also the scalar field has reached convergence together with the rest of the quantities.
<br><br>
- **MAXSTEPCH** - the maximum number of loops in TOVINIMOD for the convergence of the scalar field equation of the TOV system to a solution with a fixed metric and matter distribution. Note that the full TOV solution is obtained when also the scalar field has reached convergence together with the rest of the quantities.
<br><br>

### Parameters for convergence of the Bernoulli in HYDROEQ

- **CONVHELP** - a logical flag to activate an option in HYDROEQ to help achieve convergence. Warning: if set to .TRUE., the final central density will be slightly different from the chosen RHOINI.

### Parameters for the convergence of elliptic solvers

- **TOLCONV** - the convergence tolerance for the iterative solution of the PDEs for the conformal factor $\psi$
and the lapse $\alpha$ in XNSMAIN and for the 4-potential in HYDROEQ.
<br><br>
- **TOLCHI** - desired precision for the scalar field iterative solver in XNSMAIN.
<br><br>

### Parameters to help or stabilize convergence

- **QFACTOR** - a damping factor of the convergence loop both for solving the Bernoulli equation, used in HYDROEQ.
At the end of each sub-loop of the convergence scheme, a new set of equilibrium fluid variables is
computed. Setting QFACTOR=1 implies that these will be used, while setting QFACTOR=0 means that
the old variables $V_\mathrm{old}$ will be used (the code will never converge in this case!). A value $0 < Q_\mathrm{f} < 1$
implies that at the beginning of each loop a combination of new and old variables will be used, in the
form $V = Q_\mathrm{f}V_\mathrm{new} + (1 - Q_\mathrm{f})V_\mathrm{old}$.
Using a value less than 1 tends to give slower but more stable convergence. Values QFACTOR$ < 0.5$
are to be used only for pathological cases where the convergence is very slow or when the code fails
to converge (i.e. rotating models on the unstable branch of NS mass-radius curve).
<br><br>
- **QAPHI** - a damping factor of the convergence loop for solving the Grad-Shafranov or the Maxwell-Ampere equation. Analogous to QFACTOR but for the $\phi$-component of the 4-potential.
<br><br>
- **QRELAX** - a damping factor for the convergence loop of the scalar field in TOVINIMOD.
<br><br>
- **QFACTORMETRIC** - a damping factor for the internal convergence loop for solving the conformal factor equation in XNSMAIN.
<br><br>
- **QFACTORCONF** - a damping factor for the conformal factor in the convergence loop of XNSMAIN.
<br><br>
- **QFACTORCHI** - a damping factor for the internal convergence loop for solving the scalar field equation in XNSMAIN.
<br><br>
- **EPS** - a tolerance value. It is used in several subroutines and must be a small value. This should not
need to be changed.
<br><br>


## Grid
- **STRETCH** - a logical flag that control whether the grid is stretched or not. If .TRUE. the radial grid
is regular up to RREG with NRREG grid points and it is stretched from RREG to RMAXSTR with NR-NRREG
points. The stretching factor STRR is determined by the code determined by the code consistently with
the choices for NR, NRREG,RMAXSTR and RREG. See also Pili et al. (2015) for details.
- **NR** - the number of radial grid points. The radial grid is defined from $r=$RMIN and $r=$RMAX ($r=$RMAXSTR) is the grid is uniform (stretched).
<br><br>
- **NTH** - the number of angular grid points. The angular grid is always defined between $\theta = 0$ and $\theta = \pi$.
<br><br>
- **NRREG** - number of grid points for the regular grid if STRETCH=.TRUE..
<br><br>
- **RMIN** - the lower boundary in the radial direction. It must be always set to 0, since the metric solver
requires a compact domain and has been implemented with specific boundary conditions for RMIN = 0.
<br><br>
- **RMAX** - the maximum radius of the computational domain if the grid is uniform. This can be arbitrarily chosen. However,
one needs to guarantee that the NS is properly resolved over a sufficient number of grid points (50-
100), so this parameter and NR should be chosen consistently. In particular, the condition RMAX$ >
2R_\mathrm{TOV}$ must hold, where $R_\mathrm{TOV}$ is the NS radius of the initial TOV guess. This is because the TOV
solver is designed to converge when the ADM masses measured at RMAX and RMAX/2 (hence it must
be outside the NS) coincide within a given tolerance. If not the code will halt with a warning.
<br><br>
- **RMAXSTR** - the maximum radius of the computational domain if the grid is stretched.
<br><br>
- **RREG** - maximum radius of the regular grid.
<br><br>
- **REQMAX** - the maximum radius beyond which any NS model will be artificially truncated. This must
be set $\geq$ of the shedding mass limit. Sometimes, when working with configurations close to mass
shedding, during the convergence loop the code might get unbounded solutions or fail to converge. To
avoid this, setting a value for REQMAX will force the solution to be truncated.
<br><br>
- **MINRESREG** - desired minimum resolution of the grid, if uniform. In case the resolution is too rough, it allows to issue a warning.
<br><br>
- **MINRESSTR** - desired minimum resolution of the regular part of the grid, if stretched. In case the resolution is too rough, it allows to issue a warning.
<br><br>
- **RINI** - a very small radius used for the expansion of the TOV equations.
<br><br>

## Printing outputs and files

- **VERBOSE** - a logical flag. Setting VERBOSE=.TRUE. forces the code to output on screen all the INFOs
related to the various steps done by each subroutine (to be used only for debugging or checks). Otherwise setting VERBOSE=.FALSE. the output on screen will be produced only at the end. The latter is
the default option.
<br><br>
- **DEBUG** - additional logical flag used to print more outputs. Used for debugging purposes.
<br><br>
- **WRT** - a logical flag. Setting WRT=.TRUE. forces the code and each subroutine to write output files
at every step or substep, otherwise setting WRT=.FALSE. will prevent IO writing.
<br><br>
- **WRTF** - a logical flag. Setting WRTF=.TRUE. override WRT for the final step, and allow to write all the
files related to the final configuration. Setting WRTF=.FALSE. will prevent IO writing.
<br><br>
- **CHUP** - a logical flag. Setting CHUP=.TRUE. allows (subject to WRT, WRTF) to write the files containing the results of the metric solver and primitive solver XShiftphi.dat, Conformal.dat, Primitive.dat,
Primitive_mag.dat, Shiftphi.dat, Lapse.dat, Source.dat, Rhovec.dat, Chi.dat. Setting CHUP=.FALSE. will prevent from writing these files. The latter is the default option, unless one wishes to perform a check of the metric or primitive solvers.
<br><br>
- **WGRID** - a logical flag. Setting WGRID=.TRUE. writes the grid in a file named Grid.dat.
<br><br>
- **WCONVA** - a logical flag. Setting WCONVA=.TRUE. writes two files named Apconv.dat and Atconv.dat used to check the convergence of the $t$ and $\phi$ components of the 4-potential.
<br><br>
- **WCONVC** -  a logical flag. Setting WCONVC=.TRUE. writes a file named Chiconv.dat used to check the convergence of the scalar field.

## Rotation

**Warning: rotation in concurrency with the presence of a scalar field has not been tested!**

- **OMG** - the value of the angular velocity at the center $\Omega _\mathrm{c}$.
<br><br>
- **A2VALUE** - the value of A2. This is needed only for differentially rotating models, otherwise it should
be set to 0.
<br><br>
- **DIFFERENTIAL** - a logical flag that states whether the model is differentially rotating or not. Setting it
to .FALSE. implies uniform rotation, with $\Omega = $OMG. Setting it to .TRUE. implies differential rotation.
In this case, a value of A2VALUE must be specified.
<br><br>
- **OMGSPACE** - Logical variable to work in $\Omega$ or $J$ space. If .TRUE. the rotational law has form $J(\Omega)$, else $\Omega(J)$.
<br><br>
- **JCONSTLAW** - Rotation law: $J = A^2(\Omega_c-\Omega)$
<br><br>
- **JCMODLAW** - Rotation law: $J= A^2\Omega[(\Omega_c/\Omega)^p-1]$
<br><br>
- **URYULAW3** - Rotation law: Uryu  3, $\Omega = \Omega_c \, [1+(\dfrac{j}{B^2 \, \Omega_c})^p](1-\dfrac{j}{A^2 \, \Omega_c})$
<br><br>
- **URYULAW4** - Rotation law: Uryu  4, $\Omega = \Omega_c \, \dfrac{1+(j/B^2 \, \Omega_c)}{1+(j/A^2 \, \Omega_c)^{4}}$
<br><br>
- **PROTDIFF** - The index $p$ in either JCMODLAW or URYULAW3.
<br><br>
- **OMGMAX** - The maximum value of the rotation rate for Uryu  3 and 4 (this is given instead of A and B)
<br><br>
- **RMVALUE** - The radius at which the maximum value of the rotation rate for Uryu  3 and 4 is reached (this is given instead of A and B)

## Magnetic field

**Warning: a twisted torus configuration in concurrency with the presence of a scalar field has not been tested!**

- **IMAG** - a logical flag that states whether the model is magnetized or not. Setting it to .FALSE. implies
the non magnetized case. Setting it to .TRUE. implies the presence of a magnetic field. In the latter
case, values of parameters for the magnetic model must be specified.
<br><br>
- **ITOR** - a logical flag that must be set true only for purely toroidal configurations.
<br><br>
- **IPOL** - a logical flag that must be set true only for purely poloidal configurations.
<br><br>
- **ITWT** - a logical flag that must be set true only for mixed Twisted Torus configurations.
<br><br>
- **BCOEF** - the value of Km in the magnetic polytropic law for the case of purely toroidal field. Never
used when IMAG=.FALSE. or ITOR=.FALSE., though it is better set to 0 in this case.
<br><br>
- **MAGIND** - the value of m in the magnetic polytropic law for the case of purely toroidal field. It must be $> 1$, otherwise the magnetic energy diverges on the polar axis. It is never used when IMAG=.FALSE.
or ITOR=.FALSE..
<br><br>
- **KBPOL** - the value of $K_\mathrm{pol}$ in the magnetic law for the case of purely poloidal field. Never used when
IMAG=.FALSE. or IPOL=.FALSE., though it is better set to 0 in this case.
<br><br>
- **CSI** - the value of the parameter $\xi$ (non linear current term) in the magnetic polytropic law for the
case of purely poloidal field. It is never used when IMAG=.FALSE. or IPOL=.FALSE..
- **QNULL** - logical flag that regulates the global net charge of a rotating star with poloidal magnetic field.
If QNULL=.TRUE. the code searches for a globally uncharged star, otherwise it minimizes the electric
field at the stellar pole.
<br><br>
- **KBTT** - the value of of $K_\mathrm{pol}$ in the magnetic law for the case of mixed Twisted Torus configuration.
Never used when IMAG=.FALSE. or ITWT=.FALSE., though it is better set to 0 in this case.
<br><br>
- **ATWT** - the value of a in the magnetic law for the case of mixed Twisted Torus configuration. Never
used when IMAG=.FALSE. or ITWT=.FALSE., though it is better set to 0 in this case.
<br><br>
- **ZETA** - the value of $\zeta$ in the magnetic law for the case of mixed Twisted Torus configuration.
<br><br>
- **CUT** - the value of $\lambda$ in the magnetic law for twisted magnetosphere models. It regulates the extension
of the twist. If $\lambda = 1$ standard Twisted Torus models are recovered.
<br><br>
- **NPOL** - value of the magnetic powerlaw index. Currently it is not used.

## Initial conditions

- **RHOINI** - the central density for the starting guess. If the serial flag is used to compile the code, this is also the central density of the final solution. WARNING: if the nwtrps flag is used to compile the code, this is not the central density of the converged model, but its value for the starting guess. By default XNS will search for a solution in the range 0.8-1.2 RHOINI. If the desired solution is outside this range, XNS will output a warning, and stop.
<br><br>
- **QUOC** - the quantity of interest to which the model must converge [0 for a given central density, 1 for
a given gravitational mass, 2 for a given barionic mass].
<br><br>
- **QUCONV** - the value of the quantity of interest to which we want a model to converge. For example if one wants a model with central density 1.28$\times 10^{-3}$ (in geometrized units) set: QUOC=0, QUCONV=1.28E-3.

## Equation of state

- **K1** - the polytropic coefficient K of the EoS. The code uses a polytropic EoS. The value K1=100 is
for the standard case used for many tests of NS stability and evolution in the literature (see also the
parameter below).
<br><br>
- **GAMMA** - the adiabatic index $\gamma = 1 + 1/n$ of the polytropic EoS, where $n$ is the polytropic index. The
value GAMMA=2 ($n = $1) is for the standard case used for many tests of NS stability and evolution in
the literature (see also the parameter above).
<br><br>
- **CTP** - if set to .FALSE. the code avoid to use conservative to primitive routines. This flag is effective
only with IPOL=.TRUE..
<br><br>
- **MBARYONFC** - ratio between a reduced baryon mass and true baryon mass. If equal to 1, it has no effect of the code. In case one wants to find a solution using the SQM2 EoS (tabulated in the SQ2_resampled.dat file), it is necessary to use MBARYONFC=0.86.
<br><br>
- **EOSINT** - a logical flag to choose whether o use a tabulated or a polytropic EoS. If EOSINT=.TRUE., the EoS tabulated in the file in FILEEOS is used. Otherwise, a polytropic EoS is used.
<br><br>
- **FILEEOS** - specifies the name of the file containing the tabulated EoS that we want to use.
<br><br>
- **NPTRHO** - number of points of the tabulated EoS.
<br><br>
- **EOSJOR** - a logical flag that specifies in which frame of STTs (either Jordan or Einstein) the EoS is computed. If EOSJOR=.TRUE., it is computed in the Jordan frame; otherwise, in the Einstein frame. It should always be set to .FALSE.. Note that, in GR, the two frames coincide, and so this flag has no effect on a GR solution.

## Scalar field

- **ALPHA0** - value of the $\alpha _0$ parameter of the scalar coupling function $\mathcal{A}(\chi)=\exp [ \alpha_0 (\chi - \chi _\mathrm{inf}) + \beta _0 (\chi - \chi _\mathrm{inf})^2 ]$ by Damour, T., & Esposito-Farèse, G. (1993). Note that, if set to zero, XNS converges to a GR solution.
<br><br>
- **BETA0** - value of the $\beta _0$ parameter of the scalar coupling function $\mathcal{A}(\chi)$.
<br><br>
- **CHIINF** - value of the scalar field at infinity $\chi _\mathrm{inf}$ in the scalar coupling function $\mathcal{A}(\chi)$.  
<br><br>
- **GR** - a logical flag used to choose whether to compute the solution in GR (if .TRUE.) or in STTs (if .FALSE.). Note that setting GR=.FALSE. together with ALPHA0=0 and BETA0=0 leads to the GR solution anyway.

## Other

- **RHOSURF** - value at which the surface of the NS is set.
<br><br>
- **MLS** - number of spherical harmonics (Legendre polynomials) for spectral decomposition in $\theta$ (numbered from 0 to MLS). This should be $< $NTH. In 1D (NTH=1) it must be set to 0.
<br><br>
- **NGQ** - the number of interpolation points for the Gauss quadrature, needed to compute the integrals
over the polar direction of the source terms in the spherical harmonics decomposition. Used by all the
2D elliptical PDE solvers. It must be NGQ$$NTH. In 1D (NTH=1) it must be set to 1.
<br><br>
- **MLSL** - number of spherical harmonics used to solve the Laplace equation as described in Pili et al.
(2017). This parameter is relevant only in the cases of rotating and poloidal magnetized star.
<br><br>
- **ANALYTIC** - a logical flag. If ANALYTIC=.TRUE., the derivatives of the scalar field in TOVINIMOD are computed using their analytic form; otherwise, a numerical derivative is computed. This option may be used to achieve better convergence of the TOV solution in case one approach is numerically unstable.
<br><br>
- **DELTMU0** - variation $\Delta \mu _0$ of the central value $\mu _0$ of $\mu$ (where $-e^\mu = g _{00}$ is the $00$ element of the metric tensor in the TOV metric). It is used to compute the solution at a fixed scalar field in the Newton-Raphson method in TOVINIMOD. Its value should not be changed.
<br><br>
- **MUIN** - initial guess on $\mu _0$ in the TOV solution.
<br><br>
- **MMID** - initial guess on the ADM mass of the TOV solution computed in the middle of the grid. At convergence, it must be M$\sim$MMID.
<br><br>
- **M** - initial guess on the ADM mass of the TOV solution computed at the outer edge of the grid. At convergence, it must be M$\sim$MMID.
