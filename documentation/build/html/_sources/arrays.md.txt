# Arrays

The module system inside the file *SYSTEMXNS.f90* also contains the definitions of some arrays that
are used within the code and shared by many subroutines. We briefly describe some of them here so that the
user can have some idea of what represents what. Other arrays that are specific only to certain subroutines
are defined locally and are not discussed here. Notice that some arrays related to the poloidal components
of the velocity, shift vector, or auxiliary vectors are always zero but still defined in XNS, though related
routines are never called. This is because XNS shares the same metric solver as the full X-ECHO code.

## Grid

- **R**,**DR** - 1D arrays that store the location of the radial points, and the radial increments.
<br><br>
- **DRM**,**DRP** - additional 1D arrays that store the backward and forward radial increments.
<br><br>
- **TH**,**DTH**,**XX** - 1D arrays that store the location of the angular points, the angular increments, and the
cosine of the angle.

## Metric

- **PSI**,**PSL**,**PSS**,**PSSR**,**PSST** - 2D arrays of metric terms, respectively $\psi$, $\alpha \psi$, $X^{\phi}$, $X^{r}$, $X^{\theta}$. We have either $X^i=W^i$, or $X^i = \beta ^i$, depending on the step of the metric solver.
<br><br>
- **CURVC**,**CURVR**,**CURVT**,**CURVP** - 2D arrays containing the source terms associated with the curvature
of the metric, respectively for the two scalar Poisson equations (for $\psi$ and $\alpha$) and for the three
components of the second vector Poisson equation (that for $\beta ^i$).
<br><br>
- **MU**,**NU** - 1D arrays containing the metric terms of the radial TOV solution. The metric employed is
that for isotropic coordinates, namely $ds^2 = -e^\nu dt^2 + e^\mu (dr^2 + r^2 d\theta ^2 + r^2 \sin ^2 \theta d\phi ^2)$.

## MHD

- **RHOSRC**,**ESRC**,**PSRC**,**VPHI**,**VR**,**VTH**,**BPHI**,**SSS** - 2D arrays, respectively $\rho$, $\rho h = \rho (1 + \varepsilon) + p$, $p$, $v^\phi$, $v^r$, $v^\theta$, $B^\phi$, $S$, needed for the source terms.
<br><br>
- **RHOTV**,**PRTV**,**ETV** - 1D arrays containing the fluid variables of the radial TOV solution, respectively
$\rho$, $p$, $\rho \varepsilon$.
<br><br>
- **RHONEW**,**PNEW**,**ENEW**,**V3NEW**,**B3NEW**,**E3NEW** - 2D arrays, respectively $\rho$, $p$, $\rho \varepsilon$, $v^\phi$, $B^\phi$, $E_\phi$ computed for an equilibrium configuration on the metric at the end of each step of the convergence loop.
<br><br>
- **BPOLR**,**BPOLT**,**EPOLR**,**EPOLT**,**APHI**,**ATIM** - 2D arrays, respectively, $B^r$, $B^\theta$, $E_r$, $E_\theta$, $\tilde A ^\phi = \Phi$, $A^t = \Psi$ for the magnetic configuration with poloidal field components.
<br><br>
- **RHOTVJOR**,**PRTVJOR**,**ETVJOR** - 1D arrays, respectively $\rho$, $p$, $\rho \varepsilon$ of the TOV solution computed in the Jordan frame. Note that, in GR, these are equal to RHOTV, PRTV, ETV respectively.

## Scalar field

- **CHITV**,**DCHITV**,**DDCHITV** - 1D arrays, respectively $\chi$, $\partial _r \chi$, $\partial ^2 _r \chi$ of the TOV solution.
<br><br>
- **CHI** - 2D array containing the scalar field $\chi$.
<br><br>
- **PSCAL**,**QSCALTIM**,**QSCALR**,**QSCALT**,**QSCALP**,**QSCAL2** - 2D arrays, containing the 3+1 decomposition of the scalar field derivative. These are respectively $P$, $Q^\mu$ and $Q_\mu Q^\mu$ of Eq. 26 in SBD20.
<br><br>
- **ASCAL** - 2D array, containing the scalar coupling function $\mathcal{A}(\chi)$.

## Source terms

- **USRC**,**DSRC**,**S3SRC**,**S1SRC**,**S2SRC** - 2D arrays containing the U conservative variables needed for
the source terms, respectively $\hat E$, $\hat D$ , $\hat S _\phi$ , $\hat S _r$, $\hat S _\theta$, all multiplied by $f ^{1/2} = r^2 \sin \theta$.
<br><br>
- **USRCX**,**S3SRCX**,**S1SRCX**,**S2SRCX** - 2D arrays containing the U conservative variables associated to the scalar field, needed for the source terms, respectively $\hat E$, $\hat D$ , $\hat S _\phi$ , $\hat S _r$, $\hat S _\theta$, all multiplied by $f ^{1/2} = r^2 \sin \theta$.
<br><br>
- **ECSRC**,**ELSRC**,**ES1RC**,**ES2RC**,**ES3RC** - 2D arrays containing the source terms (the right hand side of
the equations) associated with the presence of matter in the elliptic PDEs. Respectively, the source
for the equations for $\psi$, $\alpha \psi$, $X^{\phi}$, $X^{r}$, $X^{\theta}$, where $X^i=W^i$, or $X^i = \beta ^i$.
<br><br>
- **TRACEM** - 2D array containing the trace of the energy-momentum tensor of matter, used for the source of the equation of the scalar field.
