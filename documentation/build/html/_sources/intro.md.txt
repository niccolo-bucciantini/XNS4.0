# Welcome to the XNS 4.0 documentation!

The XNS code solves for axisymmetric equilibria of polytropic magnetized and/or rotating neutron stars (NSs) using the extended conformally flat condition (XCFC) for the metric, in spherical coordinates. This is based on the metric module and the routines developed for the X-ECHO code for GRMHD in dynamical spacetimes (Bucciantini & Del Zanna 2011), which in turn is an upgrade of the Eulerian conservative high-order code (ECHO Del Zanna et al. 2007) for GRMHD in a static background metric (the so-called Cowling approximation). Like ECHO and X-ECHO, also XNS is written in the Fortran90 programming language. The reader is referred to the above cited paper for full derivation of the GRMHD equations, and for the full description of the XCFC solvers.

The 4.0 version of XNS adds the possibility of solving for the structure of a magnetised, rotating NS in a class of theories of gravity alternative to general relativity (GR), that is massless scalar tensor theories (STTs). Moreover, support for the use of realistic, tabulated equations of state (EoS) has been added.

This guide is based on the following papers, where the equations describing the approach for magnetized models are fully presented:

- Bucciantini N., & Del Zanna L., 2011, A&A, 528, A101.
- Pili A.G., Bucciantini N., & Del Zanna L., 2014, MNRAS, 439, 3541.
- Pili A.G., Bucciantini N., & Del Zanna L., 2015, MNRAS, 447, 2821.
- Pili A.G., Bucciantini N., & Del Zanna L., 2017, MNRAS, 470, 2469.
- Soldateschi, J., Bucciantini, N., & Del Zanna, L. 2020, A&A, 640, A44 (hereafter SBD20).
- Soldateschi, J., Bucciantini, N., & Del Zanna, L. 2021, A&A, 645, A39.
- Soldateschi, J., Bucciantini, N., & Del Zanna, L. 2021, A&A, arXiv:2106.00603v2 [astro-ph.HE].

If you use this software please reference at least one of the previous papers.
