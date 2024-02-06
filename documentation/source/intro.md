# Welcome to the XNS 4.0 documentation!

The XNS code solves for axisymmetric equilibria of polytropic magnetized and/or rotating neutron stars (NSs) using the extended conformally flat condition (XCFC) for the metric, in spherical coordinates. This is based on the metric module and the routines developed for the X-ECHO code for GRMHD in dynamical spacetimes (Bucciantini & Del Zanna 2011), which in turn is an upgrade of the Eulerian conservative high-order code (ECHO Del Zanna et al. 2007) for GRMHD in a static background metric (the so-called Cowling approximation). Like ECHO and X-ECHO, also XNS is written in the Fortran90 programming language. The reader is referred to the above cited paper for full derivation of the GRMHD equations, and for the full description of the XCFC solvers.

The 4.0 version of XNS adds the possibility of solving for the structure of a magnetised, rotating NS in a class of theories of gravity alternative to general relativity (GR), that is massless scalar tensor theories (STTs). Moreover, support for the use of realistic, tabulated equations of state (EoS) has been added.

You can download it via Github at (<a href="https://github.com/niccolo-bucciantini/XNS4.0">link</a>)

You can download the pdf of the guide from Github at (<a href="https://github.com/niccolo-bucciantini/XNS4.0/blob/main/documentation/build/latex/xns.pdf">link</a>)

This guide is based on the following papers, where the equations describing the approach for magnetized and rotating models are fully presented:

- Bucciantini N., & Del Zanna L., 2011, A&A, 528, A101 (<a href="https://www.aanda.org/articles/aa/abs/2011/04/aa15945-10/aa15945-10.html">link</a>).
- Pili A.G., Bucciantini N., & Del Zanna L., 2014, MNRAS, 439, 3541 (<a href="https://academic.oup.com/mnras/article/439/4/3541/1161070">link</a>).
- Pili A.G., Bucciantini N., & Del Zanna L., 2015, MNRAS, 447, 2821 (<a href="https://academic.oup.com/mnras/article/447/3/2821/2892871?login=true">link</a>).
- Pili A.G., Bucciantini N., & Del Zanna L., 2017, MNRAS, 470, 2469 (<a href="https://academic.oup.com/mnras/article/470/2/2469/3820933?login=true">link</a>).
- Soldateschi, J., Bucciantini, N., & Del Zanna, L. 2020, A&A, 640, A44 (hereafter SBD20) (<a href="https://www.aanda.org/articles/aa/abs/2020/08/aa37918-20/aa37918-20.html">link</a>).
- Soldateschi, J., Bucciantini, N., & Del Zanna, L. 2021, A&A, 645, A39 (<a href="https://www.aanda.org/articles/aa/abs/2021/01/aa38826-20/aa38826-20.html">link</a>).
- Soldateschi, J., Bucciantini, N., & Del Zanna, L. 2021, A&A, 654, A162 (<a href="https://www.aanda.org/articles/aa/abs/2021/10/aa41448-21/aa41448-21.html">link</a>).
- Franceschetti, K., Del Zanna, L., Soldateschi, J., & Bucciantini, N., 2022, Universe, 8, 172 (<a href="https://www.mdpi.com/2218-1997/8/3/172/htm">link</a>).

If you use this software please reference at least one of the previous papers.
