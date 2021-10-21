# Compiling XNS

XNS can be compiled using gfortran in several ways, as written in the file makefile:

- **make serial** - the standard way to compile XNS. In this case a single solution is found with a central density specified by the parameter RHOINI. The generated executable is named XNS-s. To compile the code this way type:
  ```shell
      make clean; make serial; ./XNS-s
  ```
- **make nwtrps** - compiles XNS using the Newton-Raphson method to converge on the quantity specified by the parameters QUOC and QUCONV. The generated executable is named XNS-nr. To compile the code this way type:
  ```shell
      make clean; make nwtrps; ./XNS-nr
  ```
- **make parspace** - compiles XNS various times with different initial conditions, in order to sample the parameter space spanned by the magnetic, rotation and density parameters, set in XNS.f90. In particular, solutions are computed with: NKB different magnetic coefficients spanning from KBMIN to KBMAX; NRHO1 different central densities spanning from RHOMIN to RHOMAX; NOMG different central angular velocities spanning from OMGMIN to OMGMAX. The generated executable is named XNS-mpi. The computation of the various models is made in parallel (only on CPUs, no GPU support) using the MPI framework, and mpif90 needs to be installed. To compile the code this way using NUMBER_OF_PROCESSES processes, type:
  ```shell
      make clean; make parspace; mpirun -n NUMBER_OF_PROCESSES ./XNS-mpi
  ```
  Note that it must be NUMBER_OF_PROCESSES$\geq$2, because one process is always only passing initial conditions to the other processes, and NUMBER_OF_PROCESSES-1$\leq$(NOMG+1)$\times$(NRHO1+1)$\times$(NKB+1), that is the number of computing processes must not be larger than the number of models to be computed.
<br><br>
- **make clean** - removes all .o and .mod files and all executables.
<br><br>
- **make cleanall** - removes all .o and .mod files, all executables and all .dat files.
