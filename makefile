#fortran compilers
FC    = gfortran-mp-7  
#FC = ifort
MPIFC = mpif90
#MPIFC = /usr/lib64/mpich/bin/mpif90

#compiler flags
CFLAGS    =  -cpp -Ofast -m64 -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 -finit-local-zero #-ffpe-trap=invalid,zero,overflow #-O3 -cpp -xCORE-AVX2 # -O3 -m64 -ffast-math
#CFLAGS    = -cpp -O3 -m64 -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 -finit-local-zero
#CFLAGS = -fpp -O3 -r8 -xHOST -ipo -no-wrap-margin
CMPI = -w -Ofast -cpp -m64 -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 -finit-local-zero #-cpp -fdefault-real-8 -O3
#-Wno-implicit -cpp -m64 -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8
# -g -qopt-report=5 -cpp -fopenmp
NWTRFLAGS = $(CFLAGS) -DNWTRPS
MPIFLAGS = $(CMPI) -DMPIV 

#DEBFLAGS  = -cpp -fcheck=all   -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8
DEBFLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 -fimplicit-none  -Wall -Wextra -Wno-tabs -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all -fbacktrace

parspace: FC=$(MPIFC)

serial:   FLAGS=$(CFLAGS)
nwtrps:   FLAGS=$(NWTRFLAGS)
debug:    FLAGS=$(DEBFLAGS)
parspace: FLAGS=$(MPIFLAGS)

#executables
EXESERIAL = XNS-s
EXENWTRPS = XNS-nr
EXEDEBUG  = XNS-s
EXEMPI    = XNS-mpi

#files to be compiled

SOURCES= SYSTEMXNS.f90 PHYSICS.f90 ROTATION.f90 FUNCTIONS.f90 XNS.f90 XNSMAIN.f90 TOVINIMOD.f90 HYDROEQ.f90 

#object files

OBJS=$(SOURCES:.f90=.o)

#suffix rule to generate .o files from .f90 files

%.o: %.f90
	$(FC) $(FLAGS) -c $< -o $@

#targets

help:
	@ echo
	@ echo "Standard Makefile for Fortran code"	
	@ echo "  make serial   : build standard XNS code"
	@ echo "  make nwtrps   : build XNS code with global Newton-Raphson"
	@ echo "  make debug    : to debug the standard XNS code"
	@ echo "  make parspace : build MPI code"
	@ echo "  make clean    : clean objects files"
	@ echo "  make cleanall : clean also DATA files"
	@ echo " "
	@ echo "Current settings are:"
	@ echo "  executable   :" $(EXESERIAL) $(EXENWTRPS) $(EXEMPI)
	@ echo "  sources      :" $(SOURCES)
	@ echo "  flags        :" $(CFLAGS)
	@ echo

.PHONY: clean help

serial:$(EXESERIAL)
nwtrps:$(EXENWTRPS)
parspace:$(EXEMPI)
debug :$(EXEDEBUG)

clean:
	rm -f *.o *.mod $(EXESERIAL) $(EXENWTRPS) $(EXEDEBUG) $(EXEMPI)
	@echo " OBJECT AND EXE FILES CLEANED"

cleanall:
	rm -f *.o *.mod $(EXESERIAL) $(EXENWTRPS) $(EXEDEBUG) $(EXEMPI)
	rm -f *.dat
	@echo " DATA, OBJECT AND EXE FILES CLEANED"

$(EXESERIAL):$(OBJS)
	$(FC) $(FLAGS) $(OBJS) -o $@

$(EXENWTRPS):$(OBJS)
	$(FC) $(FLAGS) $(OBJS) -o $@

$(EXEMPI):$(OBJS)
	$(FC) $(FLAGS) $(OBJS) -o $@
