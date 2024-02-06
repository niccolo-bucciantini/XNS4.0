#fortran compilers
FC    = gfortran
MPIFC = mpif90

#compiler flags
CFLAGS    =  -cpp -Ofast -m64 -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 -finit-local-zero

DEBFLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 -fimplicit-none  -Wall -Wextra -Wno-tabs -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all -fbacktrace

parspace: FC=$(MPIFC)

serial:   FLAGS=$(CFLAGS) 
nwtrps:   FLAGS=$(CFLAGS) -DNWTRPS
debug:    FLAGS=$(DEBFLAGS)
parspace: FLAGS=$(CFLAGS) -w -DMPIV 

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
