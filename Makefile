.SUFFIXES:	.c .f90 .o

FOPTS  = -g -traceback -132 -fpp
LOADOPTS = $(FOPTS)

MAIN    = a.out

DEPS    = metis_interface, mpi_tools, domain, telecom, main

OBJS    = $(DEPS:,=.o).o

CC = mpicc
FORTRAN	= mpiifort 
LDR	= mpiifort  

METIS_PATH = metis

LIBPATH = -L${METIS_PATH}/lib -lmetis
INCPATH	= -I${METIS_PATH}/include

default:	all	

all:	$(MAIN)

$(MAIN):	$(OBJS)
		$(LDR) $(LOADOPTS) $(OBJS) $(LIBPATH) $(LIBS) -o $(MAIN) 

.f90.o:
		$(FORTRAN) $(INCPATH) $(DEFINES) $(FOPTS) -o $*.o -c $<

clean:
		rm -rf *.o *.mod $(MAIN)

rebuild: clean all
