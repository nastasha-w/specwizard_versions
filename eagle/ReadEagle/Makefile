#
# Makefile for C/Fortran Eagle read routines
#
# For python use the setup.py script to build the python module
#
# For IDL use read_eagle.pro (doesn't need any compilation)
#

CC      = gcc
CFLAGS  = -g -DH5_USE_16_API

FC      = gfortran
FCFLAGS = -g

# Only need these for Intel compiler
#EXTRA_LIBS = -limf -lm

.SUFFIXES:
.SUFFIXES: .f90 .c .o

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $*.o

.f90.o:
	$(FC) $(FCFLAGS) -c $< -o $*.o

all:	test_c test_fortran

clean:
	\rm -rf *.o *.mod test_c test_fortran read_eagle.tar.gz build *.pyc

TEST_C_OBJS = read_eagle.o test_c.o

test_c:	$(TEST_C_OBJS)
	$(CC) $(CFLAGS) $(TEST_C_OBJS) $(LDFLAGS) $(EXTRA_LIBS) -lhdf5 -o test_c

TEST_FORTRAN_OBJS = read_eagle.o read_eagle_f.o read_eagle_fortran.o test_fortran.o

test_fortran:	$(TEST_FORTRAN_OBJS)
	$(FC) $(FCFLAGS) $(TEST_FORTRAN_OBJS) $(LDFLAGS) $(EXTRA_LIBS) -lhdf5 -o test_fortran

DIST_SRC =  Makefile _read_eagle.c read_eagle.c read_eagle_f.c read_eagle_fortran.f90 read_eagle.h read_eagle.pro read_eagle.py README setup.py test_c.c test_fortran.f90 test_py.py extract_region.py

read_eagle.tar.gz:	$(DIST_SRC)
	\rm -rf read_eagle
	\mkdir -p read_eagle
	\cp $(DIST_SRC) ./read_eagle/
	\tar cvfz read_eagle.tar.gz read_eagle
	\rm -rf read_eagle

dist:	read_eagle.tar.gz
