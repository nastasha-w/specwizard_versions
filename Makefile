include make.cosma5_2018-2.intel
#include make.cosma5.intel
#include make.COSMA4.intel
#include make.ZIJL
#include make.TITANIA
#include make.OCTOR
#include make.QUINTOR
#include make.CLUSTER
#include make.volans

FFLAGS = $(FOPTIMIZED)
#FFLAGS = $(FDEBUG)
#CFLAGS = $(CDEBUG)

#OPTIONS  = -DEAGLE
OPTIONS  = -DMPI -DEAGLE #parallel executable
#OPTIONS  = -DMPI -DOWLS -DHUBBLE=2

EXE	= specwizard
SRCS	= specwizard.F90 specwizard_subroutines.F90 specwizard_numerical.F90 specwizard_modules.F90 w4_gadget_spline_kernel_class.F90
OBJS	= ${SRCS:.F90=.o}


# hdf5 and other modules to link against
OPTHDFW += $(ISMOD)$(HDFW_MOD)
OPTHDFW += $(ISLIB)$(HDFW_LIB)
OPTHDFW += $(ISINC)$(HDFW_MOD)
OPTHDFW += $(ISRLIB)$(HDFW_LIB)

OPT += $(OPTHDFW) -lhdfwrapper

OPTHDF += $(ISMOD)$(HDF_MOD)
OPTHDF += $(ISLIB)$(HDF_LIB)
OPTHDF += $(ISINC)$(HDF_MOD)
OPTHDF += $(ISRLIB)$(HDF_LIB)

OPT += $(OPTHDF) -lhdf5
OPTIONS += $(OPT)
# Added + (JCH)
#OPTIONS +=  -DEAGLE -DH5_USE_16_API -DREADREGION
OPTIONS +=  -DEAGLE -DH5_USE_16_API -DAliVersion -DREADREGION


.SUFFIXES:
.SUFFIXES: .o .F90

.F90.o:
	$(FC) $(FFLAGS) $(OPTIONS) $(INC) -c $<

%.o: %.F90 
	$(FC) $(FFLAGS) $(OPTIONS) $(INC) -c $<

EAGLE = 
ifeq (EAGLE,$(findstring EAGLE, $(OPTIONS)))	
  EAGLE += eagle/read_eagle.o \
		    eagle/read_eagle_f.o \
		    eagle/read_eagle_fortran.o 
  OPT += ${ISMOD}eagle
endif


APPS = screen specwizard

all: $(APPS)

# Just Dummy Reporting
#=============================================================================
screen: Makefile
	@echo
	@echo "SIMTYPE=  " $(SIMTYPE)
	@echo "MACRO=    " $(MACRO)
	@echo "FFLAGS=    " $(FFLAGS)
	@echo "OPT=    " $(OPTIONS)
	@echo





specwizard: specwizard.o specwizard_subroutines.o specwizard_numerical.o specwizard_modules.o w4_gadget_spline_kernel_class.o $(EAGLE)
	$(FC) $(FFLAGS)  specwizard.o specwizard_numerical.o specwizard_modules.o w4_gadget_spline_kernel_class.o specwizard_subroutines.o $(EAGLE) $(OPTIONS) $(INC) $(HDF5_MOD) $(LIB) $(RLIB) -o specwizard

specwizard.o: specwizard_subroutines.o specwizard_numerical.o specwizard_modules.o w4_gadget_spline_kernel_class.o specwizard.F90
	$(FC) $(FFLAGS) $(OPTIONS) $(INC) $(HDF5_MOD) specwizard.F90 -c

specwizard_subroutines.o: specwizard_numerical.o specwizard_modules.o w4_gadget_spline_kernel_class.o specwizard_subroutines.F90
	$(FC) $(FFLAGS) $(OPTIONS) $(INC) $(HDF5_MOD) $(LIB) specwizard_subroutines.F90 -c

specwizard_numerical.o: specwizard_modules.o specwizard_numerical.F90
	$(FC) $(FFLAGS) $(OPTIONS) $(INC) $(HDF5_MOD) specwizard_numerical.F90 -c

specwizard_modules.o: specwizard_modules.F90
	$(FC) $(FFLAGS) $(OPTIONS) $(INC) specwizard_modules.F90 -c

w4_gadget_spline_kernel_class.o: w4_gadget_spline_kernel_class.F90
	$(FC) $(FFLAGS) $(OPTIONS) $(INC) w4_gadget_spline_kernel_class.F90 -c

clean:
	rm -f $(EXE) *.o *.mod
