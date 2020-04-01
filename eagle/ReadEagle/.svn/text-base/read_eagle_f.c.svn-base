#include <stdio.h>
#include <string.h>
#include "read_eagle.h"

/*
  Fortran wrappers for read_eagle C functions.

  Intended to be called via read_eagle_fortran.f90.

  Here we assume that strings have been null terminated
  by the fortran caller. Everything else is passed by
  reference.

  Underscores are removed from names to reduce name mangling
  problems.
*/

void geterrorf_(int *len, char *str)
{
  strncpy(str, get_error(), (size_t) *len);
}

void abortonerrorf_(int *flag)
{
  abort_on_error(*flag);
}

void opensnapshotf_(EagleSnapshot **snap, char *fname, double *boxsize,
		    long long *numpart_total, int *numfiles, int *hashbits)
{
  int i;

  *snap = open_snapshot(fname);
  if(*snap)
    {
      *boxsize = (*snap)->boxsize;
      *numfiles = (*snap)->numfiles;
      *hashbits = (*snap)->hashbits;
      for(i=0;i<6;i+=1)
	numpart_total[i] = (*snap)->numpart_total[i];
    }
}

void closesnapshotf_(EagleSnapshot **snap)
{
  close_snapshot(*snap);
}

void selectregionf_(EagleSnapshot **snap, 
		    double *xmin, double *xmax,
		    double *ymin, double *ymax,
		    double *zmin, double *zmax)
{
  select_region(*snap,
		*xmin, *xmax,
		*ymin, *ymax,
		*zmin, *zmax);
}

void setsamplingratef_(EagleSnapshot **snap, 
		       double *rate)
{
  set_sampling_rate(*snap, *rate);
}

void clearselectionf_(EagleSnapshot **snap)
{
  clear_selection(*snap);
}

void countparticlesf_(long long *n, EagleSnapshot **snap, int *itype)
{
  *n = count_particles(*snap, *itype);
}

void getparticlelocationsf_(long long *n, EagleSnapshot **snap, int *itype,
			    int *file_index, int *file_offset, int *nmax)
{
  *n = get_particle_locations(*snap, *itype, file_index, file_offset, 
			      (size_t) *nmax);
}

void readdatasetf_(long long *nread, EagleSnapshot **snap, int *itype, int *typecode, void *buf, int *n, char *name)
{
  hid_t dtype_id;
  if(*typecode==0)
    dtype_id = H5T_NATIVE_INT;
  else if(*typecode==1)
    dtype_id = H5T_NATIVE_LLONG;
  else if(*typecode==2)
    dtype_id = H5T_NATIVE_FLOAT;
  else if(*typecode==3)
    dtype_id = H5T_NATIVE_DOUBLE;
  
  *nread = read_dataset(*snap, *itype, name, dtype_id, buf, *n);
}

