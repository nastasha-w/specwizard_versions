# specwizard_versions: STANDARD FOR THIS REPO, NOT NEWEST VERSION
This is a repository for storing (and hopefully merging) different versions of specwizard and/or version control for your own copies of it. This branch (`nastasha`) is the one to use from this repository, but I am no longer maintaining this code. The newer version is at https://github.com/specwizard/specwizard .

For notes on installing SpecWizard, see the user guide.

! This Repo does not contain the most up-to-date version of specwizard !
------------------------------------------------------------------------
The newer version of this code (in python) by Andrés Arámburo García: https://github.com/specwizard/specwizard

This branch: 
-----------
version that works on cosma5 (Durham)

Includes the bug fix from the vbug2 branch (merged into this one, except for a different note in its README.md file). The MPI-involved loop for writing spectra created in parallel to a single output file fails on some systems (e.g., Leiden's quasar). The nastasha\_quasar branch has a fix for this, but note that this branch does still includes the velocity units bug, fixed in the vbug2 branch and merged into this branch. Therefore, it would be better to apply the parallel file writing fix to this branch (or a branch created from this one) instead.

When compiling, finding the right hdf5 version is an issue, addressed in the makefiles here and in the hdf5_wrapper (not included in this repository).

modules:
--------
  module load intel_comp/2018-update2
  
  module load intel_mpi/2018
  
  module load hdf5/1.8.20


