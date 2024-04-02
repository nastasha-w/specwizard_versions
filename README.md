# specwizard_versions: DO NOT USE THIS BRANCH 
This is a repository for storing (and hopefully merging) different versions of specwizard and/or version control for your own copies of it.
This branch (vbug2) was created to fix a unit error in the velocities for short spectra. (The handling of units based on the data in the EAGLE output files was ok, but the units were recorded incorrectly in those files.) **This branch** has been merged with the [nastasha](https://github.com/nastasha-w/specwizard_versions/tree/nastasha) branch, and **should not be used further**. 

For notes on installing SpecWizard, see the user guide.

this branch: 
-----------
version that works on cosma5 (Durham)

finding the right hdf5 version is an issue, addressed in the makefiles here and in the hdf5_wrapper (not included in this repository)

modules:
--------
  module load intel_comp/2018-update2
  
  module load intel_mpi/2018
  
  module load hdf5/1.8.20


