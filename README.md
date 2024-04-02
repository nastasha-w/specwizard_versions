# This branch: contains velocity units bug

This branch was set up to fix an issue on the Leiden quasar system, where writing the spectra from multiple processes to a single output file caused problems, despite MPI barriers meant to prevent the issue. However, a bug in the handling of velocity units in the EAGLE short spectrum files cropped up later. This has been fixed in the [nastasha](https://github.com/nastasha-w/specwizard_versions/tree/nastasha) branch, but not in this one. (I made an attempt in the vbugfix branch, but it stranded.) To deal with an issue arising from different processes trying to write to the same file, please apply the fix from this branch to the nastasha branch (or a different branch created from it).


Instructions for compiling SpecWizard:
--------------------------------------

  SpecWizard requires the following libraries:

      HDF5
      HDF5_Wrapper (available from us)

  Install these libraries and then create a make.$YOURCOMPUTERNAME$
  file which specifies the correct compilers and library locations for
  your machine. (Use one of the pre-existing make.* files as a template).

  Change the first line of Makefile so that it uses
  make.$YOURCOMPUTERNAME$.

  Copy your Makefile and make.* file to src/

  Fingers crossed you can then make the code with gmake.

  In case of technical problems ask Craig Booth
  (booth@strw.leidenuniv.nl) or, if you are in Durham, Tom Theuns.

Refer to ../UserGuide for detailed usage instructions
