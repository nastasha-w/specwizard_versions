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
