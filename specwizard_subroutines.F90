#define SCALE 2
! +++++++++++++++++++++++++++++++++++++++++ startup... +++++++++++++++++++++++++++++++++++++++++++ !

module read_startup_file
contains
  !
  subroutine read_parameter_file
    ! this routine sets default values for the run time parameters
    ! it then reads the parameter file to see whether any of these need changing
    ! Some variables have no default value, and allows need to be specified
    use numbers
    use runtime
    use spectra
    use modified_metallicity
    use noisedata
    implicit none
    !
    integer :: io_status
    integer, parameter    :: maxwords=80
    character(len=300)    :: inline, parm, value, out_par_file
    integer               :: first(maxwords), last(maxwords), nwords
    integer               :: i, line_len
    logical               :: in_quote, quoted
    character, parameter       :: quote = "'", space = ' ', & 
         tab = '     ', hash = '#', & 
         equals = '=', comment='%'
    character(len=1)            :: tf
    character(len=3), parameter :: bad = 'bad'
    !    
    open(unit=10,file=parameter_file,status='old', iostat=io_status)
    if(io_status /= 0) then
       write(*,*) parameter_file
       call abortrun(' error: parameter file not found ')
    endif
    !
    ! initialize required parameters to invalid
    ibdir           = invalid
    datadir         = invalid
#ifdef EAGLE
    snap_base       = invalid
#endif
    outputdir       = invalid
    do_long_spectrum = .false. ! default
    ! 
    do
      read(10,'(a300)',iostat=io_status) inline
      if(io_status > 0) then ! error in reading
        write (*,*) ' error in reading '
      else if (io_status < 0) then
        exit ! end of file
      else
        !
        i = index( inline, hash )
        if (i == 0) then
          line_len = len_trim(inline)
        else
          line_len = len_trim(inline(:i-1))
        endif
        if (line_len == 0) cycle
        !
        in_quote = .false.
        do i = 1, line_len
          if (inline(i:i) == quote) in_quote = .not. in_quote
          if ((inline(i:i) == tab) .and. (.not. in_quote))  & 
               &          inline(i:i) = space
        enddo
        !
        ! comment
        if(inline(1:1) == comment) cycle
        !
        call find_words (inline(1:line_len), first, last, maxwords, nwords)
        !
        if (nwords < 3) then
           !          call parm_syntax_error (inline)
          cycle
        endif
        ! second word not an =
        if (inline(first(2):last(2)) /= equals) then
          call parm_syntax_error (inline)
          cycle
        endif
        !
        parm = inline(first(1):last(1))
        !
        ! read required variables
        !write (*,*) trim(parm),' value = ',inline(first(3):last(3))
        if(parm == 'ibdir')then
          ibdir           = inline(first(3):last(3))
        else if (parm == 'datadir') then
          datadir         = inline(first(3):last(3))
#ifdef EAGLE
       else if (parm == 'snap_base') then
          snap_base  = inline(first(3):last(3))
#endif
        else if (parm == 'file_list') then
          file_list       = inline(first(3):last(3))
        else if (parm == 'outputdir')then
          outputdir       = inline(first(3):last(3))
        else if (parm == 'noisefile')then
          noisefile       = inline(first(3):last(3))
        else if (parm == 'SpectrumFile')then
          SpectrumFile = inline(first(3):last(3))
        else if (parm == 'do_long_spectrum') then
         do_long_spectrum = read_logical(inline(first(3):first(3)))
        else if (parm == 'use_snapshot_file') then
          use_snapshot_file = read_logical(inline(first(3):first(3)))
        else if (parm == 'use_maxdens_above_zmax') then
          use_maxdens_above_zmax = read_logical(inline(first(3):first(3)))
        else if (parm == 'output_frequency') then
          output_frequency = read_logical(inline(first(3):first(3)))
        else if (parm == 'overwrite') then
          overwrite = read_logical(inline(first(3):first(3)))
        else if (parm == 'NoPecVel') then
          NoPecVel = read_logical(inline(first(3):first(3)))
        else if (parm == 'zqso') then
          zqso = read_real(inline)
        else if (parm == 'zabsmin') then
          zabsmin = read_real(inline)
        else if (parm == 'zabsmax') then
          zabsmax = read_real(inline)
        else if (parm == 'nspec') then
          nspec = read_int(inline)
        else if (parm == 'minlambda') then
          minlambda = read_real(inline)
        else if (parm == 'maxlambda') then
          maxlambda = read_real(inline)
        else if (parm == 'fzresol') then
          fzresol = read_real(inline)
        else if (parm == 'ibfactor') then
          ibfactor = read_real(inline)
        else if (parm == 'use_fitted_ibfactor') then
          use_fitted_ibfactor = read_logical(inline(first(3):first(3)))
        else if (parm == 'ibfactor_he_reionization') then
          ibfactor_he_reionization = read_logical(inline(first(3):first(3)))
        else if (parm == 'minbother_red') then
          minbother_red = read_real(inline)
        else if (parm == 'minbother_blue') then
          minbother_blue = read_real(inline)
        else if (parm == 'sigtonoise') then
          sigtonoise = read_real(inline)
        else if (parm == 'minnoise') then
          minnoise = read_real(inline)
        else if (parm == 'generate_noise') then
          generate_noise = read_logical(inline(first(3):first(3)))
        else if (parm == 'use_noise_file') then
          use_noise_file = read_logical(inline(first(3):first(3)))
        else if (parm == 'output_realspacenionweighted_values') then
          output_realspacenionweighted_values = read_logical(inline(first(3):first(3)))
        else  if (parm == 'output_realspacemassweighted_values') then
          output_realspacemassweighted_values = read_logical(inline(first(3):first(3)))
        else  if (parm == 'output_zspaceopticaldepthweighted_values') then
          output_zspaceopticaldepthweighted_values = read_logical(inline(first(3):first(3)))
        else if (parm == 'nlyman') then
          nLyman = read_int(inline)
          if (nLyman > nlyman_all) nLyman = nlyman_all
          if (nLyman < 0) nLyman = nlyman_all
        else  if (parm == 'modify_metallicity') then
           modify_metallicity = read_logical(inline(first(3):last(3)))
        else  if (parm == 'scale_simulation_abundances') then
           scale_simulation_abundances = read_logical(inline(first(3):last(3)))
        else  if (parm == 'z_rel') then
           z_rel = read_real(inline)
        else  if (parm == 'z_mean') then
           z_mean = read_real(inline)
        else  if (parm == 'impose_z_rho_relation') then
           impose_z_rho_relation = read_logical(inline(first(3):last(3)))
        else  if (parm == 'z_index') then
           z_index = read_real(inline)
        else  if (parm == 'maxz_rel') then
           maxz_rel = read_real(inline)
        else  if (parm == 'log_normal_scatter') then
           log_normal_scatter = read_logical(inline(first(3):last(3)))
        else  if (parm == 'z_sig_bin') then
           z_sig_bin = read_real(inline)
        else  if (parm == 'z_sig_dex') then
           z_sig_dex = read_real(inline)
        else if (parm == 'doH1') then
          doH1 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doHe2') then
          doHe2 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doC2') then
          doC2 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doC3') then
          doC3 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doC4') then
          doC4 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doC5') then
          doC5 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doC6') then
          doC6 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doN2') then
          doN2 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doN3') then
          doN3 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doN4') then
          doN4 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doN5') then
          doN5 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doN6') then
          doN6 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doN7') then
          doN7 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doO1') then
          doO1 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doO3') then
          doO3 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doO4') then
          doO4 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doO5') then
          doO5 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doO6') then
          doO6 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doO7') then
          doO7 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doO8') then
          doO8 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doMg2') then
          doMg2 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doAl2') then
          doAl2 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doAl3') then
          doAl3 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doNe8') then
          doNe8 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doNe9') then
          doNe9 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doSi2') then
          doSi2 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doSi3') then
          doSi3 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doSi4') then
          doSi4 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doS5') then
          doS5 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doFe2') then
          doFe2 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doFe3') then
          doFe3 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doFe17') then
          doFe17 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doFe19') then
          doFe19 = read_logical(inline(first(3):last(3)))
        else if (parm == 'doFe21') then
          doFe21 = read_logical(inline(first(3):last(3)))
        else if (parm == 'do21cm') then
          do21cm = read_logical(inline(first(3):last(3)))
        else if (parm == 'read_part_ids_from_file') then
          read_part_ids_from_file = read_logical(inline(first(3):last(3)))
        else  if (parm == 'flagged_particle_metallicity') then
          flagged_particle_metallicity = read_real(inline)
        else if (parm == 'doall') then
          doall = read_logical(inline(first(3):last(3)))
          !
          if(doall)then
            doH1   =.true.
            doHe2  =.true.
            doC2   =.true.
            doC3   =.true.
            doC4   =.true.
            doC5   =.true.
            doC6   =.true.
            doN2   =.true.
            doN3   =.true.
            doN4   =.true.
            doN5   =.true.
            doN6   =.true.
            doN7   =.true.
            doO1   =.true.
            doO3   =.true.
            doO4   =.true.
            doO5   =.true.
            doO6   =.true.
            doO7   =.true.
            doO8   =.true.
            doMg2  =.true.
            doAl2  =.true.
            doAl3  =.true.
            doNe8  =.true.
            doNe9  =.true.
            doSi2  =.true.
            doSi3  =.true.
            doSi4  =.true.
            doS5   =.true.
            doFe2  =.true.
            doFe3  =.true.
            doFe17 =.true.
            doFe19 =.true.
            doFe21 =.true.
            do21cm =.true.
          endif
          !
        else if (parm == 'ZC_rel') then
          ZC_rel = read_real(inline)
        else if (parm == 'ZN_rel') then
          ZN_rel = read_real(inline)
        else if (parm == 'ZO_rel') then
          ZO_rel = read_real(inline)
        else if (parm == 'ZNe_rel') then
          ZNe_rel = read_real(inline)
        else if (parm == 'ZMg_rel') then
          ZMg_rel = read_real(inline)
        else if (parm == 'ZAl_rel') then
          ZAl_rel = read_real(inline)
        else if (parm == 'ZSi_rel') then
          ZSi_rel = read_real(inline)
        else if (parm == 'ZS_rel') then
          ZS_rel = read_real(inline)
        else if (parm == 'ZFe_rel') then
          ZFe_rel = read_real(inline)
        else if (parm == 'gimic') then
          gimic = read_logical(inline(first(3):last(3)))
        else if (parm == 'urchin') then
          urchin = read_logical(inline(first(3):last(3)))
        else if (parm == 'do_periodic') then
          do_periodic = read_logical(inline(first(3):last(3)))
        else if (parm == 'snap') then
          snap = read_int(inline)
        else if (parm == 'urchindir') then
          urchindir = inline(first(3):last(3))
        else if (parm == 'subtract_Hmol') then
          subtract_Hmol = read_logical(inline(first(3):last(3)))
        else if (parm == 'use_gaussian_kernel') then
          use_gaussian_kernel = read_logical(inline(first(3):last(3)))
        else if (parm == 'integrate_kernel') then
          integrate_kernel = read_logical(inline(first(3):last(3)))
        else if (parm == 'wmap7') then
          wmap7 = read_logical(inline(first(3):last(3)))
        else if (parm == 'pixsize') then
          pixsize = read_real(inline)
        else if (parm == 'do_convolve_spectrum') then
          do_convolve_spectrum = read_logical(inline(first(3):last(3)))
        else if (parm == 'fwhm') then
          fwhm = read_real(inline)
        else if (parm == 'vpixsizekms') then
          vpixsizekms = read_real(inline)
        else if (parm == 'nlyman') then
          nlyman = read_int(inline)
        else if (parm == 'use_random_los') then
          use_random_los = read_logical(inline(first(3):last(3)))
        else if (parm == 'los_coordinates_file') then
          los_coordinates_file=inline(first(3):last(3))
        else if (parm == 'impose_eos') then
          impose_eos = read_logical(inline(first(3):first(3)))
        else if (parm == 'add_turbulence') then
          add_turbulence = read_logical(inline(first(3):first(3)))
        else if (parm == 'imposed_eos_T0') then
          imposed_eos_T0 = read_real(inline)
        else if (parm == 'imposed_eos_gamma') then
          imposed_eos_gamma = read_real(inline)
        else if (parm == 'imposed_eos_maxod') then
          imposed_eos_maxod = read_real(inline)
        else if (parm == 'limsigma') then
          limsigma = read_logical(inline(first(3):last(3)))
        else if (parm == 'integrate_thermprof_exactly') then
          integrate_thermprof_exactly = read_logical(inline(first(3):last(3)))
        else if (parm == 'use_urchin_temperature') then
          use_urchin_temperature = read_logical(inline(first(3):last(3)))
        else if (parm == 'use_smoothed_abundance') then
          use_smoothed_abundance = read_logical(inline(first(3):last(3)))
        else if (parm == 'ignore_starforming') then
          ignore_starforming = read_logical(inline(first(3):last(3)))
        else if (parm == 'setmaxt4sfgas') then
          setmaxt4sfgas = read_logical(inline(first(3):last(3)))
        else if (parm == 'ionfracone') then
          ionfracone = read_logical(inline(first(3):last(3)))
        else if (parm == 'verbose') then
          verbose = read_logical(inline(first(3):last(3)))
        else
          write (0,*) ' parameter not recognised: ',inline
          call abortrun(' error reading parameter file ')
        endif
      endif
    enddo
    !
    close(10)
    !
    return
  end subroutine read_parameter_file
  !
  subroutine perform_sanity_check() !Ensure that all compulsory variables have been set in the parameterfile
    use numbers
    use runtime
    use spectra
    use modified_metallicity
    use noisedata
    implicit none
    !
    if (.not. LLT(datadir,invalid) .and. .not. LGT(datadir,invalid)) &
      
call abortrun('Variable datadir is uninitialized in parameter file. stop')
    !
    if(do_long_spectrum) then
      if (.not. LLT(SpectrumFile,invalid) .and. .not. LGT(SpectrumFile,invalid)) &
        call abortrun('Variable SpectrumFile is uninitialized in parameter file. stop')
    endif
    !
    if (.not. LLT(file_list,invalid) .and. .not. LGT(file_list,invalid)) &
      call abortrun('Variable file_list is uninitialized in parameter file. stop')
    !
    if (.not. LLT(outputdir,invalid) .and. .not. LGT(outputdir,invalid)) &
      call abortrun('Variable outputdir is uninitialized in parameter file. stop')
    !    
    if (.not. LLT(ibdir,invalid) .and. .not. LGT(ibdir,invalid)) &
      call abortrun('Variable ibdir is uninitialized in parameter file. stop')
    !
    if (.not. use_random_los .and. use_snapshot_file) then
      if (.not. LLT(los_coordinates_file,invalid) .and. .not. LGT(los_coordinates_file,invalid)) &
        call abortrun('Variable los_coordinates_file is uninitialized in parameter file. stop')
    endif
    !
    if (nspec .eq. invalid_I) &
      call abortrun('Variable nspec is uninitialized in parameter file. stop')
    !
    if (nLyman .eq. invalid_I) &
      call abortrun('Variable nLyman is uninitialized in parameter file. stop')
    !    
    if (ibfactor .eq. invalid_R) &
      call abortrun('Variable ibfactor is uninitialized in parameter file. stop')
    !
    if (z_sig_dex .eq. invalid_R) &
      call abortrun('Variable z_sig_dex is uninitialized in parameter file. stop')
    !
    if (z_sig_bin .eq. invalid_I) &
      call abortrun('Variable z_sig_bin is uninitialized in parameter file. stop')
    !
    if (z_index .eq. invalid_R) &
      call abortrun('Variable z_index is uninitialized in parameter file. stop')
    !
    if (z_rel .eq. invalid_R) &
      call abortrun('Variable z_rel is uninitialized in parameter file. stop')
    !
    if (z_mean .eq. invalid_R) &
      call abortrun('Variable z_mean is uninitialized in parameter file. stop')
    !
    if (pixsize .eq. invalid_R) &
      call abortrun('Variable pixsize is uninitialized in parameter file. stop')
    !
    if (fwhm .eq. invalid_R) &
      call abortrun('Variable fwhm is uninitialized in parameter file. stop')
    !
    if (zqso .eq. invalid_R) &
      call abortrun('Variable zqso is uninitialized in parameter file. stop')
    !
    if (zabsmin .eq. invalid_R) &
      call abortrun('Variable zabsmin is uninitialized in parameter file. stop')
    !
    if (minbother_red .eq. invalid_R) &
      call abortrun('Variable minbother_red is uninitialized in parameter file. stop')
    !
    if (minbother_blue .eq. invalid_R) &
      call abortrun('Variable minbother_blue is uninitialized in parameter file. stop')
    !
    if (zabsmax .eq. invalid_R) &
      call abortrun('Variable zabsmax is uninitialized in parameter file. stop')
    !
    if (fzresol .eq. invalid_R) &
      call abortrun('Variable fzresol is uninitialized in parameter file. stop')
      if (minlambda .eq. invalid_R) &
        call abortrun('Variable minlambda is uninitialized in parameter file. stop')
    !
    if (maxlambda .eq. invalid_R) &
      call abortrun('Variable maxlambda is uninitialized in parameter file. stop')
    !
    if (vpixsizekms .eq. invalid_R) &
      call abortrun('Variable vpixsizekms is uninitialized in parameter file. stop')
    !
    if(generate_noise) then
      if (sigtonoise .eq. invalid_R) &
        call abortrun('Variable sigtonoise is uninitialized in parameter file. stop')
      !
      if (minnoise .eq. invalid_R) &
        call abortrun('Variable minnoise is uninitialized in parameter file. stop')
    endif
    !
    if(use_noise_file) then
      if (.not. LLT(noisefile,invalid) .and. .not. LGT(noisefile,invalid)) &
        call abortrun('Variable noisefile is uninitialized in parameter file. stop')      
    endif
    if (ZC_rel .eq. invalid_R) &
      call abortrun('Variable ZC_rel is uninitialized in parameter file. stop')
    !
    if (ZN_rel .eq. invalid_R) &
      call abortrun('Variable ZN_rel is uninitialized in parameter file. stop')
    !
    if (ZO_rel .eq. invalid_R) &
      call abortrun('Variable ZO_rel is uninitialized in parameter file. stop')
    !
    if (ZNe_rel .eq. invalid_R) &
      call abortrun('Variable ZNe_rel is uninitialized in parameter file. stop')
    !
    if (ZMg_rel .eq. invalid_R) &
      call abortrun('Variable ZMg_rel is uninitialized in parameter file. stop')
    !
    if (ZAl_rel .eq. invalid_R) &
      call abortrun('Variable ZAl_rel is uninitialized in parameter file. stop')
    !
    if (ZSi_rel .eq. invalid_R) &
      call abortrun('Variable ZSi_rel is uninitialized in parameter file. stop')
    !
    if (ZS_rel .eq. invalid_R) &
      call abortrun('Variable ZS_rel is uninitialized in parameter file. stop')
    !
    if (ZFe_rel .eq. invalid_R) &
      call abortrun('Variable ZFe_rel is uninitialized in parameter file. stop')
    !
    if(.not. do_long_spectrum .and. use_noise_file) &
      call abortrun('You cannot use a noise file for a short spectrum, change noise settings! stop')
    !
#ifdef EAGLE
    if(use_snapshot_file .and. (snap_base .eq. invalid)) &
      call abortrun('As use_snapshot_file = T then need to specify snapshpt base name  in parameter file. stop')
#else
    if(use_snapshot_file .and. (snap .eq. invalid_I)) &
      call abortrun('As use_snapshot_file = T then "snap = integer" must be present in parameter file. stop')
#endif
    !
    if (use_snapshot_file .and. do_long_spectrum) &
      call abortrun('Cannot do_long_spectrum with full snaposhot. stop')
    !
    if (.not. LLT(urchindir,invalid) .and. .not. LGT(urchindir,invalid)) then
      if (use_snapshot_file .and. urchin ) &
        call abortrun('Cannot do full snapshot, with urchin, when urchin_dir is not set in parameter file. stop')
    endif
    !
    if(use_urchin_temperature .and. .not. urchin) &
      call abortrun('Cannot use_urchin_temperature given urchin is false. stop')
    !
    if (impose_eos) then
      if (imposed_eos_T0 .eq. invalid_R) &
        call abortrun('impose_eos=T but imposed_eos_T0 is uninitialized in parameter file. stop')
      if (imposed_eos_gamma .eq. invalid_R) &
        call abortrun('impose_eos=T but imposed_eos_gamma is uninitialized in parameter file. stop')
      if (imposed_eos_maxod .eq. invalid_R) &
        call abortrun('impose_eos=T but imposed_eos_od is uninitialized in parameter file. stop')
    endif
    !
    if ((integrate_kernel).and. .not.(use_gaussian_kernel)) then
      call abortrun("integrate cublic spline kernel isn't working yet. stop")
    endif
    !
    if (modify_metallicity .and. read_part_ids_from_file) then 
      write(*,*)'if YOU ARE READING METALLICITIES FROM FILE (read_part_ids_from_file) '
      write (*,*) 'YOU SHOULD NOT ALSO BE USING modify_metallicity'
      write(*,*)'modify_metallicity will just overwrite this!  Fix parameters!'
      stop
    endif
    !
  end subroutine perform_sanity_check
  !
  function read_logical(char)
    use numbers
    implicit none
    character, intent(in) :: char
    logical read_logical
    !
    if(char == 'T')then
      read_logical = .true.
    else if (char == 'F') then
      read_logical = .false.
    else
      write (*,*) ' character = ',char
      call abortrun(' logicals need to be either T or F')
    endif
    return
  end function read_logical
  !
  function read_real(string)
    use numbers
    implicit none
    character(*), intent(in) :: string
    real(kind=doubleR)       :: read_real
    !
    character(len=200)       :: dum1, dum2
    integer                  :: io_stat
    !
    read(string,*, iostat=io_stat) dum1, dum2, read_real
    if(io_stat /= 0)then
      write (*,*) ' error reading ', string
      stop
    else
      !write (*,*) ' read real value: ',read_real
    endif
    return
  end function read_real
  !
  function read_int(string)
    use numbers
    implicit none
    character(*), intent(in) :: string
    integer(kind=singleI)    :: read_int
    !
    character(len=200)       :: dum1, dum2
    integer                  :: io_stat
    !
    read(string,*, iostat=io_stat) dum1, dum2, read_int
    if(io_stat /= 0)then
      write (*,*) ' error reading ', string
      stop
    else
      !write (*,*) ' read integer value: ',read_int
    endif
    return
  end function read_int
  !
  function read_char(string)
    use numbers
    implicit none
    character(*), intent(in) :: string
    character(len=200)       :: read_char
    !
    character(len=200)       :: dum1, dum2
    integer                  :: io_stat
    !
    read(string,*, iostat=io_stat) dum1, dum2, read_char
    write (*,*) ' dum1 = ',dum1
    write (*,*) ' dum2 = ',dum2
    write (*,*) ' read_char = ',read_char
    if(io_stat /= 0)then
      write (*,*) ' error reading ', string
      stop
    else
      !write (*,*) ' read string value: ',read_char
    endif
    return
  end function read_char
  !
  subroutine parm_syntax_error(string)
    implicit none
    character(*), intent(in) :: string
    character(len=300)       :: msg
    !
    write (*,*) ' string = ',string
    msg = ' error reading: '//trim(string)
    ! call abort(msg)
  end subroutine parm_syntax_error
  !
  !  Routine:     find_words()
  !
  !  Description: Given an input string, return the character positions of the
  !               beginnings and endings of all words in the string.  Also
  !               return the number of words found.
  !
  !               What constitutes a word:
  !
  !                 <delimiter>characters<delimiter>
  !                 symbol (single character)
  !                 <delimiter>characters<end-of-line>
  !                 <beg-of-line>characters<delimiter>
  !
  !               Symbols also function as delimiters.
  !
  subroutine find_words (string, first, last, n, nwords)
    implicit none
    !
    character(len=*)   :: string
    integer            :: n, nwords, first(n), last(n)
    !
    integer, parameter :: ndel = 2, nsym = 1
    !The deliminator characters must be a space and a TAB
    ! untabifying this won't compiler
    character    :: delimiters(ndel) = (/' ', '	'/)
    character    :: symbols(nsym)    = (/'='/)
    integer              :: i, j, k, strlen
    logical              :: is_delimiter, is_symbol
    !
    strlen = len(string)
    nwords = 0
    k      = 1
    !
    do i = 1, strlen
      !
      is_delimiter = .false.
      do j = 1, ndel
        if (string(i:i) == delimiters(j)) is_delimiter = .true.
      enddo
      is_symbol = .false.
      do j = 1, nsym
        if (string(i:i) == symbols(j)) is_symbol = .true.
      enddo
      !
      if (is_delimiter .or. is_symbol) then
        if (i > k) then
          nwords = nwords + 1
          first(nwords) = k
          last(nwords)  = i - 1
        endif
        k = i + 1
      endif
      if (is_symbol) then
        nwords = nwords + 1
        first(nwords) = i
        last(nwords)  = i
      endif
      !
    enddo
    !
    ! If we ended without a delimiter or symbol, record the last word.
    if (.not. (is_delimiter .or. is_symbol)) then
      nwords = nwords + 1
      first(nwords) = k
      last(nwords)  = strlen
    endif
    !
    return
  end subroutine find_words
  !
end module read_startup_file


module get_argumentsmod
  !
  ! this module provides the subroutine get_argument(iarg,x), which
  ! returns the i'th command line argument n variable x. x can be
  ! real, integer, double precision, string or logical.
  !
  implicit none
  !
  private
  public :: get_argument, get_nargs
  !
  integer :: iargc
  character(len=500) :: usage
  logical :: usage_set = .false.
  !
  interface get_argument
    module procedure get_integer_argument
    module procedure get_real_argument
    module procedure get_double_argument
    module procedure get_string_argument
    module procedure get_logical_argument
  end interface
  !
  contains
  !
  integer function get_nargs()
    integer :: iargc    
    get_nargs = iargc()
    return
  end function get_nargs
  !
  subroutine get_integer_argument(iarg, i)
    implicit none
    !
    integer, intent(in)  :: iarg
    integer, intent(out) :: i
    integer :: ios
    character(len=256)   :: str
    integer :: iargc
    !
    if(iargc().lt.iarg)then
      write(*,*)'error in get_argument(): too few arguments'
      write(*,*)'number of command line arguments is ',iargc()
      write(*,*)'attempted to read argument ',iarg
      if(usage_set)then
        write(*,*)''
        write(*,*)'usage: ',trim(usage)
        write(*,*)''
      end if
      stop
    end if
    !
    call getarg(iarg,str)
    read(str,*,iostat=ios)i

    if(ios.ne.0)then
       write(*,*)'error in get_argument():'
       write(*,*)'unable to interpret argument ',trim(str),' as integer'
       if(usage_set)then
          write(*,*)''
          write(*,*)'usage: ',trim(usage)
          write(*,*)''
       end if
       stop
    end if
    !
    return
  end subroutine get_integer_argument
  !
  subroutine get_real_argument(iarg, r)
    implicit none
    !
    integer, intent(in)  :: iarg
    real, intent(out) :: r
    integer :: ios
    character(len=256)   :: str
    integer :: iargc
    !
    if(iargc().lt.iarg)then
      write(*,*)'error in get_argument(): too few arguments'
      write(*,*)'number of command line arguments is ',iargc()
      write(*,*)'attempted to read argument ',iarg
      if(usage_set)then
        write(*,*)''
        write(*,*)'usage: ',trim(usage)
        write(*,*)''
      end if
      stop
    end if
    !
    call getarg(iarg,str)
    read(str,*,iostat=ios)r
    !
    if(ios.ne.0)then
      write(*,*)'error in get_argument():'
      write(*,*)'unable to interpret argument ',trim(str),' as real'
      if(usage_set)then
        write(*,*)''
        write(*,*)'usage: ',trim(usage)
        write(*,*)''
      end if
      stop
    end if
    !
    return
  end subroutine get_real_argument
  !
  subroutine get_double_argument(iarg, d)
    implicit none
    !
    integer, intent(in)  :: iarg
    double precision, intent(out) :: d
    integer :: ios
    character(len=256)   :: str
    integer :: iargc
    !
    if(iargc().lt.iarg)then
      write(*,*)'error in get_argument(): too few arguments'
      write(*,*)'number of command line arguments is ',iargc()
      write(*,*)'attempted to read argument ',iarg
      if(usage_set)then
        write(*,*)''
        write(*,*)'usage: ',trim(usage)
        write(*,*)''
      end if
      stop
    end if
    !
    call getarg(iarg,str)
    read(str,*,iostat=ios)d
    !
    if(ios.ne.0)then
      write(*,*)'error in get_argument():'
      write(*,*)'unable to interpret argument ',trim(str), &
        ' as double precision'
      if(usage_set)then
        write(*,*)''
        write(*,*)'usage: ',trim(usage)
        write(*,*)''
      end if
      stop
    end if
    !
    return
  end subroutine get_double_argument
  !
  subroutine get_logical_argument(iarg, l)
    implicit none
    !
    integer, intent(in)  :: iarg
    logical, intent(out) :: l
    integer :: ios
    character(len=256)   :: str
    integer :: iargc
    !
    if(iargc().lt.iarg)then
      write(*,*)'error in get_argument(): too few arguments'
      write(*,*)'number of command line arguments is ',iargc()
      write(*,*)'attempted to read argument ',iarg
      if(usage_set)then
        write(*,*)''
        write(*,*)'usage: ',trim(usage)
        write(*,*)''
      end if
      stop
    end if
    !
    call getarg(iarg,str)
    read(str,*,iostat=ios)l
    !
    if(ios.ne.0)then
      write(*,*)'error in get_argument():'
      write(*,*)'unable to interpret argument ',trim(str), &
        ' as logical'
      if(usage_set)then
        write(*,*)''
        write(*,*)'usage: ',trim(usage)
        write(*,*)''
      end if
      stop
    end if
    !
    return
  end subroutine get_logical_argument
  !
  subroutine get_string_argument(iarg, str)
    implicit none
    !
    integer, intent(in)  :: iarg
    character(len=*), intent(out) :: str
    ! cmb bugfix, we need to declare iargc
    ! external iargc
    integer :: iargc
    !
    if(iargc().lt.iarg)then
      write(*,*)'error in get_argument(): too few arguments'
      write(*,*)'number of command line arguments is ',iargc()
      write(*,*)'attempted to read argument ',iarg
      if(usage_set)then
        write(*,*)''
        write(*,*)'usage: ',trim(usage)
        write(*,*)''
      end if
      stop
    end if
    !
    call getarg(iarg,str)
    !
    return
  end subroutine get_string_argument
  !
end module get_argumentsmod


subroutine initialize_spectral_parameters
  use numbers
  use spectra
  use hdf5_wrapper
  use physical_constants
  use runtime
  use noisedata
  use modified_metallicity
  use header, only : redshift, los_this_file
  use my_mpi
  use projection_parameters
  implicit none
  !
  ! local variables
  logical eof, file_exists
  character(len=300) :: input_files, line
  character(len=300) :: basefile, longfile
  logical            :: single_file
  integer            :: ier, nf, i, j, file_handle
  character(len=3)   :: FileNumber
  !
  ! obtain list of all available redshifts in sightline files
  simzmin   = 1.e12
  simzmax   = -1.
  nlosfiles = 0
  !
  if (use_snapshot_file) then
    nlosfiles = 1
    allocate(simfile(nlosfiles))
    allocate(nlos_in_file(nlosfiles))
    allocate(simz(nlosfiles))
    allocate(ichoose(nlosfiles))
#ifdef EAGLE
    simfile(1) = trim(snap_base)//'.0.hdf5'
#else
    write (FileNumber,'(I3.3)') snap
    simfile(1) = '/snapshot_'//trim(FileNumber)//'/snap_'//trim(FileNumber)//'.0.hdf5'
#endif
    nlos_in_file(:) = 1
    !
    x_axis = 0 
    y_axis = 1 
    z_axis = 2
    !
    longfile = trim(datadir)//'/'//trim(simfile(1))
    call hdf5_open_file(file_handle, longfile, readonly=.true.)
    call read_header(file_handle)
    call check_header
    call hdf5_close_file(file_handle)
    !
    simz(1) = redshift
    simzmin = redshift
    simzmax = redshift
    !
  else
    !
    ! file file_list (set in par file) contains names of all available spectrum files 
    input_files = trim(file_list)
    open(unit=1,file=input_files,status='old',form='formatted',iostat=ier)
    !
    if (ier .ne. 0) then
       write(*,*) 'ERROR: Problem reading file ',trim(file_list)
       write(*,*) input_files
       call abortrun('IO error')
    endif
    !
    nf = -1
    do while(ier .eq. 0)
      read (1,'(a)',iostat=ier) line
      nf = nf + 1
    enddo
    close(1)
    nlosfiles = nf
    !
    allocate(simz(nlosfiles))
    allocate(simfile(nlosfiles))
    allocate(nlos_in_file(nlosfiles))
    allocate(ichoose(nlosfiles))
    !
    open(unit=1,file=input_files,status='old',form='formatted',iostat=ier)
    do i=1,nf
      read (1,'(a)',iostat=ier) line
      simfile(i) = line
      longfile = trim(datadir)//'/'//trim(simfile(i))
      inquire(file=longfile,exist=single_file)
      if(.not. single_file)then
        ! multi-file format
        basefile = trim(longfile)
        longfile = trim(basefile)//'.0.los.hdf5'
      endif
      call hdf5_open_file(file_handle, longfile, readonly=.true.)
      call read_header(file_handle)
      call check_header
      call hdf5_close_file(file_handle)
      simz(i) = redshift
      nlos_in_file(i) = los_this_file
      if(MyPE == 0) &
           write (*,*) ' z= ',redshift,' file= ',i
    enddo
    close(1)
    !
    simzmin = minval(simz)
    simzmax = maxval(simz)
    !
    if(verbose .and. MyPE == 0) then
      write(*,*)
      write(*,'("Redshift range available in sightline files    ", &
        f7.3," - ",f7.3)') simzmin, simzmax
    endif
    !
  endif
  !
  dodo_long_spectrum: if (do_long_spectrum) then 
    !
    if(MyPE == 0)then
      !
      if (vpixsizekms .gt. pixsize/maxlambda*LightSpeed/1e5) then
        write(*,'("ERROR: vpixsizekms > min desired pixel size!", &
          f6.3,f6.3)') vpixsizekms,  &
          pixsize/maxlambda*LightSpeed/1e5
        write(*,*) pixsize, maxlambda, LightSpeed
        write(*,'("Please decrease vpixsizekms in specwizard.par.")')
        call abortrun('stop')
      endif
    endif
    ! 
    !     ------------------
    !     Derived constants:
    !     ------------------
    !
    nvpix = log(maxlambda/minlambda)/(1e5*vpixsizekms)*LightSpeed
    !
    !  Make nvpix odd (convlv expects Gaussian to be of odd size, and Gaussian is filled up to nvpix
    if (mod(nvpix,2) .eq. 1) nvpix = nvpix + 1
    !
    !     For convolution with instrumental profile we need to Fourier 
    !     transform, we thus need to increase the array so that it is a
    !     power of 2.
    nvpix   = int(2**(aint(log(dble(nvpix))/log(2.)) + 1))
    ! number of pixels for the real space weighted arrays (function of z)
    nppix = log((1.d0 + zabsmax)/(1.d0 + zabsmin))/(1e5*vpixsizekms)*LightSpeed
    !
    call allocate_spectra_long()  ! size nvpix or nppix
    !
    vpixsize = vpixsizekms * 1e5 / LightSpeed
    if(MyPE == 0) then
      write(*,'(a,e9.2,a)')  &
        "Pixel size of long spectrum before binning: ", &
        vpixsizekms," km/s"
      write(*,'("WARNING: b-values should be much greater than this!")')
    endif
    do i=1, nvpix
      voverc(i) = dble(i-1) * vpixsize
    enddo
    !
    lambda = minlambda * exp(voverc) 
    !
    if (output_realspacenionweighted_values .or. output_realspacemassweighted_values) then
      ! write(*, '("Set up voverc_realspace: ", I8, " pixels, v/c = ", E12.4)') nppix, vpixsize
      do i=1, nppix
        voverc_realspace(i) = dble(i-1) * vpixsize
      enddo
      voverc_realspace = voverc_realspace + log(1.d0 + zabsmin)
      redshift_realspace = exp(voverc_realspace) - 1.d0
    endif
    !
    !     -------------------
    !     Read in noise file.
    !     -------------------
    ! 
    if (generate_noise) then 
      if (use_noise_file) then
        !     Load table size.
        call load_noise()
        !     Check limits.
        if (n_lambda(1) .gt. minlambda .and. MyPE == 0) then
          write(*,'("ERROR: Min. lambda of noise file too large!: ", &
            f10.3)') n_lambda(1)
          call abortrun('stop')
        endif
        if (n_lambda(n_nl) .lt. maxlambda .and. MyPE == 0) then
          write(*,'("ERROR: Max. lambda of noise file too small!: ", &
            f10.3)') n_lambda(n_nl)
          call abortrun('stop')
        endif
        !
        !     Need to allow for round-off errors.
        if (n_flux(1) .ge. 0.) n_flux(1) = -0.001
        if (n_flux(n_nf) .le. 1.) n_flux(n_nf) = 1.001
      else
        if (minnoise .ge. 1d0 / sigtonoise .and. MyPE == 0) then 
          write(*,*) 'ERROR: Must have minnoise < 1/sigtonoise'
          call abortrun('stop')
        endif
      endif
    endif
    !        -----------
    !        Print info.
    !        -----------
    if(MyPE == 0)then
      write(*,*)
      write(*,'("Creating continuous spectrum...")')
      write(*,*)
      write(*,'("QSO redshift:                ",f6.3)') zqso
      write(*, &
        '("Lambda:                     ",f7.1," - ",f7.1," A")')  &
        minlambda, maxlambda
      write(*,'("FWHM:                        ",f6.3," km/s")') fwhm
      write(*,'("Pixel size:                  ",f6.3," A")') pixsize
      write(*,'("Pre-instrument pixel size:   ",f6.3," km/s")')  &
        vpixsizekms
      write(*,'("Redshift bin size dz/(1+z):  ",f6.3)') fzresol
      write(*,*)
      write(*,'("Min. max. tau to bother with:")')
      write(*,'("Redward of Ly-alpha:         ",e8.1)') minbother_red
      write(*,'("Blueward of Ly-alpha:        ",e8.1)') minbother_blue
      write(*,*)
      write(*,'("UVB multiplication factor:   ",f6.3)') ibfactor
      write(*,*)
      write(*,'("Data dir:   ",a)') trim(datadir)
      write(*,'("Output dir: ",a)') trim(outputdir)
      write(*,*)
    endif
    !
    zstart = 1.e12
    zend = zqso
    do i = 1, nion
      do j = 1, nlines(i)
        zstart = min(zstart,minlambda/lambda_rest(i,j)-1.)
      enddo
    enddo
    zstart = max(zstart,0.)
    zstart = max(zstart,zabsmin)
    zend = min(zend,zabsmax)
    if (zstart .gt. zqso .and. MyPE == 0) &
      call abortrun('ERROR: zstart > zqso!')
    !
    if(MyPE == 0)then
      write(*,'("Redshift range needed for complete coverage:    ", &
        f7.3," - ",f7.3)') zstart, zend
      write(*,'("Redshift range used:                            ", &
        f7.3," - ",f7.3)') zstart, zend
      write(*,*)
    endif
    !
    ! allocate space for info about simulation files used
    allocate(los_used(max_nsimfile_used), icshift_used(max_nsimfile_used), &
      x_physical_used(max_nsimfile_used), y_physical_used(max_nsimfile_used), ibfactor_used(max_nsimfile_used), &
      losfile_used(max_nsimfile_used), x_axis_used(max_nsimfile_used), &
      y_axis_used(max_nsimfile_used), z_axis_used(max_nsimfile_used))
  else
    !
    if(verbose .and. MyPE == 0) then
      write (*,*) ' Program will compute spectra for all sightlines in list'
    endif
  endif dodo_long_spectrum
  !
end subroutine initialize_spectral_parameters

subroutine get_los_coordinates()
  use hdf5_wrapper
  use runtime
  use header
  use my_random_numbers
  use random_numbers, only : ran_los
  use spectra, only : nspec
  implicit none
  integer :: file_handle
  !
  integer :: i, nspec_in
  real(kind=doubleR), allocatable :: proj(:)
  !
  if (use_random_los) then 
    call init_random_numbers(777)
    !
    number_of_LOS = nspec
    allocate(x_fraction_array(number_of_LOS),y_fraction_array(number_of_LOS),z_fraction_array(number_of_LOS),&
      phi_array(number_of_LOS),theta_array(number_of_LOS),ncontribute(number_of_LOS))
    ncontribute = 0
    !
    do i = 1, number_of_LOS
      x_fraction_array(i) = random(ran_los)
      y_fraction_array(i) = random(ran_los)
      z_fraction_array(i) = random(ran_los)
      phi_array(i)     = 0
      theta_array(i)   = 0
    enddo
    !
  else
#ifdef AliVersion
     ! if use_random_los=F and use_snapshot_file=T then the following code reads in
     ! coordinates from los_coordinates_file
     IF (.NOT. use_random_los .AND. use_snapshot_file ) THEN
        WRITE(*,*) '### Reading the LOS coordinates from ', los_coordinates_file,' ###'
      	OPEN(77,file=los_coordinates_file)
     	READ(77,*) number_of_LOS
     	nspec_in = number_of_LOS

     	allocate(x_fraction_array(number_of_LOS),y_fraction_array(number_of_LOS),z_fraction_array(number_of_LOS),&
		phi_array(number_of_LOS),theta_array(number_of_LOS),ncontribute(number_of_LOS),ncontribute_global(number_of_LOS))

     	DO i = 1, number_of_LOS
           READ(77,*) x_fraction_array(i), y_fraction_array(i), z_fraction_array(i)
     	ENDDO
     	WRITE(*,*) '### enforcing the LOS to be along the z axis ###'	
     	z_fraction_array = 0

     	phi_array   = 0
    	theta_array = 0
     	ncontribute = 0 ! number of particles that contributed to spectrum

     	CLOSE(77)

	IF (nspec .NE. number_of_LOS) THEN
           WRITE(*,*)'Warning, overwriting nspec with the contents of the los_coordinates_file'
     	   nspec=number_of_LOS
	ENDIF
     ENDIF

#else

     ! to avoid trying inconsistent input, read directly from hdf5 file.
     ! since current version only works for (x,y) projection // z, only allow that as input
     call hdf5_open_file(file_handle, los_coordinates_file, readonly=.true.)
     call hdf5_read_attribute(file_handle,'/Projection/nspec', number_of_LOS)
     nspec_in = number_of_LOS
     if (nspec .gt. number_of_LOS) then
        write(*,*)'Warning, overwriting nspec with the contents of the los_coordinates_file'
        nspec = number_of_LOS
     else
        number_of_LOS = nspec
     endif
     allocate(x_fraction_array(number_of_LOS),y_fraction_array(number_of_LOS),z_fraction_array(number_of_LOS),&
          phi_array(number_of_LOS),theta_array(number_of_LOS),ncontribute(number_of_LOS),ncontribute_global(number_of_LOS))
     
     ! read desired projection coordinates and copy
     allocate(proj(nspec_in))
     call hdf5_read_data(file_handle, 'Projection/x_fraction_array', proj)
     x_fraction_array = proj(1:nspec)
     call hdf5_read_data(file_handle, 'Projection/y_fraction_array', proj)
     y_fraction_array = proj(1:nspec)
     deallocate(proj)
     call hdf5_close_file(file_handle)
     z_fraction_array = 0
     phi_array   = 0
     theta_array = 0
     ncontribute = 0 ! number of particles that contributed to spectrum
    !  ! current version of specwizard now only works for projections along z-axis - so only allow those to be specified
    ! open(77,file=los_coordinates_file)
    ! read(77,*) ! comment
    ! read(77,*) number_of_LOS
    ! read(77,*) ! comment
    ! !
    ! if (nspec .gt. number_of_LOS) then
    !   write(*,*)'Warning, overwriting nspec with the contents of the los_coordinates_file'
    !   nspec = number_of_LOS
    ! else
    !   number_of_LOS = nspec
    ! endif
    ! !
    ! allocate(x_fraction_array(number_of_LOS),y_fraction_array(number_of_LOS),z_fraction_array(number_of_LOS),&
    !   phi_array(number_of_LOS),theta_array(number_of_LOS))
    ! !
    ! do i = 1, number_of_LOS                  
    !    !!      read(77,*) x_fraction_array(i), y_fraction_array(i), z_fraction_array(i), phi_array(i), theta_array(i)
    !    read(77,*) x_fraction_array(i), y_fraction_array(i)
    ! enddo
    ! z_fraction  = 0
    ! phi_array   = 0
    ! theta_array = 0
    ! !
    ! close(77)
    !
#endif 
!AliVersion

  endif
  !
  if(ANY(theta_array .ne. 0) .or. ANY(phi_array .ne. 0)) & !!!
    call abortrun('Rotation has been turned off. Non zero theta/phi found. stop') !!!
  !
end subroutine get_los_coordinates

subroutine create_spectrum_file(ifile)
  use numbers
  use hdf5_wrapper
  use runtime
  use spectra, only : simfile
  use my_mpi
  implicit none
  integer, intent(in) :: ifile
  !
  integer(kind=singleI)    :: outfile_handle, infile_handle
  character(len=200)       :: inputfile
  character(len=3)         :: FileNumber
  logical                  :: file_exists, single_file
  !
  ! generate name of output file
  if(do_long_spectrum) then
     if (.not. LLT(SpectrumFile,invalid) .and. .not. LGT(SpectrumFile,invalid)) then
        if (.not. use_snapshot_file) then
           SpectrumFile = trim(adjustl(outputdir))//'/Spectrum.'//trim(adjustl(simfile(ifile)))
        else
           write (FileNumber,'(I3.3)') snap
           SpectrumFile = 'snap_'//trim(FileNumber)//'.hdf5'
           SpectrumFile = trim(adjustl(outputdir))//'/Spectrum.'//trim(adjustl(SpectrumFile))
        endif
     else
        SpectrumFile = trim(adjustl(outputdir))//SpectrumFile
     endif
  endif
  !
  if(MyPE .eq. 0) then
    ! inquire if output file exists, otherwise stop
    inquire(file=trim(SpectrumFile),exist=file_exists)
    if(file_exists .and. .not. overwrite)then
      write (*,*) 'The required output file ',trim(adjustl(SpectrumFile)),' already exists'
      write (*,*) 'Please delete this file if you want me to continue'
      call abortrun('stop')
    endif
    !
    ! create output file
    call hdf5_create_file(outfile_handle, trim(SpectrumFile))
    !
  endif
  !
  ! find and open the first input file 
  if (.not. use_snapshot_file) then
    inputfile = trim(datadir)//'/'//trim(simfile(ifile))
  else
#ifdef EAGLE
     inputfile = trim(datadir)//trim(snap_base)//'.hdf5'
#else
    write (FileNumber,'(I3.3)') snap
    inputfile = trim(datadir)//'/snapshot_'//trim(FileNumber)//'/snap_'//trim(FileNumber)//'.hdf5'
#endif
    inquire(file=inputfile,exist=single_file)
    if(.not. single_file)then
#ifdef EAGLE
       inputfile = trim(datadir)//trim(snap_base)//'.0.hdf5'
#else
       inputfile = trim(datadir)//'/snapshot_'//trim(FileNumber)//'/snap_'//trim(FileNumber)//'.0.hdf5'
#endif
    endif
  endif
  !
  call hdf5_open_file(infile_handle,inputfile,readonly=.true.)
  !
  call read_header(infile_handle)
  call check_header()
  call read_units(infile_handle)
  call read_constants(infile_handle)
  call read_parameters(infile_handle)
  call hdf5_close_file(infile_handle)
  !
  ! write attribute groups to output file
  if(MyPE .eq. 0) then
    call write_header(outfile_handle)
    call write_units(outfile_handle)
    call write_constants(outfile_handle)
    call write_parameters(outfile_handle)
    call write_specwizard_runtime_parameters(outfile_handle)
    call hdf5_close_file(outfile_handle)
  endif
  !
end subroutine create_spectrum_file


! +++++++++++++++++++++++++++++++++++++++ ...startup +++++++++++++++++++++++++++++++++++++++++++++ !

! ++++++++++++++++++++++++++++++++++++++++++ spectra... ++++++++++++++++++++++++++++++++++++++++++ !


subroutine zero_spectra()
  use numbers
  use spectra
  implicit none
  !
  tau_long             = 0.d0
  tau_long_strongest   = 0.d0
  temp_z_ion_long      = 0.d0
  rho_z_ion_long       = 0.d0
  veloc_z_ion_long     = 0.d0
  temp_ion_long        = 0.d0
  rho_ion_long         = 0.d0
  n_ion_long           = 0.d0
  temp_long            = 0.d0
  rho_long             = 0.d0
  met_long             = 0.d0
  veloc_long           = 0.d0
  cdens_ion_integrated = 0.d0
  !
end subroutine zero_spectra


subroutine projectdata()
  ! Computes SPH estimates of the density, temperature and line of
  ! sight velocities. Then calls computespectrum to compute spectra.
  use numbers
  use header
  use particledata
  use projection_parameters
  use spectra
  use runtime
  use atomic_data
  use physical_constants
  use ionization_tables
  use modified_metallicity, only : modify_metallicity
  use w4_gadget_spline_kernel_class
  implicit none 
  !
  ! local variables
  integer(kind=singleI) :: i, ioff, iiz, j, iz, ii, ion,iz1, iz2, h1_index=-1, si2_index=-1
  real(kind=doubleR)    :: xx, yy, zz, hh, h2, b, b2, hinv2, hinv3, vr, zmingrid, zmaxgrid, dzinv, dzgrid, box, box_2
  real(kind=doubleR)    :: dzmax, zgrid, zf, zi, zedge, ztrans, dr2, qi, qf, deltaz, dvbin, dzbin, dz1, dz2
  real(kind=doubleR)    :: Density, log_dens, log_temp, DensCon, impactparameter, q, z, kernel_factor
  real(kind=doubleR)    :: RotationMatrix(3,3), dx, dy, Q1, Q2
  !
  ! read urchin flags
  do ion=1, nion
    if(trim(ions(ion)) .eq. 'h1') h1_index = ion
    if(trim(ions(ion)) .eq. 'si2') si2_index = ion
  enddo
  !
  ! Initialize arrays
  call allocate_spectra_short()
  !
  ! density-weighted real-space quantities
  rho_tot   = 0.0
  temp_tot  = 0.0
  met_tot   = 0.0
  veloc_tot = 0.0
  if(.not. do_long_spectrum) cdens_ion_integrated = 0.0
  !
  ! real-space quantities
  n_ion     = 0.0  ! number density
  temp_ion  = 0.0  ! number dens weighted temperature
  rho_ion   = 0.0  ! number dens weighted density
  veloc_ion = 0.0  ! number dens weighted velocity
  !
  ! Redshift interpolation for ionizing background
  if (nz .gt. 1) then 
    if (.not. use_maxdens_above_zmax) then
      if (zcurrent .lt. ib_redshift(1) .or. zcurrent .gt. ib_redshift(nz))then
        write(*,*) 'ERROR: z out of bounds ioniz. bal. table!'
        write(*,'("z = ",f8.5,", zmin = ",f8.5,", zmax = " &
          ,f8.5)') zcurrent, ib_redshift(1), ib_redshift(nz)
        stop
      endif
    endif
    !
    ! If z is greater than the maximim redshift then force that we use the
    ! final redshift entry of the table:
    if (use_maxdens_above_zmax .and. zcurrent .gt. ib_redshift(nz)) then
      iz1 = nz-1
      iz2 = nz
      dz1 = 0.0
      dz2 = 1.0
    else
      iz2 = 2
      do while (zcurrent .gt. ib_redshift(iz2))
        iz2 = iz2 + 1
      enddo
      iz1 = iz2-1
      dz1 = (ib_redshift(iz2) - zcurrent) / (ib_redshift(iz2) - ib_redshift(iz1))
      dz2 = 1. - dz1
    endif
  else
    write(*,'("WARNING: Ioniz. bal. table contains only", &
      " 1 redshift! (z = ",f6.3,")")') ib_redshift(1)
    iz1 = 1
    iz2 = iz1
    dz1 = 0.5
    dz2 = 1. - dz1
  endif
  !
  zmingrid  = 0.0
  zmaxgrid  = BoxSize * ExpansionFactor / HubbleParam ! size of simulation box in physical Mpc
  box       = zmaxgrid
  dzgrid    = (zmaxgrid-zmingrid) / dble(nveloc)
  dzinv     = 1. / dzgrid
  box_2     = 0.5 * box
  densscale = (ExpansionFactor/acurrent)**3 ! scale particle densities from snapshot redshift to current redshift
  !
  ! Rotate and shift particle positions relative to sightline
  !!!call SetRotationMatrix(theta_projection,phi_projection,RotationMatrix)
  !
  !!!ShiftedPosition(1,:) = Position(1,:) - x_physical
  !!!ShiftedPosition(2,:) = Position(2,:) - y_physical
  !!!ShiftedPosition(3,:) = Position(3,:) - z_physical
  !!!ShiftedPosition      = matmul(RotationMatrix,ShiftedPosition)
  !!!ShiftedVelocity      = matmul(RotationMatrix,ShiftedVelocity)
  !
  ! Periodic boundary conditions - careful, creates a truncated cube!
  !!!if(do_periodic) then
  !!!  do i=1,NGas
  !!!    do j=1,3
  !!!      if(ShiftedPosition(j,i) .gt. box_2)  then
  !!!        ShiftedPosition(j,i) = ShiftedPosition(j,i) - box
  !!!      else if(ShiftedPosition(j,i) .lt. -box_2) then
  !!!        ShiftedPosition(j,i) = ShiftedPosition(j,i) + box
  !!!      endif
  !!!    enddo
  !!!  enddo
  !!!endif
  !
    ! write (1,*) ' projection: ',x_physical, y_physical
    ! do i=1, 20
    !    write (1,100) PartID(i), Position(1,i), Position(2,i), Position(3,i)
    ! enddo
    ! write (1,*)
  ! Loop over particles in sightline.
  ncontr = 0 ! number of particles that contributes to sightline
  particle_loop: do i = 1, NGas
    !
    !!!xx = ShiftedPosition(1,i)
    !!!yy = ShiftedPosition(2,i)
    !!!zz = ShiftedPosition(3,i)
    xx = Position(1,i) !!!
    yy = Position(2,i) !!!
    zz = Position(3,i) !!!
    !
    hh = ParticleSmoothingLength(i) ! CAUTION hh is the zero point of the smoothing kernal
    h2 = hh*hh
    hinv2 = 1. / h2
    hinv3 = hinv2 / hh
    !
    dx = abs(xx - x_physical) 
    dy = abs(yy - y_physical) 
    if(dx .gt. box_2) then 
      dx = box - dx 
    endif 
    if(dy .gt. box_2) then 
      dy = box - dy 
    endif 
    b2 = dx**2 + dy**2 
    b  = sqrt(b2)
    impactparameter = b
    !
    if(impactparameter .le. hh) then
       ncontr = ncontr + 1
!        write (1,100) PartID(i), xx, yy, zz, hh
! 100 format(I16,1x,4(e12.4,1x))
       !
        ! Determine ionfrac :-
        ! Ionization Fraction         =  X_ion                     =  N_ion / N_element
        ! ParticleNeutralHFraction    =  X_H1 + X_Hmol             =  (N_H1 + 2*N_Hmol) / N_H
        ! ParticleMolecularHFraction  =  X_Hmol / (X_H1 + X_Hmol)  =  2*N_Hmol / (N_H1 + 2*N_Hmol)
        ! 
        ! Therefore:
        ! 
        ! X_H1    =  N_H1 / N_H    =  ParticleNeutralHFraction * (1 - ParticleMolecularHFraction)
        ! X_Hmol  =  N_Hmol / N_H  =  ParticleNeutralHFraction * ParticleMolecularHFraction
        ! 
        ! subtract_Hmol = F  ->  X_Si2 = X_H1 + X_Hmol
        ! subtract_Hmol = T  ->  X_Si2 = X_H1
        !
        ! Total hydrogen number density, n_H (per cc)
        Density = ParticleDensity(i)
        if (ParticleDensity(i) .lt. 0.0) then
          write(*,*)'This is not good!  Density is less than 0: ',ParticleDensity(i)
          stop
        endif
        ! Rescale density from sim redshift to current z using densscale.
        ! Note the scaling below only affects the ionization balance (it is
        ! actually not quite right if the internal ionization was used),
        ! actual densities are scaled at the end.
        Density = Density * densscale
        log_dens = log10(Density)
        !
        ! Temperature (K)
        log_temp = log10(ParticleTemperature(i)) ! T in Kelvin
        !
        call computeib(ib_redshift(iz1),ib_redshift(iz2),iz1,iz2,dz1,dz2,log_temp,log_dens,ionfrac)
        !
        ! Calculate/override ionfrac for h1 and si2 using urchin values
        if(urchin) then
          if(subtract_Hmol) then
            if(h1_index .gt. 0) &
              ionfrac(h1_index)  = ParticleNeutralHFraction(i) * (1.0d0 - ParticleMolecularHFraction(i) )
            if(si2_index .gt. 0) &
              ionfrac(si2_index) = ParticleNeutralHFraction(i) * (1.0d0 - ParticleMolecularHFraction(i) )
          else
            if(h1_index .gt. 0) &
              ionfrac(h1_index)  = ParticleNeutralHFraction(i)
            if(si2_index .gt. 0) &
              ionfrac(si2_index) = ParticleNeutralHFraction(i)
          endif
        endif
        !
        if(ionfracone) then
          totnr_ion(:) = MassFractions(ion_elnr(:),i) * Mass(i) / ElementAtomicMass(ion_elnr(:)) ! [Msun/g]
        else
          totnr_ion(:) = ionfrac(:) * MassFractions(ion_elnr(:),i) * Mass(i) / ElementAtomicMass(ion_elnr(:)) ! [Msun/g]
        endif
        !
        ! z velocity (km/s)
        if(NoPecVel) then
          vr = 0.0
        else
          !!!vr = ShiftedVelocity(3,i) ! peculiar velocity in km/s
          vr = Velocity(3,i) ! peculiar velocity in km/s !!!
        endif
        !
        ! central pixel to contribute to
        iz = (zz - zmingrid) * dzinv + 1
        !
        if(gimic)then
         if(boundary(i) == 1) spectrum_boundary(iz)=spectrum_boundary(iz)+1
        endif
        !
        ! contribute to projection segment
        dzmax = sqrt(abs(h2 - b2))
        ioff = int(dzmax * dzinv) + 2
        if (ioff*2 .gt. nveloc) stop 'ioff*2 > nveloc'
        !
        ! segment pixel loop
        do iiz = iz-ioff, iz+ioff+1
          j = iiz
          j = mod(j-1+10*nveloc,nveloc)+1
          zgrid = zmingrid + (dble(j)-0.5)*dzgrid
          deltaz = abs(zgrid - zz)
          if (deltaz .gt. box_2) deltaz = box - deltaz
          dr2 = b2 + deltaz**2
          zf = deltaz + dzgrid*0.5
          zi = deltaz - dzgrid*0.5
          !
          if(integrate_kernel) then
            if(use_gaussian_kernel) then
              ! integrate normalized truncated gaussian kernel
              kernel_factor = integrate_G3_line(hh,impactparameter,zi,zf) / dzgrid
            else
              !         integrate cubic spline M4         !
              !!!   THIS INTEGRATION ISN'T WORKING YET  !!!
              !                                           !
              zedge = sqrt(h2-b2)
              ! check for (non) intersection anywhere within pixel
              if((zf*zi .ge. 0.0) .and. (abs(zi) .ge. zedge) .and. (abs(zf) .ge. zedge)) then
                kernel_factor = 0.d0
              else
                ! limit integration to within to kernel
                ! know that zf > zi
                zf = sign(min(abs(zf),abs(zedge)),zf)
                zi = sign(min(abs(zi),abs(zedge)),zi)
                ! transition between functions
                ztrans = h2*0.25 - b2
                kernel_factor = 0.d0
                if ((abs(zf).le.ztrans) .and. (abs(zi).le.ztrans)) then
                  ! limits both in inner kernel
                  kernel_factor = Q1(zf,b,hh) - Q1(zi,b,hh)
                else if ((abs(zf).ge.ztrans) .and. (abs(zi).ge.ztrans)) then
                  ! limits both in outer kernel
                  if(zf*zi .ge. 0.0) then
                    ! just outer kernel
                    kernel_factor = Q2(zf,b,hh) - Q2(zi,b,hh)
                  else
                    ! limits straddle z=0
                    kernel_factor = Q2(zf,b,hh) - Q2(zi,b,hh)
                    kernel_factor = kernel_factor - ( Q2(ztrans,b,hh) - Q2(-ztrans,b,hh) )
                    kernel_factor = kernel_factor + ( Q1(ztrans,b,hh) - Q1(-ztrans,b,hh) )
                  endif
                else 
                  ! one inner one outer kernel limit
                  if (zf .ge. ztrans) then
                    kernel_factor = Q2(zf,b,hh) - Q2(ztrans,b,hh)
                    kernel_factor = kernel_factor + Q1(ztrans,b,hh) - Q1(zi,b,hh)
                  else if (zf .le. -ztrans) then
                    kernel_factor = Q2(zf,b,hh) - Q2(-ztrans,b,hh)
                    kernel_factor = kernel_factor + Q1(-ztrans,b,hh) - Q1(zi,b,hh)
                  else if (zi .ge. ztrans) then
                    kernel_factor = Q2(zi,b,hh) - Q2(ztrans,b,hh)
                    kernel_factor = kernel_factor + Q1(ztrans,b,hh) - Q1(zf,b,hh)
                  else if (zi .le. -ztrans) then
                    kernel_factor = Q2(zf,b,hh) - Q2(-ztrans,b,hh)
                    kernel_factor = kernel_factor + Q1(-ztrans,b,hh) - Q1(zf,b,hh)
                  endif
                endif
              endif
              kernel_factor = (kernel_factor * 8. * hinv3) / ( pi * dzgrid)
            endif
            !!!     ABOVE INTEGRATION ISN'T WORKING YET    !!!
          else
            if(use_gaussian_kernel) then
              ! normalized truncated gaussian kernel, constant over pixel
              kernel_factor = gaussian_kernel(sqrt(dr2), hh, 3)
            else
              ! cubic spline, constant over pixel
              q = sqrt(dr2 * hinv2) ! q = r/h
              if (q .le. .5) then
                kernel_factor = (1.+6.*q**2*(q-1.))
              else if (q .le. 1.) then
                kernel_factor = 2.*(1.-q)**3
              else
                kernel_factor = 0.d0
              endif
              kernel_factor = (kernel_factor * 8. * hinv3) / pi
            endif
          endif
          !
          if(gimic)then
            if(boundary(i) == 1) kernel_factor = 0.d0
          endif
          !
          if(kernel_factor .gt. 0.) then
            !
            ! quantities weighted by number of ions.  These are actually
            ! SPH estimates of A*rho, but we divide through by rho on the
            ! next loop.  This ensures that we explicitly conserve mass.
            n_ion(:,j)     = n_ion(:,j)     + kernel_factor * totnr_ion(:)
            veloc_ion(:,j) = veloc_ion(:,j) + kernel_factor * totnr_ion(:) * vr 
            temp_ion(:,j)  = temp_ion(:,j)  + kernel_factor * totnr_ion(:) * ParticleTemperature(i)
            ! Particle Density n_H = cgs, rescaled for long spectra -> Particle cgs rho
            rho_ion(:,j)   = rho_ion(:,j)   + kernel_factor * totnr_ion(:) * Density * proton_mass / MassFractions(H_index,i) 
            ! .... weighted by mass
            rho_tot(j)     = rho_tot(j)     + kernel_factor * Mass(i)
            veloc_tot(j)   = veloc_tot(j)   + kernel_factor * Mass(i) * vr
            temp_tot(j)    = temp_tot(j)    + kernel_factor * Mass(i) * ParticleTemperature(i)
            met_tot(j)     = met_tot(j)     + kernel_factor * Mass(i) * Metallicity(i)
            !
          endif ! kernel factor > 0
        enddo ! loop over contributing vertices
        !      
    endif ! b le hh
    !
  enddo particle_loop
  !
  ! Mass was computed in M_sun (was used to compute particle nr),
  ! distance was in proper Mpc, conversion factor to n (cm^-3) is
  ! thus: 
  DensCon = Msun / Mpc**3 
  ! Rescale density from sim redshift to current z using densscale.
  DensCon = DensCon * densscale
  !
  do ii = 1, nion
    do i = 1, nveloc
      if (n_ion(ii,i) .gt. 0.) then 
        veloc_ion(ii,i) = veloc_ion(ii,i) / n_ion(ii,i)
        temp_ion(ii,i)  = temp_ion(ii,i) / n_ion(ii,i)
        ! rho_ion was already cgs -> overdensity (right scaling per output time for long spectra)
        ! rhocb is at snapshot expansion factor, but densities are rescaled to acurrent 
        ! -> rescale normalizing rhocb same way
        rho_ion(ii,i)   = rho_ion(ii,i) / n_ion(ii,i) / (rhocb * densscale)
      endif
      n_ion(ii,i) = n_ion(ii,i) * DensCon ! ions/cm^3
    enddo
  enddo
  !
  do i = 1, nveloc
    if (rho_tot(i) .gt. 0.) then
      veloc_tot(i) = veloc_tot(i) / rho_tot(i)
      temp_tot(i)  = temp_tot(i)  / rho_tot(i)
      met_tot(i)   = met_tot(i)   / rho_tot(i)
      rho_tot(i)   = rho_tot(i) * DensCon / rhocb ! g/cm^3, -> overdensity (right scaling per output time for long spectra)
    endif
  enddo
  !
  return
end subroutine projectdata


function Q1(z,b,h)
  use numbers
  implicit none
  real(kind=doubleR), intent(in) :: z,b,h
  real(kind=doubleR) :: q, b_h, z_h
  real(kind=doubleR) :: Q1 
  !
  q = sqrt(z**2+b**2) / h
  b_h  = b / h
  z_h  = z / h
  !
  Q1 = z_h*( 4.0*b_h**2 + q*(2.0*q-3.0) ) + (1.0 - 3.0*b_h**2)*log(z_h+q)
  !
end function


function Q2(z,b,h)
  use numbers
  implicit none
  real(kind=doubleR), intent(in) :: z,b,h
  real(kind=doubleR) :: q, b_h, z_h
  real(kind=doubleR) :: Q2
  !
  q = sqrt(z**2+b**2) / h
  b_h  = b / h
  z_h  = z / h
  !
  Q2 = ( (9*b_h**2+6.0)*log(z_h+q) - z_h*(4.0*b_h**2 +2.0*q**2 - 9.0*q + 18.0) )/3.0
  !
end function


subroutine makespectra()
  use numbers
  use spectra
  use runtime
  use physical_constants
  use header, only : BoxSize, HubbleParam, Omega0, OmegaLambda
  use my_mpi
  use particledata, only: rhocb
  use parameters, only: SF_EOSGammaEffective, SF_THRESH_MinPhysDens_HpCM3, &
    SF_EOSEnergyAtThreshold_ERG, InitAbundance_Hydrogen
  use constants, only : GAMMA
  use atomic_data, only: massH
  implicit none
  !
  ! local variables
  integer(kind=singleI) :: i, ion
  real(kind=doubleR)    :: CurrentHubbleCt, totcdens, boxkms, dvbin, dzbin
  real(kind=doubleR)    :: GammaEOSm1, Rho0, TURBFACT, RhoH
  !
  ! Should compute vhubble as follows:
  ! dz/dr = H(z)r / c, where r is comoving radial coordinate.
  ! Integrate above to get z for each nveloc (note: up in z?)
  ! vhubble/c = ln((1+z)/(1+z_current)), where z_current is beginning of box
  ! Scale density here: dens = dens * ((1.+z)/(1.+zcurrent))^3
  ! Note that none of this is correct for fully collapsed objects
  ! which are constant in proper coordinates.
  !
  ! Compute Hubble parameter (in km/s/Mpc).
  CurrentHubbleCt = 100. * HubbleParam *  &
    sqrt(1. + omega0*(1./acurrent-1.) + OmegaLambda* &
    (acurrent**2-1.)) /acurrent
  !
  boxkms = BoxSize / HubbleParam * acurrent * CurrentHubbleCt
  dvbin = boxkms / dble(nveloc)   ! bin size (km/s)
  dzbin = BoxSize / HubbleParam * acurrent * Mpc / dble(nveloc) ! bin size (physical cm)      

#ifdef HUBBLE
  boxkms = boxkms * HUBBLE
  dvbin  = boxkms / dble(nveloc)   ! bin size (km/s)
  dzbin  = dzbin * HUBBLE
#endif

  !
  do i = 1, nveloc
    vhubble(i) = dble(i-1) * dvbin
  enddo
  !
  turbulence = 0.d0 ! real-space turbulent broadening [(cm/s)^2]
  if(add_turbulence) then
     !
     ! We add turbulent broadening to thermal line widths (b-parameters, b_T for thermal) as
     ! b^2 -> b_T^2 + 2sigma^2, where 2sigma^2 = 2 * (5/3-1) * u_0 * [(rho/rho_0)^gamma_EOS-1 -1]
     ! The imposed EOS is P = P_0 * (rho/rho_0)^gamma_EOS
     ! Note that the dimension of SF_EOSEnergyAtThreshold is (cm/s)^2, not ergs
     !
     GammaEOSm1      = SF_EOSGammaEffective - 1.0 ! gamma_EOS -1 [1]
     Rho0            = SF_THRESH_MinPhysDens_HpCM3 * massH !/ InitAbundance_Hydrogen ! rho0 [g cm^-3]
     TURBFACT        = 2.0 * (GAMMA-1.0) * SF_EOSEnergyAtThreshold_ERG ! 2sigma^2 = 2 * (5/3-1) * u_0; [(cm/s)^2]
     do i = 1,nveloc
        rhoH =  rho_tot(i) * (1.d0 - met_tot(i)) ! [g cm^-3]
        if (rhoH .ge. Rho0) then
           turbulence(i) = TURBFACT * ((rhoH / Rho0)**GammaEOSm1 - 1.0) ! [ (cm/s)^2]
        endif
     enddo
     !
  endif
  !
  do ion = 1, nion
    !
    cdens(:)    = n_ion(ion,:) * dzbin ! column density in ion / cm^2
    totcdens    = sum(cdens)
    !
    if (verbose .and. MyPE == 0) then
      write(*,'("N_",a,",tot   = ",es12.3)')  &
        trim(adjustl(ions(ion))), totcdens
    endif
    !
    cdens_ion_integrated(ion) = cdens_ion_integrated(ion) + totcdens
    !
    call computespectrum(nveloc,ion,vhubble,rho_tot,cdens,veloc_ion(ion,:),n_ion(ion,:), &
        temp_ion(ion,:),turbulence,tau,velocw,nionw,rhow,tempw)
    !
    tau_ion(ion,:)     = tau(:)
    flux_ion(ion,:)    = exp(-tau(:))
    veloc_z_ion(ion,:) = velocw(:)
    nion_z_ion(ion,:)  = nionw(:)
    rho_z_ion(ion,:)   = rhow(:) ! input density already in units of mean baryon density
    temp_z_ion(ion,:)  = tempw(:)
    !
  enddo
  !
end subroutine makespectra


subroutine computespectrum(nveloc,ion,vhubble,rho_tot,cdens,vpecul,nion,temperature,turbulence,&
     tau,velocw,nionw,rhow,tempw)
  use numbers
  use atomic_data
  use spectra, only : lambda_rest, fosc, ion_mass, integrate_thermprof_exactly, &
       minbother_red, minbother_blue, vpixsizekms, limsigma, normw, ions
  use physical_constants
  use runtime
  implicit none
  !
  integer(kind=singleI), intent(in) :: nveloc, ion
  real(kind=doubleR), intent(in)    :: vhubble(nveloc),rho_tot(nveloc),cdens(nveloc), &
       vpecul(nveloc),nion(nveloc),temperature(nveloc),turbulence(nveloc) 
  real(kind=doubleR), intent(out)   :: tau(nveloc), rhow(nveloc), tempw(nveloc), &
       velocw(nveloc), nionw(nveloc)
  !
  ! local variabales  
  real(kind=doubleR) :: lambda0, fvalue,mass
  !
  !
  ! vhubble      ! Hubble velocity (km/s)
  ! rho_tot      ! total baryonic density (g / cm^3)
  ! cdens        ! Column density (particles/cm^2)
  ! temperature  ! Temperature (K)
  ! vpecul       ! Peculiar velocity (km/s)
  ! lambda0      ! rest wavelength (A)
  ! fvalue       ! oscillator strength
  ! mass         ! mass of atom (g)
  !
  ! output:
  ! tau          ! absorption spectrum
  ! velocw       ! Optical depth peculiar velocity (km/s)
  ! nionw        ! Optical depth ion number density
  ! rhow         ! optical depth weighted density (g/cm^3)
  ! tempw        ! Optical depth weighted temperature (K)
  !
  integer(kind=singleI) ::  nveloc10, i, iveloc, j, noff, ispec
  !
  real(kind=doubleR) :: minbother
  real(kind=doubleR) :: sigma_0, tfact, cfact, vmax, vmax10, vmax2, dvbin
  real(kind=doubleR) :: dvbin_inv, taumin, bpar, bpar_inv
  real(kind=doubleR) :: tauc, slim, vz, vdiff, vpar, dtau, off
  real(kind=doubleR) :: tauint, truetauint
  real(kind=doubleR) :: taufact, taufactor, vdiff0, vdiff1, vpar0, vpar1 
  real(kind=doubleR) :: erf, derf, erfvmax2
  !
  ! parameters of this transitions
  lambda0 = lambda_rest(ion,1)  ! rest-wavelength in A
  fvalue  = fosc(ion,1)         ! oscillator strength
  mass    = ion_mass(ion)       ! mass of ion (g)
  !
  ! Cross section in cm^2:
  if(ions(ion) .eq. '21cm') then
    !For 21cm emission we use Equation 3 from Furlanetto (2008)
    !To approximate the central optical depth of a gas cloud
    !as a function of N_HI, T_S and T_K.
    !
    ! We additionally make the assumption that T_S = T_K
    !
    cfact = 0.076 * (1000.0)**1.5 / 1e21 ! cm^3/s
    taufact = sqrt(pi) * cfact / 2. ! cm^3/s	
  else !not 21 cm emission
    sigma_0 = sqrt(3.*Pi*ThomsonCross/8.) * 1.e-8 * lambda0 * fvalue
    cfact = sigma_0 * LightSpeed / sqrt(pi) ! cm^3/s
    taufact = sigma_0 * LightSpeed / 2. ! cm^3/s
  endif
  !
  tfact = 2.0 * Boltz / mass ! erg/K/g	
  !
  ! note that v = [j-1,j>*dvbin is assigned to index j = 1,2,..,nveloc
  ! --> vmax = nveloc * dvbin 
  !
  if (vhubble(1) .ne. 0.) stop 'ERROR: vhubble(1) must be zero'
  !
  dvbin = vhubble(2) - vhubble(1) ! All in km/s
  vmax = nveloc * dvbin
  vmax10 = 10. * vmax
  vmax2 = vmax / 2.
  !
  if (integrate_thermprof_exactly) erfvmax2 = erf(vmax2)
  !
  dvbin_inv = 1. / dvbin
  taufact = taufact * dvbin_inv * 1e-5
  !
  nveloc10 = 10 * nveloc
  !
  if (lambda0 .gt. 1.001 * lyalpha) then 
    minbother = minbother_red
  else
    minbother = minbother_blue
  endif
  !
  taumin = minbother / dble(nveloc)
  !
  tau   = 0.d0 ! optical depth
  normw = 0.d0 ! weighing factor
  velocw = 0.d0
  nionw = 0.d0
  rhow  = 0.d0
  tempw = 0.d0
  !  
  do i = 1,nveloc
    if(cdens(i) .gt. 0) then
      !
      bpar = sqrt(tfact * temperature(i) + turbulence(i)) ! b-parameter in cm/s
      bpar_inv = 1./bpar
      !
      tauc = cfact * cdens(i) * bpar_inv ! Central optical depth
      !
      tauint = 0.
      taufactor = taufact * cdens(i)
      !
      if(ions(ion) .eq. '21cm') then
        tauc = tauc / temperature(i)**1.5 / bpar_inv
        taufactor = taufactor / temperature(i)**1.5 / bpar_inv
      endif
      !
      if (tauc .ge. taumin) then 
        !
        !     Integrate out to tau = taumin, 
        !     v/b (tau = taumin) = sqrt(log(tauc/taumin))
        slim = int(sqrt(log(tauc/taumin)))+1
        !
        !     convert to km/s
        bpar = bpar * 1.e-5
        if (bpar .le. 3.*vpixsizekms .and.  &
          .not. integrate_thermprof_exactly) then
          write(*,'("bpar =< 3*vpixsizekms",f12.5,f12.5)')  &
            bpar, vpixsizekms 
          write(*,'("Please decrease vpixsizekms in spec_wizard_modules.F90 (spectra),")')
          write(*,'("or integrate thermal profiles exactly.")')
          stop
        endif
        bpar_inv = 1./bpar
        !     note that the line below assumes that vpecul > -vmax10
        vz = mod(vhubble(i)+vpecul(i)+vmax10,vmax) ! in km/s
        if (vz .lt. 0.) stop 'ERROR: vpecul < -10*vmax'
        !     compute velocity bin to assign to
        !     note that v = [j-1,j>*dvbin is assigned to index j = 1,2,..,nveloc
        !     --> center of bin is at velocity (j-0.5)*dvbin
        iveloc = vz * dvbin_inv + 1 
        !
        !     resolved lines
        if (.not. limsigma) then
          !     complete convolution
          do j = 1,nveloc
            !     velocity offset from line centre in km/s :
            if (integrate_thermprof_exactly) then
              vdiff1 = dble(j) * dvbin - vz
              vdiff0 = vdiff1 - dvbin
              if (vdiff1 .gt. vmax2) vdiff1 = vdiff1 - vmax
              if (vdiff1 .lt. -vmax2) vdiff1 = vdiff1 + vmax
              if (vdiff0 .gt. vmax2) vdiff0 = vdiff0 - vmax
              if (vdiff0 .lt. -vmax2) vdiff0 = vdiff0 + vmax
              !     velocity off-set in units of Doppler parameter :
              vpar1 = vdiff1 * bpar_inv 
              vpar0 = vdiff0 * bpar_inv 
              derf = erf(vpar1) - erf(vpar0)
              !     correct case for which wrap around changes order
              if (vpar0 .gt. vpar1) derf = 2.*erfvmax2 + derf
              dtau = taufactor * derf
            else
              !     velocity offset from line centre in km/s :
              vdiff = abs(dble(j-0.5) * dvbin - vz) 
              if (vdiff .gt. vmax2) vdiff = vmax - vdiff
              !     velocity off-set in units of Doppler parameter :
              vpar = vdiff * bpar_inv 
              dtau = tauc * exp(-vpar**2)
              !     integrate tau(v), used to check resolution
              tauint = tauint + dtau*dvbin
            endif
            tau(j) = tau(j) + dtau
            !
            !  optical depth weighted values
            normw(j) = normw(j) + dtau
            velocw(j) = velocw(j) + dtau * vpecul(i) 
            nionw(j) = nionw(j) + dtau * nion(i)
            rhow(j)  = rhow(j)  + dtau * rho_tot(i)
            tempw(j) = tempw(j) + dtau * temperature(i) 
            !
          enddo
        else
          off = slim * bpar * dvbin_inv
          noff = min(int(off)+1,(nveloc-1)/2)
          do j = iveloc-noff,iveloc+noff
            if (integrate_thermprof_exactly) then
              vdiff1 = dble(j) * dvbin - vz
              vdiff0 = vdiff1 - dvbin
              if (vdiff1 .gt. vmax2) vdiff1 = vdiff1 - vmax
              if (vdiff1 .lt. -vmax2) vdiff1 = vdiff1 + vmax
              if (vdiff0 .gt. vmax2) vdiff0 = vdiff0 - vmax
              if (vdiff0 .lt. -vmax2) vdiff0 = vdiff0 + vmax
              !     velocity off-set in units of Doppler parameter :
              vpar1 = vdiff1 * bpar_inv 
              vpar0 = vdiff0 * bpar_inv 
              derf = erf(vpar1) - erf(vpar0)
              !     correct case for which wrap around changes order
              if (vpar0 .gt. vpar1) derf = 2.*erfvmax2 + derf
              dtau = taufactor * derf
            else
              vdiff = abs(dble(j-0.5) * dvbin - vz) 
              if (vdiff .gt. vmax2) vdiff = vmax - vdiff
              !     velocity off-set in units of Doppler parameter :
              vpar = vdiff * bpar_inv 
              dtau = tauc * exp(-vpar**2)
              tauint = tauint + dtau*dvbin
            endif
            ispec = mod(j-1+nveloc10,nveloc) + 1
            tau(ispec) = tau(ispec) + dtau
            !
            !     optical depth weighted values
            normw(ispec) = normw(ispec) + dtau
            velocw(ispec)= velocw(ispec) + dtau * vpecul(ispec) 
            nionw(ispec) = nionw(ispec) + dtau *  nion(i)
            rhow(ispec)  = rhow(ispec)  + dtau * rho_tot(i)
            tempw(ispec) = tempw(ispec) + dtau * temperature(i) 
          enddo
        endif               ! limsigma
        !     theoretical integral of tau(v) over full thermal profile
        if (.not. integrate_thermprof_exactly) then
          truetauint = tauc * sqrt(Pi) * bpar
          if ( abs(tauint-truetauint) .gt. minbother * dvbin) then
            write(*,*)
            write(*,'("lambda0 = ",f12.3," A")') lambda0
            write(*,'("fvalue  = ",f12.3)') fvalue
            write(*,'("tauc    = ",g12.3)') tauc
            write(*,'("bpar    = ",g12.3," km/s")') bpar
            write(*,*) i,tauint,truetauint,abs(tauint-truetauint)
            write(*,*) minbother,dvbin,minbother*dvbin
            write(*,'("Line is unresolved!")')
            write(*,'(a,a)') "Try decreasing vpixkms or ", &
              "increasing minbother in pars.inc."
            stop
          endif
        endif
      endif                  ! tauc
    endif ! cdens > 0
  enddo
  !
  !     normalize density & temperature
  do i = 1,nveloc
    if (normw(i) .gt. 0.) then
      velocw(i) = velocw(i) / normw(i)
      nionw(i) = nionw(i)   / normw(i)
      rhow(i)  = rhow(i)    / normw(i)
      tempw(i) = tempw(i)   / normw(i)
    endif
  enddo
  !
end subroutine computespectrum


subroutine shift ()
  ! cyclically shift spectrum such that minimum optical depth is at start
  use spectra
  use runtime
  implicit none
  !
  integer(kind=singleI) :: iloc(1), ion, i, icyc
  !
  write(*,*) 'Warning: better not to cycle!'
  !
  ! find location of minimum optical depth
  tau(:)  = tau_ion(1,:)
  iloc    = minloc(tau)
  icshift = iloc(1)
  !
  ! shift optical depth for each ion
  do ion=1, nion
    do i=1, nveloc
      icyc = i - icshift + 1
      if (icyc .le. 0) icyc = icyc + nveloc
      tau(icyc) = tau_ion(ion,i)
    enddo
    tau_ion(ion,:) = tau(:)
    !
    do i=1, nveloc
      icyc = i - icshift
      if (icyc .le. 0) icyc = icyc + nveloc
      tau(icyc) = Rho_z_ion(ion,i)
    enddo
    Rho_z_ion(ion,:) = tau(:)
    do i=1, nveloc
      icyc = i - icshift
      if (icyc .le. 0) icyc = icyc + nveloc
      tau(icyc) = Temp_z_ion(ion,i)
    enddo
    Temp_z_ion(ion,:) = tau(:)
  enddo
  !
  return
end subroutine shift


! +++++++++++++++++++++++++++++++++++++++ ...spectra +++++++++++++++++++++++++++++++++++++++++++++ !

! +++++++++++++++++++++++++++++++++++++++ short spectra... +++++++++++++++++++++++++++++++++++++++ !


subroutine allocate_spectra_short()
  use numbers
  use spectra
  use atomic_data
  implicit none
  !
  integer, save :: nveloc_old = -1
  !
  if(nveloc .eq. nveloc_old) return
  !
  nveloc_old = nveloc
  if(allocated(rho_tot)) deallocate(rho_tot)
  if(allocated(temp_tot)) deallocate(temp_tot)
  if(allocated(met_tot)) deallocate(met_tot)
  if(allocated(veloc_tot)) deallocate(veloc_tot)
  !
  if(allocated(n_ion)) deallocate(n_ion)
  if(allocated(temp_ion)) deallocate(temp_ion)
  if(allocated(veloc_ion)) deallocate(veloc_ion)
  if(allocated(rho_ion)) deallocate(rho_ion)
  !
  if(allocated(tau)) deallocate(tau)
  if(allocated(normw)) deallocate(normw)
  if(allocated(work)) deallocate(work)
  if(allocated(work2)) deallocate(work2)
  if(allocated(voc)) deallocate(voc)
  if(allocated(vocsim)) deallocate(vocsim)
  if(allocated(cdens)) deallocate(cdens)
  if(allocated(vhubble)) deallocate(vhubble)
  if(allocated(tau_ion)) deallocate(tau_ion)
  if(allocated(flux_ion)) deallocate(flux_ion)
  if(allocated(rhow)) deallocate(rhow)
  if(allocated(tempw)) deallocate(tempw)
  if(allocated(nionw)) deallocate(nionw)
  if(allocated(velocw)) deallocate(velocw)
  if(allocated(Rho_z_ion)) deallocate(rho_z_ion)
  if(allocated(Temp_Z_ion)) deallocate(Temp_Z_ion)
  if(allocated(veloc_Z_ion)) deallocate(veloc_Z_ion)
  if(allocated(nion_Z_ion)) deallocate(nion_Z_ion)
  if(allocated(spectrum_boundary)) deallocate(spectrum_boundary)
  if(allocated(cdens_ion_integrated)) deallocate(cdens_ion_integrated)
  !
  if(allocated(turbulence)) deallocate(turbulence)
  !
  allocate(rho_tot(nveloc),temp_tot(nveloc),met_tot(nveloc),veloc_tot(nveloc))
  allocate(n_ion(nion,nveloc),temp_ion(nion,nveloc),veloc_ion(nion,nveloc),rho_ion(nion,nveloc))
  allocate(tau(nveloc),normw(nveloc),work(nveloc),work2(nveloc),vhubble(nveloc),cdens(nveloc),&
    rhow(nveloc),tempw(nveloc),voc(nveloc),vocsim(nveloc),turbulence(nveloc),velocw(nveloc),nionw(nveloc))
  allocate(tau_ion(nion,nveloc),flux_ion(nion,nveloc))
  allocate(cdens_ion_integrated(nion))
  allocate(Rho_Z_ion(nion,nveloc),Temp_Z_Ion(nion,nveloc),veloc_Z_ion(nion,nveloc),nion_Z_ion(nion,nveloc))
  allocate(spectrum_boundary(nveloc))
  !
end subroutine allocate_spectra_short



subroutine update_short_spectrum(particlefile, nlos)
  use numbers
  use spectra
  use projection_parameters
  use runtime
  use my_mpi
  use header
  use hdf5_wrapper
  use particledata, only : rhocb
  implicit none
  !
  character(*), intent(in)          :: particlefile
  character(len=400)                :: outputfile
  integer(kind=singleI), intent(in) :: nlos
  ! local variables
  logical                :: file_exists
  integer(kind=singleI)  :: file_handle, ion, number_of_transitions, line, count
  character(len=120)     :: GroupName, VarName, MassWeightedGroup, ElementGroup
  character(len=10), allocatable    :: all_ions(:)
  real(kind=doubleR), allocatable   :: all_ion_mass(:), all_lambda_rest(:), all_fosc(:)
  integer :: mpi_err

  ! calculate number of particles that contributed to each sight line

#ifdef MPI
  call mpi_reduce(ncontribute, ncontribute_global, nlos, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, MPI_err)
#else
  ncontribute_global = ncontribute 
#endif
  
  ! check for existence of output file
  if(MyPE == 0) then
     inquire(file=SpectrumFile,exist=file_exists)
     if(.not. file_exists)then
        call abortrun('The required output file '//trim(SpectrumFile)//' does not exist')
     endif
     !
     call hdf5_open_file(file_handle,trim(SpectrumFile))
     call hdf5_write_data(file_handle,'Projection/ncontr', ncontribute_global)
     call hdf5_close_file(file_handle)
  endif
end subroutine update_short_spectrum

subroutine write_short_spectrum(particlefile, los_number, nlos)
  use numbers
  use spectra
  use noisedata
  use projection_parameters
  use runtime
  use my_mpi
  use header
  use hdf5_wrapper
  use particledata, only : rhocb
  implicit none
  !
  character(*), intent(in)          :: particlefile
  character(len=400)                :: outputfile
  integer(kind=singleI), intent(in) :: los_number, nlos
  ! local variables
  logical                :: file_exists
  integer(kind=singleI)  :: file_handle, ion, number_of_transitions, line, count
  character(len=120)     :: GroupName, VarName, MassWeightedGroup, ElementGroup
  character(len=10), allocatable    :: all_ions(:)
  real(kind=doubleR), allocatable   :: all_ion_mass(:), all_lambda_rest(:), all_fosc(:)


  ! check for existence of output file
  inquire(file=SpectrumFile,exist=file_exists)
  if(.not. file_exists)then
     call abortrun('The required output file '//trim(SpectrumFile)//' does not exist')
  endif
  !
  call hdf5_open_file(file_handle,trim(SpectrumFile))
  !
  if(ispec == 0) then
    ! total number of transitions included
    number_of_transitions = sum(nlines)
    allocate(all_ions(number_of_transitions),all_lambda_rest(number_of_transitions),all_fosc(number_of_transitions),&
      all_ion_mass(number_of_transitions))
    count = 0
    !
    do ion=1, nion
      do line=1, nlines(ion)
        count = count + 1
        all_ions(count) = trim(ions(ion))
        all_lambda_rest(count) = lambda_rest(ion,line)
        all_fosc(count)        = fosc(ion,line)
        all_ion_mass(count)    = ion_mass(ion)
      enddo
    enddo

    ! write hubble velocity
    VarName = 'VHubble_KMpS'
    ! last argument: compression 1-9 normal compression, more and more, +10: more clever compression mechanisms
    call hdf5_write_data(file_handle, VarName, vhubble, gzip=16)

    ! write projection parameters
    call hdf5_write_data(file_handle,'Projection/x_fraction_array',x_fraction_array, gzip=16)
    call hdf5_write_data(file_handle,'Projection/y_fraction_array',y_fraction_array, gzip=16)
    call hdf5_write_attribute(file_handle,'Projection/nspec',nspec)
    !
  endif
  !
!!$  write (*,*) ' MyPE= ',MyPE,' rhocb= ',rhocb

  ! Create group for current spectrum
  call mylabel('/Spectrum',ispec,GroupName)
  call hdf5_create_group(file_handle, GroupName)
  !
  ! properties of sight line
  VarName = trim(GroupName)//'/X-position'
  call hdf5_write_attribute(file_handle,VarName,x_comoving)
  VarName = trim(GroupName)//'/Y-position'
  call hdf5_write_attribute(file_handle,VarName,y_comoving)
  VarName = trim(GroupName)//'/Z-position'
  call hdf5_write_attribute(file_handle,VarName,z_comoving)
  VarName = trim(GroupName)//'/theta'
  call hdf5_write_attribute(file_handle,VarName,theta_projection)
  VarName = trim(GroupName)//'/phi'
  call hdf5_write_attribute(file_handle,VarName,phi_projection)
  VarName = trim(GroupName)//'/x-axis'
  call hdf5_write_attribute(file_handle,VarName,x_axis)
  VarName = trim(GroupName)//'/y-axis'
  call hdf5_write_attribute(file_handle,VarName,y_axis)
  VarName = trim(GroupName)//'/z-axis'
  call hdf5_write_attribute(file_handle,VarName,z_axis)
  !
  ! properties of each ion
  if (output_zspaceopticaldepthweighted_values) then
    do ion=1,nion
      ElementGroup = trim(GroupName)//'/'//trim(ions(ion))//'/RedshiftSpaceOpticalDepthWeighted'
      call hdf5_create_group(file_handle, ElementGroup)
      VarName = trim(ElementGroup)//'/'//'LOSPeculiarVelocity_KMpS'
      call hdf5_write_data(file_handle, trim(VarName),veloc_z_ion(ion,:), gzip=16)
      VarName = trim(ElementGroup)//'/'//'NIon_CM3'
      call hdf5_write_data(file_handle, trim(VarName),nion_z_ion(ion,:), gzip=16)
      VarName = trim(ElementGroup)//'/'//'OverDensity'
      call hdf5_write_data(file_handle, trim(VarName),rho_z_ion(ion,:), gzip=16)
      VarName = trim(ElementGroup)//'/'//'Temperature_K'
      call hdf5_write_data(file_handle, trim(VarName),temp_z_ion(ion,:), gzip=16)
    enddo
  endif
  !
  !Real space, Nion weighted values
  if (output_realspacenionweighted_values) then
    do ion=1, nion
      ElementGroup = trim(GroupName)//'/'//trim(ions(ion))//'/RealSpaceNionWeighted/'
      call hdf5_create_group(file_handle, ElementGroup)
      ! Real space, n_ion-weighted LOS velocity
      VarName = trim(ElementGroup)//'LOSPeculiarVelocity_KMpS'
      call hdf5_write_data(file_handle, trim(VarName),veloc_ion(ion,:), gzip=16)
      ! Real space, pixel n_ion
      VarName = trim(ElementGroup)//'NIon_CM3'
      call hdf5_write_data(file_handle, trim(VarName),n_ion(ion,:), gzip=16)
      ! Real space, n_ion.-weighted overdensity
      VarName = trim(ElementGroup)//'OverDensity'
      call hdf5_write_data(file_handle,trim(VarName),rho_ion(ion,:), gzip=16) !Old version: n_ion(ion,:)*ion_mass(ion)/rhocb
      ! Real space, n_ion-weighted temperature. 
      VarName = trim(ElementGroup)//'Temperature_K'
      call hdf5_write_data(file_handle, trim(VarName),temp_ion(ion,:), gzip=16)
    enddo
  endif
  !
  ! Mass-weighted properties
  if(output_realspacemassweighted_values)then
    MassWeightedGroup = trim(GroupName)//'/RealSpaceMassWeighted'
    call hdf5_create_group(file_handle, MassWeightedGroup)
    VarName      = trim(MassWeightedGroup)//'/LOSPeculiarVelocity_KMpS'
    call hdf5_write_data(file_handle,trim(varname),veloc_tot, gzip=16)
    VarName      = trim(MassWeightedGroup)//'/OverDensity'
    call hdf5_write_data(file_handle,trim(varname),rho_tot, gzip=16)
    VarName      = trim(MassWeightedGroup)//'/Temperature_K'
    call hdf5_write_data(file_handle,trim(varname),temp_tot, gzip=16)
    !
    ! Metallicity is a mass fraction
    VarName      = trim(MassWeightedGroup)//'/MetalMassFraction'
    call hdf5_write_data(file_handle,trim(varname),met_tot, gzip=16)
  endif
  !
  do ion=1, nion
    VarName = trim(GroupName)//'/'//trim(ions(ion))//'/Flux'
    call hdf5_write_data(file_handle, trim(VarName), flux_ion(ion,1:nveloc), gzip=16)
    !
    !if(generate_noise)then
    !  VarName = trim(GroupName)//'/'//trim(ions(ion))//'/Gaussian_deviate'
    !  call hdf5_write_data(file_handle, trim(VarName), binned_noise_random)
    !  VarName = trim(GroupName)//'/'//trim(ions(ion))//'/Noise_Sigma'
    !  call hdf5_write_data(file_handle, trim(VarName), binned_noise_sigma)
    !  write(*,*) 'done writing noise info in the output'
    !endif
    !
    VarName = trim(GroupName)//'/'//trim(ions(ion))//'/OpticalDepth'
    call hdf5_write_data(file_handle, trim(VarName), tau_ion(ion,1:nveloc), gzip=16)
    VarName = trim(GroupName)//'/'//trim(ions(ion))//'/LogTotalIonColumnDensity'
    call hdf5_write_data(file_handle, trim(VarName), log10(cdens_ion_integrated(ion)))
  enddo
  !
  if(gimic)then
    VarName = trim(GroupName)//'/Boundary'
    call hdf5_write_data(file_handle, trim(VarName), spectrum_boundary(1:nveloc))
  endif
  !
  call hdf5_close_file(file_handle)
  !
end subroutine write_short_spectrum


subroutine convolve_short_spectrum()
  use numbers
  use physical_constants
  use spectra
  use runtime
  use hdf5_wrapper
  use my_mpi
  implicit none
  ! local variables
  real(kind=doubleR)    :: sigmakms, b, norm
  integer(kind=singleI) :: i, j, off, ion
  real(kind=doubleR), allocatable :: gauss(:)
  real(kind=doubleR), allocatable ::convl_flux(:),convl_flux_convolved(:)
  !
  if (verbose .or. MyPE == 0) then
    write(*,*)
    write(*,'( &
    "Convolving with instrumental profile, fwhm = ", &
    f8.3," km/s...")') fwhm
    write(*,*)
  endif
  !
  vpixsize = vhubble(2)-vhubble(1) !in km/s
  !
  ! Compute sigma in km/s
  sigmakms = fwhm / 2.35482
  ! Compute sigma in units of pixels
  b = sigmakms / vpixsize
  !
  ! For convolution with instrumental profile we need to Fourier 
  ! transform, we thus need to increase the array so that it is a power of 2.
  nvpix   = int(2**(aint(log(dble(nveloc))/log(2.)) + 1))
  !
  ! Create normalized Gaussian in wrap-around order
  allocate(gauss(nvpix))
  norm = 1d0 / (2d0 * b * b)
  do i = 0, nvpix-1
    if (i .le. nvpix-1) then 
      if (i .le. nvpix/2) then
        j = i
      else
        j = i - nvpix
      endif
      if (abs(j) .lt. 10.*b) then
        gauss(i+1) = exp(-(dble(j)**2)*norm)
      else
        gauss(i+1) = 0d0
      endif
    else
      gauss(i+1) = 0d0
    endif
  enddo
  ! Make sure Gaussian is normalized.
  gauss  = gauss / sum(gauss)
  !
  allocate(convl_flux(nvpix),convl_flux_convolved(2*nvpix))
  !
  do ion=1, nion
    convl_flux(:) = 0.0
    convl_flux(1:nveloc) = flux_ion(ion,:)
    !
    !Now copy periodic copies of the flux signal into the zero buffer
    !to avoid aliasing effects
    do i=nveloc+1,nvpix
      off = i-nveloc
      if (off .lt. (nvpix-nveloc)/2.) then
        convl_flux(i) = convl_flux(i-nveloc)
      else
        convl_flux(i) = convl_flux(i-(nvpix-nveloc))
      endif
    enddo
    convl_flux_convolved(:) = 0.0
    !
    call convlv(convl_flux,nvpix,gauss,nvpix,1,convl_flux_convolved)
    flux_ion(ion,:) = convl_flux_convolved(1:nveloc)
    !
  enddo
  !
  deallocate(convl_flux,convl_flux_convolved)
  !
end subroutine convolve_short_spectrum


! ++++++++++++++++++++++++++++++++++++ ...short spectra ++++++++++++++++++++++++++++++++++++++++++ !

! ++++++++++++++++++++++++++++++++++++++ do_long_spectra... ++++++++++++++++++++++++++++++++++++++ !


subroutine allocate_spectra_long()
  use numbers
  use spectra
  implicit none
  !
  if(allocated(lambda)) deallocate(lambda)
  !
  allocate(cdens_ion_integrated(nion))
  allocate(lambda(nvpix),voverc(nvpix))
  allocate(voverc_realspace(nppix),redshift_realspace(nppix))
  allocate(tau_long(nion,nvpix))
  allocate(tau_long_strongest(nion,nvpix))
  allocate(temp_ion_long(nion,nppix),n_ion_long(nion,nppix),rho_ion_long(nion,nppix), &
           veloc_ion_long(nion, nppix))
  allocate(temp_z_ion_long(nion,nvpix), rho_z_ion_long(nion,nvpix), &
           veloc_z_ion_long(nion,nvpix))
  allocate(rho_long(nppix),temp_long(nppix),met_long(nppix),veloc_long(nppix))
  allocate(flux(nvpix), flux_convolved(2*nvpix))
  !
  tau_long(:,:)             = 0.0
  tau_long_strongest(:,:)   = 0.0
  temp_ion_long(:,:)        = 0.0
  rho_ion_long(:,:)         = 0.0
  n_ion_long(:,:)           = 0.0
  veloc_ion_long(:,:)       = 0.0
  temp_z_ion_long(:,:)      = 0.0
  rho_z_ion_long(:,:)       = 0.0
  veloc_z_ion_long(:,:)     = 0.0
  cdens_ion_integrated(:)   = 0.0
  rho_long(:)               = 0.0
  temp_long(:)              = 0.0
  met_long(:)               = 0.0
  veloc_long(:)             = 0.0        
  flux(:)                   = 0.0
  flux_convolved(:)         = 0.0
  !
end subroutine allocate_spectra_long


subroutine insertspectra(zcurrent_next)
  use numbers
  use spectra
  use header, only : BoxSize, HubbleParam, Omega0, OmegaLambda
  use runtime
  use physical_constants
  use my_spline_interpolation
  implicit none
  ! local variables
  integer(kind=singleI) :: j, ion, i
  real(kind=doubleR), intent(out) :: zcurrent_next

  ! local variables
  real(kind=doubleR) :: boxkms, logl, minvoc, maxvoc, CurrentHubbleCt, pixvoc, term, &
       taufactor, minbother, y
  character(len=120) :: outfile
  logical :: loginterpolate
  logical, parameter :: is_positive = .true.
  !  
  CurrentHubbleCt = 100. * HubbleParam *  &
       sqrt(1. + omega0*(1./acurrent-1.) + OmegaLambda* &
       (acurrent**2-1.)) /acurrent
  boxkms = BoxSize / HubbleParam * acurrent * CurrentHubbleCt     

  !
  ! lambda = lambda0 (1+zcurent) log_e(vhubble/c)
  
  ! Compute v/c of simulation pixels when mapped onto long
  ! spectrum. Actually, this is NOT v/c, but v/c -
  ! ln(lambda_0). The rest wavelength is added later!
  pixvoc = (abs(boxkms)*1.d5) / LightSpeed  / dble(nveloc)  ! dv/c per pixel of short spectrum
  do i = 1, nveloc   ! z + 1
     voc(i) = dble(i-1) * pixvoc
  enddo
  voc = voc + log(1.+zcurrent)
  
  ! Next redshift; exp: limit case, factors of (1 + Delta v / c) for each pixel
  ! lim n-> inf (1 + x/n)^n = exp(x)
  zcurrent_next = (1.d0+zcurrent) *  exp(voc(nveloc) - voc(1) + pixvoc) - 1.

  !
  if(docycling) call shift () ! cyclically shift spectrum to have minimum HI optical depth at start/end

  ! spline interpolate real-space density, density weighted temperature 
  ! voc is log(1 + z)
  vocsim = voc
  minvoc = vocsim(1)
  maxvoc = min(vocsim(nveloc),log(1.+zqso))

  ! debug info real-space values
  !write(*,*)'Insertspectra inputs real space'
  !write(*,'("minvoc: ",f7.4," maxvoc: ", f7.4)') minvoc, maxvoc
  !write(*,'("vocsim: ",f7.4,", ", f7.4, ", ", f7.4, " ... ", f7.4)') & 
  !        vocsim(1), vocsim(2), vocsim(3), vocsim(nveloc)
  !write(*,'("voverc_realspace: ",f7.4,", ", f7.4, ", ",f7.4," ... ",f7.4)') & 
  !        voverc_realspace(1), voverc_realspace(2), voverc_realspace(3), voverc_realspace(nppix)
                  
  if(output_realspacemassweighted_values)then
     call spline_interpolate(&
          nveloc,vocsim,rho_tot &
          ,minvoc,maxvoc,1,small_rho &
          ,nppix,voverc_realspace,rho_long,is_positive=.true.)
     
     call spline_interpolate(&
          nveloc,vocsim,temp_tot &
          ,minvoc,maxvoc,1,small_temp &
          ,nppix,voverc_realspace,temp_long,is_positive=.true.)
     call spline_interpolate(&
          nveloc,vocsim,met_tot &
          ,minvoc,maxvoc,1,small_metallicity &
          ,nppix,voverc_realspace,met_long,is_positive=.true.)
     call spline_interpolate(&
          nveloc,vocsim,veloc_tot &
          ,minvoc,maxvoc,1,small_velocity &
          ,nppix,voverc_realspace,veloc_long,is_positive=.false.)
  endif
  !
  if (output_realspacenionweighted_values) then
     do ion = 1, nion
        call spline_interpolate(&
             nveloc,vocsim,n_ion(ion,:) &
             ,minvoc,maxvoc,1,small_rho &
             ,nppix,voverc_realspace,n_ion_long(ion,:),is_positive=.true.)
        
        call spline_interpolate(&
             nveloc,vocsim,temp_ion(ion,:) &
             ,minvoc,maxvoc,1,small_temp &
             ,nppix,voverc_realspace,temp_ion_long(ion,:),is_positive=.true.)
        call spline_interpolate(&
             nveloc,vocsim,rho_ion(ion,:) &
             ,minvoc,maxvoc,1,small_metallicity &
             ,nppix,voverc_realspace,rho_ion_long(ion,:),is_positive=.true.)
        call spline_interpolate(&
             nveloc,vocsim,veloc_ion(ion,:) &
             ,minvoc,maxvoc,1,small_velocity &
             ,nppix,voverc_realspace,veloc_ion_long(ion,:),is_positive=.false.)
    enddo
  endif
  !
  ion_loop: do ion = 1, nion
     taumax(ion) = maxval(tau_ion(ion,:))
     !
     lines: do j = 1, nlines(ion)
        logl   = log(lambda_rest(ion,j)/minlambda)
        minvoc = min(voc(1)      ,log(1.+zqso))+ logl
        maxvoc = min(voc(nveloc) ,log(1.+zqso))+ logl        
        !
        ! redshift-space ion-weighted density, temperature and velocity
        if (output_zspaceopticaldepthweighted_values .and. (j .eq. 1)) then
           vocsim(:)   = voc(:) + logl
           call spline_interpolate(&
                nveloc,vocsim,rho_z_ion(ion,:) &
                ,minvoc,maxvoc,ion,small_rho &
                ,nvpix,voverc,rho_z_ion_long(ion,:),is_positive=.true.)
           call spline_interpolate(&
                nveloc,vocsim,temp_z_ion(ion,:)&
                ,minvoc,maxvoc,ion,small_temp&
                ,nvpix,voverc,temp_z_ion_long(ion,:),is_positive=.true.)
           call spline_interpolate(&
                nveloc,vocsim,veloc_z_ion(ion,:)&
                ,minvoc,maxvoc,ion,small_velocity&
                ,nvpix,voverc,veloc_z_ion_long(ion,:),is_positive=.false.)
        endif 
        if (minvoc .le. voverc(nvpix) .and.  maxvoc .ge. voverc(1)) then 
           if (lambda_rest(ion,j) .gt. 1.001 * lyalpha) then
              minbother = minbother_red
           else
              minbother = minbother_blue
           endif
           ! scale factor to obtain optical depth for this transition
           taufactor=(lambda_rest(ion,j)*fosc(ion,j))/ (lambda_rest(ion,1)*fosc(ion,1))

           tau_limit: if (taumax(ion) * taufactor .gt. minbother) then
              vocsim(:)   = voc(:) + logl ! pixel velocity of this ion

              ! optical depth              
              call spline_interpolate(&
                   nveloc,vocsim,tau_ion(ion,:)*taufactor &
                   ,minvoc,maxvoc,ion,minbother &
                   ,nvpix,voverc,tau_long(ion,:) &
                   ,loginterpolate=.true.,is_positive=.true.)
              
              !If this is the strongest transition put it into tau_long_strongest
              if (j .eq. 1) then
                  call spline_interpolate(&
                  nveloc,vocsim,tau_ion(ion,:)*taufactor &
                  ,minvoc,maxvoc,ion,minbother &
                  ,nvpix,voverc,tau_long_strongest(ion,:) &
                  ,loginterpolate=.true.,is_positive=.true.)
              endif

           endif tau_limit
        endif
     enddo lines
  enddo ion_loop

end subroutine insertspectra


subroutine write_long_spectrum()
  use numbers
  use hdf5_wrapper
  use runtime
  use spectra
  use noisedata
  use physical_constants, only: pi, G,h0, LightSpeed
  use header
  use random_numbers, only: seed
  use particledata, only : rhocb
  implicit none
  ! local variables
  integer               :: file_handle, ion, infile_handle
  character(len=120)    :: outfile, GroupName, ThisSpectrum, VarName, ElementGroup, &
                           inputfile, IonWeightedGroup, MassWeightedGroup, &
                           RedshiftIonWeightedGroup
  !
  call hdf5_open_file(file_handle,trim(SpectrumFile))
  !
  if(ispec == 1) then
     !
     if (output_frequency) then
        ! write wavelengths
        call hdf5_write_data(file_handle, '/Frequency_MHz', LightSpeed / binned_lambda / 1e-8 / 1e6, gzip=16)
     else
        ! write wavelengths
        call hdf5_write_data(file_handle, '/Wavelength_Ang', binned_lambda, gzip=16)
     endif
     !
  endif
  !
  ! write current spectrum
  call mylabel('/Spectrum',ispec,ThisSpectrum)
  call hdf5_create_group(file_handle, ThisSpectrum)
  !
  VarName = trim(ThisSpectrum)//'/Flux'
  call hdf5_write_data(file_handle, trim(VarName), binned_flux, gzip=16)
  !
  if(generate_noise)then
     VarName = trim(ThisSpectrum)//'/Gaussian_deviate'
     call hdf5_write_data(file_handle, trim(VarName), binned_noise_random, gzip=16)
     VarName = trim(ThisSpectrum)//'/Noise_Sigma'
     call hdf5_write_data(file_handle, trim(VarName), binned_noise_sigma, gzip=16)
     write(*,*) 'done writing noise info in the output'
  endif
  !
  ! random seeds
  VarName = trim(ThisSpectrum)//'/RandomSeeds'
  call hdf5_write_attribute(file_handle,VarName,Seed)
  !
  ! properties of each ion
  do ion=1,nion
     ElementGroup = trim(ThisSpectrum)//'/'//trim(ions(ion))
     call hdf5_create_group(file_handle, ElementGroup)
     VarName = trim(ElementGroup)//'/'//'RedshiftSpaceOpticalDepthOfStrongestTransition'
     call hdf5_write_data(file_handle, trim(VarName),binned_tau_ion_strongest(ion,:), gzip=16)
     VarName = trim(ElementGroup)//'/'//'LogTotalIonColumnDensity'
     call hdf5_write_data(file_handle, trim(VarName),cdens_ion_integrated(ion))
  enddo
  !
  if (output_zspaceopticaldepthweighted_values) then
     do ion=1,nion
          ElementGroup = trim(ThisSpectrum)//'/'//trim(ions(ion))
          RedshiftIonWeightedGroup = trim(ElementGroup)//'/RedshiftSpaceOpticalDepthWeighted'
          call hdf5_create_group(file_handle, RedshiftIonWeightedGroup)
          VarName = trim(RedshiftIonWeightedGroup)//'/'//'OverDensity'
          call hdf5_write_data(file_handle, trim(VarName),binned_rho_z_ion(ion,:), gzip=16)
          VarName = trim(RedshiftIonWeightedGroup)//'/'//'Temperature_K'
          call hdf5_write_data(file_handle, trim(VarName),binned_temp_z_ion(ion,:), gzip=16)
          VarName = trim(RedshiftIonWeightedGroup)//'/LOSPeculiarVelocity_KMpS'
          call hdf5_write_data(file_handle, trim(VarName),binned_veloc_z_ion(ion,:), gzip=16)
      enddo
  endif
  !
  if (output_realspacenionweighted_values) then
      do ion=1,nion
          ElementGroup = trim(ThisSpectrum)//'/'//trim(ions(ion))
          IonWeightedGroup = trim(ElementGroup)//'/RealSpaceNionWeighted'
          call hdf5_create_group(file_handle, IonWeightedGroup)
          VarName = trim(ElementGroup)//'/'//'NIon_CM3'
          call hdf5_write_data(file_handle, trim(VarName),binned_n_ion(ion,:), gzip=16)
          VarName = trim(IonWeightedGroup)//'/'//'OverDensity'
          call hdf5_write_data(file_handle, trim(VarName),binned_rho_ion(ion,:), gzip=16)
          VarName = trim(IonWeightedGroup)//'/'//'Temperature_K'
          call hdf5_write_data(file_handle, trim(VarName),binned_temp_ion(ion,:), gzip=16)
          VarName = trim(IonWeightedGroup)//'/LOSPeculiarVelocity_KMpS'
          call hdf5_write_data(file_handle, trim(VarName),binned_veloc_ion(ion,:), gzip=16)
      enddo
  endif
  ! Mass-weighted properties along the los
  if (output_realspacemassweighted_values) then
      MassWeightedGroup = trim(ThisSpectrum)//'/RealSpaceMassWeighted'
      call hdf5_create_group(file_handle, MassWeightedGroup)
      VarName      = trim(MassWeightedGroup)//'/LOSPeculiarVelocity_KMpS'
      call hdf5_write_data(file_handle,trim(varname),binned_veloc, gzip=16)
      VarName      = trim(MassWeightedGroup)//'/OverDensity'
      call hdf5_write_data(file_handle,trim(varname),binned_rho, gzip=16)
      VarName      = trim(MassWeightedGroup)//'/Temperature_K'
      call hdf5_write_data(file_handle,trim(varname),binned_temp, gzip=16)
      !
      ! Metallicity is a mass fraction
      VarName      = trim(MassWeightedGroup)//'/MetalMassFraction'
      call hdf5_write_data(file_handle,trim(varname),binned_met, gzip=16)
  endif
  !
  if (output_realspacemassweighted_values .or. output_realspacenionweighted_values) then
      VarName      = trim(ThisSpectrum)//'/Redshift_RealSpace'
      call hdf5_write_data(file_handle,trim(varname),binned_redshift_realspace, gzip=16)
  endif
  !
  ! output info of small spectra used
  GroupName = trim(ThisSpectrum)//'/ShortSpectraInfo'
  call hdf5_create_group(file_handle, GroupName)
  VarName = trim(GroupName)//'/RandomSeeds'  
  call hdf5_write_data(file_handle, VarName, seed)
  VarName = trim(GroupName)//'/NumberOfShortSpectra'
  call hdf5_write_data(file_handle, trim(VarName), nsimfile_used)
  !
  !if (gimic) then ! record this for all input files
  if (.not. use_snapshot_file) then
      VarName = trim(GroupName)//'/FileUsed'
      call hdf5_write_data(file_handle,trim(VarName), losfile_used(1:nsimfile_used), gzip=16)
      VarName = trim(GroupName)//'/LosUsed'
      call hdf5_write_data(file_handle,trim(VarName), los_used(1:nsimfile_used), gzip=16)
  endif
  VarName = trim(GroupName)//'/Icshift'
  call hdf5_write_data(file_handle,trim(VarName), icshift_used(1:nsimfile_used), gzip=16)
  VarName = trim(GroupName)//'/x_simunits'
  call hdf5_write_data(file_handle,trim(VarName), x_physical_used(1:nsimfile_used), gzip=16)
  VarName = trim(GroupName)//'/y_simunits'
  call hdf5_write_data(file_handle,trim(VarName), y_physical_used(1:nsimfile_used), gzip=16)
  VarName = trim(GroupName)//'/Ibfactor'
  call hdf5_write_data(file_handle,trim(VarName), ibfactor_used(1:nsimfile_used), gzip=16)
  VarName = trim(GroupName)//'/x-axis'
  call hdf5_write_data(file_handle,trim(VarName), x_axis_used(1:nsimfile_used), gzip=16)
  VarName = trim(GroupName)//'/y-axis'
  call hdf5_write_data(file_handle,trim(VarName), y_axis_used(1:nsimfile_used), gzip=16)
  VarName = trim(GroupName)//'/z-axis'
  call hdf5_write_data(file_handle,trim(VarName), z_axis_used(1:nsimfile_used), gzip=16)
  !endif
  !
  call hdf5_close_file(file_handle)
  !
end subroutine write_long_spectrum


subroutine rebin_spectrum
  use numbers
  use runtime
  use spectra
  use my_rebin
  use hdf5_wrapper
  use noisedata
  use my_mpi
  implicit none
  ! local variables
  integer :: i,j,k,ninbin,file_handle,ion
  character(len=120) :: filename
  integer(kind=singleI), save :: n_binned_old=-1, n_binned_realspace_old=-1
  !
  if (verbose .and. MyPE == 0)  &
       write(*,'("Rebinning onto ",f7.4," A pixels...")') pixsize
  !
  ! output wavelength
  n_binned_flux = int((maxlambda-minlambda) / pixsize) + 1
  !
  if(n_binned_old /= n_binned_flux)then
     if(n_binned_old .gt. 0) then
        write (*,*) ' error: change in spectrum size? '
        stop
     endif
     n_binned_old = n_binned_flux
     !
     allocate(binned_lambda(n_binned_flux), binned_flux(n_binned_flux))
     allocate(binned_temp_z_ion(nion,n_binned_flux),binned_rho_z_ion(nion,n_binned_flux), &
              binned_veloc_z_ion(nion,n_binned_flux))
     allocate(binned_tau_ion(nion,n_binned_flux))
     allocate(binned_tau_ion_strongest(nion,n_binned_flux))
     if(generate_noise) then
        allocate(binned_noise_sigma(n_binned_flux), binned_noise_random(n_binned_flux))
     endif
     !
     do i=1, n_binned_flux
        binned_lambda(i)     = minlambda + dble(i-1)*pixsize
     enddo
  endif
  !
  if (output_realspacemassweighted_values .or. output_zspaceopticaldepthweighted_values) then
     n_binned_realspace = int((zabsmax - zabsmin) / pixsize * maxlambda) + 1
     !
     if(n_binned_realspace_old /= n_binned_realspace)then
        if(n_binned_realspace_old .gt. 0) then
           write (*,*) ' error: change in spectrum size? '
           stop
        endif
        n_binned_realspace_old = n_binned_realspace
        !
        allocate(binned_redshift_realspace(n_binned_realspace))
        allocate(binned_temp(n_binned_realspace), binned_rho(n_binned_realspace), &
                 binned_met(n_binned_realspace), binned_veloc(n_binned_realspace))
        allocate(binned_temp_ion(nion,n_binned_realspace), &
                 binned_n_ion(nion,n_binned_realspace), &
                 binned_rho_ion(nion,n_binned_realspace), &
                 binned_veloc_ion(nion,n_binned_realspace))
     endif
     !
     do i=1, n_binned_realspace
        binned_redshift_realspace(i)  = zabsmin + dble(i-1)*pixsize/maxlambda
     enddo
  endif
  !
  binned_flux(:)  = 0.d0
  call rebin(flux_convolved,binned_flux)
  !
  if (output_realspacemassweighted_values) then
     binned_temp(:)  = 0.d0
     binned_rho(:)   = 0.d0
     binned_met(:)   = 0.d0
     binned_veloc(:) = 0.d0
     temp_long = temp_long * rho_long
     met_long = met_long * rho_long
     veloc_long = veloc_long * rho_long
     call rebin_realspace(rho_long,binned_rho)
     call rebin_realspace(temp_long,binned_temp)
     call rebin_realspace(met_long,binned_met)
     call rebin_realspace(veloc_long,binned_veloc)
     binned_temp = binned_temp / binned_rho
     binned_met = binned_met / binned_rho
     binned_veloc = binned_veloc / binned_rho
  endif
  !
  if (output_zspaceopticaldepthweighted_values) then
      binned_temp_z_ion  = 0.0
      binned_rho_z_ion   = 0.0
      binned_veloc_z_ion = 0.0
      ! maintain the weighting in the rebinning process
      temp_z_ion_long(:,:) = temp_z_ion_long(:,:) * tau_long_strongest(:,:)
      rho_z_ion_long(:,:) = rho_z_ion_long(:,:) * tau_long_strongest(:,:)
      veloc_z_ion_long(:,:) = veloc_z_ion_long(:,:) * tau_long_strongest(:,:)
      do ion=1,nion
         call rebin(temp_z_ion_long(ion,:),binned_temp_z_ion(ion,:))
         call rebin(rho_z_ion_long(ion,:),binned_rho_z_ion(ion,:))    
         call rebin(veloc_z_ion_long(ion,:),binned_veloc_z_ion(ion,:))
      enddo
  endif
  !
  if (output_realspacenionweighted_values) then
     binned_temp_ion  = 0.0
     binned_rho_ion   = 0.0
     binned_veloc_ion = 0.0
     binned_n_ion     = 0.0
     temp_ion_long(:,:) = temp_ion_long(:,:) * n_ion_long(:,:)
     rho_ion_long(:,:) = rho_ion_long(:,:) * n_ion_long(:,:)
     veloc_ion_long(:,:) = veloc_ion_long(:,:) * n_ion_long(:,:)
     do ion=1,nion
        call rebin_realspace(temp_ion_long(ion,:),binned_temp_ion(ion,:))
        call rebin_realspace(rho_ion_long(ion,:),binned_rho_ion(ion,:))
        call rebin_realspace(veloc_ion_long(ion,:),binned_veloc_ion(ion,:))
        call rebin_realspace(n_ion_long(ion,:),binned_n_ion(ion,:))
     enddo
     binned_temp_ion(:,:) = binned_temp_ion(:,:) / binned_n_ion(:,:)
     binned_rho_ion(:,:) = binned_rho_ion(:,:) / binned_n_ion(:,:)
     binned_veloc_ion(:,:) = binned_veloc_ion(:,:) / binned_n_ion(:,:)
  endif
  !
  binned_tau_ion = 0.0
  binned_tau_ion_strongest = 0.0
  do ion=1,nion
     call rebin(tau_long(ion,:),binned_tau_ion(ion,:))
     call rebin(tau_long_strongest(ion,:),binned_tau_ion_strongest(ion,:))
  enddo
  !
  ! finish weighted averaging for z-space n_ion-weighted spectra
  if (output_zspaceopticaldepthweighted_values) then
     binned_temp_z_ion(:,:) = binned_temp_z_ion(:,:) / binned_tau_ion_strongest(:,:)
     binned_rho_z_ion(:,:) = binned_rho_z_ion(:,:) / binned_tau_ion_strongest(:,:)
     binned_veloc_z_ion(:,:) = binned_veloc_z_ion(:,:) / binned_tau_ion_strongest(:,:)
  endif
  !
end subroutine rebin_spectrum


subroutine convolve_long_spectrum()
  use numbers
  use physical_constants
  use spectra
  use runtime
  use hdf5_wrapper
  use my_mpi
  implicit none
  ! local variables
  real(kind=doubleR) :: sigmakms, b, norm
  integer(kind=singleI) :: i,j
  integer :: file_handle
  real(kind=doubleR), allocatable, save :: gauss(:)
  logical, save :: first_call=.true.
  character(len=120) :: filename
  !
  if (verbose .and. MyPE == 0) then
        write(*,*)
        write(*,'( &
       "Convolving with instrumental profile, fwhm = ", &
       f8.3," km/s...")') fwhm
        write(*,*)
  endif
  !
  if(first_call)then
     first_call = .false.
     !
     !     Compute sigma in km/s
     sigmakms = fwhm / 2.35482
     !     Compute sigma in units of pixels
     b = sigmakms / (vpixsize * LightSpeed / 1d5)
     !     Create normalized Gaussian in wrap-around order
     allocate(gauss(nvpix))
     norm = 1d0 / (2d0 * b * b)
     do i = 0, nvpix-1
        if (i .le. nvpix-1) then 
           if (i .le. nvpix/2) then
              j = i
           else
              j = i - nvpix
           endif
           if (abs(j) .lt. 10.*b) then
              gauss(i+1) = exp(-(dble(j)**2)*norm)
           else
              gauss(i+1) = 0d0
           endif
        else
           gauss(i+1) = 0d0
        endif
     enddo
     !     Make sure Gaussian is normalized.
     gauss  = gauss / sum(gauss)
  endif
  !
  call convlv(flux,nvpix,gauss,nvpix,1,flux_convolved)
  !
end subroutine convolve_long_spectrum


subroutine store_spectrum_info(simfile_used, los_number_used)
  use numbers
  use spectra
  use projection_parameters, only : x_comoving, y_comoving, x_axis, y_axis, z_axis
  implicit none
  character(*), intent(in)          :: simfile_used
  integer(kind=singleI), intent(in) :: los_number_used
  !
  ! store propertes of used sightline
  nsimfile_used = nsimfile_used + 1
  if(nsimfile_used .le. max_nsimfile_used)then
     los_used(nsimfile_used)      = los_number_used
     losfile_used(nsimfile_used)  = simfile_used
     !
     icshift_used(nsimfile_used)  = float(icshift)/float(nveloc)
     !
     x_physical_used(nsimfile_used)    = x_comoving
     y_physical_used(nsimfile_used)    = y_comoving
     x_axis_used(nsimfile_used)        = x_axis
     y_axis_used(nsimfile_used)        = y_axis
     z_axis_used(nsimfile_used)        = z_axis
     ibfactor_used(nsimfile_used) = ibfactor
  endif
end subroutine store_spectrum_info


subroutine add_noise()
  use numbers
  use spectra, only : binned_lambda, binned_flux, binned_noise_sigma, binned_noise_random, n_binned_flux
  use noisedata
  use runtime
  use my_mpi
  use my_random_numbers, only: random
  use random_numbers, only: ran_noise
  !
  implicit none
  ! local variables
  logical, save :: first_call = .true.
  integer(kind=singleI)  :: i, il1, il2, if1, if2, i11, i12, i21, i22
  real(kind=doubleR)     :: dl1, dl2, df1, df2, sigma
  !
  !write(*,*) 'Im here in add noise subroutine' 
  if (verbose .and. MyPE == 0) then 
    if (use_noise_file) then
      write(*,*)"Adding noise from file ",trim(noisefile)
    else 
      write(*,*)"Adding noise with S/N = minnoise + (1/sigtonoise - minnoise) * flux,"
      write(*,*)"where minnoise = ",minnoise," and sigtonoise = ",sigtonoise
    endif
  endif
  !
  do i = 1, n_binned_flux
    ! Interpolate noise level from file.
    ! Bilinear interpolation (noise is binned in binned_lambda and binned_flux).
    if (use_noise_file) then 
      il2 = 2
      do while (binned_lambda(i) .gt. n_lambda(il2))
        il2 = il2 + 1
      enddo
      il1 = il2 - 1
      dl1 = (n_lambda(il2) - binned_lambda(i)) / (n_lambda(il2) - n_lambda(il1))
      dl2 = 1. - dl1   
      if2 = 2
      do while (binned_flux(i) .gt. n_flux(if2))
        if2 = if2 + 1
      enddo
      if1 = if2 - 1
      df1 = (n_flux(if2) - binned_flux(i)) / (n_flux(if2) - n_flux(if1))
      df2 = 1. - df1   
      sigma = dl1*df1*n_sigma(il1,if1) + dl2*df1*n_sigma(il2,if1) + &
        dl1*df2*n_sigma(il1,if2) + dl2*df2*n_sigma(il2,if2)   
    else
      ! Use noise independent of wavelength.
      sigma = minnoise + (1d0 - minnoise*sigtonoise) /  sigtonoise * binned_flux(i)
    endif
    ! Actual noise realisation requires multiplying with Gaussian deviate with mean 0 and dispersion=1
    ! Array binned_noise stores sigma of noise, not actual noise
    ! Dispersion of Gaussian noise	
    binned_noise_sigma(i)  = sigma
    ! Gaussian deviate with mean 0 and dispersion 1
    binned_noise_random(i) = random(ran_noise)
  enddo
  !
end subroutine add_noise


! +++++++++++++++++++++++++++++++++++ ...do_long_spectra +++++++++++++++++++++++++++++++++++++++++ !

! ++++++++++++++++++++++++++++++++++++++++ read/write... +++++++++++++++++++++++++++++++++++++++++ !


subroutine read_header(file_handle)
  use numbers
  use header
  use runtime, only: use_snapshot_file
  USE physical_constants
  use hdf5_wrapper
  use particledata, only: rhocb
  implicit none
  integer, intent(in) ::  file_handle
  !
  call hdf5_read_attribute(file_handle,'Header/NumPart_ThisFile',NumPart_ThisFile)
  call hdf5_read_attribute(file_handle,'Header/NumPart_Total',NumPart_Total)
  call hdf5_read_attribute(file_handle,'Header/NumPart_Total_HighWord',NumPart_Total_HighWord)
  call hdf5_read_attribute(file_handle,'Header/NumFilesPerSnapshot',NumFilesPerSnapshot)
  call hdf5_read_attribute(file_handle,'Header/Flag_Sfr',Flag_Sfr)
  call hdf5_read_attribute(file_handle,'Header/Flag_Cooling',Flag_Cooling)
  call hdf5_read_attribute(file_handle,'Header/Flag_StellarAge',Flag_StellarAge)
  call hdf5_read_attribute(file_handle,'Header/Flag_Metals',Flag_Metals)
  call hdf5_read_attribute(file_handle,'Header/Flag_Feedback',Flag_Feedback)
  call hdf5_read_attribute(file_handle,'Header/MassTable',MassTable)
  call hdf5_read_attribute(file_handle,'Header/ExpansionFactor',ExpansionFactor)
  call hdf5_read_attribute(file_handle,'Header/Redshift',Redshift)
  call hdf5_read_attribute(file_handle,'Header/BoxSize',BoxSize)
  call hdf5_read_attribute(file_handle,'Header/Omega0',Omega0)
  call hdf5_read_attribute(file_handle,'Header/OmegaBaryon',OmegaBaryon)
  call hdf5_read_attribute(file_handle,'Header/OmegaLambda',OmegaLambda)
  call hdf5_read_attribute(file_handle,'Header/HubbleParam',HubbleParam)
#ifdef EAGLE
  Time =-1.
#else
  call hdf5_read_attribute(file_handle,'Header/Time_GYR',Time)
#endif
  !
  if (.not. use_snapshot_file) then
    call hdf5_read_attribute(file_handle,'Header/Number_of_sight_lines', los_this_file)
  else
    los_this_file = -1
  endif
  !
  rhoc  = 3 * (H0*HubbleParam)**2 / (8. * pi * G) ! Critical density g/cm^3 
  rhocb = rhoc / ExpansionFactor**3 * OmegaBaryon
  !
end subroutine read_header


subroutine check_header() ! check whether parameters are realistic
  use header
  implicit none
  !
  if(BoxSize .LT. 1e-5) then 
    write (*,*) ' suspicious box size: ',BoxSize
    stop
  endif
  !
  if(HubbleParam .lt. 0.1) then
     write (*,*) ' suscipicious HubbleParam: ',HubbleParam
     stop
  endif
  !
  if(Omega0 .LT. 1e-5) then
    write (*,*) ' suscipicious Omega0: ',Omega0
    stop
  endif
  !
  if(OmegaLambda .LT. 1e-5) then
    write (*,*) ' suspicious OmegaLambda: ',OmegaLambda
    stop
  endif
  !
end subroutine check_header


subroutine write_header(file_handle)
  use numbers
  use header
  use hdf5_wrapper
  implicit none
  integer, intent(in):: file_handle
  !
  call hdf5_create_group(file_handle,"/Header")
  !
  call hdf5_write_attribute(file_handle,'Header/NumPart_ThisFile',NumPart_ThisFile)
  call hdf5_write_attribute(file_handle,'Header/NumPart_Total',NumPart_Total)
  call hdf5_write_attribute(file_handle,'Header/NumPart_Total_HighWord',NumPart_Total_HighWord)
  call hdf5_write_attribute(file_handle,'Header/NumFilesPerSnapshot',NumFilesPerSnapshot)
  call hdf5_write_attribute(file_handle,'Header/Flag_Sfr',Flag_Sfr)
  call hdf5_write_attribute(file_handle,'Header/Flag_Cooling',Flag_Cooling)
  call hdf5_write_attribute(file_handle,'Header/Flag_StellarAge',Flag_StellarAge)
  call hdf5_write_attribute(file_handle,'Header/Flag_Metals',Flag_Metals)
  call hdf5_write_attribute(file_handle,'Header/Flag_Feedback',Flag_Feedback)
  call hdf5_write_attribute(file_handle,'Header/MassTable',MassTable)
  call hdf5_write_attribute(file_handle,'Header/ExpansionFactor',ExpansionFactor)
  call hdf5_write_attribute(file_handle,'Header/Redshift',Redshift)
  call hdf5_write_attribute(file_handle,'Header/BoxSize',BoxSize)
  call hdf5_write_attribute(file_handle,'Header/Omega0',Omega0)
  call hdf5_write_attribute(file_handle,'Header/OmegaBaryon',OmegaBaryon)
  call hdf5_write_attribute(file_handle,'Header/OmegaLambda',OmegaLambda)
  call hdf5_write_attribute(file_handle,'Header/HubbleParam',HubbleParam)
  call hdf5_write_attribute(file_handle,'Header/Time',Time)
  !
end subroutine write_header

#ifdef EAGLE
subroutine read_parameters(file_handle)
  use numbers
  use runtime
  use parameters
  use hdf5_wrapper
  implicit none
  !
  integer, intent(in) :: file_handle
  character(len=120)  :: GroupName
  !
  GroupName = 'Parameters/ChemicalElements'
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/BG_NELEMENTS',BG_Nelements)
  if(allocated(ElementNames)) deallocate(ElementNames)
  allocate(ElementNames(BG_Nelements))
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/ElementNames',ElementNames)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Hydrogen',InitAbundance_Hydrogen)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Helium',InitAbundance_Helium)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Carbon',InitAbundance_Carbon)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Nitrogen',InitAbundance_Nitrogen)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Oxygen',InitAbundance_Oxygen)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Neon',InitAbundance_Neon)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Magnesium',InitAbundance_Magnesium)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Silicon',InitAbundance_Silicon)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Iron',InitAbundance_Iron)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/CalciumOverSilicon',CalciumOverSilicon)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SulphurOverSilicon',SulphurOverSilicon)
  !
  !
end subroutine read_parameters
subroutine write_parameters(file_handle)
  use numbers
  use runtime
  use parameters
  use hdf5_wrapper
  implicit none
  !
  integer, intent(in):: file_handle
  character(len=120) :: GroupName
  !
  call hdf5_create_group(file_handle,'/Parameters')
  !
  GroupName = 'Parameters/ChemicalElements'
  call hdf5_create_group(file_handle,trim(GroupName))
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/BG_NELEMENTS',BG_Nelements)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ElementNames',ElementNames)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Hydrogen',InitAbundance_Hydrogen)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Helium',InitAbundance_Helium)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Carbon',InitAbundance_Carbon)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Nitrogen',InitAbundance_Nitrogen)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Oxygen',InitAbundance_Oxygen)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Neon',InitAbundance_Neon)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Magnesium',InitAbundance_Magnesium)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Silicon',InitAbundance_Silicon)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Iron',InitAbundance_Iron)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/CalciumOverSilicon',CalciumOverSilicon)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SulphurOverSilicon',SulphurOverSilicon)
  !
  !
end subroutine write_parameters

#else
subroutine read_parameters(file_handle)
  use numbers
  use runtime
  use parameters
  use hdf5_wrapper
  implicit none
  !
  integer, intent(in) :: file_handle
  character(len=120)  :: GroupName
  !
  GroupName = 'Parameters/ChemicalElements'
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/BG_NELEMENTS',BG_Nelements)
  if(allocated(ElementNames)) deallocate(ElementNames)
  allocate(ElementNames(BG_Nelements))
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/ElementNames',ElementNames)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Hydrogen',InitAbundance_Hydrogen)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Helium',InitAbundance_Helium)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Carbon',InitAbundance_Carbon)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Nitrogen',InitAbundance_Nitrogen)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Oxygen',InitAbundance_Oxygen)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Neon',InitAbundance_Neon)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Magnesium',InitAbundance_Magnesium)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Silicon',InitAbundance_Silicon)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitAbundance_Iron',InitAbundance_Iron)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/CalciumOverSilicon',CalciumOverSilicon)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SulphurOverSilicon',SulphurOverSilicon)
  !
  GroupName = 'Parameters/StellarEvolutionParameters'
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNIa_Model',SNIa_Model)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/IMF_Model',IMF_model)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/IMF_LifetimeModel',IMF_LifetimeModel)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SF_EOSEnergyAtThreshold_ERG',SF_EOSEnergyAtThreshold_ERG)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SF_EOSGammaEffective',SF_EOSGammaEffective)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SF_THRESH_MinPhysDens_HpCM3',SF_THRESH_MinPhysDens_HpCM3)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SF_THRESH_MinOverDens',SF_THRESH_MinOverDens)
  if(.not. wmap7) then ! Some attributes not present in OWLS WMAP7 runs
    SF_THRESH_MetDepExponent = 0.0 !This variable is present only in the old version of the OWLS code, it is never used in SpecWizard
    call hdf5_read_attribute(file_handle,trim(GroupName)//'/SF_THRESH_MetDepFiducialZ',SF_THRESH_MetDepFiducialZ)
    call hdf5_read_attribute(file_handle,trim(GroupName)//'/SF_THRESH_MetDepMaxPhysDens_HpCM3',SF_THRESH_MetDepMaxPhysDens_HpCM3)
  endif
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SF_THRESH_MaxTemp_K',SF_THRESH_MaxTemp_K)         
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SF_SchmidtLawCoeff_MSUNpYRpKPC2',SF_SchmidtLawCoeff_MSUNpYRpKPC2)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SF_SchmidtLawExponent',SF_SchmidtLawExponent)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SF_SchmidtLawCoeff_GpSpCM2',SF_SchmidtLawCoeff_GpSpCM2)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/IMF_MinMass_MSUN',IMF_MinMass_MSUN)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNIa_Efficiency_fracwd',SNIa_Efficiency_fracwd)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNIa_MassTransferOn',SNIa_MassTransferOn)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNIa_EnergyTransferOn',SNIa_EnergyTransferOn)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_MinMass_MSUN',SNII_MinMass_MSUN)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_MaxMass_MSUN',SNII_MaxMass_MSUN)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Hydrogen',SNII_Factor_Hydrogen)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Helium',SNII_Factor_Helium)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Carbon',SNII_Factor_Carbon)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Nitrogen',SNII_Factor_Nitrogen)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Oxygen',SNII_Factor_Oxygen)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Neon',SNII_Factor_Neon)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Silicon',SNII_Factor_Silicon)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Iron',SNII_Factor_Iron)
  if(.not. wmap7) then ! Some attributes not present in OWLS WMAP7 runs
    call hdf5_read_attribute(file_handle,trim(GroupName)//'/POPIII_Energy_ERG',POPIII_Energy_ERG)
    call hdf5_read_attribute(file_handle,trim(GroupName)//'/POPIII_NumPerMsun',POPIII_NumPerMsun)
    call hdf5_read_attribute(file_handle,trim(GroupName)//'/POPIII_MassTransferOn',POPIII_MassTransferOn)
    call hdf5_read_attribute(file_handle,trim(GroupName)//'/POPIII_EnergyTransferOn',POPIII_EnergyTransferOn)
  endif
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/AGB_MassTransferOn',AGB_MassTransferOn)
  !
  ! Wind parameters
  GroupName = 'Parameters/WindParameters'
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_WindSpeed_KMpS',SNII_WindSpeed_KMpS)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_WindOn',SNII_WindOn)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_WindIsotropicOn',SNII_WindIsotropicOn)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_WindMassLoading',SNII_WindMassLoading)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/SNII_WindDelay_YR',SNII_WindDelay_YR)
  !
  GroupName = 'Parameters/NumericalParameters'
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/ComovingIntegrationOn',ComovingIntegrationOn)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/TypeOfTimestepCriterion',TypeOfTimestepCriterion)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/TypeOfOpeningCriterion',TypeOfOpeningCriterion)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/DesNumNgb',DesNumNgb)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/MaxNumNgbDeviation',MaxNumNgbDeviation)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/BufferSize',BufferSize)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/ErrTolIntAccuracy',ErrTolIntAccuracy)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/MaxSizeTimestep',MaxSizeTimestep)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/MinSizeTimestep',MinSizeTimestep)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/ErrTolTheta',ErrTolTheta)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/ErrTolForceAcc',ErrTolForceAcc)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/TreeDomainUpdateFrequency',TreeDomainUpdateFrequency)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/ArtBulkViscConst',ArtBulkViscConst)                         
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/InitGasU_ERG',InitGasU_ERG)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/MinGasU_ERG', MinGasU_ERG)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/CourantFac',CourantFac)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/PartAllocFactor',PartAllocFactor)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/TreeAllocFactor',TreeAllocFactor)
  call hdf5_read_attribute(file_handle,trim(GroupName)//'/MinGasHsmlFractional',MinGasHsmlFractional)
  !
end subroutine read_parameters
subroutine write_parameters(file_handle)
  use numbers
  use runtime
  use parameters
  use hdf5_wrapper
  implicit none
  !
  integer, intent(in):: file_handle
  character(len=120) :: GroupName
  !
  call hdf5_create_group(file_handle,'/Parameters')
  !
  GroupName = 'Parameters/ChemicalElements'
  call hdf5_create_group(file_handle,trim(GroupName))
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/BG_NELEMENTS',BG_Nelements)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ElementNames',ElementNames)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Hydrogen',InitAbundance_Hydrogen)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Helium',InitAbundance_Helium)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Carbon',InitAbundance_Carbon)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Nitrogen',InitAbundance_Nitrogen)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Oxygen',InitAbundance_Oxygen)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Neon',InitAbundance_Neon)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Magnesium',InitAbundance_Magnesium)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Silicon',InitAbundance_Silicon)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitAbundance_Iron',InitAbundance_Iron)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/CalciumOverSilicon',CalciumOverSilicon)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SulphurOverSilicon',SulphurOverSilicon)
  !
  GroupName = 'Parameters/StellarEvolutionParameters'
  call hdf5_create_group(file_handle,trim(GroupName))
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNIa_Model',SNIa_Model)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/IMF_Model',IMF_model)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/IMF_LifetimeModel',IMF_LifetimeModel)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SF_EOSEnergyAtThreshold_ERG',SF_EOSEnergyAtThreshold_ERG)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SF_EOSGammaEffective',SF_EOSGammaEffective)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SF_THRESH_MinPhysDens_HpCM3',SF_THRESH_MinPhysDens_HpCM3)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SF_THRESH_MinOverDens',SF_THRESH_MinOverDens)
  if(.not. wmap7) then ! Some attributes not present in OWLS WMAP7 runs
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/SF_THRESH_MetDepExponent',SF_THRESH_MetDepExponent)
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/SF_THRESH_MetDepFiducialZ',SF_THRESH_MetDepFiducialZ)
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/SF_THRESH_MetDepMaxPhysDens_HpCM3',SF_THRESH_MetDepMaxPhysDens_HpCM3)
  endif
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SF_THRESH_MaxTemp_K',SF_THRESH_MaxTemp_K)         
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SF_SchmidtLawCoeff_MSUNpYRpKPC2',SF_SchmidtLawCoeff_MSUNpYRpKPC2)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SF_SchmidtLawExponent',SF_SchmidtLawExponent)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SF_SchmidtLawCoeff_GpSpCM2',SF_SchmidtLawCoeff_GpSpCM2)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/IMF_MinMass_MSUN',IMF_MinMass_MSUN)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNIa_Efficiency_fracwd',SNIa_Efficiency_fracwd)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNIa_MassTransferOn',SNIa_MassTransferOn)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNIa_EnergyTransferOn',SNIa_EnergyTransferOn)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_MinMass_MSUN',SNII_MinMass_MSUN)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_MaxMass_MSUN',SNII_MaxMass_MSUN)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Hydrogen',SNII_Factor_Hydrogen)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Helium',SNII_Factor_Helium)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Carbon',SNII_Factor_Carbon)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Nitrogen',SNII_Factor_Nitrogen)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Oxygen',SNII_Factor_Oxygen)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Neon',SNII_Factor_Neon)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Silicon',SNII_Factor_Silicon)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_Factor_Iron',SNII_Factor_Iron)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/AGB_MassTransferOn',AGB_MassTransferOn)
  if(.not. wmap7) then ! Some attributes not present in WMAP7 runs
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/POPIII_Energy_ERG',POPIII_Energy_ERG)
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/POPIII_NumPerMsun',POPIII_NumPerMsun)
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/POPIII_MassTransferOn',POPIII_MassTransferOn)
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/POPIII_EnergyTransferOn',POPIII_EnergyTransferOn)
  endif
  ! Wind parameters
  GroupName = 'Parameters/WindParameters'
  call hdf5_create_group(file_handle,trim(GroupName))
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_WindSpeed_KMpS',SNII_WindSpeed_KMpS)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_WindOn',SNII_WindOn)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_WindIsotropicOn',SNII_WindIsotropicOn)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_WindMassLoading',SNII_WindMassLoading)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SNII_WindDelay_YR',SNII_WindDelay_YR)
  !
  ! Numerical Parameters
  GroupName = 'Parameters/NumericalParameters'
  call hdf5_create_group(file_handle,trim(GroupName))
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ComovingIntegrationOn',ComovingIntegrationOn)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/TypeOfTimestepCriterion',TypeOfTimestepCriterion)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/TypeOfOpeningCriterion',TypeOfOpeningCriterion)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/DesNumNgb',DesNumNgb)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/MaxNumNgbDeviation',MaxNumNgbDeviation)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/BufferSize',BufferSize)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ErrTolIntAccuracy',ErrTolIntAccuracy)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/MaxSizeTimestep',MaxSizeTimestep)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/MinSizeTimestep',MinSizeTimestep)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ErrTolTheta',ErrTolTheta)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ErrTolForceAcc',ErrTolForceAcc)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/TreeDomainUpdateFrequency',TreeDomainUpdateFrequency)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ArtBulkViscConst',ArtBulkViscConst)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/InitGasU_ERG',InitGasU_ERG)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/MinGasU_ERG', MinGasU_ERG)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/CourantFac',CourantFac)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/PartAllocFactor',PartAllocFactor)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/TreeAllocFactor',TreeAllocFactor)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/MinGasHsmlFractional',MinGasHsmlFractional)
  !
end subroutine write_parameters
#endif

subroutine read_constants(file_handle)
  use numbers
  use constants
  use hdf5_wrapper
  implicit none
  !
  integer,intent(in):: file_handle
!!$  call hdf5_read_attribute(file_handle,'Constants/PI',pi)
  call hdf5_read_attribute(file_handle,'Constants/GAMMA',gamma)
  call hdf5_read_attribute(file_handle,'Constants/GRAVITY',gravity)
  call hdf5_read_attribute(file_handle,'Constants/SOLAR_MASS',solar_mass)
  call hdf5_read_attribute(file_handle,'Constants/SOLAR_LUM',solar_lum)
  call hdf5_read_attribute(file_handle,'Constants/BOLTZMANN',boltzmann)
  call hdf5_read_attribute(file_handle,'Constants/GAS_CONST',gas_const)
  call hdf5_read_attribute(file_handle,'Constants/CM_PER_MPC',cm_per_mpc)
  !
end subroutine read_constants


subroutine write_constants(file_handle)
  use numbers
  use constants
  use physical_constants, only : pi
  use hdf5_wrapper
  implicit none
  !
  integer,intent(in):: file_handle
  !
  call hdf5_create_group(file_handle, '/Constants')
  call hdf5_write_attribute(file_handle,'Constants/PI',pi)
  call hdf5_write_attribute(file_handle,'Constants/GAMMA',gamma)
  call hdf5_write_attribute(file_handle,'Constants/GRAVITY',gravity)
  call hdf5_write_attribute(file_handle,'Constants/SOLAR_MASS',solar_mass)
  call hdf5_write_attribute(file_handle,'Constants/SOLAR_LUM',solar_lum)
  call hdf5_write_attribute(file_handle,'Constants/BOLTZMANN',boltzmann)
  call hdf5_write_attribute(file_handle,'Constants/GAS_CONST',gas_const)
  call hdf5_write_attribute(file_handle,'Constants/CM_PER_MPC',cm_per_mpc)
  !
end subroutine write_constants


subroutine read_units(file_handle)
  use numbers
  use units
  use hdf5_wrapper
  implicit none
  !
  integer,intent(in):: file_handle
  !
  call hdf5_read_attribute(file_handle,'Units/UnitLength_in_cm',UnitLength_in_cm)
  call hdf5_read_attribute(file_handle,'Units/UnitMass_in_g',UnitMass_in_g)
  call hdf5_read_attribute(file_handle,'Units/UnitVelocity_in_cm_per_s',UnitVelocity_in_cm_per_s)
  call hdf5_read_attribute(file_handle,'Units/UnitDensity_in_cgs',UnitDensity_in_cgs)
  call hdf5_read_attribute(file_handle,'Units/UnitEnergy_in_cgs',UnitEnergy_in_cgs)
  !
end subroutine read_units


subroutine write_units(file_handle)
  use numbers
  use units
  use hdf5_wrapper
  implicit none
  !
  integer, intent(in):: file_handle
  !
  call hdf5_create_group(file_handle,"/Units")
  call hdf5_write_attribute(file_handle,'Units/UnitLength_in_cm',UnitLength_in_cm)
  call hdf5_write_attribute(file_handle,'Units/UnitMass_in_g',UnitMass_in_g)
  call hdf5_write_attribute(file_handle,'Units/UnitVelocity_in_cm_per_s',UnitVelocity_in_cm_per_s)
  call hdf5_write_attribute(file_handle,'Units/UnitDensity_in_cgs',UnitDensity_in_cgs)
  call hdf5_write_attribute(file_handle,'Units/UnitEnergy_in_cgs',UnitEnergy_in_cgs)
  !
end subroutine write_units


subroutine write_specwizard_runtime_parameters(file_handle)
  use parameters
  use runtime
  use spectra
  use noisedata
  use modified_metallicity
  use hdf5_wrapper
  implicit none
  !
  integer, intent(in):: file_handle
  character(len=120) :: GroupName 
  integer(kind=singleI)           :: number_of_transitions, line, count, ion
  character(len=10), allocatable  :: all_ions(:)
  real(kind=doubleR), allocatable :: all_lambda_rest(:), all_fosc(:), all_ion_mass(:)
  !
  ! total number of transitions included
  number_of_transitions = sum(nlines)
  allocate(all_ions(number_of_transitions),all_lambda_rest(number_of_transitions),all_fosc(number_of_transitions),&
       all_ion_mass(number_of_transitions))
  !
  count = 0
  do ion=1, nion
     do line=1, nlines(ion)
        count = count + 1
        all_ions(count) = trim(ions(ion))
        all_lambda_rest(count) = lambda_rest(ion,line)
        all_fosc(count)        = fosc(ion,line)
        all_ion_mass(count)    = ion_mass(ion)
     enddo
  enddo
  !
  GroupName = 'Parameters/SpecWizardRuntimeParameters'
  call hdf5_create_group(file_handle, trim(GroupName))
  !
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/datadir',datadir)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/outputdir',outputdir)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ibdir',ibdir)
  if(do_long_spectrum) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/do_long_spectrum','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/do_long_spectrum','FALSE')
  endif
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/SpectrumFile',SpectrumFile)
  if(generate_noise) then
     call hdf5_write_attribute(file_handle,trim(GroupName)//'/generate_noise','TRUE')
     if(use_noise_file)then
        call hdf5_write_attribute(file_handle,trim(GroupName)//'/use_noise_file','TRUE')
        call hdf5_write_attribute(file_handle,trim(GroupName)//'/noise_file',noisefile)
     else
        call hdf5_write_attribute(file_handle,trim(GroupName)//'/signal_to_noise',sigtonoise)
     endif

  else
     call hdf5_write_attribute(file_handle,trim(GroupName)//'/generate_noise','FALSE')
  endif
  if(use_snapshot_file) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/use_snapshot_file','TRUE')
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/snap',snap)
    if(use_random_los) then
      call hdf5_write_attribute(file_handle,trim(GroupName)//'/use_random_los','TRUE')
    else
      call hdf5_write_attribute(file_handle,trim(GroupName)//'/use_random_los','FALSE')
      call hdf5_write_attribute(file_handle,trim(GroupName)//'/los_coordinates_file',los_coordinates_file)
    endif
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/use_snapshot_file','FALSE')
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/file_list',file_list)
  endif
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/NumberOfSpectra',nspec)
  if(use_maxdens_above_zmax) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/use_maxdens_above_zmax','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/use_maxdens_above_zmax','FALSE')
  endif
  if(gimic) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/gimic','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/gimic','FALSE')
  endif
  if(urchin) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/urchin','TRUE')
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/urchindir',urchindir)
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/urchin','FALSE')
  endif
  if(do_periodic) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/do_periodic','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/do_periodic','FALSE')
  endif
  if(use_urchin_temperature) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/use_urchin_temperature','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/use_urchin_temperature','FALSE')
  endif
  if(use_smoothed_abundance) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/use_smoothed_abundance','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/use_smoothed_abundance','FALSE')
  endif
  if(ignore_starforming) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/ignore_starforming','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/ignore_starforming','FALSE')
  endif
  if(setmaxt4sfgas) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/setmaxt4sfgas','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/setmaxt4sfgas','FALSE')
  endif
  if(subtract_Hmol) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/subtract_Hmol','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/subtract_Hmol','FALSE')
  endif
  if(use_gaussian_kernel) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/use_gaussian_kernel','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/use_gaussian_kernel','FALSE')
  endif
  if(integrate_kernel) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/integrate_kernel','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/integrate_kernel','FALSE')
  endif
  if(wmap7) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/wmap7','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/wmap7','FALSE')
  endif
  if(read_ionbal_from_single_file) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/read_ionbal_from_single_file','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/read_ionbal_from_single_file','FALSE')
  endif
  if(output_frequency) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/output_frequency','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/output_frequency','FALSE')
  endif
  if(NoPecVel) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/NoPecVel','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/NoPecVel','FALSE')
  endif
  if(ionfracone) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/ionfracone','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/ionfracone','FALSE')
  endif
  if(docycling) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/docycling','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/docycling','FALSE')
  endif
  if(impose_eos) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/impose_eos','TRUE')
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/imposed_eos_T0',imposed_eos_T0)
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/imposed_eos_gamma',imposed_eos_gamma)
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/imposed_eos_maxod',imposed_eos_maxod)
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/impose_eos','FALSE')
  endif
  if(add_turbulence) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/add_turbulence','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/add_turbulence','FALSE')
  endif

  if(output_realspacenionweighted_values) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/output_realspacenionweighted_values','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/output_realspacenionweighted_values','FALSE')
  endif
  if(output_realspacemassweighted_values) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/output_realspacemassweighted_values','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/output_realspacemassweighted_values','FALSE')
  endif
  if(output_zspaceopticaldepthweighted_values) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/output_zspaceopticaldepthweighted_values','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/output_zspaceopticaldepthweighted_values','FALSE')
  endif
  if(do_convolve_spectrum) then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/do_convolve_spectrum','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/do_convolve_spectrum','FALSE')
  endif
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/small_rho',small_rho)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/small_temp',small_temp)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/small_metallicity',small_metallicity)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/small_velocity',small_velocity)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/small_column',small_column )
  !
  if(modify_metallicity)then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/modify_metallicity','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/modify_metallicity','FALSE')
  endif
  if(read_part_ids_from_file)then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/read_part_ids_from_file','TRUE')
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/particle_file_name',particle_file_name)
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/flagged_particle_metallicity',flagged_particle_metallicity)
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/read_part_ids_from_file','FALSE')
  endif
  if(scale_simulation_abundances)then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/scale_simulation_abundances','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/scale_simulation_abundances','FALSE')
  endif
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/z_rel',z_rel)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ZC_rel',ZC_rel)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ZN_rel',ZN_rel)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ZO_rel',ZO_rel)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ZNe_rel',ZNe_rel)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ZMg_rel',ZMg_rel)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ZAl_rel',ZAl_rel)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ZSi_rel',ZSi_rel)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ZS_rel',ZS_rel)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/ZFe_rel',ZFe_rel)
  if(impose_z_rho_relation)then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/impose_z_rho_relation','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/impose_z_rho_relation','FALSE')
  endif
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/z_index',z_index)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/z_mean',z_mean)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/maxz_rel',maxz_rel)
  if(log_normal_scatter)then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/log_normal_scatter','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/log_normal_scatter','FALSE')
  endif
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/z_sig_bin',z_sig_bin)
  call hdf5_write_attribute(file_handle,trim(GroupName)//'/z_sig_dex',z_sig_dex)
  !
  GroupName = 'Header'
  !
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/NumberOfSpectra', nspec)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Z_min', zabsmin)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Z_max', zabsmax)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Lambda_min', minlambda)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Lambda_max', maxlambda)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Z_qso', Zqso)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Z_resol', fzresol)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/PixSize_Angstrom', PiXSize)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/FWHM_kms', FWHM)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/PixSizekms_Before_convolution', vpixsizekms)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/nLyman', nLyman)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/ibfactor',ibfactor)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/minbother_red', minbother_red)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/minbother_blue', minbother_blue)
  if(use_fitted_ibfactor) then
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/use_fitted_ibfactor', 'TRUE')
  else
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/use_fitted_ibfactor', 'FALSE')
  endif
  if(ibfactor_he_reionization)then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/ibfactor_he_reionization','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/ibfactor_he_reionization','FALSE')
  endif
  if(limsigma)then
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/limsigma','TRUE')
  else
    call hdf5_write_attribute(file_handle,trim(GroupName)//'/limsigma','FALSE')
  endif
  if(integrate_thermprof_exactly)then
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/Integrate_thermprof_exactly', 'TRUE')
  else
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/Integrate_thermprof_exactly', 'FALSE')
  endif
  if(docycling)then
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/Cycle_spectra', 'TRUE')
  else
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/Cycle_spectra', 'FALSE')
  endif
  if(NoPecVel)then
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/NoPeculiarVelocities', 'TRUE')
  else
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/NoPeculiarVelocities', 'FALSE')
  endif
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/SmallRho', Small_Rho)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/SmallTemp', Small_Temp)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Small_metallicity',Small_metallicity )
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Small_velocity',Small_velocity )
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Small_column', Small_Column)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Datadir', datadir)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Ionization_directory', ibdir)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/NumberOfTransitions', nion)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Ions', all_ions)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Ion_Masses', all_ion_mass)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Transitions_Rest_Wavelength', all_lambda_rest)
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/Transitions_Oscillator_Strength', all_fosc)
  !
  GroupName = 'Header/ModifyMetallicityParameters'
  call hdf5_create_group(file_handle, trim(GroupName))
  !
  if(modify_metallicity)then
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/Imposed_Metallicity','TRUE.')
  else
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/Imposed_Metallicity','FALSE.')
  endif
  if(Modify_metallicity)then
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/Modify_metallicity','True')
  else
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/Modify_metallicity','False')
  endif
  if(scale_simulation_abundances)then
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/scale_simulation_abundances','True')
  else
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/scale_simulation_abundances','False')
  endif
  call hdf5_write_attribute(file_handle, trim(GroupName)//'/z_rel',z_rel)
  if(impose_z_rho_relation)then
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/impose_z_rho_relation','True')
  else
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/impose_z_rho_relation','False')
  endif
  if(log_normal_scatter)then
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/log_normal_scatter','True')
  else
     call hdf5_write_attribute(file_handle, trim(GroupName)//'/log_normal_scatter','False')
  endif
  !
  deallocate(all_ions,all_lambda_rest,all_fosc,all_ion_mass)
  !
end subroutine write_specwizard_runtime_parameters


subroutine load_noise
  !
  ! Reads in noise file interpolation table.
  ! Table is 2-D : lambda, flux
  use hdf5_wrapper
  use spectra
  use noisedata
  implicit none
  integer :: file_handle, il, if, rank, array_size(7)
  !
  call hdf5_open_file(file_handle, noisefile,readonly=.true.)
  !
  call hdf5_get_dimensions(file_handle,"NormalizedFlux",rank,array_size)
  if (rank .ne. 1) call abortrun('NormalizedFlux should have rank 1')
  n_nf = array_size(1)
  !
  call hdf5_get_dimensions(file_handle,"Wavelength_Angstrom",rank,array_size)
  if (rank .ne. 1) call abortrun('NormalizedFlux should have rank 1')
  n_nl = array_size(1)
  !
  allocate(n_lambda(n_nl))
  allocate(n_flux(n_nf))
  allocate(n_sigma(n_nl,n_nf))
  !
  call hdf5_read_data(file_handle,'NormalizedFlux',n_flux)  
  call hdf5_read_data(file_handle,'Wavelength_Angstrom',n_lambda)  
  call hdf5_read_data(file_handle,'NormalizedNoise',n_sigma)
  call hdf5_close_file(file_handle)
  !
end subroutine load_noise



subroutine readdata_owls(filename,los_number)
  !
  ! Read in data of OWLS simulation (hdf5): 
  ! modified -> read EAGLE or OWLS data from LOS files
  !
  use numbers
  use header
  use hdf5_wrapper
  use projection_parameters
  use particledata
  use runtime
  use physical_constants
  use atomic_data
  use constants
  use solar_data, only : zmass_solar
  use spectra
  use ionization_tables, only : nz, nt, nd, ib_logt, ib_logd
  use my_mpi
  implicit none
  !
  character(*),intent(in) :: filename
  integer, intent(in)     :: los_number
  !
  ! local variables
  character(len=300)       :: basefile, longfile
  integer                  :: file_handle
  character(len=120)       :: LosNumber, LosName, VarName
  character(len=3)         :: FileNumber
  integer                  :: nfiles
  real(kind=doubleR)       :: Pos_h_exp, Pos_aexp_exp, Vel_h_exp, Vel_aexp_exp, &
    Dens_h_exp, Dens_aexp_exp, Eint_h_exp, Eint_aexp_exp, &
    Mass_h_exp, Mass_aexp_exp,  Temp_h_exp, Temp_aexp_exp
  real(kind=doubleR)       :: Mass_cgs_unit, Dens_cgs_unit, Eint_cgs_unit, Temp_cgs_unit, Pos_cgs_unit, &
    Vel_cgs_unit
  !real(kind=doubleR)       :: proton_mass
  real(kind=singleR)       :: z_solar
  real(kind=doubleR)       :: Coordinates_conv, Velocity_conv, Density_conv, Mass_conv, Eint_conv, Temp_conv
  real(kind=doubleR)       :: tmin, tmax, dmin, dmax, metal_max, rnorm2propermpc
  integer(kind=singleI)    :: i
  logical                  :: file_exists, single_file
  logical, save            :: first_call = .true.
  integer                  :: files, NTotal, NPart_This_Sight
  real(kind=doubleR), allocatable :: tmp(:,:)
  !
  longfile = trim(datadir)//'/'//trim(filename)
  inquire(file=longfile,exist=single_file)
  if(.not. single_file)then
    ! multi-file format
    basefile = trim(longfile)
  endif
  !
  !!$ write (0,*) ' file = ',trim(longfile)
  !
  if(single_file) then
    NFiles = 1
    call hdf5_open_file(file_handle, longfile, readonly=.true.)
    call read_header(file_handle)
    call hdf5_close_file(file_handle)
  else
    NFiles = NumFilesPerSnapshot
    NGas   = 0
    ! determine total number of Gas particles
    do files=0, NFiles-1
      write(FileNumber,'(i3)') files
      longfile = trim(basefile)//'.'//trim(adjustl(FileNumber))//'.los.hdf5'        
      write (LosNumber,'(I4)') abs(los_number)
      LOSName = "LOS"//trim(adjustl(LosNumber))
      call hdf5_open_file(file_handle, longfile, readonly=.true.)
      call read_header(file_handle)
      VarName = trim(LosName)//'/Npart_this_sight'
      call hdf5_read_attribute(file_handle, VarName, NPart_This_sight)
      !
      if(files == 1) then
        call hdf5_read_attribute(file_handle,'/Constants/CM_PER_MPC',cm_per_mpc) 
        call hdf5_read_attribute(file_handle,'/Constants/PROTONMASS',proton_mass) 
        call hdf5_read_attribute(file_handle,'/Constants/SOLAR_MASS',solar_mass) 
        call hdf5_read_attribute(file_handle,'/Header/ExpansionFactor', ExpansionFactor)
        call hdf5_read_attribute(file_handle,'/Header/Redshift', Redshift)
      endif
      call hdf5_close_file(file_handle)
      NGas = Ngas + NPart_This_sight
    enddo
    call allocate_particledata()
  endif
  !
  Ntotal = 0 ! current number of particles read
  do files=0, NFiles-1 ! loop over files
    if(single_file) then
      call hdf5_open_file(file_handle, longfile, readonly=.true.)
      write (LosNumber,'(I4)') abs(los_number)
      LOSName = "LOS"//trim(adjustl(LOsNumber))
      ! read sightline # los_number
      VarName = trim(LosName)//'/Number_of_part_this_los'
      call hdf5_read_attribute(file_handle, VarName, NGas)
      call hdf5_read_attribute(file_handle,'/Constants/CM_PER_MPC',cm_per_mpc) 
      call hdf5_read_attribute(file_handle,'/Constants/PROTONMASS',proton_mass) 
      call hdf5_read_attribute(file_handle,'/Constants/SOLAR_MASS',solar_mass) 
      call hdf5_read_attribute(file_handle,'/Header/ExpansionFactor', ExpansionFactor)
      call hdf5_read_attribute(file_handle,'/Header/Redshift', Redshift)
      NPart_This_sight = NGas
      call allocate_particledata()
    else
      write(FileNumber,'(i3)') files
      longfile = trim(basefile)//'.'//trim(adjustl(FileNumber))//'.los.hdf5'        
      LOSName = "LOS"//trim(adjustl(LosNumber))
      call hdf5_open_file(file_handle, longfile, readonly=.true.)
      VarName = trim(LosName)//'/Npart_this_sight'
      call hdf5_read_attribute(file_handle, VarName, NPart_This_sight)
    endif
    !
    if(NPart_This_sight .gt. 0) then
#if defined (EAGLE) || defined (OWLS)
      VarName = trim(LosName)//'/Positions'
      call hdf5_read_data(file_handle, VarName, Position(1:3,NTotal+1:NTotal+NPart_This_sight))
      VarName = trim(LosName)//'/Positions/h-scale-exponent'
      call hdf5_read_attribute(file_handle,VarName,Pos_h_exp)
      VarName = trim(LosName)//'/Positions/aexp-scale-exponent'
      call hdf5_read_attribute(file_handle,VarName,Pos_aexp_exp)
      VarName = trim(LosName)//'/Positions/CGSConversionFactor'
      call hdf5_read_attribute(file_handle,VarName,Pos_cgs_unit)
#else
      VarName = trim(LosName)//'/Coordinates'
      call hdf5_read_data(file_handle, VarName, Position(1:3,NTotal+1:NTotal+NPart_This_sight))
      VarName = trim(LosName)//'/Coordinates/h-scale-exponent'
      call hdf5_read_attribute(file_handle,VarName,Pos_h_exp)
      VarName = trim(LosName)//'/Coordinates/aexp-scale-exponent'
      call hdf5_read_attribute(file_handle,VarName,Pos_aexp_exp)
      VarName = trim(LosName)//'/Coordinates/CGSConversionFactor'
      call hdf5_read_attribute(file_handle,VarName,Pos_cgs_unit)
#endif
      ! read projection parameters
#if defined (EAGLE) || defined(OWLS)
      VarName = trim(LosName)//'/x-position'
      call hdf5_read_attribute(file_handle,VarName,x_comoving)
      VarName = trim(LosName)//'/y-position'
      call hdf5_read_attribute(file_handle,VarName,y_comoving)
      z_comoving = 0
      theta_projection = 0
      phi_projection = 0
#else
      VarName = trim(LosName)//'/ProjectionParameters/X-position'
      call hdf5_read_attribute(file_handle,VarName,x_comoving)
      VarName = trim(LosName)//'/ProjectionParameters/Y-position'
      call hdf5_read_attribute(file_handle,VarName,y_comoving)
      VarName = trim(LosName)//'/ProjectionParameters/Z-position'
      call hdf5_read_attribute(file_handle,VarName,z_comoving)
      VarName = trim(LosName)//'/ProjectionParameters/theta'
      call hdf5_read_attribute(file_handle,VarName,theta_projection)
      VarName = trim(LosName)//'/ProjectionParameters/phi'
      call hdf5_read_attribute(file_handle,VarName,phi_projection)
#endif
      VarName = trim(LosName)//'/x-axis'
      call hdf5_read_attribute(file_handle,VarName,x_axis)
      VarName = trim(LosName)//'/y-axis'
      call hdf5_read_attribute(file_handle,VarName,y_axis)
      VarName = trim(LosName)//'/z-axis'
      call hdf5_read_attribute(file_handle,VarName,z_axis)
      !
      VarName = trim(LosName)//'/Velocity'
      call hdf5_read_data(file_handle, VarName, Velocity(1:3,NTotal+1:NTotal+NPart_This_sight))
      VarName = trim(LosName)//'/Velocity/h-scale-exponent'
      call hdf5_read_attribute(file_handle,VarName,Vel_h_exp)
      VarName = trim(LosName)//'/Velocity/aexp-scale-exponent'
      call hdf5_read_attribute(file_handle,VarName,Vel_aexp_exp)
      VarName = trim(LosName)//'/Velocity/CGSConversionFactor'
      call hdf5_read_attribute(file_handle,VarName,Vel_cgs_unit)
      !
      ! Density
      VarName = trim(LosName)//'/Density'        
      call hdf5_read_data(file_handle, VarName, ParticleDensity(NTotal+1:NTotal+NPart_This_sight))
      VarName = trim(LosName)//'/Density/h-scale-exponent'
      call hdf5_read_attribute(file_handle,VarName,Dens_h_exp)
      VarName = trim(LosName)//'/Density/aexp-scale-exponent'
      call hdf5_read_attribute(file_handle,VarName,Dens_aexp_exp)
      VarName = trim(LosName)//'/Density/CGSConversionFactor'
      call hdf5_read_attribute(file_handle,VarName,Dens_cgs_unit)
      !
      ! Star formation rate
      VarName = trim(LosName)//'/StarFormationRate'        
      call hdf5_read_data(file_handle,VarName,StarFormationRate(NTotal+1:NTotal+NPart_This_sight))
      !
      ! Same scalings as Coordinates
      VarName = trim(LosName)//'/SmoothingLength'  
      call hdf5_read_data(file_handle, VarName, ParticleSmoothingLength(NTotal+1:NTotal+NPart_This_sight))
      !
      ! Boundary type
      if(gimic)then
        VarName =  trim(LosName)//'//BoundaryType'
        call hdf5_read_data(file_handle,VarName,Boundary(NTotal+1:NTotal+NPart_This_sight))
      endif
      !
      ! Temperature
      if(use_urchin_temperature) then
        VarName = trim(LosName)//'/Temperature_urchin'
      else
        VarName = trim(LosName)//'/Temperature'
      endif
      call hdf5_read_data(file_handle, VarName, ParticleTemperature(NTotal+1:NTotal+NPart_This_sight))
      VarName = trim(LosName)//'/Temperature/h-scale-exponent'
      call hdf5_read_attribute(file_handle,VarName,Temp_h_exp)
      VarName = trim(LosName)//'/Temperature/aexp-scale-exponent'
      call hdf5_read_attribute(file_handle,VarName,Temp_aexp_exp)
      VarName = trim(LosName)//'/Temperature/CGSConversionFactor'
      call hdf5_read_attribute(file_handle,VarName,Temp_cgs_unit)
      !
      if(urchin) then
        VarName = trim(LosName)//'//HydrogenOneFraction'
        call hdf5_read_data(file_handle, VarName, ParticleNeutralHFraction(NTotal+1:NTotal+NPart_This_sight))
        VarName = trim(LosName)//'//MolecularHydrogenMassFraction'
        call hdf5_read_data(file_handle, VarName, ParticleMolecularHFraction(NTotal+1:NTotal+NPart_This_sight))
      endif
      !
      ! Element abundances
      if(use_smoothed_abundance) then
        VarName = trim(LosName)//'/SmoothedElementAbundance/'
      else
        VarName = trim(LosName)//'/ElementAbundance/'
      endif
      if (requireH) then
        call hdf5_read_data(file_handle,trim(Varname)//'Hydrogen',&
          MassFractions(H_index,NTotal+1:NTotal+NPart_This_sight))   
      endif
      if (requireHe) then
        call hdf5_read_data(file_handle,trim(Varname)//'Helium',&
          MassFractions(He_index,NTotal+1:NTotal+NPart_This_sight))
      endif
      if (requireC) then
        call hdf5_read_data(file_handle,trim(Varname)//'Carbon',&
          MassFractions(C_index,NTotal+1:NTotal+NPart_This_sight))
      endif
      if (requireN) then
        call hdf5_read_data(file_handle,trim(Varname)//'Nitrogen',&
          MassFractions(N_index,NTotal+1:NTotal+NPart_This_sight))
      endif
      if (requireO) then
        call hdf5_read_data(file_handle,trim(Varname)//'Oxygen',&
          MassFractions(O_index,NTotal+1:NTotal+NPart_This_sight))
      endif
      if (requireMg) then
        call hdf5_read_data(file_handle,trim(Varname)//'Magnesium',&
          MassFractions(Mg_index,NTotal+1:NTotal+NPart_This_sight))
      endif
      if (requireNe) then
        call hdf5_read_data(file_handle,trim(Varname)//'Neon',&
          MassFractions(Ne_index,NTotal+1:NTotal+NPart_This_sight))
      endif
      if (requireSi) then
        call hdf5_read_data(file_handle,trim(Varname)//'Silicon',&
          MassFractions(Si_index,NTotal+1:NTotal+NPart_This_sight))
      endif
      if (requireFe) then
        call hdf5_read_data(file_handle,trim(Varname)//'Iron',&
          MassFractions(Fe_index,NTotal+1:NTotal+NPart_This_sight))
      endif
      !
      ! Metallicity
      call hdf5_read_attribute(file_handle,'/Constants/Z_Solar', Z_solar)  ! assumed solar metallicity 
      call hdf5_read_data(file_handle,trim(LosName)//'/Metallicity',Metallicity(NTotal+1:NTotal+NPart_This_sight))
      !
      ! Mass
      VarName = trim(LosName)//'/Mass'
      call hdf5_read_data(file_handle, VarName, Mass(NTotal+1:NTotal+NPart_This_sight))
      VarName = trim(LosName)//'/Mass/h-scale-exponent'
      call hdf5_read_attribute(file_handle,VarName,Mass_h_exp)
      VarName = trim(LosName)//'/Mass/aexp-scale-exponent'
      call hdf5_read_attribute(file_handle,VarName,Mass_aexp_exp)
      VarName = trim(LosName)//'/Mass/CGSConversionFactor'
      call hdf5_read_attribute(file_handle,VarName,Mass_cgs_unit)
      !
      Ntotal = Ntotal + NPart_This_sight
    endif
    !
    call hdf5_close_file(file_handle)
    !
  enddo ! loop over files
  !
  x_fraction = x_comoving / BoxSize
  y_fraction = y_comoving / BoxSize
  z_fraction = z_comoving / BoxSize
  !
  x_physical    = x_fraction * BoxSize / HubbleParam * Expansionfactor ! Physical Mpc
  y_physical    = y_fraction * BoxSize / HubbleParam * Expansionfactor ! Physical Mpc
  z_physical    = z_fraction * BoxSize / HubbleParam * Expansionfactor ! Physical Mpc

  ! Eagle allows for sightlines parallel to Z (x-axis=0, y-axis=1, z_axis=2), or parallel to X (x-axis=1, y-axis=2, z_axis=0), or parallel to Y (x-axis=2, y-axis=0, z-axis=1
  ! we swap x,y and z appropriately
  if(x_axis+y_axis+z_axis .ne. 3) then
     write (*,*) ' this is not allowed'
     call abortrun(' problem with projection axes')
  endif
  if(x_axis .ne. 0 .and. y_axis .ne. 1 .and. z_axis .ne. 2) then
     allocate(tmp(3,NTotal))
     tmp = Position
     Position(1,:) = tmp(x_axis+1,:)
     Position(2,:) = tmp(y_axis+1,:)
     Position(3,:) = tmp(z_axis+1,:)
     tmp = Velocity
     Velocity(1,:) = tmp(x_axis+1,:)
     Velocity(2,:) = tmp(y_axis+1,:)
     Velocity(3,:) = tmp(z_axis+1,:)
     deallocate(tmp)
!      if(MyPE == 0) &
!           write (*,100) x_axis, y_axis, z_axis
! 100  format(' swapped axis to x_axis= ',I1,' y-axis= ',I1,' z-axis= ',I1)
  endif


  !
  ! Convert units, specwizard assumes the following units:
  ! Position:          Physical Mpc (and ParticleSmoothingLength too)
  ! Velocity:          Physical km/s
  ! ParticleTemperature:       K
  ! Density:           Hydrogen density in particles cm^-3
  !
  Coordinates_conv     = Pos_aexp_exp * log10(ExpansionFactor) + pos_h_exp * log10(HubbleParam) + log10(Pos_cgs_unit) & ! in physical cm
    - log10(cm_per_mpc) ! in physical Mpc
  !
  Coordinates_conv     = 10.d0**Coordinates_conv ! convert to proper Mpc
  Position             = Position * Coordinates_conv
  BoxPhys              = BoxSize * Coordinates_conv
  ParticleSmoothingLength  = ParticleSmoothingLength * Coordinates_conv ! CAUTION specwizard neighbours within h, OLD version: 2*h
  !
  Velocity_conv       = Vel_aexp_exp * log10(ExpansionFactor) + vel_h_exp * log10(HubbleParam) + log10(vel_cgs_unit) & 
    - log10(1.d5) ! physical km/s
  Velocity_conv       = 10.d0**Velocity_conv
  Velocity            = Velocity * Velocity_conv
  !
  Density_conv       = Dens_aexp_exp * log10(ExpansionFactor) + dens_h_exp * log10(HubbleParam) + log10(dens_cgs_unit) &
    - log10(proton_mass)
  Density_conv       = 10.d0**Density_conv
  ParticleDensity            = ParticleDensity * Density_conv
  !
  ! Convert from total density to *Hydrogen* number density
  ParticleDensity(:)         = ParticleDensity(:)*MassFractions(H_index,:)
  !
  ! metallicity in solar units
  Zmetal = Metallicity / Z_Solar  ! metallicity in solar units
  !
  Mass_conv          = Mass_aexp_exp * log10(ExpansionFactor) + mass_h_exp * log10(HubbleParam) + log10(mass_cgs_unit) &
    - log10(solar_mass)
  Mass_conv          = 10.d0**Mass_conv
  Mass               = Mass * Mass_conv
  !
  if(HubbleParam .eq. 0) call abortrun('HubbleParam = 0!')
  Temp_conv          = Temp_aexp_exp * log10(ExpansionFactor) + Temp_h_exp * log10(HubbleParam) + log10(Temp_cgs_unit) ! cgs
  Temp_Conv          = 10.d0**Temp_Conv
  ParticleTemperature        = ParticleTemperature * Temp_Conv 
  if(setmaxt4sfgas) then
    where(StarFormationRate .gt. 0) ParticleTemperature = 1.d4
  endif
  !
  nsf = COUNT(StarFormationRate .gt. 0.0)
  !
  if(nsf .GT. 0 .AND. verbose) write(*,'("WARNING: Star forming particle(s) present!")')
  !
  if(ignore_starforming) then
    where(StarFormationRate .gt. 0) Mass = 0.d0 ! Particles on the EOS do not contribute
  endif
  !
  !
  if (verbose .and. first_call .and. MyPE == 0) then 
    write(*,*) ' ++++++++++++ '
    write(*,*) 'Particle data read::'
    write(*,'("redshift      = ",f7.4)') Redshift
    write(*,'("omega_m       = ",f7.4)') omega0
    write(*,'("omega_lambda  = ",f7.4)') Omegalambda
    write(*,'("omega_b*h^2   = ",f7.4)') OmegaBaryon * HubbleParam**2
    write(*,'("h             = ",f7.4)') HubbleParam
    write(*,'("box size      = ",f8.3," h^{-1} Mpc, comoving")') &
      BoxSize
    write(*,*) ' ++++++++++++ '
    write(*,*)
  endif
  !
  ! compute extrema
  tmin        = minval(ParticleTemperature) + small_temp
  tmax        = maxval(ParticleTemperature) + small_temp
  dmin        = minval(ParticleDensity)     + small_rho
  dmax        = maxval(ParticleDensity)     + small_rho
  !
  metal_max   = maxval(Metallicity) + small_metallicity
  !
  MetallicityinSolar = Metallicity / zmass_solar ! convert metallicity to solar values
  !
  !		-------------------
  !		Print extrema info.
  !		-------------------
  !
  !		 write(*,*) ib_t(1), ib_t(ib_nt), ib_d(1), ib_d(ib_nd)
  !
  if (verbose .and. MyPE == 0) then 
    write(*,'("Redshift sim:    ",f11.4)')  Redshift
    write(*,'("# particles:     ",i11)')    Ngas
    write(*,'("# SF particles:	 ",i11)')   nsf
    write(*,'("min log(T):	 ",f11.3)') log10(tmin)
    write(*,'("max log(T):	 ",f11.3)') log10(tmax)
    write(*,'("min log(n_H):	 ",f11.3)') log10(dmin)
    write(*,'("max log(n_H):	 ",f11.3)') log10(dmax)
    write(*,'("max metallicity: ",f11.3)')  log10(metal_max)
  endif
  !
  dmin = dmin/ibfactor
  dmax = dmax/ibfactor
!  if (log10(tmin) .lt. ib_logt(1)) then
!    write(*,'("Min. temp. in ioniz. table too large!: ",f7.3, &
!      " (needed: ",f7.3,")")') ib_logt(1), log10(tmin)
!    write(*,'("Hit a key and then enter to coninue: ")')
!    !				 read(*,*) dumbo
!    !				 stop
!  endif
  if (log10(tmax) .gt. ib_logt(nt)) then
    write(*,'("Max. temp. in ioniz. table too small!: ",f7.3, &
      " (needed: ",f7.3,")")') ib_logt(nt), log10(tmax)
    write(*,'("Hit a key and then enter to coninue: ")')
    !				 read(*,*) dumbo
    !				 stop
  endif
  if (log10(dmin) .lt. ib_logd(1)) then
    write(*,'("Min. dens. in ioniz. table too large!: ",f7.3, &
      " (needed: ",f7.3,")")') ib_logd(1), log10(dmin)
    write(*,'("Hit a key and then enter to coninue: ")')
    !				 read(*,*) dumbo
    !				 stop
  endif
  if (log10(dmax) .gt. ib_logd(nd)) then
    write(*,'("Max. dens. in ioniz. table too small!: ",f7.3, &
      " (needed: ",f7.3,")")') ib_logd(nd), log10(dmax)
    write(*,'("Hit a key and then enter to coninue: ")')
    !				 read(*,*) dumbo
    !				 stop 
  endif
end subroutine readdata_owls


subroutine read_full_snapshot()
  ! Read in data of OWLS simulation (hdf5).
  use numbers
  use header
  use hdf5_wrapper
  use projection_parameters
  use particledata
  use runtime
  use physical_constants
  use atomic_data
  use constants
  use solar_data, only : zmass_solar
  use spectra
  use ionization_tables, only : nz, nt, nd, ib_logt, ib_logd
  use my_mpi
  use modified_metallicity
#ifdef READREGION
  use read_eagle
  use RegionExtent
#endif
  implicit none
  !
  ! local variables
  character(len=120)       :: basefile, urchin_basefile, longfile, urchin_longfile
  integer                  :: file_handle, urchin_file_handle
  integer                  :: part_file_handle
  character(len=120)       :: LosNumber, LosName, VarName
  character(len=3)         :: FileNumber
  integer                  :: los_number
  integer                  :: nfiles
  real(kind=singleR)       :: Pos_h_exp, Pos_aexp_exp, Vel_h_exp, Vel_aexp_exp, &
       Dens_h_exp, Dens_aexp_exp, Eint_h_exp, Eint_aexp_exp, &
       Mass_h_exp, Mass_aexp_exp,  Temp_h_exp, Temp_aexp_exp
  real(kind=doubleR)       :: Mass_cgs_unit, Dens_cgs_unit, Eint_cgs_unit, Temp_cgs_unit, &
       Pos_cgs_unit, Vel_cgs_unit
  real(kind=singleR)       :: z_solar, Bsize
  !real(kind=doubleR)       :: proton_mass
  real(kind=doubleR)       :: Coordinates_conv, Velocity_conv, Density_conv, Mass_conv, Eint_conv, Temp_conv
  real(kind=doubleR)       :: tmin, tmax, dmin, dmax, metal_max, rnorm2propermpc
  integer(kind=singleI)    :: i
  !
  logical                  :: file_exists, single_file
  logical, save :: first_call = .true.
  !
  integer                  :: files, NTotal, NPart_This_file, np
  !
#ifdef READREGION
  type (eaglesnapshot) :: snapinfo
#endif

#ifdef EAGLE
  longfile = trim(datadir)//trim(snap_base)//'.hdf5'
#else
  write (FileNumber,'(I3.3)') snap
  longfile = trim(datadir)//'/snapshot_'//trim(FileNumber)//'/snap_'//trim(FileNumber)//'.hdf5'
#endif
  if(urchin) urchin_longfile = trim(urchindir)//'/snap_999.hdf5'
  !
  inquire(file=longfile,exist=single_file)
  if(.not. single_file)then
    ! multi-file format
#ifdef EAGLE
     basefile = trim(datadir)//trim(snap_base)
#else
    basefile = trim(datadir)//'/snapshot_'//trim(FileNumber)//'/snap_'//trim(FileNumber)
#endif
    longfile = trim(basefile)//'.0.hdf5'
    if(urchin) then
      urchin_basefile = trim(urchindir)//'/snap_999.'
    endif
  else
    basefile = longfile
    urchin_basefile = urchin_longfile
  endif
  !
  if(MyPE == 0) then
     write(*,*)'Reading snapshot file(s): ',trim(basefile)
     if(urchin) write(*,*)'Reading urchin file(s): ',trim(urchin_basefile)
  endif
  !
  call hdf5_open_file(file_handle, longfile, readonly=.true.)
  call read_header(file_handle)
  call check_header
  call hdf5_close_file(file_handle)
  !
  if(single_file) then
     NFiles = 1
  else
     NFiles = NumFilesPerSnapshot
  endif
  !


  ! read all variables
#ifdef READREGION
  longfile = trim(basefile)//'.0.hdf5'        

  ! read system of units
  call hdf5_open_file(file_handle, longfile, readonly=.true.)
  call hdf5_read_attribute(file_handle,'/Constants/CM_PER_MPC',cm_per_mpc) 
  call hdf5_read_attribute(file_handle,'/Constants/PROTONMASS',proton_mass) 
  call hdf5_read_attribute(file_handle,'/Constants/SOLAR_MASS',solar_mass) 
  call hdf5_read_attribute(file_handle,'/Header/ExpansionFactor', ExpansionFactor)
  call hdf5_read_attribute(file_handle,'/Header/Redshift', Redshift)
  call read_header(file_handle)

  ! read unit conversion factors
  VarName = trim('PartType0')//'/Coordinates/h-scale-exponent'
  call hdf5_read_attribute(file_handle,VarName,Pos_h_exp)
  VarName = trim('PartType0')//'/Coordinates/aexp-scale-exponent'
  call hdf5_read_attribute(file_handle,VarName,Pos_aexp_exp)
  VarName = trim('PartType0')//'/Coordinates/CGSConversionFactor'
  call hdf5_read_attribute(file_handle,VarName,Pos_cgs_unit)

  VarName = trim('PartType0')//'/Velocity/h-scale-exponent'
  call hdf5_read_attribute(file_handle,VarName,Vel_h_exp)
  VarName = trim('PartType0')//'/Velocity/aexp-scale-exponent'
  call hdf5_read_attribute(file_handle,VarName,Vel_aexp_exp)
  VarName = trim('PartType0')//'/Velocity/CGSConversionFactor'
  call hdf5_read_attribute(file_handle,VarName,Vel_cgs_unit)
  
  VarName = trim('PartType0')//'/Density/h-scale-exponent'
  call hdf5_read_attribute(file_handle,VarName,Dens_h_exp)
  VarName = trim('PartType0')//'/Density/aexp-scale-exponent'
  call hdf5_read_attribute(file_handle,VarName,Dens_aexp_exp)
  VarName = trim('PartType0')//'/Density/CGSConversionFactor'
  call hdf5_read_attribute(file_handle,VarName,Dens_cgs_unit)
  
  if(use_urchin_temperature) then
     VarName = trim('PartType0')//'/Temperature/h-scale-exponent'
     call hdf5_read_attribute(urchin_file_handle,VarName,Temp_h_exp)
     VarName = trim('PartType0')//'/Temperature/aexp-scale-exponent'
     call hdf5_read_attribute(urchin_file_handle,VarName,Temp_aexp_exp)
     VarName = trim('PartType0')//'/Temperature/CGSConversionFactor'
     call hdf5_read_attribute(urchin_file_handle,VarName,Temp_cgs_unit)
  else
     VarName = trim('PartType0')//'/Temperature/h-scale-exponent'
     call hdf5_read_attribute(file_handle,VarName,Temp_h_exp)
     VarName = trim('PartType0')//'/Temperature/aexp-scale-exponent'
     call hdf5_read_attribute(file_handle,VarName,Temp_aexp_exp)
     VarName = trim('PartType0')//'/Temperature/CGSConversionFactor'
     call hdf5_read_attribute(file_handle,VarName,Temp_cgs_unit)
  endif
  
  VarName = trim('PartType0')//'/Mass/h-scale-exponent'
  call hdf5_read_attribute(file_handle,VarName,Mass_h_exp)
  VarName = trim('PartType0')//'/Mass/aexp-scale-exponent'
  call hdf5_read_attribute(file_handle,VarName,Mass_aexp_exp)
  VarName = trim('PartType0')//'/Mass/CGSConversionFactor'
  call hdf5_read_attribute(file_handle,VarName,Mass_cgs_unit)
  
  call hdf5_read_attribute(file_handle,'/Constants/Z_Solar', Z_solar)  ! assumed solar metallicity 

  call hdf5_close_file(file_handle)

  

  ! define selection
  snapinfo = open_snapshot(longfile)
  Bsize = BoxSize
  !!call select_region(snapinfo, 0.0, BSize, 0.0, BSize, 0.0, BSize)
  call select_region(snapinfo, RegionExtentX(1), RegionExtentX(2), RegionExtentY(1), RegionExtentY(2), RegionExtentZ(1), RegionExtentZ(2))
  !!  call select_region(snapinfo, 4.0, 5.0, 2.0, 3.0, 3.0, 4.0)
  Ngas = count_particles(snapinfo, 0)
  if(MyPE == 0) then
     write(*,*)"Boxsize   = ", snapinfo%boxsize
     write(*,*)"Numfiles  = ", snapinfo%numfiles
     write(*,*)"Hashbits  = ", snapinfo%hashbits
     write(*,'(" NumPart   = ",6(I11,1x))') snapinfo%numpart_total
     write(*,'(" Extent= ",6(f8.2,1x))') RegionExtentX, RegionExtentY, RegionExtentZ
     write (*,*)"Ngas to read = ",Ngas
  endif

  ! allocate particle data
  call allocate_particledata()

  ! read variables
  if (MyPE == 0) &
       write (*,*) ' reading position'
  np = read_dataset(snapinfo, 0, "Coordinates", Position)
  if (MyPE == 0) &
       write (*,*) ' reading velocity'
  np = read_dataset(snapinfo, 0, "Velocity", Velocity)
  if (MyPE == 0) &
       write (*,*) ' reading density'
  np = read_dataset(snapinfo, 0, "Density", ParticleDensity)
  if (MyPE == 0) &
       write (*,*) ' reading SFR'
  np = read_dataset(snapinfo, 0, "StarFormationRate", StarFormationRate)
  if (MyPE == 0) &
       write (*,*) ' reading h'
  np = read_dataset(snapinfo, 0, "SmoothingLength", ParticleSmoothingLength)
  if (MyPE == 0) &
       write (*,*) ' reading ID'
  np = read_dataset(snapinfo, 0, "ParticleIDs", PartID)

  if(urchin) then
     stop ' not implemented urchin'
  else
     ParticleNeutralHFraction = 0
     ParticleMolecularHFraction = 0
  endif
  if(use_urchin_temperature) then
     stop ' not implemented urchin temperature'
  else
     np = read_dataset(snapinfo, 0, "Temperature", ParticleTemperature)
  endif
  if (MyPE == 0) &
       write (*,*) ' reading abundances'

  ! Element abundances
  if(use_smoothed_abundance) then
     VarName = 'SmoothedElementAbundance/'
  else
     VarName = 'ElementAbundance/'
  endif

  if (requireH) &
       np = read_dataset(snapinfo, 0, trim(VarName)//'Hydrogen', MassFractions(H_index,:))
  if (requireHe) &
       np = read_dataset(snapinfo, 0, trim(VarName)//'Helium', MassFractions(He_index,:))
  if (requireC) &
       np = read_dataset(snapinfo, 0, trim(VarName)//'Carbon', MassFractions(C_index,:))
  if (requireSi) &
       np = read_dataset(snapinfo, 0, trim(VarName)//'Silicon', MassFractions(Si_index,:))
  if (requireFe) &
       np = read_dataset(snapinfo, 0, trim(VarName)//'Iron', MassFractions(Fe_index,:))
  if (requireMg) &
       np = read_dataset(snapinfo, 0, trim(VarName)//'Magnesium', MassFractions(Mg_index,:))
  if (requireN) &
       np = read_dataset(snapinfo, 0, trim(VarName)//'Nitrogen', MassFractions(N_index,:))  
  if (requireNe) &
       np = read_dataset(snapinfo, 0, trim(VarName)//'Neon', MassFractions(Ne_index,:))
  if (requireO) &
       np = read_dataset(snapinfo, 0, trim(VarName)//'Oxygen', MassFractions(O_index,:))

  ! Metallicity
  if (MyPE == 0) &
       write (*,*) ' reading Z'
  if(use_smoothed_abundance) then
     np =  read_dataset(snapinfo, 0, 'SmoothedMetallicity', Metallicity)
  else
     np =  read_dataset(snapinfo, 0, 'Metallicity', Metallicity)
  endif

  if (MyPE == 0) &
       write (*,*) ' reading Mass'
  np =  read_dataset(snapinfo, 0, 'Mass', Mass)
  if (MyPE == 0) &
       write (*,*) ' reading done'
  

  call clear_selection(snapinfo)
  call close_snapshot(snapinfo)

#else
  ! determine total number of Gas particles
  NGas   = NumPart_Total(0)
  call allocate_particledata()
  !
  Ntotal = 0 ! current number of particles read
  do files=0, NFiles-1
     if(single_file) then
        call hdf5_open_file(file_handle, longfile, readonly=.true.)
        if(urchin) call hdf5_open_file(urchin_file_handle, urchin_longfile, readonly=.true.)
     else
        write(FileNumber,'(i3)') files
        longfile = trim(basefile)//'.'//trim(adjustl(FileNumber))//'.hdf5'        
        call hdf5_open_file(file_handle, longfile, readonly=.true.)
        if(urchin) then
          urchin_longfile = trim(urchin_basefile)//trim(adjustl(FileNumber))//'.hdf5'
          call hdf5_open_file(urchin_file_handle, urchin_longfile, readonly=.true.)
        endif
     endif
     !
     call hdf5_read_attribute(file_handle,'/Constants/CM_PER_MPC',cm_per_mpc) 
     call hdf5_read_attribute(file_handle,'/Constants/PROTONMASS',proton_mass) 
     call hdf5_read_attribute(file_handle,'/Constants/SOLAR_MASS',solar_mass) 
     call hdf5_read_attribute(file_handle,'/Header/ExpansionFactor', ExpansionFactor)
     call hdf5_read_attribute(file_handle,'/Header/Redshift', Redshift)
     !
     call read_header(file_handle)
     Npart_this_file = NumPart_ThisFile(0)
     !
     VarName = trim('PartType0')//'/Coordinates'
     call hdf5_read_data(file_handle, VarName, Position(1:3,NTotal+1:NTotal+Npart_this_file))
     VarName = trim('PartType0')//'/Coordinates/h-scale-exponent'
     call hdf5_read_attribute(file_handle,VarName,Pos_h_exp)
     VarName = trim('PartType0')//'/Coordinates/aexp-scale-exponent'
     call hdf5_read_attribute(file_handle,VarName,Pos_aexp_exp)
     VarName = trim('PartType0')//'/Coordinates/CGSConversionFactor'
     call hdf5_read_attribute(file_handle,VarName,Pos_cgs_unit)
     !
     VarName = trim('PartType0')//'/Velocity'
     call hdf5_read_data(file_handle, VarName, Velocity(1:3,NTotal+1:NTotal+Npart_this_file))
     VarName = trim('PartType0')//'/Velocity/h-scale-exponent'
     call hdf5_read_attribute(file_handle,VarName,Vel_h_exp)
     VarName = trim('PartType0')//'/Velocity/aexp-scale-exponent'
     call hdf5_read_attribute(file_handle,VarName,Vel_aexp_exp)
     VarName = trim('PartType0')//'/Velocity/CGSConversionFactor'
     call hdf5_read_attribute(file_handle,VarName,Vel_cgs_unit)
     !
     VarName = trim('PartType0')//'/Density'        
     call hdf5_read_data(file_handle, VarName, ParticleDensity(NTotal+1:NTotal+Npart_this_file))
     VarName = trim('PartType0')//'/Density/h-scale-exponent'
     call hdf5_read_attribute(file_handle,VarName,Dens_h_exp)
     VarName = trim('PartType0')//'/Density/aexp-scale-exponent'
     call hdf5_read_attribute(file_handle,VarName,Dens_aexp_exp)
     VarName = trim('PartType0')//'/Density/CGSConversionFactor'
     call hdf5_read_attribute(file_handle,VarName,Dens_cgs_unit)
     !
     VarName = trim('PartType0')//'/StarFormationRate'        
     call hdf5_read_data(file_handle,VarName,StarFormationRate(NTotal+1:NTotal+Npart_this_file))
     !
     VarName = trim('PartType0')//'/SmoothingLength'  
     call hdf5_read_data(file_handle, VarName, ParticleSmoothingLength(NTotal+1:NTotal+Npart_this_file))
     !
     if(urchin)then
       VarName = trim('PartType0')//'//HydrogenOneFraction'
       call hdf5_read_data(urchin_file_handle, VarName, &
         ParticleNeutralHFraction(NTotal+1:NTotal+Npart_this_file))
       VarName = trim('PartType0')//'//MolecularHydrogenMassFraction'
       call hdf5_read_data(urchin_file_handle, VarName, &
         ParticleMolecularHFraction(NTotal+1:NTotal+Npart_this_file))
     endif
     !
     if(gimic)then
        VarName =  trim('PartType0')//'//BoundaryType'
        call hdf5_read_data(file_handle,VarName,Boundary(NTotal+1:NTotal+Npart_this_file))
     endif
     !
     if(use_urchin_temperature) then
       VarName = trim('PartType0')//'/Temperature'
       call hdf5_read_data(urchin_file_handle, VarName, ParticleTemperature(NTotal+1:NTotal+Npart_this_file))
       VarName = trim('PartType0')//'/Temperature/h-scale-exponent'
       call hdf5_read_attribute(urchin_file_handle,VarName,Temp_h_exp)
       VarName = trim('PartType0')//'/Temperature/aexp-scale-exponent'
       call hdf5_read_attribute(urchin_file_handle,VarName,Temp_aexp_exp)
       VarName = trim('PartType0')//'/Temperature/CGSConversionFactor'
       call hdf5_read_attribute(urchin_file_handle,VarName,Temp_cgs_unit)
     else
       VarName = trim('PartType0')//'/Temperature'
       call hdf5_read_data(file_handle, VarName, ParticleTemperature(NTotal+1:NTotal+Npart_this_file))
       VarName = trim('PartType0')//'/Temperature/h-scale-exponent'
       call hdf5_read_attribute(file_handle,VarName,Temp_h_exp)
       VarName = trim('PartType0')//'/Temperature/aexp-scale-exponent'
       call hdf5_read_attribute(file_handle,VarName,Temp_aexp_exp)
       VarName = trim('PartType0')//'/Temperature/CGSConversionFactor'
       call hdf5_read_attribute(file_handle,VarName,Temp_cgs_unit)
     endif
     !
     ! Element abundances
     if(use_smoothed_abundance) then
       VarName = trim('PartType0')//'/SmoothedElementAbundance/'
     else
       VarName = trim('PartType0')//'/ElementAbundance/'
     endif
     if (requireH) call hdf5_read_data(file_handle,trim(VarName)//'Hydrogen',&
                        MassFractions(H_index,NTotal+1:NTotal+Npart_this_file))
     if (requireHe) call hdf5_read_data(file_handle,trim(VarName)//'Helium',&
                        MassFractions(He_index,NTotal+1:NTotal+Npart_this_file))
     if (requireC) call hdf5_read_data(file_handle,trim(VarName)//'Carbon',&
                        MassFractions(C_index,NTotal+1:NTotal+Npart_this_file))
     if (requireSi) call hdf5_read_data(file_handle,trim(VarName)//'Silicon',&
                        MassFractions(Si_index,NTotal+1:NTotal+Npart_this_file))
     if (requireFe) call hdf5_read_data(file_handle,trim(VarName)//'Iron',&
                        MassFractions(Fe_index,NTotal+1:NTotal+Npart_this_file))
     if (requireMg) call hdf5_read_data(file_handle,trim(VarName)//'Magnesium',&
                        MassFractions(Mg_index,NTotal+1:NTotal+Npart_this_file))
     if (requireN) call hdf5_read_data(file_handle,trim(VarName)//'Nitrogen',&
                        MassFractions(N_index,NTotal+1:NTotal+Npart_this_file))
     if (requireNe) call hdf5_read_data(file_handle,trim(VarName)//'Neon',&
                        MassFractions(Ne_index,NTotal+1:NTotal+Npart_this_file))
     if (requireO) call hdf5_read_data(file_handle,trim(VarName)//'Oxygen',&
                        MassFractions(O_index,NTotal+1:NTotal+Npart_this_file))
     !
     ! Metallicity
     call hdf5_read_attribute(file_handle,'/Constants/Z_Solar', Z_solar)  ! assumed solar metallicity 
     if(use_smoothed_abundance) then
       call hdf5_read_data(file_handle,trim('PartType0')//'/SmoothedMetallicity',Metallicity(NTotal+1:NTotal+Npart_this_file))
     else
       call hdf5_read_data(file_handle,trim('PartType0')//'/Metallicity',Metallicity(NTotal+1:NTotal+Npart_this_file))
     endif
     !
     ! Mass
     VarName = trim('PartType0')//'/Mass'
     call hdf5_read_data(file_handle, VarName, Mass(NTotal+1:NTotal+Npart_this_file))
     VarName = trim('PartType0')//'/Mass/h-scale-exponent'
     call hdf5_read_attribute(file_handle,VarName,Mass_h_exp)
     VarName = trim('PartType0')//'/Mass/aexp-scale-exponent'
     call hdf5_read_attribute(file_handle,VarName,Mass_aexp_exp)
     VarName = trim('PartType0')//'/Mass/CGSConversionFactor'
     call hdf5_read_attribute(file_handle,VarName,Mass_cgs_unit)
     !
     Ntotal = Ntotal + Npart_this_file    
     call hdf5_close_file(file_handle)
     if(urchin) call hdf5_close_file(urchin_file_handle)
     if (MyPE == 0) &
          write (*,*) ' read file = ',files+1,' of ',NFiles
  enddo
#endif
  !
  ! Convert units, specwizard assumes the following units:
  ! Position:          Physical Mpc (and ParticleSmoothingLength too)
  ! Velocity:          Physical km/s
  ! ParticleTemperature:       K
  ! Density:           Hydrogen density in particles cm^-3
  !
  if(MyPE == 0) &
       write(*,*) ' read files, now converting units etc'
  !
  Coordinates_conv     = Pos_aexp_exp * log10(ExpansionFactor) + pos_h_exp * log10(HubbleParam) + log10(Pos_cgs_unit) & ! in physical cm
       - log10(cm_per_mpc) ! in physical Mpc
  !
  Coordinates_conv     = 10.d0**Coordinates_conv
  Position        = Position * Coordinates_conv ! convert to proper physical Mpc
  BoxPhys              = BoxSize * Coordinates_conv
  !
  ! note difference in definition of h between gadget (neighbours within h)
  ParticleSmoothingLength  = ParticleSmoothingLength * Coordinates_conv
  !
  Velocity_conv = Vel_aexp_exp * log10(ExpansionFactor) + vel_h_exp * log10(HubbleParam) + log10(vel_cgs_unit) & 
       - log10(1.d5) ! physical km/s
  Velocity_conv = 10.d0**Velocity_conv
  Velocity = Velocity * Velocity_conv
  !
  Density_conv = Dens_aexp_exp * log10(ExpansionFactor) + dens_h_exp * log10(HubbleParam) + log10(dens_cgs_unit) &
       - log10(proton_mass)
  !
  Density_conv       = 10.d0**Density_conv
  ParticleDensity = ParticleDensity * Density_conv
  !
  ! Convert from total density to *Hydrogen* number density
  ParticleDensity(:) = ParticleDensity(:)*MassFractions(H_index,:)
  !
  ! metallicity in solar units
  Zmetal = Metallicity / Z_Solar  ! metallicity in solar units
  !
  if (read_part_ids_from_file) then
     !
     !Open partid.hdf5:	
     call hdf5_open_file(part_file_handle, particle_file_name, readonly=.true.)
     !
     !read data from partid.hdf5 file:
     call hdf5_read_data(part_file_handle, 'MyData', partid)
     !
     call hdf5_close_file(part_file_handle)
     !
     call modify_metallicity_of_flagged_particles()
     !
  endif
  !
  Mass_conv          = Mass_aexp_exp * log10(ExpansionFactor) + mass_h_exp * log10(HubbleParam) + log10(mass_cgs_unit) &
       - log10(solar_mass)
  Mass_conv          = 10.d0**Mass_conv
  Mass               = Mass * Mass_conv
  !
  if(HubbleParam .eq. 0) call abortrun('HubbleParam = 0!')
  Temp_conv          = Temp_aexp_exp * log10(ExpansionFactor) + Temp_h_exp * log10(HubbleParam) + log10(Temp_cgs_unit) ! cgs
  Temp_Conv          = 10.d0**Temp_Conv
  ParticleTemperature        = ParticleTemperature * Temp_Conv 
  !
  if(setmaxt4sfgas) then
    where(StarFormationRate .gt. 0) ParticleTemperature = 1.d4
  endif
  !
  if(ignore_starforming) then
    where(StarFormationRate .gt. 0) Mass = 0.d0 ! Particles on the EOS do not contribute
  endif
  !
  MetallicityinSolar = Metallicity / zmass_solar ! convert metallicity to solar values
  if(MyPE == 0) &
       write(*,*) ' Snapshot read'
  !
end subroutine read_full_snapshot


! +++++++++++++++++++++++++++++++++++++ ...read/write ++++++++++++++++++++++++++++++++++++++++++++ !

! +++++++++++++++++++++++++++++++++++++++++ particles... +++++++++++++++++++++++++++++++++++++++++ !


subroutine allocate_particledata
  use particledata
  use runtime, only: gimic
  use modified_metallicity, only: read_part_ids_from_file
  implicit none
  !
  integer(kind=singleI), save :: NGas_Old=-1
  integer :: AllocateStatus
  !
  if(Ngas .ne. nGas_old)then
     if(allocated(Mass)) deallocate(Mass)
     if(allocated(Position)) deallocate(Position)
     !!!if(allocated(ShiftedPosition)) deallocate(ShiftedPosition)
     if(allocated(Velocity)) deallocate(Velocity)
     !!!if(allocated(ShiftedVelocity)) deallocate(ShiftedVelocity)
     if(allocated(ParticleDensity)) deallocate(ParticleDensity)
     if(allocated(ParticleSmoothingLength)) deallocate(ParticleSmoothingLength)
     if(allocated(ParticleTemperature)) deallocate(ParticleTemperature)
     if(allocated(ParticleNeutralHFraction)) deallocate(ParticleNeutralHFraction)
     if(allocated(ParticleMolecularHFraction)) deallocate(ParticleMolecularHFraction)
     if(allocated(Metallicity)) deallocate(Metallicity)
     if(allocated(MetallicityInSolar)) deallocate(MetallicityInSolar)
     if(allocated(MassFractions)) deallocate(MassFractions)
     if(allocated(Zmetal)) deallocate(ZMetal)
     if(allocated(Zrelat)) deallocate(Zrelat)
     if(allocated(StarFormationRate)) deallocate(StarFormationRate)
     if(allocated(Boundary)) deallocate(Boundary)
     if(allocated(partid)) deallocate(partid)
     !
     allocate(Mass(NGas),ParticleDensity(Ngas),ParticleSmoothingLength(NGas),ParticleTemperature(NGas),&
       Metallicity(NGas),MetallicityInSolar(NGas),ZMetal(NGas),ZRelat(NGas),StarFormationRate(Ngas),&
       ParticleNeutralHFraction(Ngas),ParticleMolecularHFraction(Ngas),PartID(Ngas),stat=AllocateStatus)
     if(AllocateStatus /= 0) &
          call abortrun("ran out of memory allocating particle data")


     if (gimic) allocate(Boundary(Ngas))
     !!!allocate(Position(3,NGas),Velocity(3,NGas),ShiftedPosition(3,NGas),ShiftedVelocity(3,NGas))
     allocate(Position(3,NGas),Velocity(3,NGas),stat=AllocateStatus ) !!!
     if(AllocateStatus /= 0) &
          call abortrun("ran out of memory allocating particle data")
     allocate(MassFractions(nspecies,NGas),stat=AllocateStatus)
     if(AllocateStatus /= 0) &
          call abortrun("ran out of memory allocating particle data")
     !
     if (read_part_ids_from_file) then
       allocate(partid(NGas))
     endif
     !
     NGas_old = NGas
  endif
end subroutine allocate_particledata


subroutine impose_equation_of_state()
  use numbers
  use particledata
  use header
  use my_mpi
  use runtime
  use atomic_data, only : massH
  implicit none
  !
  integer :: i
  real(kind=doubleR) :: overdens,tcalc
  !
  if (MyPE == 0) write(*,*)'Imposing EOS'
  ! !
  ! ! Equation of state is of the form T=T0*Delta^(gamma-1)
  ! ! and represents a minimum temperature
  ! ! The equation of state is applied up to a maximum
  ! ! overdensity of maxod.
  ! !
  ! do i=1, nGas
  !   overdens = ParticleDensity(i) * massH / rhocb ! ParticleDensity usually refers to hydrogen number density
  !   overdens = overdens / MassFractions(H_index,i)
  !   if (overdens .lt. imposed_eos_maxod) then
  !     tcalc = imposed_eos_T0 * overdens**(imposed_eos_gamma-1)
  !     if (ParticleTemperature(i) .lt. tcalc) then
  !       ParticleTemperature(i) = tcalc
  !     endif
  !   endif
  ! enddo
  ! !
  do i=1, NGas
     overdens = ParticleDensity(i) * massH / rhocb ! ParticleDensity usually refers to hydrogen number density
     overdens = overdens / MassFractions(H_index,i)
     tcalc    = imposed_eos_T0 * overdens**(imposed_eos_gamma-1)

     ! only impose EOS for gas in IGM
     if(ParticleTemperature(i) .lt. 2d4 .and. overdens .lt. imposed_eos_maxod) then
        ParticleTemperature(i) = tcalc
     endif
  enddo

end subroutine impose_equation_of_state


subroutine impose_metallicity()
  ! n.b. Metallicities in zmetal are relative to solar
  use numbers
  use particledata
  use header
  use solar_data
  use primordial_BBN
  use runtime
  use modified_metallicity
  use random_numbers
  use my_random_numbers
  use atomic_data
  implicit none
  !
  ! Local variables
  integer            :: i,ix,iy,iz,iczfill,ispecies
  real(kind=doubleR) :: overden
  real(kind=singleR) :: z_var_dex(z_sig_bin, z_sig_bin, z_sig_bin)
  real(kind=singleR) :: z_rvar_dex(z_sig_bin, z_sig_bin, z_sig_bin)
  real(kind=doubleR) :: z_scale, X_hydrogen
  integer, save      :: first_call=1
  logical            :: OK
  !
  if(.not. modify_metallicity) return
  !
  ! Check for sane parameter values:
  if(first_call == 1) then
    first_call = 0
    if (verbose) then 
      if(log_normal_scatter) then 
        write(*,*)"Adding ",Z_sig_dex," dex lognormal metallicity scatter"
      endif
      if(scale_simulation_abundances) then
        write(*,*)'Scaling simulation abundances by factor ',z_rel
      endif
      !
      if(impose_z_rho_relation) then
        write(*,*)'Imposing power law rho-metallicity relation with slope ',z_index,' and coefficient ',z_mean
      endif
    endif
    !
    OK = .true.
    if(z_rel .lt. 0) then
      write (*,*) 'z_rel strange value'
      OK = .false.
    endif
    !
    if(log_normal_scatter .and. z_sig_bin .le. 0) then
      write(*,*)'Cannot have negative z_sig_bin for scatter!'
      OK = .false.
    endif
    !
    if(scale_simulation_abundances .and. impose_z_rho_relation) then
      write(*,*)'impose_z_rho_relation overrides scale_simulation_abundances, can not use both'
      OK = .false.
    endif
    !
    if(.not. OK) call abortrun('unhealthy parameters in modify_metallicity')
    !
  endif
  !
  ! Scale zmetal (=Z/Z_solar) in the next bits of code. Hydrogen and Helium abundances will
  ! be scaled down to compensate in the final loop of this function:      
  !
  ! Scale simulation abundances of all metals
  if(scale_simulation_abundances) then
    zmetal(:) = z_rel * zmetal(:)
  endif
  !
  ! Divide computational volume in (z_sig_bin)^3 cells, and generate
  ! lognormally distributed random variable for assiging total metallicity
  ! (z_var_dex)
  if(log_normal_scatter)then
    do ix = 1, Z_sig_bin
      do iy = 1, Z_sig_bin
        do iz = 1, Z_sig_bin
          z_var_dex(ix,iy,iz)  = random(ran_metal)*Z_sig_dex
        enddo
      enddo
    enddo
  endif
  !
  if (impose_z_rho_relation) then
    do i=1, nGas
      ! impose relation between overdensity and metallicity
      ! z = z_rel*(rho/<rho>)^z_index        
      ! use Primordial Hydrogen mass fraction
      Overden  =  (ParticleDensity(i)/MassFractions(H_index,i))*(massH/(1.-Ymass)) / rhocb
      zmetal(i) = min(z_mean * OverDen**z_index, maxz_rel)
    enddo
  endif
  !
  if(log_normal_scatter) then
    do i=1, nGas
      ! impose log-normal scatter
      ! find sub-region particle is in
      ix = int(Position(1,i)/BoxPhys * z_sig_bin)+1
      iy = int(Position(2,i)/BoxPhys * z_sig_bin)+1
      iz = int(Position(3,i)/BoxPhys * z_sig_bin)+1
      ! impose scatter appropriate for that region
      zmetal(i) = min(10.0**(log10(zmetal(i)) + z_var_dex(ix,iy,iz)), maxz_rel)
    enddo
  endif
  !
  ! Now we update the H and He abundances, and also scale individual element abundances
  ! if required
  !
  ! ParticleDensity usually denotes *Hydrogen* density
  ! First convert from hydrogen to particle density using original Hydrogen fraction
  ! and convert back after this loop, using the new Hydrogen fraction
  ParticleDensity(:) = ParticleDensity(:) / MassFractions(H_index,:)
  !
  ! Let X=M_H/M_tot, Y=M_He/M_tot, Z=M_Z/M_tot
  !                  y=N_He/N_H
  !
  ! ***To calculate X:
  ! M_H = X M_tot   ; M_He = y * massHe/massH * M_H =  (y*massHe/massH) * X * M_tot
  ! since M_H + M_He = M_tot - M_Z = M_tot (1-Z)
  ! therefore X M_tot (1+(y*massHe/massH)) = (1-Z) M_tot
  ! => X = (1-Z)/(1+y*massHe/massH)
  !
  ! ***To calculate the mass fraction of element Zi at solar metallicity::
  ! (M_Zi/M_tot) = (M_Z/M_tot) * (N_Zi/N_H)_solar * (massZi/massH) * (M_H/M_Z)
  !              = Z * (M_Zi/M_H)_solar * (X/Z)
  !              = (M_Zi/M_H)_solar * X
  ! then scale mass relative to this with zmetal (=Z/Z_solar).
  do i=1, NGas
    ! Mass fraction of Hydrogen
    X_Hydrogen = (1.-zmetal(i)*Zmass_solar) / (1.+YNumber_solar*MassHe/MassH)
    !
    do ispecies=1, nspecies
      !
      if(ispecies .eq. H_index) then
        MassFractions(ispecies,i)  = X_Hydrogen
      else if(ispecies .eq. He_index) then
        MassFractions(ispecies,i)  = YNumber_solar * X_Hydrogen * massHe / massH
      else if(ispecies .eq. C_index) then
        MassFractions(ispecies,i)  = zmetal(i) * ZC_rel * ZC_solar   * X_hydrogen * massC/massH
      else if(ispecies .eq. N_index) then
        MassFractions(ispecies,i)  = zmetal(i) * ZN_rel * ZN_solar   * X_hydrogen * massN/massH
      else if(ispecies .eq. O_index) then
        MassFractions(ispecies,i)  = zmetal(i) * ZO_rel * ZO_solar   * X_hydrogen * massO/massH
      else if(ispecies .eq. Ne_index) then
        MassFractions(ispecies,i)  = zmetal(i) * ZNe_rel * ZNe_solar * X_hydrogen * massNe/massH
      else if(ispecies .eq. Mg_index) then
        MassFractions(ispecies,i)  = zmetal(i) * ZMg_rel * ZMg_solar * X_hydrogen * massMg/massH
      else if(ispecies .eq. Al_index) then
        MassFractions(ispecies,i)  = zmetal(i) * ZAl_rel * ZAl_solar * X_hydrogen * massAl/massH
      else if(ispecies .eq. Si_index) then
        MassFractions(ispecies,i)  = zmetal(i) * ZSi_rel * ZSi_solar * X_hydrogen * massSi/massH
      else if(ispecies .eq. S_index) then
        MassFractions(ispecies,i)  = zmetal(i) * ZS_rel * ZS_solar * X_hydrogen * massS/massH
      else if(ispecies .eq. Fe_index) then
        MassFractions(ispecies,i)  = zmetal(i) * ZFe_rel * ZFe_solar * X_hydrogen * massFe/massH
      else
        write (*,*) 'error: species: ',ispecies,' does not exist'
        call abortrun('error in modify_metallicity')
      endif
      !
    enddo
    !
  enddo
  !
  ! Rescale using new value
  ParticleDensity(:) = ParticleDensity(:) * MassFractions(H_index,:)
  !
  ! Error check
  if (minval(ParticleDensity) .lt. 0.0) then
    write(*,*)'Hydrogen density < 0 in impose metallicity, this usually happens when'
    write(*,*)'Metallicity is too high, lower z_rel or decrease metallicity scatter'
    call abortrun('stop')
  endif
  !
end subroutine impose_metallicity


subroutine modify_metallicity_of_flagged_particles()
  use numbers
  use particledata
  use solar_data
  use modified_metallicity
  use atomic_data, only : massH
  !
  implicit none 
  integer :: i,ispecies
  double precision :: X_Hydrogen
  !multiplies passed array by flagged_particle_mettalicity
  !
  ! ParticleDensity usually denotes *Hydrogen* density
  ! First convert from hydrogen to particle density using original Hydrogen fraction
  ! and convert back after this loop, using the new Hydrogen fraction
  ParticleDensity(:) = ParticleDensity(:) / MassFractions(H_index,:)
  !
  ! Let X=M_H/M_tot, Y=M_He/M_tot, Z=M_Z/M_tot
  !                  y=N_He/N_H
  !
  ! ***To calculate X:
  ! M_H = X M_tot   ; M_He = y * massHe/massH * M_H =  (y*massHe/massH) * X * M_tot
  ! since M_H + M_He = M_tot - M_Z = M_tot (1-Z)
  ! therefore X M_tot (1+(y*massHe/massH)) = (1-Z) M_tot
  ! => X = (1-Z)/(1+y*massHe/massH)
  !
  ! ***To calculate the mass fraction of element Zi at solar metallicity::
  ! (M_Zi/M_tot) = (M_Z/M_tot) * (N_Zi/N_H)_solar * (massZi/massH) * (M_H/M_Z)
  !              = Z * (M_Zi/M_H)_solar * (X/Z)
  !              = (M_Zi/M_H)_solar * X
  ! then scale mass relative to this with zmetal (=Z/Z_solar).
  do i=1, NGas
    do ispecies=1,nspecies
      if (ispecies == C_index) then
        MassFractions(ispecies,i)  = Zmetal(i) * ZC_rel * ZC_solar   * X_hydrogen * massC/massH
        cycle
      endif
      if (ispecies == N_index) then
        MassFractions(ispecies,i)  = Zmetal(i) * ZN_rel * ZN_solar   * X_hydrogen * massN/massH
        cycle
      endif
      if (ispecies == O_index) then
        MassFractions(ispecies,i)  = Zmetal(i) * ZO_rel * ZO_solar   * X_hydrogen * massO/massH
        cycle
      endif
      if (ispecies == Ne_index) then
        MassFractions(ispecies,i)  = Zmetal(i) * ZNe_rel * ZNe_solar * X_hydrogen * massNe/massH
        cycle
      endif
      if (ispecies == Mg_index) then
        MassFractions(ispecies,i)  = Zmetal(i) * ZMg_rel * ZMg_solar * X_hydrogen * massMg/massH
        cycle
      endif
      if (ispecies == Al_index) then
        MassFractions(ispecies,i)  = Zmetal(i) * ZAl_rel * ZAl_solar * X_hydrogen * massAl/massH
        cycle
      endif
      if (ispecies == Si_index) then
        MassFractions(ispecies,i)  = Zmetal(i) * ZSi_rel * ZSi_solar * X_hydrogen * massSi/massH
        cycle
      endif
      if (ispecies == S_index) then
        MassFractions(ispecies,i)  = Zmetal(i) * ZS_rel * ZS_solar * X_hydrogen * massS/massH
        cycle
      endif
      if (ispecies == Fe_index) then
        MassFractions(ispecies,i)  = Zmetal(i) * ZFe_rel * ZFe_solar * X_hydrogen * massFe/massH
        cycle
      endif
      write (*,*) 'error: species: ',ispecies,' does not exist'
      call abortrun('error in modify_metallicity')
    enddo
    !
  enddo
  !
  ! Rescale using new value
  ParticleDensity(:) = ParticleDensity(:) * MassFractions(H_index,:)
  !
end subroutine modify_metallicity_of_flagged_particles


! ++++++++++++++++++++++++++++++++++++++ ...particles ++++++++++++++++++++++++++++++++++++++++++++ !

! ++++++++++++++++++++++++++++++++++++++++ sightlines... +++++++++++++++++++++++++++++++++++++++++ !


subroutine SetRotationMatrix(theta,phi,RotationMatrix)
  use numbers
  implicit none
  real(kind=doubleR), intent(in)  :: theta, phi
  real(kind=doubleR), intent(out) :: RotationMatrix(3,3)
  ! local variables
  real(kind=doubleR) :: r1(3,3), r2(3,3)

  r1      = 0.d0
  r1(1,1) = 1.d0
  r1(2,2) = cos(theta)
  r1(2,3) = sin(theta)
  r1(3,2) = -sin(theta)
  r1(3,3) = cos(theta)
  !
  r2      = 0.d0
  r2(1,1) = cos(phi)
  r2(1,2) = sin(phi)
  r2(2,1) = -sin(phi)
  r2(2,2) = cos(phi)
  r2(3,3) = 1.d0
  !
  RotationMatrix   = matmul(r1,r2)

end subroutine SetRotationMatrix




! +++++++++++++++++++++++++++++++++++++ ...sightlines ++++++++++++++++++++++++++++++++++++++++++++ !

! ++++++++++++++++++++++++++++++++++++++++ ionization... +++++++++++++++++++++++++++++++++++++++++ !


subroutine initialize_ionization_tables
  use physical_constants
  use spectra
  use atomic_data
  use modified_metallicity
  use runtime
  use ionization_tables
  use particledata, only : H_index, He_index, C_index, N_index, &
    O_index, Ne_index, Mg_index, Al_index, Si_index, S_index, Fe_index, &
    requireH, requireHe, requireC, requireN, requireO, requireMg, requireNe, &
    requireAl, requireSi, requireS, requireFe, nspecies, ElementAtomicMass
  use spectra
  use my_mpi
  implicit none
  !
  integer(kind=singleI) ::  i, j, ib, it, id, iz, ii, ion, ioff
  real(kind=doubleR)    :: r1, ran3
  !
  ! choose which ions will be included in generation of spectra 
  ! defaults are false (see modules), except for hydrogen 
  nion = 0
  if (doH1) then
     nion = nion + 1
     requireH = .true.
  endif
  if (doHe2) then 
     nion = nion + 1
     requireHe = .true.
  endif
  if (do21cm) then 
     nion = nion + 1
     requireH = .true.
  endif
  if (doC2) then 
     nion = nion + 1
     requireC = .true.
  endif
  if (doC3) then 
     nion = nion + 1
     requireC = .true.
  endif
  if (doC4) then 
     nion = nion + 1
     requireC = .true.
  endif
  if (doC5) then 
     nion = nion + 1
     requireC = .true.
  endif
  if (doC6) then 
     nion = nion + 1
     requireC = .true.
  endif
  if (doN2) then 
     nion = nion + 1
     requireN = .true.
  endif
  if (doN3) then 
     nion = nion + 1
     requireN = .true.
  endif
  if (doN4) then 
     nion = nion + 1
     requireN = .true.
  endif
  if (doN5) then 
     nion = nion + 1
     requireN = .true.
  endif
  if (doN6) then 
     nion = nion + 1
     requireN = .true.
  endif
  if (doN7) then 
     nion = nion + 1
     requireN = .true.
  endif
  if (doO1) then 
     nion = nion + 1
     requireO = .true.
  endif
  if (doO3) then 
     nion = nion + 1
     requireO = .true.
  endif
  if (doO4) then 
     nion = nion + 1
     requireO = .true.
  endif
  if (doO5) then 
     nion = nion + 1
     requireO = .true.
  endif
  if (doO6) then 
     nion = nion + 1
     requireO = .true.
  endif
  if (doO7) then 
     nion = nion + 1
     requireO = .true.
  endif
  if (doO8) then 
     nion = nion + 1
     requireO = .true.
  endif
  if (doMg2) then 
     nion = nion + 1
     requireMg = .true.
  endif
  if (doNe8) then 
     nion = nion + 1
     requireNe = .true.
  endif
  if (doNe9) then 
     nion = nion + 1
     requireNe = .true.
  endif
  if (doAl2) then 
     nion = nion + 1
     requireAl = .true.
  endif
  if (doAl3) then 
     nion = nion + 1
     requireAl = .true.
  endif
  if (doSi2) then 
     nion = nion + 1
     requireSi = .true.
  endif
  if (doSi3) then 
     nion = nion + 1
     requireSi = .true.
  endif
  if (doSi4) then 
     nion = nion + 1
     requireSi = .true.
  endif
  if (doS5) then 
     nion = nion + 1
     requireS = .true.
  endif
  if (doFe2) then 
     nion = nion + 1
     requireFe = .true.
  endif
  if (doFe3) then 
     nion = nion + 1
     requireFe = .true.
  endif
  if (doFe17) then 
     nion = nion + 1
     requireFe = .true.
  endif
  if (doFe19) then 
     nion = nion + 1
     requireFe = .true.
  endif
  if (doFe21) then 
     nion = nion + 1
     requireFe = .true.
  endif
  !
  nspecies=0
  ioff = 1
  if(MyPE == 0) &
       write(*,*)'Elements required for this run : '
  if (requireH)  then
     if(MyPE == 0) write(*,*)' -- H'
     nspecies=nspecies+1
     H_index = ioff
     ioff = ioff + 1
  endif
  if (requireHe) then
      if(MyPE == 0) write(*,*)' -- He'
     nspecies=nspecies+1
     He_index = ioff
     ioff = ioff + 1
  endif
  if (requireC)  then
      if(MyPE == 0) write(*,*)' -- C'
     nspecies=nspecies+1
     C_index = ioff
     ioff = ioff + 1
  endif
  if (requireN)  then
      if(MyPE == 0) write(*,*)' -- N'
     nspecies=nspecies+1
     N_index = ioff
     ioff = ioff + 1
  endif
  if (requireO)  then
      if(MyPE == 0) write(*,*)' -- O'
     nspecies=nspecies+1
     O_index = ioff
     ioff = ioff + 1
  endif
  if (requireMg) then
      if(MyPE == 0) write(*,*)' -- Mg'
     nspecies=nspecies+1
     Mg_index = ioff
     ioff = ioff + 1
  endif
  if (requireNe) then 
      if(MyPE == 0) write(*,*)' -- Ne'
     nspecies=nspecies+1
     Ne_index = ioff
     ioff = ioff + 1
  endif
  if (requireAl) then
     if(MyPE == 0)  write(*,*)' -- Al'
     nspecies=nspecies+1
     Al_index = ioff
     ioff = ioff + 1
  endif
  if (requireSi) then
      if(MyPE == 0) write(*,*)' -- Si'
     nspecies=nspecies+1
     Si_index = ioff
     ioff = ioff + 1
  endif
  if (requireS)  then
      if(MyPE == 0) write(*,*)' -- S'
     nspecies=nspecies+1
     S_index = ioff
     ioff = ioff + 1
  endif
  if (requireFe) then
      if(MyPE == 0) write(*,*)' -- Fe'
     nspecies=nspecies+1
     Fe_index = ioff
     ioff = ioff + 1
  endif
  !
  allocate(ElementAtomicMass(nspecies))
  !
  if (nion .eq. 0) call abortrun('ERROR: No ions were selected!')
  if (doH1 .and. do_long_spectrum .and. nlyman .gt. nlyman_all)  &
    call abortrun('ERROR: nlyman > nlyman_all')
  if (nlyman .gt. n_lines_max) call abortrun('ERROR: nlyman > n_lines_max')
  if (Lambda_H1(1) .ne. lyalpha)  &
    call abortrun('ERROR: Lambda_H1(1) ne lyalpha')
  !
  allocate(ions(nion))  ! name of element
  allocate(ion_elnr(nion)) ! index in particle data array
  allocate(nlines(nion)) ! number of transitions for this element
  allocate(ion_mass(nion)) ! mass of ion
  allocate(ionfrac(nion))
  allocate(totnr_ion(nion))
  allocate(ion_z_abs(nion)) ! abundance by number, relative to Hydrogen
  allocate(taumax(nion))
  !
  allocate(lambda_rest(nion,n_lines_max)) ! rest wavelenght
  allocate(fosc(nion,n_lines_max)) ! oscillator strength
  !
  ii = 0
  if (doH1) then
    ii = ii + 1
    ions(ii) = 'h1  '
    ion_mass(ii) = massH
    ion_elnr(ii) = H_index
    nlines(ii)   = nlyman
    ion_H        = ii
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_H1(i)
      fosc(ii,i) = f_H1(i)
    enddo
  endif
  if (do21cm) then
    ii = ii + 1
    ions(ii) = '21cm'
    ion_mass(ii) = massH
    ion_elnr(ii) = H_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_21cm(i)
      fosc(ii,i) = f_21cm(i)
    enddo
  endif
  if (doHe2) then
    ii = ii + 1
    ions(ii) = 'he2 '
    ion_mass(ii) = massHe
    ion_elnr(ii) = He_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_He2(i)
      fosc(ii,i) = f_He2(i)
    enddo
  endif
  if (doC2) then
    ii = ii + 1
    ions(ii) = 'c2  '
    ion_mass(ii) = massC
    ion_elnr(ii) = C_index
    nlines(ii) = 2
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_C2(i)
      fosc(ii,i) = f_C2(i)
    enddo
  endif
  if (doC3) then
    ii = ii + 1
    ions(ii) = 'c3  '
    ion_mass(ii) = massC
    ion_elnr(ii) = C_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_C3(i)
      fosc(ii,i) = f_C3(i)
    enddo
  endif
  if (doC4) then
    ii = ii + 1
    ions(ii) = 'c4  '
    ion_mass(ii) = massC
    ion_elnr(ii) = C_index
    nlines(ii) = 2
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_C4(i)
      fosc(ii,i) = f_C4(i)
    enddo
  endif
  if (doC5) then
    ii = ii + 1
    ions(ii) = 'c5  '
    ion_mass(ii) = massC
    ion_elnr(ii) = C_index
    nlines(ii) = 2
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_C5(i)
      fosc(ii,i) = f_C5(i)
    enddo
  endif
  if (doC6) then
    ii = ii + 1
    ions(ii) = 'c6  '
    ion_mass(ii) = massC
    ion_elnr(ii) = C_index
    nlines(ii) = 4
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_C6(i)
      fosc(ii,i) = f_C6(i)
    enddo
  endif
  if (doN2) then
    ii = ii + 1
    ions(ii) = 'n2  '
    ion_mass(ii) = massN
    ion_elnr(ii) = N_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_N2(i)
      fosc(ii,i) = f_N2(i)
    enddo
  endif
  if (doN3) then
    ii = ii + 1
    ions(ii) = 'n3  '
    ion_mass(ii) = massN
    ion_elnr(ii) = N_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_N3(i)
      fosc(ii,i) = f_N3(i)
    enddo
  endif
  if (doN4) then
    ii = ii + 1
    ions(ii) = 'n4  '
    ion_mass(ii) = massN
    ion_elnr(ii) = N_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_N4(i)
      fosc(ii,i) = f_N4(i)
    enddo
  endif
  if (doN5) then
    ii = ii + 1
    ions(ii) = 'n5  '
    ion_mass(ii) = massN
    ion_elnr(ii) = N_index
    nlines(ii) = 2
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_N5(i)
      fosc(ii,i) = f_N5(i)
    enddo
  endif
  if (doN6) then
    ii = ii + 1
    ions(ii) = 'n6  '
    ion_mass(ii) = massN
    ion_elnr(ii) = N_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_N6(i)
      fosc(ii,i) = f_N6(i)
    enddo
  endif
  if (doN7) then
    ii = ii + 1
    ions(ii) = 'n7  '
    ion_mass(ii) = massN
    ion_elnr(ii) = N_index
    nlines(ii) = 2
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_N7(i)
      fosc(ii,i) = f_N7(i)
    enddo
  endif
  if (doO1) then
    ii = ii + 1
    ions(ii) = 'o1  '
    ion_mass(ii) = massO
    ion_elnr(ii) = O_index
    nlines(ii) = 3
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_O1(i)
      fosc(ii,i) = f_O1(i)
    enddo
  endif
  if (doO3) then
    ii = ii + 1
    ions(ii) = 'o3  '
    ion_mass(ii) = massO
    ion_elnr(ii) = O_index
    nlines(ii) = 2
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_O3(i)
      fosc(ii,i) = f_O3(i)
    enddo
  endif
  if (doO4) then
    ii = ii + 1
    ions(ii) = 'o4  '
    ion_mass(ii) = massO
    ion_elnr(ii) = O_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_O4(i)
      fosc(ii,i) = f_O4(i)
    enddo
  endif
  if (doO5) then
    ii = ii + 1
    ions(ii) = 'o5  '
    ion_mass(ii) = massO
    ion_elnr(ii) = O_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_O5(i)
      fosc(ii,i) = f_O5(i)
    enddo
  endif
  if (doO6) then
    ii = ii + 1
    ions(ii) = 'o6  '
    ion_mass(ii) = massO
    ion_elnr(ii) = O_index
    nlines(ii) = 2
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_O6(i)
      fosc(ii,i) = f_O6(i)
    enddo
  endif
  if (doO7) then
    ii = ii + 1
    ions(ii) = 'o7  '
    ion_mass(ii) = massO
    ion_elnr(ii) = O_index
    nlines(ii) = 2
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_O7(i)
      fosc(ii,i) = f_O7(i)
    enddo
  endif
  if (doO8) then
    ii = ii + 1
    ions(ii) = 'o8  '
    ion_mass(ii) = massO
    ion_elnr(ii) = O_index
    nlines(ii) = 4
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_O8(i)
      fosc(ii,i) = f_O8(i)
    enddo
  endif
  if (doMg2) then
    ii = ii + 1
    ions(ii) = 'mg2 '
    ion_mass(ii) = massMg
    ion_elnr(ii) = Mg_index
    nlines(ii) = 2
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_Mg2(i)
      fosc(ii,i) = f_Mg2(i)
    enddo
  endif
  if (doAl2) then
    ii = ii + 1
    ions(ii) = 'al2 '
    ion_mass(ii) = massAl
    ion_elnr(ii) = Al_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_Al2(i)
      fosc(ii,i) = f_Al2(i)
    enddo
  endif
  if (doAl3) then
    ii = ii + 1
    ions(ii) = 'al3 '
    ion_mass(ii) = massAl
    ion_elnr(ii) = Al_index
    nlines(ii) = 2
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_Al3(i)
      fosc(ii,i) = f_Al3(i)
    enddo
  endif
  if (doSi2) then
    ii = ii + 1
    ions(ii) = 'si2 '
    ion_mass(ii) = massSi
    ion_elnr(ii) = Si_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_Si2(i)
      fosc(ii,i) = f_Si2(i)
    enddo
  endif
  if (doSi3) then
    ii = ii + 1
    ions(ii) = 'si3 '
    ion_mass(ii) = massSi
    ion_elnr(ii) = Si_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_Si3(i)
      fosc(ii,i) = f_Si3(i)
    enddo
  endif
  if (doSi4) then
    ii = ii + 1
    ions(ii) = 'si4 '
    ion_mass(ii) = massSi
    ion_elnr(ii) = Si_index
    nlines(ii) = 2
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_Si4(i)
      fosc(ii,i) = f_Si4(i)
    enddo
  endif
  if (doS5) then
    ii = ii + 1
    ions(ii) = 's5  '
    ion_mass(ii) = massS
    ion_elnr(ii) = S_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_S5(i)
      fosc(ii,i) = f_S5(i)
    enddo
  endif
  if (doNe8) then
    ii = ii + 1
    ions(ii) = 'ne8 '
    ion_mass(ii) = massNe
    ion_elnr(ii) = Ne_index
    nlines(ii) = 2
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_Ne8(i)
      fosc(ii,i) = f_Ne8(i)
    enddo
  endif
  if (doNe9) then
    ii = ii + 1
    ions(ii) = 'ne9 '
    ion_mass(ii) = massNe
    ion_elnr(ii) = Ne_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_Ne9(i)
      fosc(ii,i) = f_Ne9(i)
    enddo
  endif
  if (doFe2) then
    ii = ii + 1
    ions(ii) = 'fe2 '
    ion_mass(ii) = massFe
    ion_elnr(ii) = Fe_index
    nlines(ii) = 9
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_Fe2(i)
      fosc(ii,i) = f_Fe2(i)
    enddo
  endif
  if (doFe3) then
    ii = ii + 1
    ions(ii) = 'fe3 '
    ion_mass(ii) = massFe
    ion_elnr(ii) = Fe_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_Fe3(i)
      fosc(ii,i) = f_Fe3(i)
    enddo
  endif
  if (doFe17) then
    ii = ii + 1
    ions(ii) = 'fe17'
    ion_mass(ii) = massFe
    ion_elnr(ii) = Fe_index
    nlines(ii) = 2
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_Fe17(i)
      fosc(ii,i) = f_Fe17(i)
    enddo
  endif
  if (doFe19) then
    ii = ii + 1
    ions(ii) = 'fe19'
    ion_mass(ii) = massFe
    ion_elnr(ii) = Fe_index
    nlines(ii) = 2
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_Fe19(i)
      fosc(ii,i) = f_Fe19(i)
    enddo
  endif
  if (doFe21) then
    ii = ii + 1
    ions(ii) = 'fe21'
    ion_mass(ii) = massFe
    ion_elnr(ii) = Fe_index
    nlines(ii) = 1
    do i = 1, nlines(ii)
      lambda_rest(ii,i) = Lambda_Fe21(i)
      fosc(ii,i) = f_Fe21(i)
    enddo
  endif
  !
  if (requireH)  ElementAtomicMass(H_index)    = massH
  if (requireHe) ElementAtomicMass(He_index)   = massHe
  if (requireC) ElementAtomicMass(C_index)     = massC
  if (requireN) ElementAtomicMass(N_index)     = massN
  if (requireO) ElementAtomicMass(O_index)     = massO
  if (requireNe) ElementAtomicMass(Ne_index)   = massNe
  if (requireMg) ElementAtomicMass(Mg_index)   = massMg
  if (requireAl) ElementAtomicMass(Al_index)   = massAl
  if (requireSi) ElementAtomicMass(Si_index)   = massSi
  if (requireS) ElementAtomicMass(S_index)     = massS
  if (requireFe) ElementAtomicMass(Fe_index)   = massFe
  !
  !     ----------------------------------
  !     Read in ionization balance tables.
  !     ----------------------------------
  if(verbose .and. mype == 0)then
    write (*,*) ' +++++++++ '
    write (*,*) ' ibdir = ',trim(ibdir)
    write (*,'('' loading ionization tables from: '',a)') trim(adjustl(ibdir))
  endif
  do ion=1, nion
    if (read_ionbal_from_single_file) then
      call load_ionbal_singlefile(ion)
    else
      call load_ionbal(ion)
    endif
    !
    j = 0
    do ib = 1, nd
      do it=1,nt
        do iz=1,nz
          j = j + 1
          ionizbal(ion,j) = ionbal(iz,it,ib)
        enddo
      enddo
    enddo
  enddo
  if(verbose .and. mype == 0)then
    write (*,*) ' +++++++++ '
    write (*,*) 
  endif
end subroutine initialize_ionization_tables


!This subroutine loads in ionization balance for JS's separate HDF5 files
subroutine load_ionbal(ion)
  ! Reads in ionization balance interpolation table.
  ! Table is 3-D : redshift, log(T) and log(n_H)
  use hdf5_wrapper
  use ionization_tables
  use runtime
  use spectra
  use my_mpi
  implicit none
  !
  integer, intent(in) :: ion
  !
  integer iz, it, id, file_handle
  character(len=200) :: longfile, VarName, filename
  character(len=50)  :: dtype
  integer            :: size, rank
  integer            :: dims(10)
  !
  logical, save :: first_call = .true.
  integer(kind=singleI), save :: nd_old, nt_old, nz_old
  ! 
  !If the ion is called '21cm' then force that we load the H1 ionization table
  if (ions(ion) .eq. '21cm') then
    filename = trim(ibdir)//trim('h1')//'.hdf5'
  else
    filename = trim(ibdir)//trim(ions(ion))//'.hdf5'
  endif
  call hdf5_open_file(file_handle, filename, readonly=.true.)
  !
  if(first_call)then
     first_call = .false.
     ! get size of each variable
     VarName = 'logd'
     call hdf5_get_dimensions(file_handle,VarName,rank,dims)
     nd = dims(1)
     VarName = 'logt'
     call hdf5_get_dimensions(file_handle,VarName,rank,dims)
     nt = dims(1)
     VarName = 'redshift'
     call hdf5_get_dimensions(file_handle,VarName,rank,dims)
     nz = dims(1)
     nd_old = nd
     nt_old = nt
     nz_old = nz
     nznt   = nz*nt
     !
     allocate(ib_redshift(nz))
     allocate(ib_logt(nt))
     allocate(ib_logd(nd))
     allocate(ionbal(nz,nt,nd))
     allocate(ionizbal(nion,nz*nd*nt))
     !
     VarName = 'logd'
     call hdf5_read_data(file_handle,VarName,ib_logd)
     VarName = 'logt'
     call hdf5_read_data(file_handle,VarName,ib_logt)
     VarName = 'redshift'
     call hdf5_read_data(file_handle,VarName,ib_redshift)
  else
     VarName = 'logd'
     call hdf5_get_dimensions(file_handle,VarName,rank,dims)
     nd = dims(1)
     VarName = 'logt'
     call hdf5_get_dimensions(file_handle,VarName,rank,dims)
     nt = dims(1)
     VarName = 'redshift'
     call hdf5_get_dimensions(file_handle,VarName,rank,dims)
     nz = dims(1)
     if(nz .ne. nz_old .or. nd .ne. nd_old .or. nt .ne. nt_old)then
        call abortrun(' error reading ionization tables')
     endif
  endif
  !
  VarName = '/ionbal'
  call hdf5_read_data(file_handle,VarName,ionbal(1:nz,1:nt,1:nd))
  call hdf5_close_file(file_handle)     
  !
  if (verbose .and. MyPE == 0) write(*,*) ions(ion), ' Ionization table loaded.'
  !  
  do id = 1, nd
     do it=1,nt
        do iz=1,nz
           if (ionbal(iz,it,id) .eq. 0.) then
              ionbal(iz,it,id) = -32.
           else
              ionbal(iz,it,id) = log10(ionbal(iz,it,id))
           endif
        enddo
     enddo
  enddo
  !
end subroutine load_ionbal


!This subroutine loads in ionization balance for TT's combined HDF5 file
subroutine load_ionbal_singlefile(ion)
  !
  ! Reads in ionization balance interpolation table.
  ! Table is 3-D : redshift, log(T) and log(n_H)
  use hdf5_wrapper
  use ionization_tables
  use runtime
  use spectra
  use my_mpi
  implicit none
  integer, intent(in) :: ion
  !
  integer iz, it, id, file_handle
  character(len=200) :: longfile, VarName, filename
  character(len=50)  :: dtype
  integer            :: size, rank
  integer            :: dims(10)
  !
  logical, save :: first_call = .true.
  integer(kind=singleI), save :: nd_old, nt_old, nz_old
  ! 
  filename = trim(ibdir)//'hm01_Q+G_tables.hdf5'
  call hdf5_open_file(file_handle, filename, readonly=.true.)
  !
  if(first_call)then
     first_call = .false.
     ! get size of each variable
     VarName = 'LogDensity'
     call hdf5_get_dimensions(file_handle,VarName,rank,dims)
     nd = dims(1)
     VarName = 'LogTemperature'
     call hdf5_get_dimensions(file_handle,VarName,rank,dims)
     nt = dims(1)
     VarName = 'Redshift'
     call hdf5_get_dimensions(file_handle,VarName,rank,dims)
     nz = dims(1)
     nd_old = nd
     nt_old = nt
     nz_old = nz
     nznt   = nz*nt
     !
     allocate(ib_redshift(nz))
     allocate(ib_logt(nt))
     allocate(ib_logd(nd))
     allocate(ionbal(nz,nt,nd))
     allocate(ionizbal(nion,nz*nd*nt))
     !
     VarName = 'LogDensity'
     call hdf5_read_data(file_handle,VarName,ib_logd)
     VarName = 'LogTemperature'
     call hdf5_read_data(file_handle,VarName,ib_logt)
     VarName = 'Redshift'
     call hdf5_read_data(file_handle,VarName,ib_redshift)
  else
     VarName = 'LogDensity'
     call hdf5_get_dimensions(file_handle,VarName,rank,dims)
     nd = dims(1)
     VarName = 'LogTemperature'
     call hdf5_get_dimensions(file_handle,VarName,rank,dims)
     nt = dims(1)
     VarName = 'Redshift'
     call hdf5_get_dimensions(file_handle,VarName,rank,dims)
     nz = dims(1)
     if(nz .ne. nz_old .or. nd .ne. nd_old .or. nt .ne. nt_old)then
        write (*,*) ' error reading ionization tables'
        stop
     endif
  endif
  !
  VarName = '/'//trim(ions(ion))//'/IonizationFraction'
  call hdf5_read_data(file_handle,VarName,ionbal(1:nz,1:nt,1:nd))
  call hdf5_close_file(file_handle)     
  !
  if (verbose .and. MyPE == 0) write(*,*) ions(ion), ' Ionization table loaded.'
  !  
  do id = 1, nd
     do it=1,nt
        do iz=1,nz
           if (ionbal(iz,it,id) .eq. 0.) then
              ionbal(iz,it,id) = -32.
           else
              ionbal(iz,it,id) = log10(ionbal(iz,it,id))
           endif
        enddo
     enddo
  enddo
  !
end subroutine load_ionbal_singlefile



function get_fitted_ibfactor(z)
  use numbers
  use spectra, only: ibfactor_he_reionization
  implicit none
  !
  !properties of x^2 e(-x) fit:
  real(kind=doubleR),parameter :: a = 0.856953
  real(kind=doubleR),parameter :: h = 1.38134
  !
  !properties of polynomial fit:
  real(kind=doubleR),parameter :: p1 = -2.07639
  real(kind=doubleR),parameter :: p2 = 2.62640
  real(kind=doubleR),parameter :: p3 = -0.734
  real(kind=doubleR),parameter :: p4 = 0.0633406
  !
  !properties for He reionization:
  real(kind=doubleR),parameter :: hea = 0.244929
  real(kind=doubleR),parameter :: heb = 3.25647
  real(kind=doubleR),parameter :: hec = 0.0966662
  !
  real(kind=doubleR) :: get_fitted_ibfactor
  real(kind=singleR) :: z
  !
  if (z .lt. 2) z = 2
  if (z .gt. 5) z = 5
  !
  !get_fitted_ibfactor = a * z**2 * exp (-z/h)
  get_fitted_ibfactor = p1 + z*p2 + z**2*p3 + z**3*p4
  !
  !Do we want to match the 'Helium reionization' peak in the data of
  !Faucher-Giguere (2007)?
  if (ibfactor_he_reionization) then
    get_fitted_ibfactor = get_fitted_ibfactor - hea * exp(-(z-heb)**2/hec**2)
  endif
  !
end function get_fitted_ibfactor


subroutine computeib(z1,z2,iz1,iz2,dz1,dz2,logtemp,logdens,fraction)
  use numbers
  use ionization_tables
  use runtime
  use spectra, only : nion, ibfactor,use_fitted_ibfactor
  implicit none
  integer(kind=singleI), intent(in) :: iz1, iz2
  real(kind=singleR), intent(in)    :: z1,z2
  real(kind=doubleR), intent(in)    :: dz1,dz2,logtemp, logdens
  real(kind=doubleR), intent(out)   :: fraction(nion)

  ! local variables
  integer i, it1, it2, id1, id2
  integer i111, i211, i121, i112, i122, i212, i221, i222 
  real(kind=doubleR) ::  w111, w211, w121, w112, w122, w212, w221, w222 
  real(kind=doubleR) ::  dd1, dd2, dt1, dt2, logd, logt,ibfactor_use
  real(kind=doubleR)  :: get_fitted_ibfactor
  external get_fitted_ibfactor
  !
  !Interpolate from CLOUDY table (trilinear interpolation).
  !Compute interpolation weights
  !
  !
  !Rescale ionizing background by factor ibfactor.
  !(instead we are actually scaling the density by 1/ibfactor)

  if (use_fitted_ibfactor) then
     ibfactor_use = 1./(get_fitted_ibfactor(z1)*dz1 + get_fitted_ibfactor(z2)*dz2)
  else
     ibfactor_use = ibfactor
  endif

  logd = logdens - log10(ibfactor_use)
  
  !    Bring temperature within range of ionization table.
  if (logtemp .lt. ib_logt(1)) then 
     logt = ib_logt(1)
  else if (logtemp .gt. ib_logt(nt)) then 
     logt = ib_logt(nt)
  else 
     logt = logtemp
  endif

!     Bring density within range of ionization table.
  if (logd .lt. ib_logd(1)) then 
     write(*,*)'WARNING! low'
     logd = ib_logd(1)
  else if (logd .gt. ib_logd(nd)) then 
    write(*,*)'WARNING! high'
     logd = ib_logd(nd)
  endif

!     If we are forcing that we use the highest density above zmax then do this here  
  if (use_maxdens_above_zmax .and. abs(z2-ib_redshift(nz)) .lt. 1e-2 .and. dz2 .eq. 1.0) then
     logd = ib_logd(nd)
  endif

  it2 = 2
  do while (logt .gt. ib_logt(it2))
     it2 = it2 + 1
  enddo
  it1 = it2 - 1
  dt1 = (ib_logt(it2) - logt) / (ib_logt(it2) - ib_logt(it1))
  dt2 = 1. - dt1
  
  id2 = 2
  do while (logd .gt. ib_logd(id2))
     id2 = id2 + 1
  enddo
  id1 = id2 - 1
  dd1 = (ib_logd(id2) - logd) / (ib_logd(id2) - ib_logd(id1))
  dd2 = 1. - dd1
  
!     Weights:
  w111 = dz1 * dt1 * dd1
  w211 = dz2 * dt1 * dd1
  w121 = dz1 * dt2 * dd1
  w221 = dz2 * dt2 * dd1
  w112 = dz1 * dt1 * dd2
  w212 = dz2 * dt1 * dd2
  w122 = dz1 * dt2 * dd2
  w222 = dz2 * dt2 * dd2
  
!     Indices:
  i111 = iz1 + (it1-1)*nz + (id1-1)*nznt
  i211 = iz2 + (it1-1)*nz + (id1-1)*nznt
  i121 = iz1 + (it2-1)*nz + (id1-1)*nznt
  i221 = iz2 + (it2-1)*nz + (id1-1)*nznt
  i112 = iz1 + (it1-1)*nz + (id2-1)*nznt
  i212 = iz2 + (it1-1)*nz + (id2-1)*nznt
  i122 = iz1 + (it2-1)*nz + (id2-1)*nznt
  i222 = iz2 + (it2-1)*nz + (id2-1)*nznt

  !
  fraction(:) =  &
          w111*ionizbal(:,i111) + w211*ionizbal(:,i211) + &
          w121*ionizbal(:,i121) + w221*ionizbal(:,i221) + &
          w112*ionizbal(:,i112) + w212*ionizbal(:,i212) + &
          w122*ionizbal(:,i122) + w222*ionizbal(:,i222)
  fraction = 10.d0**fraction

  return
end subroutine computeib


! +++++++++++++++++++++++++++++++++++++ ...ionization ++++++++++++++++++++++++++++++++++++++++++++ !

! ++++++++++++++++++++++++++++++++++++++++++ general... ++++++++++++++++++++++++++++++++++++++++++ !


subroutine mylabel(string,value,label)
  use numbers
  implicit none
  !
  character(*), intent(in)          :: string
  character(*), intent(out)         :: label
  integer(kind=singleI), intent(in) :: value
  ! local variables
  character(len=120) :: name
  !
  if(value .lt. 10)then
    write (name,'(i1)') value
  else if (value .lt. 100) then
    write (name,'(i2)') value
  else if (value .lt. 1000) then
    write (name,'(i3)') value
  else if (value .lt. 10000) then
    write (name,'(i4)') value
  else if (value .lt. 100000) then
    write (name,'(i5)') value
  else
    write (*,*) 'sorry: value too large'
    stop
  endif
  label = trim(string)//trim(adjustl(name))
  !
end subroutine mylabel


function my_cpu_time()
  use numbers
  use my_mpi
  implicit none
  !
  real(kind=doubleR) :: my_cpu_time
  !
#ifdef MPI
  my_cpu_time = mpi_wtime()
#else
  my_cpu_time = 0.0
#endif
  !
end function my_cpu_time


subroutine cpu_timer_start(timer)
  use cpu_timers
  implicit none
  !
  integer, intent(in) :: timer
  real(kind=doubleR)  :: my_cpu_time
  external my_cpu_time
  !
  if(timer .le. 0 .or. timer .gt. ntimers) &
    call abortrun(' illegal cpu timer')
  if(cpu_is_running(timer)) &
    call abortrun(' cpu timer is already running')
  cpustart(timer) = my_cpu_time()
  cpu_is_running(timer) = .true.
  !
end subroutine cpu_timer_start


subroutine cpu_timer_stop(timer)
  use cpu_timers
  implicit none
  !
  integer, intent(in) :: timer
  real(kind=doubleR)  :: my_cpu_time
  external my_cpu_time
  !
  if(timer .le. 0 .or. timer .gt. ntimers) &
       call abortrun(' illegal cpu timer')
  if(.not. cpu_is_running(timer)) &
       call abortrun(' cpu timer is not running')
  cpustop(timer) = my_cpu_time()
  cputime(timer) = cputime(timer) + cpustop(timer) - cpustart(timer)
  !
  cpu_is_running(timer) = .false.
  cpustart(timer)       = 0.0
  !
end subroutine cpu_timer_stop


subroutine cpu_timer_reset(timer)
  use cpu_timers
  implicit none
  !
  integer, intent(in) :: timer
  real(kind=doubleR) :: cpu_time
  !
  !  if(timer .le. 0 .or. timer .gt. ntimers) &
       call abortrun(' illegal cpu timer')
  cputime(timer)  = 0.0
  cpustart(timer) = 0.0
  cpustop(timer)  = 0.0
  cpu_is_running(timer) = .false.
  !
end subroutine cpu_timer_reset


subroutine initialize_cputimers
  use cpu_timers
  implicit none
  !
  cputime        = 0.0
  cpustart       = 0.0
  cpustop        = 0.0
  cpu_is_running = .false.
  !
end subroutine initialize_cputimers


! +++++++++++++++++++++++++++++++++++++++ ...general +++++++++++++++++++++++++++++++++++++++++++++ !

! ++++++++++++++++++++++++++++++++++++++++ mathematical... +++++++++++++++++++++++++++++++++++++++ !


subroutine convlv(data,n,respns,m,isign,ans)
  use numbers
  use spectra, only : fft
  implicit none
  !
  integer(kind=singleI), intent(in) :: n, m, isign
  real(kind=doubleR), intent(in)    :: data(n)
  real(kind=doubleR), intent(inout) :: respns(n)
  complex(kind=doubleR), intent(out):: ans(2*n)
  ! local variables
  INTEGER i,no2
  integer(kind=singleI), save :: nfft=-1
  !
  if(nfft .ne. n)then
    allocate(fft(n))
    nfft = n
  endif
  !  
  do  i=1,(m-1)/2
    respns(n+1-i)=respns(m+1-i)
  enddo
  do  i=(m+3)/2,n-(m-1)/2
    respns(i)=0.0
  enddo
  !
  call twofft(data,respns,fft,ans,n)
  no2=n/2
  do i=1,no2+1
    if (isign.eq.1) then
      ans(i)=fft(i)*ans(i)/no2
    else if (isign.eq.-1) then
      if (abs(ans(i)).eq.0.0) stop &
          'deconvolving at response zero in convlv'
      ans(i)=fft(i)/ans(i)/no2
    else
      stop 'no meaning for isign in convlv'
    endif
  enddo
  ans(1)=dcmplx(dble(ans(1)),dble(ans(no2+1)))
  call realft(ans,n,-1)
  !
  return
  !
end subroutine convlv


subroutine twofft(data1,data2,fft1,fft2,n)
  use numbers
  implicit none
  !
  integer(kind=singleI), intent(in) :: n
  real(kind=doubleR), intent(in)    :: data1(n), data2(n)
  complex(kind=doubleR)             :: fft1(n), fft2(n)
  !
  ! local variables
  INTEGER j,n2
  COMPLEX(kind=doubleR) ::  h1,h2,c1,c2
  !
  c1=cmplx(0.5,0.0)
  c2=cmplx(0.0,-0.5)
  do  j=1,n
    fft1(j)=cmplx(data1(j),data2(j))
  enddo
  call four1(fft1,n,1)
  fft2(1)=dcmplx(imag(fft1(1)),0d0)
  fft1(1)=dcmplx(dble(fft1(1)),0d0)
  n2=n+2
  do j=2,n/2+1
    h1=c1*(fft1(j)+conjg(fft1(n2-j)))
    h2=c2*(fft1(j)-conjg(fft1(n2-j)))
    fft1(j)=h1
    fft1(n2-j)=conjg(h1)
    fft2(j)=h2
    fft2(n2-j)=conjg(h2)
  enddo
  return
  !
end subroutine twofft


subroutine realft(data,n,isign)
  use numbers
  implicit none
  !
  integer(kind=singleI), intent(in) :: n, isign
  real(kind=doubleR), intent(inout) :: data(n)
  ! local variables
  INTEGER i,i1,i2,i3,i4,n2p3
  REAL*8 c1,c2,h1i,h1r,h2i,h2r,wis,wrs
  REAL*8 theta,wi,wpi,wpr,wr,wtemp
  theta=3.141592653589793d0/dble(n/2)
  c1=0.5
  if (isign.eq.1) then
     c2=-0.5
     call four1(data,n/2,+1)
  else
     c2=0.5
     theta=-theta
  endif
  wpr=-2.0d0*sin(0.5d0*theta)**2
  wpi=sin(theta)
  wr=1.0d0+wpr
  wi=wpi
  n2p3=n+3
  do  i=2,n/4
     i1=2*i-1
     i2=i1+1
     i3=n2p3-i2
     i4=i3+1
     wrs=sngl(wr)
     wis=sngl(wi)
     h1r=c1*(data(i1)+data(i3))
     h1i=c1*(data(i2)-data(i4))
     h2r=-c2*(data(i2)+data(i4))
     h2i=c2*(data(i1)-data(i3))
     data(i1)=h1r+wrs*h2r-wis*h2i
     data(i2)=h1i+wrs*h2i+wis*h2r
     data(i3)=h1r-wrs*h2r+wis*h2i
     data(i4)=-h1i+wrs*h2i+wis*h2r
     wtemp=wr
     wr=wr*wpr-wi*wpi+wr
     wi=wi*wpr+wtemp*wpi+wi
  enddo
  if (isign.eq.1) then
     h1r=data(1)
     data(1)=h1r+data(2)
     data(2)=h1r-data(2)
  else
     h1r=data(1)
     data(1)=c1*(h1r+data(2))
     data(2)=c1*(h1r-data(2))
     call four1(data,n/2,-1)
  endif
  return
  !
end subroutine realft


subroutine four1(data,nn,isign)
  use numbers
  implicit none
  !
  integer(kind=singleI), intent(in) :: nn, isign
  real(kind=doubleR), intent(inout) :: data(2*nn)
  ! local variables  
  INTEGER i,istep,j,m,mmax,n
  REAL(kind=doubleR) :: tempi,tempr
  REAL(kind=doubleR) ::  theta,wi,wpi,wpr,wr,wtemp
  !
  n=2*nn
  j=1
  do  i=1,n,2
     if(j.gt.i)then
        tempr=data(j)
        tempi=data(j+1)
        data(j)=data(i)
        data(j+1)=data(i+1)
        data(i)=tempr
        data(i+1)=tempi
     endif
     m=n/2
1    if ((m.ge.2).and.(j.gt.m)) then
        j=j-m
        m=m/2
        goto 1
     endif
     j=j+m
  enddo
  mmax=2
2 if (n.gt.mmax) then
     istep=2*mmax
     theta=6.28318530717959d0/(isign*mmax)
     wpr=-2.d0*sin(0.5d0*theta)**2
     wpi=sin(theta)
     wr=1.d0
     wi=0.d0
     do m=1,mmax,2
        do  i=m,n,istep
           j=i+mmax
           tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
           tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
           data(j)=data(i)-tempr
           data(j+1)=data(i+1)-tempi
           data(i)=data(i)+tempr
           data(i+1)=data(i+1)+tempi
        enddo
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
     enddo
     mmax=istep
     goto 2
  endif
  return
  !
end subroutine four1

! +++++++++++++++++++++++++++++++++++++ ...mathematical +++++++++++++++++++++++++++++++++++++++++++ !

