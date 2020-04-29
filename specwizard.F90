! to check: is ionization balance determined in terms of rho/MassH, or n_H?

program specwizard
  
  ! usage: mpirun -np 2 ./specwizard parameter_file
  ! input:       filename   : string  : optional  : contains parameter for this run
  !
  !      All simulation spectra are calculated in velocity space. For each
  !      element the velocity range is the same and can be converted to 
  !      observed wavelength space, or equivalently redshift space.
  !
  !      The long spectrum is computed in dimensionless velocity space (v/c),
  !      so that it can be convolved with the instrumental profile.
  !
  !      I make frequent use of the following relations (v_min = 0):
  !      v = c * ln(lambda/lambda_min) = c * ln[(1+z)/(1+z_min)]
  !      lambda = lambda_min * exp(v/c) = lambda_min * (1+z)/(1+z_min)
  !      1 + z = (1+z_min) * lambda/lambda_min = (1+z_min) * exp(v/c)
  !
  !      A simulated spectrum tau(v_sim) at redshift z_sim, for a transition 
  !      with rest wavelength lambda_0, can be mapped onto the long spectrum,
  !      tautot(v) using the following relation:
  !      v/c = ln(1+z_sim) + v_sim/c + ln(lambda_0/lambda_min),
  !      where we assumed v(lambda_min) = 0, v_sim(z_sim) = 0.
  !      Proof:
  !      v/c = ln(lambda/lambda_min) = ln[(1+z)*lambda_0 / lambda_min]
  !          = ln[(1+z_sim)*exp(v_sim/c) * lambda_0/lambda_min]
  !  
  use numbers
  use spectra
  use ionization_tables
  use runtime
  use my_mpi
  use particledata
  use projection_parameters
  use header
  use modified_metallicity, only: modify_metallicity
  use random_numbers, only : ran_file, ran_sight, ran_noise, ran_metal, ran_fill, seed, ran_los
  use my_random_numbers
  use cpu_timers
  use read_startup_file, only : read_parameter_file, perform_sanity_check
  use get_argumentsmod
  use noisedata
#ifdef READREGION
  use RegionExtent
#endif
  !
  implicit none 
  !
  ! local variables
  integer(kind=singleI) :: ismallspec, nlos_possible, nok, i
  integer(kind=singleI) :: ifile, los_number, ierr, iproc
  real(kind=doubleR)    :: z0, z1, CurrentHubbleCt, zcurrent_next
  integer(kind=singleI) :: this_spec, nargs, numspec, n1d
  real(kind=doubleR)    :: extent
  !
#ifdef MPI
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, MyPE, ierr)
  call mpi_comm_size(mpi_comm_world, NumPEs, ierr)
#else
  MyPE   = 0
  NumPEs = 1
#endif  
  !  
  if(MyPE == 0)then
    write(*,*)
    write(*,*) '*******************************'
    write(*,*) '*         SPECWIZARD          *'
    write(*,*) '*                             *'  
    write(*,*) '*             by              *'
    write(*,*) '*         Joop Schaye,        *'  
    write(*,*) '*         Tom Theuns,         *'     
    write(*,*) '*         Craig Booth         *'  
    write(*,*) '*                             *'  
    write(*,*) '*  schaye@strw.leidenuniv.nl  *'
    write(*,*) '*   tom.theuns@durham.ac.uk   *'
    write(*,*) '*   booth@strw.leidenuniv.nl  *'
    write(*,*) '*******************************'
    write(*,*)
    write (*,*)'* running on ',NumPEs,' processors'
  endif
  ! 
  ! input:        filename : string : optional
  nargs = get_nargs()
  if(nargs .gt. 1) then
    call abortrun('usage: ./specwizard parameterfile (optional, defaults to specwizard.par)')
  endif
  !
  if(nargs == 1) call get_argument(1,parameter_file)
  !
  ! read parameter file
  call read_parameter_file()
  !
  ! Perform sanity check on variables
  call perform_sanity_check()
  !
  ! set cpu timers
  call initialize_cputimers()
  call cpu_timer_start(main)
  !
  ! read ionization tables and set element names
  call initialize_ionization_tables()
  !
  ! set spectral parameters
  call initialize_spectral_parameters()
  !
  ! create output file and write attributes
  if(do_long_spectrum) &
       call create_spectrum_file(1)
  !
#ifdef MPI
    call mpi_barrier(mpi_comm_world, ierr)
#endif
  !
  docontspec: if (do_long_spectrum) then
    !
    do_spectra: do this_spec=1, nspec, NumPEs
      ispec = this_spec + MyPE
      ! 
      good_spectrum: if(ispec .le. nspec)then
        call cpu_timer_start(dospectra)
        !
        call zero_spectra()
        nsimfile_used = 0
        !
        call init_random_numbers(ispec)
        !
        if(verbose .and. ispec == 1)then
          write (*,*) ' +++++++++++++++++++++++++++++++++++'
          write (*,'("generating spectrum from zmin = ",F7.3," to zmax = ",F7.3)')zabsmin,zabsmax
        endif
        !     
        ismallspec      = 0
        nlos_possible   = 0
        zcurrent        = zstart
        zcurrentloop: do while (zcurrent .lt. zend)
          !
          z0 = max(zcurrent - fzresol*(1.+zcurrent)/2.,0.)
          z1 = z0 + fzresol*(1.+zcurrent)
          !
          if (verbose .and. MyPE == 0) then
            write(*,'("---------------------------------------")')
            write(*,*) 
            write(*,'("zcurrent, zmin, zmax = ",3(f6.3,x))') zcurrent, z0, z1
          endif
          !         
          nok = 0
          do i = 1, nlosfiles
            if (simz(i) .ge. z0 .and. simz(i) .le. z1) then 
              nok = nok + 1
              ichoose(nok) = i
              nlos_possible = nlos_possible + nlos_in_file(i)
            endif
          enddo
          !         
          if (nok .eq. 0) then 
            write(*,'("ERROR: There are no spectra in the range z = ", f6.3," - ",f6.3)') z0, z1
            call abortrun('stop')
          endif
          !
          ! randomly choose file between 1 and nok
          ifile = int(random(ran_file)*nok)+1
          ifile = ichoose(ifile)
          !
          ! randomly choose sightline
          los_number = int(random(ran_sight)*nlos_in_file(ifile))
          call readdata_owls(simfile(ifile),los_number)
          !
          if(modify_metallicity) call impose_metallicity()
          !
          if(impose_eos) call impose_equation_of_state()
          !
          ! Scale to zcurrent:
          acurrent = 1. / (1. + zcurrent)
          ! Compute Hubble parameter (in km/s/Mpc).
          CurrentHubbleCt = 100. * HubbleParam * sqrt(1. + omega0*(1./acurrent-1.) &
            + OmegaLambda * (acurrent**2-1.)) /acurrent
          boxkms = BoxSize / HubbleParam * acurrent * CurrentHubbleCt
          !
          ! Number of pixels in this spectrum
          nveloc = int(boxkms / vpixsizekms) + 1
          !
          if(verbose .and. MyPE == 0)then
            write(*,'("acurrent:        ",f11.4)')acurrent
            write(*,'("a_sim:           ",f11.4)')ExpansionFactor
            if (.not. use_snapshot_file) write(*,'("los_number:	 ",i11)')los_number
          endif
          !
          ! project SPH data to sightline
          call cpu_timer_start(doprojecteach)
          call projectdata()
          call cpu_timer_stop(doprojecteach)
          !
          ! Compute individual spectrum for each ion
          call cpu_timer_start(domakespectraeach)
          call makespectra()
          call cpu_timer_stop(domakespectraeach)
          !
          ! insert in resulting spectrum
          call cpu_timer_start(doinsert)
          call insertspectra(zcurrent_next)
          call cpu_timer_stop(doinsert)
          ! report single-snapshot spectrum timing
          if(verbose) &
                write (*,'('' projectdata time = '',e12.4,'' s, makespectra time = '',e12.4,&
                           '' s, insertspectra time = '',e12.4,'' s'')') & 
                cputime(doprojecteach), cputime(domakespectraeach), cputime(doinsert)
          !
          call store_spectrum_info(simfile(ifile), los_number)
          !
          zcurrent = zcurrent_next
          !
        enddo zcurrentloop
        !
        ! compute net spectrum
        if(nion == 1)then
          flux(:) = exp(-tau_long(1,:))
        else
          flux(:) = exp(-sum(tau_long,1))
        endif
        !
        ! convolve with instrumental profile
        if (do_convolve_spectrum) then
          call convolve_long_spectrum()
        else
          flux_convolved = flux
        endif
        !
        ! rebin to specified pixel size ('pixsize' angstroms)
        call rebin_spectrum()
        !
        ! add noise
        if(generate_noise) call add_noise()
        !
        call cpu_timer_stop(dospectra)
        ! report time for long spectrum
        if(verbose) &
                write (*,'('' long spectrum time = '',e12.4,'' s'')') & 
                cputime(dospectra)
        !
      endif good_spectrum
      !
      ! Output spectrum
      !
      do iproc=0, NumPEs-1
        if(iproc == MyPE .and. ispec .le. nspec)then
          write (*,*) ' MyPE = ',MyPE,' outputting spectrum = ',ispec
          call write_long_spectrum()
        endif
#ifdef MPI
        call mpi_barrier(mpi_comm_world, ierr)
#endif
      enddo
      !
    enddo do_spectra
    !
  else !docontspec
    !
    files_loop: do ifile=1,nlosfiles
      !
      if(verbose .and. MyPE == 0)  then
        write (*,*) ' file = ',ifile,' out of ',nlosfiles
        call flush(6)
      endif
      !
      if(use_snapshot_file) then
        call get_los_coordinates()
        numspec = nspec
      else
         if(nspec .gt. 0) then
            numspec = minval( (/nlos_in_file(ifile),nspec/) )
         else
            numspec = nlos_in_file(ifile) ! do all available spectra
         endif
         nspec = numspec
      endif
      !
      if(int(numspec/NumPes)*NumPEs /= numspec) then
        write (*,*) ' sorry: number of spectra to compute *MUST* be multiple of number of processors'
        write (*,*) ' this is not the case for file = ',ifile,' which has ',numspec,' spectra'
        call abortrun('stop')
      endif
      !
      if(use_snapshot_file) then
        !
#ifndef READREGION
        call read_full_snapshot()
        !
        if(modify_metallicity) call impose_metallicity()
        !
        if(impose_eos) call impose_equation_of_state()
#endif
        !
      endif
      !
      los_loop: DO los_number=MyPE, numspec-1, NumPEs
        !
        if(verbose .and. MyPE == 0) &
          write (*,*) ' sightline = ',los_number+1,' out of ', numspec
        !
        ispec = los_number
        !
        if (.not. use_snapshot_file) then
          !
          call readdata_owls(simfile(ifile), los_number)
          !
          if(modify_metallicity) call impose_metallicity()
          !
          if(impose_eos) call impose_equation_of_state()
          !
        else ! data already read above
           x_fraction = x_fraction_array(los_number+1)
           y_fraction = y_fraction_array(los_number+1)
           z_fraction = z_fraction_array(los_number+1)
           !
           x_comoving = x_fraction * BoxSize
           y_comoving = y_fraction * BoxSize
           z_comoving = z_fraction * BoxSize
           !
           x_physical = x_fraction * BoxSize / HubbleParam * Expansionfactor ! Physical Mpc
           y_physical = y_fraction * BoxSize / HubbleParam * Expansionfactor ! Physical Mpc
           z_physical = z_fraction * BoxSize / HubbleParam * Expansionfactor ! Physical Mpc
           
           
#ifdef READREGION
           ! set region extent to be read

           ! as a fraction of the box size
           N1D    = nint( float(Numpart_Total(1))**(1d0/3d0) )
           extent = 1. / float(N1d) ! seperation at mean density as fraction of box size

           RegionExtentX(1) = x_fraction - 8d0*extent
           RegionExtentX(2) = x_fraction + 8d0*extent

           RegionExtentY(1) = y_fraction - 8d0*extent
           RegionExtentY(2) = y_fraction + 8d0*extent


           ! RegionExtentX(1) = 0
           ! RegionExtentX(2) = BoxSize

           ! RegionExtentY(1) = 0
           ! RegionExtentY(2) = BoxSize

           RegionExtentZ(1) = 0
           RegionExtentZ(2) = 1d0

           ! in co-mving units
           RegionExtentX = RegionExtentX * BoxSize
           RegionExtentY = RegionExtentY * BoxSize
           RegionExtentZ = RegionExtentZ * BoxSize
           call read_full_snapshot()

           if(modify_metallicity) call impose_metallicity()
           if(impose_eos) call impose_equation_of_state()
#endif

        endif
        !
        if(los_number .le. NumPEs-1)then
          !
          !  Scale to zcurrent:
          acurrent = ExpansionFactor
          zcurrent = Redshift
          !
          ! Compute Hubble parameter (in km/s/Mpc).
          CurrentHubbleCt = 100. * HubbleParam *  &
            sqrt(1. + omega0*(1./acurrent-1.) + OmegaLambda* &
            (acurrent**2-1.)) /acurrent
          boxkms = BoxSize / HubbleParam * acurrent * CurrentHubbleCt
#ifdef HUBBLE
          write (*,*) ' changing hubble'
          boxkms = boxkms * HUBBLE
#endif
          !
          ! Number of pixels in this spectrum
          nveloc = int(boxkms / vpixsizekms) + 1
          !
        endif
        !
        ! project SPH data to sightline
        call cpu_timer_start(doprojecteach)
        call projectdata()
        call cpu_timer_stop(doprojecteach)
        ncontribute(los_number+1) = ncontr
        !
        ! compute individual spectrum for each ion
        call cpu_timer_start(domakespectraeach)
        call makespectra()
        call cpu_timer_stop(domakespectraeach)
        !
        ! convolve with instrumental profile
        if (do_convolve_spectrum) call convolve_short_spectrum()
        !
        ! add noise
        if(generate_noise) call add_noise()
        !
#ifdef MPI
        call mpi_barrier(mpi_comm_world, ierr)
#endif  
        !          
!!$        if(verbose .or. MyPE == 0) &
        if (MyPE == 0) then
             write (*,'('' MyPE = '',I2,'' file is '',i2,'' out of '',i3,'' sightline = '',i5,'' out of '',i6,'' rhocb= '',e12.4)') & 
             MyPE, ifile, nlosfiles, los_number+1, numspec, rhocb
             if(verbose) &
                write (*,'('' projectdata time = '',e12.4,'' s, makespectra time = '',e12.4,'' s'')') & 
                cputime(doprojecteach), cputime(domakespectraeach)
        endif
        !
        ! Output Spectra
        do iproc=0, numpes-1
           SpectrumFile = trim(outputdir)//'/spec.'//simfile(ifile)
           if(iproc == MyPE) then
              if(MyPE == 0 .and. los_number == 0) &
                   call create_spectrum_file(ifile)
              call write_short_spectrum(simfile(ifile), los_number, numspec)
           endif
           !
#ifdef MPI
           call mpi_barrier(mpi_comm_world, ierr)
#endif  
           !
        enddo
        !
      enddo los_loop
#ifdef MPI
      call mpi_barrier(mpi_comm_world, ierr)
#endif
      call update_short_spectrum(simfile(ifile), numspec)
    enddo files_loop
    !
  endif docontspec
  !
  call endrun(' done ')
  !
end program specwizard


subroutine endrun(message) !Shut down cleanly
  use my_mpi
  implicit none
  !
  character(*), intent(in) :: message
  integer :: ierr
  !
#ifdef MPI
  call mpi_barrier(mpi_comm_world, ierr)
  call MPI_finalize(ierr)
#endif
  !
  if(MyPE == 0) &
       write (*,*) message
  stop
  !
end subroutine endrun


subroutine abortrun(message) !Crash
  use my_mpi
  implicit none
  !
  character(*), intent(in) :: message
  integer :: ierr
  !
  write (*,*) message
  write (*,*) ' MPI abort called by processor= ',MyPE
  !
#ifdef MPI
  call mpi_abort(mpi_comm_world,324,ierr)
  call MPI_finalize(ierr)
#endif
  !
  stop
  !
end subroutine abortrun

