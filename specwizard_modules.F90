! Users should normally only modify variables in modules
! runtime: -- names of files and directories
! spectra: -- set properties of spectra to generate
! imposed_metallicity: if you want to change the element abundances from those in the simulation files

module numbers
  implicit none
  !
  ! set precision of single and double precision real and integers
  integer, parameter :: singleR = selected_real_kind(p=6,r=37)
  integer, parameter :: doubleR = selected_real_kind(p=15,r=307)
  integer, parameter :: singleI = selected_int_kind(9)
  integer, parameter :: doubleI = selected_int_kind(18)
  ! signal invalid values
  real(kind=doubleR), parameter    :: invalid_R = -1.d99
  integer(kind=singleI), parameter :: invalid_I = -1234567890
  character(len=20), parameter     :: invalid   = 'INVALID'
  !
end module numbers


module cpu_timers
  use numbers
  implicit none
  !
  integer,parameter  :: ntimers = 20
  integer, parameter :: main=1, dospectra=2, doinsert=3, dointerpolate=4
  real(kind=doubleR) :: cputime(ntimers), cpustart(ntimers), cpustop(ntimers)
  logical            :: cpu_is_running(ntimers)
  !
end module cpu_timers


module my_mpi
#ifdef MPI
  use mpi
#endif
  implicit none
  !
!#ifdef MPI
!  include   'mpif.h'
!#endif
  integer   :: MyPE, MasterPE, NumPEs
  !
end module my_mpi


module random_numbers
  use numbers
  implicit none
  !
  integer, parameter :: ntypes=6, nran_max=100000
  integer, parameter :: ran_file=1, ran_sight=2, ran_noise=3, ran_metal=4, ran_fill=5, ran_los=6
  integer, parameter :: seed0(ntypes) = (/-13,-14,-15,-16,-17,-18/)
  integer, parameter :: random_type(ntypes) = (/0,0,1,1,0,0/)  ! 0=uniform [0,1], 1 = Gaussian deviate
  integer            :: current(ntypes), cycles(ntypes)=0, seed(ntypes)
  real(kind=doubleR) :: randoms(ntypes,nran_max)
  !
end module random_numbers


module atomic_data
  use numbers
  implicit none
  !
  !    ------------	
  !    Atomic data:
  !    ------------	
  !    Atomic mass source: Atomic Weights of the Elements 2005, Pure Appl. Chem, 78, 11
  !    Atomic mass unit source: http://physics.nist.gov/cuu/Constants/index.html
  !    Rest wavelength source: Morton 2003, ApJS, 149, 205
  !    Oscillator strength source: Morton 2003, ApJS, 149, 205 (Most UV lines)
  !                                Verner et al 1994, AAPS, 108, 287
  !                                Verner et al 1996, Atomic Data Nucl. Data Tables, 64, 1 (O8 doublet)
  !                                List from Jelle Kaastra (current SPEX on 2018-02-06), X-ray lines
  !
  !    Multiplet transitions should be in order of decreasing osc strength!
  !
  !Physical constants for strength of 21cm emission:
  real(kind=doubleR), parameter :: A21=2.85e-15  !inverse lifetime of 21 cm level [s^{-1}]
  real(kind=doubleR), parameter :: nu21=1.4204e9 !21 cm frequency [s^{-1}]
  !
  real(kind=doubleR), parameter :: atom_munit = 1.66053886e-24 ! g
  !
  ! parameters for Hydrogen atom
  !  real(kind=doubleR), parameter :: massH = 1.67262171e-24 ! g (proton mass)
  real(kind=doubleR), parameter :: massH = 1.00794 * atom_munit ! g
  integer(kind=singleI), parameter ::  nlyman_all = 31
  real(kind=doubleR) ::  Lambda_H1(nlyman_all), f_H1(nlyman_all)
  data Lambda_H1 / 1215.6701, 1025.7223, 972.5368, 949.7431, &
       937.8035, 930.7483, 926.2257, 923.1504, &
       920.9631, 919.3514, 918.1294, 917.1806, &
       916.429, 915.824, 915.329, 914.919, &
       914.576, 914.286, 914.039, 913.826, &
       913.641, 913.480, 913.339, 913.215, &
       913.104, 913.006, 912.918, 912.839, &
       912.768, 912.703, 912.645 /
  data f_H1 / 0.416400, 0.079120, 0.029000, 0.013940, &
       0.007799, 0.004814, 0.003183, 0.002216, &
       0.001605, 0.00120,  0.000921, 7.226e-4, &
       0.000577, 0.000469, 0.000386, 0.000321, &
       0.000270, 0.000230, 0.000197, 0.000170, &
       0.000148, 0.000129, 0.000114, 0.000101, &
       0.000089, 0.000080, 0.000071, 0.000064, &
       0.000058, 0.000053, 0.000048 /
  !
  ! He
  real(kind=doubleR), parameter :: massHe = 4.002602 * atom_munit
  real(kind=doubleR) ::  Lambda_He2(1)
  data Lambda_He2 / 303.7822 /
  real(kind=doubleR) :: f_He2(1)
  data f_He2 / 0.416 /
  !
  ! 21cm
  real(kind=doubleR) ::  Lambda_21cm(1)
  data Lambda_21cm / 2.1e9 /
  real(kind=doubleR) :: f_21cm(1)
  data f_21cm / 0.416 /
  !
  ! C
  real(kind=doubleR), parameter :: massC  = 12.0107 * atom_munit
  real(kind=doubleR) :: Lambda_C2(2), f_C2(2)
  data Lambda_C2 / 1334.5323, 1036.3367 /
  data f_C2 / 0.127800, 0.118000 /
  real(kind=doubleR) :: Lambda_C3(1), f_C3(1)
  data Lambda_C3 / 977.0201 /      
  data f_C3 / 0.7570 /
  real(kind=doubleR) :: Lambda_C4(2), f_C4(2)
  data Lambda_C4 / 1548.2041, 1550.7812 /      
  data f_C4 / 0.189900, 0.094750 /
  real(kind=doubleR) :: Lambda_C5(2), f_C5(2)
  data Lambda_C5 / 40.2678, 34.9728 /      
  data f_C5 / 0.648, 0.141 /
  real(kind=doubleR) :: Lambda_C6(4), f_C6(4)
  data Lambda_C6 / 33.7342, 33.7396, 28.4652, 28.4663 /      
  data f_C6 / 0.277, 0.139, 0.0527, 0.0263 /
  !
  ! N  
  real(kind=doubleR), parameter :: massN  = 14.0067 * atom_munit
  real(kind=doubleR)::  Lambda_N2(1), f_N2(1)
  data Lambda_N2 / 1083.9937 /      
  data f_N2 / 0.111 /
  real(kind=doubleR) :: Lambda_N3(1), f_N3(1)
  data Lambda_N3 / 989.799 /      
  data f_N3 / 0.12287 /
  real(kind=doubleR) :: Lambda_N4(1), f_N4(1)
  data Lambda_N4 / 765.148 /      
  data f_N4 / 0.632 /
  real(kind=doubleR) :: Lambda_N5(2), f_N5(2)
  data Lambda_N5 / 1238.821, 1242.804 /      
  data f_N5 / 0.156000, 0.0770 /
  real(kind=doubleR) :: Lambda_N6(1), f_N6(1)
  data Lambda_N6 / 28.7875 /      
  data f_N6 / 0.675 /
  real(kind=doubleR) :: Lambda_N7(2), f_N7(2)
  data Lambda_N7 / 24.7792, 24.7846 /      
  data f_N7 / 0.277, 0.139 /
  !
  ! O  
  real(kind=doubleR), parameter :: massO  = 15.9994 * atom_munit
  real(kind=doubleR) :: Lambda_O1(3), f_O1(3)
  data Lambda_O1 / 1302.1685, 988.7734, 971.7382 /
  data f_O1 / 0.048000, 0.0465, 0.0116 /
  real(kind=doubleR) :: Lambda_O3(2), f_O3(2)
  data Lambda_O3 / 702.332, 832.927 /
  data f_O3 / 0.126, 0.0998 /
  real(kind=doubleR) :: Lambda_O4(1), f_O4(1)
  data Lambda_O4 / 787.711 /
  data f_O4 / 0.110 /
  real(kind=doubleR) :: Lambda_O5(1), f_O5(1)
  data Lambda_O5 / 629.730 /
  data f_O5 / 0.499 /
  real(kind=doubleR) :: Lambda_O6(2), f_O6(2)
  data Lambda_O6 / 1031.9261, 1037.6167 /      
  data f_O6 / 0.13250, 0.06580 /
  real(kind=doubleR) :: Lambda_O7(2), f_O7(2)
  data Lambda_O7 / 21.6019, 18.6284 /      
  data f_O7 / 0.696, 0.146 /
  real(kind=doubleR) :: Lambda_O8(4), f_O8(4)
  data Lambda_O8 / 18.9671, 18.9725, 16.0055, 16.0067/      
  data f_O8 / 0.277, 0.139, 0.0527, 0.0263/
  !
  ! Ne
  real(kind=doubleR), parameter :: massNe = 20.1797 * atom_munit
  real(kind=doubleR) :: Lambda_Ne8(2), f_Ne8(2)
  data Lambda_Ne8 / 770.409, 780.324 /
  data f_Ne8 / 0.103, 0.0505 /
  real(kind=doubleR) :: Lambda_Ne9(1), f_Ne9(1)
  data Lambda_Ne9 / 13.4471 /
  data f_Ne9 / 0.724 /
  !
  ! Mg
  real(kind=doubleR), parameter :: massMg = 24.3050 * atom_munit
  real(kind=doubleR)::  Lambda_Mg2(2), f_Mg2(2)
  data Lambda_Mg2 / 2796.3543, 2803.5315 /
  data f_Mg2 / 0.6155, 0.3058 /
  !
  ! Al
  real(kind=doubleR), parameter :: massAl = 26.981538 * atom_munit
  real(kind=doubleR)::  Lambda_Al2(1), f_Al2(1)
  data Lambda_Al2 / 1670.79 /
  data f_Al2 / 1.73812 /
  real(kind=doubleR)::  Lambda_Al3(2), f_Al3(2)
  data Lambda_Al3 / 1854.72, 1862.79 /
  data f_Al3 / 0.559399, 0.277866 /
  !
  ! Si
  real(kind=doubleR), parameter :: massSi = 28.0855 * atom_munit
  real(kind=doubleR) :: Lambda_Si2(1), f_Si2(1)
  data Lambda_Si2 / 1260.420 /      
  data f_Si2 / 1.17621 /
  real(kind=doubleR) :: Lambda_Si3(1), f_Si3(1)
  data Lambda_Si3 / 1206.500 /      
  data f_Si3 / 1.63 /
  real(kind=doubleR) :: Lambda_Si4(2), f_Si4(2)
  data Lambda_Si4 / 1393.76018, 1402.77291 /      
  data f_Si4 / 0.513, 0.254 /
  !
  ! S
  real(kind=doubleR), parameter :: massS = 32.065 * atom_munit
  real(kind=doubleR) :: Lambda_S5(1), f_S5(1)
  data Lambda_S5 / 786.48 /      
  data f_S5 / 1.46 /
  !
  ! Fe  
  real(kind=doubleR), parameter :: massFe = 55.8452 * atom_munit
  real(kind=doubleR) ::  Lambda_Fe2(9), f_Fe2(9)
  data Lambda_Fe2 / 1144.9379, 1608.45085, 1063.1764, 1096.8769,  &
    1260.533, 1121.9748, 1081.8748, 1143.2260,  &
    1125.4477 /      
  data f_Fe2 / 0.083, 0.0577, 0.0547, 0.032700,  &
    0.024000, 0.0290, 0.012600, 0.0192, &
    0.0156 /
  real(kind=doubleR) ::  Lambda_Fe3(1), f_Fe3(1)
  data Lambda_Fe3 / 1122.52 /      
  data f_Fe3 / 0.0544257 /
  real(kind=doubleR) ::  Lambda_Fe17(2), f_Fe17(2)
  data Lambda_Fe17 / 15.0140, 15.2610 /      
  data f_Fe17 / 2.72, 0.614 /
  real(kind=doubleR) ::  Lambda_Fe19(2), f_Fe19(2)
  data Lambda_Fe19 / 13.5180, 13.5146 /      
  data f_Fe19 / 0.717, 0.0199 /
  real(kind=doubleR) ::  Lambda_Fe21(1), f_Fe21(1)
  data Lambda_Fe21 / 12.2840 /      
  data f_Fe21 / 1.24 /
  !
end module atomic_data


module physical_constants
  use numbers
  implicit none
  !
  !    ----------------------------
  !    Natural constants and units:
  !    ----------------------------
  !    Source: http://physics.nist.gov/cuu/Constants/index.html
  real(kind=doubleR), parameter :: Boltz = 1.3806505e-16, &! erg/K 
    Planck = 6.6260693e-27, & ! erg s
    LightSpeed = 2.99792458e10, &! cm/s
    ThomsonCross = 6.65245873e-25, & ! cm^2
    G = 6.6742e-8 ! cm^3/g/s^2
  !   Source: some textbook.
  real(kind=doubleR), parameter :: Mpc = 3.085678e24, & ! cm
    yr = 3.1558e7, & ! s
    Msun = 1.989e33,  & ! g
    H0 = 1e7/Mpc, &! s^{-1} (100 km/s/Mpc)
    lyalpha = 1215.6701 !   Rest wavelength source: Morton 2003, ApJS, 149, 205
  !    --------------------
  !    Numerical constants:
  !    --------------------
  real(kind=doubleR), parameter :: tpi = 6.283185307179d0, pi = tpi/2.d0
  !
end module physical_constants


module solar_data
  use numbers
  use atomic_data
  implicit none
  !
  ! Assumed values for the Sun
  ! Solar abundances source: Hazy I, CLOUDY 94 (Anders & Grevesse 1989; Grevesse & Noels 1993)
  !
  ! Total metal mass fraction
  real(kind=doubleR), parameter :: Zmass_solar   = 0.0126637   ! M_Metal/M_tot
  real(kind=doubleR), parameter :: Ymass_solar   = 0.2466      ! M_Helium/M_total
  real(kind=doubleR), parameter :: Xmass_solar   = 1.d0-Ymass_solar-Zmass_solar ! M_Hydrogen/M_total
  real(kind=doubleR), parameter :: YNumber_solar = (Ymass_solar/Xmass_solar)*(massH/massHe) ! N_Helium/N_Hydrogen
  !
  ! Metallicities in number relative to hydrogen.	
  real(kind=doubleR),parameter  ::   ZC_solar = 3.55e-4
  real(kind=doubleR),parameter  ::   ZN_solar = 9.33e-5
  real(kind=doubleR),parameter  ::   ZO_solar = 7.41e-4
  real(kind=doubleR),parameter  ::  ZNe_solar = 1.23e-4
  real(kind=doubleR),parameter  ::  ZMg_solar = 3.80e-5
  real(kind=doubleR),parameter  ::  ZAl_solar = 2.95e-6
  real(kind=doubleR),parameter  ::  ZSi_solar = 3.55e-5
  real(kind=doubleR),parameter  ::   ZS_solar = 1.62e-5
  real(kind=doubleR),parameter  ::  ZFe_solar = 3.24e-5
  !
end module solar_data


module primordial_BBN
  use numbers
  implicit none
  !
  real(kind=doubleR), parameter :: Ymass = 0.24 ! Helium mass fraction
  real(kind=doubleR), parameter :: yNumber = Ymass / (1.-Ymass) / 4.d0 ! n_He / n_H for primordial gas
  !
end module primordial_BBN


module runtime
  use numbers
  implicit none
  !
  ! Variables here control the run. They are (mostly) read from the parameter file.
  ! Required variables are initialised to invalid
  ! Optional variables are initialized to default values
  !
  character(len=120) :: parameter_file = 'specwizard.par'  ! The parameter filename, read from command line
  character(len=300) :: datadir = invalid             ! los file directory, or data directory containing snapshots
#ifdef EAGLE
  character(len=300) :: snap_base                     ! base-name -f eagle snapshot
#endif
  character(len=300) :: outputdir = invalid           ! file output directory
  character(len=300) :: ibdir = invalid               ! directory containing ionization tables
  character(len=300) :: file_list = invalid           ! ASCII file contianing list of LOS files in directory datadir to use
  logical  :: do_long_spectrum = .true.               ! do_long_spectrum = true: combine many individual sightlines into several long spectra ...
  character(len=200) :: SpectrumFile = invalid        ! output file name, overrides automatic generation
  logical  :: use_snapshot_file = .false.             ! Generate spectra from snapshots rather than LOS files?
  integer(kind=singleI) :: snap = invalid_I           ! If working from a snapshot, then need to know the snapshot number
  logical  :: use_random_los = .true.                 ! If working from a snapshot, choose LOS randomly?
  character(len=200) :: los_coordinates_file='projections.dat' ! Otherwise need a projections file
  real(kind=doubleR), allocatable :: x_fraction_array(:), &    ! If working from a snapshot we'll need to store ...
    y_fraction_array(:), z_fraction_array(:), phi_array(:), theta_array(:)          ! ... all these LOS coordinates
  integer(kind=singleI), allocatable :: ncontribute(:), ncontribute_global(:)
  integer(kind=singleI) :: number_of_LOS = 0          ! ... and know how many los are given
  logical  :: verbose = .false.                       ! makes specwizard more wordy
  logical  :: overwrite = .false.                     ! if true, overwrite output file.
  logical  :: use_maxdens_above_zmax = .false.        ! use maximum density above the maximum redshift in the ionization tables
  logical  :: gimic = .false.                         ! gimic or owls los files?
  logical  :: urchin = .false.                        ! if true, read and use neutral fractions from post-processed snapshot files in urchin_dir
  logical  :: do_periodic = .true.                    ! if true, implements periodic boundary conditions
  character(len=300) :: urchindir = invalid           ! urchin file directory containing files for this snapshot
  logical  :: use_urchin_temperature = .false.        ! if true, read particle temperature result from urchin rather than snapshot
  logical  :: use_smoothed_abundance = .false.        ! if true, read the smoothed element abundance
  logical  :: subtract_Hmol = .false.                 ! if urchin, if T: X_Si2 = X_H1  else: X_Si2 = X_H1 + X_Hmol
  logical  :: use_gaussian_kernel = .false.           ! use a normalized, truncated gaussian kernel rather than M4
  logical  :: integrate_kernel = .false.              ! integrate kernel exactly over whole pixel, rather than using r-pixel value
  logical  :: wmap7 = .true.                         ! if wmap7, read attributes as present in the OWLS wmap7 run
  logical  :: read_ionbal_from_single_file = .false.  ! read ionization balance from TT's 'single file' format?
  logical  :: output_frequency = .false.              ! Output frequencies instead of wavelengths
  logical  :: NoPecVel=.false.                        ! if NoPecVel, then do not include peculiar velocities
  logical  :: ionfracone = .false.                    ! if true, then assume all ion fractions are 1.
  logical, parameter :: docycling=.false.             ! cycle short spectra to have minimum optical depth at ends, when generating long spectra ...
                                                      ! ... (usually a bad idea)
  logical :: ignore_starforming = .false.             ! if true, sets mass of star forming particles to 0.0
  logical :: setmaxt4sfgas = .true.                   ! if true, then the temperatures of star-forming SPH particles are set to 10^4 K
  logical :: impose_eos =.false.                      !impose EOS on IGM particles? ...
  logical :: add_turbulence =.false.                  ! add a turbulent component to the line widths
  real(kind=doubleR) :: imposed_eos_T0 = invalid_R, & ! ... if so, we'll need these
    imposed_eos_gamma = invalid_R, &                  ! ... if so, we'll need these
    imposed_eos_maxod = invalid_R                     ! ... if so, we'll need these
  logical  :: output_realspacenionweighted_values = .false. ! if true, output optical depths for each ion
  logical  :: output_realspacemassweighted_values = .false. ! if true, output density and temperature weighted by ion
  logical  :: output_zspaceopticaldepthweighted_values = .false. ! if true, output density and temperature weighted by optical depth
  logical  :: do_convolve_spectrum = .true.           ! if true, convolve final spectrum with instrumental profile of 'fwhm' km/s
  real(kind=doubleR), parameter :: &                  ! small values
    small_rho = 1.d-30, small_temp=0.1, &             ! small values
    small_metallicity = 1.d-10, small_column  = 1     ! small values
  !
end module runtime


module spectra
  use numbers
  use atomic_data, only: nlyman_all
  implicit none
  !
  ! ++++++++++++++++++++++++++++ do_long_spectrum = .true. +++++++++++++++++++++++++++++++++++
  ! Generate spectra over a large wavelength range, including transitions from many elements
  ! Spectra cover the wavelength range  [minlambda,maxlambda] and include
  ! absorption from the redshift range [zabsmin,zabsmax] or [zabsmin,zqso] if zqso<zabsmax
  ! At a given redshift z, specwizard will randomly choose a sight line file within +/- fzresol
  ! 
  ! redshift of QSO, min and max redshift of spectra to take into account
  !
  real(kind=doubleR) :: zqso = invalid_R, zabsmin = invalid_R, zabsmax = invalid_R
  real(kind=doubleR) :: fzresol = invalid_R                           ! slack in redshift for short spectra to be included in continuous spectrum
  real(kind=doubleR) :: minlambda = invalid_R, maxlambda  = invalid_R ! Angstrom
  integer(kind=singleI) :: nspec = invalid_I                          ! number of continuous spectra to generate
  !
  ! ++++++++++++++++++++++++++++ do_long_spectrum = .true. +++++++++++++++++++++++++++++++++++
  !
  ! ++++++++++++++++++++++++++++ spectral parameters +++++++++++++++++++++++++++++++++++++++++
  !
  ! Elements to be included:
  logical :: doall = .false.
  logical :: doH1= .true., &
    doHe2=.false., &
    doC2=.false., doC3=.false., doC4=.false., doC5=.false., doC6=.false., &
    doN2= .false., doN3=.false., doN4=.false., doN5=.false., doN6=.false., doN7=.false., &
    doO1=.false., doO3=.false., doO4=.false., doO5=.false., doO6=.false., doO7 = .false., doO8 = .false., &
    doNe8=.false., doNe9=.false., &
    doMg2=.false., &
    doAl2=.false., doAl3=.false., &
    doSi2=.false., doSi3=.false., doSi4 = .false., &
    doS5=.false., &
    doFe2=.false., doFe3=.false., doFe17=.false., doFe19=.false., doFe21=.false., &
    do21cm=.false.
  !
  ! number of lines to be included in Lyman-series, <= nlyman_all=31
  integer(kind=singleI) :: nLyman=invalid_I
  !
  ! amplitude of rescaling of ionising background
  real(kind=doubleR) :: ibfactor = invalid_R
  logical            :: use_fitted_ibfactor = .false.
  logical            :: ibfactor_he_reionization = .false.
  !
  ! integration of the spectra
  logical            :: limsigma = .true.! Only follow profiles until negligible?
  !
  real(kind=doubleR) :: pixsize = invalid_R, fwhm = invalid_R ! Pix size in A, FWHM in km/s
  !
  !     Pixel size (km/s) of long spectrum before instrumental processing
  !     (should be small compared to expected b-parameters, unless the
  !     thermal profiles are integrated exactly).
  !     Also pixel size of individual small spectra
  real(kind=doubleR) ::  vpixsizekms = invalid_R
  !
  ! ************************ variables needed for computing the spectra: do not modify
  !
  integer(kind=singleI)              :: ispec             ! current spectrum we are doing
  integer                            :: nion
  character(len=10), allocatable     :: ions(:)
  integer(kind=singleI), allocatable :: ion_elnr(:), nlines(:)
  real(kind=doubleR), allocatable    :: ion_mass(:), ion_z_abs(:), ionfrac(:), elementfrac(:), totnr_ion(:), taumax(:)
  integer(kind=singleI) :: ion_H=-1
  ! Max number transitions per element
  integer(kind=singleI), parameter   :: n_lines_max = max(13, nLyman_all)
  real(kind=doubleR), allocatable    :: lambda_rest(:,:), & ! rest wavelength of transition
    fosc(:,:) ! oscillator strength
  !
  real(kind=doubleR)  :: zstart,zend ! actual redshift range
  !
  !      Interpolate thermal profiles to pixel or integrate over pixel? 
  !      Integration is more accurate. 
  !      Interpolating is faster for high-resolution spectra. 
  !      There are enough built in checks to ensure that interpolation is
  !      sufficiently accurate. However, for low-resolution spectra
  !     integration may be faster. 
  logical            :: integrate_thermprof_exactly = .true.
  !
  !   Min. max. optical depth for inclusion:
  real(kind=doubleR) :: &
       minbother_blue = invalid_R, & ! lambda_0 =< Ly-alpha
       minbother_red  = invalid_R   ! lambda_0 >  Ly-alpha
  !
  !Number of star formaing particles on line of sight
  integer(kind=singleI)    :: nsf
  !
  !Number of significantly neutral particles on line of sight
  integer(kind=singleI)    :: n_neut
  !
  ! properties of simulation files
  integer(kind=singleI)  :: nlosfiles=0, nfull=0, ntotfiles=0
  real(kind=doubleR)     :: simzmin=1.e12, simzmax=1.d0
  real(kind=doubleR), allocatable    :: simz(:)
  character(len=80), allocatable     :: simfile(:)
  integer(kind=singleI), allocatable :: nlos_in_file(:), ichoose(:)
  ! info about simulation files used
  integer(kind=singleI)              :: nsimfile_used
  integer(kind=singleI), parameter   :: max_nsimfile_used = 1000
  integer(kind=singleI), allocatable :: los_used(:)
  character(len=120), allocatable    :: losfile_used(:)
  real(kind=doubleR), allocatable    :: x_physical_used(:), y_physical_used(:), ibfactor_used(:), icshift_used(:)
  !
  ! spectrum array full spectrum
  real(kind=doubleR), allocatable    :: lambda(:), voverc(:), tau_long(:,:), flux(:), tau_long_strongest(:,:)   ! size nvpix 
  !
  ! column density array
  real(kind=doubleR), allocatable    :: cdens_ion_integrated(:)
  !
  ! redshift-space quantities
  real(kind=doubleR), allocatable    :: temp_z_ion_long(:,:), rho_z_ion_long(:,:) ! size nvpix
  ! real(kind=doubleR), allocatable    :: temp_z_long(:), rho_z_long(:)
  !
  ! real-space quantities
  real(kind=doubleR), allocatable    :: temp_ion_long(:,:), n_ion_long(:,:), rho_ion_long(:,:) ! size nvpix
  real(kind=doubleR), allocatable    :: temp_long(:), rho_long(:), met_long(:)
  !
  real(kind=doubleR), allocatable    :: flux_convolved(:)  ! size 2*nvpix 
  complex(kind=doubleR), allocatable :: fft(:)
  !
  ! rebinned spectrum
  integer(kind=singleI) :: n_binned_flux
  real(kind=doubleR), allocatable    :: binned_lambda(:), binned_flux(:)! size n_binned_spectrum
  real(kind=doubleR), allocatable    :: binned_noise_sigma(:) , binned_noise_random(:)
  real(kind=doubleR), allocatable    :: binned_temp_z_ion(:,:), binned_rho_z_ion(:,:)
  real(kind=doubleR), allocatable    :: binned_temp_ion(:,:),binned_n_ion(:,:),binned_rho_ion(:,:)
  real(kind=doubleR), allocatable    :: binned_tau_ion(:,:), binned_tau_ion_strongest(:,:)
  integer(kind=singleI), allocatable :: binned_spectrum_boundary(:)
  !
  ! Computed values
  integer(kind=singleI) :: nvpix
  real(kind=doubleR)    :: vpixsize
  !
  integer(kind=singleI) :: icshift
  !
  ! spectrum arrays single spectrum
  real(kind=doubleR)    :: zcurrent, acurrent
  integer(kind=singleI) :: nveloc
  real(kind=doubleR), allocatable  :: rho_tot(:),temp_tot(:),met_tot(:),veloc_tot(:)
  !
  ! real-space quantities weigthed by number of ions
  real(kind=doubleR), allocatable  :: n_ion(:,:), temp_ion(:,:), veloc_ion(:,:), od_ion(:,:), rho_ion(:,:)
  !
  ! redshift-space quantities weigthed by number of ions
  real(kind=doubleR), allocatable  :: rho_z_ion(:,:), temp_z_ion(:,:), tau_ion(:,:), &
    veloc_z_ion(:,:), nion_z_ion(:,:), flux_ion(:,:)
  integer(kind=singleI), allocatable :: spectrum_boundary(:)
  !
  real(kind=doubleR), allocatable  :: turbulence(:)
  !
  ! temporary arrys
  real(kind=doubleR), allocatable  :: normw(:)
  real(kind=doubleR), allocatable  :: tau(:), cdens(:), &
    vhubble(:), nionw(:),rhow(:),tempw(:), velocw(:), work(:), work2(:), voc(:), vocsim(:)
  !
end module spectra


module noisedata
  use numbers
  implicit none
  !
  ! noise properties of spectra
  logical :: generate_noise = .true. 
  !.false.
  logical :: use_noise_file = .true.
  !   If either sigtonoise =< 0 or minnoise < 0, then sigma(lambda,flux) 
  !   is read in from the file "noisefile", defined below.
  !   (sigma is the standard deviation of the Gaussian random variable)
  !   If sigtonoise > 0 and minnoise >= 0, then
  !   sigma = minnoise + (1/sigtonoise - minnoise) * flux
  !   Note that we must have minnoise < 1 / sigtonoise
  real(kind=doubleR) :: sigtonoise = invalid_r, &
       minnoise = invalid_r
  !   Noise interpolation table:
  character(len=120) :: noisefile = invalid
  !
  integer(kind=singleI) :: n_nl, n_nf
  real(kind=singleR), allocatable :: n_sigma(:,:), n_lambda(:), n_flux(:)
  !
end module noisedata


module modified_metallicity
  use numbers
  implicit none
  !
  logical :: modify_metallicity = .false. ! change abundance from the one read from the snapshots
  !
  ! Variables for changing the metallicity only of flagged particles
  character(len=120) :: particle_file_name = invalid
  logical :: read_part_ids_from_file = .false.
  real(kind=doubleR) :: flagged_particle_metallicity = invalid_R
  !
  ! if true, scale abundances from simulation with factor z_rel
  ! if false, use z_rel as new metallicity (relative to solar)
  logical :: scale_simulation_abundances = .false.
  !
  ! z_rel is imposed metallicity scaling factor, relative to solar
  real(kind=doubleR) :: z_rel = invalid_R
  !
  ! Use these abundances relative to solar.
  !  .... the assumed *solar* abundances are in module atomic_data
  real(kind=doubleR) :: ZC_rel = invalid_R, ZN_rel = invalid_R, ZO_rel = invalid_R,  &
    ZNe_rel = invalid_R, ZMg_rel = invalid_R, ZAl_rel = invalid_R, ZSi_rel = invalid_R, &
    ZS_rel = invalid_R, ZFe_rel = invalid_R
  !
  ! impose metallicity z = z_mean * (density/mean density)^z_index, up to maximum metallicity maxz_rel=10.0
  logical               :: impose_z_rho_relation = .false.
  real(kind=doubleR)    :: z_index = invalid_R
  real(kind=doubleR)    :: z_mean = invalid_R
  real(kind=doubleR)    :: maxz_rel = invalid_R
  !
  ! impose log-normal metallicity
  logical               :: log_normal_scatter = .false.
  integer(kind=singleI) :: z_sig_bin = invalid_I     ! divide computational volume in (z_sig_bin)^3 cells, and impose
                                                     ! lognormal metallicity with
                                                     ! z_sig_dex  scatter for all metals
  real(kind=doubleR)    :: z_sig_dex = invalid_R
  !
!!!  ! parameters of volume filling factor of metals
!!!  logical :: impose_volume_filling_factor = .false.  !
!!!  real(kind=doubleR) :: volume_filling_factor = 1.0 ! impose metallicity in volume fraction zfilfator<1. Zero metals outside this range
  !
end module modified_metallicity

module RegionExtent
  use numbers
  implicit none
  real(kind=singleR) :: RegionExtentX(2), RegionExtentY(2), RegionExtentZ(2)
end module RegionExtent

module particledata
  use numbers
  implicit none
  !
  ! header values
  real(kind=doubleR)    :: rhocb, BoxPhys, proton_mass
  real(kind=doubleR)    :: Boxkms
  !
  integer(kind=singleI) :: NGas
#ifdef EAGLE
  integer(kind=doubleI), allocatable :: PartID(:)
#else
  integer(kind=singleI), allocatable :: PartID(:)
#endif
  !
#ifdef EAGLE
  real(kind=doubleR), allocatable :: Position(:,:)
  real(kind=doubleR), allocatable :: Velocity(:,:)!!!, ShiftedPosition(:,:), ShiftedVelocity(:,:)
  real(kind=doubleR), allocatable:: Mass(:),ParticleDensity(:),&
    ParticleSmoothingLength(:),ParticleTemperature(:), Metallicity(:), & 
    MassFractions(:,:),  &
    MetallicityInSolar(:), Zmetal(:), Zrelat(:), StarFormationRate(:),&
    ParticleNeutralHFraction(:),ParticleMolecularHFraction(:)
#else
  real(kind=singleR), allocatable :: Position(:,:)
  real(kind=singleR), allocatable :: Velocity(:,:)!!!, ShiftedPosition(:,:), ShiftedVelocity(:,:)
  real(kind=singleR), allocatable:: Mass(:),ParticleDensity(:),&
    ParticleSmoothingLength(:),ParticleTemperature(:),Metallicity(:), &
    MassFractions(:,:),  &
    MetallicityInSolar(:), Zmetal(:), Zrelat(:), StarFormationRate(:),&
    ParticleNeutralHFraction(:),ParticleMolecularHFraction(:)
#endif
  integer(kind=singleI), allocatable ::  Boundary(:)
  integer(kind=singleI), allocatable :: itype(:)
  !
  ! Indices corresponding to given elements in zmass array (order is arbitrary)
  integer :: nspecies = 0
  integer :: nspecies_max = 11
  integer(kind=singleI) :: H_index=-1, He_index=-1, C_index=-1, N_index=-1, &
    O_index=-1, Ne_index=-1, Mg_index=-1, Al_index=-1, Si_index=-1, S_index=-1, &
    Fe_index=-1
  logical :: requireH=.true., requireHe=.false., requireC=.false., requireN=.false., &
    requireO=.false., requireNe=.false., requireMg=.false., requireAl=.false., &
    requireSi=.false., requireS=.false., requireFe=.false.
    ! default: only elements of which we want ions are needed, except hydrogen, 
    !which is used for the log H number density in het ion balance tables
  real(kind=doubleR),allocatable :: ElementAtomicMass(:)
  !
  ! sightline projection parameters
  real(kind=singleR)    :: x_fraction, y_fraction, z_fraction
  real(kind=singleR)    :: x_physical, y_physical, z_physical
  integer(kind=singleI) :: ncontr
  !
  ! scaling of density from simulation output time to spectrum time
  real(kind=doubleR) :: densscale = 1.0
  !
end module particledata


module projection_parameters
  use numbers
  implicit none
  !
  integer(kind=singleI) :: x_axis, y_axis, z_axis
  real(kind=doubleR)    :: x_comoving, y_comoving, z_comoving, phi_projection, theta_projection
  !
end module projection_parameters


module ionization_tables
  use numbers
  implicit none
  !
  integer(kind=singleI), save :: nz, nd, nt, nznt
  real(kind=singleR), allocatable :: ionizbal(:,:) ! size nion, ib_nz times ib_nt times in_nd
  real(kind=singleR), allocatable :: ib_redshift(:), ib_logt(:), ib_logd(:), ionbal(:,:,:)
  real(kind=singleR), allocatable :: ion_tmp(:,:)
  character(len=80), parameter :: &
    LogDensity = 'logd',LogTemperature='logt',RedshiftLabel='redshift',IonizationFractions='ionbal'
  !
end module ionization_tables


module header
  use numbers
  implicit none
  !
  integer(kind=doubleI) :: NumPart_ThisFile(0:5), NumPart_Total(0:5), &
    NumPart_Total_HighWord(0:5), NumFilesPerSnapshot, Flag_SFR, &
    Flag_Cooling, Flag_StellarAge, Flag_Metals, Flag_Feedback, los_this_file
  real(kind=doubleR) :: MassTable(6), Time,Redshift, ExpansionFactor, &
       BoxSize, Omega0,OmegaBaryon,OmegaLambda,HubbleParam,rhoc
  !
end module header


module constants
  use numbers
  implicit none
  !
  real(kind=doubleR) :: gamma, gravity, solar_mass, solar_lum, boltzmann, &
    gas_const, cm_per_mpc
  !
end module constants


module units
  use numbers
  implicit none
  !
  real(kind=SingleR) :: UnitLength_in_cm, UnitMass_in_g, &
    UnitVelocity_in_cm_per_s, UnitDensity_in_cgs, UnitEnergy_in_cgs
  !
end module


module parameters
  use numbers
  implicit none
  !
  integer(kind=singleI) :: BG_Nelements, ComovingIntegrationOn, &
    TypeOfTimestepCriterion, TypeOfOpeningCriterion,DesNumNgb, &
    AGB_MassTransferOn, POPIII_MassTransferOn, POPIII_EnergyTransferOn, &
    SNII_WindOn,SNII_WindIsotropicOn,SNIa_MassTransferOn,SNIa_EnergyTransferOn
  real(kind=doubleR) :: InitGasU_ERG,MinGasU_ERG 
  character(len=20), allocatable     :: ElementNames(:)
  character(len=100) :: IMF_model, IMF_LifetimeModel,SNIa_Model
  real(kind=doubleR) :: BufferSize  ! mistake: should be integer
  real(kind=doubleR) :: MaxNumNgbDeviation,InitAbundance_Hydrogen, InitAbundance_Helium, InitAbundance_Carbon, &
    InitAbundance_Nitrogen, InitAbundance_Oxygen, InitAbundance_Neon, InitAbundance_Magnesium, &
    InitAbundance_Silicon, InitAbundance_Iron, CalciumOverSilicon, SulphurOverSilicon, ErrTolIntAccuracy, &
    MaxSizeTimestep, MinSizeTimestep,ErrTolTheta,ErrTolForceAcc,TreeDomainUpdateFrequency,&
    ArtBulkViscConst,CourantFac,PartAllocFactor,TreeAllocFactor, &
    MinGasHsmlFractional,SF_EOSEnergyAtThreshold_ERG, SF_EOSGammaEffective, &
    SF_THRESH_MinPhysDens_HpCM3,SF_THRESH_MinOverDens,SF_THRESH_MetDepExponent,SF_THRESH_MetDepFiducialZ, &
    SF_THRESH_MetDepMaxPhysDens_HpCM3,SF_THRESH_MaxTemp_K,SF_SchmidtLawCoeff_MSUNpYRpKPC2, &
    SF_SchmidtLawExponent,SF_SchmidtLawCoeff_GpSpCM2,IMF_MinMass_MSUN,SNIa_Efficiency_fracwd, &
    SNII_MinMass_MSUN,SNII_MaxMass_MSUN,SNII_Factor_Hydrogen, &
    SNII_Factor_Helium,SNII_Factor_Carbon,SNII_Factor_Nitrogen,SNII_Factor_Oxygen,SNII_Factor_Neon, &
    SNII_Factor_Silicon,SNII_Factor_Iron,POPIII_Energy_ERG,POPIII_NumPerMsun,SNII_WindSpeed_KMpS, &
    SNII_WindMassLoading,SNII_WindDelay_YR
  !
end module parameters

