! Copyleft 2006-2014 NASA and Blue Marble Research (http://www.bluemarble.ch)
! Author(s): Reto Stockli
! 
! This file is part of <phenoanalysis> and was started within the
! NASA Energy and Water Cycle Study (NEWS) grant No. NNG06CG42G. It
! now is a open source software project with code and documentation
! found on: http://phenoanalysis.sourceforge.net
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! IMPORTANT NOTE: Intellectual Property Right and Author statements in the
! some of the code may supersede the above Copyleft / GNU GPL license.

module driver_module

  ! MODULE storing all driver (wrapper) arrays and driver functions

  ! interface to C functions
  interface
     integer(c_int) function read_config (length, string) bind(c)
       use iso_c_binding
       implicit none
       integer(c_int) :: length
       character(c_char), dimension(*) :: string
     end function read_config
  end interface

  ! state/parameter/forcing data vectors
  integer :: nstate  ! number of model states
  integer :: nparam  ! number of model parameters
  integer :: nforce  ! number of model forcing fields
  integer :: nobser  ! number of observation fields

  integer :: nselect_state  ! number of model states for assimilation
  integer :: nselect_param  ! number of model parameters for assimilation
  integer :: nselect_force  ! number of model forcing fields to be perturbed
  integer :: nselect_obser  ! number of observation fields for assimilation

  integer :: noutput_state  ! number of model states for output
  integer :: noutput_param  ! number of model parameters for output 
  integer :: noutput_force  ! number of model forcing fields for output
  integer :: noutput_obser  ! number of observation fields for output

  real(kind=4), allocatable :: state_var(:)  ! prescribed initial variability of states
  real(kind=4), allocatable :: param_var(:,:)! prescribed initial variability of parameters
  real(kind=4), allocatable :: force_var(:)  ! prescribed run-time variability of forcing variables
  real(kind=4), allocatable :: obser_var(:)  ! prescribed run-time variability of observations

  real(kind=4), allocatable :: state(:,:,:,:,:,:)  ! states
  real(kind=4), allocatable :: param(:,:,:)        ! parameters
  real(kind=4), allocatable :: force(:,:,:,:,:)    ! forcing
  real(kind=4), allocatable :: obser(:,:,:,:)      ! observations on output grid

  real(kind=4), allocatable :: state_sim(:,:,:,:,:,:)  ! states (simulator)
  real(kind=4), allocatable :: param_sim(:,:,:)        ! parameters (simulator)
  real(kind=4), allocatable :: force_sim(:,:,:,:,:)    ! forcing (simulator)

  real(kind=4), allocatable :: state_min(:)   ! lower limit for ensemble mean of states
  real(kind=4), allocatable :: param_min(:)   ! lower limit for ensemble mean of parameters
  real(kind=4), allocatable :: force_min(:)   ! lower limit for ensemble mean of forcing
  real(kind=4), allocatable :: obser_min(:)   ! lower limit for ensemble mean of observations

  real(kind=4), allocatable :: state_max(:)   ! upper limit for ensemble mean of states
  real(kind=4), allocatable :: param_max(:)   ! upper limit for ensemble mean of parameters
  real(kind=4), allocatable :: force_max(:)   ! upper limit for ensemble mean of forcing
  real(kind=4), allocatable :: obser_max(:)   ! upper limit for ensemble mean of observations

  integer, allocatable :: select_state(:) ! index of assimilated states
  integer, allocatable :: select_param(:) ! index of assimilated parameters
  integer, allocatable :: select_force(:) ! index of perturbed forcing variables
  integer, allocatable :: select_obser(:) ! index of assimilated observations

  ! state / parameter / forcing / observation indices
  type :: stateidxtype
     integer :: fpar
     integer :: lai
     integer :: gsi
     integer :: temp_fac
     integer :: force_fac
     integer :: chill_fac
     integer :: moist_fac
     integer :: light_fac
     integer :: tmin_ave
     integer :: vpd_ave
     integer :: photo_ave
     integer :: rg_ave
     integer :: rain_ave
     integer :: chill_ave
  end type stateidxtype
  type(stateidxtype) :: stateidx

  type :: paramidxtype
     integer :: tmin_min
     integer :: tmin_max
     integer :: vpd_min
     integer :: vpd_max
     integer :: photo_min
     integer :: photo_max
     integer :: rg_min
     integer :: rg_max
     integer :: rain_min
     integer :: rain_max
     integer :: chill_min
     integer :: chill_max
     integer :: lai_min
     integer :: lai_max
     integer :: growthrate
     integer :: decayrate
     integer :: lai_sat
     integer :: fpar_sat
     integer :: fvcover
     integer :: tmin_tave
     integer :: vpd_tave
     integer :: photo_tave
     integer :: rg_tave
     integer :: rain_tave
     integer :: chill_tave
  end type paramidxtype
  type(paramidxtype) :: paramidx

  type :: forceidxtype
     integer :: tmin
     integer :: vpd
     integer :: photo
     integer :: rg
     integer :: rain
     integer :: tmean
  end type forceidxtype
  type(forceidxtype) :: forceidx

  type :: obseridxtype
     integer :: fpar
     integer :: lai
  end type obseridxtype
  type(obseridxtype) :: obseridx

  ! Naming of variables
  character(len=100), allocatable :: state_names(:)
  character(len=100), allocatable :: state_longnames(:)
  character(len=100), allocatable :: state_units(:)

  character(len=100), allocatable :: param_names(:)
  character(len=100), allocatable :: param_longnames(:)
  character(len=100), allocatable :: param_units(:)

  character(len=100), allocatable :: force_names(:)
  character(len=100), allocatable :: force_longnames(:)
  character(len=100), allocatable :: force_units(:)

  character(len=100), allocatable :: obser_names(:)
  character(len=100), allocatable :: obser_longnames(:)
  character(len=100), allocatable :: obser_units(:)

  ! Output Options
  integer, allocatable :: output_state(:) ! index of saved states
  integer, allocatable :: output_param(:) ! index of saved parameters
  integer, allocatable :: output_force(:) ! index of saved forcing variables
  integer, allocatable :: output_obser(:) ! index of saved observations

  real(kind=4), allocatable :: state_scale(:)   ! output NetCDF scale_factor of states (1 : no scale_factor)
  real(kind=4), allocatable :: param_scale(:)   ! output NetCDF scale_factor of parameters (1 : no scale_factor)
  real(kind=4), allocatable :: force_scale(:)   ! output NetCDF scale_factor of forcing (1 : no scale_factor)
  real(kind=4), allocatable :: obser_scale(:)   ! output NetCDF scale_factor of observations (1 : no scale_factor)

  real(kind=4), allocatable :: state_offset(:)   ! output NetCDF add_offset of states (0 : no add_offset)
  real(kind=4), allocatable :: param_offset(:)   ! output NetCDF add_offset of parameters (0 : no add_offset)
  real(kind=4), allocatable :: force_offset(:)   ! output NetCDF add_offset of forcing (0 : no add_offset)
  real(kind=4), allocatable :: obser_offset(:)   ! output NetCDF add_offset of observations (0 : no add_offset)

  ! output types allowed: byte, short, int, float, double
  character(len=100), allocatable :: state_type(:)   ! output NetCDF variable type of states 
  character(len=100), allocatable :: param_type(:)   ! output NetCDF variable type of parameters 
  character(len=100), allocatable :: force_type(:)   ! output NetCDF variable type of forcing
  character(len=100), allocatable :: obser_type(:)   ! output NetCDF variable type of observations 

  integer, allocatable :: ensindex(:)             ! index (1..nrens) for identification of ensemble members
  integer, allocatable :: pftindex(:)             ! index (1..npft) for identification of PFT's
  real(kind=4), allocatable :: hgtlevels(:)       ! elevation levels (minhgt .. maxhgt)

  real(kind=4) :: minsumpft                       ! minimum % pft which selected pft's cover for each observation
  integer :: nselect_pft                          ! number of pft's to use
  integer, allocatable :: select_pft(:)           ! indices of pft's to use
  logical, allocatable :: pft_is_selected(:)      ! logical array storing selected pft's

  real(kind=4) :: minsumhgt                       ! minimum % hgt which selected hgt's cover for each observation
  integer :: nselect_hgt                          ! number of hgt's to use
  integer, allocatable :: select_hgt(:)           ! indices of hgt's to use
  logical, allocatable :: hgt_is_selected(:)      ! logical array storing selected hgt's

  real(kind=4) :: maxpctwater                     ! maximum percentage (%/100) of water area needed to run prediction/assimilation

  ! Geographic arrays
  type :: projtype
     character(len=100) :: name
     real(kind=4)       :: xmin
     real(kind=4)       :: xmax
     real(kind=4)       :: ymin
     real(kind=4)       :: ymax
     real(kind=4)       :: dx
     real(kind=4)       :: dy
     integer            :: nx
     integer            :: ny
     character(len=100) :: type
     character(len=100) :: optargs
  end type projtype
 
  ! model / assimilation grid (for this process)
  type(projtype) :: geographic
  real(kind=4), allocatable :: xcoord(:,:)        ! x values of the assimilaton grid  (projection type units)
  real(kind=4), allocatable :: ycoord(:,:)        ! y values of the assimilation grid (projection type units)
  real(kind=4), allocatable :: lat(:,:)           ! latitudes of the assimilaton grid (degrees)
  real(kind=4), allocatable :: lon(:,:)           ! longitudes of the assimilation grid (degrees)
  real(kind=4), allocatable :: hgt(:,:)           ! elevations of the assimilation grid (m)
  real(kind=4), allocatable :: area(:,:)          ! area covered by grid point
!  real(kind=4), allocatable :: pft_pct(:,:,:)     ! pft distribution of the assimilation grid (%/100)
  real(kind=4), allocatable :: hgt_pct(:,:,:)     ! elevation distribution of the assimilation grid (m)

  ! output grid (for this process)
  type(projtype) :: geographic_out
  integer :: xcoord_offset                        ! offset in x direction for this process (or 1 if no offset)
  integer :: ycoord_offset                        ! offset in y direction for this process (or 1 if no offset)
  real(kind=4), allocatable :: xcoord_out(:)      ! x values of the output grid for this process (projection type units)
  real(kind=4), allocatable :: ycoord_out(:)      ! y values of the output grid for this process (projection type units)
  real(kind=4), allocatable :: pft_pct_out(:,:,:) ! pft distribution of the output grid  (%/100)
  real(kind=4), allocatable :: hgt_pct_out(:,:,:) ! elevation distribution of the output grid (%/100)
  real(kind=4), allocatable :: hgt_out(:,:)       ! mean topography of the output grid (m)

  ! full output grid (over all processes)
  type(projtype) :: geographic_full
  real(kind=4), allocatable :: xcoord_full(:)     ! x values (projection type units)
  real(kind=4), allocatable :: ycoord_full(:)     ! y values  (projection type units)

  ! grid averaged quantities
  real(kind=4), allocatable :: pft_pct_ave(:)     ! pft distribution averaged over the output grid  (%/100)
  real(kind=4), allocatable :: hgt_pct_ave(:)     ! elevation distribution averaged over the output grid (%/100)
  real(kind=4) :: hgt_ave                         ! averaged height of output grid (m)
  real(kind=4) :: hgt_stddev                      ! height standard deviation of output grid (m)

  ! DEPRECATED -- REMOVE SOONER OR LATER
  real(kind=4), allocatable :: lat_out(:,:)       ! latitudes of the output grid (degrees)
  real(kind=4), allocatable :: lon_out(:,:)       ! longitudes of the output grid (degrees)
  real(kind=4) :: pollon                          ! longitude of North Pole (default 180.)
  real(kind=4) :: pollat                          ! latitude of North Pole (default 90.)
  real(kind=4) :: latmin,latmax,lonmin,lonmax     ! geographic boundaries (degrees) of output grid
  real(kind=4) :: dlat, dlon                      ! lon/lat step (degrees) of assimilation grid
  real(kind=4) :: dlat_out, dlon_out              ! lon/lat step (degrees) of output grid
  integer :: nlat, nlon                           ! lon/lat elements of assimilation grid
  integer :: npft, nhgt                           ! number of PFT and HGT elements of output grid
  integer :: nlat_out, nlon_out                   ! lon/lat elements of output grid
  integer :: npft_out, nhgt_out                   ! number of PFT and HGT elements of output grid

  ! Array size parameters
  integer :: nrens    ! number of model ensemble members
  integer :: nptstate ! total number of state dimensions
  integer :: nptparam ! total number of parameter dimensions
  integer :: nptforce ! total number of forcing dimensions
  integer :: nptobser ! total number of observation dimensions

  ! total observations cumulated over simulation time
  ! by pft, cumulated in time and over all processes
  real(kind=4), allocatable :: cumobs(:)       


  ! MPI structure
  integer :: mpierr     ! error flag for MPI routines
  integer :: mpiid      ! process counter for MPI routines
  integer :: numprocs   ! total processes used for MPI routines
  integer :: nxprocs    ! number of processes in spatial x direction
  integer :: nyprocs    ! number of processes in spatial y direction
  integer :: xproc      ! process index (1..nxprocs)
  integer :: yproc      ! process index (1..nyprocs)
  integer :: lonoffset  ! output offset in lon direction for each process
  integer :: latoffset  ! output offset in lat direction for each process

  ! namelist
  character(len=100) :: casename        ! case name of model / analysis run

  integer :: landcover_screen           ! land use type to assimilate (MODIS UMD type)
  integer :: pft_screen                 ! plant functional type to assimilate (NCAR PFT's)
  integer :: subgrid                    ! 1 .. n, for n x n spatial subgrid points in output (only works in regional mode)
  integer :: filelength                 ! output file length: (-1) unlimited, (1) daily, (2) monthly, (3) yearly)
  logical :: restart                    ! true: restart run, false: first run
  logical :: restartdate                ! true: use start date from restart file instead of start date from namelist
  logical :: initial                    ! true: first time step
  logical :: analysis                   ! perform analysis step
  logical :: isprediction               ! perform prediction on this node
  logical :: original                   ! run original scheme with original parameters
  logical :: simulator                  ! simulate observations with prescribed parameters instaed of reading true observations
  logical :: regional                   ! regional analysis/prediction with subgrid-scale PFT/height distribution 
  logical :: hgtscreen                  ! screen observations by min/max elevation (true) or not
  character(len=100) :: package_name    ! code name
  character(len=20) :: package_version  ! code version number
  character(len=200) :: projectdir      ! project root directory
  character(len=200) :: datadir         ! data root directory
  character(len=200) :: batchdir        ! batch script directory
  character(len=200) :: outdir          ! output directory
  character(len=200) :: namelistdir     ! directory to read namelists from
  character(len=200) :: restartdir      ! restart directory
  integer :: dtphen                     ! phenology model integration time step (s)
  integer :: dtmet                      ! meteorological driver time step (s)
  integer :: dtmod                      ! model driver integration time step (s)
  integer :: dtout                      ! output write time step (s)
  character(len=8) :: date_start        ! start date (YYYY, YYYYMM or YYYYMMDD)
  character(len=8) :: date_end          ! end date (YYYY, YYYYMM or YYYYMMDD)

  logical :: use_forcing                ! use forcing (temperature-based)
  logical :: use_chilling               ! use chilling (temperature-based)
  logical :: use_moisture               ! use moisture limitation (vapor pressure deficit-based)
  logical :: use_light                  ! use light limitation (photoperiod- or global radiation-based)
  logical :: use_rainfall               ! use rainfall instead of vapor pressure deficit for moisture limitation
  logical :: use_photo                  ! use photoperiod instead of global radiation for light limitation

  ! output flags
  logical :: verbose      ! print messages
  logical :: averagepft   ! output pft-averaged fields for each grid cell
  logical :: averagehgt   ! output hgt-averaged fields for each grid cell
  logical :: averageens   ! output ensemble-averaged fields for each grid cell
  logical :: compressens  ! create output ensemble which only has two members mean + stddev
  logical :: saveobs      ! save observations on separate file in the output grid
  logical :: saveselpft   ! save only selected pft's 

  ! parameters
  real(kind=4), parameter :: nodata = -9999. ! missing value (used for flagging e.g. bad FLUXNET or MODIS observations)

contains

  ! below follow general service and IO modules needed for running the forward and assimiliation model
  ! IMPORTANT : they are self-containing and do not depend on other modules
  function get_config(string) result(retval)

    ! calls a c code which contains the information of the GNU autoconf config.h system configuration
    ! and returns the version number of this program for internal fortran use. This call is needed because
    ! we cannot directly include config.h in fortran.

    ! DEPENDENCIES: read_config.c

    use iso_c_binding

    implicit none

    ! arguments
    character(len=*), intent(inout) :: string ! contains the configuration parameter name on input
    ! and the parameter value on output (or N/A if not found)    
    ! result
    integer :: retval

    ! local
    integer :: length
    integer :: i

    length = len(string)

    ! add '\0' character to end of string (C standard string terminator)
    string = trim(string)//c_null_char

    retval = read_config(length,string)

    ! remove '\0' string termination
    i = scan(string,c_null_char)
    if (i.gt.1) string = string(1:i-1)

    return

  end function get_config

  function Translate(String,InTable,OutTable)
    implicit none

    character(len=*),intent(in) :: String,InTable,OutTable
    character(len=len(String)) :: Translate
    integer :: i , p , l

    Translate = String
    l = min(len(InTable),len(OutTable))
    do p = 1 , len(String)
       i = index(InTable(1:l),String(p:p))
       if ( i>0 ) Translate(p:p) = OutTable(i:i)
    enddo
  end function Translate

  function UpperCase(string)
    implicit none
    character(len=*),intent(in) :: string
    character(len=len(string)) :: UpperCase
    UpperCase = Translate(string,'abcdefghijklmnopqrstuvwxyz' &
         ,'ABCDEFGHIJKLMNOPQRSTUVWXYZ')
  end function UpperCase

  function LowerCase(string)
    implicit none
    character(len=*),intent(in) :: string
    character(len=len(string)) :: LowerCase
    LowerCase = Translate(string,'ABCDEFGHIJKLMNOPQRSTUVWXYZ' &
         ,'abcdefghijklmnopqrstuvwxyz')
  end function LowerCase

  subroutine init_random

    ! initializes random number generator with system clock
    implicit none

    integer :: i, s, clock
    integer, allocatable :: seed(:)

    call random_seed(size = s)
    allocate(seed(s))
    call system_clock(count=clock)
    i=1
    seed = clock + 23 * mpiid + 37 * (/ (i - 1, i = 1, s) /)

    ! uncomment line below to generate 1:1 reproducible
    ! assimilation experiments (it always initializes 
    ! the random number generator with the same numbers)

    ! seed(:) = 200

    call random_seed(put = seed)
    deallocate(seed)

  end subroutine init_random

  subroutine calc_mean(A,ave)
    ! calculates the ensemble mean of a vector A
    implicit none

    ! arguments
    real(kind=4), intent(in)  :: A(:)
    real(kind=4), intent(out) :: ave

    ! local
    integer :: nrens
    integer :: j

    nrens = size(A,1)

    ave=A(1)
    do j=2,nrens
       ave=ave+A(j)
    enddo
    ave=(1./real(nrens))*ave

  end subroutine calc_mean

  subroutine calc_variance(A,ave,var)
    ! calculates the ensemble variance of a vector A
    implicit none

    ! arguments
    real(kind=4), intent(in)  :: A(:)
    real(kind=4), intent(in)  :: ave
    real(kind=4), intent(out) :: var

    ! local
    integer :: nrens
    integer :: j

    nrens = size(A,1)

    var=0.
    do j=1,nrens
       var=var+(A(j)-ave)*(A(j)-ave)
    enddo
    var=(1./real(nrens-1))*var

  end subroutine calc_variance

  subroutine separate_minmax(A,B,optmindiff,optmaxdiff)
    ! limit the mean of vector A to be separated from the mean of vector B by a minimum and maximum difference. 
    ! If these limits are crossed, then the mean of the vector is shifted back to below and above the mean of A and B. 
    ! Individual members of the vector can therefore still be below or above the threshold. 
    ! This guarantees that the variances of vectors A and B are maintained.

    ! limit: (mean(B) - mean(A) ) > mindiff
    ! limit: (mean(B) - mean(A) ) < maxdiff

    implicit none

    ! arguments
    real(kind=4), intent(inout) :: A(:)
    real(kind=4), intent(inout) :: B(:)
    real(kind=4), intent(in), optional :: optmindiff
    real(kind=4), intent(in), optional :: optmaxdiff

    ! local variables
    real(kind=4) :: A_mean, B_mean
    real(kind=4) :: mindiff, maxdiff

    if (present(optmindiff)) then
       mindiff = optmindiff
    else
       mindiff = 1.e-3
    end if

    if (present(optmaxdiff)) then
       maxdiff = optmaxdiff
    else
       maxdiff = 1.e3
    end if

    if (maxdiff < mindiff) then
       write(*,'(A)') "separate_minmax: Maximum difference has to be above minimum difference"
    end if

    call calc_mean(A,A_mean)
    call calc_mean(B,B_mean)

    if ((B_mean - A_mean) .lt. mindiff) then
       A = A - 0.5*(A_mean - B_mean + mindiff)
       B = B + 0.5*(A_mean - B_mean + mindiff)
    endif

    if ((B_mean - A_mean) .gt. maxdiff) then
       A = A - 0.5*(A_mean - B_mean + maxdiff)
       B = B + 0.5*(A_mean - B_mean + maxdiff)
    endif

  end subroutine separate_minmax

  function calc_distance(lon1,lat1,lon2,lat2)

    ! calculates geodesic distance in units [km] between two points on a sphere
    ! note: the formula assumes a spherical earth (simplification)

    ! http://en.wikipedia.org/wiki/Great-circle_distance

    use parameter_module

    implicit none

    ! arguments
    real(kind=4) :: lon1, lat1, lon2, lat2

    ! local variables

    ! function value
    real(kind=4) :: calc_distance ! (km)

    calc_distance=2. * earthradius * asin(sqrt((sin(0.5*degrad*(lat1-lat2)))**2 + &
         cos(degrad*lat1)*cos(degrad*lat2)*(sin(0.5*degrad*(lon1-lon2)))**2))

  end function calc_distance

end module driver_module
