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

module meteorology_average_module

contains

  subroutine meteorology_average_daily(newday, dtmet, year, doy, &
       tair, mr, rg, rain, psurf, lat, hgt, hgtout, &
       tmin_daily, tmean_daily, vpd_daily, photo_daily, rg_daily, rain_daily, &
       missing)

    ! Description
    ! -----------
    ! Creates daily minimum / average / sums of meteorological forcing with or without
    ! orographic downscaling of temperature, vapor pressure deficit and global radiation.
    !
    ! The routine is written to be compatible with the phenology prediction and thus only
    ! creates daily fields of very specific meteorological variables.
    !
    ! Input: air temperature, water vapor mixing ratio, global radiation, rainfall, surface pressure.
    ! Output: minimum daily air temperature, mean daily vpd, daily photoperiod, mean daily
    !         global radiation, cumulated daily rainfall.
    !
    ! Notes
    ! -----
    ! 1) call it after each read operation of forcing data from the coupled atmospheric
    !    model or from reanalysis data
    ! 2) call it at least at each start of the day (with the newday flag set to true)
    ! 3) if no elevation downscaling is needed, set the input fields hgtout == hgt
    ! 4) the number of input and output points may differ. Two modes are accetped:
    !    a) npt = npt2        : output grid points = input grid points. No further explanation needed.
    !    b) npt = 1; npt2 > 1 : single input grid point is used to calculate daily averages for npt2 
    !                           output grid points
    ! 5) the output fields have a time history and therefore have to be treated like prognostic
    !    variables in the host model. They should be included in the restart files. If this is not 
    !    possible, they can also be initialized at restart. Set the optional missing keyword to
    !    the value that the output fields have when they are not initialized. The code will produce
    !    an automatic spin-up by taking the first instantaneous value as minimum or daily mean in 
    !    this case
    ! 6) the output fields can be calculated for n subgrid-scale elevation classes by 
    !    grid point by running this routine n times with hgtout=elevation_class(n).
    !    That however means that the following multi-day averaging and successively also
    !    the phenological prediction needs to be carried out over the elevation classes too.
    !    And the resulting FPAR and LAI fields from the phenological prediction will need
    !    to be transformed to a grid-scale value by area-weighted summation by pft and elevation
    !    class. This is how it was done in the paper
    !
    ! References
    ! ----------
    ! R. Stöckli, T. Rutishauser, I. Baker, M. Liniger, and A. S. Denning (2011): 
    ! A global reanalysis of vegetation phenology. J. Geophys. Res. - Biogeosciences, 
    ! 116 (G03020), doi: 10.1029/2010JG001545.
    !
    ! R. Stöckli, T. Rutishauser, D. Dragoni, J. O'Keefe, P. E. Thornton, M. Jolly, 
    ! L. Lu, and A. S. Denning (2008): Remote sensing data assimilation for a prognostic 
    ! phenology model. J. Geophys. Res. - Biogeosciences, 113 (G4), doi: 10.1029/2008JG000781.
    !
    ! W. M. Jolly, R. Nemani, and S. W. Running (2005). A generalized, bioclimatic 
    ! index to predict foliar phenology in response to climate. Global Change Biol., 11:619-632
    !
    ! Author
    ! ------
    ! Dr. Reto Stöckli
    ! Blue Marble Research
    ! Freiestrasse 41
    ! 3012 Bern, Switzerland
    ! Phone: +41 31 535 3989
    ! Email: reto.stockli@gmail.com
    ! Web: http://www.bluemarble.ch
    ! 
    ! Acknowledgments
    ! ---------------
    ! The NASA Energy and Water Cycle Study (NEWS) grant No. NNG06CG42G is the main funding source of this study.


    implicit none

    ! Arguments

    ! Flags & Parameters
    logical, intent(in) :: newday                    ! flag start of a new day with TRUE
    integer, intent(in) :: dtmet                     ! time step of meteorology input data (s)
    integer, intent(in) :: year                      ! Year (YYYY)
    integer, intent(in) :: doy                       ! day of year (1..366)

    ! Input: Instantaneous Forcing
    real(kind=4), intent(in) :: tair            ! air temperature (K)
    real(kind=4), intent(in) :: mr              ! water vapor mixing ratio (g/kg)
    real(kind=4), intent(in) :: rg              ! global radiation (W/m2)
    real(kind=4), intent(in) :: rain            ! rainfall (mm)
    real(kind=4), intent(in) :: psurf           ! surface pressure (Pa) (at hgt, not at sea-level!)

    ! Geographic information of instantaneous forcing
    real(kind=4), intent(in) :: lat             ! latitude (degrees N)
    real(kind=4), intent(in) :: hgt             ! input elevation (meters above sea level)
    real(kind=4), intent(in) :: hgtout          ! output elevation (meters above sea level)

    ! Output: Daily Forcing
    real(kind=4), intent(inout) :: tmin_daily(:)  ! minimum daily temperature (K)
    real(kind=4), intent(inout) :: tmean_daily(:) ! mean daily temperature (K)
    real(kind=4), intent(inout) :: vpd_daily(:)   ! mean daily VPD (mb)
    real(kind=4), intent(inout) :: photo_daily(:) ! cumulative daily Photoperiod (hours)
    real(kind=4), intent(inout) :: rg_daily(:)    ! mean daily global radiation (W/m2)
    real(kind=4), intent(inout) :: rain_daily(:)  ! cumulative daily rainfall (mm)

    ! Initial value for uninitialized Multi-Day Averaged Forcing
    ! set this optional keyword if these variables are not initialized
    real(kind=4), intent(in), optional :: missing    ! optional missing data value

    ! Local Variables
    integer :: npt                     ! number of daily averaged forcing data points

    integer :: n
    real(kind=4) :: weightdaily                      ! weighting (1) daily .. (0) instantaneous values
    real(kind=4) :: doy_adj                          ! adjusted astronomical day of year
    real(kind=4), allocatable :: psurf_local(:)      ! surface pressure (Pa) at output elevation

    ! Parameters
    real(kind=4), parameter :: pi  = 3.141592         ! the constant pi
    real(kind=4), parameter :: decmax = 23.441        ! maximum declination
    real(kind=4), parameter :: degrad = 0.0174532935  ! pi / 180

    real(kind=4), parameter :: tair_lapse = -0.0065   ! temperature lapse rate: -0.65 K / 100 m
    real(kind=4), parameter :: rg_lapse = 0.003       ! global radiation lapse rate: 0.3 W/m2 / 100m

    ! diagnose array sizes
    npt = size(tmin_daily)

    ! allocate local arrays
    allocate(psurf_local(npt))

    ! weight factor for daily means versus instantanous values
    weightdaily  = exp(-1. / (1. * 86400. / real(dtmet)))

    ! reset time-integrated / averaged properties each day at the first time step
    if (newday) then
       rain_daily(:)=0.
       tmin_daily(:)=9999.
    end if

    ! check if we have the absolutely first initial time (all forcings are set to missing)
    ! then initialize running mean variables with instantaneous variables
    if (present(missing)) then
       if (vpd_daily(1).eq.missing) then
          weightdaily  = 0.
       endif
    endif

    ! calculate daily averages/sums/minimums

    ! daily minimum temperature [K]
    do n=1,npt
       tmin_daily(n) = min(tmin_daily(n), tair+(hgtout - hgt)*tair_lapse) 
    enddo

    ! daily mean temperature [K]
    tmean_daily(:) = weightdaily * tmean_daily(:) + (1.-weightdaily) * (tair+(hgtout - hgt)*tair_lapse)

    ! calculate local surface pressure
    psurf_local(:) = psurf * exp( -9.80665 * 0.0289644 * (hgtout - hgt) / (8.31432 * tair) )

    ! daily mean vapor pressure deficit [mb] fromm mixing ratio [g/kg], 
    ! temperature [K] and air pressure [Pa]
    ! elevation downscaling: assume constant mixing ratio but temperature varying with lapse rate 
    vpd_daily(:) = weightdaily * vpd_daily(:) + (1.-weightdaily) * &
         ( 6.11*10**(7.5*(tair+(hgtout - hgt)*tair_lapse-273.16) / &
         (237.7+(tair+(hgtout - hgt)*tair_lapse-273.16))) - &
         mr * psurf_local(:) / 100. / (mr+622.) )

    ! daily mean RG [W/m2]
    rg_daily(:)  = weightdaily * rg_daily(:) + (1.-weightdaily) * &
         max((rg +(hgtout - hgt)*rg_lapse),0.)

    ! daily photoperiod [h]
    ! formula from http://herbert.gandraxa.com/length_of_day.xml
    ! astronomical doy shifts by 0.25 days per year after each "leap" year
    ! this is an approximate calculation and ignores years 100 and 400 for simplicity
    ! photoperiod is defined here as the integrated time over the day where the sun is
    ! has elevation of 5 or more degrees above the horizon
    doy_adj = real(doy) + real(mod(year,4))*0.25
    photo_daily(:) = acos(1. - min(max(1. - tan(degrad*lat) * &
         tan(cos(pi * doy_adj / 182.625) * decmax * degrad) - &
         tan(degrad*5.) / cos(degrad*lat),0.),2.)) / pi * 24.

    ! daily cumulative rainfall [mm]
    rain_daily(:) = rain_daily(:) + rain

    ! deallocate local arrays
    deallocate(psurf_local)

  end subroutine meteorology_average_daily

  subroutine meteorology_average_multiday(lat, tmin_daily, tmean_daily, &
       vpd_daily, photo_daily, rg_daily, rain_daily, &
       tmin_ave, vpd_ave, photo_ave, rg_ave, rain_ave, chill_ave, &
       tmin_tave, vpd_tave, photo_tave, rg_tave, rain_tave, chill_tave, &
       dtphen, doy, missing)

    ! Description
    ! -----------
    ! Creates multi-day averages of daily meteorology input.
    !
    ! The routine is written to be compatible with the phenology prediction and thus only
    ! creates output fields of very specific daily meteorological variables.
    !
    ! Input: minimum daily air temperature, mean daily vpd, daily photoperiod, mean daily
    !         global radiation, cumulated daily rainfall.
    ! Output: averaged minimum daily air temperature, averaged mean daily vpd, 
    !         averaged daily photoperiod, averaged mean daily
    !         global radiation, averaged cumulated daily rainfall.
    !
    ! Notes
    ! -----
    ! 1) call it before each call to the phenological prediction (e.g. daily)
    ! 2) the number of input points can be different from the number of output points.
    !    Two modes of operation are possible:
    !    npt = npt2 : average each input grid point with its own parameter set 
    !                 (externally loop over grid points and optionally elevation classes)
    !    npt > 1, npt2 = 1 : average several grid points with a single parameter set 
    !                        (externally loop over pft's and optionally elevation classes)
    ! 3) the time averaging parameters tave_xyz are phenological parameters and vary by plant
    !    functional type. Please make sure that the pft-specific values are correctly supplied.
    ! 4) the output fields have a time history and therefore have to be treated like prognostic
    !    variables in the host model. They should be included in the restart files. If this is not 
    !    possible, they can also be initialized at restart. Set the optional missing keyword to
    !    the value that the output fields have when they are not initialized. The code will produce
    !    an automatic spin-up by taking the first daily mean / minimum value as first guess. The 
    !    length of spin up will vary depending on the magnitude of the time averaging parameter.
    !
    ! References
    ! ----------
    ! R. Stöckli, T. Rutishauser, I. Baker, M. Liniger, and A. S. Denning (2011): 
    ! A global reanalysis of vegetation phenology. J. Geophys. Res. - Biogeosciences, 
    ! 116 (G03020), doi: 10.1029/2010JG001545.
    !
    ! R. Stöckli, T. Rutishauser, D. Dragoni, J. O'Keefe, P. E. Thornton, M. Jolly, 
    ! L. Lu, and A. S. Denning (2008): Remote sensing data assimilation for a prognostic 
    ! phenology model. J. Geophys. Res. - Biogeosciences, 113 (G4), doi: 10.1029/2008JG000781.
    !
    ! W. M. Jolly, R. Nemani, and S. W. Running (2005). A generalized, bioclimatic 
    ! index to predict foliar phenology in response to climate. Global Change Biol., 11:619-632
    !
    ! Author
    ! ------
    ! Dr. Reto Stöckli
    ! Blue Marble Research
    ! Freiestrasse 41
    ! 3012 Bern, Switzerland
    ! Phone: +41 31 535 3989
    ! Email: reto.stockli@gmail.com
    ! Web: http://www.bluemarble.ch
    ! 
    ! Acknowledgments
    ! ---------------
    ! The NASA Energy and Water Cycle Study (NEWS) grant No. NNG06CG42G is the main funding source of this study.

    implicit none

    ! Arguments

    ! Geographic Boundary Conditions
    real(kind=4), intent(in) :: lat(:)             ! latitude (degrees north)

    ! Daily Averaged Forcing
    real(kind=4), intent(in) :: tmin_daily(:)      ! minimum daily temperature (K)
    real(kind=4), intent(in) :: tmean_daily(:)     ! mean daily temperature (K)
    real(kind=4), intent(in) :: vpd_daily(:)       ! mean daily VPD (mb)
    real(kind=4), intent(in) :: photo_daily(:)     ! cumulative daily Photoperiod (hours)
    real(kind=4), intent(in) :: rg_daily(:)        ! mean daily global radiation (W/m2)
    real(kind=4), intent(in) :: rain_daily(:)      ! cumulative daily rainfall (mm)

    ! Multi-Day Averaged Forcing
    real(kind=4), intent(inout) :: tmin_ave(:)     ! average minimum daily temperature (K)
    real(kind=4), intent(inout) :: vpd_ave(:)      ! average mean daily VPD (mb)
    real(kind=4), intent(inout) :: photo_ave(:)    ! average cumulative daily Photoperiod (hours)
    real(kind=4), intent(inout) :: rg_ave(:)       ! average mean daily global radiation (W/m2)
    real(kind=4), intent(inout) :: rain_ave(:)     ! average cumulative daily rainfall (mm)
    real(kind=4), intent(inout) :: chill_ave(:)    ! average cumulative chilling state (K*days)

    ! Time Averaging Parameters
    real(kind=4), intent(in) :: tmin_tave(:)      ! n day averaging interval for temperature (days)
    real(kind=4), intent(in) :: vpd_tave(:)       ! n day averaging interval for vpd (days)
    real(kind=4), intent(in) :: photo_tave(:)     ! n day averaging interval for photoperiod (days)
    real(kind=4), intent(in) :: rg_tave(:)        ! n day averaging interval for rg (days)
    real(kind=4), intent(in) :: rain_tave(:)      ! n day averaging interval for rain (days)
    real(kind=4), intent(in) :: chill_tave(:)     ! n day averaging interval for chilling state (days)

    ! Time Step Parameter
    integer, intent(in) :: dtphen                    ! prediction time step [s]
    integer, intent(in) :: doy                       ! day of year [days]

    ! Initial value for uninitialized Multi-Day Averaged Forcing
    ! set this optional keyword if these variables are not initialized
    real(kind=4), intent(in), optional :: missing

    ! Local variables

    ! Array size
    integer :: npt                       ! number of forcing points
    integer :: npt2                      ! number of parameter points

    ! Diagnostics
    real(kind=4), allocatable :: weight_tmin(:)
    real(kind=4), allocatable :: weight_vpd(:)
    real(kind=4), allocatable :: weight_photo(:)
    real(kind=4), allocatable :: weight_rg(:)
    real(kind=4), allocatable :: weight_rain(:)
    real(kind=4), allocatable :: weight_chill(:)

    ! diagnose array sizes
    npt = size(lat)
    npt2 = size(tmin_tave)

    if ((npt.ne.npt2).and.(npt2.ne.1)) then
       write(*,'(A)') "Forcing and Parameter need either to be of same dimension (one parameter set per grid point)"
       write(*,'(A)') "Or Parameters have to be scalar (same parameters for all grid points)"
       stop
    endif

    ! allocate local arrays
    allocate(weight_tmin(npt))
    allocate(weight_vpd(npt))
    allocate(weight_photo(npt))
    allocate(weight_rg(npt))
    allocate(weight_rain(npt))
    allocate(weight_chill(npt))

    ! factors to average daily mean / sums / minimum values
    ! to multi-day averages 

    if (npt.ne.npt2) then
       weight_tmin  = exp(-real(dtphen) / 86400. / min(max(tmin_tave(1),1.0),1000.))
       weight_vpd   = exp(-real(dtphen) / 86400. / min(max(vpd_tave(1),1.0),1000.))
       weight_photo = exp(-real(dtphen) / 86400. / min(max(photo_tave(1),1.0),1000.))
       weight_rg    = exp(-real(dtphen) / 86400. / min(max(rg_tave(1),1.0),1000.))
       weight_rain  = exp(-real(dtphen) / 86400. / min(max(rain_tave(1),1.0),1000.))
       weight_chill = exp(-real(dtphen) / 86400. / min(max(chill_tave(1),1.0),1000.))
    else
       weight_tmin  = exp(-real(dtphen) / 86400. / min(max(tmin_tave,1.0),1000.))
       weight_vpd   = exp(-real(dtphen) / 86400. / min(max(vpd_tave,1.0),1000.))
       weight_photo = exp(-real(dtphen) / 86400. / min(max(photo_tave,1.0),1000.))
       weight_rg    = exp(-real(dtphen) / 86400. / min(max(rg_tave,1.0),1000.))
       weight_rain  = exp(-real(dtphen) / 86400. / min(max(rain_tave,1.0),1000.))
       weight_chill = exp(-real(dtphen) / 86400. / min(max(chill_tave,1.0),1000.))
    endif

    ! check if multi-day averaged forcing is uninitialized
    ! if yes, initialize with daily mean / sums / minimum values
    if (present(missing)) then 
       if (tmin_ave(1).eq.missing) tmin_ave = tmin_daily
       if (vpd_ave(1).eq.missing) vpd_ave = vpd_daily
       if (photo_ave(1).eq.missing) photo_ave = photo_daily
       if (rg_ave(1).eq.missing) rg_ave = rg_daily
       if (rain_ave(1).eq.missing) rain_ave = rain_daily
       if (chill_ave(1).eq.missing) chill_ave = 150.
    endif

    ! calculate multi-day averages from daily averages (except photoperiod, which is not averaged over days)
    tmin_ave  = weight_tmin * tmin_ave  + (1.-weight_tmin) * tmin_daily
    vpd_ave   = weight_vpd * vpd_ave   + (1.-weight_vpd) * vpd_daily
    photo_ave = weight_photo * photo_ave + (1.-weight_photo) * photo_daily
    rg_ave    = weight_rg * rg_ave    + (1.-weight_rg) * rg_daily
    rain_ave  = weight_rain * rain_ave  + (1.-weight_rain) * rain_daily

    ! reset chill state between
    ! 1st October to 1st December (NH) or 1st April to 1 June (SH) 
    ! under the constraint that daylength is below 10 hours
    if ((doy.ge.274).and.(doy.le.335)) then
       where ((lat.gt.0.).and.(photo_daily.lt.10.).and.(chill_ave.gt.100.))
          chill_ave = 0.
       end where
    end if

    if ((doy.ge.91).and.(doy.le.152)) then
       where ((lat.le.0.).and.(photo_daily.lt.10.).and.(chill_ave.gt.100.))
          chill_ave = 0.
       end where
    end if

    ! Dormphot chill model by Caffarra 2011
    !  chill_ave = chill_ave  + 1. / (1. + exp(0.03 * (tmean_daily - 273.16 - 13.89)**2 + (tmean_daily - 273.16 - 13.89)))

    ! General chill model
    chill_ave = chill_ave + (1. - min(1.,max(0.,(tmean_daily - 278.16)/5.)))

    ! deallocate local arrays
    deallocate(weight_tmin)
    deallocate(weight_vpd)
    deallocate(weight_photo)
    deallocate(weight_rg)
    deallocate(weight_rain)
    deallocate(weight_chill)
    
  end subroutine meteorology_average_multiday

end module meteorology_average_module
