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

module phenology_model_module

contains

  subroutine phenology_prognostic_gsi( fpar, lai, gsi, &
       temp_fac, moist_fac, light_fac, force_fac, chill_fac, &
       temp_ave, moist_ave, light_ave, chill_ave, &
       temp_min, temp_max, moist_min, moist_max, light_min, light_max, chill_min, chill_max, &
       lai_min, lai_max, growthrate, decayrate, lai_sat, fpar_sat, fvcover, &
       use_forcing, use_chilling, use_moisture, use_light, use_rainfall, use_photo)

    ! Description
    ! -----------
    ! Creates phenological prediction from averaged meteorological forcing, previous phenological
    ! state and plant-functional-type specific phenological parameters.
    !
    ! Input:  previous leaf area index, averaged minimum daily air temperature, 
    !         averaged mean daily vpd, averaged daily photoperiod, averaged mean daily
    !         global radiation, averaged cumulated daily rainfall, phenological parameters (17).
    ! Output: updated leaf area index, diagnosed FPAR, Growing Season Index GSI, 
    !         temperature limitation factor, moisture limitation factor and light limitation factor
    !
    ! Notes
    ! -----
    ! 1) call it after each calculation of multi-day averaged meteorological drivers
    ! 2) the number of input points can be different from the number of output points.
    !    Two modes of operation are possible:
    !    npt = npt2 : predict phenology for each input grid point with its own parameter set 
    !                 (externally loop over grid points and optionally elevation classes)
    !    npt > 1, npt2 = 1 : predict phenology for several grid points with a single parameter set 
    !                        (externally loop over pft's and optionally elevation classes)
    ! 3) the phenological parameters like tminmin, decayrate or fpar_sat vary by plant
    !    functional type. Please make sure that the pft-specific values are correctly supplied.
    ! 4) LAI is a prognostic and therefore it has to be treated like prognostic
    !    variables in the host model. It has to be included in the host model's restart files. 
    !    There is no spin up procedure included in this routine here. Either LAI is set to
    !    a arbitrary value globally (e.g. 2.0) externally to this routine in the first year,
    !    or a climatology is used as first guess. In either case the spin up will be around
    !    one year to get the model's LAI into a usable state.
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

    ! Prognostic States
    real(kind=4), intent(inout) :: lai(:)          ! Phenological state LAI (m2/m2)

    ! Diagnostic States
    real(kind=4), intent(inout) :: fpar(:)         ! Phenological state FPAR (-)
    real(kind=4), intent(inout) :: gsi(:)          ! Growing Season Index (-)
    real(kind=4), intent(inout) :: temp_fac(:)     ! temperature stress factor (forcing and chilling) (0..1)
    real(kind=4), intent(inout) :: moist_fac(:)    ! moisture stress factor(0..1)
    real(kind=4), intent(inout) :: light_fac(:)    ! light stress factor (0..1)
    real(kind=4), intent(inout) :: force_fac(:)    ! temperature forcing stress factor(0..1)
    real(kind=4), intent(inout) :: chill_fac(:)    ! temperature chilling stress factor (0..1)

    ! Multi-Day Averaged Forcing
    real(kind=4), intent(in) :: temp_ave(:)        ! average minimum daily temperature (K)
    real(kind=4), intent(in) :: moist_ave(:)       ! average mean daily VPD (mb) or cumulative rainfall (mm)
    real(kind=4), intent(in) :: light_ave(:)       ! average cumulative daily photoperiod (hours) or global radiation (W/m2)
    real(kind=4), intent(in) :: chill_ave(:)       ! average cumulative chilling state (K*days)

    ! Climate Control Parameters
    real(kind=4), intent(in) :: temp_min(:)       ! minimum of minimum daily temperature (K) parameter
    real(kind=4), intent(in) :: temp_max(:)       ! maximum of minimum daily temperature (K) parameter
    real(kind=4), intent(in) :: moist_min(:)      ! minimum of daily mean VPD (mb) or cumulative precipitation (mm) parameter
    real(kind=4), intent(in) :: moist_max(:)      ! maximum of daily mean VPD (mb) or cumulative precipitation (mm) parameter
    real(kind=4), intent(in) :: light_min(:)      ! minimum of daily cumulative photoperiod (hours) or global radiation (W/m2) parameter 
    real(kind=4), intent(in) :: light_max(:)      ! maximum of daily cumulative photoperiod (hours) or global radiation (W/m2) parameter      
    real(kind=4), intent(in) :: chill_min(:)      ! minimum of daily cumulative chilling state (K*days) parameter 
    real(kind=4), intent(in) :: chill_max(:)      ! maximum of daily cumulative chilling state (K*days) parameter            

    ! Structural Parameters
    real(kind=4), intent(in) :: lai_min(:)        ! minimum phenological state parameter (m2/m2)
    real(kind=4), intent(in) :: lai_max(:)        ! maximum phenological state parameter (m2/m2)
    real(kind=4), intent(in) :: growthrate(:)     ! growth rate parameter (1/day)
    real(kind=4), intent(in) :: decayrate(:)      ! decay rate parameter (1/day)
    real(kind=4), intent(in) :: lai_sat(:)        ! maximum LAI at saturated FPAR parameter (m2/m2)
    real(kind=4), intent(in) :: fpar_sat(:)       ! asymptotic FPAR for LAI calculation parameter (-)
    real(kind=4), intent(in) :: fvcover(:)        ! fraction of vegetation cover parameter (-)

    ! Flags
    logical, intent(in) :: use_forcing               ! use forcing (temperature-based)
    logical, intent(in) :: use_chilling              ! use chilling (temperature-based)
    logical, intent(in) :: use_moisture              ! use moisture limiting (vapor pressure deficit-based)
    logical, intent(in) :: use_light                 ! use light limiting (photoperiod- or global radiation-based)
    logical, intent(in) :: use_rainfall              ! use rainfall instead of vapor pressure deficit for moisture limitation
    logical, intent(in) :: use_photo           ! use photoperiod instead of global radiation for light limitation

    ! Local Variables

    ! Array size
    integer :: npt                       ! number of forcing/state points
    integer :: npt2                      ! number of parameter points

    ! Diagnostics
    real(kind=4), allocatable :: dgsi(:)

    real(kind=4), allocatable :: dlai(:)
    real(kind=4), allocatable :: rlai(:)

    ! diagnose state and parameter array sizes
    npt = size(lai)
    npt2 = size(temp_min)

    if ((npt.ne.npt2).and.(npt2.ne.1)) then
       write(*,'(A)') "States/Forcing and Parameters need either to be of same dimension (one parameter set per grid point)"
       write(*,'(A)') "Or Parameters have to be scalar (same parameters for all grid points)"
       stop
    endif

    ! allocate local arrays
    allocate(dgsi(npt))
    allocate(dlai(npt))
    allocate(rlai(npt))

    ! Climate control factors calculated with multi-day averaged forcing: 
    ! Temperature, Photoperiod, VPD (or, alternatively Rainfall and Global Radiation)
    if (npt.ne.npt2) then
       if (use_forcing) then
          force_fac=(temp_ave-temp_min(1))/max((temp_max(1)-temp_min(1)),1.e-3)
       else
          force_fac=1.0
       endif
       if (use_chilling) then
          chill_fac=(chill_ave-chill_min(1))/max((chill_max(1)-chill_min(1)),1.e-3)
       else
          chill_fac=1.0
       endif
       if (use_moisture) then
          if (use_rainfall) then
             moist_fac=(moist_ave-moist_min(1))/max((moist_max(1)-moist_min(1)),1.e-3)
          else
             moist_fac=1.-(moist_ave-moist_min(1))/max((moist_max(1)-moist_min(1)),1.e-3)
          endif
       else
          moist_fac=1.0
       endif
       if (use_light) then
          if (use_photo) then
             light_fac=(light_ave-light_min(1))/max((light_max(1)-light_min(1)),1.e-3)
          else
             light_fac=(light_ave-light_min(1))/max((light_max(1)-light_min(1)),1.e-3)
          endif
       else
          light_fac=1.0
       endif
    else
       if (use_forcing) then
          force_fac=(temp_ave-temp_min)/max((temp_max-temp_min),1.e-3)
       else
          force_fac=1.0
       endif
       if (use_chilling) then
          chill_fac=(chill_ave-chill_min)/max((chill_max-chill_min),1.e-3)
       else
          chill_fac=1.0
       endif
       if (use_moisture) then
          if (use_rainfall) then
             moist_fac=(moist_ave-moist_min)/max((moist_max-moist_min),1.e-3)
          else
             moist_fac=1.-(moist_ave-moist_min)/max((moist_max-moist_min),1.e-3)
          endif
       else
          moist_fac=1.0
       endif
       if (use_light) then
          if (use_photo) then
             light_fac=(light_ave-light_min)/max((light_max-light_min),1.e-3)
          else
             light_fac=(light_ave-light_min)/max((light_max-light_min),1.e-3)
          endif
       else
          light_fac=1.0
       endif
    endif

    ! unite temperature forcing and chilling additively
    temp_fac = force_fac - min(1.0,max(0.0,(1.0 - chill_fac)))

    ! bind limitation factors to the range 0..1
    temp_fac=min(1.0,max(0.0,temp_fac))
    moist_fac=min(1.0,max(0.0,moist_fac))
    light_fac=min(1.0,max(0.0,light_fac))
    force_fac=min(1.0,max(0.0,force_fac))
    chill_fac=min(1.0,max(0.0,chill_fac))

    ! 1. calculate new GSI
    gsi = temp_fac * moist_fac * light_fac

    ! 2. calculate prognostic LAI

    ! calculate dLAI/dt
    if (npt.ne.npt2) then
       rlai = (lai-lai_min(1))/max((lai_max(1)-lai_min(1)),1.e-3)
    else
       rlai = (lai-lai_min)/max((lai_max-lai_min),1.e-3)
    endif

    dgsi= min(max(gsi - rlai,-1.),1.)

    if (npt.ne.npt2) then
       where (dgsi.gt.0.)
          dlai = min(100.,min(5.,max(5.e-1,growthrate(1)))) * dgsi * min(max(rlai*(1. - rlai),0.05),0.95)
       elsewhere
          dlai = min(100.,min(5.,max(5.e-1,decayrate(1))))  * dgsi * min(max(rlai*(1. - rlai),0.05),0.95)
       end where
    else
       where (dgsi.gt.0.)
          dlai = min(100.,min(5.,max(5.e-1,growthrate))) * dgsi * min(max(rlai*(1. - rlai),0.05),0.95)
       elsewhere
          dlai = min(100.,min(5.,max(5.e-1,decayrate)))  * dgsi * min(max(rlai*(1. - rlai),0.05),0.95)
       end where
    endif

    ! update prognostic LAI
    lai = lai + dlai

    ! limit LAI
    if (npt.ne.npt2) then
       lai=min(max(lai_min(1),lai),lai_max(1))
    else
       lai=min(max(lai_min,lai),lai_max)
    endif

    ! 3. calculated diagnostic FPAR from LAI

    ! update FPAR
    if (npt.ne.npt2) then
       fpar =  max(fvcover(1),1.e-3) * (1. - exp ( max(min(lai, max(lai_sat(1)*fvcover(1), 1.e-3 )),0.) * &
            log(max(1.-real(fpar_sat(1)),1.e-3)) / max( lai_sat(1) * fvcover(1), 1.e-3 ) ) )
    else
       fpar =  max(fvcover,1.e-3) * (1. - exp ( max(min(lai, max(lai_sat*fvcover, 1.e-3 )),0.) * &
            log(max(1.-real(fpar_sat),1.e-3)) / max( lai_sat * fvcover, 1.e-3 ) ) )
    endif

    ! deallocate local arrays
    deallocate(dgsi)
    deallocate(dlai)
    deallocate(rlai)

  end subroutine phenology_prognostic_gsi

subroutine phenology_diagnostic_gsi (fpar, lai, &
     gsi, temp_fac, moist_fac, light_fac, &
     temp, vpd, photo, dtphen)

  implicit none   
  
  ! Arguments

  ! Diagnostic States
  real(kind=4), intent(inout) :: fpar(:)
  real(kind=4), intent(inout) :: lai(:)
  real(kind=4), intent(inout) :: gsi(:)
  real(kind=4), intent(inout) :: temp_fac(:)
  real(kind=4), intent(inout) :: moist_fac(:)
  real(kind=4), intent(inout) :: light_fac(:)

  ! Multi-Day Averaged Forcing
  real(kind=4), intent(in) :: temp(:)
  real(kind=4), intent(in) :: vpd(:)
  real(kind=4), intent(in) :: photo(:)

  ! Time Step Parameters
  integer, intent(in) :: dtphen ! prediction step [s]

  ! Local Variables

  ! Diagnostics
  real(kind=4) :: weight
  
  ! Climate Control Parameters
  real(kind=4) :: tempmin
  real(kind=4) :: tempmax
  real(kind=4) :: vpdmin
  real(kind=4) :: vpdmax
  real(kind=4) :: photomin
  real(kind=4) :: photomax

  ! Structural Parameters
  real(kind=4) :: laimin
  real(kind=4) :: laimax
  real(kind=4) :: laisat
  real(kind=4) :: fparsat
  real(kind=4) :: fvcover

  ! evaluate number of grid points

  ! old (non-prognostic) GSI model after Jolly et al. 2005
  
  ! 21 day moving average GSI (original formulation)
  weight = exp(-real(dtphen) / 86400. / 21.0) 
  
  ! pre-defined environmental controlling factors (after Jolly et al. 2005)
  ! structural parameters defined after common values valid in LSM's and from MODIS
  laimin = 0.5
  laimax = 6.5
  laisat = 6.5
  fvcover = 0.9     
  fparsat = 0.98
  
  tempmax = 278.16 
  tempmin = 271.16 
  vpdmax = 9.
  vpdmin = 41.
  photomax = 39.6d3/3600.
  photomin = 36.d3/3600.

  ! control factors calculated with daily mean values of T, VPD, and Photoperiod
  temp_fac=(temp-tempmin)/(tempmax-tempmin)
  temp_fac=min(1.,max(0.,temp_fac))
  moist_fac=1.-(vpd-vpdmin)/(vpdmax-vpdmin)
  moist_fac=min(1.,max(0.,moist_fac))
  light_fac=(photo-photomin)/(photomax-photomin)
  light_fac=min(1.,max(0.,light_fac))
  
  ! then calculate 21 day running mean GSI state
  gsi = weight * gsi + (1.-weight) * temp_fac * moist_fac * light_fac
  
  lai = gsi * (laimax - laimin) + laimin
  
  fpar = fvcover * ( 1.-exp( lai * log(1.-fparsat) / (laisat * fvcover) ) )
  
end subroutine phenology_diagnostic_gsi

end module phenology_model_module
