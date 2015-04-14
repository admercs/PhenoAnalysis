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

module phenology_module

  ! The diagnostic / prognostic phenology model wrapper for MODIS FPAR/LAI data assimilation.

contains

  subroutine predict_phenology

    use parameter_module
    use driver_module
    use time_module
    use meteorology_average_module
    use phenology_model_module

    implicit none

    integer :: h, p, j, i, pidx
    real(kind=4) :: temp_ave(nrens)
    real(kind=4) :: moist_ave(nrens)
    real(kind=4) :: light_ave(nrens)
    real(kind=4) :: chill_ave(nrens)
    real(kind=4) :: temp_min(nrens)
    real(kind=4) :: temp_max(nrens)
    real(kind=4) :: moist_min(nrens)
    real(kind=4) :: moist_max(nrens)
    real(kind=4) :: light_min(nrens)
    real(kind=4) :: light_max(nrens)
    real(kind=4) :: chill_min(nrens)
    real(kind=4) :: chill_max(nrens)
    real(kind=4) :: lat_ens(nrens)

    real(kind=4) :: temp_ave_sim(1)
    real(kind=4) :: moist_ave_sim(1)
    real(kind=4) :: light_ave_sim(1)
    real(kind=4) :: chill_ave_sim(1)
    real(kind=4) :: temp_min_sim(1)
    real(kind=4) :: temp_max_sim(1)
    real(kind=4) :: moist_min_sim(1)
    real(kind=4) :: moist_max_sim(1)
    real(kind=4) :: light_min_sim(1)
    real(kind=4) :: light_max_sim(1)
    real(kind=4) :: chill_min_sim(1)
    real(kind=4) :: chill_max_sim(1)
    real(kind=4) :: lat_ens_sim(1)

    if (isprediction) then

       do h=1,nselect_hgt
          do p=1,nselect_pft
             pidx = select_pft(p)
             do j=1,geographic%ny
                do i=1,geographic%nx
                   if ((lon(i,j) - nodata).gt.eps) then
                      if (original) then
                         call phenology_diagnostic_gsi( &
                              state(:,i,j,p,h,stateidx%fpar), &
                              state(:,i,j,p,h,stateidx%lai), &
                              state(:,i,j,p,h,stateidx%gsi), &
                              state(:,i,j,p,h,stateidx%temp_fac), &
                              state(:,i,j,p,h,stateidx%moist_fac), &
                              state(:,i,j,p,h,stateidx%light_fac), &
                              force(:,i,j,h,forceidx%tmin), &
                              force(:,i,j,h,forceidx%vpd), &
                              force(:,i,j,h,forceidx%photo), &
                              dtphen)
                      else
                         lat_ens(:) = lat(i,j)

                         call meteorology_average_multiday(&
                              lat_ens, &
                              force(:,i,j,h,forceidx%tmin), &
                              force(:,i,j,h,forceidx%tmean), &
                              force(:,i,j,h,forceidx%vpd), &
                              force(:,i,j,h,forceidx%photo), &
                              force(:,i,j,h,forceidx%rg), &
                              force(:,i,j,h,forceidx%rain), &
                              state(:,i,j,p,h,stateidx%tmin_ave), &
                              state(:,i,j,p,h,stateidx%vpd_ave), &
                              state(:,i,j,p,h,stateidx%photo_ave), &
                              state(:,i,j,p,h,stateidx%rg_ave), &
                              state(:,i,j,p,h,stateidx%rain_ave), &
                              state(:,i,j,p,h,stateidx%chill_ave), &
                              param(:,pidx,paramidx%tmin_tave), &
                              param(:,pidx,paramidx%vpd_tave), &
                              param(:,pidx,paramidx%photo_tave), &
                              param(:,pidx,paramidx%rg_tave), &
                              param(:,pidx,paramidx%rain_tave), &
                              param(:,pidx,paramidx%chill_tave), &
                              dtphen, doy, nodata)

                         temp_ave(:) = state(:,i,j,p,h,stateidx%tmin_ave)
                         temp_min(:) = param(:,pidx,paramidx%tmin_min)
                         temp_max(:) = param(:,pidx,paramidx%tmin_max)
                         if (use_rainfall) then
                            moist_ave(:) = state(:,i,j,p,h,stateidx%rain_ave)
                            moist_min(:) = param(:,pidx,paramidx%rain_min)
                            moist_max(:) = param(:,pidx,paramidx%rain_max)
                         else
                            moist_ave(:) = state(:,i,j,p,h,stateidx%vpd_ave)
                            moist_min(:) = param(:,pidx,paramidx%vpd_min)
                            moist_max(:) = param(:,pidx,paramidx%vpd_max)
                         endif
                         if (use_photo) then
                            light_ave(:) = state(:,i,j,p,h,stateidx%photo_ave)
                            light_min(:) = param(:,pidx,paramidx%photo_min)
                            light_max(:) = param(:,pidx,paramidx%photo_max)
                         else
                            light_ave(:) = state(:,i,j,p,h,stateidx%rg_ave)
                            light_min(:) = param(:,pidx,paramidx%rg_min)
                            light_max(:) = param(:,pidx,paramidx%rg_max)
                         endif
                         chill_ave(:) = state(:,i,j,p,h,stateidx%chill_ave)
                         chill_min(:) = param(:,pidx,paramidx%chill_min)
                         chill_max(:) = param(:,pidx,paramidx%chill_max)

                         call phenology_prognostic_gsi( &
                              state(:,i,j,p,h,stateidx%fpar), &
                              state(:,i,j,p,h,stateidx%lai), &
                              state(:,i,j,p,h,stateidx%gsi), &
                              state(:,i,j,p,h,stateidx%temp_fac), &
                              state(:,i,j,p,h,stateidx%moist_fac), &
                              state(:,i,j,p,h,stateidx%light_fac), &
                              state(:,i,j,p,h,stateidx%force_fac), &
                              state(:,i,j,p,h,stateidx%chill_fac), &
                              temp_ave(:), moist_ave(:), light_ave(:), chill_ave(:), &
                              temp_min(:), temp_max(:), &
                              moist_min(:), moist_max(:), &
                              light_min(:), light_max(:), &
                              chill_min(:), chill_max(:), &
                              param(:,pidx,paramidx%lai_min), &
                              param(:,pidx,paramidx%lai_max), &
                              param(:,pidx,paramidx%growthrate), &
                              param(:,pidx,paramidx%decayrate), &
                              param(:,pidx,paramidx%lai_sat), &
                              param(:,pidx,paramidx%fpar_sat), &
                              param(:,pidx,paramidx%fvcover), &
                              use_forcing, use_chilling, use_moisture, use_light, use_rainfall, use_photo) 

                      endif

                      if (simulator) then 
                         if (original) then
                            call phenology_diagnostic_gsi( &
                                 state_sim(:,i,j,p,h,stateidx%fpar), &
                                 state_sim(:,i,j,p,h,stateidx%lai), &
                                 state_sim(:,i,j,p,h,stateidx%gsi), &
                                 state_sim(:,i,j,p,h,stateidx%temp_fac), &
                                 state_sim(:,i,j,p,h,stateidx%moist_fac), &
                                 state_sim(:,i,j,p,h,stateidx%light_fac), &
                                 force_sim(:,i,j,h,forceidx%tmin), &
                                 force_sim(:,i,j,h,forceidx%vpd), &
                                 force_sim(:,i,j,h,forceidx%photo), &
                                 dtphen)
                         else

                            lat_ens_sim(:) = lat(i,j)

                            call meteorology_average_multiday(&
                                 lat_ens_sim, &
                                 force_sim(:,i,j,h,forceidx%tmin), &
                                 force_sim(:,i,j,h,forceidx%tmean), &
                                 force_sim(:,i,j,h,forceidx%vpd), &
                                 force_sim(:,i,j,h,forceidx%photo), &
                                 force_sim(:,i,j,h,forceidx%rg), &
                                 force_sim(:,i,j,h,forceidx%rain), &
                                 state_sim(:,i,j,p,h,stateidx%tmin_ave), &
                                 state_sim(:,i,j,p,h,stateidx%vpd_ave), &
                                 state_sim(:,i,j,p,h,stateidx%photo_ave), &
                                 state_sim(:,i,j,p,h,stateidx%rg_ave), &
                                 state_sim(:,i,j,p,h,stateidx%rain_ave), &
                                 state_sim(:,i,j,p,h,stateidx%chill_ave), &
                                 param_sim(:,pidx,paramidx%tmin_tave), &
                                 param_sim(:,pidx,paramidx%vpd_tave), &
                                 param_sim(:,pidx,paramidx%photo_tave), &
                                 param_sim(:,pidx,paramidx%rg_tave), &
                                 param_sim(:,pidx,paramidx%rain_tave), &
                                 param_sim(:,pidx,paramidx%chill_tave), &
                                 dtphen, doy, nodata)

                            temp_ave_sim(:) = state_sim(:,i,j,p,h,stateidx%tmin_ave)
                            temp_min_sim(:) = param_sim(:,pidx,paramidx%tmin_min)
                            temp_max_sim(:) = param_sim(:,pidx,paramidx%tmin_max)
                            if (use_rainfall) then
                               moist_ave_sim(:) = state_sim(:,i,j,p,h,stateidx%rain_ave)
                               moist_min_sim(:) = param_sim(:,pidx,paramidx%rain_min)
                               moist_max_sim(:) = param_sim(:,pidx,paramidx%rain_max)
                            else
                               moist_ave_sim(:) = state_sim(:,i,j,p,h,stateidx%vpd_ave)
                               moist_min_sim(:) = param_sim(:,pidx,paramidx%vpd_min)
                               moist_max_sim(:) = param_sim(:,pidx,paramidx%vpd_max)
                            endif
                            if (use_photo) then
                               light_ave_sim(:) = state_sim(:,i,j,p,h,stateidx%photo_ave)
                               light_min_sim(:) = param_sim(:,pidx,paramidx%photo_min)
                               light_max_sim(:) = param_sim(:,pidx,paramidx%photo_max)
                            else
                               light_ave_sim(:) = state_sim(:,i,j,p,h,stateidx%rg_ave)
                               light_min_sim(:) = param_sim(:,pidx,paramidx%rg_min)
                               light_max_sim(:) = param_sim(:,pidx,paramidx%rg_max)
                            endif
                            chill_ave_sim(:) = state_sim(:,i,j,p,h,stateidx%chill_ave)
                            chill_min_sim(:) = param_sim(:,pidx,paramidx%chill_min)
                            chill_max_sim(:) = param_sim(:,pidx,paramidx%chill_max)

                            call phenology_prognostic_gsi( &
                                 state_sim(:,i,j,p,h,stateidx%fpar), &
                                 state_sim(:,i,j,p,h,stateidx%lai), &
                                 state_sim(:,i,j,p,h,stateidx%gsi), &
                                 state_sim(:,i,j,p,h,stateidx%temp_fac), &
                                 state_sim(:,i,j,p,h,stateidx%moist_fac), &
                                 state_sim(:,i,j,p,h,stateidx%light_fac), &
                                 state_sim(:,i,j,p,h,stateidx%force_fac), &
                                 state_sim(:,i,j,p,h,stateidx%chill_fac), &
                                 temp_ave_sim(:), moist_ave_sim(:), light_ave_sim(:), chill_ave_sim(:), &
                                 temp_min_sim(:), temp_max_sim(:), &
                                 moist_min_sim(:), moist_max_sim(:), &
                                 light_min_sim(:), light_max_sim(:), &
                                 chill_min_sim(:), chill_max_sim(:), &
                                 param_sim(:,pidx,paramidx%lai_min), &
                                 param_sim(:,pidx,paramidx%lai_max), &
                                 param_sim(:,pidx,paramidx%growthrate), &
                                 param_sim(:,pidx,paramidx%decayrate), &
                                 param_sim(:,pidx,paramidx%lai_sat), &
                                 param_sim(:,pidx,paramidx%fpar_sat), &
                                 param_sim(:,pidx,paramidx%fvcover), &
                                 use_forcing, use_chilling, use_moisture, use_light, use_rainfall, use_photo) 

                         end if
                      end if
                   else 
                      state(:,i,j,:,:,stateidx%fpar) = nodata
                      state(:,i,j,:,:,stateidx%lai) = nodata
                      state(:,i,j,:,:,stateidx%gsi) = nodata
                      state(:,i,j,:,:,stateidx%temp_fac) = nodata
                      state(:,i,j,:,:,stateidx%moist_fac) = nodata
                      state(:,i,j,:,:,stateidx%light_fac) = nodata
                      state(:,i,j,:,:,stateidx%force_fac) = nodata
                      state(:,i,j,:,:,stateidx%chill_fac) = nodata
                      state(:,i,j,:,:,stateidx%tmin_ave) = nodata
                      state(:,i,j,:,:,stateidx%vpd_ave) = nodata
                      state(:,i,j,:,:,stateidx%photo_ave) = nodata
                      state(:,i,j,:,:,stateidx%rg_ave) = nodata
                      state(:,i,j,:,:,stateidx%rain_ave) = nodata
                      state(:,i,j,:,:,stateidx%chill_ave) = nodata
                      if (simulator) then
                         state_sim(:,i,j,:,:,stateidx%fpar) = nodata
                         state_sim(:,i,j,:,:,stateidx%lai) = nodata
                         state_sim(:,i,j,:,:,stateidx%gsi) = nodata
                         state_sim(:,i,j,:,:,stateidx%temp_fac) = nodata
                         state_sim(:,i,j,:,:,stateidx%moist_fac) = nodata
                         state_sim(:,i,j,:,:,stateidx%light_fac) = nodata
                         state_sim(:,i,j,:,:,stateidx%force_fac) = nodata
                         state_sim(:,i,j,:,:,stateidx%chill_fac) = nodata
                         state_sim(:,i,j,:,:,stateidx%tmin_ave) = nodata
                         state_sim(:,i,j,:,:,stateidx%vpd_ave) = nodata
                         state_sim(:,i,j,:,:,stateidx%photo_ave) = nodata
                         state_sim(:,i,j,:,:,stateidx%rg_ave) = nodata
                         state_sim(:,i,j,:,:,stateidx%rain_ave) = nodata
                         state_sim(:,i,j,:,:,stateidx%chill_ave) = nodata
                      end if
                   end if
                end do
             end do
          end do
       end do

    end if

  end subroutine predict_phenology

end module phenology_module
