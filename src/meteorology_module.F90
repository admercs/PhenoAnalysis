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

module meteorology_module

  ! This module reads meteorological data from stations (ALMA compliant NetCDF),
  ! ECMWF ERA Interim and Operational Analysis and COSMO

  character(len=100) :: meteorologytype     ! string of meteorological driver data name
  character(len=100) :: meteorologydir      ! directory where Meteorological driver data is located

  ! Met Forcing arrays
  type :: mettype
     real(kind=4), allocatable :: ta(:,:)
     real(kind=4), allocatable :: rh(:,:)
     real(kind=4), allocatable :: ws(:,:)
     real(kind=4), allocatable :: sdr(:,:)
     real(kind=4), allocatable :: ldr(:,:)
     real(kind=4), allocatable :: ps(:,:)
     real(kind=4), allocatable :: lsp(:,:)
     real(kind=4), allocatable :: cup(:,:)
     real(kind=4), allocatable :: mr(:,:)
     real(kind=4), allocatable :: vpd(:,:)
     real(kind=4), allocatable :: hgt(:,:)
  end type mettype
  type(mettype) :: met

contains

  subroutine init_meteorology

    use driver_module
    use exit_module
    use time_module
    use ecmwf_module
    use cosmo_module
    use alma_module

    implicit none

    if (isprediction) then

       ! allocate Met Forcing arrays
       allocate(met%ta(geographic%nx,geographic%ny))
       allocate(met%rh(geographic%nx,geographic%ny))
       allocate(met%ws(geographic%nx,geographic%ny))
       allocate(met%sdr(geographic%nx,geographic%ny))
       allocate(met%ldr(geographic%nx,geographic%ny))
       allocate(met%ps(geographic%nx,geographic%ny))
       allocate(met%lsp(geographic%nx,geographic%ny))
       allocate(met%cup(geographic%nx,geographic%ny))
       allocate(met%mr(geographic%nx,geographic%ny))
       allocate(met%vpd(geographic%nx,geographic%ny))
       allocate(met%hgt(geographic%nx,geographic%ny))

       ! only open met driver files if we are within integration period
       if ((timestep+timestep_start-1).le.timestep_end) then

          ! initialize met driver
          select case (trim(meteorologytype))
          case ("ecmwf") 
             call init_ecmwf(meteorologydir, dtmet)
          case ("cosmo") 
             call init_cosmo(meteorologydir, dtmet)
          case ("alma") 
             call init_alma(meteorologydir, geographic%name, dtmet)
          case default
             write(*,*)
             write(*,'(A,A)') 'Unknown meteorological driver data: ',trim(meteorologytype)
             call model_exit(-1323)
          end select

       endif

    endif

  end subroutine init_meteorology

  subroutine exit_meteorology
    
    use driver_module
    
    implicit none
    
    if (isprediction) then
       deallocate(met%ta,met%rh,met%ws,met%sdr,met%ps,met%ldr,met%lsp,met%cup,met%mr,met%vpd,met%hgt)
    endif
    
  end subroutine exit_meteorology

  subroutine read_meteorology

    ! read meteorology

    use parameter_module
    use driver_module
    use exit_module
    use time_module
    use ecmwf_module
    use cosmo_module
    use alma_module
    use meteorology_average_module
    use phenology_module

    implicit none

    integer :: h, j, i, hidx

    if (isprediction) then

       if (mod((timestep-1)*dtmod,dtmet).eq.0) then

          select case (trim(meteorologytype))
          case ("ecmwf") ! ecmwf met driver
             call read_ecmwf(geographic%type, geographic%xmin, geographic%xmax,&
                  geographic%ymin, geographic%ymax, geographic%dx, geographic%dy, &
                  geographic%optargs, nodata, dtmod, dtmet, &
                  met%ta, met%rh, met%ws, met%sdr, met%ldr, &
                  met%ps, met%lsp, met%cup, met%mr, met%vpd, &
                  met%hgt)   
          case ("cosmo") ! cosmo met driver
             call read_cosmo(geographic%type, geographic%xmin, geographic%xmax,&
                  geographic%ymin, geographic%ymax, geographic%dx, geographic%dy, &
                  geographic%optargs, nodata, dtmod, dtmet, &
                  met%ta, met%rh, met%ws, met%sdr, met%ldr, &
                  met%ps, met%lsp, met%cup, met%mr, met%vpd, &
                  met%hgt)   
          case ("alma")  ! tower met driver 
             call read_alma(geographic%type, geographic%xmin, geographic%xmax,&
                  geographic%ymin, geographic%ymax, geographic%dx, geographic%dy, &
                  geographic%optargs, nodata, dtmod, dtmet, &
                  met%ta, met%rh, met%ws, met%sdr, met%ldr, &
                  met%ps, met%lsp, met%cup, met%mr, met%vpd, &
                  met%hgt)   
          case default
             write(*,*)
             write(*,'(A,A)') 'Unknown meteorological driver data: ',trim(meteorologytype)
             call model_exit(-1324)
          end select

          ! for local analysis, no topography information is used
          ! and thus meteorological data cannot be downscaled
          if (.not.regional) then
             met%hgt(:,:) = 0.
          end if

          ! calculate daily averages
          do h=1,nselect_hgt
             hidx = select_hgt(h)
             do j=1,geographic%ny
                do i=1,geographic%nx
                   if ((lon(i,j) - nodata).gt.eps) then
                      if (size(hgtlevels).eq.1) then
                         ! single elevation level: downscale to output grid elevation
                         call meteorology_average_daily(newday, dtmet, year, doy, &
                              met%ta(i,j), met%mr(i,j), met%sdr(i,j), met%lsp(i,j), met%ps(i,j), &
                              lat(i,j), met%hgt(i,j), hgt(i,j), &
                              force(:,i,j,h,forceidx%tmin), &
                              force(:,i,j,h,forceidx%tmean), &
                              force(:,i,j,h,forceidx%vpd), &
                              force(:,i,j,h,forceidx%photo), &
                              force(:,i,j,h,forceidx%rg), &
                              force(:,i,j,h,forceidx%rain), &
                              nodata)

                         if (simulator) then
                            call meteorology_average_daily(newday, dtmet, year, doy, &
                                 met%ta(i,j), met%mr(i,j), met%sdr(i,j), met%lsp(i,j), met%ps(i,j), &
                                 lat(i,j), met%hgt(i,j), hgt(i,j), &
                                 force_sim(:,i,j,h,forceidx%tmin), &
                                 force_sim(:,i,j,h,forceidx%tmean), &
                                 force_sim(:,i,j,h,forceidx%vpd), &
                                 force_sim(:,i,j,h,forceidx%photo), &
                                 force_sim(:,i,j,h,forceidx%rg), &
                                 force_sim(:,i,j,h,forceidx%rain), &
                                 nodata)
                         endif
                      else
                         ! multiple elevation levels: downscale to discrete elevation bins
                         ! hgtlevels != met%hgt != hgt
                         call meteorology_average_daily(newday, dtmet, year, doy, &
                              met%ta(i,j), met%mr(i,j), met%sdr(i,j), met%lsp(i,j), met%ps(i,j), &
                              lat(i,j), met%hgt(i,j), hgtlevels(hidx), &
                              force(:,i,j,h,forceidx%tmin), &
                              force(:,i,j,h,forceidx%tmean), &
                              force(:,i,j,h,forceidx%vpd), &
                              force(:,i,j,h,forceidx%photo), &
                              force(:,i,j,h,forceidx%rg), &
                              force(:,i,j,h,forceidx%rain), &
                              nodata)
                         if (simulator) then
                            call meteorology_average_daily(newday, dtmet, year, doy, &
                                 met%ta(i,j), met%mr(i,j), met%sdr(i,j), met%lsp(i,j), met%ps(i,j), &
                                 lat(i,j), met%hgt(i,j), hgtlevels(hidx), &
                                 force_sim(:,i,j,h,forceidx%tmin), &
                                 force_sim(:,i,j,h,forceidx%tmean), &
                                 force_sim(:,i,j,h,forceidx%vpd), &
                                 force_sim(:,i,j,h,forceidx%photo), &
                                 force_sim(:,i,j,h,forceidx%rg), &
                                 force_sim(:,i,j,h,forceidx%rain), &
                                 nodata)
                         end if
                      end if
                   else
                      force(:,i,j,:,forceidx%tmin) = nodata
                      force(:,i,j,:,forceidx%tmean) = nodata
                      force(:,i,j,:,forceidx%vpd) = nodata
                      force(:,i,j,:,forceidx%photo) = nodata
                      force(:,i,j,:,forceidx%rg) = nodata
                      force(:,i,j,:,forceidx%rain) = nodata
                      if (simulator) then
                         force_sim(:,i,j,:,forceidx%tmin) = nodata
                         force_sim(:,i,j,:,forceidx%tmean) = nodata
                         force_sim(:,i,j,:,forceidx%vpd) = nodata
                         force_sim(:,i,j,:,forceidx%photo) = nodata
                         force_sim(:,i,j,:,forceidx%rg) = nodata
                         force_sim(:,i,j,:,forceidx%rain) = nodata
                      end if
                   end if
                end do
             end do
          end do

       endif

    endif

  end subroutine read_meteorology

end module meteorology_module
