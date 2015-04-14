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

module ecmwf_module

  ! 2009/02/20 Reto Stockli (Blue Marble Research)

  ! This module handles ECMWF global analysis and reanalysis data IO
  ! The ECMWF files have been processed from GRIB to NetCDF on the ECMWF archive server by scripts
  ! provided with the data

  ! Dependencies:
  ! driver_module.F90, meteo_module.F90, bilinear_module.F90,
  ! a working NetCDF library with F90 interface

  ! updates:
  ! 2009/02/01 Reto Stockli: Inclusion of ERA-Interim, new formatting for forecast fields:
  ! 4 times a day reading of cumulative fields. 00FC: 06-12, 00FC: 12-18, 12FC: 18-24, 12FC: 24-30
  ! those are shifted by one time step when transforming from GRIB to NetCDF so that we receive 4 fields a day: 
  ! 00-06, 06-12, 12-18, 18-24 etc.

  ! elevation of ECMWF data in output projection
  real(kind=4), allocatable :: hgt_save(:,:)
  logical, allocatable :: hgt_valid(:,:)

  ! ECMWF geographic boundaries
  integer :: xmin_ecmwf, xmax_ecmwf, ymin_ecmwf, ymax_ecmwf
  integer :: nx_ecmwf, ny_ecmwf

  ! ECMWF time / record variables
  integer :: ntime_ecmwf, ntime_fc
  integer :: rec1_ecmwf, rec1_fc(2)
  integer :: nrec_ecmwf, nrec_fc

  ! lon/lat spacing of ECMWF driver
  real(kind=4), allocatable :: lons_ecmwf(:), lats_ecmwf(:)
  real(kind=4) :: dlon_ecmwf, dlat_ecmwf

  integer :: maxfiles_ecmwf 
  integer, allocatable :: ncid_ecmwf(:)
  integer, allocatable :: ncid_fc(:)
  character(len=4), allocatable :: filename_ecmwf(:)
  character(len=4), allocatable :: ncname_ecmwf(:)
  logical, allocatable :: select_ecmwf(:)
  logical, allocatable :: is_forecast(:)

  character(len=100) :: directory_ecmwf 

  ! forecast dataset helpers
  integer :: dtmetcum  ! cumulative time step for forecast fields [s]
  integer :: n_fc ! number of forecast steps per day [1-4]

  ! new forecast variables are stored like this:
  ! fc time: 00:00 has 6 12 18 24 and fc time 12:00 has  18 24 30 36 h cumulated forecast values
  ! we calculate from these: 18-6 = 12h value, 24-12 = 18h value, 30-18 = 24h value, 36-24 = 06h value (next day)
  ! shift all values by 12h, repeat first 12h value in first month

contains

  subroutine init_ecmwf(meteorologydir, dtmet)

    use time_module

    implicit none

    ! arguments
    character(len=*), intent(in) :: meteorologydir
    integer, intent(out) :: dtmet

    ! ERA-INTERIM / Operational Analysis/ Forecast specific parameters

    ! Time step
    dtmet = 21600

    ! Define file structure
    maxfiles_ecmwf = 9
    allocate(ncid_ecmwf(maxfiles_ecmwf))
    allocate(filename_ecmwf(maxfiles_ecmwf))
    allocate(ncname_ecmwf(maxfiles_ecmwf))
    allocate(select_ecmwf(maxfiles_ecmwf))
    allocate(is_forecast(maxfiles_ecmwf))

    ncid_ecmwf(:)  = -1
    filename_ecmwf = (/'2T  ', 'MSL ', '2D  ', '10U ', '10V ', 'SSRD','STRD', 'LSP ','CP  '/)
    ncname_ecmwf   = (/'2T  ', 'MSL ', '2D  ', '10U ', '10V ', 'SSRD','STRD', 'LSP ','CP  '/)
    select_ecmwf   = (/.TRUE., .TRUE., .TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE./)
    is_forecast    = (/.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.TRUE.,.TRUE., .TRUE.,.TRUE./)

    directory_ecmwf = trim(meteorologydir)
    dtmetcum = dtmet
    n_fc = 4
    ncid_ecmwf(:) = -1

  end subroutine init_ecmwf

  subroutine read_ecmwf(type, xmin, xmax, ymin, ymax, dx, dy, optargs, missing, &
       dtmod, dtmet, ta, rh, ws, sdr, ldr, ps, lsp, cup, mr, vpd, hgt)

    ! This subroutines reads met data from the ECMWF ERA-40 and Operational Archive dataset 
    ! attention: we have averaged and cumulated fields in the dataset

    use exit_module
    use time_module
    use proj_module
    use reproject_module
    use parameter_module

    use netcdf
    use typeSizes


    implicit none

    ! arguments
    character(len=*), intent(in) :: type
    real(kind=4), intent(in) :: xmin
    real(kind=4), intent(in) :: xmax
    real(kind=4), intent(in) :: ymin
    real(kind=4), intent(in) :: ymax
    real(kind=4), intent(in) :: dx
    real(kind=4), intent(in) :: dy
    character(len=*), intent(in) :: optargs
    real(kind=4), intent(in) :: missing
    integer, intent(in) :: dtmod                    ! time step of model (s)
    integer, intent(in) :: dtmet                    ! time step of meteorology data (s)
    real(kind=4), intent(inout) :: ta(:,:)       ! air temperature (K)    
    real(kind=4), intent(inout) :: rh(:,:)       ! air humidity (%)    
    real(kind=4), intent(inout) :: ws(:,:)       ! wind speed (m/s)    
    real(kind=4), intent(inout) :: sdr(:,:)      ! global radiation (W/m2)    
    real(kind=4), intent(inout) :: ldr(:,:)      ! long wave downward radiation (W/m2)    
    real(kind=4), intent(inout) :: ps(:,:)       ! long wave downward radiation (W/m2)    
    real(kind=4), intent(inout) :: lsp(:,:)      ! large scale precipitation (mm)    
    real(kind=4), intent(inout) :: cup(:,:)      ! convective precipitation (mm)    
    real(kind=4), intent(inout) :: mr(:,:)       ! mixing ratio (g/kg)    
    real(kind=4), intent(inout) :: vpd(:,:)      ! vapor pressure deficit (mb)    
    real(kind=4), intent(inout) :: hgt(:,:)      ! elevation (m a.s.l.)    

    ! local variables
    character(len=6) :: curmon
    integer :: i,f,k
    integer :: status, irec, irec_fc
    integer :: ncid, dimid, varid
    character(len=1) :: string

    real(kind=4) :: lonmin, lonmax, latmin, latmax
    real(kind=4) :: lonmin_ecmwf, lonmax_ecmwf, latmin_ecmwf, latmax_ecmwf
    integer :: nlon_ecmwf, nlat_ecmwf

    ! physical variables in analysis grid (mean + cumulated)
    real(kind=4), allocatable ::  data_ecmwf(:,:), cumdata_ecmwf(:,:,:)

    ! physical variables in output (model) grid
    real(kind=4), allocatable ::  data(:,:)

    ! flags
    logical, parameter :: verbose = .true.

    ! parameters
    real(kind=4),parameter :: gph = 9.806555 ! conversion factor gph -> m (ecmwf)

    ! close the last month's driver files
    if (newmonth) then
       do i=1,maxfiles_ecmwf
          if (ncid_ecmwf(i).ge.0) then
             status = nf90_close (ncid_ecmwf(i))
             if (status /= nf90_noerr) call handle_err (status,routine="ecmwf_module: error closing NetCDF files")
             ncid_ecmwf(i) = -1
          endif
       end do
    endif

    ! open new driver files  
    if (newmonth.or.(timestep.eq.1)) then
       write(unit=curmon,fmt='(i6)') 100*year + month
       if (verbose) write(*,'(A,A)') 'Opening ECMWF meteo driver files for: ',curmon
       do i=1,maxfiles_ecmwf
          status = nf90_open(trim(directory_ecmwf)//'/'//trim(filename_ecmwf(i))//'.'//curmon//'.nc', &
               nf90_nowrite,ncid_ecmwf(i) )
          if (status /= nf90_noerr) call handle_err(status,routine='ecmwf_module: error opening '// &
               trim(directory_ecmwf)//'/'//trim(filename_ecmwf(i))//'.'//curmon//'.nc NetCDF file')
       end do

       ! analysis, averaged fields
       status = nf90_inq_dimid ( ncid_ecmwf(1), 'time', dimid )
       if (status /= nf90_noerr) call handle_err(status,routine='ecmwf_module: error reading dimension id: time')
       status = nf90_inquire_dimension ( ncid_ecmwf(1), dimid, string, ntime_ecmwf )
       if (status /= nf90_noerr) call handle_err(status,routine='ecmwf_module: error reading dimension size: time')

       ! forecast, cumulated fields
       status = nf90_inq_dimid ( ncid_ecmwf(6), 'time', dimid )
       if (status /= nf90_noerr) call handle_err(status,routine='ecmwf_module: error reading dimension: time')
       status = nf90_inquire_dimension ( ncid_ecmwf(6), dimid, string, ntime_fc )
       if (status /= nf90_noerr) call handle_err(status,routine='ecmwf_module: error reading dimension size: time')

    endif

    ! first read: get ECMWF boundaries and read ECMWF topography
    if (.not.allocated(lons_ecmwf)) then

       ! check how far the requested area extends in the ECMWF domain
       call getxybounds(type,xmin,xmax,ymin,ymax,dx,dy,optargs,"latitude_longitude",lonmin,lonmax,latmin,latmax,"",missing)

       ! get ECMWF dimensions
       status = nf90_inq_dimid ( ncid_ecmwf(1), 'lon', dimid )
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_inquire_dimension ( ncid_ecmwf(1), dimid, string, nlon_ecmwf )
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_inq_dimid ( ncid_ecmwf(1), 'lat', dimid )
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_inquire_dimension ( ncid_ecmwf(1), dimid, string, nlat_ecmwf )
       if (status /= nf90_noerr) call handle_err(status)

       ! ECMWF lat/lon bounds
       allocate(lons_ecmwf(nlon_ecmwf))
       allocate(lats_ecmwf(nlat_ecmwf))
       status = nf90_inq_varid(ncid_ecmwf(1), "lon", varid)
       if(status /= nf90_NoErr) call handle_err(status)
       status = nf90_get_var(ncid_ecmwf(1), varid, lons_ecmwf )
       if(status /= nf90_NoErr) call handle_err(status)    
       status = nf90_inq_varid(ncid_ecmwf(1), "lat", varid)
       if(status /= nf90_NoErr) call handle_err(status)
       status = nf90_get_var(ncid_ecmwf(1), varid, lats_ecmwf )
       if(status /= nf90_NoErr) call handle_err(status)   

       ! assume ECMWF is a regular grid
       dlon_ecmwf = lons_ecmwf(nlon_ecmwf/2+1) - lons_ecmwf(nlon_ecmwf/2)
       dlat_ecmwf = lats_ecmwf(nlat_ecmwf/2+1) - lats_ecmwf(nlat_ecmwf/2)
       dlon_ecmwf = 1./real(nint(1./dlon_ecmwf),kind=4)
       dlat_ecmwf = 1./real(nint(1./dlon_ecmwf),kind=4)

       lonmin_ecmwf = minval(lons_ecmwf) - 0.5*dlon_ecmwf
       lonmax_ecmwf = maxval(lons_ecmwf) + 0.5*dlon_ecmwf
       latmin_ecmwf = minval(lats_ecmwf) - 0.5*dlat_ecmwf
       latmax_ecmwf = maxval(lats_ecmwf) + 0.5*dlat_ecmwf

       xmin_ecmwf = -1
       xmax_ecmwf = -1
       ymin_ecmwf = -1
       ymax_ecmwf = -1

       ! the multiplication and nint() is a fix for numerical imprecision in comparisons
       ! do not extrapolate over driver data bounds, only interpolate within bounds
       ! driver data dimension 1 goes west -> east
       ! driver data dimension 2 goes south -> north
       do k=nlon_ecmwf,1,-1
          if ((nint((lons_ecmwf(k)-0.5*dlon_ecmwf)/eps).le.nint(lonmin/eps)).and.(xmin_ecmwf.eq.-1)) xmin_ecmwf=k
       enddo
       do k=nlat_ecmwf,1,-1
          if ((nint((lats_ecmwf(k)-0.5*dlat_ecmwf)/eps).le.nint(latmin/eps)).and.(ymin_ecmwf.eq.-1)) ymin_ecmwf=k
       enddo
       do k=1,nlon_ecmwf
          if ((nint((lons_ecmwf(k)+0.5*dlon_ecmwf)/eps).ge.nint(lonmax/eps)).and.(xmax_ecmwf.eq.-1)) xmax_ecmwf=k
       enddo
       do k=1,nlat_ecmwf
          if ((nint((lats_ecmwf(k)+0.5*dlat_ecmwf)/eps).ge.nint(latmax/eps)).and.(ymax_ecmwf.eq.-1)) ymax_ecmwf=k
       enddo

       if ((xmin_ecmwf.eq.-1).or.(ymin_ecmwf.eq.-1).or.(xmax_ecmwf.eq.-1).or.(ymax_ecmwf.eq.-1)) then
          write(*,'(A)') 'ECMWF file does not cover the requested area: '
          write(*,'(A,4F10.5)') 'We need: ',&
               real(nint(lonmin/eps),kind=4)*eps, &
               real(nint(lonmax/eps),kind=4)*eps, &
               real(nint(latmin/eps),kind=4)*eps, &
               real(nint(latmax/eps),kind=4)*eps
          write(*,'(A,4F10.5)') 'but get: ',&
               real(nint(lonmin_ECMWF/eps),kind=4)*eps, &
               real(nint(lonmax_ECMWF/eps),kind=4)*eps, &
               real(nint(latmin_ECMWF/eps),kind=4)*eps, &
               real(nint(latmax_ECMWF/eps),kind=4)*eps
          write(*,'(A,4I8)') 'ECMWF File array bounds: ',xmin_ECMWF, xmax_ECMWF, ymin_ECMWF, ymax_ECMWF
          stop
       end if

       ! extend area for correct bilinear interpolation
       ! TODO: check for what to do at outer edges!
       xmin_ecmwf = max(xmin_ecmwf - 1,1)
       xmax_ecmwf = min(xmax_ecmwf + 1,nlon_ecmwf)
       ymin_ecmwf = max(ymin_ecmwf - 1,1)
       ymax_ecmwf = min(ymax_ecmwf + 1,nlat_ecmwf)

!       write(*,'(A,4I8)') 'ECMWF File array bounds: ',xmin_ECMWF, xmax_ECMWF, ymin_ECMWF, ymax_ECMWF

       nx_ecmwf = xmax_ecmwf - xmin_ecmwf + 1
       ny_ecmwf = ymax_ecmwf - ymin_ecmwf + 1

       allocate(data_ecmwf(nx_ecmwf,ny_ecmwf))

       status = nf90_open(trim(directory_ecmwf)//'/Z.nc', nf90_nowrite,ncid )
       if (status /= nf90_noerr) call handle_err(status,routine='ecmwf_module: error opening ' &
            //trim(directory_ecmwf)//'/Z.nc NetCDF file')

       status = nf90_inq_varid(ncid, "Z", varid)
       if(status /= nf90_NoErr) call handle_err(status)

       status = nf90_get_var(ncid, varid, data_ecmwf, start = (/xmin_ecmwf,ymin_ecmwf,1/), count = (/nx_ecmwf,ny_ecmwf,1/) )
       if(status /= nf90_NoErr) call handle_err(status)  

       status = nf90_close(ncid)
       if(status /= nf90_NoErr) call handle_err(status)  

       call reproject("latitude_longitude",lons_ecmwf(xmin_ecmwf)-0.5*dlon_ecmwf,lons_ecmwf(xmax_ecmwf)+0.5*dlon_ecmwf,&
            lats_ecmwf(ymin_ecmwf)-0.5*dlat_ecmwf,lats_ecmwf(ymax_ecmwf)+0.5*dlat_ecmwf,dlon_ecmwf,dlat_ecmwf,"",&
            type,xmin,xmax,ymin,ymax,dx,dy,optargs,missing,purge=.false.,idataf4=data_ecmwf,odataf4=hgt_save, &
            missingf4=missing)

       deallocate(data_ecmwf)

       ! make missing data mask
       allocate(hgt_valid(size(hgt_save,1),size(hgt_save,2)))
       hgt_valid = .false.

       ! geopotential hgt -> metric hgt
       where ((hgt_save-missing).gt.eps)
          hgt_valid = .true.
          hgt_save = hgt_save / gph
       end where

    end if ! perform subsetting and topo reading only once

    hgt = hgt_save

    allocate(data_ecmwf(nx_ecmwf,ny_ecmwf))
    allocate(cumdata_ecmwf(nx_ecmwf,ny_ecmwf,n_fc))

    ! record pointer: we read ecmwf era-40 4 x per day
    irec = (timestep-timestep_monthstart) * dtmod / dtmet
    irec_fc = irec/4*4 + 1
    irec = irec + 1

    ! read data 
    ws = 0.
    lsp = 0.

    do f=1,maxfiles_ecmwf

       print*,trim(ncname_ecmwf(f))

       if (is_forecast(f).and.select_ecmwf(f)) then
          ! cumulated fields:  00 and 12h forecasts (Operational) or 6/18h Forecast composites (ERA-40) (2 per day)

          status = nf90_inq_varid(ncid_ecmwf(f), trim(ncname_ecmwf(f)), varid)
          if(status /= nf90_NoErr) call handle_err(status)

          status = nf90_get_var(ncid_ecmwf(f), varid, cumdata_ecmwf,  &
               start = (/xmin_ecmwf,ymin_ecmwf,irec_fc/), count = (/nx_ecmwf,ny_ecmwf,n_fc/) )           
          if(status /= nf90_NoErr) call handle_err(status)  

          data_ecmwf = (cumdata_ecmwf(:,:,1) + cumdata_ecmwf(:,:,2) + cumdata_ecmwf(:,:,3) + cumdata_ecmwf(:,:,4)) * 0.25

       endif

       if ((.not.is_forecast(f)).and.select_ecmwf(f)) then
          ! instantaneous averaged fields  6, 12, 18, 24 ERA-40/Operational Analyses (4 per day)

          status = nf90_inq_varid(ncid_ecmwf(f), trim(ncname_ecmwf(f)), varid)
          if(status /= nf90_NoErr) call handle_err(status)

          status = nf90_get_var(ncid_ecmwf(f), varid, data_ecmwf,  &
               start = (/xmin_ecmwf,ymin_ecmwf,irec/), count = (/nx_ecmwf,ny_ecmwf,1/) )
          if(status /= nf90_NoErr) call handle_err(status)  

       endif

       if (select_ecmwf(f)) then
          call reproject("latitude_longitude",lons_ecmwf(xmin_ecmwf)-0.5*dlon_ecmwf,lons_ecmwf(xmax_ecmwf)+0.5*dlon_ecmwf,&
               lats_ecmwf(ymin_ecmwf)-0.5*dlat_ecmwf,lats_ecmwf(ymax_ecmwf)+0.5*dlat_ecmwf,dlon_ecmwf,dlat_ecmwf,"",&
               type,xmin,xmax,ymin,ymax,dx,dy,optargs,missing,purge=.false.,idataf4=data_ecmwf,odataf4=data, &
               missingf4=missing)

          select case (f)
          case (1) 
             ta = data  ! air temperature [K]       
          case (2)
             where (hgt_valid)
                ps = data*exp((-1.)*0.029*gph*hgt/(8.314*ta))   ! Pa MSL -> [Pa local pressure]
             elsewhere
                ps = missing
             end where
          case (3)
             where (hgt_valid)
                ! mixing ratio [g/kg]
                mr = 6.11*10**(7.5*(data-273.16)/(237.7+(data-273.16))) * 62200. / &
                     ( ps - 100. * 6.11*10**(7.5*(data-273.16)/(237.7+(data-273.16))) )

                ! vpd [mb]
                vpd =   6.11*10**(7.5*(ta-273.16)/(237.7+(ta-273.16))) - &
                     6.11*10**(7.5*(data-273.16)/(237.7+(data-273.16)))

                ! relative humidity [%]
                rh = 100. * 6.11*10**(7.5*(data-273.16)/(237.7+(data-273.16))) / &
                     (6.11*10**(7.5*(ta-273.16)/(237.7+(ta-273.16))))
             elsewhere
                mr = missing
                vpd = missing
                rh = missing
             end where
          case (4)
             where (hgt_valid)
                ws = ws + data**2  ! add quadratic wind speed in u direction [m/s]
             end where
          case (5)
             where (hgt_valid)
                ws = ws + data**2  ! add quadratic wind speed in v direction [m/s]
             end where
          case (6)
             where (hgt_valid)
                sdr = data / real(dtmetcum,kind=4)      ! cumulated [W/m2*s] over dtmet
             elsewhere
                sdr = missing
             end where
          case (7)
             where (hgt_valid)
                ldr = data / real(dtmetcum,kind=4)      ! cumulated [W/m2*s] over dtmet
             elsewhere
                ldr = missing
             end where
          case (8)
             where (hgt_valid)
                lsp = lsp + data*1.e3*real(dtmet,kind=4)/real(dtmetcum,kind=4)  ! cumulated [mm] over dtmet
             elsewhere
                lsp = missing
             end where
          case (9)
             where (hgt_valid)
                lsp = lsp + data*1.e3*real(dtmet,kind=4)/real(dtmetcum,kind=4)  ! cumulated [mm] over dtmet
                cup = data*1.e3*real(dtmet,kind=4)/real(dtmetcum,kind=4)
             elsewhere
                cup = missing
             end where
          end select

       end if
    end do

    where (hgt_valid)
       ws = sqrt(ws)
    elsewhere
       ws = missing
    end where

    !    write(*,*) 'ECMWF DATA: ',year,month,day, irec,irec_fc, ntime_ecmwf, ntime_fc
    !    write(*,*) 'hgt: ',hgt(:,:)
    !    write(*,*) 'ta:  ',ta(:,:)
    !    write(*,*) 'rh:  ',rh(:,:)
    !    write(*,*) 'ws:  ',ws(:,:)
    !    write(*,*) 'sdr: ',sdr(:,:)
    !    write(*,*) 'ps:  ',ps(:,:)
    !    write(*,*) 'ldr: ',ldr(:,:)
    !    write(*,*) 'lsp: ',lsp(:,:)
    !    write(*,*) 'cup: ',cup(:,:)
    !    write(*,*) 'vpd: ',vpd(:,:)

    deallocate(data_ecmwf)
    deallocate(cumdata_ecmwf)

  end subroutine read_ecmwf

end module ecmwf_module

