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

module cosmo_module

  ! 2009/05/04 Reto Stockli (Blue Marble Research)

  ! This module handles COSMO model meteorological driver data read
  ! The COSMO files have been processed from GRIB to monthly NetCDF by Jean Marie Bettems by use of fieldextra

  ! Dependencies:
  ! driver_module.F90, bilinear_module.F90,
  ! a working NetCDF library with F90 interface

  ! updates:
  ! none

  ! elevation of COSMO data in output projection
  real(kind=4), allocatable :: hgt_save(:,:)
  logical, allocatable :: hgt_valid(:,:)

  ! COSMO geographic boundaries arrays
  integer :: xmin_cosmo, xmax_cosmo, ymin_cosmo, ymax_cosmo
  integer :: nx_cosmo, ny_cosmo

  ! COSMO record / time variables
  integer :: ntime_cosmo
  integer :: rec1_cosmo
  integer :: nrec_cosmo
  integer :: nvar_cosmo 

  ! COSMO projection information
  character(len=100) :: type_cosmo
  character(len=100) :: optargs_cosmo
  real(kind=4) :: pollon_cosmo
  real(kind=4) :: pollat_cosmo

  ! lon/lat array on cosmo grid
  real(kind=4), allocatable :: lons_cosmo(:), lats_cosmo(:)
  real(kind=4) :: dlon_cosmo, dlat_cosmo

  ! NetCDF variables
  integer :: ncid_cosmo
  character(len=100) :: filename_cosmo
  character(len=9), allocatable :: ncname_cosmo(:)
  logical, allocatable :: select_cosmo(:)  ! select variable for reading
 
  character(len=100) :: directory_cosmo 

contains

  subroutine init_cosmo(meteorologydir, dtmet)

    use time_module
    use exit_module
    
    use netcdf
    use typeSizes

    implicit none

    ! arguments
    character(len=*), intent(in) :: meteorologydir
    integer, intent(out) :: dtmet

    ! local variables

    ! cosmo met driver time step
    dtmet = 3600

    ! initialize COSMO-specific parameters

    nvar_cosmo = 9
    allocate(ncname_cosmo(nvar_cosmo))
    allocate(select_cosmo(nvar_cosmo))

    ncname_cosmo   = (/'T_2M     ','PMSL     ','TD_2M    ','U_10M    ','V_10M    ','ASWDIR_S ','ASWDIFD_S','ALWD_S   ','TOT_PREC '/)
    select_cosmo   = (/.TRUE.,.TRUE.,.TRUE. ,.TRUE.,.TRUE.,.TRUE.,.TRUE. ,.TRUE.,.TRUE./)

    directory_cosmo = trim(meteorologydir)
    filename_cosmo = 'cosmo7_pheno_analysis.'

    ncid_cosmo = -1

  end subroutine init_cosmo

  subroutine read_cosmo(type, xmin, xmax, ymin, ymax, dx, dy, optargs, missing, &
       dtmod, dtmet, ta, rh, ws, sdr, ldr, ps, lsp, cup, mr, vpd, hgt)

    ! This subroutines reads met data from the COSMO model NetCDF data 

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
    character(len=10) :: sdate ! YYYYMMDD00
    integer :: f, k
    integer :: status, irec
    integer :: dimid, varid
    character(len=1) :: string

    real(kind=4) :: lonmin, lonmax, latmin, latmax
    real(kind=4) :: lonmin_cosmo, lonmax_cosmo, latmin_cosmo, latmax_cosmo
    integer :: nlon_cosmo, nlat_cosmo
    
    real(kind=4) :: ceps

    ! physical variables in analysis grid (mean + cumulated)
    real(kind=4), allocatable ::  data_cosmo(:,:)

    ! physical variables in output (model) grid
    real(kind=4), allocatable ::  data(:,:)

    ! flags
    logical, parameter :: verbose = .false.

    ! parameters
    real(kind=4),parameter :: gph = 9.806555 ! conversion factor gph -> m (ecmwf)

    ceps = 1.e-4

    ! close the last month's driver files
    if (newday.and.(ncid_cosmo.ge.0)) then
       if (verbose) write(*,'(A)') 'Closing COSMO meteo driver file.'
       status = nf90_close (ncid_cosmo)
       if (status /= nf90_noerr) call handle_err (status,routine="cosmo_module: error closing NetCDF file")
       ncid_cosmo = -1
    endif

    ! open new driver files  
    if (newday.or.(timestep.eq.1)) then
       !write(*,*) 'Record at Init: ',nt,irec
       write(unit=sdate,fmt='(i10)') 1000000*year + 10000*month + 100*day
       if (verbose) write(*,'(A,A)') 'Opening COSMO meteo driver file for: ',sdate
       status = nf90_open(trim(directory_cosmo)//'/'//trim(filename_cosmo)//sdate//'.nc', &
            nf90_nowrite,ncid_cosmo )
       if (status /= nf90_noerr) call handle_err(status,routine='cosmo_module: error opening '// &
            trim(directory_cosmo)//'/'//trim(filename_cosmo)//sdate//'.nc NetCDF file')

       ! fetch time axis
       status = nf90_inq_dimid ( ncid_cosmo, 'time', dimid )
       if (status /= nf90_noerr) call handle_err(status,routine='cosmo_module: error reading dimension id: time')
       status = nf90_inquire_dimension ( ncid_cosmo, dimid, string, ntime_cosmo )
       if (status /= nf90_noerr) call handle_err(status,routine='cosmo_module: error reading dimension size: time') 
    endif

    ! first read: get COSMO boundaries and read COSMO topography
    if (.not.allocated(lons_cosmo)) then

       ! get COSMO projection information
       status = nf90_inq_varid(ncid_cosmo, "grid_mapping_1", varid)
       if(status /= nf90_NoErr) call handle_err(status)
       status = nf90_get_att(ncid_cosmo, varid, "grid_mapping_name", type_cosmo)
       if(status /= nf90_NoErr) call handle_err(status)   
       status = nf90_get_att(ncid_cosmo, varid, "grid_north_pole_latitude", pollat_cosmo)
       if(status /= nf90_NoErr) call handle_err(status)   
       status = nf90_get_att(ncid_cosmo, varid, "grid_north_pole_longitude", pollon_cosmo)
       if(status /= nf90_NoErr) call handle_err(status)   

       if (pollon_cosmo.gt.180.0) pollon_cosmo = pollon_cosmo - 360.0
       if (pollon_cosmo.lt.-180.0) pollat_cosmo = pollon_cosmo + 360.0

       write(optargs_cosmo,'(A,F0.1,A,F0.1)') "grid_north_pole_longitude=",pollon_cosmo,";grid_north_pole_latitude=",pollat_cosmo

       ! check how far the requested area extends in the COSMO domain
       call getxybounds(type,xmin,xmax,ymin,ymax,dx,dy,optargs,type_cosmo,lonmin,lonmax,latmin,latmax,optargs_cosmo,missing)

       ! get COSMO dimensions
       status = nf90_inq_dimid ( ncid_cosmo, 'x_1', dimid )
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_inquire_dimension ( ncid_cosmo, dimid, string, nlon_cosmo )
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_inq_dimid ( ncid_cosmo, 'y_1', dimid )
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_inquire_dimension ( ncid_cosmo, dimid, string, nlat_cosmo )
       if (status /= nf90_noerr) call handle_err(status)     

       ! COSMO lat/lon bounds (rotated pole)
       allocate(lons_cosmo(nlon_cosmo))
       allocate(lats_cosmo(nlat_cosmo))
       status = nf90_inq_varid(ncid_cosmo, "x_1", varid)
       if(status /= nf90_NoErr) call handle_err(status)
       status = nf90_get_var(ncid_cosmo, varid, lons_cosmo )
       if(status /= nf90_NoErr) call handle_err(status)    
       status = nf90_inq_varid(ncid_cosmo, "y_1", varid)
       if(status /= nf90_NoErr) call handle_err(status)
       status = nf90_get_var(ncid_cosmo, varid, lats_cosmo )
       if(status /= nf90_NoErr) call handle_err(status)   

       ! assume COSMO is a regular grid
       dlon_cosmo = lons_cosmo(nlon_cosmo/2+1) - lons_cosmo(nlon_cosmo/2)
       dlat_cosmo = lats_cosmo(nlat_cosmo/2+1) - lats_cosmo(nlat_cosmo/2)
       dlon_cosmo = 3./real(nint(3./dlon_cosmo),kind=4)
       dlat_cosmo = 3./real(nint(3./dlon_cosmo),kind=4)

       lonmin_cosmo = minval(lons_cosmo) - 0.5*dlon_cosmo
       lonmax_cosmo = maxval(lons_cosmo) + 0.5*dlon_cosmo
       latmin_cosmo = minval(lats_cosmo) - 0.5*dlat_cosmo
       latmax_cosmo = maxval(lats_cosmo) + 0.5*dlat_cosmo

       xmin_cosmo = -1
       xmax_cosmo = -1
       ymin_cosmo = -1
       ymax_cosmo = -1

       ! the multiplication and nint() is a fix for numerical imprecision in comparisons
       ! do not extrapolate over driver data bounds, only interpolate within bounds
       ! driver data dimension 1 goes west -> east
       ! driver data dimension 2 goes south -> north
       do k=nlon_cosmo,1,-1
          if ((nint((lons_cosmo(k)-0.5*dlon_cosmo)/ceps).le.nint(lonmin/ceps)).and.(xmin_cosmo.eq.-1)) xmin_cosmo=k
       enddo
       do k=nlat_cosmo,1,-1
          if ((nint((lats_cosmo(k)-0.5*dlat_cosmo)/ceps).le.nint(latmin/ceps)).and.(ymin_cosmo.eq.-1)) ymin_cosmo=k
       enddo
       do k=1,nlon_cosmo
          if ((nint((lons_cosmo(k)+0.5*dlon_cosmo)/ceps).ge.nint(lonmax/ceps)).and.(xmax_cosmo.eq.-1)) xmax_cosmo=k
       enddo
       do k=1,nlat_cosmo
          if ((nint((lats_cosmo(k)+0.5*dlat_cosmo)/ceps).ge.nint(latmax/ceps)).and.(ymax_cosmo.eq.-1)) ymax_cosmo=k
       enddo

       if ((xmin_cosmo.eq.-1).or.(ymin_cosmo.eq.-1).or.(xmax_cosmo.eq.-1).or.(ymax_cosmo.eq.-1)) then
          write(*,'(A)') 'COSMO file does not cover the requested area: '
          write(*,'(A,4F10.5)') 'We need: ',&
               real(nint(lonmin/ceps),kind=4)*ceps, &
               real(nint(lonmax/ceps),kind=4)*ceps, &
               real(nint(latmin/ceps),kind=4)*ceps, &
               real(nint(latmax/ceps),kind=4)*ceps
          write(*,'(A,4F10.5)') 'but get: ',&
               real(nint(lonmin_COSMO/ceps),kind=4)*ceps, &
               real(nint(lonmax_COSMO/ceps),kind=4)*ceps, &
               real(nint(latmin_COSMO/ceps),kind=4)*ceps, &
               real(nint(latmax_COSMO/ceps),kind=4)*ceps
          write(*,'(A,4I8)') 'COSMO File array bounds: ',xmin_COSMO, xmax_COSMO, ymin_COSMO, ymax_COSMO
          stop
       end if

       ! extend area for correct bilinear interpolation
       ! TODO: check for what to do at outer edges!
       xmin_cosmo = max(xmin_cosmo - 1,1)
       xmax_cosmo = min(xmax_cosmo + 1,nlon_cosmo)
       ymin_cosmo = max(ymin_cosmo - 1,1)
       ymax_cosmo = min(ymax_cosmo + 1,nlat_cosmo)

!       write(*,'(A,4I8)') 'COSMO File array bounds: ',xmin_COSMO, xmax_COSMO, ymin_COSMO, ymax_COSMO

       nx_cosmo = xmax_cosmo - xmin_cosmo + 1
       ny_cosmo = ymax_cosmo - ymin_cosmo + 1

       allocate(data_cosmo(nx_cosmo,ny_cosmo))

       status = nf90_inq_varid(ncid_cosmo, "HSURF", varid)
       if(status /= nf90_NoErr) call handle_err(status)
       status = nf90_get_var(ncid_cosmo, varid, data_cosmo, start = (/xmin_cosmo,ymin_cosmo,1/), count = (/nx_cosmo,ny_cosmo,1/) )
       if(status /= nf90_NoErr) call handle_err(status)  

       call reproject(type_cosmo,lons_cosmo(xmin_cosmo)-0.5*dlon_cosmo,lons_cosmo(xmax_cosmo)+0.5*dlon_cosmo, &
            lats_cosmo(ymin_cosmo)-0.5*dlat_cosmo,lats_cosmo(ymax_cosmo)+0.5*dlat_cosmo,dlon_cosmo,dlat_cosmo,optargs_cosmo, &
            type,xmin,xmax,ymin,ymax,dx,dy,optargs,missing,purge=.false.,idataf4=data_cosmo,odataf4=hgt_save, &
            missingf4=missing)

       deallocate(data_cosmo)

       ! make missing data mask
       allocate(hgt_valid(size(hgt_save,1),size(hgt_save,2)))
       hgt_valid = .false.

       ! geopotential hgt -> metric hgt
       where ((hgt_save-missing).gt.eps)
          hgt_valid = .true.
       end where

    end if ! perform subsetting and topo reading only once

    hgt = hgt_save

    allocate(data_cosmo(nx_cosmo,ny_cosmo))

    ! record pointer
    irec = (timestep-timestep_daystart) * dtmod / dtmet + 1

    ! read data 
    ws = 0.
    sdr = 0.

    do f = 1,nvar_cosmo 

       if (select_cosmo(f)) then

          status = nf90_inq_varid(ncid_cosmo, trim(ncname_cosmo(f)), varid)
          if(status /= nf90_NoErr) call handle_err(status)

          status = nf90_get_var(ncid_cosmo, varid, data_cosmo,  &
               start = (/xmin_cosmo,ymin_cosmo,irec/), count = (/nx_cosmo,ny_cosmo,1/) )           
          if(status /= nf90_NoErr) call handle_err(status)  

          call reproject(type_cosmo,lons_cosmo(xmin_cosmo)-0.5*dlon_cosmo,lons_cosmo(xmax_cosmo)+0.5*dlon_cosmo, &
               lats_cosmo(ymin_cosmo)-0.5*dlat_cosmo,lats_cosmo(ymax_cosmo)+0.5*dlat_cosmo,dlon_cosmo,dlat_cosmo,optargs_cosmo, &
               type,xmin,xmax,ymin,ymax,dx,dy,optargs,missing,purge=.false.,idataf4=data_cosmo,odataf4=data, &
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
                sdr = sdr + data      !  Direct Irradiance [W/m2]
             elsewhere
                sdr = missing
             end where
          case (7)
             where (hgt_valid)
                sdr = sdr + data  ! Diffuse irradiance[mm] over dtmet
             elsewhere
                sdr = missing
             end where
          case (8)
             where (hgt_valid)
                ldr = data      ! [W/m2] 
             elsewhere
                ldr = missing
             end where
          case (9)
             where (hgt_valid)
                lsp = data  ! cumulated [mm] over dtmet
                cup = 0.
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

    !    write(*,*) 'COSMO DATA: ',year,month,day, hour, irec, ntime_cosmo
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

    deallocate(data_cosmo)

  end subroutine read_cosmo

end module cosmo_module

