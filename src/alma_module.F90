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

module alma_module

  character(len=100) :: directory_alma
  character(len=100) :: filename_alma
  integer :: ncid_alma
  integer :: ntime_alma ! number of time steps in a year

contains

  subroutine init_alma(meteorologydir, sitename, dtmet)

    use exit_module
    use time_module

    use netcdf
    use typeSizes

    implicit none

    ! arguments
    character(len=*), intent(in) :: meteorologydir
    character(len=*), intent(in) :: sitename
    integer, intent(out) :: dtmet

    ! local variables
    character(len=4) :: curyear ! string of the year
    integer :: status ! file operation status
    integer :: varid ! netcdf variable id
    character(len=1) :: string ! temporary character array

    directory_alma = trim(meteorologydir)//trim(sitename)//'/'
    filename_alma = trim(sitename)//'_'

    ! open new alma driver file 
    write(unit=curyear,fmt='(i4)') year
    write(*,'(A,A)') 'Opening ALMA driver files for: ',curyear
    status = nf90_open(trim(directory_alma)//trim(filename_alma)//curyear//'.nc', &
         nf90_nowrite,ncid_alma )
    if (status /= nf90_noerr) call handle_err(status,routine='alma_module: error opening '// &
         trim(directory_alma)//trim(filename_alma)//curyear//'.nc NetCDF file')

    ! get time step length
    status = nf90_inq_dimid ( ncid_alma, 'time', varid )
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inquire_dimension ( ncid_alma, varid, string, ntime_alma )
    if (status /= nf90_noerr) call handle_err(status)

    if (isleap) then 
       dtmet = 366*86400 / ntime_alma 
    else
       dtmet = 365*86400 / ntime_alma
    endif

    status = nf90_close(ncid_alma)
    if(status /= nf90_NoErr) call handle_err(status) 
    ncid_alma = -1

  end subroutine init_alma

  subroutine read_alma(type, xmin, xmax, ymin, ymax, dx, dy, optargs, missing, &
       dtmod, dtmet, ta, rh, ws, sdr, ldr, ps, lsp, cup, mr, vpd, hgt)

    ! this routine reads met forcing from a yearly ALMA compatible file

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
 
    ! parameters
    integer, parameter :: nvar_alma = 8
    character(len=6), parameter :: alma_name(nvar_alma) = (/'Tair  ','Qair  ','Wind  ','Rainf ', &
         'PSurf ','SWdown','LWdown','CO2air'/)

    ! local variables
    character(len=4) :: curyear ! string of the year
    integer :: i ! counter
    integer :: irec ! record pointer
    integer :: status ! file operation status
    integer :: varid ! netcdf variable id
    real(kind=4) :: alma_var(nvar_alma), e, esat
    real(kind=4) :: lon_alma,lat_alma
    real(kind=4) :: missing_alma

    write(unit=curyear,fmt='(i4)') year

    ! open yearly alma file
    if (newyear.or.(timestep.eq.1)) then
       ! close previous file if existent

       write(*,'(A,A)') 'Opening ALMA driver files for: ',curyear
       status = nf90_open(trim(directory_alma)//trim(filename_alma)//curyear//'.nc', &
            nf90_nowrite,ncid_alma )
       if (status /= nf90_noerr) call handle_err(status,routine='alma_module: error opening '// &
            trim(directory_alma)//trim(filename_alma)//curyear//'.nc NetCDF file')
    endif

    ! determine site location
    if (timestep.eq.1) then

       if ((trim(type).ne."latitude_longitude").and.(trim(optargs).ne."")) then
          write(*,'(A)') "Please choose latitude_longitude grid when reading ALMA data"
          stop
       end if

       ! assign height/lat/lon
       status = nf90_inq_varid(ncid_alma, 'lon', varid)
       if(status /= nf90_NoErr) call handle_err(status)
       status = nf90_get_var(ncid_alma, varid, lon_alma )
       if(status /= nf90_NoErr) call handle_err(status)  

       status = nf90_inq_varid(ncid_alma, 'lat', varid)
       if(status /= nf90_NoErr) call handle_err(status)
       status = nf90_get_var(ncid_alma, varid, lat_alma )
       if(status /= nf90_NoErr) call handle_err(status)  

       hgt(:,:) = 0.

       if ((lon_alma.lt.(xmin-0.0*dx)).or.(lon_alma.gt.(xmax+0.0*dx))) then
          write(*,'(A,F12.8,F12.8)') 'Station and model longitudes differ: ',lon_alma,0.5*(xmax+xmin)
       endif
       if ((lat_alma.lt.(ymin-0.0*dy)).or.(lat_alma.gt.(ymax+0.0*dy))) then
          write(*,'(A,F12.8,F12.8)') 'Station and model latitudes differ: ',lat_alma,0.5*(xmax+xmin)
       endif
    endif

    ! calculate alma time step
    irec = (timestep-timestep_yearstart)*dtmod/dtmet+1

    ! read the data for a single time step
    do i=1,nvar_alma
       status = nf90_inq_varid(ncid_alma, trim(alma_name(i)), varid)
       if(status /= nf90_NoErr) call handle_err(status)
       status = nf90_get_var(ncid_alma, varid, alma_var(i:i), &
            start = (/1,1,1,irec/), count = (/1,1,1,1/) )
       if(status /= nf90_NoErr) call handle_err(status)  
       status = nf90_get_att(ncid_alma, varid, "missing_value", missing_alma)     
    enddo

    where (alma_var.eq.missing_alma)
       alma_var = missing
    end where

    ! we spatially aggregate met forcing to the ensemble domain
    ! for flux tower, this is a 1:1 operation

    ! derive vapor pressure (mb) from specific humidity (kg/kg)
    e = alma_var(2) * alma_var(5) * 0.01 /(0.62197 + 0.378030 * alma_var(2))

    ! derive saturated vapor pressure (mb) from temperature (K)
    esat = 6.11*10**(7.5*(alma_var(1)-273.16)/(237.7+(alma_var(1)-273.16)))

    ! create mixing ratio (g/kg) from specific humidity (kg/kg)
    mr(:,:) = 1.e3 * alma_var(2) / (1.-alma_var(2))

    ta(:,:)=alma_var(1) ! K
    rh(:,:)=e/esat*100. ! %
    ws(:,:)=alma_var(3) ! m/s
    lsp(:,:)=alma_var(4)*dtmet ! mm 
    ps(:,:)=alma_var(5) ! Pa
    sdr(:,:)=alma_var(6) ! W/m2
    ldr(:,:)=alma_var(7) ! W/m2
    cup(:,:)=0. ! convective precip not in file

    ! create vp (mb) and vp_sat (mb) and calculate vpd (mb)
    vpd(:,:) = 6.11*10**(7.5*(ta-273.16)/(237.7+(ta-273.16))) - &
         mr * ps/(mr+622.)/100.


!    write(*,*) 'ALMA DATA: ',year,month,day, hour, irec, ntime_alma
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

    if (oldyear) then
       write(*,'(A,A)') 'Closing ALMA driver file for: ',curyear       
       status = nf90_close(ncid_alma)
       if(status /= nf90_NoErr) call handle_err(status)  
       ncid_alma = -1
    endif

  end subroutine read_alma

end module alma_module
