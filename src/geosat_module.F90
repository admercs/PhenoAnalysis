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

module geosat_module
  
  private
  
  public :: geosat2ll
  public :: ll2geosat
  
contains
  
  ! handles forward coordinate transformation from lon/lat to geostationary satellite column/line
  subroutine ll2geosat(lon,lat,x,y,satname,sensor,channel,satlon,satlat,sathgt,missing)
    
    use msg_seviri_module
    use mfg_mviri_module

    implicit none

    ! arguments
    real(kind=4), intent(in)        :: lon(:,:), lat(:,:)
    real(kind=4), intent(out)       :: x(:,:), y(:,:)
    character(len=*), intent(in)    :: satname
    character(len=*), intent(in)    :: sensor
    character(len=*), intent(in)    :: channel
    real(kind=4), intent(in)        :: satlon
    real(kind=4), intent(in)        :: satlat
    real(kind=4), intent(in)        :: sathgt
    real(kind=4), intent(in)        :: missing
    
    ! local
    integer :: i,j,ni,nj
    logical :: ischannel

    ni = size(lon,1)
    nj = size(lat,2)

    select case (trim(satname))
    case ("mfg")
       select case (trim(sensor))
       case ("mviri")
          if (trim(channel) == "wvir") then
             ischannel = .true.
          else
             ischannel = .false.
          end if
          do j=1,nj
             do i=1,ni
                call ll2mfg_mviri(lon(i,j),lat(i,j),x(i,j),y(i,j),satlon,satlat,sathgt,missing,ischannel)
             end do
          end do
       case default
          write(*,'(A,A,A)') "Geographic Projection for satellite sensor ",trim(sensor)," not implemented."
          x = missing
          y = missing
       end select
    case ("msg")
       select case (trim(sensor))
       case ("seviri")
          if (trim(channel) == "hrvis") then
             ischannel = .true.
          else
             ischannel = .false.
          end if
          do j=1,nj
             do i=1,ni
                call ll2msg_seviri(lon(i,j),lat(i,j),x(i,j),y(i,j),satlon,satlat,sathgt,missing,ischannel)
             end do
          end do
       case default
          write(*,'(A,A,A)') "Geographic Projection for satellite sensor ",trim(sensor)," not implemented."
          x = missing
          y = missing
       end select
    case default
       write(*,'(A,A,A)') "Geographic Projection for satellite name ",trim(satname)," not implemented."
       x = missing
       y = missing
    end select

  end subroutine ll2geosat

  ! handles inverse coordinate transformation from geostationary satellite column/line to lon/lat
  subroutine geosat2ll(x,y,lon,lat,satname,sensor,channel,satlon,satlat,sathgt,missing)

    use msg_seviri_module
    use mfg_mviri_module

    implicit none

    ! arguments
    real(kind=4), intent(in)        :: x(:,:), y(:,:)
    real(kind=4), intent(out)       :: lon(:,:), lat(:,:)
    character(len=*), intent(in)    :: satname
    character(len=*), intent(in)    :: sensor
    character(len=*), intent(in)    :: channel
    real(kind=4), intent(in)        :: satlon
    real(kind=4), intent(in)        :: satlat
    real(kind=4), intent(in)        :: sathgt
    real(kind=4), intent(in)        :: missing
    
    ! local
    integer :: i,j,ni,nj
    logical :: ischannel

    ni = size(x,1)
    nj = size(x,2)

    select case (trim(satname))
    case ("mfg")
       select case (trim(sensor))
       case ("mviri")
          if (trim(channel) == "wvir") then
             ischannel = .true.
          else
             ischannel = .false.
          end if
          do j=1,nj
             do i=1,ni
                call mfg_mviri2ll(x(i,j),y(i,j),lon(i,j),lat(i,j),satlon,satlat,sathgt,missing,ischannel)
             end do
          end do
       case default
          write(*,'(A,A,A)') "Geographic Projection for satellite sensor ",trim(sensor)," not implemented."
          lon = missing
          lat = missing
       end select
    case ("msg")
       select case (trim(sensor))
       case ("seviri")
          if (trim(channel) == "hrvis") then
             ischannel = .true.
          else
             ischannel = .false.
          end if
          do j=1,nj
             do i=1,ni
                call msg_seviri2ll(x(i,j),y(i,j),lon(i,j),lat(i,j),satlon,satlat,sathgt,missing,ischannel)
             end do
          end do
       case default
          write(*,'(A,A,A)') "Geographic Projection for satellite sensor ",trim(sensor)," not implemented."
          lon = missing
          lat = missing
       end select
    case default
       write(*,'(A,A,A)') "Geographic Projection for satellite name ",trim(satname)," not implemented."
       lon = missing
       lat = missing
    end select

  end subroutine geosat2ll


end module geosat_module
