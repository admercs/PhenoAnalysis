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

module photoperiod_module

contains

function fff(x)

  ! calculate earth's orbital parameters
  implicit none

  ! arguments
  real(kind=4), intent(in) :: x

  ! local
  real(kind=4), parameter :: perhl  = 102.7
  real(kind=4), parameter :: eccn   = 0.016715
  real(kind=4), parameter :: pi  = 3.141592

  real(kind=4) :: perhlr, reccn, fff

  PERHLR = PERHL * (pi/180.)
  RECCN  = 1./(1.-ECCN*ECCN)**1.5
  FFF = RECCN*(1.-ECCN*COS(X-PERHLR))**2.

end function fff

function photoperiod(npt, year, doy, lat)

  ! calculate photoperiod length (hours) for given day and latitude
  ! daylength is calculated as the time where the sun is higher than
  ! the maximum 85 degrees of the solar zenith angle
  implicit none

  ! arguments
  integer, intent(in) :: npt            ! number of points in latitude/photoperiod vector
  integer, intent(in) :: year           ! Year YYYY
  integer, intent(in) :: doy            ! day of year (1..365 (366))
  real(kind=4), intent(in) :: lat(npt)  ! latitude (-90.0 .. 90.0)

  ! local parameters
  real(kind=4), parameter :: pi  = 3.141592
  real(kind=4), parameter :: eqnx    = 80.
  real(kind=4), parameter :: decmax = 23.441

  integer, parameter :: ntime =  1440        ! number of time steps per day (10 minute time step)

  ! local variables
  real(kind=4) :: photoperiod(npt)   ! photoperiod in [hours]
  logical :: is_leap
  real(kind=4) :: totdays, h, y
  real(kind=4) :: iiday
  real(kind=4) :: t1, t2, t3, t4
  real(kind=4) :: sindd, cosdd
  real(kind=4) :: con, hrangle
  real(kind=4) :: sinll(npt), cosll(npt), zenith(npt)
  integer :: t

  is_leap = ( ((mod(year,4).eq.0).and.(mod(year,100).ne.0)).or.(mod(year,400).eq.0) )

  if (is_leap) then
     totdays = 366.
  else
     totdays = 365.
  endif
  
  h = (pi*2.)/totdays
  y = 0.
  
  iiday  = real(doy) - 1. - eqnx
  if (iiday.ne.0.) then
     if (iiday.lt.0.) iiday = iiday + totdays + 0.1
     do while (iiday.gt.0.)
        iiday = iiday - 1.
        t1 = fff(y)*h
        t2 = fff(y+t1*0.5)*h
        t3 = fff(y+t2*0.5)*h
        t4 = fff(y+t3)*h
        y = y + (t1+2.*(t2+t3)+t4)/6.
     end do
  endif
    
  sindd   = sin(decmax*(pi/180.)) * sin(y)
  cosdd   = sqrt(1. -sindd*sindd)
  
  photoperiod(:) = 0.
  
  do t=1,ntime
     
     con = 2.*real(t-1)/real(ntime)

     ! Hour angle: Number of minutes before noon divided by 4
     hrangle = cos( pi*(1.-con))
     sinll = sin(pi*lat/180.)

     cosll = sqrt(1.-sinll**2.)
     zenith = acos( max(sinll*sindd + cosll*cosdd*hrangle,0.))*180./pi  ! solar zenith angle in degrees
     
     where (zenith.lt.85.) 
        photoperiod = photoperiod + 24./real(ntime)  ! augment photoperiod length by finite time step
     endwhere     

  enddo
  
end function photoperiod

end module photoperiod_module
