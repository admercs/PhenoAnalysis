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

module parameter_module

  implicit none

  real(kind=4) :: eps
  real(kind=8) :: deps

  real(kind=4), parameter :: earthradius = 6371.  ! earth's radius (km)
  real(kind=4), parameter :: pi = 3.141592653589793       ! pi (single precision)
  real(kind=4), parameter :: degrad = 0.017453292519943   ! pi / 180
  real(kind=4), parameter :: raddeg = 57.295779513082323  ! 180 / pi
  real(kind=8), parameter :: dpi = 3.141592653589793d0      ! pi (double precision)
  real(kind=8), parameter :: ddegrad = 0.017453292519943d0  ! pi / 180
  real(kind=8), parameter :: draddeg = 57.295779513082323d0 ! 180 / pi

  real(kind=4), parameter :: miss_real = -9999.0

contains

  subroutine init_parameters

    implicit none

    ! determine smallest difference in 32 bit and 64bit precision
    ! and add some safe range
    eps = 10.0*epsilon(0.0)
    deps = 10.d0*epsilon(0.d0)


!    print*,eps
!    print*,deps

  end subroutine init_parameters

end module parameter_module
