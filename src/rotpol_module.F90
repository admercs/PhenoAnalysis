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

! Copyright 2014 Blue Marble Research
! 
! This file is part of the <phenoanalysis> project.
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

module rotpol_module

! Conversion routines between regular cylindrical and rotated pole coordinate system
! 2009/06/01 Reto Stockli (Blue Marble Research), based on COSMO code by Oliver Fuhrer (MeteoSwiss)

  implicit none

  private

  public :: ll2rpol
  public :: rpol2ll

contains

  subroutine ll2rpol (lon,lat,rlon,rlat,pollon,pollat,missing)

    ! Converts longitude / latitude array from regular cylindrical to rotated pole coordinates

    implicit none

    ! Arguments
    real (kind=4), intent (in) :: &
         pollon,  & ! longitude of the rotated north pole
         pollat,  & ! latitude of the rotated north pole
         lon,     & ! longitude in the non-rotated system
         lat        ! latitude in the non-rotated system
    real(kind=4), intent (inout) :: &
         rlon,    & ! longitude in the rotated system
         rlat       ! latitutde in the rotated system
    real(kind=4), intent(in) :: missing

    ! Local variables
    real (kind=4) :: &
         zsinpol, zcospol, zlonpol, zlat, zlon, zarg1, zarg2

    ! pi
    real (kind=4), parameter :: pi=3.14159265358979

    ! conversion factors
    real (kind=4), parameter :: deg2rad=pi/180.0
    real (kind=4), parameter :: rad2deg=180.0/pi

    real(kind=4)            :: macheps

    macheps = 10.0*epsilon(0.0)

    if ((abs(lon-missing).gt.macheps).and.(abs(lat-missing).gt.macheps).and. &
         (lon.ge.-180.0).and.(lon.le.180.0).and.(lat.ge.-90.0).and.(lat.le.90.0).and. &
         (pollon.ge.-180.0).and.(pollon.le.180.0).and.(pollat.ge.-90.0).and.(pollat.le.90.0)) then

       zsinpol  = sin (deg2rad * pollat)
       zcospol  = cos (deg2rad * pollat)
       zlonpol  =      deg2rad * pollon

       zlat = deg2rad * lat
       zlon = deg2rad * lon

       ! rotated longitude
       zarg1  = - sin (zlon-zlonpol) * cos(zlat)
       zarg2  = - zsinpol * cos(zlat) * cos(zlon-zlonpol) + zcospol * sin(zlat)

       if (abs(zarg2 - 0.0) < macheps) zarg2 = macheps

       rlon = rad2deg * atan2 (zarg1,zarg2)

       ! rotated latitude
       zarg1 = sin (zlat) * zsinpol
       zarg2 = cos (zlat) * zcospol * cos (zlon - zlonpol)

       rlat = rad2deg * asin (zarg1 + zarg2)

    else
       rlon = missing
       rlat = missing
    end if

  end subroutine ll2rpol

  subroutine rpol2ll (rlon,rlat,lon,lat,pollon,pollat,missing)
    ! Converts longitude / latitude array from rotated pole to regular cylindrical coordinates

    implicit none

    ! Arguments
    real (kind=4), intent (in) :: &
         pollon,  & ! longitude of the rotated north pole
         pollat,  & ! latitude of the rotated north pole
         rlon, &    ! longitude in the rotated system
         rlat       ! latitude in the rotated system
    real(kind=4), intent (inout) :: &
         lon, &     ! longitude in the non-rotated system
         lat        ! latitude in the non-rotated system
    real(kind=4), intent(in) :: missing

    ! Local variables
    real (kind=4) :: &
         zsinpol, zcospol, zlonpol, zlats, zlons, zarg, zarg1, zarg2

    ! pi
    real (kind=4), parameter :: pi=3.14159265358979

    ! conversion factors
    real (kind=4), parameter :: deg2rad=pi/180.0
    real (kind=4), parameter :: rad2deg=180.0/pi

    real(kind=4)            :: macheps

    macheps = 10.0*epsilon(0.0)

    if ((abs(rlon-missing).gt.macheps).and.(abs(rlat-missing).gt.macheps).and. &
         (rlon.ge.-180.0).and.(rlon.le.180.0).and.(rlat.ge.-90.0).and.(rlat.le.90.0).and. &
         (pollon.ge.-180.0).and.(pollon.le.180.0).and.(pollat.ge.-90.0).and.(pollat.le.90.0))  then

       zsinpol = sin (deg2rad * pollat)
       zcospol = cos (deg2rad * pollat)
       zlonpol =      deg2rad * pollon

       zlats = deg2rad * rlat
       zlons = deg2rad * rlon

       ! calculate longitude
       zarg1 = sin (zlonpol) * (-zsinpol * cos(zlons) * cos(zlats) + &
            zcospol *  sin(zlats)) - cos (zlonpol) * sin(zlons) * cos(zlats)
       zarg2 = cos (zlonpol) * (-zsinpol * cos(zlons) * cos(zlats)  + &
            zcospol * sin(zlats)) +  sin (zlonpol) * sin(zlons) * cos(zlats)

       if (abs(zarg2 - 0.0) < macheps) zarg2 = macheps
       lon = rad2deg * atan2(zarg1,zarg2)

       ! calculate latitude
       zarg  = zcospol * cos (zlats) * cos (zlons) + zsinpol * sin (zlats)
       lat = rad2deg * asin (zarg)

    else
       lon = missing
       lat = missing
    end if

  end subroutine rpol2ll

end module rotpol_module
