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

module modsin_module
!   Transforms MODIS sinusoidal projection into longitude/latitude
!   Formulas implemented by Reto Stockli (Blue Marble Research) 
!   after tilemap3 (v4.0) by Robert Wolfe, November 13, 2003

!   x is the global horizontal pixel coordinate. UL pixel center of UL tile (h00v00)
!   has x coordinate 0. Longitude/Latitude -180.0/0.0 has x coordinate -0.5.
!   x increases from east to east and reaches 43199 for the easternmost pixel of 
!   the easternmost tile (h35vyy) with the 1 km grid.
!   y is the global vertical pixel coordinate. UL pixel center of UL tile (h00v00)
!   has x coordinate 0. Longitude/Latitude -180.0/90.0 has y coordinate -0.5.
!   y increases from north to south and reaches 21599 for the southernmost pixel of
!   the southernmost tile (hxxv17) with the 1 km grid.

  implicit none

  private

  public :: ll2modsin
  public :: modsin2ll
  public :: modsinxy2tile
  public :: modsintile2xy

contains

  subroutine ll2modsin(lon, lat, x, y, res, missing)

    implicit none

    ! arguments
    real(kind=4), intent(in)      :: lon,lat    ! longitude / latitude (degrees E / degrees N)
    real(kind=4), intent(inout)   :: x,y        ! global MODIS pixel coordinates
    real(kind=4), intent(in)      :: res        ! spatial resolution of grid (1=1km, 0.5=500m, 0.25=250m)
    real(kind=4), intent(in)      :: missing    ! missing data value for input/output

    ! local variables
    real(kind=4), parameter :: pi = 3.14159265358979 ! Pi
    integer, parameter      :: nhtile = 36 ! number of horizontal global tiles
    integer, parameter      :: nvtile = 18 ! number of vertical global tiles
    integer                 :: nsample, nline
    real(kind=4)            :: macheps

    macheps = 10.0*epsilon(0.0)

    if (res.gt.macheps) then
       nsample = nint(1200.0 / res) ! number of horizontal samples per tile 
       nline = nint(1200.0 / res) ! number of vertical lines per tile

       if ((abs(lon-missing).gt.macheps).and.(abs(lat-missing).gt.macheps).and.&
            (lon.ge.-180.0).and.(lon.le.180.0).and.(lat.ge.-90.0).and.(lat.le.90.0))  then
          x = (cos(lat / 90.0 * 0.5 * pi) * lon / 180.0 * 0.5 + 0.5) * real(nhtile) * real(nsample) - 0.5
          y = (-lat / 90.0 * 0.5 + 0.5) * real(nvtile) * real(nline) - 0.5       
       else
          x = missing
          y = missing
       end if
    else
       x = missing
       y = missing
    end if

  end subroutine ll2modsin

  subroutine modsin2ll(x, y, lon, lat, res, missing)

    implicit none

    ! arguments
    real(kind=4), intent(in)      :: x,y        ! global MODIS pixel coordinates
    real(kind=4), intent(inout)   :: lon,lat    ! longitude / latitude (degrees E / degrees N)
    real(kind=4), intent(in)      :: res        ! spatial resolution of grid (1=1km, 0.5=500m, 0.25=250m)
    real(kind=4), intent(in)      :: missing    ! missing data value for input/output

    ! local variables
    real(kind=4), parameter :: pi = 3.14159265358979 ! Pi
    integer, parameter      :: nhtile = 36 ! number of horizontal global tiles
    integer, parameter      :: nvtile = 18 ! number of vertical global tiles
    integer                 :: nsample, nline
    real(kind=4)            :: macheps
    
    macheps = 10.0*epsilon(0.0)

    if (res.gt.macheps) then
       nsample = nint(1200.0 / res) ! number of horizontal samples per tile 
       nline = nint(1200.0 / res) ! number of vertical lines per tile
       
       if ((abs(x-missing).gt.macheps).and.(abs(y-missing).gt.macheps)) then
          lat = -((y + 0.5) / (real(nvtile) * real(nline)) - 0.5) * 90.0 / 0.5
          lon = ((x + 0.5) / (real(nhtile) * real(nsample)) - 0.5) * 180.0 / 0.5 / cos(lat / 90.0 * 0.5 * pi)
          
          if ((lon.lt.-180.0).or.(lon.gt.180.0).or.(lat.lt.-90.0).or.(lat.gt.90.0)) then
             lon = missing
             lat = missing
          end if
       else 
          lon = missing
          lat = missing
       end if
    else
       lon = missing
       lat = missing
    end if
 
  end subroutine modsin2ll

  subroutine modsinxy2tile(x, y, htile, vtile, sample, line, res, missing)

    implicit none

    ! arguments
    real(kind=4), intent(in)      :: x,y          ! input global MODIS pixel coordinates
    integer, intent(out)          :: htile, vtile ! output horizontal/vertical tile number (0-based)
    real(kind=4), intent(out)     :: sample, line ! output tile sample/line number (0-based)
    real(kind=4), intent(in)      :: res          ! spatial resolution of grid (1=1km, 0.5=500m, 0.25=250m)
    real(kind=4), intent(in)      :: missing      ! missing data value for input/output

    ! local variables
    integer, parameter      :: nhtile = 36 ! number of horizontal global tiles
    integer, parameter      :: nvtile = 18 ! number of vertical global tiles
    integer                 :: nsample, nline
    real(kind=4)            :: macheps

    macheps = 10.0*epsilon(0.0)

    if (res.gt.macheps) then
       nsample = nint(1200.0 / res) ! number of horizontal samples per tile 
       nline = nint(1200.0 / res) ! number of vertical lines per tile

       ! calculate modis tile number and tile sample/line number
       htile = min(max(int((x + 0.5) / real(nsample)),0),nhtile-1)
       vtile = min(max(int((y + 0.5) / real(nline)),0),nvtile-1)

       if ((abs(x-missing).gt.macheps).and.(abs(y-missing).gt.macheps)) then
          sample = x - real(htile)*real(nsample)
          line = y - real(vtile)*real(nline)
       else
          htile = 0
          vtile = 0
          sample = missing
          line = missing
       end if
    else
       htile = 0
       vtile = 0
       sample = missing
       line = missing
    end if

  end subroutine modsinxy2tile

  subroutine modsintile2xy(htile, vtile, sample, line, x, y, res, missing)

    implicit none

    ! arguments
    integer, intent(in)           :: htile, vtile ! input horizontal/vertical tile number (0-based)
    real(kind=4), intent(in)      :: sample, line ! input tile sample/line number (0-based)
    real(kind=4), intent(out)     :: x,y          ! output global MODIS pixel coordinates
    real(kind=4), intent(in)      :: res          ! spatial resolution of grid (1=1km, 0.5=500m, 0.25=250m)
    real(kind=4), intent(in)      :: missing      ! missing data value for input/output

    ! local variables
    integer                 :: nsample, nline
    real(kind=4)            :: macheps

    macheps = 10.0*epsilon(0.0)

    if (res.gt.macheps) then
       nsample = nint(1200.0 / res) ! number of horizontal samples per tile 
       nline = nint(1200.0 / res) ! number of vertical lines per tile
       
       if ((abs(sample-missing).gt.macheps).and.(abs(line-missing).gt.macheps)) then
          x = sample + real(htile)*real(nsample)
          y = line + real(vtile)*real(nline)
       else
          x = missing
          y = missing
       end if
    else
       x = missing
       y = missing
    end if

  end subroutine modsintile2xy

end module modsin_module
