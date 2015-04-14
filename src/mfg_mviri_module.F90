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

module mfg_mviri_module

  private

  real(kind=8), parameter :: pi = 3.141592653589793d0   ! pi
  real(kind=8), parameter :: deg2rad  = pi/180.d0       ! convert from degrees to radians
  real(kind=8), parameter :: rad2deg  = 180.d0/pi       ! convert from radians to degrees 

  real(kind=8), parameter :: coff_vis  = 2500.d0
  real(kind=8), parameter :: loff_vis  = 2500.d0
  real(kind=8), parameter :: coff_wvir = 1250.d0
  real(kind=8), parameter :: loff_wvir = 1250.d0

  real(kind=8), parameter :: r_eq  = 6378.140d0      ! radius from Earth centre to equator [km]
  real(kind=8), parameter :: r_pol = 6356.755d0      ! radius from Earth centre to pole [km]
  real(kind=8), parameter :: oblate = 1.d0 / 298.257d0

  public :: mfg_mviri2ll
  public :: ll2mfg_mviri

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! This subroutine converts digital to geographical co-ordinates.
  ! taken from: Meteosat First Generation User Handbook
  !             EUM/OPS/USR/10/1537
  !             http://www.eumetsat.int/idcplg?IdcService=GET_FILE&dDocName=PDF_TD06_MARF&RevisionSelectionMethod=LatestReleased
  ! adapted by: Rebekka Posselt (MeteoSwiss)
  !
  ! Input parameters:
  ! lin          - line number, measured from southern end of frame
  ! col          - column number, measured from eastern end of frame
  ! satlon      - longitude of the subsatellite point (usually 0°)
  !
  ! Line and column values are real numbers to enable sub-pixel
  ! accuracy. Integer values correspond to the middle of the pixel,
  ! e.g.(500, 800) would correspond to the middle of the pixel with
  ! corners (499.5, 799.5), (499.5, 800.5), (500.5, 799.5), (500.5, 800.5).
  !
  ! Output parameters
  ! lat          - latitude of this pixel (degrees North from Equator)
  ! lon          - longitude of this pixel (degrees East from Greenwich)
  !
  ! (c) EUMETSAT 1997
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mfg_mviri2ll (col, lin, lon, lat, satlon, satlat, sathgt, missing, wvir)

    implicit none

    ! arguments
    real(kind=4), intent(in)    :: lin, col
    real(kind=4), intent(out)   :: lat, lon
    real(kind=4), intent(in)    :: satlon
    real(kind=4), intent(in)    :: satlat
    real(kind=4), intent(in)    :: sathgt
    real(kind=4), intent(in)    :: missing
    logical, intent(in)         :: wvir

    ! local
    real(kind=8) :: aline, asamp, tanal, tanas, det, k, cenlat, lontmp
    real(kind=8) :: step, x, y, z, a, b, c, p, q, r
    real(kind=8) :: coff, loff

    if (wvir) then
       coff = coff_wvir
       loff = loff_wvir
    else
       coff = coff_vis
       loff = loff_vis
    end if

    ! Step is the radiometer step as seen by the spacecraft,
    ! in degrees. The image represents an 18° x 18° field
    ! of view divided up on an equi-angular basis. For this
    ! program an IR channel of 2500 x 2500 is assumed but
    ! in the real code the size of each channel must be accounted for.
    step = 18.d0 / (2.d0 * coff) 

    ! Convert lin/col values to angular offsets from centre point
    asamp = - (real(col,kind=8) - 0.5d0 - coff) * step 
    aline =   (real(lin,kind=8) - 0.5d0 - loff) * step

    asamp = asamp * deg2rad
    aline = aline * deg2rad

    ! Calculate tangents of angles
    tanal = tan(aline)
    tanas = tan(asamp)

    ! Calculate components of an arbitrary vector from the spacecraft
    ! in the viewing direction.
    p = -1.d0
    q = tanas  
    r = tanal * sqrt (1.d0+ q*q)

    ! The location of the point on the earth can be identified by
    ! solving a quadratic equation for the intersection between
    ! the earths surface and the viewing line from the spacecraft.
    ! If this equation has no real roots then there is no intersection;
    ! otherwise the r_equired root is the one nearer to the spacecraft
    ! (on the visible side of the earth).
    a = q*q + (r*r_eq/r_pol)**2.d0 + p*p
    b = 2.d0 * real(sathgt,kind=8) * p
    c = sathgt * real(sathgt,kind=8) - r_eq*r_eq

    ! Calculate determinant. If it is negative (no real roots to
    ! quadratic equation) there is no intersection between the
    ! line of sight and the disc and so the pixel does not correspond
    ! to visible data.
    det = b*b - 4.d0 * a * c
    if (det .le. 0.d0) then
       lat = missing
       lon = missing
    else
       k = (- b - sqrt(det)) / (2.d0 * a)
       x = sathgt + k * p
       y = k * q
       z = k * r
       lontmp = atan (y/x)
       cenlat = atan (z * cos(lontmp) / x)

       ! This is the geocentric latitude. Convert it to the geodetic
       ! (or geographic) latitude before returning it to the calling program
       lat = real(atan ( tan(cenlat) / ((1.d0 - oblate)**2d0) ) * rad2deg,kind=4)
       lon = real(lontmp * rad2deg,kind=4) + satlon
    end if

    return

  end subroutine mfg_mviri2ll

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! This subroutine converts pixel position from geographical
  ! (lat / long) co-ordinates to digital (lin / col) co-ordinates.
  ! taken from: Meteosat First Generation User Handbook
  !             EUM/OPS/USR/10/1537
  !             http://www.eumetsat.int/idcplg?IdcService=GET_FILE&dDocName=PDF_TD06_MARF&RevisionSelectionMethod=LatestReleased
  ! adapted by: Rebekka Posselt (MeteoSwiss)
  !
  ! Input parameters:
  ! lat         - latitude of pixel (North is pos, South is neg)
  ! lon         - longitude of pixel (East is pos, West is neg)
  ! satlon      - longitude of the subsatellite point (usually 0°)
  !
  ! Note that these are standard geographic co-ordinates as would
  ! be found in an atlas.
  !
  ! Output parameters
  ! lin         - line number, measured from southern end of frame
  ! col         - column number, measured from eastern end of frame
  !
  ! (c) EUMETSAT 1997
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ll2mfg_mviri (lon, lat, col, lin, satlon, satlat, sathgt, missing, wvir)

    implicit none

    ! arguments
    real(kind=4), intent (in)     :: lat, lon
    real(kind=4), intent (out)    :: lin, col
    real(kind=4), intent (in)     :: satlon
    real(kind=4), intent (in)     :: satlat
    real(kind=4), intent (in)     :: sathgt
    real(kind=4), intent (in)     :: missing
    logical, intent(in)           :: wvir

    ! local
    real(kind=8) :: latrad, lonrad
    real(kind=8) :: x, y, z, rtheta, aline, asamp, dotprod
    real(kind=8) :: nlines, nsamps
    real(kind=8) :: coff, loff

    if (wvir) then
       coff = coff_wvir
       loff = loff_wvir
    else
       coff = coff_vis
       loff = loff_vis
    end if

    ! check if the values are sane, otherwise return error value
    if ((lat .lt. -90.0) .or. (lat .gt. 90.0) .or. &
         (lon .lt. -180.0) .or. (lon .gt. 180.0) ) then
       lin = missing
       col = missing
    else
       ! convert inputs to radians
       latrad    = real(lat,kind=8) * deg2rad
       lonrad    = real(lon-satlon,kind=8) * deg2rad

       !  Convert geodetic latitudes (as input) to geocentric latitudes
       !  for use within the algorithm
       latrad = atan( ((1.d0 - oblate)**2.d0) * tan(latrad) )

       !  Calculate rtheta. This is the distance from the earth centre to
       !  a point on the surface at latitude 'lat'.
       rtheta = (r_eq*r_pol) / sqrt (r_pol**2.d0*cos(latrad)**2.d0 + r_eq**2*sin(latrad)**2.d0)

       !  Calculate Cartesian co-ordinates of target point. This is
       !  basic geometry. The co-ordinate system is geocentric with
       !  the x-axis towards the spacecraft, the y-axis to the East
       !  and the x-axis towards the N pole.
       x = rtheta * cos(latrad) * cos(lonrad)
       y = rtheta * cos(latrad) * sin(lonrad)
       z = rtheta * sin(latrad)

       !  Check for invisibility. This is done using the basic geometric
       !  theorem that the dot product of two vectors A and B is equal
       !  to
       !  |A||B| cos (theta)
       ! 
       !  where theta is the angle between them. In this case, the test
       !  is simple. The horizon is defined as the locus of points where
       !  the local normal is perpendicular to the spacecraft sightline
       !  vector. All visible points have (theta) less than 90°
       !  and all invisible points have (theta) greater than 90°.
       !  The test therefore reduces to whether the sign of the dot
       !  product is +ve or -ve; if it is -ve the point is invisible.
       !  The vector from the point to the spacecraft has components
       !  Rs-x, -y, -z where Rs is the distance from the origin to the
       !  satellite. The vector for the normal has components
       !  x y z(Re/Rp)^2
       dotprod = (sathgt-x)*x - y*y - z*z*((r_eq/r_pol)**2.d0)

       if (dotprod .le. 0.d0) then
          lin = missing
          col = missing
       else
          !  In this co-ordinate system the spacecraft (S) is at position
          !  (sathgt,0,0), the earth centre (O) at (0,0,0) and the point (P)
          !  at (x,y,z). Two additional points need to be defined, so that the
          !  angles from the reference planes to the target point (i.e. the
          !  position of the point in the sensor FOV) can be extracted.
          !  These points are defined by dropping lines perpendicularly from P
          !  onto the equatorial plane and the Greenwich meridian plane.
          !  Their co-ordinates are defined as:
          ! 
          !  O' = (x, y, 0) and O'' = (x, 0, z).
          ! 
          !  With these points, right-angled triangles can be defined SO'P
          !  and SO''P which can be used directly to determine the angular
          !  co-ordinates (aline, asamp) of P in the FOV.
          asamp = atan (y / (real(sathgt,kind=8) - x))
          aline = atan (z / sqrt (y**2.d0 + (real(sathgt,kind=8) - x)**2.d0))

          !  Convert back to degrees
          asamp = asamp * rad2deg
          aline = aline * rad2deg

          !  Calculate lin, col. Note that since pixels are measured from
          !  the right of the image, and the angular conversion was measured in
          !  the x (east) direction, a sign correction has to be included for
          !  pixels. The image represents an 18° x 18° field of view
          !  divided up on an equi-angular basis.
          nlines = 2.d0*loff 
          nsamps = 2.d0*coff
          asamp = asamp / (18.d0 / nsamps)
          aline = aline / (18.d0 / nlines)
          col = real(nsamps/2.d0 + 0.5d0 - asamp,kind=4)
          lin = real(nlines/2.d0 + 0.5d0 + aline,kind=4)
       end if
    end if

  end subroutine ll2mfg_mviri

end module mfg_mviri_module
