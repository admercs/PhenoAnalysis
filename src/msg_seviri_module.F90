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

module msg_seviri_module

  private

  real(kind=8), parameter :: pi = 3.141592653589793d0   ! pi
  real(kind=8), parameter :: deg2rad  = pi/180.d0       ! convert from degrees to radians
  real(kind=8), parameter :: rad2deg  = 180.d0/pi       ! convert from radians to degrees 

  real(kind=8), parameter :: coff_visir  = 1856.d0
  real(kind=8), parameter :: loff_visir  = 1856.d0
  real(kind=8), parameter :: cfac_visir  = -781648343.d0
  real(kind=8), parameter :: lfac_visir  = -781648343.d0

  real(kind=8), parameter :: coff_hrvis  = 5566.d0
  real(kind=8), parameter :: loff_hrvis  = 5566.d0
  real(kind=8), parameter :: cfac_hrvis  = -2344944937.d0
  real(kind=8), parameter :: lfac_hrvis  = -2344944937.d0

  real(kind=8), parameter    :: r_eq  = 6378.1690d0      ! radius from Earth centre to equator [km]
  real(kind=8), parameter    :: r_pol = 6356.5838d0      ! radius from Earth centre to pole [km]

  public :: msg_seviri2ll
  public :: ll2msg_seviri

contains
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  ! subroutine msg_seviri2ll 
  !  
  !  PURPOSE: 
  !  return the geograhic latitude and longitude of an MSG image   
  !    for a given pair of latitude/longitude.             
  !    (based on the formulas given in Ref. [1])                
  !                                                        
  !                                                        
  !  DEPENDENCIES:                                         
  !    none                                                 
  !                                                        
  !                                                        
  !  REFERENCE:                                            
  !  [1] LRIT/HRIT Global Specification                     
  !      (CGMS 03, Issue 2.6, 12.08.1999)                  
  !      for the parameters used in the program           
  !  [2] MSG Ground Segment LRIT/HRIT Mission Specific 
  !      Implementation, EUMETSAT Document, 
  !      (EUM/MSG/SPE/057, Issue 6, 21. June 2006).
  !                                                        
  !  MODIFICATION HISTORY:                                 
  !    Version 1.01
  !  08.08.2008 removed a bug in longi = atan(s2/s1) statement
  !    Copyright(c) EUMETSAT 2005, 2009
  !                                                        
  !  INPUT:                                                
  !    lin   (int) lin-value of the pixel
  !    col   (int) colum-value of the pixel
  !                                                        
  !  OUTPUT:                                               
  !    lat (double) geographic Latitude of the wanted pixel [Degrees]
  !    lon (double) geographic Longitude of the wanted pixel [Degrees]
  !                                                        
  !                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine msg_seviri2ll(col, lin, lon, lat, satlon, satlat, sathgt, missing, hrvis)

    implicit none

    ! arguments
    real(kind=4), intent (IN)    :: col, lin 
    real(kind=4), intent (OUT)   :: lat, lon
    real(kind=4), intent (IN)    :: satlon
    real(kind=4), intent (IN)    :: satlat
    real(kind=4), intent (IN)    :: sathgt
    real(kind=4), intent (IN)    :: missing
    logical, intent (IN) :: hrvis

    ! local
    real(kind=8) :: s1, s2, s3, sn, sd, sxy
    real(kind=8) :: x, y
    real(kind=8) :: longi, lati

    real(kind=8) :: sa

    real(kind=8) :: cfac,lfac,coff,loff

    if (hrvis) then
       cfac = cfac_hrvis
       lfac = lfac_hrvis
       coff = coff_hrvis
       loff = loff_hrvis
    else
       cfac = cfac_visir
       lfac = lfac_visir
       coff = coff_visir
       loff = loff_visir
    endif

    !  calculate viewing angle of the satellite by use of the equation 
    !  on page 28, Ref [1]. 
    x = (2.d0**16.d0 * ( col - coff )) / cfac
    y = (2.d0**16.d0 * ( lin - loff )) / lfac

    !  now calculate the inverse projection using equations on page 25, Ref. [1]  

    !  first check for visibility, whether the pixel is located on the Earth 
    !  surface or in space. 
    !  To do this calculate the argument to sqrt of "sd", which is named "sa". 
    !  If it is negative then the sqrt will return NaN and the pixel will be 
    !  located in space, otherwise all is fine and the pixel is located on the 
    !  Earth surface.
    sa =  (sathgt * cos(x) * cos(y) )**2.d0 - (cos(y)*cos(y) + &
         1.006803d0 * sin(y)*sin(y)) * 1737121856.d0

    ! take care if the pixel is in space, that an error code will be returned
    if ( sa .le. 0.0 ) then
       lat = missing
       lon = missing
       return 
    end if

    ! now calculate the rest of the formulas using eq. on page 25 Ref [1]
    sd = sqrt(sa)
    sn = (sathgt * cos(x) * cos(y) - sd) / ( cos(y)*cos(y) + 1.006803d0 * sin(y)*sin(y) )

    s1 = sathgt - sn * cos(x) * cos(y)
    s2 = sn * sin(x) * cos(y)
    s3 = -sn * sin(y)

    sxy = sqrt( s1*s1 + s2*s2 )

    ! using the previous calculations now the inverse projection can be
    ! calculated, which means calculating the lat./long. from the pixel
    ! row and column by equations on page 25, Ref [1].
    longi = atan(s2/s1)
    lati  = atan((1.006803d0*s3)/sxy)

    ! convert from radians into degrees
    lat = real(lati*rad2deg,kind=4)
    lon = real(longi*rad2deg,kind=4) + satlon

  end subroutine msg_seviri2ll


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! subroutine ll2msg_seviri                                 
  !                                                       
  ! PURPOSE:                                              
  !   return the pixel column and line of an MSG image 
  !   for a given pair of geographic latitude/longitude.                   
  !   (based on the formulas given in Ref. [1])                
  !                                                       
  !                                                       
  ! DEPENDENCIES:                                         
  !   none                                       
  !                                                       
  !                                                       
  ! REFERENCE:                                            
  ! [1] LRIT/HRIT Global Specification                     
  !     (CGMS 03, Issue 2.6, 12.08.1999)                  
  !     for the parameters used in the program.
  ! [2] MSG Ground Segment LRIT/HRIT Mission Specific 
  !     Implementation, EUMETSAT Document, 
  !     (EUM/MSG/SPE/057, Issue 6, 21. June 2006).
  !                                                       
  !                                                       
  ! MODIFICATION HISTORY:
  !   Version 1.01
  !    Copyright(c) EUMETSAT 2005, 2009
  !                                                       
  !                                                       
  ! INPUT:                                                
  !   lat   (double)     geographic Latitude of a point [Degrees] 
  !   lon   (double)     geographic Longitude of a point [Degrees]
  !                                                       
  !                                                       
  ! OUTPUT:                                               
  !   lin (int)          lin-value of the pixel
  !   col (int)          col-value of the pixel
  !                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ll2msg_seviri( lon, lat, col, lin, satlon, satlat, sathgt, missing, hrvis)

    implicit none

    ! arguments
    real(kind=4), intent (IN)    :: lat, lon
    real(kind=4), intent (OUT)   :: col, lin
    real(kind=4), intent (IN)    :: satlon
    real(kind=4), intent (IN)    :: satlat
    real(kind=4), intent (IN)    :: sathgt
    real(kind=4), intent (IN)    :: missing
    logical, intent (IN) :: hrvis

    ! local
    real(kind=8) :: latrad, lonrad
    real(kind=8) :: cenlat
    real(kind=8) :: r1, r2, r3, rn, re, rl
    real(kind=8) :: xx, yy
    real(kind=8) :: dotprod

    real(kind=8) :: cfac,lfac,coff,loff

    if (hrvis) then
       cfac = cfac_hrvis
       lfac = lfac_hrvis
       coff = coff_hrvis
       loff = loff_hrvis
    else
       cfac = cfac_visir
       lfac = lfac_visir
       coff = coff_visir
       loff = loff_visir
    endif

    ! check if the values are sane, otherwise return error value
    if ((lat .lt. -90.0) .or. (lat .gt. 90.0) .or. &
         (lon .lt. -180.0) .or. (lon .gt. 180.0) ) then
       lin = missing
       col = missing
       return
    end if

    ! convert them to radians 
    latrad    = real(lat,kind=8)*deg2rad
    lonrad   = real(lon-satlon,kind=8)*deg2rad

    ! calculate the geocentric latitude from the       
    ! geographic one using equations on page 24, Ref. [1] 
    cenlat = atan ( (0.993243d0*(sin(latrad)/cos(latrad)) ))

    ! using cenlat calculate the length from the Earth 
    ! centre to the surface of the Earth ellipsoid    
    ! equations on page 23, Ref [1]
    re = R_POL / sqrt( (1.d0 - 0.00675701d0 * cos(cenlat) * cos(cenlat) ) )

    ! calculate the forward projection using equations on page 24, Ref. [1]
    rl = re
    r1 = sathgt - rl * cos(cenlat) * cos(lonrad)
    r2 = - rl *  cos(cenlat) * sin(lonrad)
    r3 = rl * sin(cenlat)
    rn = sqrt( r1*r1 + r2*r2 +r3*r3 )

    ! check for visibility, whether the point on the Earth given by the
    ! latitude/lonradtude pair is visible from the satellte or not. This 
    ! is given by the dot product between the vectors of:
    ! 1) the point to the spacecraft,
    ! 2) the point to the centre of the Earth.
    ! If the dot product is positive the point is visible otherwise it
    ! is invisible.
    dotprod = r1*(rl * cos(cenlat) * cos(lonrad)) - r2*r2 - r3*r3*((r_EQ/R_POL)**2.d0)

    if (dotprod .le. 0.d0 ) then
       col = missing
       lin = missing
       return
    end if

    ! the forward projection is x and y 
    xx = atan( (-r2/r1) )
    yy = asin( (-r3/rn) )

    ! convert to pixel column and row using the scaling functions on 
    ! page 28, Ref. [1]. And finding nearest integer value for them. 
    col = real(coff + xx *  2.d0**(-16.d0) * cfac,kind=4)
    lin = real(loff + yy *  2.d0**(-16.d0) * lfac,kind=4)

  end subroutine ll2msg_seviri

end module msg_seviri_module
