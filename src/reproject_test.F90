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

program reproject_test

  use reproject_module
  use proj_module
  use modsin_module
  use rotpol_module
  use histogram_module
  use quicksort_module

  implicit none

  real(kind=4)              :: missing = -9999.0
  integer(kind=4), parameter:: missingi4 = -9999  
  real(kind=4), parameter   :: missingf4 = -9999.0  
  integer, parameter        :: type = 5
  integer, parameter        :: csize = 2

  character(len=100)        :: itype, otype
  character(len=100)        :: ioptargs, ooptargs
  real(kind=4)              :: ixmin, ixmax, iymin, iymax, idx, idy
  real(kind=4)              :: oxmin, oxmax, oymin, oymax, odx, ody
  integer                   :: inx, iny, onx, ony
  integer                   :: i,j
  real(kind=4)              :: x(1,1),y(1,1),lon(1,1),lat(1,1)
  real(kind=4), allocatable :: idataf4(:,:), odataf4(:,:)
  integer(kind=4), allocatable :: idatai4(:,:), odatai4(:,:)
  real(kind=4), allocatable :: ixgrid(:,:), iygrid(:,:), oxgrid(:,:), oygrid(:,:)

  if (type.eq.1) then
     itype = "latitude_longitude"
     ioptargs = ""
     ixmin = 6.0
     ixmax = 11.0
     iymin = 45.0
     iymax = 48.0
     idx = 0.01
     idy = 0.01
     
     otype = "swisscors"
     ooptargs = ""
     oxmin = 500000.0
     oxmax = 800000.0
     oymin = 100000.0
     oymax = 400000.0
     odx = 2000.0
     ody = 2000.0
  end if

  if (type.eq.2) then
     itype = "latitude_longitude"
     ioptargs = ""
     ixmin = 6.0
     ixmax = 11.0
     iymin = 45.0
     iymax = 48.0
     idx = 0.5
     idy = 0.5
     
     otype = "modsin"
     ooptargs = 'res=1.0'
     oxmin = 22100.0
     oxmax = 22450.0
     oymin = 5000.0
     oymax = 5300.0
     odx = 10.0
     ody = 10.0
  end if

 
  if (type.eq.3) then
     itype = "latitude_longitude"
     ioptargs = ""
     ixmin = 6.0
     ixmax = 9.0
     iymin = 45.0
     iymax = 48.0
     idx = 0.1
     idy = 0.1
     
     otype = "rotated_latitude_longitude"
     ooptargs = 'pollon=-172.5;pollat=43.5'
     oxmin = -1.0
     oxmax = 1.0
     oymin = -1.5
     oymax = 0.5
     odx = 0.2
     ody = 0.2
  end if

  if (type.eq.4) then
     itype = "rotated_latitude_longitude"
     ioptargs = "pollon=-172.5;pollat=43.5"
     ixmin = 4.0
     ixmax = 6.0
     iymin = -47.0
     iymax = -45.0
     idx = 0.1
     idy = 0.1

     otype = "rotated_latitude_longitude"
     ooptargs = "pollon=-172.5;pollat=43.5"
     oxmin = 6.0
     oxmax = 6.0
     oymin = -45.0
     oymax = -45.0
     odx = 0.0
     ody = 0.0     
  end if

  if (type.eq.5) then
     itype = "geosat"
     ioptargs = "satname=msg;sensor=seviri;channel=visir;satlon=0.0;satlat=0.0;sathgt=42164.0"
     ixmin = 1688.0
     ixmax = 1690.0
     iymin = 3306.0
     iymax = 3309.0
     idx = 1.0
     idy = 1.0

     otype = "latitude_longitude"
     ooptargs = ""
     oxmin = 6.9
     oxmax = 7.0
     oymin = 46.75
     oymax = 46.85
     odx = 0.02
     ody = 0.02     
  end if

  if (type.eq.6) then
     itype = "geosat"
     ioptargs = "satname=mfg;sensor=mviri;channel=vis;satlon=0.0;satlat=0.0;sathgt=42164.0"
     ixmin = 1138.0*2.0
     ixmax = 1140.0*2.0
     iymin = 2221.0*2.0
     iymax = 2223.0*2.0
     idx = 1.0
     idy = 1.0

     otype = "latitude_longitude"
     ooptargs = ""
     oxmin = 6.9
     oxmax = 7.0
     oymin = 46.75
     oymax = 46.85
     odx = 0.02
     ody = 0.02     
  end if

  
  ! generate input data with checkerboard pattern
  if ((idx.ne.0.0).and.(idy.ne.0.0)) then
     inx = nint((ixmax - ixmin)/idx)
     iny = nint((iymax - iymin)/idy)
  else
     inx = 1
     iny = 1
  end if

  allocate(idataf4(inx,iny))
  do j=1,iny
     do i=1,inx
        idataf4(i,j) = 9*mod((i-1) / csize,2) * mod((j-1) / csize,2) + &
                     9*mod((i-1) / csize+1,2) * mod((j-1) / csize+1,2) + 1
     end do
!     write(*,'(20F8.3)') idataf4(:,y)
  end do

  allocate(idatai4(inx,iny))
  idatai4 = nint(idataf4)

!  idataf4(inx/2,iny/2) = missingf4
!  idatai4(inx/2,iny/2) = missingi4

  x = ixmin
  y = iymin
  call proj_inv(itype,x,y,lon,lat,ioptargs,missing)

  print*,x,y
  print*,lon,lat
  call proj_fwd(itype,lon,lat,x,y,ioptargs,missing)
  print*,x,y

  stop

  call reproject(itype, ixmin, ixmax, iymin, iymax, idx, idy, ioptargs, &
       otype, oxmin, oxmax, oymin, oymax, odx, ody, ooptargs, &
       missing, &
       idatai4=idatai4, odatai4=odatai4, missingi4=missingi4, &
       ixgrid=ixgrid, iygrid=iygrid, oxgrid=oxgrid, oygrid=oygrid)

  call reproject(itype, ixmin, ixmax, iymin, iymax, idx, idy, ioptargs, &
       otype, oxmin, oxmax, oymin, oymax, odx, ody, ooptargs, &
       missing, &
       idataf4=idataf4, odataf4=odataf4, missingf4=missingf4, &
       ixgrid=ixgrid, iygrid=iygrid, oxgrid=oxgrid, oygrid=oygrid)

!  onx = size(odataf4,1)
!  ony = size(odataf4,2)
  onx = size(odatai4,1)
  ony = size(odatai4,2)

  do j=1,ony
     write(*,'(20I10)') odatai4(:,j)     
  end do
 
 do j=1,ony
     write(*,'(20F10.3)') odataf4(:,j)     
  end do

end program reproject_test
