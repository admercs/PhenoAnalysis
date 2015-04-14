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

module proj_module

  private

  public :: proj_fwd
  public :: proj_inv
  public :: getoptargs
  public :: getxybounds
  
contains

  subroutine proj_fwd(type,lon,lat,x,y,optargs,missing)

    use geosat_module
    use modsin_module
    use rotpol_module

    implicit none

    ! arguments
    character(len=*), intent(in)    :: type
    real(kind=4), intent(in)        :: lon(:,:), lat(:,:)
    real(kind=4), intent(out)       :: x(:,:), y(:,:)
    character(len=*), intent(in)    :: optargs
    real(kind=4), intent(in)        :: missing

    ! local variables
    integer :: i,j,ni,nj
    integer :: nargs
    character(len=100), allocatable :: arguments(:)
    real(kind=4), allocatable :: fvalues(:)
    character(len=100), allocatable :: cvalues(:)

    ni = size(lon,1)
    nj = size(lon,2)

    select case (type)
    case ("rotated_latitude_longitude")
       nargs = 2
       allocate(arguments(nargs))
       arguments(1) = "grid_north_pole_longitude"
       arguments(2) = "grid_north_pole_latitude"
       nargs = getoptargs(optargs,arguments,fvalues)
       if (nargs.ne.2) then
          write(*,'(A)') "rotated pole projection requires two arguments: grid_north_pole_longitude"//&
               " and grid_north_pole_latitude"
          x = missing
          y = missing
       else 
          do j=1,nj
             do i=1,ni
                call ll2rpol(lon(i,j),lat(i,j),x(i,j),y(i,j),fvalues(1),fvalues(2),missing)
             end do
          end do
       end if
       deallocate(arguments,fvalues)
    case ("modsin")
       nargs = 1
       allocate(arguments(nargs))
       arguments(1) = "resolution"
       nargs = getoptargs(optargs,arguments,fvalues)
       if (nargs.ne.1) then
          write(*,'(A)') "modsin projection requires one argument: resolution"
          x = missing
          y = missing
       else 
          do j=1,nj
             do i=1,ni
                call ll2modsin(lon(i,j),lat(i,j),x(i,j),y(i,j),fvalues(1),missing)
             end do
          end do
       end if
       deallocate(arguments,fvalues)
    case ("ortho")
    case ("swisscors")
    case ("geosat")
       nargs = 6
       allocate(arguments(nargs))
       arguments(1) = "satname"
       arguments(2) = "sensor"
       arguments(3) = "channel"
       arguments(4) = "satlon"
       arguments(5) = "satlat"
       arguments(6) = "sathgt"
       nargs = getoptargs(optargs,arguments,fvalues,cvalues=cvalues)
       if (nargs.ne.6) then
          write(*,'(A)') "geostationary satellite projection requires six arguments: "
          write(*,'(A)') "satname : short name of satellite (character)"
          write(*,'(A)') "sensor  : name of sensor (character)"
          write(*,'(A)') "channel : name of channel, such as VIS, HRV, IR10.8 (character)"
          write(*,'(A)') "satlon  : sub-satellite longitude (float, degrees east)"
          write(*,'(A)') "satlat  : sub-satellite latitude (float, degrees north)"
          write(*,'(A)') "sathgt  : satellite height above earth center (float, kilometers)"
          x = missing
          y = missing
       else 
          call ll2geosat(lon,lat,x,y,cvalues(1),cvalues(2),cvalues(3),fvalues(4),fvalues(5),fvalues(6),missing)
       end if
       deallocate(arguments,fvalues,cvalues)
    case ("latitude_longitude")
       x = lon
       y = lat
    case default
       write(*,'(A,A,A)') "Geographic Projection Type ",trim(type)," not implemented."
       x = missing
       y = missing
    end select

  end subroutine proj_fwd

  subroutine proj_inv(type,x,y,lon,lat,optargs,missing)

    use geosat_module
    use modsin_module
    use rotpol_module

    implicit none

    ! arguments
    character(len=*), intent(in)    :: type
    real(kind=4), intent(in)        :: x(:,:), y(:,:)
    real(kind=4), intent(out)       :: lon(:,:), lat(:,:)
    character(len=*), intent(in)    :: optargs
    real(kind=4), intent(in)        :: missing

    ! local variables
    integer :: i,j,ni,nj
    integer :: nargs
    character(len=100), allocatable :: arguments(:)
    real(kind=4), allocatable :: fvalues(:)
    character(len=100), allocatable :: cvalues(:)

    ni = size(x,1)
    nj = size(x,2)

    select case (type)
    case ("rotated_latitude_longitude")
       nargs = 2
       allocate(arguments(nargs))
       arguments(1) = "grid_north_pole_longitude"
       arguments(2) = "grid_north_pole_latitude"
       nargs = getoptargs(optargs,arguments,fvalues)
       if (nargs.ne.2) then
          write(*,'(A)') "rotated pole projection requires two arguments: grid_north_pole_longitude"//&
               " and grid_north_pole_latitude"
          lon = missing
          lat = missing
       else 
          do j=1,nj
             do i=1,ni
                call rpol2ll(x(i,j),y(i,j),lon(i,j),lat(i,j),fvalues(1),fvalues(2),missing)
             end do
          end do
       end if
       deallocate(arguments,fvalues)
    case ("modsin")
       nargs = 1
       allocate(arguments(nargs))
       arguments(1) = "resolution"
       nargs = getoptargs(optargs,arguments,fvalues)
       if (nargs.ne.1) then
          write(*,'(A)') "modsin projection requires one argument: resolution"
          lon = missing
          lat = missing
       else 
          do j=1,nj
             do i=1,ni
                call modsin2ll(x(i,j),y(i,j),lon(i,j),lat(i,j),fvalues(1),missing)
             end do
          end do
       end if
       deallocate(arguments,fvalues)
    case ("ortho")
    case ("swisscors")
    case ("geosat")
       nargs = 6
       allocate(arguments(nargs))
       arguments(1) = "satname"
       arguments(2) = "sensor"
       arguments(3) = "channel"
       arguments(4) = "satlon"
       arguments(5) = "satlat"
       arguments(6) = "sathgt"
       nargs = getoptargs(optargs,arguments,fvalues,cvalues=cvalues)
       if (nargs.ne.6) then
          write(*,'(A)') "geostationary satellite projection requires six arguments: "
          write(*,'(A)') "satname : short name of satellite (character)"
          write(*,'(A)') "sensor  : name of sensor (character)"
          write(*,'(A)') "channel : name of channel, such as VIS, HRV, IR10.8 (character)"
          write(*,'(A)') "satlon  : sub-satellite longitude (float, degrees east)"
          write(*,'(A)') "satlat  : sub-satellite latitude (float, degrees north)"
          write(*,'(A)') "sathgt  : satellite height above earth center (float, kilometers)"
          lon = missing
          lat = missing
       else 
          call geosat2ll(x,y,lon,lat,cvalues(1),cvalues(2),cvalues(3),fvalues(4),fvalues(5),fvalues(6),missing)
       end if
       deallocate(arguments,fvalues,cvalues)       
    case ("latitude_longitude")
       lon = x
       lat = y
    case default
       write(*,'(A,A,A)') "Geographic Projection Type ",trim(type)," not implemented."
       lon = missing
       lat = missing
    end select

  end subroutine proj_inv

  function getoptargs(string,arguments,fvalues,ivalues,cvalues) result (nargs)

    ! this subroutine separates a string composed of arguments and values separated by semicolons:
    ! arg1=0.24;arg2=222.0 into a character array of arguments and a real array of values
    ! integer and string values are optionally returned in separate arrays

    ! unallocated arguments variable : all arguments and values are returned
    ! allocated arguments variable : values are returned for arguments provided in the arguments array

    implicit none

    ! arguments
    character(len=*), intent(in) :: string
    character(len=*), allocatable, intent(inout) :: arguments(:)
    real(kind=4), allocatable, intent(inout) :: fvalues(:)
    integer, allocatable, intent(inout), optional :: ivalues(:)
    character(len=*), allocatable, intent(inout), optional :: cvalues(:)

    ! result
    integer :: nargs

    ! local variables
    integer :: a, pos1, pos2, pos3, nchars
    logical :: break
    logical :: all
    logical :: found

    nchars = len(string)

    if (allocated(arguments)) then
       ! read only values for arguments provided in 'arguments'
       all = .false.

       nargs = size(arguments)

       if (nargs.ne.0) then
          if (allocated(fvalues)) allocate(fvalues(nargs))
          allocate(fvalues(nargs))
          fvalues(:) = 0.0
          if (present(ivalues)) then
             if (allocated(ivalues)) deallocate(ivalues)
             allocate(ivalues(nargs))
             ivalues(:) = 0
          end if
          if (present(cvalues)) then
             if (allocated(cvalues)) deallocate(cvalues)
             allocate(cvalues(nargs))
             cvalues(:) = ""
          end if
       end if
    else
       ! read all values to all arguments in string
       all = .true.

       ! find number of arguments with values
       ! tokenize read string by semicolon-separated columns
       pos1 = 1
       a = 0
       break = .false.
       do while (.not.break)         
          pos2=index(string(pos1:nchars),';')
          if (pos2.gt.0) then
             pos2 = pos2 + pos1 - 2
          else
             pos2 = nchars
             break = .true.
          endif
          if (pos2.ge.pos1) then
             pos3=index(string(pos1:pos2),'=')
             if ((pos3.gt.1).and.(pos3.le.(pos2-pos1))) then
                a = a + 1
             end if
          end if
          pos1 = pos2+2
       end do

       nargs = a

       if (nargs.ne.0) then
          allocate(arguments(nargs))
          if (allocated(fvalues)) deallocate(fvalues)
          allocate(fvalues(nargs))
          fvalues(:) = 0.0
          if (present(ivalues)) then
             if (allocated(ivalues)) deallocate(ivalues)
             allocate(ivalues(nargs))
             ivalues(:) = 0
          end if
          if (present(cvalues)) then
             if (allocated(cvalues)) deallocate(cvalues)
             allocate(cvalues(nargs))
             cvalues(:) = ""
          end if
       end if
    end if

    if (nargs.ne.0) then
       ! read values from argument1=value1; argument2=value2 ... string
       ! tokenize read string by semicolon-separated columns
       pos1 = 1
       if (all) a = 0
       break = .false.
       do while (.not.break)         
          pos2=index(string(pos1:nchars),';')
          if (pos2.gt.0) then
             pos2 = pos2 + pos1 - 2
          else
             pos2 = nchars
             break = .true.
          endif
          if (pos2.ge.pos1) then
             pos3=index(string(pos1:pos2),'=')
             if ((pos3.gt.1).and.(pos3.le.(pos2-pos1))) then
                found = .false.
                if (all) then
                   a = a + 1
                   found = .true.
                else
                   a = size(arguments)
                   do while ((a.ne.0).and.(.not.found))
                      if (trim(arguments(a)).eq.trim(string(pos1:pos1+pos3-2))) then
                         found=.true.
                      else
                         a = a - 1
                      end if
                   end do
                end if
                if (found) then
                   arguments(a)=string(pos1:pos1+pos3-2)
                   if (scan(string(pos1+pos3:pos2),"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ") > 0) then
                      if (present(cvalues)) then
                         cvalues(a) = string(pos1+pos3:pos2)
                      end if
                   else
                      if (present(ivalues)) then
                         if (scan(string(pos1+pos3:pos2),".") > 0) then
                            read(string(pos1+pos3:pos2),*) fvalues(a)
                         else
                            read(string(pos1+pos3:pos2),*) ivalues(a)
                         end if
                      else
                         read(string(pos1+pos3:pos2),*) fvalues(a)
                      end if
                   end if
                end if
             end if
          end if
          pos1 = pos2+2
       end do
    end if

  end function getoptargs

  subroutine getxybounds(itype,ixmin,ixmax,iymin,iymax,idx,idy,ioptargs,otype,oxmin,oxmax,oymin,oymax,ooptargs,missing)

    ! this routine finds the outermost bounds of the input grid with given input projection in the given output projection

    implicit none

    ! arguments
    character(len=*), intent(in) :: itype
    real(kind=4), intent(in) :: ixmin
    real(kind=4), intent(in) :: ixmax
    real(kind=4), intent(in) :: iymin
    real(kind=4), intent(in) :: iymax
    real(kind=4), intent(in) :: idx
    real(kind=4), intent(in) :: idy
    character(len=*), intent(in) :: ioptargs
    character(len=*), intent(in) :: otype
    real(kind=4), intent(out) :: oxmin
    real(kind=4), intent(out) :: oxmax
    real(kind=4), intent(out) :: oymin
    real(kind=4), intent(out) :: oymax
    character(len=*), intent(in) :: ooptargs
    real(kind=4), intent(in) :: missing

    ! local variables
    integer :: nx, ny
    integer :: i,j
    real(kind=4), allocatable :: x(:,:), y(:,:)
    real(kind=4), allocatable :: lon(:,:), lat(:,:)
    logical, allocatable :: mask(:,:)
    real(kind=4) :: eps

    if ((trim(itype).ne.trim(otype)).or.(trim(ioptargs).ne.trim(ooptargs))) then

       ! determine smallest difference in 32 bit precision
       eps = 10.0*epsilon(0.0)

       nx = nint((ixmax-ixmin)/idx)
       ny = nint((iymax-iymin)/idy)

       allocate(x(nx*2,ny*2))
       allocate(y(nx*2,ny*2))
       allocate(lon(nx*2,ny*2))
       allocate(lat(nx*2,ny*2))
       allocate(mask(nx*2,ny*2))

       do i=1,nx
          x(i,:) = ixmin + idx*real(i-1)
          x(i+nx,:) = ixmin + idx*real(i)
       end do
       do j=1,ny
          y(:,j) = iymin + idy*real(j-1)
          y(:,j+ny) = iymin + idy*real(j)
       end do
       call proj_inv(itype,x,y,lon,lat,ioptargs,missing)
       call proj_fwd(otype,lon,lat,x,y,ooptargs,missing)

       mask = ((abs(x-missing).gt.eps).and.(abs(y-missing).gt.eps))

       if (count(mask).ne.0) then
          oxmin = minval(x,mask=mask)
          oxmax = maxval(x,mask=mask)
          oymin = minval(y,mask=mask)
          oymax = maxval(y,mask=mask)
       else
          oxmin = missing
          oxmax = missing
          oymin = missing
          oymax = missing
       end if

       deallocate(x,y,lon,lat,mask)
    else
       oxmin = ixmin
       oxmax = ixmax
       oymin = iymin
       oymax = iymax       
    end if

  end subroutine getxybounds

end module proj_module
