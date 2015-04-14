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

module modis_module

  ! 2009/02/18 Reto Stockli (Blue Marble Research)

  ! This module contains the reading and gridding routines to read
  ! MODIS HDF4-eos tiles

  ! MODIS-ASCII support removed 2008/12/29
  ! MODIS binary tiles support removed 2008/12/29 (only HDF4 allowed at this point)

  ! DEPENDENCIES:
  ! driver_module.F90, hdf_module.F90, topo_module.F90, pft_module.F90

  implicit none

  private

  integer, parameter :: xtilekm = 1200  ! number of horizontal km per L3 MODIS tile
  integer, parameter :: ytilekm = 1200  ! number of vertical km per L3 MODIS tile

  type :: modis_type
     ! geographic
     real(kind=4)           :: xmin
     real(kind=4)           :: xmax
     real(kind=4)           :: ymin
     real(kind=4)           :: ymax
     real(kind=4)           :: dx
     real(kind=4)           :: dy
     integer                :: nx
     integer                :: ny
     character(len=100)     :: type
     character(len=100)     :: optargs
     real(kind=4)               :: resolution     ! pixel resolution [km]
     integer                    :: np        ! number of pixels
     integer                    :: hmin      ! minimum horizontal tile number
     integer                    :: hmax      ! maximum horizontal tile number
     integer                    :: vmin      ! minimum vertical tile number
     integer                    :: vmax      ! maximum vertical tile number
     real(kind=4)               :: smin      ! minimum horizontal sample number (at hmin)
     real(kind=4)               :: smax      ! maximum horizontal sample number (at hmax)
     real(kind=4)               :: lmin      ! minimum vertical line number (at vmin)
     real(kind=4)               :: lmax      ! maximum vertical line number (at vmax)
     integer, allocatable       :: h(:)      ! horizontal tile number
     integer, allocatable       :: v(:)      ! vertical tile number
     integer, allocatable       :: x(:)      ! horizontal pixel number within tile
     integer, allocatable       :: y(:)      ! vertical pixel number within tile
     ! physical
     character(len=100)         :: dir       ! HDF file directory
     character(len=50)          :: prefix    ! Product Name / Prefix
     character(len=50)          :: version   ! Product Version
     character(len=50)          :: name      ! variable name
     character(len=50)          :: hdfname   ! variable name (in HDF file)
     character(len=50)          :: qc1name   ! Quality Layer #1 name (in HDF file)
     character(len=50)          :: qc2name   ! Quality Layer #2 name (in HDF file)
     integer                    :: qc1bits   ! number of QC bits (Layer #1)
     integer                    :: qc2bits   ! number of QC bits (Layer #2)
     integer, allocatable       :: qc1error(:) ! add this value times minerror if QC bit is set (Layer #1)
     integer, allocatable       :: qc2error(:) ! add this value times minerror if QC bit is set (Layer #2)
     character(len=100)         :: longname  ! variable full name
     character(len=50)          :: units     ! variable units
     real(kind=4)               :: scale     ! scale factor (to convert MODIS DN to floating point value)
     real(kind=4)               :: offset    ! add offset factor (to convert MODIS DN to floating point value)
     real(kind=4)               :: magnitude ! order of magnitude for variable
     integer                    :: range(2)  ! valid range for MODIS product (DN)
     integer                    :: fill      ! fill value for modis product (DN)
     real(kind=4)               :: minerror  ! minimum error of specified MODIS product (arbitrary value)
     ! temporal
     integer                    :: nt     ! number of observations per year
     integer                    :: t0     ! first observation time step in year (seconds)
     integer                    :: dt     ! observation time step (seconds)
  end type modis_type

  type(modis_type), allocatable :: modis(:)

  integer :: nmodis               ! number of MODIS products employed in data assimilation

  public :: init_modis
  public :: exit_modis
  public :: read_modis

contains

  subroutine init_modis(directory,name,type,xmin,xmax,ymin,ymax,dx,dy,optargs,missing, &
           otype,oxmin,oxmax,oymin,oymax,odx,ody,ooptargs,onx,ony,ont,ot0,odt)

    use proj_module
    use modsin_module
    use parameter_module

    implicit none

    ! arguments
    character(len=*), intent(in) :: directory
    character(len=*), intent(in) :: name(:)
    character(len=*), intent(in) :: type
    real(kind=4), intent(in)     :: xmin
    real(kind=4), intent(in)     :: xmax
    real(kind=4), intent(in)     :: ymin
    real(kind=4), intent(in)     :: ymax
    real(kind=4), intent(in)     :: dx
    real(kind=4), intent(in)     :: dy
    character(len=*), intent(in) :: optargs
    real(kind=4), intent(in)     :: missing
    character(len=*), intent(inout) :: otype(:)
    real(kind=4), intent(inout)     :: oxmin(:)
    real(kind=4), intent(inout)     :: oxmax(:)
    real(kind=4), intent(inout)     :: oymin(:)
    real(kind=4), intent(inout)     :: oymax(:)
    real(kind=4), intent(inout)     :: odx(:)
    real(kind=4), intent(inout)     :: ody(:)
    character(len=*), intent(inout) :: ooptargs(:)
    integer, intent(inout)       :: onx(:)
    integer, intent(inout)       :: ony(:)
    integer, intent(inout)       :: ont(:)
    integer, intent(inout)       :: ot0(:)
    integer, intent(inout)       :: odt(:)
    

    ! local variables
    integer :: m
    
    nmodis = size(name)
    allocate(modis(nmodis))

    do m = 1,nmodis

       ! set static parameters for selected MODIS dataset
       select case (trim(name(m)))
       case ("FPAR")
          modis(m)%prefix = "MOD15A2"
          modis(m)%name = "FPAR"
          modis(m)%version = "005"
          modis(m)%resolution = 1.0
          modis(m)%type = "modsin"
          modis(m)%dx = 1.0
          modis(m)%dy = 1.0
          write(modis(m)%optargs,'(A,F0.3)') "resolution=",modis(m)%resolution
          modis(m)%hdfname = "Fpar_1km"
          modis(m)%qc1name = "FparLai_QC"
          modis(m)%qc2name = "FparExtra_QC"
          modis(m)%longname = "Fraction of Photosynthetically Absorbed Radiation"
          modis(m)%units = '-'
          modis(m)%qc1bits = 8
          modis(m)%qc2bits = 8
          allocate(modis(m)%qc1error(modis(m)%qc1bits))
          allocate(modis(m)%qc2error(modis(m)%qc2bits))
          modis(m)%qc1error = (/9,4,2,9,9,9,0,1/)
          modis(m)%qc2error = (/0,9,9,9,5,9,9,9/)
          modis(m)%scale = 0.01
          modis(m)%offset = 0.0
          modis(m)%magnitude = 1.0
          modis(m)%range = (/0,100/)
          modis(m)%fill = 255
          modis(m)%minerror = 0.05
          modis(m)%nt = 46
          modis(m)%dt = 8*86400
          modis(m)%t0 = 86400
       case ("LAI")
          modis(m)%prefix = "MOD15A2"
          modis(m)%name = "LAI"
          modis(m)%version = "005"
          modis(m)%resolution = 1.0
          modis(m)%type = "modsin"
          modis(m)%dx = 1.0
          modis(m)%dy = 1.0
          write(modis(m)%optargs,'(A,F0.3)') "resolution=",modis(m)%resolution
          modis(m)%hdfname = "Lai_1km"
          modis(m)%qc1name = "FparLai_QC"
          modis(m)%qc2name = "FparExtra_QC"
          modis(m)%longname = "Leaf Area Index"
          modis(m)%units = 'm2/m2'
          modis(m)%qc1bits = 8
          modis(m)%qc2bits = 8
          allocate(modis(m)%qc1error(modis(m)%qc1bits))
          allocate(modis(m)%qc2error(modis(m)%qc2bits))
          modis(m)%qc1error = (/9,4,2,9,9,9,0,1/)
          modis(m)%qc2error = (/0,9,9,9,5,9,9,9/)
          modis(m)%scale = 0.1
          modis(m)%offset = 0.0
          modis(m)%magnitude = 10.0
          modis(m)%range = (/0,100/)
          modis(m)%fill = 255
          modis(m)%minerror = 0.25
          modis(m)%nt = 46
          modis(m)%dt = 8*86400
          modis(m)%t0 = 86400
       case default
          write(*,'(A,A)') "Undefined MODIS dataset: ",trim(name(m))
          stop
       end select

       ! pre-calculate projection bounds
       call getxybounds(type,xmin,xmax,ymin,ymax,dx,dy,optargs,&
            modis(m)%type,modis(m)%xmin,modis(m)%xmax,modis(m)%ymin,modis(m)%ymax,modis(m)%optargs,missing)
       
       ! MODIS sample/line are integer numbers, their bounds are +/- 0.5 of the respective pixel centers.
       ! thus set the absolute bounds to +/- 0.5 of pixel centers
       if (((modis(m)%xmin-missing).gt.eps).and. &
            ((modis(m)%xmax-missing).gt.eps).and.&
            ((modis(m)%ymin-missing).gt.eps).and. &
            ((modis(m)%ymax-missing).gt.eps)) then
          modis(m)%xmin = nint(modis(m)%xmin-0.5*modis(m)%dx)+0.5*modis(m)%dx
          modis(m)%xmax = nint(modis(m)%xmax+0.5*modis(m)%dx)-0.5*modis(m)%dx
          modis(m)%ymin = nint(modis(m)%ymin-0.5*modis(m)%dy)+0.5*modis(m)%dy
          modis(m)%ymax = nint(modis(m)%ymax+0.5*modis(m)%dy)-0.5*modis(m)%dy
          
          if ((modis(m)%xmax-modis(m)%xmin).lt.modis(m)%dx) then
             modis(m)%xmax = modis(m)%xmin
             modis(m)%xmax = modis(m)%xmax + 0.5*modis(m)%dx
             modis(m)%xmin = modis(m)%xmin - 0.5*modis(m)%dx
          end if
          if ((modis(m)%ymax-modis(m)%ymin).lt.modis(m)%dy) then
             modis(m)%ymax = modis(m)%ymin
             modis(m)%ymax = modis(m)%ymax + 0.5*modis(m)%dy
             modis(m)%ymin = modis(m)%ymin - 0.5*modis(m)%dy
          end if
          
          modis(m)%nx = nint((modis(m)%xmax - modis(m)%xmin)/modis(m)%dx)
          modis(m)%ny = nint((modis(m)%ymax - modis(m)%ymin)/modis(m)%dy)
          
          call modsinxy2tile(modis(m)%xmin, modis(m)%ymin, modis(m)%hmin , modis(m)%vmin, &
               modis(m)%smin, modis(m)%lmin, modis(m)%resolution, missing)
          call modsinxy2tile(modis(m)%xmax, modis(m)%ymax, modis(m)%hmax , modis(m)%vmax, &
               modis(m)%smax, modis(m)%lmax, modis(m)%resolution, missing)
       
          modis(m)%dir = trim(directory)//trim(modis(m)%prefix)//"."//trim(modis(m)%version)
       else
          modis(m)%xmin = missing
          modis(m)%xmax = missing
          modis(m)%ymin = missing
          modis(m)%ymax = missing
          modis(m)%nx = 0
          modis(m)%ny = 0
          modis(m)%dir = ""
       end if

    end do

    otype = modis%type
    oxmin = modis%xmin
    oxmax = modis%xmax
    oymin = modis%ymin
    oymax = modis%ymax
    odx = modis%dx
    ody = modis%dy
    onx = modis%nx
    ony = modis%ny
    ooptargs = modis%optargs
    ont = modis%nt
    ot0 = modis%t0
    odt = modis%dt

  end subroutine init_modis

  subroutine exit_modis

    ! clean up modis variables

    implicit none

    ! local variables
    integer :: m

    do m=1,nmodis

       if (allocated(modis(m)%h)) then
          deallocate(modis(m)%h)
          deallocate(modis(m)%v)
          deallocate(modis(m)%x)
          deallocate(modis(m)%y)
          deallocate(modis(m)%qc1error)
          deallocate(modis(m)%qc2error)
       endif

    enddo

    deallocate(modis)

  end subroutine exit_modis

  subroutine read_modis(name,data,error,year,doy,missing,hasdata)

    use hdf_module
    use file_module

    implicit none

    ! Arguments
    character(len=*), intent(in) :: name      ! observation parameter name
    real(kind=4), intent(inout) :: data(:,:)  ! observation data  
    real(kind=4), intent(inout) :: error(:,:) ! observation error
    integer, intent(in) :: year               ! observation year
    integer, intent(in) :: doy                ! observation day-of-year
    real(kind=4), intent(in) :: missing        ! missing data value
    logical, intent(out) :: hasdata           ! check if any data was found at all

    ! Local Variables
    integer :: m,h,v,k,l,b
    integer :: x0,x1,y0,y1,nx,ny
    integer :: k0,k1,l0,l1

    character(len=256) :: ftemp
    integer :: sd_id

    character(len=1), allocatable :: val(:,:)
    character(len=1), allocatable :: qc1(:,:)
    character(len=1), allocatable :: qc2(:,:)
    integer(kind=4), allocatable :: byte(:,:)

    logical, parameter :: verbose = .true.

    character(len=7) :: datestring
    character(len=6) :: tile

    hasdata = .false.

    ! search and loop through modis tile h/v files
    do m=1,nmodis
       if (trim(modis(m)%name).eq.trim(name)) then

          if (size(data,1).ne.modis(m)%nx) then
             write(*,'(A,I6,I6)') "Requested and initialized MODIS domain size do not match: ",size(data,1),modis(m)%nx
             stop
          end if
          if (size(data,2).ne.modis(m)%ny) then
             write(*,'(A,I6,I6)') "Requested and initialized MODIS domain size do not match: ",size(data,2),modis(m)%ny
             stop
          end if

          allocate(val(modis(m)%nx,modis(m)%ny))
          allocate(qc1(modis(m)%nx,modis(m)%ny))
          allocate(qc2(modis(m)%nx,modis(m)%ny))

          val = achar(modis(m)%fill)
          qc1 = achar(modis(m)%fill)
          qc2 = achar(modis(m)%fill)

          do v=modis(m)%vmin,modis(m)%vmax
             do h=modis(m)%hmin,modis(m)%hmax

                write(unit=datestring,fmt='(i7)') 1000*year + doy

                write(unit=tile,fmt='(i6)') 1000*h + v
                if (h.lt.10) tile(2:2)='0'
                tile(1:1)='h'
                tile(4:4)='v'

                ftemp = trim(modis(m)%dir)//'/'//datestring//'/'//trim(modis(m)%prefix)//  &
                     '.A'//datestring//'.'//tile//'.*.hdf'

                ! search file based on above pattern
                if (find_file(ftemp)) then    

                   hasdata = .true.

!                   print*,modis(m)%smin+0.5*modis(m)%dx,modis(m)%smax-0.5*modis(m)%dx,&
!                        modis(m)%lmin+0.5*modis(m)%dy,modis(m)%lmax-0.5*modis(m)%dy

                   if (h.eq.modis(m)%hmin) then
                      x0 = nint(modis(m)%smin+0.5*modis(m)%dx)
                      k0 = 1
                   else
                      x0 = 0
                      k0 = xtilekm/nint(modis(m)%resolution) - nint(modis(m)%smin+0.5*modis(m)%dx) + &
                           (h-modis(m)%hmin-1)*xtilekm/nint(modis(m)%resolution) + 1
                   end if
                   if (h.eq.modis(m)%hmax) then
                      x1 = nint(modis(m)%smax-0.5*modis(m)%dx)
                      k1 = modis(m)%nx
                   else
                      x1 = xtilekm/nint(modis(m)%resolution) - 1
                      k1 = modis(m)%nx - nint(modis(m)%smax-0.5*modis(m)%dx) - 1 - &
                           (modis(m)%hmax-h-1)*xtilekm/nint(modis(m)%resolution)
                   end if

                   if (v.eq.modis(m)%vmin) then
                      y0 = nint(modis(m)%lmin+0.5*modis(m)%dy)
                      l0 = 1
                   else
                      y0 = 0
                      l0 = ytilekm/nint(modis(m)%resolution) - nint(modis(m)%lmin+0.5*modis(m)%dy) + &
                           (v-modis(m)%vmin-1)*ytilekm/nint(modis(m)%resolution) + 1
                   end if
                   if (v.eq.modis(m)%vmax) then
                      y1 = nint(modis(m)%lmax-0.5*modis(m)%dy)
                      l1 = modis(m)%ny
                   else
                      y1 = ytilekm/nint(modis(m)%resolution) - 1
                      l1 = modis(m)%ny - nint(modis(m)%lmax-0.5*modis(m)%dy) - 1 - &
                           (modis(m)%vmax-v-1)*ytilekm/nint(modis(m)%resolution)
                   end if

                   nx = x1-x0+1
                   ny = y1-y0+1

!                   print*,trim(modis(m)%hdfname)
!                   print*,h,v
!                   print*, nint(modis(m)%smin+0.5*modis(m)%dx),nint(modis(m)%smax-0.5*modis(m)%dx), &
!                        nint(modis(m)%lmin+0.5*modis(m)%dy),nint(modis(m)%lmax-0.5*modis(m)%dy)
!                   print*,x0,x1,y0,y1,nx,ny
!                   print*,k0,k1,l0,l1,modis(m)%nx,modis(m)%ny
!                   print*,

                   if ((x0.lt.0).or.(y0.lt.0).or.(x1.ge.xtilekm/nint(modis(m)%resolution)).or.&
                        (y1.ge.ytilekm/nint(modis(m)%resolution))) then
                      write(*,'(A,4I6)') "MODIS sample/line out of range: ",x0,x1,y0,y1
                      stop
                   end if
                   if ((k0.lt.1).or.(l0.lt.1).or.(k1.gt.modis(m)%nx).or.(l1.gt.modis(m)%ny)) then
                      write(*,'(A,4I6)') "OBS x/y out of range: ",k0,k1,l0,l1
                      stop
                   end if
                   if (nx.ne.(k1-k0+1)) then
                      write(*,'(A,I6,I6)') "MODIS subsetting in x direction does not match: ",nx,k1-k0+1
                      stop
                   end if
                   if (ny.ne.(l1-l0+1)) then
                      write(*,'(A,I6,I6)') "MODIS subsetting in y direction does not match: ",ny,l1-l0+1
                      stop
                   end if

                   call open_hdf_file(ftemp,sd_id)
                   call read_hdf_file(sd_id,modis(m)%hdfname,x0,y0,nx,ny,val_char=val(k0:k1,l0:l1))
                   call close_hdf_file(sd_id)
                   if (modis(m)%qc1name.ne.'') then
                      call open_hdf_file(ftemp,sd_id)
                      call read_hdf_file(sd_id,modis(m)%qc1name,x0,y0,nx,ny,val_char=qc1(k0:k1,l0:l1))
                      call close_hdf_file(sd_id)
                   endif
                   if (modis(m)%qc2name.ne.'') then
                      call open_hdf_file(ftemp,sd_id)
                      call read_hdf_file(sd_id,modis(m)%qc2name,x0,y0,nx,ny,val_char=qc2(k0:k1,l0:l1))
                      call close_hdf_file(sd_id)
                   endif
                else
                   if (verbose) write(*,'(A,A,A,I4,A,I3,A,I2.2,A,I2.2)') &
                        'No ',trim(modis(m)%prefix),' files found for (YYYY, DDD, H, V): ', &
                        year,' ',doy,' h',h,' v',v
                endif
             enddo
          enddo

          if (hasdata) then

             allocate(byte(modis(m)%nx,modis(m)%ny))

             error = modis(m)%minerror ! minimum error (10% of max error)

             ! bitwise screening QA flags (btests starts from LSB -> MSB with argument 0 .. modis(m)%qc1bits-1
             ! above in the bitmask, MSB is left and LSB is right
             if (modis(m)%qc1name.ne.'') then
                byte = iachar(qc1)
                do l=1,modis(m)%ny
                   do k=1,modis(m)%nx
                      do b=1,modis(m)%qc1bits
                         ! add punishment to error, error can be come more than 100% if many quality flags are set
                         if ((error(k,l).ne.missing).and.(modis(m)%qc1error(b).ne.-1)) then
                            if ( btest(byte(k,l),modis(m)%qc1bits-b) ) &
                                 error(k,l)=error(k,l)+float(modis(m)%qc1error(b))*modis(m)%minerror
                            if ((modis(m)%qc1error(b).eq.9).and.btest(byte(k,l),modis(m)%qc1bits-b) ) error(k,l)=missing
                         end if
                      end do
                      ! hard-code flagging value 65 (backup algorithm and other quality) as very bad quality
                      if (byte(k,l).eq.65) error(k,l) = missing
                    end do
                end do
             end if

             if (modis(m)%qc2name.ne.'') then
                byte = iachar(qc2)
                do l=1,modis(m)%ny
                   do k=1,modis(m)%nx
                      do b=1,modis(m)%qc2bits
                         ! add punishment to error, error can be come more than 100% if many quality flags are set
                         if ((error(k,l).ne.missing).and.(modis(m)%qc2error(b).ne.-1)) then
                            if ( btest(byte(k,l),modis(m)%qc2bits-b) ) &
                                 error(k,l)=error(k,l)+float(modis(m)%qc2error(b))*modis(m)%minerror
                            if ((modis(m)%qc2error(b).eq.9).and.btest(byte(k,l),modis(m)%qc2bits-b) ) error(k,l)=missing
                         endif
                      enddo
                   enddo
                end do
             endif

             ! process the data and fill observation vector
             byte = iachar(val)

             where ((byte.ge.modis(m)%range(1)).and.(byte.le.modis(m)%range(2)).and. &
                  (byte.ne.modis(m)%fill).and.(error.ne.missing))                
                data = (real(byte,kind=4) - modis(m)%offset)*modis(m)%scale
             elsewhere
                data = missing
                error = missing
             end where

!             print*,modis(m)%nx,modis(m)%ny
!             print*,data(:,1)
!             print*,
!             print*,data(:,modis(m)%ny)

!             stop

!             do l=1,modis(m)%ny
!                do k=1,modis(m)%nx
!                   print*,x0+k-1,y0+l-1,data(k,l),error(k,l)
!                end do
!             end do

             deallocate(byte)

          endif

          deallocate(val)
          deallocate(qc1)
          deallocate(qc2)

       endif

    enddo

  end subroutine read_modis

end module modis_module
