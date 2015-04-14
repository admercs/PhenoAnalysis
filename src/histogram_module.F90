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

! Code by Erin Sheldon
! https://github.com/esheldon/misc/blob/master/fortran/objshear/histogram.f90
module histogram_module

  use quicksort_module

  implicit none

  private

  public :: histf8
  public :: histf4
  public :: histi8
  public :: histi4
  public :: histi2
  public :: histi1

contains
    !   usage:
    !       histf8/f4/i8/i4/i2/i1(array, binsize, h, rev [, binmin, binmax, omin, omax])
    !   inputs:
    !       array: A 1-d array
    !       binsize: The bin size for the histogram
    !   optional inputs:
    !       binmin, binmax: The min,max values to consider. Default is
    !           the min and max of the input array
    !   outputs:
    !       h: The histogram
    !   optional outputs:
    !       rev: The reverse indices.
    !       omin, omax: the min,max values considered in histogram

  subroutine histf8(array, binsize, h, rev, binmin, binmax, omin, omax)

    ! arguments
    real(kind=8), intent(in) :: array(:)
    real(kind=8), intent(in) :: binsize

    integer, intent(inout), allocatable            :: h(:)
    integer, intent(inout), allocatable, optional  :: rev(:)

    real(kind=8), intent(in), optional :: binmin
    real(kind=8), intent(in), optional :: binmax
    real(kind=8), intent(out), optional :: omin
    real(kind=8), intent(out), optional :: omax

    ! local variables
    real(kind=8) :: mbinmin, mbinmax
    real(kind=8), allocatable :: temp_array(:)
    integer, allocatable :: sort_index(:)
    real(kind=8) :: bininv

    integer :: nbin
    logical :: dorev

    integer :: binnum, binnum_old, array_index, tbin, offset
    integer :: i, na

    na = size(array)

    ! generate sort indices
    allocate(temp_array(na))
    allocate(sort_index(na))
    temp_array = array
    sort_index = (/ (i, i=1,na) /)
    call quicksortf8(temp_array,sort_index,na)
    deallocate(temp_array)

    dorev = .false. 
    if (present(rev)) then
       dorev=.true.
    endif

    if (present(binmin)) then
       mbinmin = binmin
    else
       mbinmin = array(sort_index(1))
    endif
    if (present(binmax)) then
       mbinmax = binmax
    else
       mbinmax = array(sort_index(na))
    endif

    if (present(omin)) omin = mbinmin
    if (present(omax)) omax = mbinmax

    bininv = 1.0/binsize
    nbin = int( (mbinmax-mbinmin)*bininv) + 1

    ! allocate the outputs
    allocate(h(nbin))
    h=0
    if (dorev) then
       allocate(rev(na + nbin + 1))
       rev=size(rev)+1
    endif

    binnum_old = 0
    do i=1,na
       array_index = sort_index(i)

       ! offset into rev
       offset = i + nbin + 1
       if (dorev) then
          rev(offset) = array_index
       endif

       binnum = int( (array(array_index)-mbinmin)*bininv ) + 1

       if ( (binnum >= 1) .and. (binnum <= nbin) ) then

          ! should we update the reverse indices?
          if (dorev .and. (binnum > binnum_old)) then
             tbin = binnum_old + 1
             do while(tbin <= binnum)
                rev(tbin) = offset
                tbin = tbin + 1
             end do
          endif

          ! update the histogram
          h(binnum) = h(binnum) + 1
          binnum_old = binnum

       endif
    end do

    deallocate(sort_index)

  end subroutine histf8

  subroutine histf4(array, binsize, h, rev, binmin, binmax, omin, omax)

    ! arguments
    real(kind=4), intent(in) :: array(:)
    real(kind=4), intent(in) :: binsize

    integer, intent(inout), allocatable            :: h(:)
    integer, intent(inout), allocatable, optional  :: rev(:)

    real(kind=4), intent(in), optional :: binmin
    real(kind=4), intent(in), optional :: binmax
    real(kind=4), intent(out), optional :: omin
    real(kind=4), intent(out), optional :: omax

    ! local variables
    real(kind=4) :: mbinmin, mbinmax
    real(kind=4), allocatable :: temp_array(:)
    integer, allocatable :: sort_index(:)
    real(kind=4) :: bininv

    integer :: nbin
    logical :: dorev

    integer :: binnum, binnum_old, array_index, tbin, offset
    integer :: i, na

    na = size(array)

    ! generate sort indices
    allocate(temp_array(na))
    allocate(sort_index(na))
    temp_array = array
    sort_index = (/ (i, i=1,na) /)
    call quicksortf4(temp_array,sort_index,na)
    deallocate(temp_array)

    dorev = .false. 
    if (present(rev)) then
       dorev=.true.
    endif

    if (present(binmin)) then
       mbinmin = binmin
    else
       mbinmin = array(sort_index(1))
    endif
    if (present(binmax)) then
       mbinmax = binmax
    else
       mbinmax = array(sort_index(na))
    endif

    if (present(omin)) omin = mbinmin
    if (present(omax)) omax = mbinmax

    bininv = 1.0/binsize
    nbin = int( (mbinmax-mbinmin)*bininv) + 1

    ! allocate the outputs
    allocate(h(nbin))
    h=0
    if (dorev) then
       allocate(rev(na + nbin + 1))
       rev=size(rev)+1
    endif

    binnum_old = 0
    do i=1,na
       array_index = sort_index(i)

       ! offset into rev
       offset = i + nbin + 1
       if (dorev) then
          rev(offset) = array_index
       endif

       binnum = int( (array(array_index)-mbinmin)*bininv ) + 1

       if ( (binnum >= 1) .and. (binnum <= nbin) ) then

          ! should we update the reverse indices?
          if (dorev .and. (binnum > binnum_old)) then
             tbin = binnum_old + 1
             do while(tbin <= binnum)
                rev(tbin) = offset
                tbin = tbin + 1
             end do
          endif

          ! update the histogram
          h(binnum) = h(binnum) + 1
          binnum_old = binnum

       endif
    end do

    deallocate(sort_index)

  end subroutine histf4

  subroutine histi8(array, binsize, h, rev, binmin, binmax, omin, omax)

    ! arguments
    integer(kind=8), intent(in) :: array(:)
    integer(kind=8), intent(in) :: binsize

    integer, intent(inout), allocatable            :: h(:)
    integer, intent(inout), allocatable, optional  :: rev(:)

    integer(kind=8), intent(in), optional :: binmin
    integer(kind=8), intent(in), optional :: binmax
    integer(kind=8), intent(out), optional :: omin
    integer(kind=8), intent(out), optional :: omax

    ! local variables
    integer(kind=8) :: mbinmin, mbinmax
    integer(kind=8), allocatable :: temp_array(:)
    integer, allocatable :: sort_index(:)
    real(kind=8) :: bininv

    integer :: nbin
    logical :: dorev

    integer :: binnum, binnum_old, array_index, tbin, offset
    integer :: i, na

    na = size(array)

    ! generate sort indices
    allocate(temp_array(na))
    allocate(sort_index(na))
    temp_array = array
    sort_index = (/ (i, i=1,na) /)
    call quicksorti8(temp_array,sort_index,na)
    deallocate(temp_array)

    dorev = .false. 
    if (present(rev)) then
       dorev=.true.
    endif

    if (present(binmin)) then
       mbinmin = binmin
    else
       mbinmin = array(sort_index(1))
    endif
    if (present(binmax)) then
       mbinmax = binmax
    else
       mbinmax = array(sort_index(na))
    endif

    if (present(omin)) omin = mbinmin
    if (present(omax)) omax = mbinmax

    bininv = 1.d0/real(binsize,kind=8)
    nbin = nint( real(mbinmax-mbinmin,kind=8)*bininv) + 1

    ! allocate the outputs
    allocate(h(nbin))
    h=0
    if (dorev) then
       allocate(rev(na + nbin + 1))
       rev=size(rev)+1
    endif

    binnum_old = 0
    do i=1,na
       array_index = sort_index(i)

       ! offset into rev
       offset = i + nbin + 1
       if (dorev) then
          rev(offset) = array_index
       endif

       binnum = int( real(array(array_index)-mbinmin)*bininv ) + 1

       if ( (binnum >= 1) .and. (binnum <= nbin) ) then

          ! should we update the reverse indices?
          if (dorev .and. (binnum > binnum_old)) then
             tbin = binnum_old + 1
             do while(tbin <= binnum)
                rev(tbin) = offset
                tbin = tbin + 1
             end do
          endif

          ! update the histogram
          h(binnum) = h(binnum) + 1
          binnum_old = binnum

       endif
    end do

    deallocate(sort_index)

  end subroutine histi8

  subroutine histi4(array, binsize, h, rev, binmin, binmax, omin, omax)

    ! arguments
    integer(kind=4), intent(in), dimension(:) :: array
    integer(kind=4), intent(in)               :: binsize

    integer, intent(inout), dimension(:), allocatable            :: h
    integer, intent(inout), dimension(:), allocatable, optional  :: rev

    integer(kind=4), intent(in), optional :: binmin
    integer(kind=4), intent(in), optional :: binmax
    integer(kind=4), intent(out), optional :: omin
    integer(kind=4), intent(out), optional :: omax

    ! local variables
    integer(kind=4) :: mbinmin, mbinmax
    integer(kind=4), allocatable :: temp_array(:)
    integer, allocatable :: sort_index(:)
    real(kind=8) :: bininv
 
    integer :: nbin
    logical :: dorev

    integer :: binnum, binnum_old, array_index, tbin, offset
    integer :: i, na

    na = size(array)
    ! generate sort indices
    allocate(temp_array(na))
    allocate(sort_index(na))
    temp_array = array
    sort_index = (/ (i, i=1,na) /)
    call quicksorti4(temp_array,sort_index,na)
    deallocate(temp_array)

    dorev = .false. 
    if (present(rev)) then
       dorev=.true.
    endif

    if (present(binmin)) then
       mbinmin = binmin
    else
       mbinmin = array(sort_index(1))
    endif
    if (present(binmax)) then
       mbinmax = binmax
    else
       mbinmax = array(sort_index(na))
    endif

    if (present(omin)) omin = mbinmin
    if (present(omax)) omax = mbinmax

    bininv = 1.d0/real(binsize,kind=8)
    nbin = nint( real(mbinmax-mbinmin,kind=8)*bininv ) + 1

    ! allocate the outputs
    allocate(h(nbin))
    h=0
    if (dorev) then
       allocate(rev(na + nbin + 1))
       rev=size(rev)+1
    endif

    binnum_old = 0
    do i=1,na
       array_index = sort_index(i)

       ! offset into rev
       offset = i + nbin + 1
       if (dorev) then
          rev(offset) = array_index
       endif

       binnum = int( real(array(array_index)-mbinmin)*bininv ) + 1

       if ( (binnum >= 1) .and. (binnum <= nbin) ) then

          ! should we update the reverse indices?
          if (dorev .and. (binnum > binnum_old)) then
             tbin = binnum_old + 1
             do while(tbin <= binnum)
                rev(tbin) = offset
                tbin = tbin + 1
             end do
          endif

          ! update the histogram
          h(binnum) = h(binnum) + 1
          binnum_old = binnum

       endif
    end do

  end subroutine histi4

  subroutine histi2(array, binsize, h, rev, binmin, binmax, omin, omax)

    ! arguments
    integer(kind=2), intent(in), dimension(:) :: array
    integer(kind=2), intent(in)               :: binsize

    integer, intent(inout), dimension(:), allocatable            :: h
    integer, intent(inout), dimension(:), allocatable, optional  :: rev

    integer(kind=2), intent(in), optional :: binmin
    integer(kind=2), intent(in), optional :: binmax
    integer(kind=2), intent(out), optional :: omin
    integer(kind=2), intent(out), optional :: omax

    ! local variables
    integer(kind=2) :: mbinmin, mbinmax
    integer(kind=2), allocatable :: temp_array(:)
    integer, allocatable :: sort_index(:)
    real(kind=8) :: bininv
 
    integer :: nbin
    logical :: dorev

    integer :: binnum, binnum_old, array_index, tbin, offset
    integer :: i, na

    na = size(array)
    ! generate sort indices
    allocate(temp_array(na))
    allocate(sort_index(na))
    temp_array = array
    sort_index = (/ (i, i=1,na) /)
    call quicksorti2(temp_array,sort_index,na)
    deallocate(temp_array)

    dorev = .false. 
    if (present(rev)) then
       dorev=.true.
    endif

    if (present(binmin)) then
       mbinmin = binmin
    else
       mbinmin = array(sort_index(1))
    endif
    if (present(binmax)) then
       mbinmax = binmax
    else
       mbinmax = array(sort_index(na))
    endif

    if (present(omin)) omin = mbinmin
    if (present(omax)) omax = mbinmax

    bininv = 1.d0/real(binsize,kind=8)
    nbin = nint( real(mbinmax-mbinmin,kind=8)*bininv ) + 1

    ! allocate the outputs
    allocate(h(nbin))
    h=0
    if (dorev) then
       allocate(rev(na + nbin + 1))
       rev=size(rev)+1
    endif

    binnum_old = 0
    do i=1,na
       array_index = sort_index(i)

       ! offset into rev
       offset = i + nbin + 1
       if (dorev) then
          rev(offset) = array_index
       endif

       binnum = int( real(array(array_index)-mbinmin)*bininv) + 1

       if ( (binnum >= 1) .and. (binnum <= nbin) ) then

          ! should we update the reverse indices?
          if (dorev .and. (binnum > binnum_old)) then
             tbin = binnum_old + 1
             do while(tbin <= binnum)
                rev(tbin) = offset
                tbin = tbin + 1
             end do
          endif

          ! update the histogram
          h(binnum) = h(binnum) + 1
          binnum_old = binnum

       endif
    end do

  end subroutine histi2

  subroutine histi1(array, binsize, h, rev, binmin, binmax, omin, omax)

    ! arguments
    integer(kind=1), intent(in), dimension(:) :: array
    integer(kind=1), intent(in)               :: binsize

    integer, intent(inout), dimension(:), allocatable            :: h
    integer, intent(inout), dimension(:), allocatable, optional  :: rev

    integer(kind=1), intent(in), optional :: binmin
    integer(kind=1), intent(in), optional :: binmax
    integer(kind=1), intent(out), optional :: omin
    integer(kind=1), intent(out), optional :: omax

    ! local variables
    integer(kind=1) :: mbinmin, mbinmax
    integer(kind=1), allocatable :: temp_array(:)
    integer, allocatable :: sort_index(:)
    real(kind=8) :: bininv
 
    integer :: nbin
    logical :: dorev

    integer :: binnum, binnum_old, array_index, tbin, offset
    integer :: i, na

    na = size(array)
    ! generate sort indices
    allocate(temp_array(na))
    allocate(sort_index(na))
    temp_array = array
    sort_index = (/ (i, i=1,na) /)
    call quicksorti1(temp_array,sort_index,na)
    deallocate(temp_array)

    dorev = .false. 
    if (present(rev)) then
       dorev=.true.
    endif

    if (present(binmin)) then
       mbinmin = binmin
    else
       mbinmin = array(sort_index(1))
    endif
    if (present(binmax)) then
       mbinmax = binmax
    else
       mbinmax = array(sort_index(na))
    endif

    if (present(omin)) omin = mbinmin
    if (present(omax)) omax = mbinmax

    bininv = 1.d0/real(binsize,kind=8)
    nbin = nint( real(mbinmax-mbinmin,kind=8)*bininv ) + 1

    ! allocate the outputs
    allocate(h(nbin))
    h=0
    if (dorev) then
       allocate(rev(na + nbin + 1))
       rev=size(rev)+1
    endif

    binnum_old = 0
    do i=1,na
       array_index = sort_index(i)

       ! offset into rev
       offset = i + nbin + 1
       if (dorev) then
          rev(offset) = array_index
       endif

       binnum = int( real(array(array_index)-mbinmin)*bininv) + 1

       if ( (binnum >= 1) .and. (binnum <= nbin) ) then

          ! should we update the reverse indices?
          if (dorev .and. (binnum > binnum_old)) then
             tbin = binnum_old + 1
             do while(tbin <= binnum)
                rev(tbin) = offset
                tbin = tbin + 1
             end do
          endif

          ! update the histogram
          h(binnum) = h(binnum) + 1
          binnum_old = binnum

       endif
    end do

  end subroutine histi1

end module histogram_module
