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

module quicksort_module

! Code based on:
! A Comparative Study of Programming Languages in Rosetta Code
! Sebastian Nanz and Carlo A. Furia
! arXiv:1409.0252, 2014. 

  private

  public :: quicksortf8
  public :: quicksortf4
  public :: quicksorti8
  public :: quicksorti4
  public :: quicksorti2
  public :: quicksorti1

contains

  recursive subroutine quicksortf8(a,i,na)

    ! arguments
    integer, intent(in) :: na
    real(kind=8), intent(inout) :: a(na)
    integer, intent(inout) :: i(na)

    ! local variables
    integer :: left, right
    real(kind=4) :: random
    integer :: marker
    integer :: itemp

    real(kind=8) :: pivot
    real(kind=8) :: atemp

    if (na > 1) then

       call random_number(random)
       pivot = a(int(random*real(na-1))+1)   ! random pivot (not best performance, but avoids worst-case)
       left = 0
       right = na + 1

       do while (left < right)
          right = right - 1
          do while (a(right) > pivot)
             right = right - 1
          end do
          left = left + 1
          do while (a(left) < pivot)
             left = left + 1
          end do
          if (left < right) then
             atemp = a(left)
             a(left) = a(right)
             a(right) = atemp

             itemp = i(left)
             i(left) = i(right)
             i(right) = itemp
          end if
       end do

       if (left == right) then
          marker = left + 1
       else
          marker = left
       end if

       call quicksortf8(a(:marker-1),i(:marker-1),marker-1)
       call quicksortf8(a(marker:),i(marker:),na-marker+1)

    end if

  end subroutine quicksortf8

  recursive subroutine quicksortf4(a,i,na)

    ! arguments
    integer, intent(in) :: na
    real(kind=4), intent(inout) :: a(na)
    integer, intent(inout) :: i(na)

    ! local variables
    integer :: left, right
    real(kind=4) :: random
    integer :: marker
    integer :: itemp

    real(kind=4) :: pivot
    real(kind=4) :: atemp

    if (na > 1) then

       call random_number(random)
       pivot = a(int(random*real(na-1))+1)   ! random pivot (not best performance, but avoids worst-case)
       left = 0
       right = na + 1

       do while (left < right)
          right = right - 1
          do while (a(right) > pivot)
             right = right - 1
          end do
          left = left + 1
          do while (a(left) < pivot)
             left = left + 1
          end do
          if (left < right) then
             atemp = a(left)
             a(left) = a(right)
             a(right) = atemp

             itemp = i(left)
             i(left) = i(right)
             i(right) = itemp
          end if
       end do

       if (left == right) then
          marker = left + 1
       else
          marker = left
       end if

       call quicksortf4(a(:marker-1),i(:marker-1),marker-1)
       call quicksortf4(a(marker:),i(marker:),na-marker+1)

    end if

  end subroutine quicksortf4

  recursive subroutine quicksorti8(a,i,na)

    ! arguments
    integer, intent(in) :: na
    integer(kind=8), intent(inout) :: a(na)
    integer, intent(inout) :: i(na)

    ! local variables
    integer :: left, right
    real(kind=4) :: random
    integer :: marker
    integer :: itemp

    integer(kind=8) :: pivot
    integer(kind=8) :: atemp

    if (na > 1) then

       call random_number(random)
       pivot = a(int(random*real(na-1))+1)   ! random pivot (not best performance, but avoids worst-case)
       left = 0
       right = na + 1

       do while (left < right)
          right = right - 1
          do while (a(right) > pivot)
             right = right - 1
          end do
          left = left + 1
          do while (a(left) < pivot)
             left = left + 1
          end do
          if (left < right) then
             atemp = a(left)
             a(left) = a(right)
             a(right) = atemp

             itemp = i(left)
             i(left) = i(right)
             i(right) = itemp
          end if
       end do

       if (left == right) then
          marker = left + 1
       else
          marker = left
       end if

       call quicksorti8(a(:marker-1),i(:marker-1),marker-1)
       call quicksorti8(a(marker:),i(marker:),na-marker+1)

    end if

  end subroutine quicksorti8

  recursive subroutine quicksorti4(a,i,na)

    ! arguments
    integer, intent(in) :: na
    integer(kind=4), intent(inout) :: a(na)
    integer, intent(inout) :: i(na)

    ! local variables
    integer :: left, right
    real(kind=4) :: random
    integer :: marker
    integer :: itemp

    integer(kind=4) :: pivot
    integer(kind=4) :: atemp

    if (na > 1) then

       call random_number(random)
       pivot = a(int(random*real(na-1))+1)   ! random pivot (not best performance, but avoids worst-case)
       left = 0
       right = na + 1

       do while (left < right)
          right = right - 1
          do while (a(right) > pivot)
             right = right - 1
          end do
          left = left + 1
          do while (a(left) < pivot)
             left = left + 1
          end do
          if (left < right) then
             atemp = a(left)
             a(left) = a(right)
             a(right) = atemp

             itemp = i(left)
             i(left) = i(right)
             i(right) = itemp
          end if
       end do

       if (left == right) then
          marker = left + 1
       else
          marker = left
       end if

       call quicksorti4(a(:marker-1),i(:marker-1),marker-1)
       call quicksorti4(a(marker:),i(marker:),na-marker+1)

    end if

  end subroutine quicksorti4

  recursive subroutine quicksorti2(a,i,na)

    ! arguments
    integer, intent(in) :: na
    integer(kind=2), intent(inout) :: a(na)
    integer, intent(inout) :: i(na)

    ! local variables
    integer :: left, right
    real(kind=4) :: random
    integer :: marker
    integer :: itemp

    integer(kind=2) :: pivot
    integer(kind=2) :: atemp

    if (na > 1) then

       call random_number(random)
       pivot = a(int(random*real(na-1))+1)   ! random pivot (not best performance, but avoids worst-case)
       left = 0
       right = na + 1

       do while (left < right)
          right = right - 1
          do while (a(right) > pivot)
             right = right - 1
          end do
          left = left + 1
          do while (a(left) < pivot)
             left = left + 1
          end do
          if (left < right) then
             atemp = a(left)
             a(left) = a(right)
             a(right) = atemp

             itemp = i(left)
             i(left) = i(right)
             i(right) = itemp
          end if
       end do

       if (left == right) then
          marker = left + 1
       else
          marker = left
       end if

       call quicksorti2(a(:marker-1),i(:marker-1),marker-1)
       call quicksorti2(a(marker:),i(marker:),na-marker+1)

    end if

  end subroutine quicksorti2

  recursive subroutine quicksorti1(a,i,na)

    ! arguments
    integer, intent(in) :: na
    integer(kind=1), intent(inout) :: a(na)
    integer, intent(inout) :: i(na)

    ! local variables
    integer :: left, right
    real(kind=4) :: random
    integer :: marker
    integer :: itemp

    integer(kind=1) :: pivot
    integer(kind=1) :: atemp

    if (na > 1) then

       call random_number(random)
       pivot = a(int(random*real(na-1))+1)   ! random pivot (not best performance, but avoids worst-case)
       left = 0
       right = na + 1

       do while (left < right)
          right = right - 1
          do while (a(right) > pivot)
             right = right - 1
          end do
          left = left + 1
          do while (a(left) < pivot)
             left = left + 1
          end do
          if (left < right) then
             atemp = a(left)
             a(left) = a(right)
             a(right) = atemp

             itemp = i(left)
             i(left) = i(right)
             i(right) = itemp
          end if
       end do

       if (left == right) then
          marker = left + 1
       else
          marker = left
       end if

       call quicksorti1(a(:marker-1),i(:marker-1),marker-1)
       call quicksorti1(a(marker:),i(marker:),na-marker+1)

    end if

  end subroutine quicksorti1

end module quicksort_module
