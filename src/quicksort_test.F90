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

program quicksort_Test

  use quicksort_module

  implicit none

  ! ------------------------------------------------------------------------------
  ! Local variables
  ! ------------------------------------------------------------------------------

  integer, parameter :: gnx = 20
  integer, dimension(gnx) :: indarr
  real(kind=8), dimension(gnx) :: randvals = 0.0, randcopy = 0.0

  integer :: i

  indarr = (/ (i,i=1,gnx) /) ! indexical array

  call random_number(randvals)
  randcopy = randvals

  print *,' org_arr index'
  do i = 1,gnx
     print '(f10.4,1x,i3)',randvals(i),indarr(i)
  enddo

  call quicksortf8(randvals,indarr,gnx)

  print *,'       sorted      O_index'
  do i = 1,gnx
     print '(i8,1x,f10.4,1x,i8,1x,f10.4)',i,randvals(i),indarr(i),randcopy(indarr(i)) 
  enddo

  print *,' sorted - indexed original = ',sum(randvals - randcopy(indarr))

end program quicksort_Test
