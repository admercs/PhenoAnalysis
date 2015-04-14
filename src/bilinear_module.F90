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

module bilinear_module

  ! 2009/02/19 Reto Stockli (Blue Marble Research)

  ! This module interpolates a regular grid bilinearly and its algorithm is largely based on
  ! the ESMF/CCSM bilinear interpolation routines

  ! 2009/12/01 all real variables are now dble (kind=8)

  ! 2014/12/22 separate interpolation for integer1/2/4/8 and float4/8

  ! Dependencies: 
  ! none

  private

  public :: bilineari1
  public :: bilineari2
  public :: bilineari4
  public :: bilineari8
  public :: bilinearf4
  public :: bilinearf8
  public :: find4corners
  public :: getweights
  public :: linearweights

contains

  subroutine bilineari1(input,output,pti,ptj,ptw,missing)

    implicit none

    ! arguments
    integer(kind=1),  intent(in) :: input(:,:)
    integer(kind=1),  intent(out) :: output(:,:)
    integer, intent(in) :: pti(:,:,:)
    integer, intent(in) :: ptj(:,:,:)
    real(kind=4),  intent(in) :: ptw(:,:,:)
    integer(kind=1), intent(in) :: missing

    ! local variables
    integer i,j
    integer ni,nj,np

    ni = size(pti,1)
    nj = size(pti,2)
    np = size(pti,3)

    do j = 1,nj
       do i = 1,ni
          if ((pti(i,j,1).ne.0).and.(pti(i,j,2).ne.0).and. &
               (pti(i,j,3).ne.0).and.(pti(i,j,4).ne.0)) then
             if ((input(pti(i,j,1),ptj(i,j,1)).ne.missing).and. &
                  (input(pti(i,j,2),ptj(i,j,2)).ne.missing).and. &
                  (input(pti(i,j,3),ptj(i,j,3)).ne.missing).and. &
                  (input(pti(i,j,4),ptj(i,j,4)).ne.missing)) then
                output(i,j) = nint(&
                     real(input(pti(i,j,1),ptj(i,j,1)),kind=4) * ptw(i,j,1) + &
                     real(input(pti(i,j,2),ptj(i,j,2)),kind=4) * ptw(i,j,2) + &
                     real(input(pti(i,j,3),ptj(i,j,3)),kind=4) * ptw(i,j,3) + &
                     real(input(pti(i,j,4),ptj(i,j,4)),kind=4) * ptw(i,j,4),kind=1)
             else
                output(i,j) = missing
             end if
          else
             output(i,j) = missing
          end if
       end do
    end do

  end subroutine bilineari1

  subroutine bilineari2(input,output,pti,ptj,ptw,missing)

    implicit none

    ! arguments
    integer(kind=2),  intent(in) :: input(:,:)
    integer(kind=2),  intent(out) :: output(:,:)
    integer, intent(in) :: pti(:,:,:)
    integer, intent(in) :: ptj(:,:,:)
    real(kind=4),  intent(in) :: ptw(:,:,:)
    integer(kind=2), intent(in) :: missing

    ! local variables
    integer i,j
    integer ni,nj,np

    ni = size(pti,1)
    nj = size(pti,2)
    np = size(pti,3)

    do j = 1,nj
       do i = 1,ni
          if ((pti(i,j,1).ne.0).and.(pti(i,j,2).ne.0).and. &
               (pti(i,j,3).ne.0).and.(pti(i,j,4).ne.0)) then
             if ((input(pti(i,j,1),ptj(i,j,1)).ne.missing).and. &
                  (input(pti(i,j,2),ptj(i,j,2)).ne.missing).and. &
                  (input(pti(i,j,3),ptj(i,j,3)).ne.missing).and. &
                  (input(pti(i,j,4),ptj(i,j,4)).ne.missing)) then
                output(i,j) = nint(&
                     real(input(pti(i,j,1),ptj(i,j,1)),kind=4) * ptw(i,j,1) + &
                     real(input(pti(i,j,2),ptj(i,j,2)),kind=4) * ptw(i,j,2) + &
                     real(input(pti(i,j,3),ptj(i,j,3)),kind=4) * ptw(i,j,3) + &
                     real(input(pti(i,j,4),ptj(i,j,4)),kind=4) * ptw(i,j,4),kind=2)
             else
                output(i,j) = missing
             end if
          else
             output(i,j) = missing
          end if
       end do
    end do

  end subroutine bilineari2

  subroutine bilineari4(input,output,pti,ptj,ptw,missing)

    implicit none

    ! arguments
    integer(kind=4),  intent(in) :: input(:,:)
    integer(kind=4),  intent(out) :: output(:,:)
    integer, intent(in) :: pti(:,:,:)
    integer, intent(in) :: ptj(:,:,:)
    real(kind=4),  intent(in) :: ptw(:,:,:)
    integer(kind=4), intent(in) :: missing

    ! local variables
    integer i,j
    integer ni,nj,np

    ni = size(pti,1)
    nj = size(pti,2)
    np = size(pti,3)

    do j = 1,nj
       do i = 1,ni
          if ((pti(i,j,1).ne.0).and.(pti(i,j,2).ne.0).and. &
               (pti(i,j,3).ne.0).and.(pti(i,j,4).ne.0)) then
             if ((input(pti(i,j,1),ptj(i,j,1)).ne.missing).and. &
                  (input(pti(i,j,2),ptj(i,j,2)).ne.missing).and. &
                  (input(pti(i,j,3),ptj(i,j,3)).ne.missing).and. &
                  (input(pti(i,j,4),ptj(i,j,4)).ne.missing)) then
                output(i,j) = nint(&
                     real(input(pti(i,j,1),ptj(i,j,1)),kind=4) * ptw(i,j,1) + &
                     real(input(pti(i,j,2),ptj(i,j,2)),kind=4) * ptw(i,j,2) + &
                     real(input(pti(i,j,3),ptj(i,j,3)),kind=4) * ptw(i,j,3) + &
                     real(input(pti(i,j,4),ptj(i,j,4)),kind=4) * ptw(i,j,4),kind=4)
             else
                output(i,j) = missing
             end if
          else
             output(i,j) = missing
          end if
       end do
    end do

  end subroutine bilineari4

  subroutine bilineari8(input,output,pti,ptj,ptw,missing)

    implicit none

    ! arguments
    integer(kind=8),  intent(in) :: input(:,:)
    integer(kind=8),  intent(out) :: output(:,:)
    integer, intent(in) :: pti(:,:,:)
    integer, intent(in) :: ptj(:,:,:)
    real(kind=4),  intent(in) :: ptw(:,:,:)
    integer(kind=8), intent(in) :: missing

    ! local variables
    integer i,j
    integer ni,nj,np

    ni = size(pti,1)
    nj = size(pti,2)
    np = size(pti,3)

    do j = 1,nj
       do i = 1,ni
          if ((pti(i,j,1).ne.0).and.(pti(i,j,2).ne.0).and. &
               (pti(i,j,3).ne.0).and.(pti(i,j,4).ne.0)) then
             if ((input(pti(i,j,1),ptj(i,j,1)).ne.missing).and. &
                  (input(pti(i,j,2),ptj(i,j,2)).ne.missing).and. &
                  (input(pti(i,j,3),ptj(i,j,3)).ne.missing).and. &
                  (input(pti(i,j,4),ptj(i,j,4)).ne.missing)) then
                output(i,j) = nint(&
                     real(input(pti(i,j,1),ptj(i,j,1)),kind=4) * ptw(i,j,1) + &
                     real(input(pti(i,j,2),ptj(i,j,2)),kind=4) * ptw(i,j,2) + &
                     real(input(pti(i,j,3),ptj(i,j,3)),kind=4) * ptw(i,j,3) + &
                     real(input(pti(i,j,4),ptj(i,j,4)),kind=4) * ptw(i,j,4),kind=8)
             else
                output(i,j) = missing
             end if
          else
             output(i,j) = missing
          end if
       end do
    end do

  end subroutine bilineari8

  subroutine bilinearf4(input,output,pti,ptj,ptw,missing)

    implicit none

    ! arguments
    real(kind=4),  intent(in) :: input(:,:)
    real(kind=4),  intent(out) :: output(:,:)
    integer, intent(in) :: pti(:,:,:)
    integer, intent(in) :: ptj(:,:,:)
    real(kind=4),  intent(in) :: ptw(:,:,:)
    real(kind=4), intent(in) :: missing

    ! local variables
    integer i,j
    integer ni,nj,np
    real(kind=4) :: eps

    ni = size(pti,1)
    nj = size(pti,2)
    np = size(pti,3)

    ! determine smallest number
    eps = 10.0*epsilon(0.0)

    do j = 1,nj
       do i = 1,ni
          if ((pti(i,j,1).ne.0).and.(pti(i,j,2).ne.0).and. &
               (pti(i,j,3).ne.0).and.(pti(i,j,4).ne.0)) then
             if ((abs(input(pti(i,j,1),ptj(i,j,1))-missing).gt.eps).and. &
                  (abs(input(pti(i,j,2),ptj(i,j,2))-missing).gt.eps).and. &
                  (abs(input(pti(i,j,3),ptj(i,j,3))-missing).gt.eps).and. &
                  (abs(input(pti(i,j,4),ptj(i,j,4))-missing).gt.eps)) then
                output(i,j) = &
                     input(pti(i,j,1),ptj(i,j,1)) * ptw(i,j,1) + &
                     input(pti(i,j,2),ptj(i,j,2)) * ptw(i,j,2) + &
                     input(pti(i,j,3),ptj(i,j,3)) * ptw(i,j,3) + &
                     input(pti(i,j,4),ptj(i,j,4)) * ptw(i,j,4)
             else
                output(i,j) = missing
             end if
          else
             output(i,j) = missing
          end if
       end do
    end do

  end subroutine bilinearf4

  subroutine bilinearf8(input,output,pti,ptj,ptw,missing)

    implicit none

    ! arguments
    real(kind=8),  intent(in) :: input(:,:)
    real(kind=8),  intent(out) :: output(:,:)
    integer, intent(in) :: pti(:,:,:)
    integer, intent(in) :: ptj(:,:,:)
    real(kind=4),  intent(in) :: ptw(:,:,:)
    real(kind=8), intent(in) :: missing

    ! local variables
    integer i,j
    integer ni,nj,np
    real(kind=8) :: eps

    ni = size(pti,1)
    nj = size(pti,2)
    np = size(pti,3)

    ! determine smallest number
    eps = 10.d0*epsilon(0.d0)

    do j = 1,nj
       do i = 1,ni
          if ((pti(i,j,1).ne.0).and.(pti(i,j,2).ne.0).and. &
               (pti(i,j,3).ne.0).and.(pti(i,j,4).ne.0)) then
             if ((abs(input(pti(i,j,1),ptj(i,j,1))-missing).gt.eps).and. &
                  (abs(input(pti(i,j,2),ptj(i,j,2))-missing).gt.eps).and. &
                  (abs(input(pti(i,j,3),ptj(i,j,3))-missing).gt.eps).and. &
                  (abs(input(pti(i,j,4),ptj(i,j,4))-missing).gt.eps)) then
                output(i,j) = &
                     input(pti(i,j,1),ptj(i,j,1)) * real(ptw(i,j,1),kind=8) + &
                     input(pti(i,j,2),ptj(i,j,2)) * real(ptw(i,j,2),kind=8) + &
                     input(pti(i,j,3),ptj(i,j,3)) * real(ptw(i,j,3),kind=8) + &
                     input(pti(i,j,4),ptj(i,j,4)) * real(ptw(i,j,4),kind=8)
             else
                output(i,j) = missing
             end if
          else
             output(i,j) = missing
          end if
       end do
    end do

  end subroutine bilinearf8

  ! Next Subroutines are taken from the CCSM software repository
  ! Bilinear and NN interpolation

  !===============================================================================
  !XXBOP ===========================================================================
  !
  ! !IROUTINE: shr_map_getWts -- local code that sets weights for a point
  !
  ! !DESCRIPTION:
  !     Local code that sets weights for a point.  Executes searches
  !     and computes weights.  For bilinear remap for example.
  !
  ! !REMARKS:
  !     Assumes Xsrc,Ysrc are regular lat/lon grids, monotonicallly increasing
  !        on constant latitude and longitude lines.  
  !     Assumes Xdst,Ydst,Xsrc,Ysrc are all either radians or degrees
  !
  ! !REVISION HISTORY:
  !     2005-Mar-27 - T. Craig - first version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine getweights(Xdst,Ydst,Xsrc,Ysrc,pti,ptj,ptw)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    !XXEOP

    real(kind=4)   ,intent(in) :: Xdst,Ydst
    real(kind=4)   ,intent(in) :: Xsrc(:,:),Ysrc(:,:)
    integer  ,intent(out):: pti(:),ptj(:)
    real(kind=4)   ,intent(out):: ptw(:)

    !--- local ---
    integer :: isize,jsize   ! array sizes
    integer :: il,ir         ! index of i left/mid/right
    integer :: jl,ju         ! index of j lower/mid/upper
    real(kind=4)   :: xsl,xsr       ! value of Xsrc, left/right
    real(kind=4)   :: ysl,ysu       ! value of Ysrc, left/right
    real(kind=4)   :: dx,dy,dx1,dy1 ! some d_lengths for weights calc


    !-------------------------------------------------------------------------------

    isize = size(Xsrc,1)
    jsize = size(Xsrc,2)

    if ((Xdst < Xsrc(1,1)).or.(Xdst > Xsrc(isize,1)).or.&
         (Ydst < Ysrc(1,1)).or.(Ydst > Ysrc(1,jsize))) then
       ! destination pixel is outside of input grid
       pti(:) = 0
       ptj(:) = 0
       ptw(:) = 0.0
    else
       ! destination pixel is within input grid

       call find4corners(Xdst,Ydst,Xsrc,Ysrc,il,ir,jl,ju)
       
       !--- bilinear ---
       xsl = Xsrc(il,1)
       xsr = Xsrc(ir,1)
       ysl = Ysrc(1,jl)
       ysu = Ysrc(1,ju)
       
       !--- compute dx1,dy1; distance from src(1) to dst
       dx  = (xsr-xsl)
       dy  = (ysu-ysl)
       dx1 = (Xdst-xsl)
       dy1 = (Ydst-Ysl)
       
       if ((dx < 0.0) .or. (dy < 0.0) .or. (dx1 > dx) .or. (dy1 > dy)) then
          write(*,*) ' '
          write(*,*) 'ERROR in dx,dy: ',dx1,dx,dy1,dy
          write(*,*) '   dst: ',Xdst,Ydst
          write(*,*) '   ind: ',il,ir,jl,ju
          write(*,*) '   dis: ',dx1,dx,dy1,dy
          write(*,*) '   x3 : ',xsl,Xdst,xsr
          write(*,*) '   y3 : ',ysl,Ydst,ysu
          write(*,*) ' '
          write(*,*) ' ERROR in dx,dy calc'
          stop
          return
       endif
       
       if (dx.gt.0.) dx1 = dx1 / dx
       if (dy.gt.0.) dy1 = dy1 / dy

       pti(1) = il
       pti(2) = ir
       pti(3) = il
       pti(4) = ir
       
       ptj(1) = jl
       ptj(2) = jl
       ptj(3) = ju
       ptj(4) = ju
       
       ptw(1) = (1.0-dx1)*(1.0-dy1)
       ptw(2) = (    dx1)*(1.0-dy1)
       ptw(3) = (1.0-dx1)*(    dy1)
       ptw(4) = (    dx1)*(    dy1)

    end if

  end subroutine getweights

  !===============================================================================

  subroutine find4corners(Xdst,Ydst,Xsrc,Ysrc,il,ir,jl,ju)

    ! finds 4 corner points surrounding dst in src
    ! returns left, right, lower, and upper i and j index

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    real(kind=4)  ,intent(in) :: Xdst,Ydst
    real(kind=4)  ,intent(in) :: Xsrc(:,:),Ysrc(:,:)
    integer,intent(out):: il,ir,jl,ju

    !--- local ---
    integer :: isize,jsize
    integer :: im,jm

    !-------------------------------------------------------------------------------

    isize = size(Xsrc,1)
    jsize = size(Xsrc,2)

    !--- find i index where Xsrc(i) <=  Xdst < Xsrc(i+1) ---
    il = 1
    ir = isize
    do while (ir-il > 1)
       im = (ir+il)/2
       if (Xdst >= Xsrc(im,1)) then
          il = im
       else
          ir = im
       endif
    enddo

    !--- find j index where Ysrc(j) <=  Ydst < Ysrc(j+1) ---
    jl = 1
    ju = jsize
    do while (ju-jl > 1)
       jm = (ju+jl)/2
       if (Ydst >= Ysrc(1,jm)) then
          jl = jm
       else
          ju = jm
       endif
    enddo
 
  end subroutine find4corners

  subroutine linearweights ( a, lo, hi, wgts )

    !*****************************************************************************
    ! this routine creates linear weights for array a across a number (wnum) of bins
    ! bounded by lo and hi

    implicit none

    ! arguments
    real(kind=4), intent(in) :: a(:)
    real(kind=4), intent(in) :: lo
    real(kind=4), intent(in) :: hi
    real(kind=4), intent(inout) :: wgts(:)

    ! local
    integer n
    integer i, j, i0, i1
    integer :: wnum
    real(kind=4) delta, fac
    real(kind=4), allocatable :: wbins(:)

    n = size(a)
    wnum = size(wgts)

    allocate(wbins(wnum))

    wgts(:) = 0.0

    delta = ( hi - lo ) / float(2 * wnum)

    do j = 1, wnum
       wbins(j) = delta + 2.0*delta*(j-1) + lo
    enddo

    do i = 1, n

       if ( a(i) < (lo + delta) ) then

          wgts(1) = wgts(1) + 1.0

       else if ( a(i) <= (hi - delta) ) then

          ! find left/right bounds
          j = wnum
          do while (wbins(j).gt.a(i))
             j=j-1
          enddo
          i0=j

          j = 1
          do while (wbins(j).lt.a(i))
             j=j+1
          enddo
          i1=j

          if (i0.eq.i1) then
             wgts(i0) = wgts(i0) + 1.0
          else
             fac = (a(i) - wbins(i0)) / (wbins(i1) - wbins(i0))
             wgts(i0) = wgts(i0) + 1.0-fac
             wgts(i1) = wgts(i1) + fac

          endif

       else if ( (hi - delta) < a(i) ) then

          wgts(wnum) = wgts(wnum) + 1.0

       end if

    end do

    wgts = wgts / sum(wgts)

    deallocate(wbins)

    return
  end subroutine linearweights

end module bilinear_module
