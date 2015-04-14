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

module reproject_module

  implicit none

  private

  integer, parameter :: maxprojsave = 100
  type :: projtype
     character(len=100)         :: itype, otype
     real(kind=4)               :: ixmin, ixmax, iymin, iymax, idx, idy 
     real(kind=4)               :: oxmin, oxmax, oymin, oymax, odx, ody
     character(len=100)         :: ioptargs, ooptargs
     integer                    :: inx, iny, onx, ony
     logical                    :: ibounds, obounds
     real(kind=4), allocatable  :: ixg(:,:), iyg(:,:)
     real(kind=4), allocatable  :: oxg(:,:), oyg(:,:)
     real(kind=4), allocatable  :: ixo(:,:), iyo(:,:)
     real(kind=4), allocatable  :: oxi(:,:), oyi(:,:)
     integer, allocatable       :: x_nn(:,:), y_nn(:,:)
     integer, allocatable       :: x_bil(:,:,:), y_bil(:,:,:)
     real(kind=4), allocatable  :: w_bil(:,:,:)
     integer, allocatable       :: rev(:)
     integer                    :: omin, omax
     real(kind=4)               :: scale
  end type projtype
  type(projtype) :: projsave(maxprojsave)

  public :: reproject

contains

  subroutine reproject(itype, ixmin, ixmax, iymin, iymax, idx, idy, ioptargs, &
       otype, oxmin, oxmax, oymin, oymax, odx, ody, ooptargs,  &
       missing, purge, bincount, &
       idataf4, odataf4, idataf8, odataf8, &
       idatai1, odatai1, idatai2, odatai2, &
       idatai4, odatai4, idatai8, odatai8, &
       missingf4, missingf8, missingi1, missingi2, missingi4, missingi8, &
       ixgrid, iygrid, oxgrid, oygrid, &
       ixoval, iyoval, oxival, oyival)

    use proj_module
    use histogram_module
    use bilinear_module

    implicit none

    ! arguments
    character(len=*), intent(in)           :: itype, otype
    real(kind=4), intent(in)               :: ixmin, ixmax, iymin, iymax, idx, idy 
    real(kind=4), intent(in)               :: oxmin, oxmax, oymin, oymax, odx, ody 
    character(len=*), intent(in)           :: ioptargs, ooptargs
    real(kind=4), intent(in)               :: missing
    logical, intent(in), optional          :: purge
    real(kind=4), intent(in), optional  :: idataf4(:,:)
    real(kind=8), intent(in), optional  :: idataf8(:,:)
    integer(kind=1), intent(in), optional  :: idatai1(:,:)
    integer(kind=2), intent(in), optional  :: idatai2(:,:)
    integer(kind=4), intent(in), optional  :: idatai4(:,:)
    integer(kind=8), intent(in), optional  :: idatai8(:,:)
    real(kind=4), intent(inout), allocatable, optional  :: odataf4(:,:)
    real(kind=8), intent(inout), allocatable, optional  :: odataf8(:,:)
    integer(kind=1), intent(inout), allocatable, optional  :: odatai1(:,:)
    integer(kind=2), intent(inout), allocatable, optional  :: odatai2(:,:)
    integer(kind=4), intent(inout), allocatable, optional  :: odatai4(:,:)
    integer(kind=8), intent(inout), allocatable, optional  :: odatai8(:,:)
    real(kind=4), intent(in), optional :: missingf4
    real(kind=8), intent(in), optional :: missingf8
    integer(kind=1), intent(in), optional :: missingi1
    integer(kind=2), intent(in), optional :: missingi2
    integer(kind=4), intent(in), optional :: missingi4
    integer(kind=8), intent(in), optional :: missingi8
    real(kind=4), intent(inout), allocatable, optional  :: ixgrid(:,:), iygrid(:,:)
    real(kind=4), intent(inout), allocatable, optional  :: oxgrid(:,:), oygrid(:,:)
    real(kind=4), intent(inout), allocatable, optional  :: ixoval(:,:), iyoval(:,:)
    real(kind=4), intent(inout), allocatable, optional  :: oxival(:,:), oyival(:,:)
    integer(kind=4), intent(inout), allocatable, optional :: bincount(:,:)

    ! local variables
    integer :: inx, iny
    integer :: onx, ony
    logical :: ibounds, obounds
    logical :: opt_purge
    logical :: found
    real(kind=4) :: eps
    real(kind=8) :: deps
    integer :: c,i,j,k,n,p

    real(kind=4) :: xbounds(4,1), ybounds(4,1)
    real(kind=4) :: lonbounds(4,1), latbounds(4,1)
    real(kind=4) :: oxibounds(4,1), oyibounds(4,1)
    real(kind=4) :: ixobounds(4,1), iyobounds(4,1)
    real(kind=4), allocatable :: lon(:,:), lat(:,:)

    real(kind=4), allocatable :: ixg(:,:), iyg(:,:)
    real(kind=4), allocatable :: oxg(:,:), oyg(:,:)
    real(kind=4), allocatable :: ixo(:,:), iyo(:,:)
    real(kind=4), allocatable :: oxi(:,:), oyi(:,:)

    integer, allocatable :: x_nn(:,:), y_nn(:,:)
    integer, allocatable :: x_bil(:,:,:), y_bil(:,:,:)
    real(kind=4), allocatable :: w_bil(:,:,:)

    integer, allocatable :: x_1d(:,:), y_1d(:,:), xy_2d(:,:)
    integer, allocatable :: rev(:), hist(:)
    integer :: omin, omax
    integer, allocatable :: ind(:)
    integer :: nind

    real(kind=4) :: scalex, scaley, scale
    real(kind=4), parameter :: scale_default = 1.0 ! use histogram binning (aggregating) and NN filling

    real(kind=4) :: dataf4
    real(kind=8) :: dataf8
    integer(kind=8) :: datai
    integer :: datatype


    ! determine smallest difference in 32 bit and 64bit precision
    eps = 10.0*epsilon(0.0)
    deps = 10.d0*epsilon(0.d0)

    ! assume that we want to keep the projection information
    if (present(purge)) then
       opt_purge = purge
    else
       opt_purge = .false.
    end if

    ! check for already stored projection information
    found = .false.
    p = 1
    do while ((p.le.maxprojsave).and.(.not.found))
       if (allocated(projsave(p)%ixg)) then
          if ((abs(ixmin-projsave(p)%ixmin).lt.eps).and.&
               (abs(ixmax-projsave(p)%ixmax).lt.eps).and.&
               (abs(iymin-projsave(p)%iymin).lt.eps).and.&
               (abs(iymax-projsave(p)%iymax).lt.eps).and.&
               (abs(idx-projsave(p)%idx).lt.eps).and.&
               (abs(idy-projsave(p)%idy).lt.eps)) then
             if ((abs(oxmin-projsave(p)%oxmin).lt.eps).and.&
                  (abs(oxmax-projsave(p)%oxmax).lt.eps).and.&
                  (abs(oymin-projsave(p)%oymin).lt.eps).and.&
                  (abs(oymax-projsave(p)%oymax).lt.eps).and.&
                  (abs(odx-projsave(p)%odx).lt.eps).and.&
                  (abs(ody-projsave(p)%ody).lt.eps)) then
                if ((itype.eq.projsave(p)%itype).and.&
                     (otype.eq.projsave(p)%otype).and.&
                     (ioptargs.eq.projsave(p)%ioptargs).and.&
                     (ooptargs.eq.projsave(p)%ooptargs)) then
                   found = .true.
                end if
             end if
          end if
       end if
       if (.not.found) then
          p = p + 1
       end if
    end do

    if (found) then

!       write(*,'(A)') "Stored projection information found."

       ibounds = projsave(p)%ibounds
       obounds = projsave(p)%obounds
       inx = projsave(p)%inx
       iny = projsave(p)%iny
       onx = projsave(p)%onx
       ony = projsave(p)%ony
       allocate(ixg(inx,iny))
       allocate(iyg(inx,iny))
       ixg = projsave(p)%ixg
       iyg = projsave(p)%iyg
       allocate(oxg(onx,ony))
       allocate(oyg(onx,ony))
       oxg = projsave(p)%oxg
       oyg = projsave(p)%oyg
       allocate(ixo(inx,iny))
       allocate(iyo(inx,iny))
       ixo = projsave(p)%ixo
       iyo = projsave(p)%iyo
       allocate(oxi(onx,ony))
       allocate(oyi(onx,ony))
       oxi = projsave(p)%oxi
       oyi = projsave(p)%oyi

       if (ibounds) then
          allocate(x_nn(onx,ony))
          allocate(y_nn(onx,ony))
          x_nn = projsave(p)%x_nn
          y_nn = projsave(p)%y_nn
          allocate(x_bil(onx,ony,4))
          allocate(y_bil(onx,ony,4))
          allocate(w_bil(onx,ony,4))
          x_bil = projsave(p)%x_bil
          y_bil = projsave(p)%y_bil
          w_bil = projsave(p)%w_bil
       end if

       if (obounds) then
          allocate(rev(size(projsave(p)%rev)))
          rev = projsave(p)%rev
          omin = projsave(p)%omin
          omax = projsave(p)%omax
       end if

       scale = projsave(p)%scale
    else

!       write(*,'(A)') "Creating new projection information."

       ! check if input bounds are specified
       if ((abs(ixmin-missing).gt.eps).and.(abs(ixmax-missing).gt.eps).and.&
            (abs(iymin-missing).gt.eps).and.(abs(iymax-missing).gt.eps)) then 
          ibounds = .true. 
       else 
          ibounds = .false.
          if (.not.present(ixgrid).or..not.present(iygrid)) then
             write(*,'(A)') "please supply ixgrid and iygrid if input bounds are set to missing"
             stop
          end if
       end if

       ! check if output bounds are specified
       if ((abs(oxmin-missing).gt.eps).and.(abs(oxmax-missing).gt.eps).and.&
            (abs(oymin-missing).gt.eps).and.(abs(oymax-missing).gt.eps)) then 
          obounds = .true. 
       else 
          obounds = .false.
          if (.not.present(oxgrid).or..not.present(oygrid)) then
             write(*,'(A)') "please supply oxgrid and oygrid if output bounds are set to missing"
             stop
          end if
       end if

       ! number of input grid elements
       if (ibounds) then
          if (idx.lt.eps) then 
             inx = 1 
          else 
             inx = max(nint((ixmax-ixmin)/idx),1)
          end if
          if (idy.lt.eps) then 
             iny = 1 
          else 
             iny = max(nint((iymax-iymin)/idy),1)
          end if
       else
          opt_purge = .true.
          inx = size(ixgrid,1)
          iny = size(ixgrid,2)
       end if

       ! number of output grid elements
       if (obounds) then
          if (odx.lt.eps) then 
             onx = 1 
          else 
             onx = max(nint((oxmax-oxmin)/odx),1)
          end if
          if (ody.lt.eps) then 
             ony = 1 
          else 
             ony = max(nint((oymax-oymin)/ody),1)
          end if
       else
          opt_purge = .true.
          onx = size(oxgrid,1)
          ony = size(oxgrid,2)
       end if

       allocate(ixg(inx,iny))
       allocate(iyg(inx,iny))
       if (ibounds) then
          ! generate regular input grid coordinates
          do i=1,inx
             ixg(i,:) = ixmin + (real(i-1,kind=4) + 0.5)*idx
          end do
          do j=1,iny
             iyg(:,j) = iymin + (real(j-1,kind=4) + 0.5)*idy
          end do
       else 
          ! irregular input coordinates are provided in ixgrid and iygrid
          ixg = ixgrid
          iyg = iygrid
       end if

       allocate(oxg(onx,ony))
       allocate(oyg(onx,ony))
       if (obounds) then
          ! generate regular output grid coordinates
          do i=1,onx
             oxg(i,:) = oxmin + (real(i-1,kind=4) + 0.5)*odx
          end do
          do j=1,ony
             oyg(:,j) = oymin + (real(j-1,kind=4) + 0.5)*ody
          end do
       else
          ! irregular output coordinates are provided in oxgrid and oygrid
          oxg = oxgrid
          oyg = oygrid
       end if

       allocate(ixo(inx,iny))
       allocate(iyo(inx,iny))
       allocate(oxi(onx,ony))
       allocate(oyi(onx,ony))

       ! 1: reprojection of input grid to input grid in lon/lat coordinates
       ! 2: reprojection of input grid in lon/lat to input grid in output coordinates
       ! 3: reprojection of output grid to output grid in lon/lat coordinates
       ! 4: reprojection of output grid in lon/lat coordinates to output grid in input coordinates
       ! projection: input   -> lon/lat
       allocate(lon(inx,iny))
       allocate(lat(inx,iny))
       if (itype.eq."latitude_longitude") then
          lon = ixg
          lat = iyg
          if (ibounds) then
             lonbounds(:,1) = (/ixmin,ixmin,ixmax,ixmax/)
             latbounds(:,1) = (/iymin,iymax,iymin,iymax/)
          else
             lonbounds = missing
             latbounds = missing
          end if
       else
          call proj_inv(itype,ixg,iyg,lon,lat,ioptargs,missing)
          if (ibounds) then
             xbounds(:,1) = (/ixmin,ixmin,ixmax,ixmax/)
             ybounds(:,1) = (/iymin,iymax,iymin,iymax/)
             call proj_inv(itype,xbounds,ybounds,lonbounds,latbounds,ioptargs,missing)
          else
             lonbounds = missing
             latbounds = missing
          end if
       end if

       ! projection: lon/lat -> output
       if (otype.eq."latitude_longitude") then
          ixo = lon
          iyo = lat
          ixobounds = lonbounds
          iyobounds = latbounds
       else
          call proj_fwd(otype,lon,lat,ixo,iyo,ooptargs,missing)
          call proj_fwd(otype,lonbounds,latbounds,ixobounds,iyobounds,ooptargs,missing)
       end if

       ! projection: output  -> lon/lat
       deallocate(lon,lat)
       allocate(lon(onx,ony))
       allocate(lat(onx,ony))
       if (otype.eq."latitude_longitude") then
          lon = oxg
          lat = oyg
          if (obounds) then
             lonbounds(:,1) = (/oxmin,oxmin,oxmax,oxmax/)
             latbounds(:,1) = (/oymin,oymax,oymin,oymax/)
          else
             lonbounds = missing
             latbounds = missing
          end if
       else
          call proj_inv(otype,oxg,oyg,lon,lat,ooptargs,missing)
          if (obounds) then
             xbounds(:,1) = (/oxmin,oxmin,oxmax,oxmax/)
             ybounds(:,1) = (/oymin,oymax,oymin,oymax/)
             call proj_inv(otype,xbounds,ybounds,lonbounds,latbounds,ooptargs,missing)
          else
             lonbounds = missing
             latbounds = missing
          end if
       end if

       ! projection: lon/lat -> input
       if (itype.eq."latitude_longitude") then
          oxi = lon
          oyi = lat
          oxibounds = lonbounds
          oyibounds = latbounds
       else
          call proj_fwd(itype,lon,lat,oxi,oyi,ioptargs,missing)
          call proj_fwd(itype,lonbounds,latbounds,oxibounds,oyibounds,ooptargs,missing)
       end if

       deallocate(lon,lat)

       ! generate nearest neighbor indices
       if (ibounds) then
          allocate(x_nn(onx,ony))
          allocate(y_nn(onx,ony))
          if (idx.lt.eps) then
             x_nn = int(oxi - ixmin) + 1
          else
             x_nn = int((oxi - ixmin) / idx) + 1
          end if

          where ((x_nn < 1).or.(x_nn > inx))
             x_nn = 0
          end where

          if (idy.lt.eps) then
             y_nn = int(oyi - iymin) + 1
          else
             y_nn = int((oyi - iymin) / idy) + 1
          end if

          where ((y_nn < 1).or.(y_nn > iny))
             y_nn = 0
          end where
       end if

       ! generate bilinear interpolation indices and weights
       if (ibounds) then
          allocate(x_bil(onx,ony,4))
          allocate(y_bil(onx,ony,4))
          allocate(w_bil(onx,ony,4))
          do j=1,ony
             do i=1,onx
                call getweights(oxi(i,j),oyi(i,j),ixg,iyg,x_bil(i,j,:),y_bil(i,j,:),w_bil(i,j,:))
             end do
          end do
       end if

       ! generate area-aggregating histogram bins
       if (obounds) then
          allocate(x_1d(inx,iny))
          allocate(y_1d(inx,iny))
          allocate(xy_2d(inx,iny))
          if (odx.lt.eps) then
             x_1d = int(ixo - oxmin) + 1
          else
             x_1d = int((ixo - oxmin) / odx) + 1
          end if

          where ((x_1d < 1).or.(x_1d > onx))
             x_1d = 0
          end where

          if (ody.lt.eps) then
             y_1d = int(iyo - oymin) + 1
          else
             y_1d = int((iyo - oymin) / ody) + 1
          end if

          where ((y_1d < 1).or.(y_1d > ony))
             y_1d = 0
          end where

          where ((x_1d.ne.0).and.(y_1d.ne.0))
             xy_2d = (y_1d-1) * onx + x_1d
          elsewhere
             xy_2d = 0
          end where

          deallocate(x_1d,y_1d)

          call histi4(reshape(xy_2d,(/inx*iny/)), 1, hist, rev, binmin=1, omin=omin, omax=omax)

          deallocate(xy_2d)
       end if

       ! generate triangulation indices

       ! evaluate if the output grid is subgrid-scale to the input
       ! grid. This will affect the choice of reprojection
       ! method. We use bilinear interpolation for scale < 0.5.
       ! We use triangulation for scale 0.5 <= scale <= 1.5 and
       ! we use histogram binning with NN filling for scale >
       ! 1.5. If the output grid or the input grids are singular, then
       ! bilinear interpolation is used.

       ! update: triangulation is not yet implemented. We use histogram
       ! binning and NN filling for scale > 0.5

       if ((odx.lt.eps).or.(ody.lt.eps).or.(idx.lt.eps).or.(idy.lt.eps)) then
          scale = scale_default
       else
          if ((abs(oxibounds(1,1)-missing).gt.eps).and.(abs(oxibounds(2,1)-missing).gt.eps).and.&
               (abs(oxibounds(3,1)-missing).gt.eps).and.(abs(oxibounds(4,1)-missing).gt.eps).and.&                
               (abs(ixmin-missing).gt.eps).and.(abs(ixmax-missing).gt.eps)) then
             scalex = 0.5 * (oxibounds(3,1) - oxibounds(1,1) + oxibounds(4,1) - oxibounds(2,1)) / real(onx,kind=4) / &
                  ((ixmax - ixmin) / real(inx,kind=4))
          else
             scalex = scale_default
          end if
          if ((abs(oyibounds(1,1)-missing).gt.eps).and.(abs(oyibounds(2,1)-missing).gt.eps).and.&
               (abs(oyibounds(3,1)-missing).gt.eps).and.(abs(oyibounds(4,1)-missing).gt.eps).and.&                
               (abs(iymin-missing).gt.eps).and.(abs(iymax-missing).gt.eps)) then
             scaley = 0.5 * (oyibounds(2,1) - oyibounds(1,1) + oyibounds(4,1) - oyibounds(3,1)) / real(ony,kind=4) / &
                  ((iymax - iymin) / real(iny,kind=4))
          else
             scaley = scale_default
          end if
          scale = 0.5*(abs(scalex) + abs(scaley))
       end if

!       print*,"SCALE: ",scale

    end if

    ! reproject input data to output projection
    if ((present(idataf4).and.present(odataf4)).or. &
         (present(idataf8).and.present(odataf8)).or. &
         (present(idatai1).and.present(odatai1)).or. &
         (present(idatai2).and.present(odatai2)).or. &
         (present(idatai4).and.present(odatai4)).or. &
         (present(idatai8).and.present(odatai8))) then

       if (present(idataf4)) datatype = 1
       if (present(idataf8)) datatype = 2
       if (present(idatai1)) datatype = 3
       if (present(idatai2)) datatype = 4
       if (present(idatai4)) datatype = 5
       if (present(idatai8)) datatype = 6

       if (present(bincount)) then
          if (allocated(bincount)) then
             if ((size(bincount,1).ne.onx).or.(size(bincount,2).ne.ony)) then 
                deallocate(bincount)
                allocate(bincount(onx,ony))
             end if
          else
             allocate(bincount(onx,ony))
          end if
          bincount = 0
       end if

       select case (datatype)
       case (1)
          if (allocated(odataf4)) then
             if ((size(odataf4,1).ne.onx).or.(size(odataf4,2).ne.ony)) then 
                deallocate(odataf4)
                allocate(odataf4(onx,ony))
             end if
          else
             allocate(odataf4(onx,ony))
          end if
          odataf4 = missingf4
       case (2)
          if (allocated(odataf8)) then
             if ((size(odataf8,1).ne.onx).or.(size(odataf8,2).ne.ony)) then 
                deallocate(odataf8)
                allocate(odataf8(onx,ony))
             end if
          else
             allocate(odataf8(onx,ony))
          end if
          odataf8 = missingf8
       case (3)
          if (allocated(odatai1)) then
             if ((size(odatai1,1).ne.onx).or.(size(odatai1,2).ne.ony)) then 
                deallocate(odatai1)
                allocate(odatai1(onx,ony))
             end if
          else
             allocate(odatai1(onx,ony))
          end if
          odatai1 = missingi1
       case (4)
          if (allocated(odatai2)) then
             if ((size(odatai2,1).ne.onx).or.(size(odatai2,2).ne.ony)) then 
                deallocate(odatai2)
                allocate(odatai2(onx,ony))
             end if
          else
             allocate(odatai2(onx,ony))
          end if
          odatai2 = missingi2
       case (5)
          if (allocated(odatai4)) then
             if ((size(odatai4,1).ne.onx).or.(size(odatai4,2).ne.ony)) then 
                deallocate(odatai4)
                allocate(odatai4(onx,ony))
             end if
          else
             allocate(odatai4(onx,ony))
          end if
          odatai4 = missingi4
       case (6)
          if (allocated(odatai8)) then
             if ((size(odatai8,1).ne.onx).or.(size(odatai8,2).ne.ony)) then 
                deallocate(odatai8)
                allocate(odatai8(onx,ony))
             end if
          else
             allocate(odatai8(onx,ony))
          end if
          odatai8 = missingi8
       end select

       if (scale.gt.0.5) then

          ! fill result with NN values
          if (allocated(x_nn).and.allocated(y_nn)) then
             do j=1,ony
                do i=1,onx
                   if ((x_nn(i,j).ne.0).and.(y_nn(i,j).ne.0)) then
                      select case (datatype)
                      case (1) 
                         odataf4(i,j) = idataf4(x_nn(i,j),y_nn(i,j))
                         if (present(bincount)) then
                            if (abs(odataf4(i,j)-missingf4).gt.eps) bincount(i,j) = 1
                         end if
                      case (2) 
                         odataf8(i,j) = idataf8(x_nn(i,j),y_nn(i,j))
                         if (present(bincount)) then
                            if (abs(odataf8(i,j)-missingf8).gt.deps) bincount(i,j) = 1
                         end if
                      case (3) 
                         odatai1(i,j) = idatai1(x_nn(i,j),y_nn(i,j))
                         if (present(bincount)) then
                            if (odatai1(i,j).ne.missingi1) bincount(i,j) = 1
                         end if
                      case (4) 
                         odatai2(i,j) = idatai2(x_nn(i,j),y_nn(i,j))
                         if (present(bincount)) then
                            if (odatai2(i,j).ne.missingi2) bincount(i,j) = 1
                         end if
                       case (5) 
                         odatai4(i,j) = idatai4(x_nn(i,j),y_nn(i,j))
                         if (present(bincount)) then
                            if (odatai4(i,j).ne.missingi4) bincount(i,j) = 1
                         end if
                      case (6) 
                         odatai8(i,j) = idatai8(x_nn(i,j),y_nn(i,j))
                         if (present(bincount)) then
                            if (odatai8(i,j).ne.missingi8) bincount(i,j) = 1
                         end if
                      end select
                   end if
                end do
             end do
          end if

          ! sort through indices and find predominant pixel among selected
          ! indices for each output pixels. Fill predominant pixels
          ! into result where available (e.g. in areas with dense
          ! coverage of input data).
          if (allocated(rev)) then
!             print*,omin,omax,inx,iny
             do k=1,omax-omin+1
                if (rev(k).ne.rev(k+1)) then
                   nind=rev(k+1)-rev(k)
                   allocate(ind(nind))
                   ind=rev(rev(k):rev(k+1)-1)
                   c=0
                   dataf4 = 0.0
                   dataf8 = 0.d0
                   datai = 0
                   do n=1,nind
                      j = (ind(n)-1)/inx + 1
                      i = ind(n) - (j-1)*inx  
                      select case (datatype)
                      case (1)
                         if (abs(idataf4(i,j)-missingf4).gt.eps) then
                            dataf4 = dataf4 + idataf4(i,j)
                            c = c + 1
                         end if
                      case (2)
                         if (abs(idataf8(i,j)-missingf8).gt.deps) then
                            dataf8 = dataf8 + idataf8(i,j)
                            c = c + 1
                         end if
                      case (3)
                         if (idatai1(i,j).ne.missingi1) then
                            datai = datai + int(idatai1(i,j),kind=8)
                            c = c + 1
                         end if
                      case (4)
                         if (idatai2(i,j).ne.missingi2) then
                            datai = datai + int(idatai2(i,j),kind=8)
                            c = c + 1
                         end if
                      case (5)
                         if (idatai4(i,j).ne.missingi4) then
                            datai = datai + int(idatai4(i,j),kind=8)
                            c = c + 1
                         end if
                      case (6)
                         if (idatai8(i,j).ne.missingi8) then
                            datai = datai + int(idatai8(i,j),kind=8)
                            c = c + 1
                         end if
                      end select
                   end do
                   deallocate(ind)
                   if (c.ne.0) then
                      j = (k-2+omin)/onx + 1
                      i = k-1+omin - (j-1)*onx
                      select case (datatype)
                      case (1)
                         odataf4(i,j) = dataf4 / real(c,kind=4)
                      case (2)
                         odataf8(i,j) = dataf8 / real(c,kind=8)
                      case (3)
                         odatai1(i,j) = int(datai/c,kind=1)
                      case (4)
                         odatai2(i,j) = int(datai/c,kind=2)
                      case (5)
                         odatai4(i,j) = int(datai/c,kind=4)
                      case (6)
                         odatai8(i,j) = int(datai/c,kind=8)
                      end select
                      if (present(bincount)) bincount(i,j) = c
                   end if
                end if
             end do
          end if

       else

          ! do bilinear interpolation for downscaling
          if (allocated(x_bil).and.allocated(y_bil).and.allocated(w_bil)) then
             select case (datatype)
             case (1) 
                call bilinearf4(idataf4,odataf4,x_bil,y_bil,w_bil,missingf4)
                if (present(bincount)) then
                   where(abs(odataf4-missingf4).gt.eps) bincount = 4
                end if
             case (2) 
                call bilinearf8(idataf8,odataf8,x_bil,y_bil,w_bil,missingf8)
                if (present(bincount)) then
                   where(abs(odataf8-missingf8).gt.deps) bincount = 4
                end if
             case (3) 
                call bilineari1(idatai1,odatai1,x_bil,y_bil,w_bil,missingi1)
                if (present(bincount)) then
                   where(odatai1.ne.missingi1) bincount = 4
                end if
             case (4) 
                call bilineari2(idatai2,odatai2,x_bil,y_bil,w_bil,missingi2)
                if (present(bincount)) then
                   where(odatai2.ne.missingi2) bincount = 4
                end if
             case (5) 
                call bilineari4(idatai4,odatai4,x_bil,y_bil,w_bil,missingi4)
                if (present(bincount)) then
                   where(odatai4.ne.missingi4) bincount = 4
                end if
             case (6) 
                call bilineari8(idatai8,odatai8,x_bil,y_bil,w_bil,missingi8)
                if (present(bincount)) then
                   where(odatai8.ne.missingi8) bincount = 4
                end if
             end select
          end if

          ! quick fix on boundaries (where bilinear interpolation would require to 
          ! extend boundaries by 1 grid point : fill result with NN values
          if (allocated(x_nn).and.allocated(y_nn)) then
             do j=1,ony
                do i=1,onx
                   if ((x_nn(i,j).ne.0).and.(y_nn(i,j).ne.0)) then
                      select case (datatype)
                      case (1) 
                         if (abs(odataf4(i,j)-missingf4).le.eps) then
                            odataf4(i,j) = idataf4(x_nn(i,j),y_nn(i,j))
                            if (present(bincount)) then
                               if (abs(odataf4(i,j)-missingf4).gt.eps) bincount(i,j) = 1
                            end if
                         end if
                      case (2) 
                         if (abs(odataf4(i,j)-missingf8).le.deps) then
                            odataf8(i,j) = idataf8(x_nn(i,j),y_nn(i,j))
                            if (present(bincount)) then
                               if (abs(odataf8(i,j)-missingf8).gt.deps) bincount(i,j) = 1
                            end if
                         end if
                      case (3) 
                         if (odatai1(i,j).ne.missingi1) then
                            odatai1(i,j) = idatai1(x_nn(i,j),y_nn(i,j))
                            if (present(bincount)) then
                               if (odatai1(i,j).ne.missingi1) bincount(i,j) = 1
                            end if
                         end if
                      case (4) 
                         if (odatai1(i,j).ne.missingi2) then
                            odatai2(i,j) = idatai2(x_nn(i,j),y_nn(i,j))
                            if (present(bincount)) then
                               if (odatai2(i,j).ne.missingi2) bincount(i,j) = 1
                            end if
                         end if
                       case (5) 
                         if (odatai1(i,j).ne.missingi4) then
                            odatai4(i,j) = idatai4(x_nn(i,j),y_nn(i,j))
                            if (present(bincount)) then
                               if (odatai4(i,j).ne.missingi4) bincount(i,j) = 1
                            end if
                         end if
                      case (6) 
                         if (odatai1(i,j).ne.missingi8) then
                            odatai8(i,j) = idatai8(x_nn(i,j),y_nn(i,j))
                            if (present(bincount)) then
                               if (odatai8(i,j).ne.missingi8) bincount(i,j) = 1
                            end if
                         end if
                      end select
                   end if
                end do
             end do
          end if

       end if

    end if

    ! store projection information for later re-use
    if ((.not.opt_purge).and.(.not.found)) then
       p = 1
       do while ((p.le.maxprojsave).and.(.not.found))
          if (.not.allocated(projsave(p)%ixg)) then
             found = .true.
          else
             p = p + 1
          end if
       end do
       
       if (p.gt.maxprojsave) then
          write(*,'(A,I4)') "Exceeded maximum # of storable projections: ",maxprojsave
       else
          projsave(p)%itype = itype
          projsave(p)%ixmin = ixmin
          projsave(p)%ixmax = ixmax
          projsave(p)%iymin = iymin
          projsave(p)%iymax = iymax
          projsave(p)%idx = idx
          projsave(p)%idy = idy
          projsave(p)%ioptargs = ioptargs
          projsave(p)%otype = otype
          projsave(p)%oxmin = oxmin
          projsave(p)%oxmax = oxmax
          projsave(p)%oymin = oymin
          projsave(p)%oymax = oymax
          projsave(p)%odx = odx
          projsave(p)%ody = ody
          projsave(p)%ooptargs = ooptargs
          
          projsave(p)%ibounds = ibounds
          projsave(p)%obounds = obounds
          projsave(p)%inx = inx
          projsave(p)%iny = iny
          projsave(p)%onx = onx
          projsave(p)%ony = ony
          allocate(projsave(p)%ixg(inx,iny))
          allocate(projsave(p)%iyg(inx,iny))
          projsave(p)%ixg = ixg
          projsave(p)%iyg = iyg
          allocate(projsave(p)%oxg(onx,ony))
          allocate(projsave(p)%oyg(onx,ony))
          projsave(p)%oxg = oxg
          projsave(p)%oyg = oyg
          allocate(projsave(p)%ixo(inx,iny))
          allocate(projsave(p)%iyo(inx,iny))
          projsave(p)%ixo = ixo
          projsave(p)%iyo = iyo
          allocate(projsave(p)%oxi(onx,ony))
          allocate(projsave(p)%oyi(onx,ony))
          projsave(p)%oxi = oxi
          projsave(p)%oyi = oyi

          if (ibounds) then
             allocate(projsave(p)%x_nn(onx,ony))
             allocate(projsave(p)%y_nn(onx,ony))
             projsave(p)%x_nn = x_nn
             projsave(p)%y_nn = y_nn
             allocate(projsave(p)%x_bil(onx,ony,4))
             allocate(projsave(p)%y_bil(onx,ony,4))
             allocate(projsave(p)%w_bil(onx,ony,4))
             projsave(p)%x_bil = x_bil
             projsave(p)%y_bil = y_bil
             projsave(p)%w_bil = w_bil
          end if

          if (obounds) then
             allocate(projsave(p)%rev(size(rev)))
             projsave(p)%rev = rev
             projsave(p)%omin = omin
             projsave(p)%omax = omax
          end if

          projsave(p)%scale = scale
          
       end if
    end if

    ! set optional output arrays
    if (present(ixgrid)) then
       if (allocated(ixgrid)) deallocate(ixgrid)
       allocate(ixgrid(inx,iny))
       ixgrid = ixg
    end if
    if (present(iygrid)) then
       if (allocated(iygrid)) deallocate(iygrid)
       allocate(iygrid(inx,iny))
       iygrid = iyg
    end if
    if (present(oxgrid)) then
       if (allocated(oxgrid)) deallocate(oxgrid)
       allocate(oxgrid(onx,ony))
       oxgrid = oxg
    end if
    if (present(oygrid)) then
       if (allocated(oygrid)) deallocate(oygrid)
       allocate(oygrid(onx,ony))
       oygrid = oyg
    end if
    if (present(ixoval)) then
       if (allocated(ixoval)) deallocate(ixoval)
       allocate(ixoval(inx,iny))
       ixoval = ixo
    end if
    if (present(iyoval)) then
       if (allocated(iyoval)) deallocate(iyoval)
       allocate(iyoval(inx,iny))
       iyoval = iyo
    end if
    if (present(oxival)) then
       if (allocated(oxival)) deallocate(oxival)
       allocate(oxival(onx,ony))
       oxival = oxi
    end if
    if (present(oyival)) then
       if (allocated(oyival)) deallocate(oyival)
       allocate(oyival(onx,ony))
       oyival = oyi
    end if

    deallocate(ixg,iyg,oxg,oyg)
    deallocate(ixo,iyo,oxi,oyi)

    if (ibounds) then
       deallocate(x_nn,y_nn)
       deallocate(x_bil,y_bil,w_bil)
    end if
    if (obounds) then
       deallocate(rev)
    end if

  end subroutine reproject

end module reproject_module
