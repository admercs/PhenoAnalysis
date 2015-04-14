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

module topo_module

  ! 2009/05/02 Reto Stockli (Blue Marble Research)

  ! This Module extracts sub-regions of the 1km topography grid based on CGIAR SRTM and GTOPO30.
  ! 2009/05/02 : Added support for rotated pole grid
  ! 2009/12/29 : speed up

  ! DEPENDENCIES:
  ! driver_module.F90, bilinear_module.F90, working NetCDF library
  
  character(len=100) :: topotype  ! topography data name
  character(len=100) :: topodir   ! topography data directory
  character(len=256) :: topofile  ! topography file name
  character(len=10), parameter :: toponames = 'Z'  ! NetCDF short name for Elevation
  character(len=10), parameter :: topounits = 'm'  ! NetCDF units for Elevation
  character(len=100), parameter :: topolongnames = 'Elevation above sea level' ! NetCDF long name

contains

  subroutine init_topo

    use file_module
    use exit_module

    implicit none

    topofile = trim(topodir)//'topo.'//trim(topotype)//'.nc'

    ! search file based on above pattern
    call handle_error(.not.find_file(topofile),"Topography file not found: "//topofile)

  end subroutine init_topo

  subroutine read_topo(type,xmin,xmax,ymin,ymax,dx,dy,optargs,missing,topo,topohist,topolevels)

    use file_module
    use exit_module
    use proj_module
    use reproject_module

    use netcdf
    use typeSizes

    implicit none

    ! arguments
    character(len=*), intent(in) :: type
    real(kind=4), intent(in) :: xmin
    real(kind=4), intent(in) :: xmax
    real(kind=4), intent(in) :: ymin
    real(kind=4), intent(in) :: ymax
    real(kind=4), intent(in) :: dx
    real(kind=4), intent(in) :: dy
    character(len=*), intent(in) :: optargs
    real(kind=4), intent(in) :: missing
    real(kind=4), intent(out), optional :: topo(:,:)
    real(kind=4), intent(out), optional :: topohist(:,:,:)
    real(kind=4), intent(in), optional :: topolevels(:)

    ! local
    integer :: k, h, nh

    real(kind=4) :: lonmin, lonmax, latmin, latmax
    real(kind=4) :: lonmin_topo, lonmax_topo, latmin_topo, latmax_topo

    integer :: xmin_topo, xmax_topo, ymin_topo, ymax_topo
    integer :: nx_topo, ny_topo

    integer :: nlon_topo ! W-E direction global points
    integer :: nlat_topo ! N-S direction global points
    real(kind=4), allocatable :: lons_topo(:), lats_topo(:)
    real(kind=4) :: dlon_topo, dlat_topo

    integer(kind=2), allocatable :: indata(:,:)  ! Topography stored in NetCDF
    integer(kind=2), allocatable :: inlevel(:,:) ! Topography stored in NetCDF (level by level for histogram)
    integer(kind=2), allocatable :: outdata(:,:) ! Topography transformed to output projection
    integer(kind=2) :: missing_topo ! NetCDF topography fill value
    integer(kind=2) :: topomin,topomax
    real(kind=4) :: dh

    integer(kind=4), allocatable :: outbins(:,:) ! Topography level count transformed to output projection

    ! lon/lat geolocation precision 
    ! (used to avoid out-of-bounds indices due to rounding errors)
    real(kind=4), parameter :: pr = 1.e-6

    ! NetCDF
    integer :: ncid_topo, dimid, varid
    integer :: status
    character(len=1) :: string

    ! check how far the requested area extends in the topo domain
    call getxybounds(type,xmin,xmax,ymin,ymax,dx,dy,optargs,"latitude_longitude",lonmin,lonmax,latmin,latmax,"",missing)

!    print*,xmin,xmax,ymin,ymax
!    print*,lonmin,lonmax,latmin,latmax

    ! open NetCDF file
    status = nf90_open(topofile,nf90_nowrite,ncid_topo )
    if (status /= nf90_noerr) call handle_err(status,routine='topo_module: error opening '//trim(topofile))

    ! get lon/lat number of elements
    status = nf90_inq_dimid ( ncid_topo, 'lon', dimid )
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inquire_dimension ( ncid_topo, dimid, string, nlon_topo )
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_dimid ( ncid_topo, 'lat', dimid )
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inquire_dimension ( ncid_topo, dimid, string, nlat_topo )
    if (status /= nf90_noerr) call handle_err(status)

    ! get geographic extents
    allocate(lons_topo(nlon_topo))
    allocate(lats_topo(nlat_topo))
    status = nf90_inq_varid(ncid_topo, "lon", varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_get_var(ncid_topo, varid, lons_topo )
    if(status /= nf90_NoErr) call handle_err(status)    
    status = nf90_inq_varid(ncid_topo, "lat", varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_get_var(ncid_topo, varid, lats_topo )
    if(status /= nf90_NoErr) call handle_err(status)   

    dlon_topo = lons_topo(nlon_topo/2+1) - lons_topo(nlon_topo/2)
    dlat_topo = lats_topo(nlat_topo/2+1) - lats_topo(nlat_topo/2)    
    dlon_topo = 1./real(nint(1./dlon_topo),kind=4)
    dlat_topo = 1./real(nint(1./dlon_topo),kind=4)

    lonmin_topo = minval(lons_topo) - 0.5*dlon_topo
    lonmax_topo = maxval(lons_topo) + 0.5*dlon_topo
    latmin_topo = minval(lats_topo) - 0.5*dlat_topo
    latmax_topo = maxval(lats_topo) + 0.5*dlat_topo

    xmin_topo = -1
    xmax_topo = -1
    ymin_topo = -1
    ymax_topo = -1

    ! the multiplication and nint() is a fix for numerical imprecision in comparisons
    ! south -> north
    ! west -> east

    do k=nlon_topo,1,-1
       if ((nint((lons_topo(k)-0.5*dlon_topo)/pr).le.nint(lonmin/pr)).and.(xmin_topo.eq.-1)) xmin_topo=k
    enddo
    do k=nlat_topo,1,-1
       if ((nint((lats_topo(k)-0.5*dlat_topo)/pr).le.nint(latmin/pr)).and.(ymin_topo.eq.-1)) ymin_topo=k
    enddo
    do k=1,nlon_topo
       if ((nint((lons_topo(k)+0.5*dlon_topo)/pr).ge.nint(lonmax/pr)).and.(xmax_topo.eq.-1)) xmax_topo=k
    enddo
    do k=1,nlat_topo
       if ((nint((lats_topo(k)+0.5*dlat_topo)/pr).ge.nint(latmax/pr)).and.(ymax_topo.eq.-1)) ymax_topo=k
    enddo

    if ((xmin_topo.eq.-1).or.(ymin_topo.eq.-1).or.(xmax_topo.eq.-1).or.(ymax_topo.eq.-1)) then

       write(*,'(A)') 'TOPO file does not cover the requested area: '
       write(*,'(A,4F10.5)') 'We need: ',real(nint(lonmin/pr),kind=4)*pr, real(nint(lonmax/pr))*pr, &
            real(nint(latmin/pr))*pr, real(nint(latmax/pr))*pr
       write(*,'(A,4F10.5)') 'but get: ',real(nint(lonmin_topo/pr))*pr, real(nint(lonmax_topo/pr))*pr, &
            real(nint(latmin_topo/pr))*pr, real(nint(latmax_topo/pr))*pr
       write(*,'(A,4I8)') 'TOPO File array bounds: ',xmin_topo, xmax_topo, ymin_topo, ymax_topo

       if (present(topo)) topo = missing
       if (present(topohist)) topohist = missing

    else

!       write(*,'(A,4I8)') 'TOPO File array bounds: ',xmin_topo, xmax_topo, ymin_topo, ymax_topo

       nx_topo = xmax_topo - xmin_topo + 1
       ny_topo = ymax_topo - ymin_topo + 1

       ! allocate input array
       allocate(indata(nx_topo,ny_topo))

       ! read TOPO from NetCDF 
       status = nf90_inq_varid(ncid_topo, toponames, varid)
       if(status /= nf90_NoErr) call handle_err(status)
       status = nf90_get_att(ncid_topo, varid, "_FillValue", missing_topo)       
       if(status /= nf90_NoErr) call handle_err(status)
       status = nf90_get_var(ncid_topo, varid, indata, start = (/xmin_topo,ymin_topo/), count = (/nx_topo,ny_topo/) )
       if(status /= nf90_NoErr) call handle_err(status)   

       ! close files and clean up
       status = nf90_close (ncid_topo)
       if (status /= nf90_noerr) call handle_err (status)

       if (present(topo)) then
!          print*,(lons_topo(xmax_topo)+0.5*dlon_topo-lons_topo(xmin_topo)+0.5*dlon_topo)/dlon_topo,&
!               (lats_topo(ymax_topo)+0.5*dlat_topo-lats_topo(ymin_topo)+0.5*dlat_topo)/dlat_topo
!          print*,nx_topo,ny_topo
!          print*,lons_topo(xmin_topo)-0.5*dlon_topo,lons_topo(xmax_topo)+0.5*dlon_topo,&
!               lats_topo(ymin_topo)-0.5*dlat_topo,lats_topo(ymax_topo)+0.5*dlat_topo,dlon_topo,dlat_topo
!          print*,xmin,xmax,ymin,ymax
          call reproject("latitude_longitude",lons_topo(xmin_topo)-0.5*dlon_topo,lons_topo(xmax_topo)+0.5*dlon_topo,&
               lats_topo(ymin_topo)-0.5*dlat_topo,lats_topo(ymax_topo)+0.5*dlat_topo,dlon_topo,dlat_topo,"",&
               type,xmin,xmax,ymin,ymax,dx,dy,optargs,missing,purge=.true.,idatai2=indata,odatai2=outdata,&
               missingi2=missing_topo)

          where(outdata.ne.missing_topo)
             topo = real(outdata,kind=4)
          elsewhere
             topo = missing
          end where

          deallocate(outdata)
       end if

       if (present(topohist)) then
          if (present(topolevels)) then
             nh = size(topolevels)
             if (nh.gt.1) then
                allocate(inlevel(nx_topo,ny_topo))
                dh = (topolevels(nh)-topolevels(1))/real(nh-1)
                do h=1,nh
                   topomin = nint(topolevels(h)-0.5*dh,kind=2)
                   topomax = nint(topolevels(h)+0.5*dh,kind=2)
                   if (h.eq.1) topomin = -9000
                   if (h.eq.nh) topomax = 9000
                   where ((indata.ge.topomin).and.(indata.lt.topomax))
                      inlevel = indata
                   elsewhere
                      inlevel = missing_topo
                   end where

                   call reproject("latitude_longitude",lons_topo(xmin_topo)-0.5*dlon_topo,lons_topo(xmax_topo)+0.5*dlon_topo,&
                        lats_topo(ymin_topo)-0.5*dlat_topo,lats_topo(ymax_topo)+0.5*dlat_topo,dlon_topo,dlat_topo,"",&
                        type,xmin,xmax,ymin,ymax,dx,dy,optargs,missing,purge=.true.,&
                        idatai2=inlevel,odatai2=outdata,missingi2=missing_topo,bincount=outbins)

                   topohist(:,:,h) = real(outbins,kind=4)
                end do
                deallocate(inlevel)
                deallocate(outbins)
                deallocate(outdata)

                do h=1,nh
                   where(sum(topohist,dim=3).lt.pr)
                      topohist(:,:,h) = missing
                   end where
                end do
             else
                topohist(:,:,:) = 1.0
             end if
          else
             topohist = missing
          end if
       end if

       deallocate(indata)

    endif

    deallocate(lons_topo,lats_topo)

  end subroutine read_topo

end module topo_module
