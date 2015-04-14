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


module pft_module
  
  ! 2009/05/02 Reto Stockli (Blue Marble Research)

  ! This module reads % of PFT's for each 1km pixel for a chosen area from the fractional vegetation cover map
  ! It currently supports the full NCAR CLM-type 17 PFT classification.
  ! 2009/05/02 : Added support for rotated pole grid
  ! 2009/12/29 : speed up

  ! Dependencies
  ! bilinear_module.F90, working NetCDF Library with F90 interface
  ! driver_module: handle_err()

  implicit none

  character(len=100) :: pfttype    ! name of pft dataset
  character(len=100) :: pftdir     ! directory of pft dataset
  character(len=256) :: pftfile    ! full name of pft file
  logical :: manyfiles_pft ! per-pft NetCDF files (global 1km resolution) or single NetCDF file (by region)

  ! 35 PFT's
  integer, parameter :: numpft = 35
  character(len=10), parameter :: pftunits(numpft) = (/ &
       '%','%','%','%','%','%','%','%','%','%','%','%','%','%','%','%','%', &
       '%','%','%','%','%','%','%','%','%','%','%','%','%','%','%','%','%','%'/)

  character(len=10), parameter :: pftnames(numpft) = (/ &
       'bar_all',&
       'enf_tem',&
       'enf_bor',&
       'dnf_bor',&
       'ebf_tro',&
       'ebf_tem',&
       'dbf_tro',&
       'dbf_tem',&
       'dbf_bor',&
       'ebs_all',&
       'dbs_tem',&
       'dbs_bor',&
       'c3g_arc',&
       'c3g_nar',&
       'c4g_all',&
       'cro_brl',&
       'cro_cas',&
       'cro_cot',&
       'cro_grn',&
       'cro_mze',&
       'cro_mil',&
       'cro_oil',&
       'cro_oth',&
       'cro_pot',&
       'cro_pul',&
       'cro_rap',&
       'cro_ric',&
       'cro_rye',&
       'cro_sor',&
       'cro_soy',&
       'cro_sgb',&
       'cro_sgc',&
       'cro_sun',&
       'cro_wht',&
       'wat_all'/)

  character(len=50), parameter :: pftlongnames(numpft) = (/ & 
       'Bare Soil (-)                         ',&
       'Evergreen Needleleaf Trees (Temperate)',&
       'Evergreen Needleleaf Trees (Boreal)   ',&
       'Deciduous Needleleaf Trees (Boreal)   ',&
       'Evergreen Broadleaf Trees  (Tropical) ',&
       'Evergreen Broadleaf Trees  (Temperate)',&
       'Deciduous Broadleaf Trees  (Tropical) ',&
       'Deciduous Broadleaf Trees  (Temperate)',&
       'Deciduous Broadleaf Trees  (Boreal)   ',&
       'Evergreen Broadleaf Shrubs (All)      ',&
       'Deciduous Broadleaf Shrubs (Temperate)',&
       'Deciduous Broadleaf Shrubs (Boreal)   ',&
       'C3 Grass (Arctic)                     ',&
       'C3 Grass (Non-Arctic)                 ',&
       'C4 Grass (All)                        ',&
       'Crops (Barley)                        ',&
       'Crops (Cassava)                       ',&
       'Crops (Cotton)                        ',&
       'Crops (Groundnuts)                    ',&
       'Crops (Maize)                         ',&
       'Crops (Millet)                        ',&
       'Crops (Oilpalm)                       ',&
       'Crops (Other)                         ',&
       'Crops (Potato)                        ',&
       'Crops (Pulses)                        ',&
       'Crops (Rape)                          ',&
       'Crops (Rice)                          ',&
       'Crops (Rye)                           ',&
       'Crops (Sorghum)                       ',&
       'Crops (Soy)                           ',&
       'Crops (Sugarbeets)                    ',&
       'Crops (Sugarcane)                     ',&
       'Crops (Sunflower)                     ',&
       'Crops (Wheat)                         ',&
       'Water (-)                             '/)

contains

  subroutine init_pft

    ! call this before calling read_pft

    ! modules
    use file_module
    use exit_module
    use netcdf
    use typeSizes

    implicit none

    ! local
    integer :: nDimensions, nVariables, nAttributes, unlimdimid, formatnum

    ! NetCDF
    integer :: ncid_pft
    integer :: status
   
    integer :: n

    ! check whether PFT's are stored in a single file
    pftfile = trim(pftdir)//'pft.'//trim(pfttype)//'.nc'

    if (find_file(pftfile)) then
       ! PFT's contained in single file

       manyfiles_pft = .false. 

       ! check for number of PFT's within this file
       status = nf90_open(pftfile, nf90_nowrite,ncid_pft )
       if (status /= nf90_noerr) call handle_err(status,routine='pft_module: error opening '// &
            trim(pftdir)//'pft.'//trim(pfttype)//'.nc'//' NetCDF file')

       status = nf90_inquire(ncid_pft, nDimensions, nVariables, nAttributes, unlimdimid, formatnum)
       if (status /= nf90_noerr) call handle_err(status,routine='pft_module: error inquiring '// &
            trim(pftdir)//'pft.'//trim(pfttype)//'.nc'//' NetCDF file')

       ! fix for xlf (it doesn't get ndimensions read correctly from file):
       nDimensions = 2

       ! nDimension variables are used to specify dimensions (e.g. lon,lat), the rest are pft's
       call handle_error(numpft.ne.(nVariables - nDimensions),"PFT file contains more or less PFT's than defined ")

       status = nf90_close (ncid_pft)
       if (status /= nf90_noerr) call handle_err (status)

    else
       ! PFT's spread over several files: one file per PFT
       manyfiles_pft = .true.

       do n=1,numpft
          pftfile = trim(pftdir)//trim(pftnames(n))//'.'//trim(pfttype)//'.nc'
          call handle_error(.not.find_file(pftfile),"PFT file not found: "//trim(pftfile))
       end do
    endif

  end subroutine init_pft

  subroutine read_pft(type,xmin,xmax,ymin,ymax,dx,dy,optargs,missing,pft3d,pft2d,pftsel)

    use file_module
    use exit_module
    use netcdf
    use typeSizes
    use proj_module
    use reproject_module

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
    real(kind=4), intent(out), optional :: pft2d(:,:)
    real(kind=4), intent(out), optional :: pft3d(:,:,:)
    integer, intent(in), optional :: pftsel

    ! local
    integer :: k, n, n0, n1

    real(kind=4) :: lonmin, lonmax, latmin, latmax
    real(kind=4) :: lonmin_pft, lonmax_pft, latmin_pft, latmax_pft

    integer :: xmin_pft, xmax_pft, ymin_pft, ymax_pft
    integer :: nx_pft, ny_pft

    integer :: nlon_pft ! W-E direction global points
    integer :: nlat_pft ! N-S direction global points
    real(kind=4), allocatable :: lons_pft(:), lats_pft(:)
    real(kind=4) :: dlon_pft, dlat_pft

    integer(kind=1), allocatable :: indatai1(:,:)  ! PFT stored in NetCDF
    real(kind=4), allocatable :: indata(:,:)  ! PFT stored in NetCDF
    real(kind=4), allocatable :: outdata(:,:) ! PFT transformed to output projection
    integer(kind=1) :: missing_pft

    ! lon/lat geolocation precision 
    ! (used to avoid out-of-bounds indices due to rounding errors)
    real(kind=4), parameter :: pr = 1.e-6
    real(kind=4) :: eps

    ! NetCDF
    integer :: ncid_pft, dimid, varid
    integer :: status
    character(len=1) :: string

   ! determine smallest difference in 32 bit
    eps = 10.0*epsilon(0.0)
 
    if ((.not.present(pft2d)).and.(.not.present(pft3d))) then
       write(*,'(A)') "Please provide either pft2d or pft3d as arguments"
       stop
    end if
    if (present(pft2d).and.(.not.present(pftsel))) then
       write(*,'(A)') "Please provide pftsel argument when setting pft2d argument"
       stop
    end if
    if (present(pftsel)) then
       if ((pftsel.lt.1).or.(pftsel.gt.numpft)) then
          write(*,'(A)') "pftsel has to be within the range 1..numpft"
          stop
       end if
    end if

    ! check how far the requested area extends in the pft domain
    call getxybounds(type,xmin,xmax,ymin,ymax,dx,dy,optargs,"latitude_longitude",lonmin,lonmax,latmin,latmax,"",missing)

    status = nf90_open(pftfile,nf90_nowrite,ncid_pft)
    if (status /= nf90_noerr) call handle_err(status,routine='pft_module: error opening '//trim(pftfile))

    ! get lon/lat number of elements
    status = nf90_inq_dimid ( ncid_pft, 'lon', dimid )
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inquire_dimension ( ncid_pft, dimid, string, nlon_pft )
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_dimid ( ncid_pft, 'lat', dimid )
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inquire_dimension ( ncid_pft, dimid, string, nlat_pft )
    if (status /= nf90_noerr) call handle_err(status)

    ! get geographic extents
    allocate(lons_pft(nlon_pft))
    allocate(lats_pft(nlat_pft))
    status = nf90_inq_varid(ncid_pft, "lon", varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_get_var(ncid_pft, varid, lons_pft )
    if(status /= nf90_NoErr) call handle_err(status)    

    status = nf90_inq_varid(ncid_pft, "lat", varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_get_var(ncid_pft, varid, lats_pft )
    if(status /= nf90_NoErr) call handle_err(status)   

    ! close first PFT file
    status = nf90_close (ncid_pft)
    if (status /= nf90_noerr) call handle_err (status)

    dlon_pft = lons_pft(nlon_pft/2+1) - lons_pft(nlon_pft/2)
    dlat_pft = lats_pft(nlat_pft/2+1) - lats_pft(nlat_pft/2)
    dlon_pft = 1./real(nint(1./dlon_pft),kind=4)
    dlat_pft = 1./real(nint(1./dlon_pft),kind=4)

    lonmin_pft = minval(lons_pft) - 0.5*dlon_pft
    lonmax_pft = maxval(lons_pft) + 0.5*dlon_pft
    latmin_pft = minval(lats_pft) - 0.5*dlat_pft
    latmax_pft = maxval(lats_pft) + 0.5*dlat_pft

    xmin_pft = -1
    xmax_pft = -1
    ymin_pft = -1
    ymax_pft = -1

    ! the multiplication and nint() is a fix for numerical imprecision in comparisons
    ! south -> north
    ! west -> east

    do k=nlon_pft,1,-1
       if ((nint((lons_pft(k)-0.5*dlon_pft)/pr).le.nint(lonmin/pr)).and.(xmin_pft.eq.-1)) xmin_pft=k
    enddo
    do k=nlat_pft,1,-1
       if ((nint((lats_pft(k)-0.5*dlat_pft)/pr).le.nint(latmin/pr)).and.(ymin_pft.eq.-1)) ymin_pft=k
    enddo
    do k=1,nlon_pft
       if ((nint((lons_pft(k)+0.5*dlon_pft)/pr).ge.nint(lonmax/pr)).and.(xmax_pft.eq.-1)) xmax_pft=k
    enddo
    do k=1,nlat_pft
       if ((nint((lats_pft(k)+0.5*dlat_pft)/pr).ge.nint(latmax/pr)).and.(ymax_pft.eq.-1)) ymax_pft=k
    enddo

    if ((xmin_pft.eq.-1).or.(ymin_pft.eq.-1).or.(xmax_pft.eq.-1).or.(ymax_pft.eq.-1)) then

       write(*,'(A)') 'PFT file does not cover the requested area: '
       write(*,'(A,4F10.5)') 'We need: ',real(nint(lonmin/pr))*pr, real(nint(lonmax/pr))*pr, &
            real(nint(latmin/pr))*pr, real(nint(latmax/pr))*pr
       write(*,'(A,4F10.5)') 'but get: ',real(nint(lonmin_pft/pr))*pr, real(nint(lonmax_pft/pr))*pr, &
            real(nint(latmin_pft/pr))*pr, real(nint(latmax_pft/pr))*pr
       write(*,'(A,4I8)') 'PFT File array bounds: ',xmin_pft, xmax_pft, ymin_pft, ymax_pft

       if (present(pft2d)) pft2d = missing
       if (present(pft3d)) pft3d = missing

    else

!       write(*,'(A,4I8)') 'PFT File array bounds: ',xmin_pft, xmax_pft, ymin_pft, ymax_pft

       nx_pft = xmax_pft - xmin_pft + 1
       ny_pft = ymax_pft - ymin_pft + 1

       ! allocate input array
       allocate(indatai1(nx_pft,ny_pft))
       allocate(indata(nx_pft,ny_pft))

       ! get PFT distribution 
       if (.not.manyfiles_pft) then
          status = nf90_open(trim(pftfile), nf90_nowrite,ncid_pft )
          if (status /= nf90_noerr) call handle_err(status,routine='pft_module: error opening '//trim(pftfile))
       endif

       if (present(pft3d)) then
          n0 = 1
          n1 = numpft
       else
          n0 = pftsel
          n1 = pftsel
       end if

       do n=n0,n1
          if (manyfiles_pft) then
             pftfile = trim(pftdir)//trim(pftnames(n))//'.'//trim(pfttype)//'.nc'
             status = nf90_open(trim(pftfile), nf90_nowrite,ncid_pft )
             if (status /= nf90_noerr) call handle_err(status,routine='pft_module: error opening '//trim(pftfile))
          endif

          status = nf90_inq_varid(ncid_pft, pftnames(n), varid)
          if(status /= nf90_NoErr) call handle_err(status)
          status = nf90_get_att(ncid_pft, varid, "_FillValue", missing_pft)       
          if(status /= nf90_NoErr) call handle_err(status)
          status = nf90_get_var(ncid_pft, varid, indatai1, start = (/xmin_pft,ymin_pft,1,1/), count = (/nx_pft,ny_pft,1,1/) )
          if(status /= nf90_NoErr) call handle_err(status)   

          where (indatai1.ne.missing_pft)
             indata = real(indatai1,kind=4)/100.0
          elsewhere
             indata = missing
          end where

          call reproject("latitude_longitude",lons_pft(xmin_pft)-0.5*dlon_pft,lons_pft(xmax_pft)+0.5*dlon_pft,&
               lats_pft(ymin_pft)-0.5*dlat_pft,lats_pft(ymax_pft)+0.5*dlat_pft,dlon_pft,dlat_pft,"",&
               type,xmin,xmax,ymin,ymax,dx,dy,optargs,missing,purge=.true.,idataf4=indata,odataf4=outdata,&
               missingf4=missing)

          if (present(pft3d)) then
             where (abs(outdata-missing).gt.eps)
                pft3d(:,:,n) = outdata
             elsewhere
                pft3d(:,:,n) = missing
             end where
          else
             where (abs(outdata-missing).gt.eps)
                pft2d = outdata
             elsewhere
                pft2d = missing
             end where
          end if

          if (manyfiles_pft) then
             status = nf90_close (ncid_pft)
             if (status /= nf90_noerr) call handle_err (status)
          endif

       enddo

       if (.not.manyfiles_pft) then
          status = nf90_close (ncid_pft)
          if (status /= nf90_noerr) call handle_err (status)
       endif

       if (allocated(outdata)) deallocate(outdata)
       deallocate(indata)

    endif

    deallocate(lons_pft,lats_pft)

  end subroutine read_pft

end module pft_module
