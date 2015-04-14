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

module netcdf_module

  ! 2009/02/20 Reto Stockli (Blue Marble Research) 

  ! MODULE for writing NetCDF files with aggregation possibilities in all of the 5 dimensions
  ! specifically created for subgrid-scale pft and hgt dimensions in the phenological data assimilation system
  
  ! July 2008: added: if a file exists from previous simulation we will append to it instead of re-creating it
  ! 2008/12/30: added: support for multiple output file streams (e.g. model, simulation, observation)

  ! DEPENDENCIES:
  ! A working NetCDF library with F90 support

  ! NetCDF arrays 
  character(len=255), allocatable :: NCfilenames(:)      ! names of NetCDF files

  character(len=100), allocatable :: NCnames(:,:)         ! names of variables to be stored in NetCDF file
  character(len=100), allocatable :: NClongnames(:,:)    ! long names of variables to be stored in NetCDF file
  character(len=100), allocatable :: NCunits(:,:)         ! units of variables to be stored in NetCDF file
  real(kind=4), allocatable :: NCscale_factor(:,:)       ! scale_factor of variables to be stored in NetCDF file
  real(kind=4), allocatable :: NCadd_offset(:,:)         ! add_offset of variables to be stored in NetCDF file
  character(len=100), allocatable :: NCtype(:,:)          ! type of variables to be stored in NetCDF file

  character(len=100), allocatable :: NCdimnames(:,:)      ! names of dimensions to be stored in NetCDF file
  character(len=100), allocatable :: NCdimvarnames(:,:)   ! names of dimension variables to be stored in NetCDF file
  character(len=100), allocatable :: NCdimlongnames(:,:) ! long names of dimensions to be stored in NetCDF file
  character(len=100), allocatable :: NCdimunits(:,:)      ! units of dimensions to be stored in NetCDF file

  integer, allocatable :: NCdimid(:)                     ! NetCDF dimension ID's
  integer, allocatable :: NCdimvarid(:)                  ! NetCDF dimension variable ID's
  integer, allocatable :: NCvarid(:,:)                   ! NetCDF variable ID's
  integer, allocatable :: NCdim(:,:,:)                   ! chosen dimension per NCvar
  logical, allocatable :: NCdimave(:,:)                  ! average output across a given dimension or not
  logical, allocatable :: NCdimcompress(:,:)             ! compress values across output dimension by mean and sd

  logical, allocatable :: NChasfile(:)         ! whether to write/open/read this file or not
  integer, allocatable :: Ncid_output(:)       ! NetCDF file ID
  logical :: NCparallel                        ! whether Parallel I/O is used (only works with NetCDF-4)

  integer NCnvar     ! number of variables to store in NetCDF files
  integer NCnauxvar  ! number of auxiliary variables (not dimension and regular model variables)
  integer NCnfiles   ! Number of Output Files

  integer NCcomm     ! MPI communicator for selected I/O processes

  real(kind=4) :: NCmissing           ! missing value

  integer, parameter :: NCmaxdim = 6  ! maximum NetCDF dimensions

contains

subroutine netcdf_allocate

  implicit none

  ! allocate NetCDF output arrays
  allocate(NCfilenames(NCnfiles))

  allocate(NCnames(NCnvar,NCnfiles))
  allocate(NClongnames(NCnvar,NCnfiles))
  allocate(NCunits(NCnvar,NCnfiles))
  allocate(NCscale_factor(NCnvar,NCnfiles))
  allocate(NCadd_offset(NCnvar,NCnfiles))
  allocate(NCtype(NCnvar,NCnfiles))

  allocate(NCdimnames(NCmaxdim,NCnfiles))
  allocate(NCdimvarnames(NCmaxdim,NCnfiles))
  allocate(NCdimlongnames(NCmaxdim,NCnfiles))
  allocate(NCdimunits(NCmaxdim,NCnfiles))
  allocate(NCdimid(NCmaxdim))
  allocate(NCdimvarid(NCmaxdim))
  allocate(NCvarid(NCnvar,NCnfiles))
  allocate(NCdim(NCnvar,NCmaxdim,NCnfiles))
  allocate(NCdimave(NCmaxdim,NCnfiles))
  allocate(NCdimcompress(NCmaxdim,NCnfiles))

  allocate(NChasfile(NCnfiles))
  allocate(ncid_output(NCnfiles))

end subroutine netcdf_allocate

subroutine netcdf_deallocate

  implicit none

  ! deallocate NetCDF output arrays
  deallocate(NCfilenames)

  deallocate(NCnames)
  deallocate(NClongnames)
  deallocate(NCunits)
  deallocate(NCscale_factor)
  deallocate(NCadd_offset)
  deallocate(NCtype)

  deallocate(NCdimnames)
  deallocate(NCdimvarnames)
  deallocate(NCdimlongnames)
  deallocate(NCdimunits)
  deallocate(NCdimid)
  deallocate(NCdimvarid)
  deallocate(NCvarid)
  deallocate(NCdim)
  deallocate(NCdimave)
  deallocate(NCdimcompress)

  deallocate(NChasfile)
  deallocate(ncid_output)
  
end subroutine netcdf_deallocate

subroutine netcdf_init(file, ens, xcoord, ycoord, pft, hgt, year, name, version, nodata, ntime, type, optargs)

  use netcdf
  use typeSizes

  use proj_module

#ifndef HIDE_MPI
  use mpi
#endif

  implicit none

  ! arguments
  integer, intent(in)      :: file
  integer, intent(in)      :: ens(:)
  real(kind=4), intent(in) :: xcoord(:)
  real(kind=4), intent(in) :: ycoord(:)
  integer, intent(in)      :: pft(:)
  real(kind=4), intent(in) :: hgt(:)
  integer, intent(in)      :: year
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: version
  real(kind=4), intent(in) :: nodata
  integer, intent(in)      :: ntime
  character(len=*), intent(in) :: type
  character(len=*), intent(in) :: optargs

  ! local
  character(len=4) :: syear
  integer :: a, v, d, k, ndims
  logical :: append
  character(len=8) :: datestr
  character(len=10) :: timestr

  ! NC dimensions
  integer, allocatable :: dimvec(:)

  ! dimension ids (NetCDF)
  integer :: dimsize(NCmaxdim)

  ! File id's (NetCDF)
  integer :: status
  character(len=1) :: string
  integer :: ival
  integer :: mode_flag

  ! projection arguments
  character(len=100), allocatable :: arguments(:)
  real(kind=4), allocatable :: fvalues(:)
  integer, allocatable :: ivalues(:)
  character(len=100), allocatable :: cvalues(:)
  integer :: narguments
  
  ! set missing value
  NCmissing = nodata

  ! determine dimension sizes
  dimsize(1) = size(ens)
  dimsize(2) = size(xcoord)
  dimsize(3) = size(ycoord)
  dimsize(4) = size(pft)
  dimsize(5) = size(hgt)
!  dimsize(6) = NF90_UNLIMITED    !! does not work with independent access of parallel NetCDF 4 I/O
  dimsize(6) = ntime  
  
  ! check if dimension d gets averaged in output
  do d = 1,NCmaxdim
     if (NCdimave(d,file)) then
!        dimsize(d) = 1       ! becomes singular: n=1 (ensemble mean)
        NCdim(NCnauxvar+1:NCnvar,d,file) = 0  ! and is therefore removed from dimension vector
     endif
  end do

  ! check if dimension d gets compressed in output
  do d = 1,NCmaxdim
     if (NCdimcompress(d,file)) then
!        dimsize(d) = 1       ! becomes singular: n=1 (ensemble mean)
        NCdim(NCnauxvar+1:NCnvar,d,file) = 0  ! and is therefore removed from the dimension vector
     endif
  end do

  ! create file and enter define mode
  write(unit=syear,fmt='(i4)') year

  append = .false.

  ! check if there already is a file with all dimensions and variables defined

  mode_flag = nf90_write 
  if (NCparallel) then
#ifndef HIDE_MPI
#ifndef HIDE_NETCDF4_MPI
     mode_flag = ior(mode_flag, nf90_netcdf4)
     mode_flag = ior(mode_flag, nf90_classic_model) 
     mode_flag = ior(mode_flag, nf90_mpiio) 
     status = nf90_open(trim(NCfilenames(file)), mode_flag, &
          ncid_output(file), comm = NCcomm, info = MPI_INFO_NULL)

#endif
#endif
  else 
     status = nf90_open(trim(NCfilenames(file)), mode_flag, ncid_output(file) )
  endif

  if (status == nf90_noerr) then

     append = .true. ! assume file is ok

     ! inquire dimension NetCDF id's
     do d=1,NCmaxdim
        status = nf90_inq_dimid ( ncid_output(file),  NCdimnames(d,file), NCdimid(d) )
        if (status .ne. nf90_noerr) write(*,'(A)')  nf90_strerror(status)
        status = nf90_inquire_dimension ( ncid_output(file), NCdimid(d), string, ival)
        if (status .ne. nf90_noerr) write(*,'(A)')  nf90_strerror(status)

        if (d.ne.NCmaxdim) then
           ! check compatibility dimensions (not applied for time dimension)
           if (ival.ne.dimsize(d)) append = .false. ! sorry, file is not ok
        endif

     enddo

     ! inquire dimension variable NetCDF id's 
     do d=1,NCmaxdim
        status = nf90_inq_varid ( ncid_output(file),  trim(NCdimvarnames(d,file)), NCdimvarid(d) )
        if (status .ne. nf90_noerr) append = .false.
     enddo

     ! inquire variable NetCDF id's 
     do v=1,NCnvar
        if (trim(NCnames(v,file)).ne.'') then
           status = nf90_inq_varid ( ncid_output(file),  trim(NCnames(v,file)), ncvarid(v,file) )
           if (status .ne. nf90_noerr) append = .false.
        endif
     enddo

     if (.not.append) then
        ! close the file since it is not complete
        status = nf90_close(ncid_output(file))
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     endif

  endif

  if (.not.append) then
     ! open new file and create dimensions and variables first

     mode_flag = nf90_clobber
     if (NCparallel) then
#ifndef HIDE_MPI
#ifndef HIDE_NETCDF4_MPI
        mode_flag = ior(mode_flag, nf90_netcdf4)
        mode_flag = ior(mode_flag, nf90_classic_model) 
        mode_flag = ior(mode_flag, nf90_mpiio) 
        status = nf90_create(trim(NCfilenames(file)), mode_flag, &
             ncid_output(file), comm = NCcomm, info = MPI_INFO_NULL)

#endif
#endif
     else 
        status = nf90_create(trim(NCfilenames(file)), mode_flag, ncid_output(file) )
     endif
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)//&
          '. Netcdf_module: error opening file: '//trim(NCfilenames(file))

     ! define dimensions
     do d=1,NCmaxdim
        status = nf90_def_dim(ncid_output(file), NCdimnames(d,file),dimsize(d), NCdimid(d))
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     enddo

     d=1 ! ensemble members
     status = nf90_def_var(ncid_output(file), NCdimvarnames(d,file), nf90_short, NCdimid(d), NCdimvarid(d))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     status = nf90_put_att(ncid_output(file), NCdimvarid(d), 'long_name',trim(NCdimlongnames(d,file)))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     status = nf90_put_att(ncid_output(file), NCdimvarid(d), 'units',trim(NCdimunits(d,file)))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     d=2 ! lon or x coordinate
     status = nf90_def_var(ncid_output(file), NCdimvarnames(d,file), nf90_real, NCdimid(d), NCdimvarid(d))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     status = nf90_put_att(ncid_output(file), NCdimvarid(d), 'long_name',trim(NCdimlongnames(d,file)))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     status = nf90_put_att(ncid_output(file), NCdimvarid(d), 'units',trim(NCdimunits(d,file)))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     d=3 ! lat or y coordinate
     status = nf90_def_var(ncid_output(file), NCdimvarnames(d,file), nf90_real, NCdimid(d), NCdimvarid(d))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     status = nf90_put_att(ncid_output(file), NCdimvarid(d), 'long_name',trim(NCdimlongnames(d,file)))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     status = nf90_put_att(ncid_output(file), NCdimvarid(d), 'units',trim(NCdimunits(d,file)))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     d=4 ! pfts
     status = nf90_def_var(ncid_output(file), NCdimvarnames(d,file), nf90_short, &
          (/NCdimid(d)/), NCdimvarid(d))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     status = nf90_put_att(ncid_output(file), NCdimvarid(d), 'long_name',trim(NCdimlongnames(d,file)))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     status = nf90_put_att(ncid_output(file), NCdimvarid(d), 'units',trim(NCdimunits(d,file)))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     d=5 ! elevation
     status = nf90_def_var(ncid_output(file), NCdimvarnames(d,file), nf90_float, &
          (/NCdimid(d)/), NCdimvarid(d))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     status = nf90_put_att(ncid_output(file), NCdimvarid(d), 'long_name',trim(NCdimlongnames(d,file)))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     status = nf90_put_att(ncid_output(file), NCdimvarid(d), 'units',trim(NCdimunits(d,file)))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     d=6 ! time
     status = nf90_def_var(ncid_output(file), NCdimvarnames(d,file), nf90_float, NCdimid(d), NCdimvarid(d))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     status = nf90_put_att(ncid_output(file), NCdimvarid(d), 'long_name',trim(NCdimlongnames(d,file)))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     status = nf90_put_att(ncid_output(file), NCdimvarid(d), 'units',trim(NCdimunits(d,file))//' '//syear//'-01-01 00:00:00')
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     status = nf90_put_att(ncid_output(file), NCdimvarid(d), 'calendar','gregorian')
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)  

     ! get date and time
     call date_and_time(datestr,timestr)

     ! define Variables
     do v=1,NCnvar

        if (trim(NCnames(v,file)).ne.'') then

           ! check if this variable is present in file (any dimensions?)
           ndims = sum(NCdim(v,:,file))

           allocate(dimvec(ndims))

           k=1
           do d = 1, NCmaxdim
              if (NCdim(v,d,file).eq.1) then
                 dimvec(k) = NCdimid(d)
                 k=k+1
              endif
           enddo

           select case (NCtype(v,file))
           case ("byte")
              status = nf90_def_var(ncid_output(file), trim(NCnames(v,file)), NF90_BYTE, dimvec, NCvarid(v,file))
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
              status = nf90_put_att(ncid_output(file), NCvarid(v,file), '_FillValue',NF90_FILL_BYTE)
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
           case ("short")
              status = nf90_def_var(ncid_output(file), trim(NCnames(v,file)), NF90_SHORT, dimvec, NCvarid(v,file))
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
              status = nf90_put_att(ncid_output(file), NCvarid(v,file), '_FillValue',NF90_FILL_SHORT)
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
           case ("int")
              status = nf90_def_var(ncid_output(file), trim(NCnames(v,file)), NF90_INT, dimvec, NCvarid(v,file))
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
              status = nf90_put_att(ncid_output(file), NCvarid(v,file), '_FillValue',NF90_FILL_INT)
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
           case ("float")
              status = nf90_def_var(ncid_output(file), trim(NCnames(v,file)), NF90_FLOAT, dimvec, NCvarid(v,file))
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
              status = nf90_put_att(ncid_output(file), NCvarid(v,file), '_FillValue',NF90_FILL_FLOAT)
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
           case ("double")
              status = nf90_def_var(ncid_output(file), trim(NCnames(v,file)), NF90_DOUBLE, dimvec, NCvarid(v,file))
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
              status = nf90_put_att(ncid_output(file), NCvarid(v,file), '_FillValue',NF90_FILL_DOUBLE)
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
           case ("char")
              status = nf90_def_var(ncid_output(file), trim(NCnames(v,file)), NF90_CHAR, dimvec, NCvarid(v,file))
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
           end select

           if ((NCscale_factor(v,file).ne.1.).or.(NCadd_offset(v,file).ne.0.)) then
              status = nf90_put_att(ncid_output(file), NCvarid(v,file), 'scale_factor',NCscale_factor(v,file))
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
              status = nf90_put_att(ncid_output(file), NCvarid(v,file), 'add_offset',NCadd_offset(v,file))
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
           endif

           if (trim(NCnames(v,file)).eq."grid_mapping") then
              status = nf90_put_att(ncid_output(file), NCvarid(v,file), 'grid_mapping_name',trim(type))
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
              if (trim(type).eq."latitude_longitude") then
                 status = nf90_put_att(ncid_output(file), NCvarid(v,file), 'longitude_of_prime_meridian',0.d0)
                 if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
                 status = nf90_put_att(ncid_output(file), NCvarid(v,file), 'semi_major_axis',6378137.d0)
                 if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
                 status = nf90_put_att(ncid_output(file), NCvarid(v,file), 'inverse_flattening',298.257223563d0)
                 if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
              else
                 narguments = getoptargs(optargs,arguments,fvalues,ivalues=ivalues,cvalues=cvalues)

                 do a=1,narguments
                    if (trim(cvalues(a)).ne."") then
                       status = nf90_put_att(ncid_output(file), NCvarid(v,file), trim(arguments(a)),cvalues(a))
                       if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
                    else
                       if (ivalues(a).ne.0) then
                          status = nf90_put_att(ncid_output(file), NCvarid(v,file), trim(arguments(a)),ivalues(a))
                          if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
                       else
                          status = nf90_put_att(ncid_output(file), NCvarid(v,file), trim(arguments(a)),fvalues(a))
                          if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
                       end if
                    end if
                 end do

                 if (allocated(arguments)) deallocate(arguments)
                 if (allocated(fvalues)) deallocate(fvalues)
                 if (allocated(ivalues)) deallocate(ivalues)
                 if (allocated(cvalues)) deallocate(cvalues)
              end if
           else
              status = nf90_put_att(ncid_output(file), NCvarid(v,file), 'long_name',trim(NClongnames(v,file)))
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
              status = nf90_put_att(ncid_output(file), NCvarid(v,file), 'units',trim(NCunits(v,file)))
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
              status = nf90_put_att(ncid_output(file), NCvarid(v,file), 'version',trim(version))
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
              status = nf90_put_att(ncid_output(file), NCvarid(v,file), &
                   'prod_date',datestr(1:4)//'-'//datestr(5:6)// &
                   '-'//datestr(7:8)//' '//timestr(1:2)//':'//timestr(3:4)//':'//timestr(5:6))
              if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

              ! add reference to coordinate system for geographically mapped variables
              if ((NCdim(v,2,file).eq.1).and.(NCdim(v,3,file).eq.1)) then
                 status = nf90_put_att(ncid_output(file), NCvarid(v,file),'grid_mapping','grid_mapping')
                 if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
              end if
            end if

           deallocate(dimvec)

        endif

     enddo

     ! put global attributes
     status = nf90_put_att(ncid_output(file), nf90_global, 'Title',trim(name))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     status = nf90_put_att(ncid_output(file), nf90_global, 'Summary', &
          'Global Vegetation State Estimation with MODIS observations using the Ensemble Kalman Filter')
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     status = nf90_put_att(ncid_output(file), nf90_global, 'Institution','Blue Marble Research')
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     status = nf90_put_att(ncid_output(file), nf90_global, 'References', &
          'R. Stockli, T. Rutishauser, I. Baker, M. Liniger, and A. S. Denning.'// &
          'A global reanalysis of vegetation phenology. J. Geophys. Res. -- '// &
          'Biogeosciences, 116(G03020), 2011. doi: 10.1029/2010JG001545.')
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     status = nf90_put_att(ncid_output(file), nf90_global, 'Conventions','CF-1.6')
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     status = nf90_put_att(ncid_output(file), nf90_global, 'Digital_Object_Identifier','not assigned')
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     status = nf90_put_att(ncid_output(file), nf90_global, 'Creator_Name','Reto Stockli')
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     status = nf90_put_att(ncid_output(file), nf90_global, 'Creator_Email','reto.stockli@gmail.com')
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     status = nf90_put_att(ncid_output(file), nf90_global, 'Creator_Url','http://www.bluemarble.ch')
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     ! leave define mode
     status = NF90_ENDDEF(ncid_output(file))
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

     ! write dimension variables (except time)
     d=1 
     status = nf90_put_var( ncid_output(file), NCdimvarid(d), int(ens(1:dimsize(d)),2), (/1/), (/dimsize(d)/) )
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     d=2
     status = nf90_put_var( ncid_output(file), NCdimvarid(d), xcoord, (/1/), (/dimsize(d)/) )
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     d=3
     status = nf90_put_var( ncid_output(file), NCdimvarid(d), ycoord, (/1/), (/dimsize(d)/) )
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     d=4
     status = nf90_put_var( ncid_output(file), NCdimvarid(d), int(pft(1:dimsize(d)),2), (/1/), (/dimsize(d)/) )
     IF (status .NE. NF90_NOERR) WRITE(*,'(A)')  NF90_STRERROR(status)
     d=5
     status = nf90_put_var( ncid_output(file), NCdimvarid(d), hgt(1:dimsize(d)), (/1/), (/dimsize(d)/) )
     IF (status .NE. NF90_NOERR) WRITE(*,'(A)')  NF90_STRERROR(status)

  endif

  if (NCparallel) then
#ifndef HIDE_MPI
#ifndef HIDE_NETCDF4_MPI
     do v=1,NCnvar
        if (trim(NCnames(v,file)).ne.'') then
           ! force collective parallel access (independent access with unlimited time dimension does not work)
           status = nf90_var_par_access(ncid_output(file), NCvarid(v,file), nf90_collective)
           !        status = nf90_var_par_access(ncid_output(file), NCvarid(v,file), nf90_independent)
           if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
        end if
     end do
#endif
#endif
  end if

end subroutine netcdf_init

subroutine netcdf_write_2d(file, var, v, ncoffset, nccount)

  ! write static 2D data to NetCDF file
  use netcdf
  use typeSizes
  implicit none

  ! arguments
  integer, intent(in) :: file
  real(kind=4), intent(in) :: var(:,:)
  integer, intent(in) :: v
  integer,intent(in) :: ncoffset(:)
  integer,intent(in) :: nccount(:)

  ! local
  integer :: status
  real(kind=4), allocatable :: varin(:,:)

  ! generate temporary local data array
  allocate(varin(nccount(1),nccount(2)))
  varin = var

  ! add scale/offset to output variable if selected
  if ((NCscale_factor(v,file).ne.1.).or.(NCadd_offset(v,file).ne.0.)) then
     where(varin.ne.NCmissing)
        varin = (varin - NCadd_offset(v,file)) / NCscale_factor(v,file)
     endwhere
  endif

  ! limit variables to output variable range (probably only needed for byte/short)
  select case (NCtype(v,file))
  case ("byte")
     where(varin.ne.NCmissing)
        varin = max(min(varin,127.),-126.)
     end where
  case("short")
     where(varin.ne.NCmissing)
        varin = max(min(varin,32767.),-32766.)
     end where
  case ("int")
  case ("float")
  case ("double")
  end select

  ! add standard NetCDF missing values
  select case (NCtype(v,file))
  case ("byte")
     where(varin.eq.NCmissing)
        varin = real(NF90_FILL_BYTE,kind=4)
     end where
  case("short")
     where(varin.eq.NCmissing)
        varin = real(NF90_FILL_SHORT,kind=4)
     end where
  case ("int")
     where(varin.eq.NCmissing)
        varin = real(NF90_FILL_INT,kind=4)
     end where
  case ("float")
     where(varin.eq.NCmissing)
        varin = NF90_FILL_FLOAT
     end where
  case ("double")
     where(varin.eq.NCmissing)
        varin = real(NF90_FILL_DOUBLE,kind=4)
     end where
  end select

  ! write 2d data
  select case (NCtype(v,file))
  case ("byte")
     status = nf90_put_var( ncid_output(file), NCvarid(v,file), int(varin,kind=1), ncoffset, nccount )
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
  case("short")
     status = nf90_put_var( ncid_output(file), NCvarid(v,file), int(varin,kind=2), ncoffset, nccount )
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
  case ("int")
     status = nf90_put_var( ncid_output(file), NCvarid(v,file), int(varin,kind=4), ncoffset, nccount )
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
  case ("float")
     status = nf90_put_var( ncid_output(file), NCvarid(v,file), varin, ncoffset, nccount )
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
  case ("double")
     status = nf90_put_var( ncid_output(file), NCvarid(v,file), real(varin,kind=8), ncoffset, nccount )
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
  end select

  deallocate(varin)

end subroutine netcdf_write_2d

subroutine netcdf_write_3d(file, var, v, ncoffset, nccount)

  ! write static 3D data to NetCDF file
  use netcdf
  use typeSizes
  implicit none

  ! arguments
  integer, intent(in) :: file
  real(kind=4), intent(in) :: var(:,:,:)
  integer, intent(in) :: v
  integer,intent(in) :: ncoffset(:)
  integer,intent(in) :: nccount(:)

  ! local
  integer :: status
  real(kind=4), allocatable :: varin(:,:,:)

  ! generate temporary local data array
  allocate(varin(nccount(1),nccount(2),nccount(3)))
  varin = var

  ! add scale/offset to output variable if selected
  if ((NCscale_factor(v,file).ne.1.).or.(NCadd_offset(v,file).ne.0.)) then
     where(varin.ne.NCmissing)
        varin = (varin - NCadd_offset(v,file)) / NCscale_factor(v,file)
     endwhere
  endif

  ! limit variables to output variable range (probably only needed for byte/short)
  select case (NCtype(v,file))
  case ("byte")
     where(varin.ne.NCmissing)
        varin = max(min(varin,127.),-126.)
     end where
  case("short")
     where(varin.ne.NCmissing)
        varin = max(min(varin,32767.),-32766.)
     end where
  case ("int")
  case ("float")
  case ("double")
  end select

  ! add standard NetCDF missing values
  select case (NCtype(v,file))
  case ("byte")
     where(varin.eq.NCmissing)
        varin = real(NF90_FILL_BYTE,kind=4)
     end where
  case("short")
     where(varin.eq.NCmissing)
        varin = real(NF90_FILL_SHORT,kind=4)
     end where
  case ("int")
     where(varin.eq.NCmissing)
        varin = real(NF90_FILL_INT,kind=4)
     end where
  case ("float")
     where(varin.eq.NCmissing)
        varin = NF90_FILL_FLOAT
     end where
  case ("double")
     where(varin.eq.NCmissing)
        varin = real(NF90_FILL_DOUBLE,kind=4)
     end where
  end select

  ! write 3d data
  select case (NCtype(v,file))
  case ("byte")
     status = nf90_put_var( ncid_output(file), NCvarid(v,file), int(varin,kind=1), ncoffset, nccount )
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
  case("short")
     status = nf90_put_var( ncid_output(file), NCvarid(v,file), int(varin,kind=2), ncoffset, nccount )
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
  case ("int")
     status = nf90_put_var( ncid_output(file), NCvarid(v,file), int(varin,kind=4), ncoffset, nccount )
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
  case ("float")
     status = nf90_put_var( ncid_output(file), NCvarid(v,file), varin, ncoffset, nccount )
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
  case ("double")
     status = nf90_put_var( ncid_output(file), NCvarid(v,file), real(varin,kind=8), ncoffset, nccount )
     if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
  end select

  deallocate(varin)

end subroutine netcdf_write_3d

subroutine netcdf_write_time(file, timestep, timeval)

  ! collective call in parallel access to write unlimited time dimension. HDF5 layer
  ! does not allow to extend dimension in independent calls.

  use netcdf
  use typeSizes
  implicit none

  ! arguments
  integer, intent(in) :: file
  integer, intent(in) :: timestep
  real(kind=4), intent(in) :: timeval

  ! local
  integer status

  ! write time step
  status = nf90_put_var( ncid_output(file), NCdimvarid(6) , timeval, (/timestep/) )
  if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

end subroutine netcdf_write_time

subroutine netcdf_write(file, var, v, timestep, pft_pct, hgt_pct, ncoffset, nccount, &
     select_pft, select_hgt)

  ! writes NetCDF output. For PFTs and HGTs we need % subgrid-scale distribution to spatially redistribute
  ! can be independent call in parallel access with unlimited time dimension if time dimension has been previously
  ! extended (see above call) by collective call

  use netcdf
  use typeSizes
  implicit none

  ! arguments
  integer, intent(in) :: file
  real(kind=4), intent(in) :: var(:,:,:,:,:)
  integer, intent(in) :: v
  integer, intent(in) :: timestep
  real(kind=4), intent(in) :: pft_pct(:,:,:)
  real(kind=4), intent(in) :: hgt_pct(:,:,:)
  integer,intent(in) :: ncoffset(:)
  integer,intent(in) :: nccount(:)
  integer,intent(in), optional :: select_pft(:)
  integer,intent(in), optional :: select_hgt(:)

  ! local
  integer, parameter :: varndims = 5
  integer :: indimsize (varndims)
  integer :: outdimsize (varndims)
  integer status
  integer ndims
  integer k, d, e, i, j, p, h
  real(kind=4) dinv, dinv_minus_one
  real(kind=4) stddev
  real(kind=4) variance
  real(kind=4) mean
  real(kind=4) isum

  ! averaging quantitites
  real(kind=4), allocatable :: tempmean(:)
  real(kind=4), allocatable :: tempvar(:)

  ! averaged quantities
  real(kind=4), allocatable :: varin(:,:,:,:,:)
  real(kind=4), allocatable :: varout(:,:,:,:,:)

  integer, allocatable :: startvec(:)
  integer, allocatable :: countvec(:)

  ! determine input variable size
  indimsize = shape(var)

  ! determine output dimension sizes
  outdimsize(1) = nccount(1)
  outdimsize(2) = nccount(2)
  outdimsize(3) = nccount(3)
  outdimsize(4) = nccount(4)
  outdimsize(5) = nccount(5)

  ! check if dimension d gets averaged in output
  do d = 1,NCmaxdim
     if (NCdimave(d,file)) then
        outdimsize(d) = 1 ! becomes singular: n=1 (ensemble mean)
     endif
  end do

  ! check if dimension d gets compressed to two separate variables in output
  ! e.g. ensemble mean and standard deviation
  do d = 1,NCmaxdim
     if (NCdimcompress(d,file)) then
        outdimsize(d) = 2 ! becomes ensemble mean and ensemble standard deviation
     endif
  end do

  ! average along dimensions, start with ensemble members first, then pft and height if needed
  ! the code below is specific to the application here since we have weighted averaging for pft and hgt

  ! local copy of data array where we can work
  allocate(varin(indimsize(1),indimsize(2),indimsize(3),indimsize(4),indimsize(5)))
  varin = var

  ! ensemble averaging
  d = 1 
  if ((NCdimave(d,file).or.NCdimcompress(d,file)).and.(indimsize(d).ne.outdimsize(d))) then

     allocate(varout(outdimsize(1),indimsize(2),indimsize(3),indimsize(4),indimsize(5)))

     dinv = 1. / real(indimsize(d))
     dinv_minus_one = 1. / real(indimsize(d)-1)
     do h = 1,indimsize(5)
        do p = 1,indimsize(4)
           do j = 1, indimsize(3)
              do i = 1, indimsize(2)
                 if (varin(1,i,j,p,h).eq.NCmissing) then
                    varout(1,i,j,p,h) = NCmissing
                    if (NCdimcompress(d,file)) varout(2,i,j,p,h) = NCmissing
                 else
                    mean =  sum(varin(:,i,j,p,h)) * dinv
                    if (NCdimcompress(d,file)) then
                       stddev = (sum((varin(:,i,j,p,h)-mean)**2)*dinv_minus_one)**0.5
                       varout(2,i,j,p,h) = stddev
                    endif
                    varout(1,i,j,p,h) = mean
                 endif
              enddo
           enddo
        enddo
     enddo

     deallocate(varin)
     allocate(varin(outdimsize(1),indimsize(2),indimsize(3),indimsize(4),indimsize(5)))
     varin = varout
     deallocate(varout)

  endif

  ! x coordinate redistribution to high res grid (or grid averaging) goes here!
  d = 2
  if (indimsize(d).ne.outdimsize(d)) then

     allocate(varout(outdimsize(1),outdimsize(2),indimsize(3),indimsize(4),indimsize(5)))  

     if (indimsize(d).lt.outdimsize(d)) then
        ! subdivide grid
        do h = 1,indimsize(5)
           do p = 1,indimsize(4)
              do j = 1, indimsize(3)
                 do e = 1, outdimsize(1)

                    varout(e,:,j,p,h) = pack(spread(varin(e,:,j,p,h),1,outdimsize(d)/indimsize(d)),.true.)

                 enddo
              enddo
           enddo
        enddo
     endif

     if ((NCdimave(d,file).or.NCdimcompress(d,file)).and.(indimsize(d).gt.outdimsize(d))) then
        ! average grid
        dinv = 1. / real(indimsize(d))
        dinv_minus_one = 1. / real(indimsize(d)-1)

        do h = 1,indimsize(5)
           do p = 1,indimsize(4)
              do j = 1, indimsize(3)
                 do e = 1, outdimsize(1)
                    if (varin(e,1,j,p,h).eq.NCmissing) then
                       varout(e,1,j,p,h) = NCmissing
                       if (NCdimcompress(d,file)) varout(e,2,j,p,h) = NCmissing
                    else
                       mean =  sum(varin(e,:,j,p,h)) * dinv
                       if (NCdimcompress(d,file)) then
                          stddev = (sum((varin(e,:,j,p,h)-mean)**2)*dinv_minus_one)**0.5
                          varout(e,2,j,p,h) = stddev
                       endif
                       varout(e,1,j,p,h) = mean
                    endif
                 enddo
              enddo
           enddo
        enddo
     endif

     deallocate(varin)
     allocate(varin(outdimsize(1),outdimsize(2),indimsize(3),indimsize(4),indimsize(5)))
     varin = varout
     deallocate(varout)

  endif

  ! y coordinate redistribution to high res grid (or grid averaging) goes here!
  d = 3
  if (indimsize(d).ne.outdimsize(d)) then

     allocate(varout(outdimsize(1),outdimsize(2),outdimsize(3),indimsize(4),indimsize(5)))  

     if (indimsize(d).lt.outdimsize(d)) then
        ! subdivide grid
        do h = 1,indimsize(5)
           do p = 1,indimsize(4)
              do i = 1, outdimsize(2)
                 do e = 1, outdimsize(1)

                    varout(e,i,:,p,h) = pack(spread(varin(e,i,:,p,h),1,outdimsize(d)/indimsize(d)),.true.)

                 enddo
              enddo
           enddo
        enddo
     endif

     if ((NCdimave(d,file).or.NCdimcompress(d,file)).and.(indimsize(d).gt.outdimsize(d))) then
        ! average grid
        dinv = 1. / real(indimsize(d))
        dinv_minus_one = 1. / real(indimsize(d)-1)

        do h = 1,indimsize(5)
           do p = 1,indimsize(4)
              do i = 1, outdimsize(2)
                 do e = 1, outdimsize(1)
                    if (varin(e,i,1,p,h).eq.NCmissing) then
                       varout(e,i,1,p,h) = NCmissing
                       if (NCdimcompress(d,file)) varout(e,i,2,p,h) = NCmissing
                    else
                       mean =  sum(varin(e,i,:,p,h)) * dinv
                       if (NCdimcompress(d,file)) then
                          stddev = (sum((varin(e,i,:,p,h)-mean)**2)*dinv_minus_one)**0.5
                          varout(e,i,2,p,h) = stddev
                       endif
                       varout(e,i,1,p,h) = mean
                    endif
                 enddo
              enddo
           enddo
        enddo
     endif

     deallocate(varin)
     allocate(varin(outdimsize(1),outdimsize(2),outdimsize(3),indimsize(4),indimsize(5)))
     varin = varout
     deallocate(varout)
  endif

  ! pft averaging
  d = 4
  if (indimsize(d).ne.outdimsize(d)) then

     allocate(varout(outdimsize(1),outdimsize(2),outdimsize(3),outdimsize(4),indimsize(5)))  

     if (NCdimave(d,file).or.NCdimcompress(d,file)) then
        ! average or compress

        allocate(tempmean(indimsize(d)))
        allocate(tempvar(indimsize(d)))

        dinv = 1. / real(outdimsize(1))
        dinv_minus_one = 1. / max(real(outdimsize(1)-1),1.)

        do h = 1,indimsize(5)
           do j = 1, outdimsize(3)
              do i = 1, outdimsize(2)

                 if (varin(1,i,j,1,h).eq.NCmissing) then
                    ! input grid points are missing
                    varout(:,i,j,:,h) = NCmissing
                 else
                    if (sum(pft_pct(i,j,:)).lt.1.e-2) then
                       ! none of the selected PFT's present
                       varout(:,i,j,:,h) = NCmissing
                    else
                       ! selected PFT's are present
                       isum = 1./max(sum(pft_pct(i,j,:)),1.e-2)

                       if (outdimsize(1).gt.2) then
                          ! full ensemble output:
                          ! subtract mean of ensembles to have weighted centered ensembles
                          do p=1,indimsize(d)
                             tempmean(p) =  sum(varin(:,i,j,p,h)) * dinv
                             tempvar(p) =  sum((varin(:,i,j,p,h)-tempmean(p))**2)*dinv_minus_one
                          enddo

                          do e = 1, outdimsize(1)
                             varout(e,i,j,1,h) = sum((varin(e,i,j,:,h) - tempmean(:)) * pft_pct(i,j,:))            
                          enddo

                          ! rescale variance after adding weighted variances
                          mean = sum(varout(:,i,j,1,h)) * dinv
                          variance = sum((varout(:,i,j,1,h) - mean)**2)*dinv_minus_one

                          varout(:,i,j,1,h) = varout(:,i,j,1,h) * sqrt(sum(tempvar * pft_pct(i,j,:)) / &
                               max(variance,1.e-8))

                          ! add back mean weighted ensemble
                          varout(:,i,j,1,h) = varout(:,i,j,1,h) + sum(tempmean * pft_pct(i,j,:))

                       else if (outdimsize(1).eq.2) then
                          ! compressed ensemble output (2)
                          do p=1,indimsize(d)
                             varout(1,i,j,1,h) = sum(varin(1,i,j,:,h)*pft_pct(i,j,:))*isum
                             varout(2,i,j,1,h) = sum(varin(2,i,j,:,h)*pft_pct(i,j,:))*isum
                          enddo
                       else if (outdimsize(1).eq.1) then
                          ! mean output (1)
                          do p=1,indimsize(d)
                             varout(1,i,j,1,h) = sum(varin(1,i,j,:,h)*pft_pct(i,j,:))*isum
                          enddo
                       endif
                    endif
                 endif

              enddo
           enddo
        enddo

        deallocate(tempmean)
        deallocate(tempvar)

     else
        ! redistribute
        varout = NCmissing
        varout(:,:,:,select_pft,:) = varin(:,:,:,:,:)
     endif

     deallocate(varin)
     allocate(varin(outdimsize(1),outdimsize(2),outdimsize(3),outdimsize(4),indimsize(5)))
     varin = varout
     deallocate(varout)
  endif

  ! hgt averaging
  d = 5
  if (indimsize(d).ne.outdimsize(d)) then

     allocate(varout(outdimsize(1),outdimsize(2),outdimsize(3),outdimsize(4),outdimsize(5)))  

     if (NCdimave(d,file).or.NCdimcompress(d,file)) then
        ! average or compress

        allocate(tempmean(indimsize(d)))
        allocate(tempvar(indimsize(d)))

        dinv = 1. / real(outdimsize(1))
        dinv_minus_one = 1. / max(real(outdimsize(1)-1),1.)

        do p = 1,outdimsize(4)
           do j = 1, outdimsize(3)
              do i = 1, outdimsize(2)

                 if (varin(1,i,j,p,1).eq.NCmissing) then
                    ! input grid points are missing
                    varout(:,i,j,p,:) = NCmissing
                 else
                    if (sum(hgt_pct(i,j,:)).lt.1.e-2) then
                       ! none of the selected HGT's present
                       varout(:,i,j,p,:) = NCmissing
                    else
                       ! selected HGT's are present
                       isum = 1./max(sum(hgt_pct(i,j,:)),1.e-2)

                       if (outdimsize(1).gt.2) then
                          ! full ensemble output
                          ! subtract mean of ensembles to have weighted centered ensembles
                          do h=1,indimsize(d)
                             tempmean(h) =  sum(varin(:,i,j,p,h)) * dinv
                             tempvar(h) =  sum((varin(:,i,j,p,h)-tempmean(h))**2)*dinv_minus_one
                          enddo

                          do e = 1, outdimsize(1)
                             varout(e,i,j,p,1) = sum((varin(e,i,j,p,:) - tempmean(:)) * hgt_pct(i,j,:))            
                          enddo

                          ! rescale variance after adding weighted variances
                          mean = sum(varout(:,i,j,p,1)) * dinv
                          variance = sum((varout(:,i,j,p,1) - mean)**2)*dinv_minus_one

                          varout(:,i,j,p,1) = varout(:,i,j,p,1) * sqrt(sum(tempvar * hgt_pct(i,j,:)) / &
                               max(variance,1.e-8))

                          ! add back mean weighted ensemble
                          varout(:,i,j,p,1) = varout(:,i,j,p,1) + sum(tempmean * hgt_pct(i,j,:))

                       else if (outdimsize(1).eq.2) then
                          ! compressed ensemble output
                          varout(1,i,j,p,1) = sum(varin(1,i,j,p,:)*hgt_pct(i,j,:))*isum
                          varout(2,i,j,p,1) = sum(varin(2,i,j,p,:)*hgt_pct(i,j,:))*isum
                       else if (outdimsize(1).eq.1) then
                          ! single ensemble output
                          varout(1,i,j,p,1) = sum(varin(1,i,j,p,:)*hgt_pct(i,j,:))*isum
                       endif
                    endif
                 endif

              enddo
           enddo
        enddo

        deallocate(tempmean)
        deallocate(tempvar)

     else
        ! redistribute
        varout = NCmissing
        varout(:,:,:,:,select_hgt) = varin(:,:,:,:,:)        
     endif

     deallocate(varin)
     allocate(varin(outdimsize(1),outdimsize(2),outdimsize(3),outdimsize(4),outdimsize(5)))
     varin = varout
     deallocate(varout)

  endif

  ndims = sum(NCdim(v,:,file))
  allocate(startvec(ndims))
  allocate(countvec(ndims))

  ! space dimensions start/count
  k=1
  do d = 1, varndims
     if (NCdim(v,d,file).eq.1) then
        startvec(k) = ncoffset(d)
        countvec(k) = outdimsize(d)
        k=k+1
     endif
  enddo

  ! time start/count (last dimension)
  startvec(k) = timestep
  countvec(k) = 1

  ! add scale/offset to output variable if selected
  if ((NCscale_factor(v,file).ne.1.).or.(NCadd_offset(v,file).ne.0.)) then
     where(varin.ne.NCmissing)
        varin = (varin - NCadd_offset(v,file)) / NCscale_factor(v,file)
     endwhere
  endif

  ! limit variables to output variable range (probably only needed for byte/short)
  select case (NCtype(v,file))
  case ("byte")
     where(varin.ne.NCmissing)
        varin = max(min(varin,127.),-126.)
     end where
  case("short")
     where(varin.ne.NCmissing)
        varin = max(min(varin,32767.),-32766.)
     end where
  case ("int")
  case ("float")
  case ("double")
  end select

  ! add standard NetCDF missing values
  select case (NCtype(v,file))
  case ("byte")
     where(varin.eq.NCmissing)
        varin = real(NF90_FILL_BYTE,kind=4)
     end where
  case("short")
     where(varin.eq.NCmissing)
        varin = real(NF90_FILL_SHORT,kind=4)
     end where
  case ("int")
     where(varin.eq.NCmissing)
        varin = real(NF90_FILL_INT,kind=4)
     end where
  case ("float")
     where(varin.eq.NCmissing)
        varin = NF90_FILL_FLOAT
     end where
  case ("double")
     where(varin.eq.NCmissing)
        varin = real(NF90_FILL_DOUBLE,kind=4)
     end where
  end select

  if (NCdimcompress(1,file)) then
     ! assume that only first dimension can get compressed
     select case (NCtype(v,file))
     case ("byte")
        status = nf90_put_var( ncid_output(file), NCvarid(v,file), nint(varin(1,:,:,:,:),1), startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
        status = nf90_put_var( ncid_output(file), NCvarid((NCnvar-NCnauxvar)/2+v,file), &
             nint(varin(2,:,:,:,:),1), startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     case("short")
        status = nf90_put_var( ncid_output(file), NCvarid(v,file), &
             nint(varin(1,:,:,:,:),2), startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
        status = nf90_put_var( ncid_output(file), NCvarid((NCnvar-NCnauxvar)/2+v,file), &
             nint(varin(2,:,:,:,:),2), startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     case ("int")
        status = nf90_put_var( ncid_output(file), NCvarid(v,file), &
             nint(varin(1,:,:,:,:),4), startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
        status = nf90_put_var( ncid_output(file), NCvarid((NCnvar-NCnauxvar)/2+v,file), &
             nint(varin(2,:,:,:,:),4), startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     case ("float")
        status = nf90_put_var( ncid_output(file), NCvarid(v,file), &
             varin(1,:,:,:,:), startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
        status = nf90_put_var( ncid_output(file), NCvarid((NCnvar-NCnauxvar)/2+v,file), &
             varin(2,:,:,:,:), startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     case ("double")
        status = nf90_put_var( ncid_output(file), NCvarid(v,file), &
             real(varin(1,:,:,:,:),kind=8), startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
        status = nf90_put_var( ncid_output(file), NCvarid((NCnvar-NCnauxvar)/2+v,file), &
             real(varin(2,:,:,:,:),kind=8), startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     end select
  else
     select case (NCtype(v,file))
     case ("byte")
        status = nf90_put_var( ncid_output(file), NCvarid(v,file), int(varin,1), startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     case("short")
        status = nf90_put_var( ncid_output(file), NCvarid(v,file), int(varin,2), startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     case ("int")
        status = nf90_put_var( ncid_output(file), NCvarid(v,file), int(varin,4), startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     case ("float")
        status = nf90_put_var( ncid_output(file), NCvarid(v,file), varin, startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     case ("double")
        status = nf90_put_var( ncid_output(file), NCvarid(v,file), real(varin,kind=8), startvec, countvec )
        if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)
     end select
  endif

  deallocate(startvec)
  deallocate(countvec)

  deallocate(varin)

end subroutine netcdf_write

subroutine netcdf_close(file)

  use netcdf
  use typeSizes

  implicit none
  
  ! arguments
  integer, intent(in) :: file

  ! local
  integer status

  ! close file
  status = nf90_close(ncid_output(file))
  if (status .ne. NF90_NOERR) write(*,'(A)')  NF90_STRERROR(status)

end subroutine netcdf_close

end module netcdf_module
