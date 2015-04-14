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

module output_module

contains

  subroutine init_output

    use netcdf
    use typeSizes

    use driver_module
    use exit_module
    use time_module
    use netcdf_module

#ifndef HIDE_MPI
  use mpi
#endif

    implicit none

    ! local variables
    integer :: f,i
#ifndef HIDE_MPI
#ifndef HIDE_NETCDF4_MPI
    integer :: mpistatus(MPI_STATUS_SIZE)
    integer :: group_world, group_io
    integer :: count_io
    integer :: c, proc
    logical :: l
    integer, allocatable :: procs_io(:)
#endif
#endif

    ! define total number of NetCDF files
    NCnfiles = 3

    ! define total NetCDF variables
    NCnvar = noutput_state + noutput_param + noutput_force + noutput_obser ! number of variables for NetCDF file

    if (trim(geographic%type).eq."latitude_longitude") then
       NCnauxvar = 4   ! PFT_PCT + HGT_PCT + Z + GRID_MAPPING
    else
       NCnauxvar = 6   ! LON + LAT +  PFT_PCT + HGT_PCT + Z + GRID_MAPPING
    end if

    ! double standard model variables to account with standard deviation 
    ! for variables where ensemble compression or
    ! the observation error is chosen
    NCnvar = NCnvar * 2 + NCnauxvar

    ! construct output NetCDF names + dimensions
    call netcdf_allocate

    if (verbose) write(*,'(A)') 'NetCDF arrays allocated.'

    ! check whether we have Parallel I/O with NetCDF-4 / HDF5
    if ((nxprocs.gt.1).or.(nyprocs.gt.1)) then
       NCparallel = .true.
#ifdef HIDE_MPI
       write(*,'(A)') 'MPI needs to be used for parallel I/O on a distributed grid'
      call model_exit(-7934)
#endif
#ifdef HIDE_NETCDF4_MPI
      write(*,'(A)') 'NetCDF-4 needs to be linked for parallel I/O on a distributed grid'
      call model_exit(-7935)
#endif
    else
       NCparallel = .false.
    endif


    ! generate I/O communicator for selected output processes
#ifndef HIDE_MPI
#ifndef HIDE_NETCDF4_MPI
    if (NCparallel) then
       if (isprediction) then
          c = 1
       else
          c = 0
       end if
       call MPI_ALLREDUCE(c, count_io, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
       allocate(procs_io(count_io))
       ! now gather all process numbers which run a prediction (and thus output data)
       if (mpiid.eq.0) then
          c = 1
          if (isprediction) then
             procs_io(c) = mpiid
             c = c + 1
          end if
          do proc=2,numprocs
             call MPI_RECV(l,1,MPI_LOGICAL,proc-1,proc,MPI_COMM_WORLD,mpistatus,mpierr)
             if (l) then
                procs_io(c) = proc-1
                c = c + 1
             end if
          end do
       else
          call MPI_SEND(isprediction,1,MPI_LOGICAL,0,mpiid+1,MPI_COMM_WORLD,mpierr)
       end if

       call MPI_BCAST(procs_io,count_io,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

       call MPI_Comm_group(MPI_COMM_WORLD, group_world, mpierr)
       call MPI_Group_incl(group_world, count_io, procs_io, group_io, mpierr)  
       call MPI_Comm_create(MPI_COMM_WORLD, group_io, NCcomm, mpierr)
       deallocate(procs_io)
    end if
#endif
#endif

    ! find out which dimensions are averaged (scalar)
    NCdimave(:,:) = .false.
    if (averageens) NCdimave(1,:) = .true.
    if (averagepft) NCdimave(4,:) = .true.
    if (averagehgt) NCdimave(5,:) = .true.

    ! find out which dimensions are compressed (mean+stddev)
    NCdimcompress(:,:) = .false.

    ! compress ensemble dimension to mean and standard deviation on output
    if (compressens) NCdimcompress(1,1) = .true.

    ! simulator only has averaged quantities
    NCdimcompress(1,2) = .false.

    ! observations have mean and standard deviation by default on output
    NCdimcompress(1,3) = .true.

    ! reset variable names, long names and units
    NCnames(:,:) = ''
    NClongnames(:,:) = ''
    NCunits(:,:) = ''

    ! set file presence to off
    NChasfile(:) = .false.

    ! reset dimension assignment
    NCdim(:,:,:)=0

    ! reset variable ID's
    NCvarid(:,:) = -1

    ! define auxiliary variables (the same for all files)
    do f=1,NCnfiles
       if (trim(geographic%type).eq."latitude_longitude") then
          ! PFT_PCT + HGT_PCT + Z
          NCdim(1,:,f) = (/0,1,1,1,0,0/)
          NCdim(2,:,f) = (/0,1,1,0,1,0/)
          NCdim(3,:,f) = (/0,1,1,0,0,0/)
          NCdim(4,:,f) = (/0,0,0,0,0,0/)
       else
          ! LON + LAT +  PFT_PCT + HGT_PCT + Z
          NCdim(1,:,f) = (/0,1,1,0,0,0/)
          NCdim(2,:,f) = (/0,1,1,0,0,0/)
          NCdim(3,:,f) = (/0,1,1,1,0,0/)
          NCdim(4,:,f) = (/0,1,1,0,1,0/)
          NCdim(5,:,f) = (/0,1,1,0,0,0/)
          NCdim(6,:,f) = (/0,0,0,0,0,0/)
       end if
       if (trim(geographic%type).eq."latitude_longitude") then
          ! PFT_PCT + HGT_PCT + Z
          NCnames(1,f) = "PFT_PCT"
          NCnames(2,f) = "HGT_PCT"
          NCnames(3,f) = "Z"
          NCnames(4,f) = "grid_mapping"
          NClongnames(1,f) = "Percent Cover of Plant Functional Type"
          NClongnames(2,f) = "Percent Cover of Elevation Class"
          NClongnames(3,f) = "Elevation above sea level"
          NClongnames(4,f) = ""
          NCunits(1,f) = "%"
          NCunits(2,f) = "%"
          NCunits(3,f) = "m"
          NCunits(4,f) = ""
          NCscale_factor(1,f) = 0.01
          NCscale_factor(2,f) = 0.01
          NCscale_factor(3,f) = 1.0
          NCscale_factor(4,f) = 1.0
          NCadd_offset(1,f) = 0.0
          NCadd_offset(2,f) = 0.0
          NCadd_offset(3,f) = 0.0
          NCadd_offset(4,f) = 0.0
          NCtype(1,f) = "byte"
          NCtype(2,f) = "byte"
          NCtype(3,f) = "float"
          NCtype(4,f) = "char"
       else
          ! LON + LAT +  PFT_PCT + HGT_PCT + Z
          NCnames(1,f) = "lon"
          NCnames(2,f) = "lat"
          NCnames(3,f) = "PFT_PCT"
          NCnames(4,f) = "HGT_PCT"
          NCnames(5,f) = "Z"
          NCnames(6,f) = "grid_mapping"
          NClongnames(1,f) = "longitude"
          NClongnames(2,f) = "latitude"
          NClongnames(3,f) = "Percent Cover of Plant Functional Type"
          NClongnames(4,f) = "Percent Cover of Elevation Class"
          NClongnames(5,f) = "Elevation above sea level"
          NClongnames(6,f) = ""
          NCunits(1,f) = "degrees_east"
          NCunits(2,f) = "degrees_north"
          NCunits(3,f) = "%"
          NCunits(4,f) = "%"
          NCunits(5,f) = "m"
          NCunits(6,f) = ""
          NCscale_factor(1,f) = 1.0
          NCscale_factor(2,f) = 1.0
          NCscale_factor(3,f) = 0.01
          NCscale_factor(4,f) = 0.01
          NCscale_factor(5,f) = 1.0
          NCscale_factor(6,f) = 1.0
          NCadd_offset(1,f) = 0.0
          NCadd_offset(2,f) = 0.0
          NCadd_offset(3,f) = 0.0
          NCadd_offset(4,f) = 0.0
          NCadd_offset(5,f) = 0.0
          NCadd_offset(6,f) = 0.0
          NCtype(1,f) = "float"
          NCtype(2,f) = "float"
          NCtype(3,f) = "byte"
          NCtype(4,f) = "byte"
          NCtype(5,f) = "float"
          NCtype(6,f) = "char"
       end if
    end do

    ! define initialization/open/close file per write step to off
    f = 1

    NChasfile(f) = ((noutput_state+noutput_force+noutput_param).gt.0)

    if (NChasfile(f)) then
       if (noutput_state.gt.0) then
          NCnames(1+NCnauxvar:noutput_state+NCnauxvar,f) = state_names(output_state)
          NClongnames(1+NCnauxvar:noutput_state+NCnauxvar,f) = state_longnames(output_state)
          NCunits(1+NCnauxvar:noutput_state+NCnauxvar,f) = state_units(output_state)
          NCscale_factor(1+NCnauxvar:noutput_state+NCnauxvar,f) = state_scale(output_state)
          NCadd_offset(1+NCnauxvar:noutput_state+NCnauxvar,f) = state_offset(output_state)
          NCtype(1+NCnauxvar:noutput_state+NCnauxvar,f) = state_type(output_state)
       endif
       if (noutput_param.gt.0) then
          NCnames(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,f) = param_names(output_param)
          NClongnames(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,f) = param_longnames(output_param)
          NCunits(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,f) = param_units(output_param)
          NCscale_factor(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,f) = param_scale(output_param)
          NCadd_offset(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,f) = param_offset(output_param)
          NCtype(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,f) = param_type(output_param)
       endif
       if (noutput_force.gt.0) then
          NCnames(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,f) = &
               force_names(output_force)
          NClongnames(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,f) = &
               force_longnames(output_force)
          NCunits(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,f) = &
               force_units(output_force)
          NCscale_factor(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,f) = &
               force_scale(output_force)
          NCadd_offset(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,f) = &
               force_offset(output_force)
          NCtype(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,f) = &
               force_type(output_force)
       endif

       if (compressens) then
          do i=1,(NCnvar-NCnauxvar)/2
             if (trim(NCnames(i+NCnauxvar,f)).ne.'') then 
                NCnames((NCnvar-NCnauxvar)/2+i+NCnauxvar,f) = trim(NCnames(i+NCnauxvar,f))//'_stddev'
                NClongnames((NCnvar-NCnauxvar)/2+i+NCnauxvar,f) = &
                     trim(NClongnames(i+NCnauxvar,f))//' (ensemble standard deviation)'
                NCunits((NCnvar-NCnauxvar)/2+i+NCnauxvar:NCnvar,f) = NCunits(i+NCnauxvar,f)
                NCscale_factor((NCnvar-NCnauxvar)/2+i+NCnauxvar:NCnvar,f) = NCscale_factor(i+NCnauxvar,f)
                NCadd_offset((NCnvar-NCnauxvar)/2+i+NCnauxvar:NCnvar,f) = NCadd_offset(i+NCnauxvar,f)
                NCtype((NCnvar-NCnauxvar)/2+i+NCnauxvar:NCnvar,f) = NCtype(i+NCnauxvar,f)
             endif
          end do
       endif

    endif

    ! simulator has no ensemble thus no compressed ensemble
    if (simulator) then
       f = 2

       NChasfile(f) = ((noutput_state+noutput_force+noutput_param).gt.0)

       if (NChasfile(f)) then

          if (noutput_state.gt.0) then
             NCnames(1+NCnauxvar:noutput_state+NCnauxvar,f) = state_names(output_state)
             NClongnames(1+NCnauxvar:noutput_state+NCnauxvar,f) = state_longnames(output_state)
             NCunits(1+NCnauxvar:noutput_state+NCnauxvar,f) = state_units(output_state)
             NCscale_factor(1+NCnauxvar:noutput_state+NCnauxvar,f) = state_scale(output_state)
             NCadd_offset(1+NCnauxvar:noutput_state+NCnauxvar,f) = state_offset(output_state)
             NCtype(1+NCnauxvar:noutput_state+NCnauxvar,f) = state_type(output_state)
          endif
          if (noutput_param.gt.0) then
             NCnames(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,f) = param_names(output_param)
             NClongnames(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,f) = param_longnames(output_param)
             NCunits(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,f) = param_units(output_param)
             NCscale_factor(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,f) = param_scale(output_param)
             NCadd_offset(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,f) = param_offset(output_param)
             NCtype(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,f) = param_type(output_param)
          endif
          if (noutput_force.gt.0) then
             NCnames(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,f) = &
                  force_names(output_force)
             NClongnames(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,f) = &
                  force_longnames(output_force)
             NCunits(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,f) = &
                  force_units(output_force)
             NCscale_factor(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,f) = &
                  force_scale(output_force)
             NCadd_offset(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,f) = &
                  force_offset(output_force)
             NCtype(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,f) = &
                  force_type(output_force)
          endif

       endif

    endif

    if (saveobs) then
       f = 3

       NChasfile(f) = (noutput_obser.gt.0)

       if (NChasfile(f)) then

          if (noutput_obser.gt.0) then
             NCnames(1+NCnauxvar:noutput_obser+NCnauxvar,f) = obser_names(output_obser)
             NClongnames(1+NCnauxvar:noutput_obser+NCnauxvar,f) = obser_longnames(output_obser)
             NCunits(1+NCnauxvar:noutput_obser+NCnauxvar,f) =  obser_units(output_obser)
             NCscale_factor(1+NCnauxvar:noutput_obser+NCnauxvar,f) =  obser_scale(output_obser)
             NCadd_offset(1+NCnauxvar:noutput_obser+NCnauxvar,f) =  obser_offset(output_obser)
             NCtype(1+NCnauxvar:noutput_obser+NCnauxvar,f) =  obser_type(output_obser)
          endif

          do i=1,(NCnvar-NCnauxvar)/2
             if (trim(NCnames(i+NCnauxvar,f)).ne.'') then 
                NCnames((NCnvar-NCnauxvar)/2+i+NCnauxvar,f) = trim(NCnames(i+NCnauxvar,f))//'_stddev'
                NClongnames((NCnvar-NCnauxvar)/2+i+NCnauxvar,f) = &
                     trim(NClongnames(i+NCnauxvar,f))//' (observation uncertainty)'
                NCunits((NCnvar-NCnauxvar)/2+i+NCnauxvar:NCnvar,f) = NCunits(i+NCnauxvar,f)
                NCscale_factor((NCnvar-NCnauxvar)/2+i+NCnauxvar:NCnvar,f) = NCscale_factor(i+NCnauxvar,f)
                NCadd_offset((NCnvar-NCnauxvar)/2+i+NCnauxvar:NCnvar,f) = NCadd_offset(i+NCnauxvar,f)
                NCtype((NCnvar-NCnauxvar)/2+i+NCnauxvar:NCnvar,f) = NCtype(i+NCnauxvar,f)
             endif
          end do

       endif

    endif

    NCdimnames(1,:) = 'ens'
    if (trim(geographic%type).eq."latitude_longitude") then
       NCdimnames(2,:) = 'lon'
       NCdimnames(3,:) = 'lat'
    else
       NCdimnames(2,:) = 'x'
       NCdimnames(3,:) = 'y'
    end if
    NCdimnames(4,:) = 'pft'
    NCdimnames(5,:) = 'hgt'
    NCdimnames(6,:) = 'time'

    NCdimvarnames(1,:) = 'ens'
    if (trim(geographic%type).eq."latitude_longitude") then
       NCdimvarnames(2,:) = 'lon'
       NCdimvarnames(3,:) = 'lat'
    else
       NCdimvarnames(2,:) = 'x'
       NCdimvarnames(3,:) = 'y'
    end if
    NCdimvarnames(4,:) = 'pft'
    NCdimvarnames(5,:) = 'hgt'
    NCdimvarnames(6,:) = 'time'

    NCdimlongnames(1,:) = 'ensemble members'
    if (trim(geographic%type).eq."latitude_longitude") then
       NCdimlongnames(2,:) = 'longitude'
       NCdimlongnames(3,:) = 'latitude'
    else
       NCdimlongnames(2,:) = 'x coordinate of projection'
       NCdimlongnames(3,:) = 'y coordinate of projection'
    end if
    NCdimlongnames(4,:) = 'plant functional type'
    NCdimlongnames(5,:) = 'elevation above sea level'
    NCdimlongnames(6,:) = 'time'

    NCdimunits(1,:) = '1'
    if (trim(geographic%type).eq."latitude_longitude") then
       NCdimunits(2,:) = 'degrees_east'
       NCdimunits(3,:) = 'degrees_north'
    else
       NCdimunits(2,:) = 'x coordinate units of projection'
       NCdimunits(3,:) = 'y coordinate units of projection'
    end if
    NCdimunits(4,:) = '1'
    NCdimunits(5,:) = 'm'
    NCdimunits(6,:) = 'days since'

    f = 1
    if (NChasfile(f)) then
       ! states (nrens,geographic%nx,geographic%ny,npft,nhgt,time)
       if (noutput_state.gt.0) then
          NCdim(1+NCnauxvar:noutput_state+NCnauxvar,:,f) = 1
       endif
       ! parameters (nrens,npft,time)
       if (noutput_param.gt.0) then     
          NCdim(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,1,f)=1
          NCdim(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,4,f)=1
          NCdim(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,6,f)=1
       endif
       ! forcings (nrens,geographic%nx,geographic%ny,nhgt,time)
       if (noutput_force.gt.0) then
          NCdim(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,:,f)=1
          NCdim(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,4,f)=0
       endif

       if (compressens) then
          NCdim((NCnvar-NCnauxvar)/2+1+NCnauxvar:NCnvar,:,f) = NCdim(1+NCnauxvar:(NCnvar-NCnauxvar)/2+NCnauxvar,:,f)
       endif
    endif

    if (simulator) then
       f=2
       if (NChasfile(f)) then
          ! states (nrens,geographic%nx,geographic%ny,npft,nhgt,time)
          if (noutput_state.gt.0) then
             NCdim(1+NCnauxvar:noutput_state+NCnauxvar,:,f) = 1
          endif
          ! parameters (nrens,npft,time)
          if (noutput_param.gt.0) then     
             NCdim(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,1,f)=1
             NCdim(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,4,f)=1
             NCdim(noutput_state+1+NCnauxvar:noutput_state+noutput_param+NCnauxvar,6,f)=1
          endif
          ! forcings (nrens,geographic%nx,geographic%ny,nhgt,time)
          if (noutput_force.gt.0) then
             NCdim(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,:,f)=1
             NCdim(noutput_state+noutput_param+1+NCnauxvar:noutput_state+noutput_param+noutput_force+NCnauxvar,4,f)=0
          endif
       endif
    endif

    if (saveobs) then
       f=3
       if (NChasfile(f)) then
          ! observations (2,geographic%nx,geographic%ny,time)
          if (noutput_obser.gt.0) then
             NCdim(1+NCnauxvar:noutput_obser+NCnauxvar,1,f)=1
             NCdim(1+NCnauxvar:noutput_obser+NCnauxvar,2,f)=1
             NCdim(1+NCnauxvar:noutput_obser+NCnauxvar,3,f)=1
             NCdim(1+NCnauxvar:noutput_obser+NCnauxvar,6,f)=1
          endif

          NCdim((NCnvar-NCnauxvar)/2+1+NCnauxvar:NCnvar,:,f) = NCdim(1+NCnauxvar:(NCnvar-NCnauxvar)/2+NCnauxvar,:,f)
       endif
    endif

  end subroutine init_output

  subroutine open_output

    use driver_module
    use time_module
    use netcdf_module

    implicit none

    ! local variables
    integer :: f, n
    character(len=8) :: sdate
    integer :: ntout

    if (isprediction) then

       if (((timestep+timestep_start-1).eq.timestep_start).or. &
            (newyear.and.(filelength.eq.3)).or.(newmonth.and.(filelength.eq.2)).or. &
            (newday.and.(filelength.eq.1))) then

          ! create file and enter define mode
          select case (filelength)
          case (1)  
             write(unit=sdate,fmt='(i8)') year*10000+month*100+day
             ntout = 86400 / dtout
          case (2)
             write(unit=sdate,fmt='(i6)') year*100+month
             ntout = calc_daysofmonth(year,month) * 86400 / dtout
          case (3)
             write(unit=sdate,fmt='(i4)') year
             ntout = calc_daysofyear(year) * 86400 / dtout
          case default
             stop
          end select

          ! define filenames
          NCfilenames(1) = trim(outdir)//trim(geographic%name)//'.analysis.'//trim(sdate)//'.nc'
          NCfilenames(2) = trim(outdir)//trim(geographic%name)//'.simulator.'//trim(sdate)//'.nc'
          NCfilenames(3) = trim(outdir)//trim(geographic%name)//'.observations.'//trim(sdate)//'.nc'

          ! initialize files
          do f=1,3

             if (NChasfile(f)) then

                ! initialize NetCDF output file

                ! initialize only sub-domains that have prediction, or collective initialization
                ! from all processes if parallel output is used
                if (NCparallel) then
                   if (verbose) write(*,'(A,I5,A,A)') 'Parallel Initialize output file ',mpiid,' : ',trim(NCfilenames(f))
                   ! single file for multiple processes
                   ! initialize full domain at once
                   call netcdf_init(f, ensindex, xcoord_full, ycoord_full, pftindex, hgtlevels, &
                        year, package_name, package_version, nodata, ntout, geographic%type, geographic%optargs)
                else
                   write(*,'(A,I5,A,A)') 'Serial Initialize output file ',mpiid,' : ',trim(NCfilenames(f))
                   ! one file per process
                   ! initalize one file per region
                   call netcdf_init(f, ensindex, xcoord_out, ycoord_out, pftindex, hgtlevels, &
                        year, package_name, package_version, nodata, ntout, geographic%type, geographic%optargs)
                endif

                ! write static pft and hgt distribution and elevation, and optionally irregular lon/lat fields
                n = 1
                if (trim(geographic%type).ne."latitude_longitude") then
                   call netcdf_write_2d(f, lon, n, &
                        (/xcoord_offset,ycoord_offset/),(/geographic_out%nx,geographic_out%ny/))
                   n=n+1
                   call netcdf_write_2d(f, lat, n, &
                        (/xcoord_offset,ycoord_offset/),(/geographic_out%nx,geographic_out%ny/))
                   n=n+1
                end if
                call netcdf_write_3d(f, pft_pct_out, n, &
                     (/xcoord_offset,ycoord_offset,1/),(/geographic_out%nx,geographic_out%ny,npft/))
                n=n+1
                call netcdf_write_3d(f, hgt_pct_out, n, &
                     (/xcoord_offset,ycoord_offset,1/),(/geographic_out%nx,geographic_out%ny,nhgt/))
                n=n+1
                call netcdf_write_2d(f, hgt_out, n, &
                     (/xcoord_offset,ycoord_offset/),(/geographic_out%nx,geographic_out%ny/))

             endif

          enddo

       endif

    endif

  end subroutine open_output

  subroutine write_output

    use driver_module
    use time_module
    use netcdf_module

#ifndef HIDE_MPI
    use mpi
#endif

    implicit none

    ! arguments

    ! local variables
    integer :: f, n, t

    if (isprediction) then

       if (mod(timestep*dtmod,dtout).eq.0) then

          select case (filelength)
          case (1)  
             t = (timestep-timestep_daystart)*dtmod/dtout + 1
          case (2)
             t = (timestep-timestep_monthstart)*dtmod/dtout + 1
          case (3)
             t = (timestep-timestep_yearstart)*dtmod/dtout + 1
          case default
             t = timestep*dtmod/dtout
          end select

          do f=1,3

             ! only save stuff if there are selected spatially distributed variables to save, except
             ! for the first process in distributed runs or during multiple region assimilation experiments.

             ! in parallel mode saving happens independently, so not all processes need to be involved
             if (NChasfile(f)) then

                ! write time step (extend unlimited time dimension) and time value (collective)
                call netcdf_write_time(f,t,real(t-1))

                ! write states & parameters & forcings to output NetCDF file

                ! states (spatially distributed)
                do n=1,noutput_state 

                   if (averagehgt) then
                      if (averagepft) then

                         ! average only selected hgt's and selected pft's
                         if (f.eq.1) call netcdf_write(f, &
                              state(:,:,:,:,:,output_state(n)), &
                              n+NCnauxvar, t, &
                              pft_pct_out(:,:,select_pft), &
                              hgt_pct_out(:,:,select_hgt), &
                              (/1,xcoord_offset,ycoord_offset,1,1/), &
                              (/nrens,geographic_out%nx,geographic_out%ny,npft,nhgt/), &
                              select_pft=select_pft, &
                              select_hgt=select_hgt) 
                         if (f.eq.2) call netcdf_write(f, &
                              state_sim(:,:,:,:,:,output_state(n)), &
                              n+NCnauxvar, t, &
                              pft_pct_out(:,:,select_pft), &
                              hgt_pct_out(:,:,select_hgt), &
                              (/1,xcoord_offset,ycoord_offset,1,1/), &
                              (/1,geographic_out%nx,geographic_out%ny,npft,nhgt/))     
                      else
                         ! average only selected hgt's and output all pft's
                         if (f.eq.1) call netcdf_write(f, &
                              state(:,:,:,:,:,output_state(n)), &
                              n+NCnauxvar, t, &
                              pft_pct_out,hgt_pct_out(:,:,select_hgt), &
                              (/1,xcoord_offset,ycoord_offset,1,1/), &
                              (/nrens,geographic_out%nx,geographic_out%ny,npft,nhgt/), &
                              select_pft=select_pft) 
                         if (f.eq.2) call netcdf_write(f, &
                              state_sim(:,:,:,:,:,output_state(n)), &
                              n+NCnauxvar, t, &
                              pft_pct_out,hgt_pct_out(:,:,select_hgt), &
                              (/1,xcoord_offset,ycoord_offset,1,1/), &
                              (/1,geographic_out%nx,geographic_out%ny,npft,nhgt/), &
                              select_pft=select_pft)     
                      endif
                   else
                      if (averagepft) then
                         ! average only selected pft's and output all hgt's
                         if (f.eq.1) call netcdf_write(f, &
                              state(:,:,:,:,:,output_state(n)), &
                              n+NCnauxvar, t, &
                              pft_pct_out(:,:,select_pft),&
                              hgt_pct_out, &
                              (/1,xcoord_offset,ycoord_offset,1,1/), &
                              (/nrens,geographic_out%nx,geographic_out%ny,nselect_pft,nhgt/),&
                              select_hgt=select_hgt) 
                         if (f.eq.2) call netcdf_write(f, &
                              state_sim(:,:,:,:,:,output_state(n)), &
                              n+NCnauxvar, t, &
                              pft_pct_out(:,:,select_pft), &
                              hgt_pct_out,&
                              (/1,xcoord_offset,ycoord_offset,1,1/), &
                              (/1,geographic_out%nx,geographic_out%ny,npft,nhgt/), &
                              select_hgt=select_hgt)     
                      else
                         ! output all pft's and hgt's
                         if (f.eq.1) call netcdf_write(f,&
                              state(:,:,:,:,:,output_state(n)), &
                              n+NCnauxvar, t, &
                              pft_pct_out,hgt_pct_out, &
                              (/1,xcoord_offset,ycoord_offset,1,1/), &
                              (/nrens,geographic_out%nx,geographic_out%ny,npft,nhgt/),&
                              select_pft=select_pft, & 
                              select_hgt=select_hgt) 
                         if (f.eq.2) call netcdf_write(f,&
                              state_sim(:,:,:,:,:,output_state(n)),&
                              n+NCnauxvar, t, &
                              pft_pct_out,hgt_pct_out, &
                              (/1,xcoord_offset,ycoord_offset,1,1/), &
                              (/1,geographic_out%nx,geographic_out%ny,npft,nhgt/), &
                              select_pft=select_pft, &
                              select_hgt=select_hgt)     
                      endif
                   endif
                enddo

                ! parameters (always save all pft's since this is a global parameter set 
                ! not necessarily reflecting the local pft distribution)
                do n=1,noutput_param 
                   if (averagepft) then
                      if (f.eq.1) call netcdf_write(f, &
                           reshape( param(:,select_pft,output_param(n)),(/nrens,1,1,nselect_pft,1/) ), &
                           n+noutput_state+NCnauxvar, t, &
                           reshape(pft_pct_ave(select_pft),(/1,1,nselect_pft/)), &
                           reshape(hgt_pct_ave(select_hgt),(/1,1,nselect_hgt/)),&
                           (/1,1,1,1,1/), &
                           (/nrens,1,1,npft,1/))
                      if (f.eq.2) call netcdf_write(f, &
                           reshape( param_sim(:,select_pft,output_param(n)),(/1,1,1,nselect_pft,1/) ), &
                           n+noutput_state+NCnauxvar, t, &
                           reshape(pft_pct_ave(select_pft),(/1,1,nselect_pft/)), &
                           reshape(hgt_pct_ave(select_hgt),(/1,1,nselect_hgt/)),&
                           (/1,1,1,1,1/), &
                           (/1,1,1,npft,1/))
                   else
                      if (f.eq.1) call netcdf_write(f, &
                           reshape( param(:,select_pft,output_param(n)),(/nrens,1,1,nselect_pft,1/) ), &
                           n+noutput_state+NCnauxvar, t, &
                           reshape(pft_pct_ave,(/1,1,npft/)), &
                           reshape(hgt_pct_ave,(/1,1,nhgt/)),&
                           (/1,1,1,1,1/),&
                           (/nrens,1,1,npft,1/),&
                           select_pft=select_pft)
                      if (f.eq.2) call netcdf_write(f, &
                           reshape( param_sim(:,select_pft,output_param(n)),(/1,1,1,nselect_pft,1/) ), &
                           n+noutput_state+NCnauxvar, t, &
                           reshape(pft_pct_ave,(/1,1,npft/)), &
                           reshape(hgt_pct_ave,(/1,1,nhgt/)),&
                           (/1,1,1,1,1/), &
                           (/1,1,1,npft,1/), &
                           select_pft=select_pft)
                   endif
                enddo

                ! forcings (spatially distributed)
                do n=1,noutput_force 
                   if (averagehgt) then
                      ! average only selected hgt's
                      if (f.eq.1) call netcdf_write(f, &
                           reshape( force(:,:,:,:,output_force(n)),&
                           (/nrens,geographic%nx,geographic%ny,1,nselect_hgt/) ), &
                           n+noutput_param+noutput_state+NCnauxvar, t, &
                           pft_pct_out, &
                           hgt_pct_out(:,:,select_hgt), &
                           (/1,xcoord_offset,ycoord_offset,1,1/), &
                           (/nrens,geographic_out%nx,geographic_out%ny,1,nhgt/))
                      if (f.eq.2) call netcdf_write(f, &
                           reshape( force_sim(:,:,:,:,output_force(n)),&
                           (/1,geographic%nx,geographic%ny,1,nselect_hgt/) ), &
                           n+noutput_param+noutput_state+NCnauxvar, t, &
                           pft_pct_out, &
                           hgt_pct_out(:,:,select_hgt), &
                           (/1,xcoord_offset,ycoord_offset,1,1/), &
                           (/1,geographic_out%nx,geographic_out%ny,1,nhgt/))
                   else
                      ! output all hgt's
                      if (f.eq.1) call netcdf_write(f, &
                           reshape( force(:,:,:,:,output_force(n)),&
                           (/nrens,geographic%nx,geographic%ny,1,nselect_hgt/) ), &
                           n+noutput_param+noutput_state+NCnauxvar, t, &
                           pft_pct_out, &
                           hgt_pct_out, &
                           (/1,xcoord_offset,ycoord_offset,1,1/), &
                           (/nrens,geographic_out%nx,geographic_out%ny,1,nhgt/), &
                           select_hgt=select_hgt)
                      if (f.eq.2) call netcdf_write(f, &
                           reshape( force_sim(:,:,:,:,output_force(n)),&
                           (/1,geographic%nx,geographic%ny,1,nselect_hgt/) ), &
                           n+noutput_param+noutput_state+NCnauxvar, t, &
                           pft_pct_out, &
                           hgt_pct_out, &
                           (/1,xcoord_offset,ycoord_offset,1,1/), &
                           (/1,geographic_out%nx,geographic_out%ny,1,nhgt/), &
                           select_hgt=select_hgt)
                   endif
                enddo

                ! observations (gridded observations)
                do n=1,noutput_obser
                   if (f.eq.3) call netcdf_write(f, &
                        reshape( obser(:,:,:,output_obser(n)),(/2,geographic%nx,geographic%ny,1,1/) ), &
                        n+NCnauxvar, t, &
                        pft_pct_out(:,:,:), &
                        hgt_pct_out(:,:,:),&
                        (/1,xcoord_offset,ycoord_offset,1,1/), &
                        (/2,geographic_out%nx,geographic_out%ny,1,1/))
                enddo

             endif

          enddo

       endif

    endif

  end subroutine write_output

  subroutine close_output

    use driver_module
    use time_module
    use netcdf_module

    ! local variables
    integer :: f

    if (isprediction) then

       if (((timestep+timestep_start-1).eq.timestep_end).or. &
            (oldyear.and.(filelength.eq.3)).or.(oldmonth.and.(filelength.eq.2)).or. &
            (oldday.and.(filelength.eq.1))) then

          do f=1,3

             ! only save stuff if there are selected spatially distributed variables to save, except
             ! for the first process in distributed runs or during multiple region assimilation experiments.
             if (NChasfile(f)) then

                ! close NetCDF file
                if (verbose) write(*,'(A,I5,A,A)') 'Close output file ',mpiid,' : ',trim(NCfilenames(f))
                call netcdf_close(f)

             endif

          enddo

       endif

    endif

  end subroutine close_output

  subroutine exit_output

    use driver_module
    use netcdf_module

    implicit none

    ! destroy output NetCDF names + dimensions
    call netcdf_deallocate

  end subroutine exit_output

end module output_module





