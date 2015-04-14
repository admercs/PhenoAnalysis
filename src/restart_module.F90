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

module restart_module

  contains

    subroutine init_restart

      use driver_module
      use netcdf
      use typeSizes

      implicit none

      ! arguments

      ! local variables

    end subroutine init_restart
    
    subroutine open_restart

      use driver_module
      use netcdf
      use typeSizes

      implicit none

      ! arguments

      ! local variables

    end subroutine open_restart
    
    subroutine write_restart

      use driver_module
      use file_module

      use netcdf
      use typeSizes

      implicit none

      ! arguments

      ! local variables
      integer :: unit
      character(len=4) :: sid
      character(len=200) :: filename

      write(unit=sid,fmt='(i4.4)') mpiid

      unit = getunit()

      ! write restart file
      filename = trim(restartdir)//trim(geographic%name)//'.'//sid//'.restart.bin'
      if (verbose) write(*,'(A,A)') 'Writing model states/parameters/forcings to file : ',trim(filename)
      open(unit=unit, file=trim(filename), form='unformatted',action='write') 
      if (isprediction) then
         write(unit) param,state,force,cumobs
      else
         write(unit) param
      end if
      close(unit)

      if (simulator) then
         filename  = trim(restartdir)//trim(geographic%name)//'.'//sid//'.restart_sim.bin'
         if (verbose) write(*,'(A,A)') 'Writing simulator states/forcings to file : ',trim(filename)
         open(unit=unit, file=trim(filename), form='unformatted',action='write') 
         if (isprediction) then
            write(unit) param_sim,state_sim,force_sim
         else
            write(unit) param_sim
         end if
         close(unit)
      end if

    end subroutine write_restart
    
    subroutine read_restart

      use driver_module
      use file_module
      use exit_module

      use netcdf
      use typeSizes

      implicit none

      ! arguments

      ! local variables
      integer :: unit
      character(len=4) :: sid
      character(len=200) :: filename
      logical :: file_exists

      write(unit=sid,fmt='(i4.4)') mpiid

      ! read initialized state/parameters/forcings from previous model run
      filename = trim(restartdir)//trim(geographic%name)//'.'//sid//'.restart.bin'
      inquire(file=filename, exist=file_exists)
      if (.not.file_exists) then
         write(*,'(A,A)') 'Warning: restart file does not exist: ',trim(filename)
      else
         if (verbose) write(*,'(A,A)') 'Reading saved states/parameters/forcings from file : ',trim(filename)
         
         unit = getunit()
         open(unit=unit, file=trim(filename), form='unformatted',action='read') 
         if (isprediction) then
            read(unit) param,state,force,cumobs
         else
            read(unit) param
         end if         
         close(unit)
      end if

      if (simulator) then
         ! do the same for the simulator stream
         filename = trim(restartdir)//trim(geographic%name)//'.'//sid//'.restart_sim.bin'
         inquire(file=filename, exist=file_exists)
         if (.not.file_exists) then
            write(*,'(A,A)') 'Warning: simulator restart file does not exist: ',trim(filename)
         else
            if (verbose) write(*,'(A,A)') 'Reading saved simulator states/forcings from file : ',trim(filename)
            
            open(unit=unit, file=trim(filename), form='unformatted',action='read') 
            if (isprediction) then
               read(unit) param_sim,state_sim,force_sim
            else
               read(unit) param_sim
            endif
            close(unit)
         end if
      end if

    end subroutine read_restart

    subroutine close_restart

      use driver_module
      use netcdf
      use typeSizes

      implicit none

      ! arguments

      ! local variables

    end subroutine close_restart

    subroutine exit_restart

      use driver_module
      use netcdf
      use typeSizes

      implicit none

      ! arguments

      ! local variables

    end subroutine exit_restart

    subroutine write_parameters

      ! write mean parameters to ascii file

      use driver_module
      use file_module
      use pft_module

      implicit none

      ! arguments

      ! local variables
      integer :: p, n
      real(kind=4) :: ensemble_mean
      real(kind=4) :: ensemble_variance
      integer :: unit
      character(len=200) :: filename

      unit = getunit()

      if ( ((mpiid.eq.0).and.((nxprocs.ne.1).or.(nyprocs.ne.1))) &
           .or.(isprediction.and.(nxprocs.eq.1).and.(nyprocs.eq.1)) ) then

         if (regional) then
            filename = trim(outdir)//trim(geographic%name)//'.parameters_regional.dat'
            open(unit=unit, file=trim(filename), form='formatted',action='write')
            write(unit,'(I6,A,I4,A)') nparam,' Parameters and ',npft,' Plant Functional Types.'
         else
            filename = trim(outdir)//trim(geographic%name)//'.parameters_local.dat'
            open(unit=unit, file=trim(filename), form='formatted',action='write')
            write(unit,'(I6,A,I4,A)') nparam,' Parameters.'
         endif

         if (verbose) write(*,'(A,A)') 'Writing model parameters to file : ',trim(filename)

         do p=1,npft
            if (regional) then
               write(unit,'(A,A,F8.2,A,I14)') 'pft: ',pftnames(p), pft_pct_ave(p)*100., &
                    ' % of area. Total obs used: ',nint(cumobs(p),8)
            else
               write(unit,'(A,I14)') 'Local-scale parameters. Total obs used: ',nint(cumobs(p),8)
            endif
            do n=1,nparam
               call calc_mean(param(:,p,n),ensemble_mean)
               call calc_variance(param(:,p,n),ensemble_mean,ensemble_variance)
               write(unit,'(G12.5,A8,G12.5,A12,A20,A2,A10,A1)') ensemble_mean,' (mean) ',&
                    ensemble_variance,' (variance) ', &
                    trim(param_names(n)),' (', trim(param_units(n)),')'
            enddo
         enddo

         close(unit)

      endif

    end subroutine write_parameters

    subroutine read_parameters
      ! read mean and variance of parameters (ascii)

      use driver_module
      use file_module
      use exit_module
      use pft_module

      implicit none

      ! arguments

      ! local
      integer :: p, n
      integer :: unit
      character(len=1) :: cval
      real(kind=4) :: fmean, fvar
      character(len=4) :: sid
      character(len=200) :: filename
      logical :: file_exists

      logical, parameter :: nochill = .FALSE.

      ! analysis and NOT simulator mode: do nothing
      ! analysis and simulator mode: read all simulator parameters from simulator ascii file, 
      ! prediction mode: read mean/variance of pft-specific parameters from ascii file 
      if ((.not.analysis).or.simulator) then

         write(unit=sid,fmt='(i4.4)') mpiid

         unit = getunit()

         if (regional) then
            if (simulator) then
               filename = trim(namelistdir)//'simulator_parameters_regional.dat'
            else
               filename = trim(namelistdir)//'prediction_parameters_regional.dat'
            endif
         else
            if (simulator) then
               filename = trim(namelistdir)//'simulator_parameters_local.dat'
            else
               filename = trim(namelistdir)//'prediction_parameters_local.dat'
            endif
         endif

         inquire(file=filename, exist=file_exists)
         if (.not.file_exists) then
            write(*,'(A,A)') 'Parameter file does not exist: ',trim(filename)
            call model_exit(-1)
         end if

         if (verbose) write(*,'(A,A)') &
              'Reading fixed (non-assimilated) model parameters from file : ',trim(filename)
         open(unit=unit, file=trim(filename), form='formatted',action='read') 

         read(unit,*) cval
         do p=1,npft
            read(unit,*) cval
            do n=1,nparam

               if ((nochill).and.((n.eq.11).or.(n.eq.12).or.(n.eq.25))) then
                  fmean = 0.0
                  fvar = 0.0
               else
                  read(unit,'(G12.5,A8,G12.5)') fmean,cval,fvar
               end if

               ! prediction: prescribe parameter mean and variance
               if (.not.analysis) then
                  param(:,p,n) = fmean
                  param_var(p,n) = fvar
               endif
               ! simulator: prescribe parameter mean
               if (simulator) param_sim(:,p,n) = fmean
            enddo
         enddo

         close(unit)

      endif

    end subroutine read_parameters
 
  end module restart_module
