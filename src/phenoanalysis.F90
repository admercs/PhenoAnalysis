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

program phenoanalysis

  ! Main driver routine for the phenological data assimilation system which control both 
  ! the assimilation and the prognostic routines

  ! DEPENDENCIES: 
  use parameter_module
  use driver_module
  use init_module
  use exit_module
  use phenology_module
  use meteorology_module
  use output_module
  use assimilation_module
  use restart_module
  use time_module

#ifndef HIDE_MPI
    use mpi
#endif

  implicit NONE

  ! local variables

  ! cpu time profiling variables
  integer :: count0,count1,count_tot0,count_tot1,count_rate,count_max
  real(kind=4) :: time_tot
  real(kind=4) :: time_init
  real(kind=4) :: time_met
  real(kind=4) :: time_obs
  real(kind=4) :: time_predict
  real(kind=4) :: time_enkf
  real(kind=4) :: time_sync
  real(kind=4) :: time_write
  real(kind=4) :: time_exit
#ifndef HIDE_MPI
  real(kind=4) :: time_temp
#endif

  ! program arguments
  character(len=200) :: namelistfile
  character(len=200) :: sitefile
  character(len=200) :: sitelist
  character(len=8)   :: date_curr

  !---------------------------------------------
  ! INITIALIZATION
  !---------------------------------------------

  ! initialize wall clock variables
  time_tot = 0.
  time_init = 0.
  time_met = 0.
  time_obs = 0.
  time_predict = 0.
  time_enkf = 0.
  time_sync = 0.
  time_write = 0.
  time_exit = 0.

  call system_clock(count_tot0,count_rate,count_max)
  call system_clock(count0,count_rate,count_max)

    ! initialize number of processes and process id
#ifndef HIDE_MPI
    ! case: parallel processing using MPI
    call MPI_INIT(mpierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,mpiid,mpierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,mpierr)
    if (mpiid.eq.0) then
       write(*,'(A,I5,A)') 'MPI-support enabled with ',numprocs,' processes participating.'
       write(*,*)
    endif
#else
    ! case: serial processing
    write(*,'(A)') 'MPI-support not compiled. Force single process run.'
    write(*,*)
    mpiid = 0
    numprocs = 1
#endif

  ! get program arguments
  call getarg(1,namelistfile)
  call getarg(2,sitefile)
  call getarg(3,sitelist)
  call getarg(4,date_curr) ! optional for resubmission: current date YYYYMMDD, YYYYMM or YYYY

  call handle_error(len_trim(namelistfile) == 0, &
     "Please supply a namelist file name as the first argument")
  call handle_error(len_trim(sitefile) == 0, &
     "Please supply a site file name as the second argument")
  call handle_error(len_trim(sitelist) == 0, &
     "Please supply a site or a list of sites (e.g. 2..5) as the third argument")

  ! initialize random number generator with fixed value (if exact replicate is needed)
  call init_random

  ! initialize physical parameters
  call init_parameters

  ! read namelist and allocate + initialize driver arrays
  call init_driver(namelistfile,sitefile,sitelist,date_curr)

  ! read default parameter mean/variances
  call read_parameters
  
  ! perturb parameters for either analysis or prediction
  ! analysis: include parameter uncertainty for assimilation
  ! prediction: include parameter uncertainty in prediction
  if (nrens.gt.1) call perturb_parameters
  ! perturb states for analysis
  if (analysis) call perturb_states

  if (restart) then
     ! read restart file (states, parameters, forcings) 
     ! will overwrite default initial states if file is available
     call read_restart
  end if
    
  ! initialize observation stream
  if (analysis) call analysis_init

  ! initialize time
  call update_time(date_start, date_end, dtmod, .true.)

  ! allocate meteorological arrays
  call init_meteorology

  ! initialize output arrays
  call init_output

  ! bring all computation nodes to the same model day to avoid timeouts
  call synchronize_processes

  call system_clock(count1,count_rate,count_max)
  if (count0.gt.count1) count0 = count0 - count_max
  time_init = time_init + real(count1-count0,kind=4)/real(count_rate,kind=4)

  call system_clock(count_tot1,count_rate,count_max)
  if (count_tot0.gt.count_tot1) count_tot0 = count_tot0 - count_max
  time_tot = time_tot + real(count_tot1-count_tot0,kind=4)/real(count_rate,kind=4)

  !---------------------------------------------
  ! TIME STEPPING
  !---------------------------------------------

  ! run until defined ending period reached
  do while ((timestep+timestep_start-1).le.timestep_end)

     call system_clock(count_tot0,count_rate,count_max)

     ! open output files if needed
     call open_output

     !---------------------------------------------
     ! READ METEOROLOGICAL DATA
     !---------------------------------------------
     call system_clock(count0,count_rate,count_max)

     ! read meteo data and calculate daily averages
     call read_meteorology

     call system_clock(count1,count_rate,count_max)
     if (count0.gt.count1) count0 = count0 - count_max
     time_met = time_met + real(count1-count0,kind=4)/real(count_rate,kind=4)

     ! run the phenology model at dedicated time steps
     ! for instance: end of hour, end of day
     if (mod(timestep*dtmod,dtphen).eq.0) then

        if (mpiid.eq.0) write(*,'(A,I4,A,I3)') 'Date : ',year,' / ',doy

        !---------------------------------------------
        ! PREDICTION (insert your model here)
        !---------------------------------------------

        call system_clock(count0,count_rate,count_max)

        ! perturb met forcing at each integration time step
        if (analysis) call perturb_forcings

        ! predict phenological states
        call predict_phenology

        call system_clock(count1,count_rate,count_max)
        if (count0.gt.count1) count0 = count0 - count_max
        time_predict = time_predict + real(count1-count0,kind=4)/real(count_rate,kind=4)

        ! Analyze States and Parameters by use of EnKF
        if (analysis) then

           !---------------------------------------------
           ! READ SATELLITE DATA (MODIS L3 DATASETS)
           !---------------------------------------------

           call system_clock(count0,count_rate,count_max)

           ! read observations
           call analysis_observation

           call system_clock(count1,count_rate,count_max)
           if (count0.gt.count1) count0 = count0 - count_max
           time_obs = time_obs + real(count1-count0,kind=4)/real(count_rate,kind=4)

           !---------------------------------------------
           ! DATA ASSIMILATION (ANALYSIS)
           !---------------------------------------------

           call system_clock(count0,count_rate,count_max)
           
           ! global analysis of parameters
           call analysis_global

           ! local analysis of states
           call analysis_local

           call system_clock(count1,count_rate,count_max)
           if (count0.gt.count1) count0 = count0 - count_max
           time_enkf = time_enkf + real(count1-count0,kind=4)/real(count_rate,kind=4)

        endif  ! if we do analysis in this simulation

     endif ! if phenology model time step

     !---------------------------------------------
     ! OUTPUT
     !---------------------------------------------
     
     call system_clock(count0,count_rate,count_max)

     ! write states/parameters/forcings to output files
     call write_output

     ! close output files if needed
     call close_output

     call system_clock(count1,count_rate,count_max)
     if (count0.gt.count1) count0 = count0 - count_max
     time_write = time_write + real(count1-count0,kind=4)/real(count_rate,kind=4)

     !---------------------------------------------
     ! MPI SYNC
     !---------------------------------------------

     call system_clock(count0,count_rate,count_max)

     ! bring all computation nodes to the same model day to avoid timeouts
     call synchronize_processes

     call system_clock(count1,count_rate,count_max)
     if (count0.gt.count1) count0 = count0 - count_max
     time_sync = time_sync + real(count1-count0,kind=4)/real(count_rate,kind=4)

     call update_time(date_start, date_end, dtmod, .false.)
     initial = .false.

     call system_clock(count_tot1,count_rate,count_max)
     if (count_tot0.gt.count_tot1) count_tot0 = count_tot0 - count_max
     time_tot = time_tot + real(count_tot1-count_tot0,kind=4)/real(count_rate,kind=4)

  enddo   ! main loop

  call system_clock(count_tot0,count_rate,count_max)
  call system_clock(count0,count_rate,count_max)

  ! deallocate output arrays
  call exit_output

  ! write restart data at the end of model integration for prediction mode
  ! or in analysis mode if the restart date is not used
  if ((.not.analysis).or.(.not.restartdate)) then
     call write_restart
  endif

  ! clean up and exit
  call exit_meteorology

  if (analysis) call analysis_exit

  ! bring all computation nodes to the same model day have nicer final diagnostics
  call synchronize_processes

  call system_clock(count1,count_rate,count_max)
  if (count0.gt.count1) count0 = count0 - count_max
  time_exit = time_exit + real(count1-count0,kind=4)/real(count_rate,kind=4)

  call system_clock(count_tot1,count_rate,count_max)
  if (count_tot0.gt.count_tot1) count_tot0 = count_tot0 - count_max
  time_tot = time_tot + real(count_tot1-count_tot0,kind=4)/real(count_rate,kind=4)

#ifndef HIDE_MPI
  call MPI_ALLREDUCE(time_tot, time_temp, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, mpierr)
  time_tot = time_temp / real(numprocs,kind=4)
  call MPI_ALLREDUCE(time_init, time_temp, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, mpierr)
  time_init = time_temp / real(numprocs,kind=4)
  call MPI_ALLREDUCE(time_met, time_temp, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, mpierr)
  time_met = time_temp / real(numprocs,kind=4)
  call MPI_ALLREDUCE(time_obs, time_temp, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, mpierr)
  time_obs = time_temp / real(numprocs,kind=4)
  call MPI_ALLREDUCE(time_predict, time_temp, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, mpierr)
  time_predict = time_temp / real(numprocs,kind=4)
  call MPI_ALLREDUCE(time_enkf, time_temp, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, mpierr)
  time_enkf = time_temp / real(numprocs,kind=4)
  call MPI_ALLREDUCE(time_sync, time_temp, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, mpierr)
  time_sync = time_temp / real(numprocs,kind=4)
  call MPI_ALLREDUCE(time_write, time_temp, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, mpierr)
  time_write = time_temp / real(numprocs,kind=4)
  call MPI_ALLREDUCE(time_exit, time_temp, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, mpierr)
  time_exit = time_temp / real(numprocs,kind=4)
#endif

  if (mpiid.eq.0) then
     write(*,*)
     write(*,'(A,F12.2)') 'Total wall clock      [s]: ',time_tot
     write(*,'(A,F12.2)') 'Initialization        [s]: ',time_init
     write(*,'(A,F12.2)') 'Reading Met data      [s]: ',time_met
     write(*,'(A,F12.2)') 'Reading Observations  [s]: ',time_obs
     write(*,'(A,F12.2)') 'Predicting Phenology  [s]: ',time_predict
     write(*,'(A,F12.2)') 'EnKF Analysis         [s]: ',time_enkf
     write(*,'(A,F12.2)') 'MPI synchronizing     [s]: ',time_sync
     write(*,'(A,F12.2)') 'Writing output        [s]: ',time_write
     write(*,'(A,F12.2)') 'Exiting               [s]: ',time_exit
  endif

  call model_exit(0)

END PROGRAM phenoanalysis
