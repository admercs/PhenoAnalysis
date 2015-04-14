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

module assimilation_module

  ! This module contains the interface to fill the EnKF matrices with observation and states 
  ! and retrieve analyzed parameter and states from the updated matrices back into model space.

  ! counters
  integer, allocatable :: gnrobs(:)            ! global observation counter per cpu

  ! flags
  logical :: regridobs                         ! regrid observations before comparing to state vector
  logical :: useobspct                         ! model states represent observed pft/hgt pct distribution
  logical :: yearlyanalysis                    ! perform yearly analysis or instantaneous analysis
  logical :: globalanalysis                    ! perform global analysis of parameter (and states)
  logical :: localanalysis                     ! perform local analysis of states (and NOT parameters)
  logical :: weightmag                         ! normalize HA-D and observation errors to observation magnitude
  logical :: weightdist                        ! normalize HA-D and observation errors to observation distance
  logical :: weightarea                        ! normalize HA-D and observation errors to observation area
  logical :: weightgrid                        ! weight gridding by observation uncertainty
  logical :: sepminmax                         ! mean of MIN-parameter always lower than of mean MAX-parameter

  ! parameters
  integer :: nrobsmax                          ! upper limit for # observations
  integer :: nrobsmin                          ! lower limit for # observations
  real(kind=4) :: deflation                    ! deflation factor 0.0-1.0
  real(kind=4) :: inflation                    ! inflation factor 0.0-1.0
  real(kind=4) :: influence                    ! HW of influence distance relative to grid point extent

  ! directories
  character(len=100) :: satellitetype          ! name of assimilation data
  character(len=100) :: satellitedir           ! directory of assimilation data
  character(len=100) :: landcovertype          ! name of landcover data
  character(len=100) :: landcoverdir           ! directory of landcover data

  type :: obs_type
     ! names/units
     character(len=100) :: name                ! name of observation
     real(kind=4) :: magnitude                 ! magnitude (units of observations)
     ! geographic extent
     character(len=100) :: type                ! observation geographic projection type
     real(kind=4)       :: xmin                ! minimum x coordinate bound (edge) of observation grid
     real(kind=4)       :: xmax                ! maximum x coordinate bound (edge) of observation grid
     real(kind=4)       :: ymin                ! minimum y coordinate bound (edge) of observation grid
     real(kind=4)       :: ymax                ! maximum y coordinate bound (edge) of observation grid
     real(kind=4)       :: dx                  ! coordinate x distance of observation grid
     real(kind=4)       :: dy                  ! coordinate y distance of observation grid
     integer            :: nx                  ! number of x coordinate points of observation grid
     integer            :: ny                  ! number of y coordinate points of observation grid
     character(len=100) :: optargs             ! optional arguments for projection type
     ! data and error
     real(kind=4), allocatable :: dat(:,:)       ! observation nx,ny (units of observation)
     real(kind=4), allocatable :: err(:,:)       ! observation error nx,ny  (units of observation)
     ! absolute 4D location (relative to global geographic coordinates)
     real(kind=4), allocatable :: x(:,:)         ! x coordinate of observation
     real(kind=4), allocatable :: y(:,:)         ! y coordinate of observation
     real(kind=4), allocatable :: hgt(:,:)       ! elevation nx,ny (meters a.s.l.)
     real(kind=4), allocatable :: pft(:,:)       ! predominant pft nx,ny (% cover)
     real(kind=4), allocatable :: lon(:,:)       ! longitude nx,ny (degrees east)
     real(kind=4), allocatable :: lat(:,:)       ! latitude nx,ny (degrees north)
     ! integer 4D location (relative to local model grid) (we might not need them!)
     integer, allocatable :: xid(:,:)          ! index of horizontal grid point nx,ny (1..nlon)
     integer, allocatable :: yid(:,:)          ! index of vertical grid point nx,ny (1..nlat)
     integer, allocatable :: hgtid(:,:)        ! index of elevation class nx,ny (1..nhgt)
     integer, allocatable :: pftid(:,:)        ! index of pft class nx,ny (1..npft)
     ! geodesic distance of observation to model grid points
     real(kind=4), allocatable :: dist(:,:,:,:)  ! nx x ny x nlon x nlat (km)
     ! area covered by each observation point
     real(kind=4), allocatable :: area(:,:)      ! nx,ny (km^2)
     ! observation subgrid parameters (elevation and pft classes)
     real(kind=4), allocatable :: hgt_sum(:,:)   ! sum of selected elevation classes nx x ny (%)
     real(kind=4), allocatable :: pft_sum(:,:)   ! sum of selected pft classes nx x ny (%)
     real(kind=4), allocatable :: hgt_pct(:,:,:) ! elevation classes covered by observation nx x ny x nhgt (%)
     real(kind=4), allocatable :: pft_pct(:,:,:) ! pft classes covered by observation nx x ny x npft (%)
     ! size
     integer :: dt                             ! observation time step (seconds)
     integer :: t0                             ! first observation time step in year (seconds)
     integer :: nt                             ! maximum number of observation time steps (per year)
     integer :: nv                             ! number of valid observations
  end type obs_type

  type(obs_type), allocatable :: obs(:)

  type :: ana_type
     real(kind=4), allocatable :: obsdat(:,:)  ! observation nx x ny x nt (units of observation)
     real(kind=4), allocatable :: obserr(:,:)  ! observation error nx x ny x nt (units of observation)
     real(kind=4), allocatable :: HA(:,:,:)    ! model states in observation space np x nt x nrens (units of observation)
     integer, allocatable :: xid(:)          ! index of W-E point np (1..nlon)
     integer, allocatable :: yid(:)          ! index of N-S grid point np (1..nlat)
     integer :: np
     integer :: nt
  end type ana_type

  type(ana_type), allocatable :: ana(:)

  ! global analysis matrix 
  real(kind=8), allocatable :: Xglobal(:,:)

  ! random state matrix
  integer :: nrand
  real(kind=4), allocatable :: Arand(:,:)

contains

  subroutine perturb_parameters

    ! Add gaussian noise with initial parameter variance to parameters 

#ifndef HIDE_MPI
    use mpi    
#endif

    use driver_module
    use enkf_module_r8

    IMPLICIT NONE

    ! local variables
    integer :: n, p

    ! perturb parameters on all processes in order to keep random number generation synchronized
    do n=1,nparam
       do p=1,npft
          call ensemble_perturb(param(:,p,n),param_var(p,n))
       end do
    end do

#ifndef HIDE_MPI
    ! then redistribute single parameter set among all processes
    do n=1,nparam          
       do p = 1, npft
          call MPI_BCAST(param(:,p,n),nrens,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
       end do ! pfts
    end do ! parameters
#endif

  end subroutine perturb_parameters

  subroutine perturb_states

    ! Add gaussian noise with initial state variance to states

    use driver_module
    use enkf_module_r8

    IMPLICIT NONE

    ! local variables
    integer :: n, h, p, i, j

    if (isprediction) then

       do n=1,nselect_state
          do h=1,nselect_hgt
             do p=1,nselect_pft
                do j=1,geographic%ny
                   do i=1,geographic%nx
                      call ensemble_perturb(state(:,i,j,p,h,select_state(n)),state_var(select_state(n)))
                   end do
                end do
             end do
          end do
       end do

    endif

  end subroutine perturb_states

  subroutine perturb_forcings

    ! Add gaussian noise with initial forcing variance to forcings

    use driver_module
    use enkf_module_r8

    IMPLICIT NONE

    ! local variables
    integer :: n, h, i, j

    if (isprediction) then

       do n=1,nselect_force
          !       do n=1,nforce
          do h=1,nselect_hgt
             do j=1,geographic%ny
                do i=1,geographic%nx
                   call ensemble_perturb(force(:,i,j,h,select_force(n)),force_var(select_force(n)))
                   !                   call ensemble_perturb(force(:,i,j,h,n),force_var(n))
                end do
             end do
          end do
       end do

    endif

  end subroutine perturb_forcings

  subroutine analysis_init

    ! this subroutine initializes the observation streams

    use driver_module
    use time_module
    use modis_module
    use enkf_module_r8
    use topo_module
    use pft_module
    use file_module
    use exit_module
    use proj_module
#ifndef HIDE_MPI
    use mpi
#endif

    implicit none

    ! local variables
    integer :: i,j,k,l,q,s
    integer :: yyyy, mm, dd, hh, min, sec
    integer :: unit
    character(len=200) :: filename
    integer :: ar

    ! allocate and initialize random matrix
    nrand = 10
    allocate(Arand(nrand,nrens))
    Arand(:,:) = 0.
    do ar=1,nrand
       call ensemble_perturb(Arand(ar,:),1.)
    end do

    ! allocate global observation counter
    if (mpiid.eq.0) then
       allocate(gnrobs(numprocs))
       gnrobs(:) = 0
    endif

    ! allocate global analysis matrix
    if (globalanalysis) then
       allocate(Xglobal(nrens,nrens))
    endif

    if (isprediction) then

       ! initialize observation stream
       allocate(obs(nselect_obser))
       obs%name = obser_names(select_obser)
       obs%nt = 0
       obs%dt = -1

       ! allocate observation arrays
       select case (trim(satellitetype))
       case ("modis")
          call init_modis(satellitedir,obs%name,geographic%type,geographic%xmin,geographic%xmax, &
               geographic%ymin,geographic%ymax,geographic%dx,geographic%dy,geographic%optargs,nodata, &
               obs%type,obs%xmin,obs%xmax,obs%ymin,obs%ymax,obs%dx,obs%dy, obs%optargs, &
               obs%nx, obs%ny, obs%nt, obs%t0, obs%dt)
       case default
          write(*,'(A,A)') "Assimilation Data Type not implemented: ",trim(satellitetype)
          call model_exit(-1)
       end select

       do s=1,nselect_obser

          if (.not.yearlyanalysis) then
             obs(s)%nt = 1
          endif

          if ((obs(s)%nx.ne.0).and.(obs(s)%nt.ne.0)) then
             allocate(obs(s)%dat(obs(s)%nx,obs(s)%ny))
             allocate(obs(s)%err(obs(s)%nx,obs(s)%ny))

             allocate(obs(s)%x(obs(s)%nx,obs(s)%ny))
             allocate(obs(s)%y(obs(s)%nx,obs(s)%ny))
             allocate(obs(s)%lon(obs(s)%nx,obs(s)%ny))
             allocate(obs(s)%lat(obs(s)%nx,obs(s)%ny))
             allocate(obs(s)%hgt(obs(s)%nx,obs(s)%ny))
             allocate(obs(s)%pft(obs(s)%nx,obs(s)%ny))

             allocate(obs(s)%dist(obs(s)%nx,obs(s)%ny,geographic%nx,geographic%ny))
             allocate(obs(s)%area(obs(s)%nx,obs(s)%ny))

             allocate(obs(s)%xid(obs(s)%nx,obs(s)%ny))
             allocate(obs(s)%yid(obs(s)%nx,obs(s)%ny))
             allocate(obs(s)%hgtid(obs(s)%nx,obs(s)%ny))
             allocate(obs(s)%pftid(obs(s)%nx,obs(s)%ny))

             allocate(obs(s)%hgt_sum(obs(s)%nx,obs(s)%ny))
             allocate(obs(s)%pft_sum(obs(s)%nx,obs(s)%ny))
             allocate(obs(s)%hgt_pct(obs(s)%nx,obs(s)%ny,nhgt))
             allocate(obs(s)%pft_pct(obs(s)%nx,obs(s)%ny,npft))

          endif

          ! 3. precalculate geographic location and transformation grid
          do k=1,obs(s)%nx
             obs(s)%x(k,:) = obs(s)%xmin + real(k-1,kind=4) + 0.5*obs(s)%dx
          end do
          do l=1,obs(s)%ny
             obs(s)%y(:,l) = obs(s)%ymin + real(l-1,kind=4) + 0.5*obs(s)%dy
          end do

          call proj_inv(obs(s)%type,obs(s)%x,obs(s)%y,obs(s)%lon,obs(s)%lat, obs(s)%optargs,nodata)

          if (regional) then

             ! 4. get observation elevation and elevation distribution
             call read_topo(obs(s)%type,obs(s)%xmin,obs(s)%xmax,obs(s)%ymin,obs(s)%ymax, &
                  obs(s)%dx,obs(s)%dy,obs(s)%optargs,nodata,topo=obs(s)%hgt)    
             call read_topo(obs(s)%type,obs(s)%xmin,obs(s)%xmax,obs(s)%ymin,obs(s)%ymax, &
                  obs(s)%dx,obs(s)%dy,obs(s)%optargs,nodata,topohist=obs(s)%hgt_pct,topolevels=hgtlevels) 
             
             ! generate % out of histogram counts
             do l=1,obs(s)%ny
                do k=1,obs(s)%nx
                   obs(s)%hgt_pct(k,l,:) = obs(s)%hgt_pct(k,l,:) / max(sum(obs(s)%hgt_pct(k,l,:)),1.0)
                end do
             end do
             obs(s)%hgt_sum = sum(obs(s)%hgt_pct(:,:,select_hgt),3)
             obs(s)%hgtid = maxloc(obs(s)%hgt_pct,3)

             ! 5. get observation pft distribution
             call read_pft(obs(s)%type,obs(s)%xmin,obs(s)%xmax,obs(s)%ymin,obs(s)%ymax, &
                  obs(s)%dx,obs(s)%dy,obs(s)%optargs,nodata,pft3d=obs(s)%pft_pct)    

             ! remove observations with too high % water
             do l=1,obs(s)%ny
                do k=1,obs(s)%nx
                   if (obs(s)%pft_pct(k,l,npft).gt.0.01) then
                      obs(s)%pft_pct(k,l,1:npft-1) = 0.0
                      obs(s)%pft_pct(k,l,npft) = 1.0
                   end if
                end do
             end do

             obs(s)%pft_sum = sum(obs(s)%pft_pct(:,:,select_pft),3)
             obs(s)%pftid = maxloc(obs(s)%pft_pct,3)
             obs(s)%pft = maxval(obs(s)%pft_pct,3)

          else
             ! local-scale analysis: no elevation and pft distribution used
             obs(s)%hgt = nodata
             obs(s)%hgt_pct = 1.0
             obs(s)%hgt_sum = 1.0
             obs(s)%pft = nodata
             obs(s)%pft_pct = 1.0
             obs(s)%pft_sum = 1.0
          endif

          ! calculate distance from observation locations to each model grid point (units: km)
          do j=1,geographic%ny
             do i=1,geographic%nx
                do l=1,obs(s)%ny
                   do k=1,obs(s)%nx
                      obs(s)%dist(k,l,i,j) = calc_distance(obs(s)%lon(k,l),obs(s)%lat(k,l),lon(i,j),lat(i,j))
                   end do
                end do
             end do
          end do

       end do

       ! initialize analysis stream
       allocate(ana(nselect_obser))
       ana%np = 0
       ana%nt = obs%nt

       ! count maximum number of model-observation couples
       do s=1,nselect_obser
          if (regridobs) then
             ana(s)%np = geographic%nx * geographic%ny
          else
             do j=1,geographic%ny
                do i=1,geographic%nx
                   do l=1,obs(s)%ny
                      do k=1,obs(s)%nx
                         if ((obs(s)%dist(k,l,i,j).le.(0.5*sqrt(area(i,j)))).and. &
                              (obs(s)%pft_sum(k,l).ge.minsumpft).and. &
                              (obs(s)%hgt_sum(k,l).ge.minsumhgt)) then
                            ana(s)%np = ana(s)%np + 1
                         endif
                      end do
                   end do
                end do
             end do
          endif
       end do

       ! allocate analysis stream
       do s=1,nselect_obser
          if ((ana(s)%np.ne.0).and.(ana(s)%nt.ne.0)) then
             allocate(ana(s)%obsdat(ana(s)%np,ana(s)%nt))
             allocate(ana(s)%obserr(ana(s)%np,ana(s)%nt))
             allocate(ana(s)%HA(ana(s)%np,ana(s)%nt,nrens))
             allocate(ana(s)%xid(ana(s)%np))
             allocate(ana(s)%yid(ana(s)%np))
          endif
       end do

       ! assign static grid point id's to analysis stream
       do s=1,nselect_obser
          ana(s)%obsdat = nodata
          ana(s)%obserr = nodata
          ana(s)%HA = nodata

          q=0
          do j=1,geographic%ny
             do i=1,geographic%nx
                if (regridobs) then
                   q = q + 1
                   ana(s)%xid(q) = i
                   ana(s)%yid(q) = j
                else
                   do l=1,obs(s)%ny
                      do k=1,obs(s)%nx
                         if ((obs(s)%dist(k,l,i,j).le.(0.5*sqrt(area(i,j)))).and. &
                              (obs(s)%pft_sum(k,l).ge.minsumpft).and. &
                              (obs(s)%hgt_sum(k,l).ge.minsumhgt)) then
                            q = q + 1
                            ana(s)%xid(q) = i
                            ana(s)%yid(q) = j
                         endif
                      end do
                   end do
                endif
             end do
          end do
       end do

       if (verbose) then
          write(*,*)
          write(*,'(A,I8)') 'Process ID:                                                   ',mpiid
          write(*,'(A,I8)') 'Total Number of potential raw observations in Analysis :      ',sum(obs%nx*obs%ny)
          write(*,'(A,I8)') 'Total Number of potential gridded observations in Analysis :  ',sum(ana%np)
          write(*,'(A,I8)') 'Total Number of observation integration timesteps :           ',maxval(obs%nt)
          write(*,*)
       end if
    else
       if (verbose) then
          write(*,*)
          write(*,'(A,I8)') 'Process ID:                                                   ',mpiid
          write(*,'(A)') 'No Prediction / Observation / Analysis '
          write(*,*)
       end if
    endif

    ! read analysis date
    if ((restart).and.(restartdate)) then
       if (mpiid.eq.0) then
          unit = getunit()
          filename  = trim(outdir)//trim(geographic%name)//'.analysis-date.dat'
          if (verbose) write(*,'(A,A)') 'Read analysis date from file : ',trim(filename)
          open(unit=unit, file=trim(filename), form='formatted',action='read')
          read(unit,'(I4,1X,I2,1X,I2)') yyyy,mm,dd
          close(unit)

          ! set start date to last date of analysis plus one day
          call calc_date(calc_timestep(yyyy,mm,dd,0,0,0,dtmod) + 86400/dtmod,yyyy,mm,dd,hh,min,sec,dtmod)
          write(date_start, *) yyyy*10000 + mm*100 + dd
       endif

#ifndef HIDE_MPI
       ! distribute namelist entries among processes
       call MPI_BCAST(date_start,20,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
#endif

    endif

  end subroutine analysis_init

  subroutine analysis_exit

    use driver_module
    use modis_module

    IMPLICIT NONE

    ! local variables
    integer :: s

    if (isprediction) then

       ! deallocate observation arrays
       do s=1,nselect_obser
          if ((obs(s)%nx.ne.0).and.(obs(s)%ny.ne.0).and.(obs(s)%nt.ne.0)) then
             deallocate(obs(s)%dat)
             deallocate(obs(s)%err)

             deallocate(obs(s)%x)
             deallocate(obs(s)%y)
             deallocate(obs(s)%lon)
             deallocate(obs(s)%lat)
             deallocate(obs(s)%hgt)
             deallocate(obs(s)%pft)

             deallocate(obs(s)%dist)
             deallocate(obs(s)%area)

             deallocate(obs(s)%xid)
             deallocate(obs(s)%yid)
             deallocate(obs(s)%hgtid)
             deallocate(obs(s)%pftid)

             deallocate(obs(s)%hgt_sum)
             deallocate(obs(s)%pft_sum)
             deallocate(obs(s)%hgt_pct)
             deallocate(obs(s)%pft_pct)

          endif
       end do

       deallocate(obs)

       ! deallocate analysis stream
       do s=1,nselect_obser
          if ((ana(s)%np.gt.0).and.(ana(s)%nt.gt.0)) then
             deallocate(ana(s)%obsdat)
             deallocate(ana(s)%obserr)
             deallocate(ana(s)%HA)
             deallocate(ana(s)%xid)
             deallocate(ana(s)%yid)
          endif
       end do
       deallocate(ana)

       ! dealocate global analysis matrix
       if (globalanalysis) then
          deallocate(Xglobal)
       endif

       ! deallocate random matrix
       deallocate(Arand)

       ! deallocate MODIS tile arrays
       call exit_modis

    endif

  end subroutine analysis_exit

  subroutine analysis_observation

    ! this subroutine reads the observation streams
    ! and derives the number of observations needed for assimilation 
    ! this subroutine is very specific to the chosen observation stream

    use driver_module
    use time_module
    use modis_module
    use enkf_module_r8

#ifndef HIDE_MPI
    use mpi
#endif

    implicit none

    ! local variables
    integer :: nrobs                      ! number of observations for this process
    real(kind=4) :: totobs(npft)          ! number of observations for this process weighted by pft
    integer :: i, j, h, k, l, n, p, q, r, s, t
    logical :: hasdata

    real(kind=4), allocatable :: tempmean(:,:), tempwgts(:,:)
    real(kind=4) :: w

#ifndef HIDE_MPI
    integer :: proc
    integer :: mpistatus(MPI_STATUS_SIZE)
#endif

    ! gridding to save observation stream
    type :: obsgrid_type
       real(kind=4), allocatable :: dat(:,:)
       real(kind=4), allocatable :: err(:,:)
       real(kind=4), allocatable :: totdat(:,:)
       real(kind=4), allocatable :: toterr(:,:)
       real(kind=4), allocatable :: pft_sum(:,:)
       real(kind=4), allocatable :: pft_pct(:,:,:)
       real(kind=4), allocatable :: hgt_sum(:,:)
       real(kind=4), allocatable :: hgt_pct(:,:,:)
       integer, allocatable      :: np(:,:)
       integer, allocatable      :: npmax(:,:)
    end type obsgrid_type

    type(obsgrid_type) :: obsgrid

    ! reset global observation counters and observation / analysis data
    if ((.not.yearlyanalysis).or.(doy.eq.1)) then
       if (mpiid.eq.0) then
          gnrobs(:) = 0
       endif
       if (isprediction) then
          do s=1,nselect_obser
             ana(s)%obsdat(:,:) = nodata
             ana(s)%obserr(:,:) = nodata
             ana(s)%HA(:,:,:) = nodata
          end do
       endif
    endif

    ! reset local observation counters
    nrobs = 0
    totobs(:) = 0.

    if (isprediction) then

       allocate(tempmean(nselect_pft, nselect_hgt))
       allocate(tempwgts(nselect_pft, nselect_hgt))

       if (saveobs.or.regridobs) then
          allocate(obsgrid%dat(geographic%nx,geographic%ny))
          allocate(obsgrid%err(geographic%nx,geographic%ny))
          allocate(obsgrid%totdat(geographic%nx,geographic%ny))
          allocate(obsgrid%toterr(geographic%nx,geographic%ny))
          allocate(obsgrid%pft_sum(geographic%nx,geographic%ny))
          allocate(obsgrid%pft_pct(geographic%nx,geographic%ny,npft))
          allocate(obsgrid%hgt_sum(geographic%nx,geographic%ny))
          allocate(obsgrid%hgt_pct(geographic%nx,geographic%ny,nhgt))
          allocate(obsgrid%np(geographic%nx,geographic%ny))
          allocate(obsgrid%npmax(geographic%nx,geographic%ny))
       endif

       if (saveobs) then
          obser(:,:,:,:) = nodata
       endif

       do s=1,nselect_obser

          if (mod(doy*86400-obs(s)%t0,obs(s)%dt).eq.0) then

             ! reset observations
             obs(s)%dat = nodata
             obs(s)%err = nodata

             if (yearlyanalysis) then
                t = (doy*86400-obs(s)%t0)/obs(s)%dt + 1
             else
                t = 1
             endif

             if (simulator) then

                ! assign simulator states to observation stream
                hasdata = .false.

                do n=1,nstate

                   if (trim(state_names(n)) == trim(obs(s)%name)) then

                      do l=1,obs(s)%ny

                         do k=1,obs(s)%nx

                            do j=1,geographic%ny

                               do i=1,geographic%nx

                                  if (obs(s)%dist(k,l,i,j).le.(0.5*sqrt(area(i,j)))) then

                                     if ((obs(s)%pft_sum(k,l).ge.minsumpft).and.(obs(s)%hgt_sum(k,l).ge.minsumhgt)) then

                                        ! weighted summation of simulated obs by %pft & %hgt
                                        do h=1,nselect_hgt
                                           do p=1,nselect_pft
                                              tempmean(p,h) = state_sim(1,i,j,p,h,n)
                                              tempwgts(p,h) = obs(s)%pft_pct(k,l,select_pft(p)) * &
                                                   obs(s)%hgt_pct(k,l,select_hgt(h)) / &
                                                   (obs(s)%pft_sum(k,l) * obs(s)%hgt_sum(k,l))
                                           end do
                                        end do

                                        obs(s)%dat(k,l) = sum(tempmean*tempwgts)
                                        obs(s)%err(k,l) = obser_var(select_obser(s))

                                        hasdata = .true.

                                     endif

                                  end if

                               end do

                            end do

                         end do

                      end do

                   endif

                end do

             else
                ! read MODIS stream

                call read_modis(obs(s)%name,obs(s)%dat,obs(s)%err,year,doy,nodata,hasdata)

             endif

             if (hasdata) then

                ! regrid observations to output grid
                if (saveobs.or.regridobs) then

                   ! reset counters and sums
                   obsgrid%dat(:,:) = 0.
                   obsgrid%err(:,:) = 0.
                   obsgrid%totdat(:,:) = 0.
                   obsgrid%toterr(:,:) = 0.
                   obsgrid%pft_sum(:,:) = 0.
                   obsgrid%pft_pct(:,:,:) = 0.
                   obsgrid%hgt_sum(:,:) = 0.
                   obsgrid%hgt_pct(:,:,:) = 0.
                   obsgrid%np(:,:) = 0
                   obsgrid%npmax(:,:) = 0

                   do l=1,obs(s)%ny

                      do k=1,obs(s)%nx

                         do j=1,geographic%ny

                            do i=1,geographic%nx

                               if (obs(s)%dist(k,l,i,j).le.(0.5*sqrt(area(i,j)))) then

                                  if ((obs(s)%dat(k,l).ne.nodata).and.(obs(s)%err(k,l).ne.nodata).and. &
                                       (obs(s)%pft_sum(k,l).ge.minsumpft).and.(obs(s)%hgt_sum(k,l).ge.minsumhgt)) then

                                     w = 1.

                                     if (weightdist) w = w / max((obs(s)%dist(k,l,i,j) / sqrt(area(i,j))),1.)
                                     if (weightarea) w = w / max(sqrt(area(i,j) / obs(s)%area(k,l)),1.)
                                     if (weightgrid) w = w / obs(s)%err(k,l)

                                     ! weighted addition: 
                                     ! x_bar = sum(w*x)/sum(w)
                                     ! e_bar = sum(w^2*e)/sum(w^2)
                                     obsgrid%dat(i,j) = obsgrid%dat(i,j) + obs(s)%dat(k,l) * w
                                     obsgrid%err(i,j) = obsgrid%err(i,j) + obs(s)%err(k,l) * w !**2 
                                     obsgrid%totdat(i,j) = obsgrid%totdat(i,j) + w
                                     obsgrid%toterr(i,j) = obsgrid%toterr(i,j) + w !**2
                                     obsgrid%pft_sum(i,j) = obsgrid%pft_sum(i,j) + w * obs(s)%pft_sum(k,l)                         
                                     obsgrid%pft_pct(i,j,:) = obsgrid%pft_pct(i,j,:) + w * obs(s)%pft_pct(k,l,:) 
                                     obsgrid%hgt_sum(i,j) = obsgrid%hgt_sum(i,j) + w * obs(s)%hgt_sum(k,l)                         
                                     obsgrid%hgt_pct(i,j,:) = obsgrid%hgt_pct(i,j,:) + w * obs(s)%hgt_pct(k,l,:)

                                     obsgrid%np(i,j) = obsgrid%np(i,j) + 1
                                  endif

                                  if ((obs(s)%pft_sum(k,l).ge.minsumpft).and.(obs(s)%hgt_sum(k,l).ge.minsumhgt)) then
                                     obsgrid%npmax(i,j) = obsgrid%npmax(i,j) + 1
                                  endif

                               end if

                            end do

                         end do

                      end do

                   end do

                   ! compute gridded mean and error values
                   do j=1,geographic%ny
                      do i=1,geographic%nx

                         if (obsgrid%np(i,j).gt.int(0.1*obsgrid%npmax(i,j))) then
                            obsgrid%dat(i,j) = obsgrid%dat(i,j)/obsgrid%totdat(i,j)
                            obsgrid%err(i,j) = obsgrid%err(i,j)/obsgrid%toterr(i,j)
                            obsgrid%pft_sum(i,j) = obsgrid%pft_sum(i,j)/obsgrid%totdat(i,j)
                            obsgrid%pft_pct(i,j,:) = obsgrid%pft_pct(i,j,:)/obsgrid%totdat(i,j)
                            obsgrid%hgt_sum(i,j) = obsgrid%hgt_sum(i,j)/obsgrid%totdat(i,j)
                            obsgrid%hgt_pct(i,j,:) = obsgrid%hgt_pct(i,j,:)/obsgrid%totdat(i,j)
                         else
                            obsgrid%dat(i,j) = nodata
                            obsgrid%err(i,j) = nodata
                            obsgrid%pft_sum(i,j) = 1.0
                            obsgrid%pft_pct(i,j,:) = 0.0
                            obsgrid%hgt_sum(i,j) = 1.0
                            obsgrid%hgt_pct(i,j,:) = 0.0
                         endif

                      end do ! lon
                   end do ! lat

                end if

                ! fill output observation grid
                if (saveobs) then
                   do n=1,nobser
                      if (trim(obser_names(n)) == trim(obs(s)%name)) then
                         do j=1,geographic%ny
                            do i=1,geographic%nx
                               obser(1,i,j,n) = obsgrid%dat(i,j)
                               obser(2,i,j,n) = obsgrid%err(i,j)
                            end do ! lon
                         end do ! lat  
                      endif
                   end do
                endif

                ! fill observations and corresponding model states into analysis stream
                ! this is where the application-specific model-observation operator is formulated
                ! also screen observations by pft and hgt sum

                ! for FPAR and LAI the model states need to be aggregated by grid point in order to match
                ! the subgrid-scale hgt and pft distribution of each observation

                do n=1,nstate

                   if (trim(state_names(n)) == trim(obs(s)%name)) then

                      if (regridobs) then
                         ! analyze gridded observations

                         q = 0

                         do j = 1,geographic%ny

                            do i = 1,geographic%nx

                               ! aggregation of subgrid-scale gridded model states to gridded observations

                               q = q + 1

                               if ((obsgrid%dat(i,j).ne.nodata).and.(obsgrid%err(i,j).ne.nodata).and. &
                                    (obsgrid%pft_sum(i,j).ge.minsumpft).and. &
                                    (obsgrid%hgt_sum(i,j).ge.minsumhgt)) then

                                  ! assign gridded observations
                                  ana(s)%obsdat(q,t) = obsgrid%dat(i,j)
                                  ana(s)%obserr(q,t) = obsgrid%err(i,j)

                                  ! aggregate gridded model states

                                  ! weighted linear mixing of each observation with state vector by %pft and %hgt
                                  ! average ensemble state vector by pft/hgt.
                                  do h=1,nselect_hgt
                                     do p=1,nselect_pft
                                        tempwgts(p,h) = obsgrid%pft_pct(i,j,select_pft(p)) * &
                                             obsgrid%hgt_pct(i,j,select_hgt(h)) / &
                                             (obsgrid%pft_sum(i,j) * obsgrid%hgt_sum(i,j))
                                     enddo
                                  enddo

                                  do r=1,nrens
                                     ana(s)%HA(q,t,r) = sum(state(r,i,j,:,:,n) * tempwgts)
                                  enddo

                                  ! update observation counter
                                  nrobs = nrobs + 1

                                  ! update total (weighted observation counter by pft and hgt)
                                  ! do not count observed pft's which are not assimilated
                                  totobs(select_pft) = totobs(select_pft) + obsgrid%pft_pct(i,j,select_pft) / obsgrid%pft_sum(i,j)

                                  ! match output pft/hgt coverage to observation pft/hgt coverage
                                  ! for validation purposes (only during assimilation)
                                  if (useobspct) then
                                     pft_pct_out(i,j,:) =  obsgrid%pft_pct(i,j,:)
                                     hgt_pct_out(i,j,:) =  obsgrid%hgt_pct(i,j,:)
                                  endif

                               endif ! valid observation

                            end do ! lon
                         end do ! lat

                      else
                         ! analyze raw observations

                         q = 0

                         do j=1,geographic%ny

                            do i=1,geographic%nx

                               do l=1,obs(s)%ny

                                  do k=1,obs(s)%nx

                                     if ((obs(s)%dist(k,l,i,j).le.(0.5*sqrt(area(i,j)))).and. &
                                          (obs(s)%pft_sum(k,l).ge.minsumpft).and. &
                                          (obs(s)%hgt_sum(k,l).ge.minsumhgt)) then

                                        q = q + 1

                                        if ((obs(s)%dat(k,l).ne.nodata).and.(obs(s)%err(k,l).ne.nodata)) then

                                           ! assign raw observations
                                           ana(s)%obsdat(q,t) = obs(s)%dat(k,l)
                                           ana(s)%obserr(q,t) = obs(s)%err(k,l)

                                           ! aggregate subgrid-scale model states

                                           ! weighted linear mixing of each observation with state vector by %pft and %hgt
                                           ! average ensemble state vector by pft/hgt.
                                           do h=1,nselect_hgt
                                              do p=1,nselect_pft
                                                 tempwgts(p,h) = obs(s)%pft_pct(k,l,select_pft(p)) * &
                                                      obs(s)%hgt_pct(k,l,select_hgt(h)) / &
                                                      (obs(s)%pft_sum(k,l) * obs(s)%hgt_sum(k,l))
                                              enddo
                                           enddo

                                           do r=1,nrens
                                              ana(s)%HA(q,t,r) = sum(state(r,i,j,:,:,n) * tempwgts)
                                           enddo

                                           ! update observation counter
                                           nrobs = nrobs + 1

                                           ! update total (weighted observation counter by pft and hgt)
                                           ! do not count observed pft's which are not assimilated
                                           totobs(select_pft) = totobs(select_pft) + &
                                                obs(s)%pft_pct(k,l,select_pft) / obs(s)%pft_sum(k,l)

                                        end if ! valid observation

                                     end if ! observation inside gridpoint i,j

                                  end do ! observation grid points 1..nx

                               end do ! observation grid points 1..ny

                            end do ! model grid points 1..nx

                         end do ! model grid points 1..ny

                      end if ! gridded or raw observations

                   end if ! state name == observation name

                end do ! loop over state names

             end if ! has observation data

          end if ! observation time

       end do ! observation streams

       deallocate(tempmean)
       deallocate(tempwgts)

       if (saveobs.or.regridobs) then
          deallocate(obsgrid%dat)
          deallocate(obsgrid%err)
          deallocate(obsgrid%totdat)
          deallocate(obsgrid%toterr)
          deallocate(obsgrid%pft_sum)
          deallocate(obsgrid%pft_pct)
          deallocate(obsgrid%hgt_sum)
          deallocate(obsgrid%hgt_pct)
          deallocate(obsgrid%np)
          deallocate(obsgrid%npmax)
       endif

    endif ! is prediction

    ! gather nrobs from all processes
    if (mpiid.eq.0) then
       gnrobs(1) = gnrobs(1) + nrobs
       cumobs = cumobs + totobs

#ifndef HIDE_MPI
       do proc=2,numprocs
          call MPI_RECV(nrobs,1,MPI_INTEGER,proc-1, proc,MPI_COMM_WORLD,mpistatus,mpierr)
          gnrobs(proc) = gnrobs(proc) + nrobs
          call MPI_RECV(totobs,npft,MPI_REAL,proc-1, proc,MPI_COMM_WORLD,mpistatus,mpierr)
          cumobs = cumobs + totobs
       end do

#endif

    else
#ifndef HIDE_MPI
       call MPI_SEND(nrobs,1,MPI_INTEGER,0,mpiid+1,MPI_COMM_WORLD,mpierr)
       call MPI_SEND(totobs,npft,MPI_REAL,0,mpiid+1,MPI_COMM_WORLD,mpierr)
#endif
    endif

#ifndef HIDE_MPI
    ! redistribute cumulative (for both processes and time) to all processes
    call MPI_BCAST(cumobs,npft,MPI_REAL,0,MPI_COMM_WORLD,mpierr)    
#endif

  end subroutine analysis_observation

  subroutine analysis_global

    ! global analysis of model parameters
    ! first process performs the analysis
    ! update is performed for all parameters globally (and if selected, for all states globally)

    use driver_module
    use time_module
    use enkf_module_r8
    use file_module
    use restart_module
    use pft_module

#ifndef HIDE_MPI
    use mpi
#endif

    implicit none

    ! local variables
    integer :: nrobs                      ! number of observations for this process
    logical :: doanalysis                 ! perform global analysis or not
    logical :: subset                     ! subset observations or not
    integer :: k, k0, k1
    integer :: h, i, j, l, n, o, s, p, q, r, t
    integer :: nrobstot
    integer, allocatable :: obsindextot(:)
    integer, allocatable :: obsindex(:)

    real(kind=8), allocatable :: HA(:,:)
    real(kind=8), allocatable :: obsdat(:)
    real(kind=8), allocatable :: obserr(:)
    real(kind=8), allocatable :: A(:,:)
    real(kind=8), allocatable :: ATEMP(:,:)

    real(kind=8) :: val_var, val_mean

    real(kind=8) :: varscale
    integer :: ar

    character(len=200) :: filename
    integer :: unit

#ifndef HIDE_MPI
    integer :: proc
    integer :: mpistatus(MPI_STATUS_SIZE)
    integer, allocatable :: tempindex(:)
    real(kind=8), allocatable :: tempval(:)
    integer :: e, nr
#endif

    if (((.not.yearlyanalysis).or.(oldyear)).and.globalanalysis) then

       doanalysis = .false.

       ! check if we have enough observations for analysis
       if (mpiid.eq.0) then
          nrobstot = sum(gnrobs)
          if (nrobstot.ge.nrobsmin) doanalysis = .true.
       endif

#ifndef HIDE_MPI
       call MPI_BCAST(doanalysis,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)    
#endif

       if (doanalysis) then

          ! check if we have to subset observation vector
          if (mpiid.eq.0) then

             ! first process: distribute nrobs and observation index
             if (nrobstot.gt.nrobsmax) then
                allocate(obsindextot(nrobsmax))
                do k=1,nrobsmax
                   obsindextot(k) = int(int(k,8)*int(nrobstot,8)/int(nrobsmax,8),4)
                end do
                nrobstot = min(nrobsmax,nrobstot)
                subset = .true.
             else                
                subset = .false.
             endif
          endif

#ifndef HIDE_MPI
          call MPI_BCAST(subset,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)    
#endif

          if (subset) then
             ! distribute subsetted observations among processes that contribute to global analysis
             if (mpiid.eq.0) then
#ifndef HIDE_MPI

                do proc=2,numprocs

                   if (gnrobs(proc).gt.0) then
                      k0 = -1
                      k1 = -1
                      do k=1,nrobstot 
                         if ((obsindextot(k).ge.(sum(gnrobs(1:proc-1))+1)).and.(k0.eq.-1)) k0=k
                      end do
                      do k=nrobstot,1,-1
                         if ((obsindextot(k).le.sum(gnrobs(1:proc))).and.(k1.eq.-1)) k1=k
                      end do

                      nrobs = k1-k0+1

                      if ((nrobs.lt.0).or.(k0.lt.0).or.(k1.lt.0)) nrobs=0 
                   else
                      nrobs = 0
                   endif

                   ! send new nrobs 
                   call MPI_SEND(nrobs,1,MPI_INTEGER,proc-1,proc,MPI_COMM_WORLD,mpierr)

                   if (nrobs.gt.0) then
                      ! send observations selected from this process
                      allocate(tempindex(nrobs))
                      tempindex = obsindextot(k0:k1) - sum(gnrobs(1:proc-1))
                      call MPI_SEND(tempindex,nrobs,MPI_INTEGER,proc-1,proc,MPI_COMM_WORLD,mpierr)
                      deallocate(tempindex)
                   endif

                end do
#endif


                ! first process
                if (gnrobs(1).gt.0) then
                   k0 = -1
                   k1 = -1
                   do k=1,nrobstot 
                      if ((obsindextot(k).ge.1).and.(k0.eq.-1)) k0=k
                   end do
                   do k=nrobstot,1,-1
                      if ((obsindextot(k).le.gnrobs(1)).and.(k1.eq.-1)) k1=k
                   end do

                   nrobs = k1-k0+1

                   if ((nrobs.lt.0).or.(k0.lt.0).or.(k1.lt.0)) nrobs=0 

                   if (nrobs.gt.0) then
                      allocate(obsindex(nrobs))
                      obsindex(:) = obsindextot(k0:k1)
                   endif
                else
                   ! reset nrobs of first process to saved nrobs
                   nrobs = gnrobs(1)
                endif

             else
#ifndef HIDE_MPI
                ! other processes: receive new nrobs and observation index
                call MPI_RECV(nrobs,1,MPI_INTEGER,0, mpiid+1,MPI_COMM_WORLD,mpistatus,mpierr)
                if (nrobs.gt.0) then
                   ! observations for global analysis selected for this process
                   allocate(obsindex(nrobs))
                   call MPI_RECV(obsindex,nrobs,MPI_INTEGER,0, mpiid+1,MPI_COMM_WORLD,mpistatus,mpierr)
                endif
#endif
             endif

          else
             ! no subsetting: send observation count to processes
             if (mpiid.eq.0) then
                nrobs = gnrobs(1)
#ifndef HIDE_MPI
                do proc=2,numprocs
                   call MPI_SEND(gnrobs(proc),1,MPI_INTEGER,proc-1,proc,MPI_COMM_WORLD,mpierr)
                end do
#endif
             else
#ifndef HIDE_MPI
                call MPI_RECV(nrobs,1,MPI_INTEGER,0, mpiid+1,MPI_COMM_WORLD,mpistatus,mpierr)
#endif
             endif

          endif ! subset

          if (mpiid.eq.0) then
             ! allocate analysis matrices
             allocate(HA(nrobstot,nrens))
             allocate(obsdat(nrobstot))
             allocate(obserr(nrobstot))
          else
             if (nrobs.gt.0) then
                allocate(HA(nrobs,nrens))
                allocate(obsdat(nrobs))
                allocate(obserr(nrobs))
             endif
          endif

          allocate(A(1,nrens))
          allocate(ATEMP(1,nrens))

          ! fill analysis matrices
          r = 1
          o = 1
          if (nrobs.gt.0) then
             do s=1,nselect_obser
                do t=1,ana(s)%nt
                   do q=1,ana(s)%np
                      if (ana(s)%obsdat(q,t).ne.nodata) then
                         if (subset) then
                            if (obsindex(min(r,nrobs)).eq.o) then
                               HA(r,:) = real(ana(s)%HA(q,t,:),kind=8)
                               obsdat(r) = real(ana(s)%obsdat(q,t),kind=8)
                               obserr(r) = real(ana(s)%obserr(q,t),kind=8)
                               r=r+1
                            endif
                         else
                            HA(r,:) = real(ana(s)%HA(q,t,:),kind=8)
                            obsdat(r) = real(ana(s)%obsdat(q,t),kind=8)
                            obserr(r) = real(ana(s)%obserr(q,t),kind=8)
                            r=r+1
                         endif
                         o=o+1                         
                      endif
                   end do
                end do
             end do
          endif

#ifndef HIDE_MPI
          do proc=2,numprocs
             if (mpiid.eq.0) then
                call MPI_RECV(nr,1,MPI_INTEGER,proc-1, proc,MPI_COMM_WORLD,mpistatus,mpierr)
                if (nr.gt.0) then
                   allocate(tempval(nr))
                   do e=1,nrens
                      call MPI_RECV(tempval,nr,MPI_DOUBLE_PRECISION,proc-1, &
                           proc,MPI_COMM_WORLD,mpistatus,mpierr)
                      HA(r:r+nr-1,e) = tempval
                   end do
                   call MPI_RECV(tempval,nr,MPI_DOUBLE_PRECISION,proc-1, &
                        proc,MPI_COMM_WORLD,mpistatus,mpierr)
                   obsdat(r:r+nr-1) = tempval
                   call MPI_RECV(tempval,nr,MPI_DOUBLE_PRECISION,proc-1, &
                        proc,MPI_COMM_WORLD,mpistatus,mpierr)
                   obserr(r:r+nr-1) = tempval
                   deallocate(tempval)
                   r = r + nr
                endif
             else
                if ((proc-1).eq.mpiid) then
                   call MPI_SEND(nrobs,1,MPI_INTEGER,0,mpiid+1,MPI_COMM_WORLD,mpierr)
                   if (nrobs.gt.0) then
                      do e=1,nrens
                         call MPI_SEND(HA(:,e),nrobs,MPI_DOUBLE_PRECISION,0,mpiid+1,MPI_COMM_WORLD,mpierr)
                      end do
                      call MPI_SEND(obsdat,nrobs,MPI_DOUBLE_PRECISION,0,mpiid+1,MPI_COMM_WORLD,mpierr)
                      call MPI_SEND(obserr,nrobs,MPI_DOUBLE_PRECISION,0,mpiid+1,MPI_COMM_WORLD,mpierr)
                   endif
                endif
             endif

             ! only proceed to exchange of HA, obsdat, obserr of next process when current 
             ! process has finished in order to avoid MPI buffer overflow (problem was before:
             ! all processes send HA, obsdat and obserr to process 0)
             call synchronize_processes

          end do
#endif

          ! run the analysis (first process only
          if (mpiid.eq.0) then
             allocate(ROT(nrens,nrens))
             call calc_rot(nrens)
             call enkf(HA, obsdat, obserr, 23, X=Xglobal)
             deallocate(ROT)
          endif

          ! clean up
          if (subset) then
             if (nrobs.gt.0) then
                deallocate(obsindex)
             endif
             if (mpiid.eq.0) then
                deallocate(obsindextot)
             endif
          endif

          if (mpiid.eq.0) then
             deallocate(HA)
             deallocate(obsdat)
             deallocate(obserr)
          else
             if (nrobs.gt.0) then
                deallocate(HA)
                deallocate(obsdat)
                deallocate(obserr)
             endif
          endif

#ifndef HIDE_MPI
          ! send analysis matrix to all processes
          call MPI_BCAST(Xglobal,nrens*nrens,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif

          ! retrieve intrinsic analysis variance deflation factor
          varscale = 0.d0
          do ar = 1, nrand

             A(1,:) = real(Arand(ar,:),kind=8)
             ATEMP = A

             call dgemm('n','n', 1, nrens, nrens, &
                  1.d0, ATEMP, 1, Xglobal, nrens, 0.d0, A, 1)

             ! posterior ensemble variance
             call ensmean(A(1,:),val_mean)
             call ensvar(A(1,:),val_mean,val_var)

             varscale = varscale + val_var

          end do
          varscale = real(nrand,kind=8) / varscale

          l=0 ! parameter counter

          ! update global parameters with global analysis matrix
          if (nptparam*nselect_param.gt.0) then

             do n=1,nselect_param
                ! update all parameters for all pft's (not only selected pft's)
                ! in order to obtain a global parameter set on each process
                do p=1,npft
                   ! assign parameter ensemble
                   A(1,:) = real(param(:,p,select_param(n)),kind=8)

                   ! prior ensemble variance
                   call ensmean(A(1,:),val_mean)
                   call ensvar(A(1,:),val_mean,val_var)

                   ! update ensemble with analysis
                   ATEMP = A
                   call dgemm('n','n', 1, nrens, nrens, &
                        1.d0, ATEMP, 1, Xglobal, nrens, 0.d0, A, 1)

                   ! rescale, inflate & deflate posterior ensemble variance
                   call kernel_rescale(A(1,:),varscale)
                   call kernel_inflate(A(1,:),val_var,real(inflation,kind=8))
                   call kernel_deflate(A(1,:),val_var,real(deflation,kind=8))

                   ! limit posterior ensemble mean to min/max values
                   !call kernel_limit(A(1,:),real(param_min(select_param(n)),kind=8), &
                   !     real(param_max(select_param(n)),kind=8))

                   ! return updated parameter ensemble
                   param(:,p,select_param(n)) = real(A(1,:),kind=4)

                   l=l+1

                end do
             end do

             if (sepminmax) then
                do n=1,nselect_param-1
                   ! separate parameters so that the mean of the xyzMIN is lower than the
                   ! mean of xyzMAX. The xyzMIN and xyzMAX need to be listed one after the other
                   ! in the parameters file.

                   if ((len(trim(param_names(select_param(n)))).gt.3).and. &
                        (len(trim(param_names(select_param(n)))).eq. &
                        len(trim(param_names(select_param(n+1))))).and. &
                        (select_param(n).eq.(select_param(n+1)-1))) then

                      if (index(param_names(select_param(n)),'MIN',.true.).eq. &
                           index(param_names(select_param(n+1)),'MAX',.true.)) then

                         if (verbose) write(*,'(A,A,A,A)') 'Separate ensemble mean of Parameters : ', &
                              trim(param_names(select_param(n))),' / ',trim(param_names(select_param(n+1)))

                         do p=1,npft
                           call separate_minmax(param(:,p,select_param(n)),&
                                 param(:,p,select_param(n+1)), &
                                 optmindiff=0.1*param_var(p,select_param(n))) 
                         end do

                      endif

                   endif

                end do
             endif

          endif


          k=0 ! state counter

          if (isprediction) then

             ! update local states with global analysis matrix
             if (nptstate*nselect_state.gt.0) then

                do n=1,nselect_state
                   do h=1,nselect_hgt
                      do p=1,nselect_pft
                         do j=1,geographic%ny
                            do i=1,geographic%nx
                               A(1,:) = real(state(:,i,j,p,h,select_state(n)),kind=8)

                               ! prior ensemble variance
                               call ensmean(A(1,:),val_mean)
                               call ensvar(A(1,:),val_mean,val_var)

                               ! update ensemble with analysis
                               ATEMP = A
                               call dgemm('n','n', 1, nrens, nrens, &
                                    1.d0, ATEMP, 1, Xglobal, nrens, &
                                    0.d0, A, 1)

                               ! rescale,  inflate & deflate posterior ensemble variance
                               call kernel_rescale(A(1,:),varscale)
                               call kernel_inflate(A(1,:),val_var,real(inflation,kind=8))
                               call kernel_deflate(A(1,:),val_var,real(deflation,kind=8))

                               ! limit posterior ensemble mean to min/max values
                               ! call kernel_limit(A(1,:),real(state_min(select_state(n)),kind=8, &
                               !     real(state_max(select_state(n)),kind=8))

                               ! return updated parameter ensemble
                               state(:,i,j,p,h,select_state(n)) = real(A(1,:),kind=4)

                               k=k+1

                            end do
                         end do
                      end do
                   end do
                end do

             endif

          endif

          ! clean up
          deallocate(A)
          deallocate(ATEMP)

          if ((mpiid.eq.0).and.(verbose)) then
             write(*,*)
             write(*,'(A)')    'Global Analysis:'
             write(*,'(A,I8)') 'Total number of observations used:    ',sum(gnrobs)
             write(*,'(A,I8)') 'Total number of innovations analyzed: ',nrobstot
             write(*,'(A,I8)') 'Total number of parameters updated:   ',l
             write(*,'(A,I8)') 'Total number of states updated:       ',k
          endif

          if (mpiid.eq.0) write(*,*)

          ! write analyzed parameters to ascii file
          call write_parameters

          ! if we use the restart date, then only write restart after each analysis
          if (restartdate) then

             ! write restart with analyzed states/forcings/parameters
             call write_restart

             ! write analysis date
             if (mpiid.eq.0) then
                unit = getunit()
                filename  = trim(outdir)//trim(geographic%name)//'.analysis-date.dat'
                if (verbose) write(*,'(A,A)') 'Writing analysis date to file : ',trim(filename)
                open(unit=unit, file=trim(filename), form='formatted',action='write')
                write(unit,'(I4,1X,I2,1X,I2)') year,month,day
                close(unit)
             endif

          endif

       endif ! global analysis

    endif

  end subroutine analysis_global

  subroutine analysis_local

    ! local analysis of model states
    ! update is performed by grid point with observations localized by grid point
    ! the same observation can be used by several grid point if within the influence area

    use driver_module
    use time_module
    use enkf_module_r8

#ifndef HIDE_MPI
    use mpi
#endif

    implicit none

    ! local variables
    integer :: nrobs
    integer :: i, j, h, k, n, o, p, q, r, s, t

    real(kind=8), allocatable :: HA(:,:)
    real(kind=8), allocatable :: obsdat(:)
    real(kind=8), allocatable :: obserr(:)
    real(kind=8), allocatable :: Xlocal(:,:)
    real(kind=8), allocatable :: A(:,:)
    real(kind=8), allocatable :: ATEMP(:,:)

    real(kind=8) :: val_var, val_mean

    real(kind=8) :: varscale
    integer :: ar

#ifndef HIDE_MPI
    integer :: proc
    integer :: mpistatus(MPI_STATUS_SIZE)
#endif

    if (((.not.yearlyanalysis).or.(oldyear)).and.localanalysis) then

       ! send observation count to processes
       if (mpiid.eq.0) then
          nrobs = gnrobs(1)
#ifndef HIDE_MPI
          do proc=2,numprocs
             call MPI_SEND(gnrobs(proc),1,MPI_INTEGER,proc-1,proc,MPI_COMM_WORLD,mpierr)
          end do
#endif
       else
#ifndef HIDE_MPI
          call MPI_RECV(nrobs,1,MPI_INTEGER,0, mpiid+1,MPI_COMM_WORLD,mpistatus,mpierr)
#endif
       endif

       ! perform ensemble rotation
       allocate(ROT(nrens,nrens))
       if (mpiid.eq.0) then
          if (sum(gnrobs).ge.nrobsmin) then
             call calc_rot(nrens)
          endif
       endif

#ifndef HIDE_MPI
       ! send analysis matrix to all processes
       call MPI_BCAST(ROT,nrens*nrens,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif

       ! check if we have enough observations for analysis
       if ((nrobs.ge.nrobsmin).and.(isprediction)) then

          k=0  ! state counter
          o=0  ! observation counter

          ! cycle through model grid points and perform local analysis per grid point
          do j=1,geographic%ny
             do i=1,geographic%nx

                nrobs = 0
                do s=1,nselect_obser
                   do t=1,ana(s)%nt
                      nrobs = nrobs + &
                           count((ana(s)%xid.eq.i).and.(ana(s)%yid.eq.j).and.(ana(s)%obsdat(:,t).ne.nodata))
                   end do
                end do

                if (nrobs.ge.nrobsmin) then

                   ! allocate analysis arrays
                   allocate(HA(nrobs,nrens))
                   allocate(obsdat(nrobs))
                   allocate(obserr(nrobs))
                   allocate(A(1,nrens))
                   allocate(ATEMP(1,nrens))
                   allocate(Xlocal(nrens,nrens))

                   ! fill analysis arrays
                   r = 1
                   do s=1,nselect_obser
                      do t=1,ana(s)%nt
                         do q=1,ana(s)%np
                            if ((ana(s)%xid(q).eq.i).and.&
                                 (ana(s)%yid(q).eq.j).and.&
                                 (ana(s)%obsdat(q,t).ne.nodata)) then
                               obsdat(r) = real(ana(s)%obsdat(q,t),kind=8)
                               obserr(r) = real(ana(s)%obserr(q,t),kind=8)
                               HA(r,:) = real(ana(s)%HA(q,t,:),kind=8)
                               r = r + 1
                            endif
                         end do
                      end do
                   end do

                   o = o + nrobs

                   ! perform analysis
                   call enkf(HA, obsdat, obserr, 23, X=Xlocal)

                   ! retrieve intrinsic analysis variance deflation factor
                   varscale = 0.d0
                   do ar = 1, nrand

                      A(1,:) = real(Arand(ar,:),kind=8)
                      ATEMP = A

                      call dgemm('n','n', 1, nrens, nrens, &
                           1.d0, ATEMP, 1, Xlocal, nrens, 0.d0, A, 1)

                      ! posterior ensemble variance
                      call ensmean(A(1,:),val_mean)
                      call ensvar(A(1,:),val_mean,val_var)

                      varscale = varscale + val_var

                   end do
                   varscale = real(nrand) / varscale

                   ! update only the difference between global and local analysis
                   ! if global analysis was performed before
                   if (globalanalysis) Xlocal = Xlocal - Xglobal

                   ! update states
                   do n=1,nselect_state
                      do h=1,nselect_hgt
                         do p=1,nselect_pft

                            A(1,:) = real(state(:,i,j,p,h,select_state(n)),kind=8)

                            ! prior ensemble variance
                            call ensmean(A(1,:),val_mean)
                            call ensvar(A(1,:),val_mean,val_var)

                            ! update ensemble with analysis
                            ATEMP = A
                            call dgemm('n','n', 1, nrens, nrens, &
                                 1.d0, ATEMP, 1, Xlocal, nrens, &
                                 0.d0, A, 1)

                            ! add differences of local analysis back to global analysis
                            if (globalanalysis) A = ATEMP + A

                            ! rescale, inflate & deflate posterior ensemble variance
                            call kernel_rescale(A(1,:),varscale)
                            call kernel_inflate(A(1,:),val_var,real(inflation,kind=8))
                            call kernel_deflate(A(1,:),val_var,real(deflation,kind=8))

                            ! limit posterior ensemble mean to min/max values
                            call kernel_limit(A(1,:),real(state_min(select_state(n)),kind=8), &
                                 real(state_max(select_state(n)),kind=8))

                            ! return updated parameter ensemble
                            state(:,i,j,p,h,select_state(n)) = real(A(1,:),kind=4)

                            k=k+1

                         end do
                      end do
                   end do


                   ! clean up
                   deallocate(Xlocal)
                   deallocate(ATEMP)
                   deallocate(A)
                   deallocate(HA)
                   deallocate(obsdat)
                   deallocate(obserr)

                endif

             end do
          end do

          if (verbose) then
             write(*,*)
             write(*,'(A,I8)') 'Local Analysis: ',mpiid
             write(*,'(A,I8)') 'Total number of observations used:    ',nrobs
             write(*,'(A,I8)') 'Total number of innovations analyzed: ',o
             write(*,'(A,I8)') 'Total number of states updated:       ',k
          end if
       endif

       deallocate(ROT)

    endif

  end subroutine analysis_local

  subroutine synchronize_processes
    ! this routine catches up with all processes on MPI runs
    ! should be called regularly in order to avoid timeouts
    ! on certain MPI implementations

#ifndef HIDE_MPI

    use mpi
    use driver_module

    call MPI_Barrier(MPI_COMM_WORLD,mpierr)

#endif

  end subroutine synchronize_processes

  subroutine ensemble_perturb(vector,variance)
    ! perturbs an ensemble vector with given variance

    use enkf_module_r8

    implicit none

    ! arguments
    real(kind=4), intent(inout) :: vector(:)
    real(kind=4), intent(in) :: variance

    ! local variables
    real(kind=8), allocatable :: e_rand(:)
    real(kind=8) :: e_mean
    real(kind=8) :: e_var
    integer :: nrens

    nrens=size(vector,1)

    allocate(e_rand(nrens))

    call random(e_rand,nrens)

    ! evaluate 'real' mean/variance of random numbers (should be close to 0.0/1.0)
    call ensmean(e_rand,e_mean)
    call ensvar(e_rand,e_mean,e_var)

    ! normalize to variance of 1.0 and mean of 0.0
    e_rand = (e_rand - e_mean)/sqrt(e_var)

    ! scale variance and add to mean
    vector = vector + real(e_rand,kind=4) * sqrt(variance)

    deallocate(e_rand)

  end subroutine ensemble_perturb

end module assimilation_module
