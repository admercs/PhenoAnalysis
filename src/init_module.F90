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

module init_module

contains

  subroutine init_driver(namelistfile,sitefile,sitelist,date_curr)

    ! Subroutine initializing the driver arrays and reading the namelist

    use driver_module
    use pft_module
    use topo_module
    use phenology_module
    use modis_module
    use meteorology_module
    use assimilation_module
    use time_module
    use file_module
    use exit_module
    use proj_module
    use quicksort_module
    use parameter_module

#ifndef HIDE_MPI
    use mpi
#endif

    implicit none

    ! arguments
    character(len=*), intent(in) :: namelistfile
    character(len=*), intent(in) :: sitefile
    character(len=*), intent(in) :: sitelist
    character(len=*), intent(in), optional :: date_curr

    ! local variables
    integer :: i, f, h, j, k, n, o, p, r, s, v
    integer :: i0,i1,j0,j1
    character(len=100) :: filename
    integer c, pos1, pos2, nchars, unit, status
    real(kind=4) :: fval
    integer :: nlin, ncol
    integer, allocatable :: sel(:), out(:)
    real(kind=4) :: minselpft
    real(kind=4) :: minselhgt
    integer :: proc
    integer :: minpft, maxpft
    integer :: site0, site1, nsites, site
    integer, allocatable :: sites(:)
    character(100), allocatable :: table(:,:)
  
    real(kind=4)  :: maxhgt
    real(kind=4)  :: minhgt

    integer, allocatable :: indexnpt(:)
    integer :: npt
    real(kind=4) :: xbounds(4,1), ybounds(4,1), lonbounds(4,1), latbounds(4,1)
    real(kind=4), allocatable :: xmin(:),xmax(:),ymin(:),ymax(:)

    real(kind=4), allocatable :: waterfrac(:,:)     ! water fraction (1: water, 0: land)
    real(kind=4), allocatable :: waterfrac1d(:,:)
    logical, allocatable :: landmask(:,:)           ! land mask (true=land, false=water)
    logical, allocatable :: validmask(:,:)          
    integer :: nlandproc                            ! land tile counter

#ifndef HIDE_MPI
    ! MPI variables
    integer :: mpistatus(MPI_STATUS_SIZE)
    real(kind=4), allocatable :: pft_tmp(:)
#endif

    ! define model namelist
    namelist /model/ &
         regional, &
         original, &
         nsites, &
         nrens, &
         nhgt, &
         use_forcing, &
         use_chilling, &
         use_moisture, &
         use_light, &
         use_rainfall, &
         use_photo, &
         averagepft, &
         averagehgt, &
         averageens, &
         compressens, &
         restart, &
         restartdate, &
         filelength, &
         verbose

    namelist /assimilation/ &
         analysis, &
         saveobs, &
         regridobs, &
         useobspct, &
         simulator, &       
         yearlyanalysis, &
         globalanalysis, &
         localanalysis, &
         influence, &
         weightmag, &
         weightdist, &
         weightarea, &
         weightgrid, &
         sepminmax, &
         nrobsmax, &
         nrobsmin, &
         deflation, &
         inflation

    namelist /spatial/ &
         landcover_screen, &
         pft_screen, &
         nxprocs, &
         nyprocs, &
         subgrid

    namelist /temporal/ &
         date_start, &
         date_end, &
         dtmod, &
         dtphen, &
         dtout

    namelist /data/ &
         meteorologytype, &
         pfttype, &
         landcovertype, &
         topotype, &
         satellitetype

    namelist /paths/ &
         projectdir, &
         datadir, &
         meteorologydir, &
         pftdir, &
         landcoverdir, &
         topodir, &
         satellitedir

    ! retrieve package name and version
    package_name = "PACKAGE_NAME"
    status = get_config(package_name)
    package_version = "PACKAGE_VERSION"
    status = get_config(package_version)

    ! MANUALLY SET SOME FLAGS WHICH HAVE NOT YET PROPAGATED INTO THE NAMELIST
    hgtscreen = .false.  ! screen observations by elevation or not

    ! set initial flag
    initial = .true.

    ! minimum % coverage (e.g. minimum land area fraction if water pixels are not assimilated)
    ! to make it into the state and observation vector (only pixels are assimilated that are 
    ! well represented within grid, this minimizes noise contribution to the assimilation)
    minsumpft = 0.75
    minsumhgt = 0.75

    if (mpiid.eq.0) then
       write(*,'(A,A,A,A)') 'Package ',trim(package_name),' Version ',trim(package_version)
       write(*,*)
       write(*,'(A)') 'Copyright (C) 2006 - 2014 Reto Stockli (Blue Marble Research)'
       write(*,'(A)') 'This program comes with ABSOLUTELY NO WARRANTY'
       write(*,'(A)') 'This is free software, and you are welcome to redistribute it'
       write(*,'(A)') 'under the terms of the GNU General Public License as published by'
       write(*,'(A)') 'the Free Software Foundation, either version 3 of the License, or'    
       write(*,'(A)') 'any later version. For details see http://www.gnu.org/licenses/gpl.html.'
       write(*,*)
    endif

    !--------------------------
    ! NAMELISTS
    !--------------------------

    !! evaluate project case from namelist filename
    pos2 = index(namelistfile,"/",back=.true.)
    namelistdir = namelistfile(1:pos2)
    pos1 = index(namelistfile(1:pos2-1),"/",back=.true.)
    casename = namelistfile(pos1+1:pos2-1)

    ! read namelist

    unit = getunit()
    ! read namelists
    open(unit=unit,file=namelistfile,form='formatted')
    read(unit,model)
    read(unit,assimilation)
    read(unit,spatial)
    read(unit,temporal)
    read(unit,data)
    read(unit,paths)
    close(unit)

    if (present(date_curr)) then
       ! for resubmission, replace full date range with the date range
       ! of current resubission
       if (trim(date_curr).ne."") then
          date_start = trim(date_curr)
          date_end = trim(date_curr)
       end if
    end if

    ! only print stuff on first process
    if (mpiid.ne.0) verbose = .false.

    ! cannot save observations (since not reading observations) in prediction-only mode
    if (.not.analysis) saveobs = .false.

    if (nrens.eq.1) then
       ! scalar simulation (no ensemble member dimensions on output)
       averageens = .true. 
       compressens = .false.
    endif
    
    if (compressens) then
       ! if ensemble compression is chosen, by default also average the ensemble
       ! members to the ensemble mean. The compressed ensemble will create a second
       ! output variable for each variable with the ensemble standard deviation
       averageens = .true.
    endif

    ! define data and output directories, add absolute path and trailing slash
    pos1 = index(trim(datadir),"/",back=.true.)
    if (pos1.ne.len_trim(datadir)) datadir = trim(datadir)//"/"

    pos1 = index(meteorologydir,"/")
    if (pos1.ne.1) meteorologydir = trim(datadir)//trim(meteorologydir)
    pos1 = index(trim(meteorologydir),"/",back=.true.)
    if (pos1.ne.len_trim(meteorologydir)) meteorologydir = trim(meteorologydir)//"/"

    pos1 = index(topodir,"/")
    if (pos1.ne.1) topodir = trim(datadir)//trim(topodir)
    pos1 = index(trim(topodir),"/",back=.true.)
    if (pos1.ne.len_trim(topodir)) topodir = trim(topodir)//"/"

    pos1 = index(pftdir,"/")
    if (pos1.ne.1) pftdir = trim(datadir)//trim(pftdir)
    pos1 = index(trim(pftdir),"/",back=.true.)
    if (pos1.ne.len_trim(pftdir)) pftdir = trim(pftdir)//"/"

    pos1 = index(satellitedir,"/")
    if (pos1.ne.1) satellitedir = trim(datadir)//trim(satellitedir)
    pos1 = index(trim(satellitedir),"/",back=.true.)
    if (pos1.ne.len_trim(satellitedir)) satellitedir = trim(satellitedir)//"/"

    pos1 = index(landcoverdir,"/")
    if (pos1.ne.1) landcoverdir = trim(datadir)//trim(landcoverdir)
    pos1 = index(trim(landcoverdir),"/",back=.true.)
    if (pos1.ne.len_trim(landcoverdir)) landcoverdir = trim(landcoverdir)//"/"

    pos1 = index(trim(projectdir),"/",back=.true.)
    if (pos1.ne.len_trim(projectdir)) projectdir = trim(projectdir)//"/"

    batchdir = trim(projectdir)//trim(casename)//'/batch/'
    outdir = trim(projectdir)//trim(casename)//'/output/'
    restartdir = trim(projectdir)//trim(casename)//'/restart/'

    ! define number of subgrid-scale divisions and pft and height classes
    if (regional) then
       if (mpiid.eq.0) write(*,'(A)') 'Model running in Regional Mode'

       ! check for topo and pft files and set total pft number (numpft)
       call init_topo
       call init_pft
       npft = numpft
    else
       if (mpiid.eq.0) write(*,'(A)') 'Model running in Local Mode'
       npft = 1
       nhgt = 1
    endif

    ! evaluate site string and read site information
    ! first syntax option range of sites: a..c
    pos1 = index(sitelist,"..")
    if ((pos1.gt.1).and.(pos1.le.(len(sitelist)-3))) then
       read(sitelist(1:pos1-1),*) site0
       read(sitelist(pos1+2:len(sitelist)),*) site1
       nsites = site1-site0+1
       allocate(sites(nsites))
       do c=1,nsites
          sites(c) = site0 + c - 1
       end do
    else
       ! second syntax option list of sites: a,b,c (including single site)
       nchars = len(sitelist)
       nsites = 0
       do c=1,len(sitelist)
          if (sitelist(c:c).eq.',') nsites = nsites + 1
       end do
       nsites = nsites + 1
       allocate(sites(nsites))

       ! tokenize read string by comma-separated columns
       c=1
       pos1=1
       pos2=index(sitelist(pos1:nchars),',')
       if (pos2.gt.0) then 
          pos2=pos2+pos1-2
       else 
          pos2=nchars
       end if

       ! loop through columns 
       do while (pos1.gt.0)
          if (pos2.ge.pos1) then
             if (c.le.nsites) read(sitelist(pos1:pos2),*) sites(c)
          else
             pos2=pos1
          endif

          pos1=index(sitelist(pos2:nchars),',')
          if (pos1.gt.0) then 
             pos1=pos1+pos2 
             pos2=index(sitelist(pos1:nchars),',')
             if (pos2.gt.0) then
                pos2=pos2+pos1-2 
             else 
                pos2=nchars
             end if
          endif
          c=c+1
       end do
    end if

    call handle_error(nsites.gt.numprocs,'# processes need to be equal or larger than # of regions (sites)')

    ! read site case file
    call read_csv_table(sitefile,table,',',2)
      
    if (nsites.eq.1) then
       site = sites(1)
    else
       if (mpiid.lt.nsites) then
          site = sites(mpiid+1)
       else
          site = 0
       end if
    end if

    if (site.gt.0) then
       do s=1,size(table,2)
          read(table(1,s),*) site0
          if (site0.eq.site) then
             geographic%name = trim(table(3,s))
             read(table(4,s),*) geographic%xmin
             read(table(5,s),*) geographic%xmax
             read(table(6,s),*) geographic%ymin
             read(table(7,s),*) geographic%ymax
             read(table(8,s),*) geographic%dx
             read(table(9,s),*) geographic%dy
             geographic%type = trim(table(11,s))
             geographic%optargs = trim(table(12,s))
          endif
       end do

       if (geographic%dx.lt.eps) then
          if (geographic%type.eq."latitude_longitude") geographic%dx = 0.05
          if (geographic%type.eq."rotated_latitude_longitude") geographic%dx = 0.05
          if (geographic%type.eq."swisscors") geographic%dx = 25.0
          if (geographic%type.eq."geosat") geographic%dx = 1.0
          if (geographic%dx.lt.eps) stop
          geographic%xmin = geographic%xmin - 0.5*geographic%dx
          geographic%xmax = geographic%xmax + 0.5*geographic%dx
       end if
       if (geographic%dy.lt.eps) then
          if (geographic%type.eq."latitude_longitude") geographic%dy = 0.05
          if (geographic%type.eq."rotated_latitude_longitude") geographic%dy = 0.05
          if (geographic%type.eq."swisscors") geographic%dy = 25.0
          if (geographic%type.eq."geosat") geographic%dy = 1.0
          if (geographic%dy.lt.eps) stop
          geographic%ymin = geographic%ymin - 0.5*geographic%dy
          geographic%ymax = geographic%ymax + 0.5*geographic%dy
       end if

       geographic%nx = nint((geographic%xmax - geographic%xmin)/geographic%dx)
       geographic%ny = nint((geographic%ymax - geographic%ymin)/geographic%dy)

       ! derive full geographic grid before subdividing into MPI processes
       geographic_full = geographic
       geographic_full%dx = geographic_full%dx/real(subgrid,kind=4)
       geographic_full%dy = geographic_full%dy/real(subgrid,kind=4)
       geographic_full%nx = geographic_full%nx*subgrid
       geographic_full%ny = geographic_full%ny*subgrid

       ! create full domain output grid on each process
       allocate(xcoord_full(geographic_full%nx))
       allocate(ycoord_full(geographic_full%ny))
       do i=1,geographic_full%nx
          xcoord_full(i) = geographic_full%xmin + geographic_full%dx*(real(i)-0.5)
       end do
       do j=1,geographic_full%ny
          ycoord_full(j) = geographic_full%ymin + geographic_full%dy*(real(j)-0.5)
       end do
    else
       geographic%name = ""
       geographic%xmin = nodata
       geographic%xmax = nodata
       geographic%ymin = nodata
       geographic%ymax = nodata
       geographic%dx = nodata
       geographic%dy = nodata
       geographic%nx = 0
       geographic%ny = 0
       geographic%type = ""
       geographic%optargs = ""

       geographic_full = geographic
    end if

    deallocate(table)
    deallocate(sites)

    ! We currently support two options for MPI Parallelization
    ! a) multiple region assimilation or single region 
    ! b) single region assimilation with no subdivision

    if ((nsites.gt.1).or.((nsites.eq.1).and.(nxprocs.eq.1).and.(nyprocs.eq.1).and.(numprocs.eq.1))) then
       ! a) number of selected regions >= 1 and nxprocs/nyprocs == 1

       if (site.ne.0) then
          isprediction = .true.
       else
          isprediction = .false.
       end if

    else
       ! b) single region distributed on many processes
       ! to be used to assimilate a large region on many cpu's
       ! region will be subdivided and land areas will be filled into available process space

       call handle_error((geographic%nx.lt.nxprocs).or.(geographic%ny.lt.nyprocs),&
            'nx/ny need to exceed nxprocs/nyprocs')

       if (mpiid.eq.0) then

          ! read global grid-cell mean pft distribution
          allocate(waterfrac(geographic%nx,geographic%ny))

          ! distribute grid cells among processes
          allocate(indexnpt(numprocs))
          do proc=1,numprocs
             indexnpt(proc) = proc*geographic%nx*geographic%ny/numprocs
          end do

          npt = indexnpt(1)

#ifndef HIDE_MPI
          do proc=2,numprocs
             call MPI_SEND(indexnpt(proc)-indexnpt(proc-1),1,MPI_INTEGER,proc-1,proc,MPI_COMM_WORLD,mpierr)
          enddo
#endif

          if (npt.gt.0) then
             allocate(xmin(npt))
             allocate(xmax(npt))
             allocate(ymin(npt))
             allocate(ymax(npt))
          endif

          proc = 1
          k = 1
          do j = 1,geographic%ny
             do i = 1,geographic%nx
                do while ((k.gt.indexnpt(proc)).and.(proc.lt.numprocs)) 
                   proc = proc + 1
                end do
                if (proc.eq.1) then
                   xmin(k) = geographic%xmin + geographic%dx*real(i-1,kind=4)
                   xmax(k) = geographic%xmin + geographic%dx*real(i,kind=4)
                   ymin(k) = geographic%ymin + geographic%dy*real(j-1,kind=4)
                   ymax(k) = geographic%ymin + geographic%dy*real(j,kind=4)
                 else
#ifndef HIDE_MPI
                   call MPI_SEND(geographic%xmin + geographic%dx*real(i-1,kind=4),1,MPI_REAL,proc-1,proc,MPI_COMM_WORLD,mpierr)
                   call MPI_SEND(geographic%xmin + geographic%dx*real(i,kind=4),1,MPI_REAL,proc-1,proc,MPI_COMM_WORLD,mpierr)
                   call MPI_SEND(geographic%ymin + geographic%dy*real(j-1,kind=4),1,MPI_REAL,proc-1,proc,MPI_COMM_WORLD,mpierr)
                   call MPI_SEND(geographic%ymin + geographic%dy*real(j,kind=4),1,MPI_REAL,proc-1,proc,MPI_COMM_WORLD,mpierr)
#endif        
                endif
                k = k + 1
             end do
          end do

       else
#ifndef HIDE_MPI
          call MPI_RECV(npt,1,MPI_INTEGER,0,mpiid+1,MPI_COMM_WORLD,mpistatus,mpierr)
          if (npt.gt.0) then
             allocate(xmin(npt))
             allocate(xmax(npt))
             allocate(ymin(npt))
             allocate(ymax(npt))
             do k = 1,npt
                call MPI_RECV(xmin(k),1,MPI_REAL,0,mpiid+1,MPI_COMM_WORLD,mpistatus,mpierr)
                call MPI_RECV(xmax(k),1,MPI_REAL,0,mpiid+1,MPI_COMM_WORLD,mpistatus,mpierr)              
                call MPI_RECV(ymin(k),1,MPI_REAL,0,mpiid+1,MPI_COMM_WORLD,mpistatus,mpierr)
                call MPI_RECV(ymax(k),1,MPI_REAL,0,mpiid+1,MPI_COMM_WORLD,mpistatus,mpierr)              
             end do
          endif
#endif        
       endif
 
       if (npt.gt.0) then
          allocate(waterfrac1d(npt,1))
          if (regional) then
             do k=1,npt
                call read_pft(geographic%type,xmin(k),xmax(k),ymin(k),ymax(k), &
                     geographic%dx,geographic%dy,geographic%optargs, &
                     nodata,pft2d=waterfrac1d(k:k,:),pftsel=npft)  
             end do
          else
             waterfrac1d = 0.0
          endif

#ifndef HIDE_MPI
          if (mpiid.ne.0) then
             do k=1,npt
                call MPI_SEND(waterfrac1d(k,1),1,MPI_REAL,0,mpiid+1,MPI_COMM_WORLD,mpierr)
             end do
          endif
#endif        

          deallocate(xmin,xmax,ymin,ymax)
       endif
 
       if (mpiid.eq.0) then

          proc = 1
          k = 1
          do j = 1,geographic%ny
             do i = 1,geographic%nx
                do while ((k.gt.indexnpt(proc)).and.(proc.lt.numprocs)) 
                   proc = proc + 1
                end do
                if (proc.eq.1) then
                   waterfrac(i,j) = waterfrac1d(k,1)
!                   waterfrac(i,j) = 0.8 ! test for water tile
                else
#ifndef HIDE_MPI
                   call MPI_RECV(waterfrac(i,j),1,MPI_REAL,proc-1,proc,MPI_COMM_WORLD,mpistatus,mpierr)
#endif        
                endif
                k = k + 1
             end do
          end do

          deallocate(indexnpt)

       end if

       if (npt.gt.0) deallocate(waterfrac1d)

       ! create land mask
       allocate(landmask(nxprocs,nyprocs))

       if (mpiid.eq.0) then

          if (analysis) then
             maxpctwater = 0.5
          else
             maxpctwater = 1.0
          endif
          
          ! only process grid cells in domain with at least 5% land
          ! read in pft distribution and create water mask
          do yproc=1,nyprocs
             do xproc=1,nxprocs

                i0 = geographic%nx*(xproc-1)/nxprocs + 1
                j0 = geographic%ny*(yproc-1)/nyprocs + 1
                i1 = geographic%nx*(xproc)/nxprocs
                j1 = geographic%ny*(yproc)/nyprocs

                if (sum(waterfrac(i0:i1,j0:j1))/real((i1-i0+1)*(j1-j0+1)).gt.maxpctwater) then
                   landmask(xproc,yproc)=.false.
                else
                   landmask(xproc,yproc)=.true.
                endif

             enddo
          enddo

          deallocate(waterfrac)

          ! count number of land cells
          nlandproc = count(landmask)

          write(*,*)
          write(*,'(A,I6,I6)') "# nx/ny points (model):    ",geographic%nx,geographic%ny
          write(*,'(A,I6,I6)') "# nx/ny points (output):   ",geographic%nx*subgrid,geographic%ny*subgrid
          write(*,'(A,I6,I6)') "# x / y tiles:             ",nxprocs,nyprocs
          write(*,'(A,I6)')    "# regions with land area:  ",nlandproc
          write(*,'(A,I6)')    "# processes allocated:     ",numprocs
          write(*,*)

       endif

#ifndef HIDE_MPI
       ! distribute information on full domain grid among processes
       call MPI_BCAST(nlandproc,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
       call MPI_BCAST(landmask,nxprocs*nyprocs,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
#endif

       call handle_error(nlandproc.gt.numprocs,"# processes need to be equal or larger than # regions with land coverage.")

       call handle_error(numprocs.gt.(nxprocs*nyprocs),"# processes need to be equal or lower than nxproc x nyproc.")

       call handle_error(nlandproc.eq.0,"# regions with land coverage (>5% land) is 0. Water only area was chosen.")

       isprediction = .false.

       ! fill up land points
       proc = 0
       do yproc=1,nyprocs
          do xproc=1,nxprocs
             if (landmask(xproc,yproc)) then
                if ((proc+numprocs-nlandproc).eq.mpiid) then    
                   geographic%xmax = geographic%xmin+real(geographic%nx*(xproc)/nxprocs)*geographic%dx
                   geographic%xmin = geographic%xmin+real(geographic%nx*(xproc-1)/nxprocs)*geographic%dx
                   geographic%ymax = geographic%ymin+real(geographic%ny*(yproc)/nyprocs)*geographic%dy
                   geographic%ymin = geographic%ymin+real(geographic%ny*(yproc-1)/nyprocs)*geographic%dy
                   isprediction = .true.
                end if
                proc = proc + 1                                 
             endif
          enddo
       enddo

       ! fill up water points for unused processes
       proc = 0
       do yproc=1,nyprocs
          do xproc=1,nxprocs
             if (.not.landmask(xproc,yproc)) then
                if ((mpiid.eq.proc).and.(mpiid.lt.(numprocs-nlandproc))) then
                   geographic%xmax = geographic%xmin+real(geographic%nx*(xproc)/nxprocs)*geographic%dx
                   geographic%xmin = geographic%xmin+real(geographic%nx*(xproc-1)/nxprocs)*geographic%dx
                   geographic%ymax = geographic%ymin+real(geographic%ny*(yproc)/nyprocs)*geographic%dy
                   geographic%ymin = geographic%ymin+real(geographic%ny*(yproc-1)/nyprocs)*geographic%dy
                end if
                proc = proc + 1                 
             endif
          enddo
       enddo

       deallocate(landmask)

       ! generate grid boundaries for local tile
       geographic%nx = nint((geographic%xmax - geographic%xmin)/geographic%dx)
       geographic%ny = nint((geographic%ymax - geographic%ymin)/geographic%dy)

    endif ! **** DISTRIBUTED OR MULTISITE ****

    !--------------------------
    ! PFT/HGT/MODEL GRIDS
    !--------------------------

    ! output pft/hgt distribution
    if (averagepft) then
       npft_out = 1 
    else 
       npft_out = npft
    endif

    if (averagehgt) then
       nhgt_out = 1 
    else 
       nhgt_out = nhgt
    endif

    if (isprediction) then
       ! generate assimilation / model grid for this process
       allocate(xcoord(geographic%nx,geographic%ny))
       allocate(ycoord(geographic%nx,geographic%ny))
       allocate(lon(geographic%nx,geographic%ny))
       allocate(lat(geographic%nx,geographic%ny))
       allocate(area(geographic%nx,geographic%ny))
       do i=1,geographic%nx
          xcoord(i,:) = geographic%xmin + geographic%dx*(real(i)-0.5)
       end do
       do j=1,geographic%ny
          ycoord(:,j) = geographic%ymin + geographic%dy*(real(j)-0.5)
       end do

       call proj_inv(geographic%type,xcoord,ycoord,lon,lat,geographic%optargs,nodata)

       ! if whole region is outside of the geographical domain, do not perform prediction
       if (count((lon-nodata).gt.eps).eq.0) isprediction = .false.

       do j=1,geographic%ny
          do i=1,geographic%nx
       
             xbounds(:,1) = (/xcoord(i,j)-0.5*geographic%dx, xcoord(i,j)+0.5*geographic%dx, xcoord(i,j), xcoord(i,j)/)
             ybounds(:,1) = (/ycoord(i,j), ycoord(i,j), ycoord(i,j)-0.5*geographic%dy, ycoord(i,j)+0.5*geographic%dy/)
             
             call proj_inv(geographic%type,xbounds,ybounds,lonbounds,latbounds,geographic%optargs,nodata)
                
             if (((lonbounds(1,1)-nodata).gt.eps).and.&
                  ((lonbounds(2,1)-nodata).gt.eps).and.&
                  ((lonbounds(3,1)-nodata).gt.eps).and.&
                  ((lonbounds(4,1)-nodata).gt.eps).and.&
                  ((latbounds(1,1)-nodata).gt.eps).and.&
                  ((latbounds(2,1)-nodata).gt.eps).and.&
                  ((latbounds(3,1)-nodata).gt.eps).and.&
                  ((latbounds(4,1)-nodata).gt.eps)) then
                area(i,j) =  calc_distance(lonbounds(1,1),latbounds(1,1),lonbounds(2,1),latbounds(2,1)) * &
                     calc_distance(lonbounds(3,1),latbounds(3,1),lonbounds(4,1),latbounds(4,1))
             else
                area(i,j) = nodata
             end if
          end do
       end do

       ! generate output grid for this process
       geographic_out = geographic
       geographic_out%dx = geographic_out%dx/real(subgrid,kind=4)
       geographic_out%dy = geographic_out%dy/real(subgrid,kind=4)
       geographic_out%nx = geographic_out%nx*subgrid
       geographic_out%ny = geographic_out%ny*subgrid 
       allocate(xcoord_out(geographic_out%nx))
       allocate(ycoord_out(geographic_out%ny))
       do i=1,geographic_out%nx
          xcoord_out(i) = geographic_out%xmin + geographic_out%dx*(real(i)-0.5)
       end do
       do j=1,geographic_out%ny
          ycoord_out(j) = geographic_out%ymin + geographic_out%dy*(real(j)-0.5)
       end do

       ! generate x/y output grid offsets for this process with respect to full output grid
       xcoord_offset = nint((geographic_out%xmin - geographic_full%xmin)/geographic_full%dx) + 1
       ycoord_offset = nint((geographic_out%ymin - geographic_full%ymin)/geographic_full%dy) + 1

    end if

    if (regional) then

       ! regional analysis: read pft and hgt percentage distribution of output grid

       allocate(pft_pct_ave(npft))
       allocate(hgt_pct_ave(nhgt))

       if (mpiid.eq.0) write(*,'(A)') 'Reading and deriving subgrid-scale topography.'

       if (isprediction) then

          allocate(pft_pct_out(geographic_out%nx,geographic_out%ny,npft))
          allocate(hgt_pct_out(geographic_out%nx,geographic_out%ny,nhgt))
          allocate(hgt_out(geographic_out%nx,geographic_out%ny))
          allocate(validmask(geographic_out%nx,geographic_out%ny))
          allocate(hgt(geographic%nx,geographic%ny))

          ! read topography on model grid
          call read_topo(geographic%type,geographic%xmin,geographic%xmax,geographic%ymin,geographic%ymax, &
               geographic%dx,geographic%dy,geographic%optargs,nodata,topo=hgt)    

          ! read topography on output grid
          call read_topo(geographic_out%type,geographic_out%xmin,geographic_out%xmax,geographic_out%ymin,geographic_out%ymax, &
               geographic_out%dx,geographic_out%dy,geographic_out%optargs,nodata,topo=hgt_out)    

          ! derive minimum and maximum hgt from output topography
          validmask = (hgt_out-nodata).gt.eps
          if (count(validmask).ne.0) then
             hgt_ave = sum(hgt_out,mask=validmask)/real(count(validmask))
             hgt_stddev = (sum((hgt_out-hgt_ave)**2,mask=validmask)/real(count(validmask)))**0.5

             if (hgt_stddev.lt.0.1) then
                ! output grid too coarse or single point
                minhgt = max(hgt_ave - 500.,0.)
                maxhgt = min(hgt_ave + 500.,4000.)
             else
                minhgt = max(hgt_ave - 3. * hgt_stddev,0.)
                maxhgt = min(hgt_ave + 3. * hgt_stddev,4000.)
             endif
          else
             hgt_ave = nodata
             hgt_stddev = nodata
             minhgt = nodata
             maxhgt = nodata
          end if
          deallocate(validmask)

       else
          hgt_ave = nodata
          hgt_stddev = nodata
          minhgt = nodata
          maxhgt = nodata
       endif

#ifndef HIDE_MPI
       if ((nxprocs.gt.1).or.(nyprocs.gt.1)) then
          ! for single distributed domain, create a global set of elevation bins
          ! take minimum of the minimum, or maximum of the maximum, of all processes

          if (mpiid.eq.0) then
             if ((minhgt-nodata).le.eps) minhgt = 99999.
             if ((maxhgt-nodata).le.eps) maxhgt = -99999.

             do proc=2,numprocs
                call MPI_RECV(fval,1,MPI_REAL,proc-1,proc,MPI_COMM_WORLD,mpistatus,mpierr)
                if (((fval-nodata).gt.eps).and.(fval.lt.minhgt)) minhgt = fval
                call MPI_RECV(fval,1,MPI_REAL,proc-1,proc,MPI_COMM_WORLD,mpistatus,mpierr)
                if (((fval-nodata).gt.eps).and.(fval.gt.maxhgt)) maxhgt = fval
             enddo
          else
             call MPI_SEND(minhgt,1,MPI_REAL,0,mpiid+1,MPI_COMM_WORLD,mpierr)
             call MPI_SEND(maxhgt,1,MPI_REAL,0,mpiid+1,MPI_COMM_WORLD,mpierr)
          endif

          call MPI_BCAST(minhgt,1,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
          call MPI_BCAST(maxhgt,1,MPI_REAL,0,MPI_COMM_WORLD,mpierr)

       endif
#endif
       
       allocate(hgtlevels(nhgt))
       do i=1,nhgt
          hgtlevels(i) = (maxhgt-minhgt)/real(nhgt,kind=4) * (real(i,kind=4)-0.5) + minhgt
       enddo

       if (isprediction) then

          allocate(validmask(geographic_out%nx,geographic_out%ny))
          validmask = (hgt_out-nodata).gt.eps

          call read_topo(geographic_out%type,geographic_out%xmin,geographic_out%xmax,geographic_out%ymin,geographic_out%ymax, &
               geographic_out%dx,geographic_out%dy,geographic_out%optargs,nodata,topohist=hgt_pct_out,topolevels=hgtlevels)    

          ! generate % out of histogram counts
          do j=1,geographic_out%ny
             do i=1,geographic_out%nx
                if ((hgt_pct_out(i,j,1)-nodata).gt.eps) then
                   hgt_pct_out(i,j,:) = hgt_pct_out(i,j,:) / max(sum(hgt_pct_out(i,j,:)),1.0)
                end if
             end do
          end do

          ! average hgt's over output grid
          if (count(validmask).ne.0) then
             do h = 1, nhgt
                hgt_pct_ave(h) = sum(hgt_pct_out(:,:,h),mask=validmask)/real(count(validmask),kind=4)
             enddo
          else
             hgt_pct_ave = nodata
          end if

          if (mpiid.eq.0) write(*,'(A)') 'Reading and deriving subgrid-scale pft distribution.'

          call read_pft(geographic_out%type,geographic_out%xmin,geographic_out%xmax,geographic_out%ymin,geographic_out%ymax, &
               geographic_out%dx,geographic_out%dy,geographic_out%optargs,nodata,pft3d=pft_pct_out)    

          ! average pft's over output grid
          if (count(validmask).ne.0) then
             do p = 1, npft
                pft_pct_ave(p) = sum(pft_pct_out(:,:,p),mask=validmask)/real(count(validmask),kind=4)
             enddo
          else
             pft_pct_ave = nodata
          end if
          deallocate(validmask)

          ! regional analysis, select pft's and hgt's which need to be assimilated

          ! select all pft's that cover more than x% of (100%-minsumpft) averaged over entire grid
          
          if (analysis) then
             minselpft = 0.1 * (1.-minsumpft) * (1.-pft_pct_ave(npft))
          else
             minselpft = 0.
          endif

          if (analysis.or.averagepft) then
             ! all, but no water pft's.
             minpft = 1
             maxpft = npft-1

             ! all, but no bare soil and water pft's.
             ! That may not be wise since MODIS can see soil so MODIS pixel integrate soil reflectance.
             ! minpft = 2
             ! maxpft = npft-1

             ! only tall vegetation and shrub pft's.
             ! minpft = 1
             ! maxpft = 12
             
             ! only tall vegetation pft's.
!             minpft = 1
!             maxpft = 9
          else
             ! all pft's, only if pft's are output separately and in pure prediction mode
             minpft = 1
             maxpft = npft
          endif

          nselect_pft = count(pft_pct_ave(minpft:maxpft).ge.minselpft)

          allocate(select_pft(nselect_pft))
          allocate(pft_is_selected(npft))

          pft_is_selected(:) = .false.
          r=1
          do p = minpft,maxpft
             if (pft_pct_ave(p).ge.minselpft) then
                pft_is_selected(p) = .true.
                select_pft(r) = p
                r=r+1
             else
                pft_is_selected(p) = .false.
             endif
          enddo

          if (analysis) then
             minselhgt = 0.1 * (1.-minsumhgt)
          else
             minselhgt = 0.0
          endif

          ! count all hgt's that substantially contribute to this grid
          nselect_hgt = count(hgt_pct_ave(1:nhgt).ge.minselhgt)

          allocate(select_hgt(nselect_hgt))
          allocate(hgt_is_selected(nhgt))

          hgt_is_selected(:) = .false.
          r=1
          do h = 1,nhgt
             if (hgt_pct_ave(h).ge.minselhgt) then
                hgt_is_selected(h) = .true.
                select_hgt(r) = h
                r=r+1
             else
                hgt_is_selected(h) = .false.
             endif
          enddo

          ! do not perform prediction on water only areas
          if ((nselect_pft.eq.0).or.(nselect_hgt.eq.0)) then
             isprediction = .false.
          else
             isprediction = .true.
          endif

       else
          ! no prediction on this node
          nselect_hgt = 1
          nselect_pft = 1
          allocate(select_pft(nselect_pft))
          allocate(select_hgt(nselect_hgt))
          select_pft(:) = npft
          select_hgt(:) = nhgt
          pft_pct_ave(:) = nodata
          hgt_pct_ave(:) = nodata
          allocate(pft_pct_out(geographic_out%nx,geographic_out%ny,npft))
          allocate(hgt_pct_out(geographic_out%nx,geographic_out%ny,nhgt))
          allocate(hgt_out(geographic_out%nx,geographic_out%ny))
          pft_pct_out(:,:,:) = nodata
          hgt_pct_out(:,:,:) = nodata
          hgt_out(:,:) = nodata
       endif

       ! generate averaged pft distribution over all prediction processes in distributed mode
       if ((nxprocs.gt.1).or.(nyprocs.gt.1)) then

          if (mpiid.eq.0) then
#ifndef HIDE_MPI
             if (sum(pft_pct_ave).gt.eps) then
                n = 1
             else
                n = 0
                pft_pct_ave = 0.0
             end if
             allocate(pft_tmp(npft))
             do proc=2,numprocs
                call MPI_RECV(pft_tmp,npft,MPI_REAL,proc-1,proc,MPI_COMM_WORLD,mpistatus,mpierr)
                if (sum(pft_tmp).gt.eps) then
                   pft_pct_ave = pft_pct_ave + pft_tmp
                   n = n + 1
                end if
             enddo
             deallocate(pft_tmp)
             if (n.eq.0) then
                pft_pct_ave = nodata
             else
                pft_pct_ave = pft_pct_ave / real(n,kind=4)
             end if
#endif
          else
#ifndef HIDE_MPI
             call MPI_SEND(pft_pct_ave,npft,MPI_REAL,0,mpiid+1,MPI_COMM_WORLD,mpierr)
#endif 
          endif

#ifndef HIDE_MPI
          ! send this pft distribution to all processes
          call MPI_BCAST(pft_pct_ave,npft,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
#endif 

       endif

    else

       ! local analysis: single pft class, single elevation class

       allocate(pft_pct_ave(npft))
       allocate(hgt_pct_ave(nhgt))

       if (isprediction) then

          allocate(pft_pct_out(geographic_out%nx,geographic_out%ny,npft))
          allocate(hgt_pct_out(geographic_out%nx,geographic_out%ny,nhgt))
          allocate(hgt_out(geographic_out%nx,geographic_out%ny))
          allocate(hgt(geographic%nx,geographic%ny))
          hgt_ave = 0.

          minhgt = hgt_ave - 500.0
          maxhgt = hgt_ave + 500.0

          hgt_stddev = 0.0
          pft_pct_out(:,:,:) = 1.0
          hgt_pct_out(:,:,:) = 1.0
          hgt_out(:,:) = hgt_ave
          pft_pct_ave(:) = 1.0
          hgt_pct_ave(:) = 1.0
          allocate(hgtlevels(nhgt))
          hgtlevels(:) = hgt_ave
          hgt(:,:) = hgt_ave

          ! select pft's and hgt's for assimilation
          ! local mode: not pft-specific and not elevation-specific
          nselect_pft = 1
          nselect_hgt = 1
          allocate(select_pft(nselect_pft))
          allocate(select_hgt(nselect_hgt))
          allocate(pft_is_selected(1))
          allocate(hgt_is_selected(1))
          select_pft(1) = 1
          select_hgt(1) = 1
          pft_is_selected(1) = .true.
          hgt_is_selected(1) = .true.
       else
          ! no prediction on this node
!          nselect_hgt = 1
!          nselect_pft = 1
!          allocate(select_pft(nselect_pft))
!          allocate(select_hgt(nselect_hgt))
!          select_pft(:) = npft
!          select_hgt(:) = nhgt
!          pft_pct_ave(:) = 0.
!          hgt_pct_ave(:) = 0.
!          allocate(pft_pct_out(geographic%nx_out,geographic%ny_out,npft))
!          allocate(hgt_pct_out(geographic%nx_out,geographic%ny_out,nhgt))
!          pft_pct_out(:,:,:) = 0.
!          hgt_pct_out(:,:,:) = 0.
       endif
    endif

    !--------------------------
    ! COUNTERS + TIMINGS
    !--------------------------

    ! allocate observation counter
    allocate(cumobs(npft))
    cumobs(:) = 0.0

    !--------------------------
    ! SELECTION OF VARIABLES
    !--------------------------

    ! total number of points in assimilation vector (including ensembles in domain area geographic%ny*geographic%nx)
    nptstate = geographic%nx*geographic%ny*nselect_pft*nselect_hgt
    nptparam = npft
    nptforce = geographic%nx*geographic%ny*nselect_hgt
    nptobser = geographic%nx*geographic%ny

    ! read state/parameter/force definition files

    do f=1,4

       select case (f)
       case (1)  
          filename = trim(namelistdir)//'states.dat'
       case (2)  
          filename = trim(namelistdir)//'parameters.dat'
       case (3)  
          filename = trim(namelistdir)//'forcings.dat'
       case (4)  
          filename = trim(namelistdir)//'observations.dat'
       end select

       if (mpiid.eq.0) write(*,'(A,A)') 'Parsing definition file: ',trim(filename)

       ! read site case file
       call read_csv_table(filename,table,',',2)
       ncol = size(table,1)
       nlin = size(table,2)

       ! allocate arrays for variables that are going to be read
       select case (f)
       case (1)
          nstate = nlin
          allocate(state_var(nstate))
          allocate(state_names(nstate))
          allocate(state_longnames(nstate))
          allocate(state_units(nstate))
          allocate(state_min(nstate))
          allocate(state_max(nstate))
          allocate(state_scale(nstate))
          allocate(state_offset(nstate))
          allocate(state_type(nstate))
          allocate(state(nrens,geographic%nx,geographic%ny,nselect_pft,nselect_hgt,nstate))
          if (simulator) allocate(state_sim(1,geographic%nx,geographic%ny,nselect_pft,nselect_hgt,nstate))
          allocate(sel(nstate))
          allocate(out(nstate))
       case (2)
          nparam = nlin
          allocate(param_var(npft,nparam))
          allocate(param_names(nparam))
          allocate(param_longnames(nparam))
          allocate(param_units(nparam))
          allocate(param_min(nparam))
          allocate(param_max(nparam))
          allocate(param_scale(nparam))
          allocate(param_offset(nparam))
          allocate(param_type(nparam))
          allocate(param(nrens,npft,nparam))
          if (simulator) allocate(param_sim(1,npft,nparam))
          allocate(sel(nparam))
          allocate(out(nparam))
       case (3)
          nforce = nlin
          allocate(force_var(nforce))
          allocate(force_names(nforce))
          allocate(force_longnames(nforce))
          allocate(force_units(nforce))
          allocate(force_min(nforce))
          allocate(force_max(nforce))
          allocate(force_scale(nforce))
          allocate(force_offset(nforce))
          allocate(force_type(nforce))
          allocate(force(nrens,geographic%nx,geographic%ny,nselect_hgt,nforce))
          if (simulator) allocate(force_sim(1,geographic%nx,geographic%ny,nselect_hgt,nforce))
          allocate(sel(nforce))
          allocate(out(nforce))
       case (4)
          nobser = nlin
          allocate(obser_var(nobser))
          allocate(obser_names(nobser))
          allocate(obser_longnames(nobser))
          allocate(obser_units(nobser))
          allocate(obser_min(nobser))
          allocate(obser_max(nobser))
          allocate(obser_scale(nobser))
          allocate(obser_offset(nobser))
          allocate(obser_type(nobser))
          allocate(obser(2,geographic%nx,geographic%ny,nobser))
          allocate(sel(nobser))
          allocate(out(nobser))
       end select

       do v = 1, nlin

          ! assign variables, initial conditions, names and units
          select case (f)
          case (1)
             read(table(2,v),'(I4)') sel(v)
             read(table(3,v),'(I4)') out(v)
             read(table(4,v),'(g14.6)') fval
             if (isprediction) then
                state(:,:,:,:,:,v)=fval
                if (simulator) state_sim(:,:,:,:,:,v)=fval
             else
                state(:,:,:,:,:,v)=nodata
                if (simulator) state_sim(:,:,:,:,:,v)=nodata
             endif
             read(table(5,v),'(g14.6)') state_var(v)
             read(table(6,v),'(g14.6)') state_min(v)
             read(table(7,v),'(g14.6)') state_max(v)
             state_names(v) = table(8,v)
             state_longnames(v) = table(9,v)
             state_units(v) = table(10,v)
             read(table(11,v),'(g14.6)') state_scale(v)
             read(table(12,v),'(g14.6)') state_offset(v)
             state_type(v) = table(13,v)
          case (2)
             read(table(2,v),'(I4)') sel(v)
             read(table(3,v),'(I4)') out(v)
             read(table(4,v),'(g14.6)') fval
             param(:,:,v)=fval
             if (simulator) param_sim(:,:,v)=fval
             read(table(5,v),'(g14.6)') fval
             if (sel(v).eq.1) then
                param_var(:,v)=fval
             else
                param_var(:,v)=0.
             endif
             read(table(6,v),'(g14.6)') param_min(v) 
             read(table(7,v),'(g14.6)') param_max(v) 
             param_names(v) = table(8,v)
             param_longnames(v) = table(9,v)
             param_units(v) = table(10,v)
             read(table(11,v),'(g14.6)') param_scale(v)
             read(table(12,v),'(g14.6)') param_offset(v)
             param_type(v) = table(13,v)
          case (3)
             read(table(2,v),'(I4)') sel(v)
             read(table(3,v),'(I4)') out(v)
             read(table(4,v),'(g14.6)') fval
             if (isprediction) then
                force(:,:,:,:,v)=fval
                if (simulator) force_sim(:,:,:,:,v)=fval
             else
                force(:,:,:,:,v)=nodata
                if (simulator) force_sim(:,:,:,:,v)=nodata
             endif
             read(table(5,v),'(g14.6)') force_var(v)
             read(table(6,v),'(g14.6)') force_min(v)
             read(table(7,v),'(g14.6)') force_max(v)
             force_names(v) = table(8,v)
             force_longnames(v) = table(9,v)
             force_units(v) = table(10,v)
             read(table(11,v),'(g14.6)') force_scale(v)
             read(table(12,v),'(g14.6)') force_offset(v)
             force_type(v) = table(13,v)
          case (4)
             read(table(2,v),'(I4)') sel(v)
             read(table(3,v),'(I4)') out(v)
             read(table(4,v),'(g14.6)') fval
             if (isprediction) then
                obser(1,:,:,v)=fval
             else
                obser(1,:,:,v)=nodata
             endif
             read(table(5,v),'(g14.6)') obser_var(v)
             if (isprediction) then
                obser(2,:,:,v)=obser_var(v)
             else
                obser(2,:,:,v)=nodata
             endif
             read(table(6,v),'(g14.6)') obser_min(v)
             read(table(7,v),'(g14.6)') obser_max(v)
             obser_names(v) = table(8,v)
             obser_longnames(v) = table(9,v)
             obser_units(v) = table(10,v)
             read(table(11,v),'(g14.6)') obser_scale(v)
             read(table(12,v),'(g14.6)') obser_offset(v)
             obser_type(v) = table(13,v)
          end select

       end do

       deallocate(table)

       ! in distributed mode: do not write parameters since they would vary
       ! spatially depending on which pft distribution is present by area
       if (((nxprocs.gt.1).or.(nyprocs.gt.1)).and.(f.eq.2)) then
          out(:) = 0
       endif

       ! check which variables are selected for assimilation/perturbing/saving
       select case (f)
       case (1)
          s = 0
          o = 0
          nselect_state = sum(sel)
          noutput_state = sum(out)
          if (nselect_state.gt.0) allocate(select_state(nselect_state))
          if (noutput_state.gt.0) allocate(output_state(noutput_state))
          do v = 1, nlin
             if (sel(v)==1) then
                s = s + 1
                select_state(s) = v
             endif
             if (out(v)==1) then
                o = o + 1
                output_state(o) = v
             endif
          enddo
       case (2)
          s = 0
          o = 0
          nselect_param = sum(sel)
          noutput_param = sum(out)
          if (nselect_param.gt.0) allocate(select_param(nselect_param))
          if (noutput_param.gt.0) allocate(output_param(noutput_param))
          do v = 1, nlin
             if (sel(v)==1) then
                s = s + 1
                select_param(s) = v
             endif
             if (out(v)==1) then
                o = o + 1
                output_param(o) = v
             endif
          enddo
       case (3)
          s = 0
          o = 0
          nselect_force = sum(sel)
          noutput_force = sum(out)
          if (nselect_force.gt.0) allocate(select_force(nselect_force))
          if (noutput_force.gt.0) allocate(output_force(noutput_force))
          do v = 1, nlin
             if (sel(v)==1) then
                s = s + 1
                select_force(s) = v
             endif
             if (out(v)==1) then
                o = o + 1
                output_force(o) = v
             endif
          enddo
       case (4)
          s = 0
          o = 0
          nselect_obser = sum(sel)
          noutput_obser = sum(out)
          if (nselect_obser.gt.0) allocate(select_obser(nselect_obser))
          if (noutput_obser.gt.0) allocate(output_obser(noutput_obser))
          do v = 1, nlin
             if (sel(v)==1) then
                s = s + 1
                select_obser(s) = v
             endif
             if (out(v)==1) then
                o = o + 1
                output_obser(o) = v
             endif
          enddo
       end select

       deallocate(sel,out)

    end do

    ! assign state / parameter / force / observation indices
    do n=1,nstate
       select case (trim(state_names(n)))
       case ("FPAR") 
          stateidx%fpar = n
       case ("LAI") 
          stateidx%lai = n
       case ("GSI") 
          stateidx%gsi = n
       case ("TEMP_FAC") 
          stateidx%temp_fac = n
       case ("FORCE_FAC") 
          stateidx%force_fac = n
       case ("CHILL_FAC") 
          stateidx%chill_fac = n
       case ("MOIST_FAC") 
          stateidx%moist_fac = n
       case ("LIGHT_FAC") 
          stateidx%light_fac = n
       case ("TMIN_AVE") 
          stateidx%tmin_ave = n
       case ("VPD_AVE") 
          stateidx%vpd_ave = n
       case ("PHOTO_AVE") 
          stateidx%photo_ave = n
       case ("RG_AVE") 
          stateidx%rg_ave = n
       case ("RAIN_AVE") 
          stateidx%rain_ave = n
       case ("CHILL_AVE") 
          stateidx%chill_ave = n
       case default
          write(*,'(A,A)') "Unknown state in state namelist: ",trim(state_names(n))
          call model_exit(-1)
       end select
    end do

    do p=1,nparam
       select case (trim(param_names(p)))
       case ("TMIN_MIN") 
          paramidx%tmin_min = p
       case ("TMIN_MAX") 
          paramidx%tmin_max = p
       case ("VPD_MIN") 
          paramidx%vpd_min = p
       case ("VPD_MAX") 
          paramidx%vpd_max = p
       case ("PHOTO_MIN") 
          paramidx%photo_min = p
       case ("PHOTO_MAX") 
          paramidx%photo_max = p
       case ("RG_MIN") 
          paramidx%rg_min = p
       case ("RG_MAX") 
          paramidx%rg_max = p
       case ("RAIN_MIN") 
          paramidx%rain_min = p
       case ("RAIN_MAX") 
          paramidx%rain_max = p
       case ("CHILL_MIN") 
          paramidx%chill_min = p
       case ("CHILL_MAX") 
          paramidx%chill_max = p
       case ("LAI_MIN") 
          paramidx%lai_min = p
       case ("LAI_MAX") 
          paramidx%lai_max = p
       case ("GROWTHRATE") 
          paramidx%growthrate = p
       case ("DECAYRATE") 
          paramidx%decayrate = p
       case ("LAI_SAT") 
          paramidx%lai_sat = p
       case ("FPAR_SAT") 
          paramidx%fpar_sat = p
       case ("FVCOVER") 
          paramidx%fvcover = p
       case ("TMIN_TAVE") 
          paramidx%tmin_tave = p
       case ("VPD_TAVE") 
          paramidx%vpd_tave = p
       case ("PHOTO_TAVE") 
          paramidx%photo_tave = p
       case ("RG_TAVE") 
          paramidx%rg_tave = p
       case ("RAIN_TAVE") 
          paramidx%rain_tave = p
       case ("CHILL_TAVE") 
          paramidx%chill_tave = p
       case default
          write(*,'(A,A)') "Unknown parameter in parameter namelist: ",trim(param_names(p))
          call model_exit(-1)
       end select
    end do

    do f=1,nforce
       select case (trim(force_names(f)))
       case ("TMIN") 
          forceidx%tmin = f
       case ("VPD") 
          forceidx%vpd = f
       case ("PHOTO") 
          forceidx%photo = f
       case ("RG") 
          forceidx%rg = f
       case ("RAIN") 
          forceidx%rain = f
       case ("TMEAN") 
          forceidx%tmean = f
       case default
          write(*,'(A,A)') "Unknown forcing in forcing namelist: ",trim(force_names(f))
          call model_exit(-1)
       end select
    end do

    do o=1,nobser
       select case (trim(obser_names(o)))
       case ("FPAR") 
          obseridx%fpar = o
       case ("LAI") 
          obseridx%lai = o
       case default
          write(*,'(A,A)') "Unknown observation in observation namelist: ",trim(obser_names(o))
          call model_exit(-1)
       end select
    end do

    if (verbose) then

       ! print diagnostics
       write(*,*)
       write(*,'(A,I6)') 'Process ID:               ',mpiid
       write(*,'(A,L6)') 'Prediction:               ',isprediction
       write(*,'(A,L6)') 'Assimilation:             ',analysis
       write(*,'(A,I6)') 'Ensemble members:         ',nrens
       write(*,*)

       if (isprediction) then
          write(*,'(A,A,A,A)') 'Projection Type: ',trim(geographic%type),' / Parameters: ',trim(geographic%optargs)
          write(*,'(A,F8.3,A,F8.3,A,F8.3,A,F8.3)') 'Coordinate bounds: ',&
               geographic%xmin,' to ',geographic%xmax,'; ',geographic%ymin,' to ',geographic%ymax
          write(*,'(A,4I6)') 'Model  grid: nx,ny,npft,nhgt: ',geographic%nx,geographic%ny,npft,nhgt
          write(*,'(A,4I6)') 'Output grid: nx,ny,npft,nhgt: ',geographic_out%nx,geographic_out%ny,npft_out,nhgt_out

          write(*,*)
          write(*,'(A,I4)') 'Subgrid-scale PFT distribution on output grid: '
          do p=1,npft
             write(*,'(A,I4,A,A,A,F6.1,A,L1)') 'PFT ',p,' type: ',pftnames(p),': ',pft_pct_ave(p)*100., &
                  ' % of area. Selected: ',pft_is_selected(p)
          enddo
          write(*,'(A,I6)') 'Subgrid-scale PFT selected #: ',nselect_pft
          write(*,*)
          write(*,'(A,F6.1,A,F6.1,A)') 'Mean HGT of output grid: ',hgt_ave,' m, Stddev: ',hgt_stddev,' m'
          write(*,'(A,I4)') 'Subgrid-scale HGT distribution (nhgt=): ',nhgt
          do h=1,nhgt
             write(*,'(A,I4,A,F6.1,A,F6.1,A,L1)') 'level ',h,' height: ',hgtlevels(h), ' m, ',hgt_pct_ave(h)*100.,&
                  ' % of area. Selected: ',hgt_is_selected(h)
          enddo
          write(*,'(A,I6)') 'Subgrid-scale HGT selected #: ',nselect_hgt
          if (hgtscreen) then
             write(*,'(A,F6.1,A,F6.1,A)') 'Min/Max HGT used to screen observations: ',minhgt,' m / ',maxhgt,' m'
          else
             write(*,'(A)') 'No HGT screening of observations applied'
          endif
          write(*,*)
       endif

       if (analysis) then
          write(*,'(A,I6)') 'Assimilated Parameters:   ',nselect_param * nptparam
          write(*,'(A,I6)') 'Assimilated States:       ',nselect_state * nptstate
       endif

    endif

    if ((analysis).and.(nselect_obser.eq.0)) then
       write(*,'(A)') 'Please select at least one observation stream'
       write(*,'(A)') 'if running in data assimilation (analysis) mode'
       call model_exit(-7400)
    endif

    !--------------------------
    ! MODEL ARRAY ALLOCATION
    !--------------------------

    allocate(ensindex(nrens))
    allocate(pftindex(npft))

    do r = 1, nrens
       ensindex(r) = r
    enddo

    do p = 1, npft
       pftindex(p) = p
    enddo


    ! adjust parameters for some pft's in regional mode
    ! choosing the right start parameter set enables a more
    ! precise parameter adjustment

    if (regional) then

       do p = 1,npft

          ! reset bare soil pft structural parameters and variability to 0
          if (trim(pftnames(p)).eq.'bar_all') then

             if (verbose) write(*,'(A,A)') &
                  'Adjusting parameters of pft ', trim(pftnames(p))

             param(:,p,paramidx%lai_min) = 0. 
             param(:,p,paramidx%lai_max) = 0. 
             param_var(p,paramidx%lai_min) = 0.
             param_var(p,paramidx%lai_max) = 0.

          endif

          ! set temporal averaging parameters of evergreen species to beyond a year
          ! in remove seasonality in predicted states
!          if ((trim(pftnames(p)).eq.'enf_bor').or.(trim(pftnames(p)).eq.'enf_tem').or. &
!               (trim(pftnames(p)).eq.'ebf_tro').or.(trim(pftnames(p)).eq.'ebf_tem').or. &
!               (trim(pftnames(p)).eq.'ebs_all')) then

!             if (verbose) write(*,'(A,A)') &
!                  'Adjusting parameters of pft ', trim(pftnames(p))

!             param(:,p,paramidx%tmin_tave) = 1000.0
!             param(:,p,paramidx%vpd_tave) = 1000.0
!             param(:,p,paramidx%photo_tave) = 1000.0
!             param_var(p,paramidx%tmin_tave) = 0.01
!             param_var(p,paramidx%vpd_tave) = 0.01
!             param_var(p,paramidx%photo_tave) = 0.01

!          end if

!          if (trim(pftnames(p)).eq.'ebf_tro') then

!             if (verbose) write(*,'(A,A)') &
!                  'Adjusting parameters of pft ', trim(pftnames(p))

!             param(:,p,paramidx%lai_min) = 5.0  
!             param(:,p,paramidx%lai_max) = 8.0  

!          endif

!          if (trim(pftnames(p)).eq.'enf_tem') then

!             if (verbose) write(*,'(A,A)') &
!                  'Adjusting parameters of pft ', trim(pftnames(p))

!             param(:,p,paramidx%lai_min) = 2.0 
!             param(:,p,paramidx%lai_max) = 7.0 

!          endif

!          if (trim(pftnames(p)).eq.'enf_bor') then

!             if (verbose) write(*,'(A,A)') &
!                  'Adjusting parameters of pft ', trim(pftnames(p))

!             param(:,p,paramidx%lai_min) = 2.0 
!             param(:,p,paramidx%lai_max) = 7.0 

!          endif

!          if (trim(pftnames(p)).eq.'dbf_tem') then

!             if (verbose) write(*,'(A,A)') &
!                  'Adjusting parameters of pft ', trim(pftnames(p))

!             param(:,p,paramidx%tmin_min) = 271.0
!             param(:,p,paramidx%tmin_max) = 278.0
!             param_var(p,paramidx%tmin_min) = 1.0
!             param_var(p,paramidx%tmin_max) = 1.0

!             param(:,p,paramidx%photo_min) = 10.0
!             param(:,p,paramidx%photo_max) = 11.0
!             param_var(p,paramidx%photo_min) = 1.0e-2
!             param_var(p,paramidx%photo_max) = 1.0e-2

!             param(:,p,paramidx%lai_max) = 5.5

!          endif

!          if ((trim(pftnames(p)).eq.'enf_bor').or.(trim(pftnames(p)).eq.'dnf_bor').or. &
!               (trim(pftnames(p)).eq.'dbf_bor')) then

!             if (verbose) write(*,'(A,A)') &
!                  'Adjusting parameters of pft ', trim(pftnames(p))

!             param(:,p,paramidx%photo_min) = param(:,p,paramidx%photo_min) + 1.
!             param(:,p,paramidx%photo_max) = param(:,p,paramidx%photo_max) + 3.

!          endif

!          if (trim(pftnames(p)).eq.'dbs_bor') then

!             if (verbose) write(*,'(A,A)') &
!                  'Adjusting parameters of pft ', trim(pftnames(p))

!             param(:,p,paramidx%photo_min) = param(:,p,paramidx%photo_min) + 1.
!             param(:,p,paramidx%photo_max) = param(:,p,paramidx%photo_max) + 3.

!             param(:,p,paramidx%lai_max) = 2.0 

!          endif

!          if (trim(pftnames(p)).eq.'c3g_arc') then

!             if (verbose) write(*,'(A,A)') &
!                  'Adjusting parameters of pft ', trim(pftnames(p))

!             param(:,p,paramidx%photo_min) = param(:,p,paramidx%photo_min) + 1.
!             param(:,p,paramidx%photo_max) = param(:,p,paramidx%photo_max) + 3.

!             param(:,p,paramidx%lai_max) = 2.0 

!          endif

!          if (trim(pftnames(p)).eq.'c3g_nar') then

!             if (verbose) write(*,'(A,A)') &
!                  'Adjusting parameters of pft ', trim(pftnames(p))

!             param(:,p,paramidx%lai_max) = 3.0

!          endif

!          if (trim(pftnames(p)).eq.'c4g_all') then

!             if (verbose) write(*,'(A,A)') &
!                  'Adjusting parameters of pft ', trim(pftnames(p))

!             param(:,p,paramidx%lai_max) = 4.0 

!          endif

          if (trim(pftnames(p)).eq.'wat_all') then

             if (verbose) write(*,'(A,A)') &
                  'Adjusting parameters of pft ', trim(pftnames(p))

             param(:,p,paramidx%lai_min) = 0. 
             param(:,p,paramidx%lai_max) = 0. 
             param_var(p,paramidx%lai_min) = 0.
             param_var(p,paramidx%lai_max) = 0.

          endif

       end do

    endif

    ! wait with model time stepping until all processes have caught up
#ifndef HIDE_MPI
!    write(*,'(A,I6)') "MPIID: ",mpiid
    call MPI_Barrier(MPI_COMM_WORLD,mpierr)
#endif

    if (mpiid.eq.0) then
       write(*,'(A)') 'Finalized Initialization. Proceeding to time stepping.'
       write(*,*)
    endif

  end subroutine init_driver

end module init_module
