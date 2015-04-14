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

module exit_module

  private

  public :: handle_err
  public :: handle_error
  public :: model_exit

contains

  subroutine handle_err( status, routine, number )
    ! NetCDF error handling with error reporting

    use netcdf
    use typeSizes

    implicit none

    ! arguments
    integer, intent(in) :: status                       ! error status
    character(len=*), intent(in), optional :: routine   ! routine where err occurred
    integer, intent(in), optional :: number             ! error "line" number

    ! used to pinpoint which line caused the error
    if(status/=nf90_noerr) then
       write(*,'(A,I6,1X,A)') 'error', status, trim(nf90_strerror(status))
       if ( present(routine) ) write(*,'(A)')  trim(routine)
       if ( present(number) ) write(*,'(I8)')  number
    endif

    call model_exit(status)

  end subroutine handle_err

  subroutine handle_error(error, message)
    ! this subroutine checks for a true  error status and prints the 
    ! error message, handles the error with correct termination of
    ! parallel code.

    ! DO NOT CALL THIS ROUTINE ASYNCHRONOUSLY:
    ! make sure this code gets called for error evaluation from all processes

#ifndef HIDE_MPI
    use mpi
#endif

    implicit none

    ! arguments
    logical, intent(in) :: error
    character(len=*), intent(in) :: message

    ! local
    integer :: status
    logical :: allerror

    ! bring all processes to same error level.
#ifndef HIDE_MPI
    integer :: mpierr
    integer :: mpiid
    call MPI_COMM_RANK(MPI_COMM_WORLD,mpiid,mpierr)
    call MPI_ALLREDUCE(error, allerror, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, mpierr)
#else
    allerror = error
#endif     

    if (allerror) then
       
       if (error) then
          status = -2
#ifndef HIDE_MPI
          write(*,'(A,I4,A,A)') "MPI Process ",mpiid," with Error: ",trim(message)
#else
          write(*,'(A,A)') "Error: ",trim(message)
#endif     
       else
          status = -1
       end if

       call model_exit(status)

    end if

  end subroutine handle_error

  subroutine model_exit(status)
    ! this subroutine handles the exit
    ! terminates MPI processes gracefully

    ! DO NOT CALL THIS ROUTINE ASYNCHRONOUSLY:
    ! make sure this code gets called for exit from all processes
    ! alternatively, call handle_error() from all processes, which synchronizes
    ! error handling, and exits all processes, if one of the processes signals an error

    use file_module

#ifndef HIDE_MPI
    use mpi
#endif

    implicit none

    ! arguments
    integer, intent(in) :: status

    ! local
    integer :: unit
    integer :: mpiid
!    character(len=100) :: c

#ifndef HIDE_MPI
    integer :: mpierr

    call MPI_COMM_RANK(MPI_COMM_WORLD,mpiid,mpierr)
#else
    mpiid = 0
#endif

    if (status.ne.0) then
       ! print error status
#ifndef HIDE_MPI
       write(*,'(A,I4,A,I6)') 'Terminating MPI process ',mpiid,' with status: ',status
#else
       write(*,'(A,I6)') 'Terminating with status: ',status
#endif
    else
       if (mpiid.eq.0) then
          write(*,'(A)') 'Successfully exiting processing'

!          call getcwd(c)
!          write(*,'(A)') trim(c)

          ! write success flag to file
          unit = getunit()
          open(unit=unit, file='model_successful', form='formatted',action='write')
          write(unit,'(A)') 'OK'
          close(unit)
       endif

    endif

#ifndef HIDE_MPI
    ! clean up MPI if needed
    call mpi_finalize(mpierr)
#endif

    stop

  end subroutine model_exit

end module exit_module
