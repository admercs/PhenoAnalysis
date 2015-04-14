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

module file_module

  private

  public :: find_file
  public :: getunit
  public :: read_csv_table

  ! interface to C functions
  interface
     integer(c_int) function match_file_init (length,string) bind(c)
       use iso_c_binding
       implicit none
       integer(c_int) :: length
       character(c_char), dimension(*) :: string
     end function match_file_init
  end interface

  interface
     subroutine match_file_next (length,string) bind(c)
       use iso_c_binding
       implicit none
       integer(c_int) :: length
       character(c_char), dimension(*) :: string
     end subroutine match_file_next
  end interface

  interface
     subroutine match_file_exit () bind(c)
       use iso_c_binding
       implicit none
     end subroutine match_file_exit
  end interface

contains

  function find_file(ftemp) result(found)

    ! returns actual file name matched by ftemp pattern
    ! only returns the FIRST name found if many are found

    ! DEPENDENCIES: glob.c

    implicit none

    ! arguments
    character(len=*), intent(inout) :: ftemp

    ! result
    logical :: found

    ! local
    integer :: nfile, nf, length

    ! supply 0-trimmed length of string 
    length = len_trim(ftemp)

    ! set up file search
    nfile = match_file_init(length,ftemp)

    if (nfile.gt.0) then
       found = .true.

       ! supply full length of string to hold matched file
       length = len(ftemp)

       do nf = 1,nfile
          ! pick last file in list only
          call match_file_next(length,ftemp)
       enddo
    else
       ! initialize return value
       found = .false.
    endif

    call match_file_exit()

    return

  end function find_file

  integer function getunit()
 
    ! searches for and returns the next free file unit

    implicit none

    integer :: unit
    logical :: opened

    unit=10
    inquire(unit,opened=opened)
    DO WHILE (opened)
       unit=unit+1
       inquire(unit,opened=opened)
    END DO

    getunit = unit

    return

  end function getunit

  subroutine read_csv_table(filename,table,separator,headerlength)

    implicit none

    ! arguments
    character(*), intent(in) :: filename                   ! name of the file
    character(*), intent(out), allocatable :: table(:,:)   ! table entries (string array ncol x nlin)
    character(1), intent(in) :: separator                  ! separator (e.g. "," or ";"). do NOT use white space
    integer, intent(in) :: headerlength                    ! number of header lines

    ! local variables
    logical                        :: exist
    integer                        :: unit
    integer                        :: h, l, c, pos1, pos2, nchars, status
    integer                        :: nlin
    integer                        :: ncol 
    integer                        :: ncol_temp
    character(len=1000)            :: strtemp
    character(len=1), parameter    :: tab = char(9)

    inquire(file=filename, exist=exist)
    if (.not. exist) then
       write(*,'(A,A)') "Table file does not exist: ",filename
       stop
    end if

    unit = getunit()

    !! get number of table lines and columns
    l = 0
    open(unit=unit,file=trim(filename),form='formatted',action='read',status='old')
    do h=1,headerlength
       read(unit,'(A)',ADVANCE="no", SIZE=nchars, IOSTAT=status) strtemp ! read header
       l = l + 1
    end do

    nlin = 0
    ncol = 0
    do while (status.ne.-1)
       read(unit,'(A)',ADVANCE="no", SIZE=nchars, IOSTAT=status) strtemp ! read data
       l = l + 1
       if (status.ne.-1) then
          nlin = nlin + 1
          ncol_temp = 0
          do c=1,nchars
             if (strtemp(c:c).eq.separator) ncol_temp = ncol_temp + 1
          end do
          ncol_temp = ncol_temp + 1
          if (ncol.eq.0) then
             ncol = ncol_temp
          else
             if (ncol.ne.ncol_temp) then
                write(*,'(A,I4,A,I4,A,I4,A,A)') "Variable number of columns ",ncol_temp," instead of ",&
                     ncol," at line ",l," in table file: ",trim(filename)
                stop
             endif
          endif
       endif
    end do

    if ((nlin.gt.0).and.(ncol.gt.0)) then

       !! allocate table
       allocate(table(ncol,nlin))

       !! get table entries
       rewind(unit)  
       do h=1,headerlength
          read(unit,'(A)',ADVANCE="no", SIZE=nchars, IOSTAT=status) strtemp ! read header
       end do

       do l = 1, nlin 

          read(unit,'(A)',ADVANCE="no", SIZE=nchars, IOSTAT=status) strtemp ! read data column

          ! exit code -2 is ok since it marks an end-of-record during non-advancing read
          ! exit code -1 is not ok since it marks an end-of-file
          ! exit code 0 is always fine
          if (status.eq.-1) then
             write(*,'(A,A,A,I6)') 'Error in Output Namelist. End of file: ',trim(filename), &
                  ' encountered at line: ',2+l
             stop
          endif

          ! remove tab characters from string
          do c=1,nchars
             if (strtemp(c:c).eq.tab) strtemp(c:c) = ' '
          end do

          ! tokenize read string by comma-separated columns
          c=1
          pos1=1
          pos2=index(strtemp(pos1:nchars),separator)
          if (pos2.gt.0) then 
             pos2=pos2+pos1-2
          else 
             pos2=nchars
          end if

          ! loop through columns of data file
          do while (pos1.gt.0)
             if (pos2.ge.pos1) then
                if (c.le.ncol) table(c,l)=trim(adjustl(strtemp(pos1:pos2)))
             else
                if (c.le.ncol) table(c,l)=""
                pos2=pos1
             endif

             pos1=index(strtemp(pos2:nchars),separator)
             if (pos1.gt.0) then 
                pos1=pos1+pos2 
                pos2=index(strtemp(pos1:nchars),separator)
                if (pos2.gt.0) then
                   pos2=pos2+pos1-2 
                else 
                   pos2=nchars
                end if
             endif
             c=c+1
          end do
          c=c-1

          if (c.lt.ncol) then
             write(*,'(A,A,A,I6)') 'Error in Table file. Too few columns in file: ', &
                  trim(filename), ' encountered at line: ',2+l
             stop
          endif
          if (c.gt.ncol) then
             write(*,'(A,A,A,I6)') 'Error in Table file. Too many columns in file: ', &
                  trim(filename), ' encountered at line: ',2+l
             stop
          endif

       end do

    else
       write(*,'(A,A)') "Error in Table file. No columns or lines: ",trim(filename)
       stop
    end if

    close(unit)

  end subroutine read_csv_table

end module file_module
