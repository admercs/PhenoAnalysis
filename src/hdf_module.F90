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

module hdf_module

  ! 2009/02/20 Reto Stockli (Blue Marble Research)

  ! This module is the interface to the F77-style HDF4 library
  ! It can be used to access the MODIS HDF4 EOS files from within Fortran 90

  ! Dependencies:
  ! a compiled HDF4 library with the Fortran 77 interface enabled and NetCDF disabled
  ! IMPORTANT: NetCDF needs to be disabled in the HDF4 library, since else we get multiple symbols errors
  ! during linking here when the separate NetCDF 3.x or 4.x library is linked.

contains

  subroutine open_hdf_file(filename, sd_id)
    implicit none

    ! inout
    character(len=*), intent(in) :: filename
    integer, intent(out) :: sd_id

    ! HDF4
    integer :: sfstart
    integer, parameter :: DFACC_READ = 1

    sd_id = sfstart(filename,DFACC_READ)

  end subroutine open_hdf_file


  subroutine close_hdf_file(sd_id)
    implicit none

    ! inout
    integer, intent(in) :: sd_id

    ! HDF4
    integer :: sfend

    ! local
    integer :: status

    status = sfend(sd_id)

  end subroutine close_hdf_file

  subroutine read_hdf_attribute(sds_id, name, val)
    ! Only reads single attribute now (no attributes with 2 or more values)
    ! returns 32 bit floating point
    implicit none

    ! inout
    integer, intent(in) :: sds_id
    character(len=*), intent(in) :: name
    real(kind=4), intent(out) :: val

    ! HDF4
    integer :: sffattr, sfrattr, sfgainfo

    ! local
    character(len=64) :: attr_name
    integer :: ret
    integer :: cnt
    integer :: attr_index
    integer :: num_type
    integer(kind=2) :: att_int16
    character(len=1) :: att_char
    integer(kind=4) :: att_int32
    real(kind=4) :: att_real
    real(kind=4) :: att_double

    attr_index = sffattr (sds_id,name)
    if (attr_index.ge.0) then
       ret = sfgainfo (sds_id, attr_index, attr_name, num_type, cnt)

       if (cnt.eq.1) then
          !        print*,attr_name,num_type,cnt,attr_index
          if ((num_type.eq.3).or.(num_type.eq.4.).or.(num_type.eq.20).or.(num_type.eq.21)) then
             ret = sfrattr (sds_id, attr_index, att_char)
             val = real(iachar(att_char),kind=4)
          else if ((num_type.eq.22).or.(num_type.eq.23)) then
             ret = sfrattr (sds_id, attr_index, att_int16)
             val = real(att_int16,kind=4)
          else if ((num_type.eq.24).or.(num_type.eq.25)) then
             ret = sfrattr (sds_id, attr_index, att_int32)
             val = real(att_int32,kind=4)
          else if (num_type.eq.5) then
             ret = sfrattr (sds_id, attr_index, att_real)
             val = att_real
          else if (num_type.eq.6) then
             ret = sfrattr (sds_id, attr_index, att_double)
             val = real(att_double,kind=4)
         else
             val = 0.0
          endif

       else
          val = 0.0
       endif

    else
       val = 0.0
    endif

  end subroutine read_hdf_attribute

  subroutine read_hdf_file(sd_id,varname,x,y,nx,ny,val_char,val_int16,val_int32,val_real,val_double)

    implicit none

    ! inout
    character(len=*), intent(in) :: varname
    integer, intent(in) :: sd_id
    integer, intent(in) :: x
    integer, intent(in) :: y
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    character(len=1), intent(inout), optional :: val_char(:,:)
    integer(kind=2), intent(inout), optional :: val_int16(:,:)
    integer(kind=4), intent(inout), optional :: val_int32(:,:)
    real(kind=4), intent(inout), optional :: val_real(:,:)
    real(kind=4), intent(inout), optional :: val_double(:,:)

    ! local
    integer :: sds_id, sds_index, ret
!    real(kind=4) :: fill_value, scale_factor, add_offset

    integer :: rank
    integer :: dim_sizes(2)
    integer :: num_type
    integer :: num_attrs
    character(len=64) :: sds_name

    integer :: start(2), stride(2), edge(2)

    ! HDF4
    integer :: sfn2index, sfselect, sfendacc
    integer :: sfginfo, sfrdata, sfrcdata

    ! select SDS
    sds_index = sfn2index(sd_id,varname)
    sds_id = sfselect (sd_id, sds_index)

    if (sds_id.lt.0) write(*,'(A,A,A)') 'ERROR: SDS NAME ',varname,' NOT FOUND IN FILE'

    ! read attributes (for some reason scale factor is always 0.0 ?!
!    call read_hdf_attribute(sds_id,'_FillValue', fill_value)  
!    write(*,'(A,F15.8)') 'Fill Value: ',fill_value
!    call read_hdf_attribute(sds_id,'scale_factor', scale_factor)  
!    write(*,'(A,F15.8)') 'Scale:      ',scale_factor
!    call read_hdf_attribute(sds_id,'add_offset', add_offset)  
!    write(*,'(A,F15.8)') 'Offset:     ',add_offset

    ! read data
    start(1)=x   ! 0-based
    start(2)=y   ! 0-based
    stride(1)=1
    stride(2)=1
    edge(1)=nx
    edge(2)=ny

    ret = sfginfo(sds_id,sds_name,rank,dim_sizes,num_type,num_attrs)
    !write(*,'(I8,A,5I6) sds_id,trim(sds_name),rank,dim_sizes,num_type,num_attrs

!    print*,start
!    print*,edge

!    print*,size(val_char,1)
!    print*,size(val_char,2)

    ! Attention: we deliver integer values as output (may not be ok for non-integer MODIS data)
    if ((num_type.eq.3).or.(num_type.eq.4.).or.(num_type.eq.20).or.(num_type.eq.21)) then
       if (present(val_char)) ret = sfrcdata (sds_id, start, stride , edge, val_char)
    else if ((num_type.eq.22).or.(num_type.eq.23)) then
       if (present(val_int16)) ret = sfrdata (sds_id, start, stride, edge, val_int16)
    else if ((num_type.eq.24).or.(num_type.eq.25)) then
       if (present(val_int32)) ret = sfrdata (sds_id, start, stride, edge, val_int32)
    else if (num_type.eq.5) then
       if (present(val_real)) ret = sfrdata (sds_id, start, stride, edge, val_real)
    else if (num_type.eq.6) then
       if (present(val_double)) ret = sfrdata (sds_id, start, stride, edge, val_double)
    endif

    ! stop reading SDS
    ret = sfendacc(sds_id)

  end subroutine read_hdf_file

end module hdf_module
