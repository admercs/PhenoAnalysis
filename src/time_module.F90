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

module time_module

  ! time variables
  integer :: year                       ! simulated year
  integer :: month                      ! simulated month
  integer :: day                        ! simulated day
  integer :: hour                       ! simulated hour
  integer :: minute                     ! simulated minute
  integer :: second                     ! simulated second
  integer :: doy                        ! simulated day-of-year
  logical :: newyear, newmonth, newday  ! start of year, month, day
  logical :: oldyear, oldmonth, oldday  ! end of year, momth, day
  logical :: isleap                     ! is a leap year (366 days) or not (365 days)
  integer :: timestep                   ! model time step
  integer :: timestep_start             ! model start time step
  integer :: timestep_end               ! model end time step
  integer :: timestep_yearstart         ! model time step at start of year
  integer :: timestep_monthstart        ! model time step at start of month
  integer :: timestep_daystart          ! model time step at start of day
  integer :: timestep_yearend           ! model time step at end of year
  integer :: timestep_monthend          ! model time step at end of month
  integer :: timestep_dayend            ! model time step at end of day

contains

  subroutine update_time(date_start, date_end, dt, initialize)

    ! calculates model time/date variables
    ! attention: timestep is relative to timestep_start (and timestep_end)

    implicit none

    ! parameters
    character(len=*), intent(in) :: date_start
    character(len=*), intent(in) :: date_end
    integer, intent(in) :: dt
    logical, intent(in) :: initialize

    ! local variables
    integer :: yyyy, mm, dd

    if (initialize) then

       call date_string_to_integer(date_start,yyyy,mm,dd,.false.)
       timestep_start = calc_timestep(yyyy,mm,dd,0,0,0,dt)

       call date_string_to_integer(date_end,yyyy,mm,dd,.true.)
       timestep_end = calc_timestep(yyyy,mm,dd,0,0,0,dt) + 86400/dt - 1

       timestep = 1

    else

       timestep = timestep + 1

    endif

    ! calculate current date and doy
    call calc_date(timestep+timestep_start-1,year,month,day,hour,minute,second,dt)
    doy = calc_doy(year,month,day)
    isleap = leap_year(year)

    ! calculate timesteps at start and end of year/month/day
    timestep_yearstart = calc_timestep(year,1,1,0,0,0,dt) - timestep_start + 1
    timestep_monthstart = calc_timestep(year,month,1,0,0,0,dt) - timestep_start + 1
    timestep_daystart = calc_timestep(year,month,day,0,0,0,dt) - timestep_start + 1
    timestep_yearend = calc_timestep(year+1,1,1,0,0,0,dt) - 1 - timestep_start + 1
    timestep_monthend = calc_timestep(year,month,1,0,0,0,dt) + &
         calc_daysofmonth(year,month)*86400/dt - 1 - timestep_start + 1
    timestep_dayend = calc_timestep(year,month,day,0,0,0,dt) + 86400/dt - 1 - timestep_start + 1

    ! determine whether we have a start or end of year/month/day
    if (timestep.eq.timestep_yearstart) then
       newyear = .true.
    else
       newyear = .false.
    endif
    if (timestep.eq.timestep_monthstart) then
       newmonth = .true.
    else
       newmonth = .false.
    endif
    if (timestep.eq.timestep_daystart) then
       newday = .true.
    else
       newday = .false.
    endif
    if (timestep.eq.timestep_yearend) then
       oldyear = .true.
    else
       oldyear = .false.
    endif
    if (timestep.eq.timestep_monthend) then
       oldmonth = .true.
    else
       oldmonth = .false.
    endif
    if (timestep.eq.timestep_dayend) then
       oldday = .true.
    else
       oldday = .false.
    endif

!    print*,timestep,timestep_start,timestep_end
!    print*,timestep_yearstart,timestep_monthstart,timestep_daystart
!    print*,newyear,newmonth,newday
!    print*,timestep_yearend,timestep_monthend,timestep_dayend
!    print*,oldyear,oldmonth,oldday
!    print*,year,month,day,hour,minute,second
!    print*,doy
!    print*,dt

  end subroutine update_time

  integer function calc_timestep(year,month,day,hour,minute,second,dt,calendar)
    ! calculates an integer time step from the calendar date by use of the 
    ! time step length dt

    implicit none

    ! arguments
    integer, intent(in) :: year                         ! year (starting at -4713 ,-4714 for gregorian) 
    integer, intent(in) :: month                        ! month 1-12
    integer, intent(in) :: day                          ! day 1-31
    integer, intent(in) :: hour                         ! hour 0-23 (utc)
    integer, intent(in) :: minute                       ! minute 0-59 (utc)
    integer, intent(in) :: second                       ! second 0-59 (utc)
    integer, intent(in) :: dt                           ! time step length (seconds)
    character(len=*), intent(in), optional :: calendar  ! type of calendar: proleptic gregorian, gregorian or julian

    ! local variables
    real(kind=8) :: jd

    ! calculate the julian date from the calendar date
    jd = calc_julian_date(year,month,day,hour,minute,second,calendar)

    ! calculate the time step
    calc_timestep = nint(jd / real(dt,kind=8) * 86400.d0)

    return
  end function calc_timestep

  subroutine date_string_to_integer(sdate,yyyy,mm,dd,last,calendar)
    ! interprets string and returns integer year, month and day

    implicit none

    ! arguments
    character(*), intent(in) :: sdate       ! string date YYYY, YYYYMM or YYYYMMDD
    integer, intent(out) :: yyyy            ! year
    integer, intent(out) :: mm              ! month
    integer, intent(out) :: dd              ! day
    logical, intent(in)  :: last            ! if true, set the output date to last day in year or month
    character(len=*), intent(in), optional :: calendar  ! type of calendar: proleptic gregorian, gregorian or julian

    ! local variables
    integer :: date

    read (sdate, *) date
    select case (len_trim(sdate))
    case(4)
       yyyy = date
       if (last) then
          mm = 12
          dd = calc_daysofmonth(yyyy,mm,calendar)
       else
          mm = 1
          dd = 1
       endif
    case(6)
       yyyy = date / 100
       mm = date - yyyy*100
       if (last) then
          dd = calc_daysofmonth(yyyy,mm,calendar)
       else
          dd = 1
       endif
    case(8)
       yyyy = date / 10000
       mm = (date - yyyy*10000) / 100
       dd = date - yyyy*10000 - mm*100
    case DEFAULT
       write(*,'(A,A)') "Wrong date format: ",trim(sdate)
    end select
    
  end subroutine date_string_to_integer

  subroutine calc_date(timestep,year,month,day,hour,minute,second,dt,calendar)
    ! calculates the calendar date from the integer time step with time step length dt

    implicit none

    ! arguments
    integer, intent(in) :: timestep                     ! integer time step
    integer, intent(out) :: year                         ! year (starting at -4713 ,-4714 for gregorian) 
    integer, intent(out) :: month                        ! month 1-12
    integer, intent(out) :: day                          ! day 1-31
    integer, intent(out) :: hour                         ! hour 0-23 (utc)
    integer, intent(out) :: minute                       ! minute 0-59 (utc)
    integer, intent(out) :: second                       ! second 0-59 (utc)
    integer, intent(in) :: dt                           ! time step length (seconds)
    character(len=*), intent(in), optional :: calendar  ! type of calendar: proleptic gregorian, gregorian or julian

    ! local variables
    real(kind=8) :: jd

    ! calculate the julian date
    jd = real(timestep,kind=8) / 86400.d0 * real(dt,kind=8)

    ! calculate the calendar date from the julian date
    call calc_calendar_date(jd,year,month,day,hour,minute,second,calendar)

    return
  end subroutine calc_date

  integer function calc_doy(year, month, day, calendar)
    ! returns the day of year (1-366) from year, month and day

    implicit none

    ! arguments
    integer, intent(in) :: year                         ! year  
    integer, intent(in) :: month                        ! month 1-12
    integer, intent(in) :: day                          ! day 1-31
    character(len=*), intent(in), optional :: calendar  ! type of calendar: proleptic gregorian, gregorian or julian

    ! local variables
    integer, parameter :: monthsum_noleap(12) = (/0,31,59,90,120,151,181,212,243,273,304,334/)
    integer, parameter :: monthsum_leap(12) = (/0,31,60,91,121,152,182,213,244,274,305,335/)

    if (leap_year(year,calendar)) then
       calc_doy = monthsum_leap(month) + day
    else
       calc_doy = monthsum_noleap(month) + day
    endif

    return

  end function calc_doy

  integer function calc_daysofmonth(year, month, calendar)
    ! returns the total number of days of current month 

    implicit none

    ! arguments
    integer, intent(in) :: year                         ! year  
    integer, intent(in) :: month                        ! month 1-12
    character(len=*), intent(in), optional :: calendar  ! type of calendar: proleptic gregorian, gregorian or julian

    ! local variables
    integer, parameter :: monthdays_noleap(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    integer, parameter :: monthdays_leap(12) = (/31,29,31,30,31,30,31,31,30,31,30,31/)

    if (leap_year(year,calendar)) then
       calc_daysofmonth = monthdays_leap(month)
    else
       calc_daysofmonth = monthdays_noleap(month)
    endif

    return

  end function calc_daysofmonth

  integer function calc_daysofyear(year, calendar)
    ! returns the total number of days of current month 

    implicit none

    ! arguments
    integer, intent(in) :: year                         ! year  
    character(len=*), intent(in), optional :: calendar  ! type of calendar: proleptic gregorian, gregorian or julian

    ! local variables
    integer, parameter :: yeardays_noleap = 365
    integer, parameter :: yeardays_leap = 366

    if (leap_year(year,calendar)) then
       calc_daysofyear = yeardays_leap
    else
       calc_daysofyear = yeardays_noleap
    endif

    return

  end function calc_daysofyear

  subroutine calc_month_day(doy, year, month, day, calendar)
    ! returns the year, month and day from day of year (1-366)
    ! only works for gregorian calendar

    implicit none

    ! arguments
    integer, intent(in) :: doy                           ! day of year (1..366)
    integer, intent(in) :: year                          ! year  
    integer, intent(out) :: month                        ! month 1-12
    integer, intent(out) :: day                          ! day 1-31
    character(len=*), intent(in), optional :: calendar   ! type of calendar: proleptic gregorian, gregorian or julian

    ! local variables
    integer, parameter :: monthsum_noleap(12) = (/31,59,90,120,151,181,212,243,273,304,334,365/)
    integer, parameter :: monthsum_leap(12) = (/31,60,91,121,152,182,213,244,274,305,335,366/)

    month = 1
    if (leap_year(year,calendar)) then
       do while ((month < 12).and.(monthsum_leap(month).lt.doy))
          month = month + 1
       end do
       if (month.gt.1) then
          day = doy - monthsum_leap(month-1)
       else
          day = doy
       endif
    else
       do while ((month < 12).and.(monthsum_noleap(month).lt.doy))
          month = month + 1
       end do
       if (month.gt.1) then
          day = doy - monthsum_noleap(month-1)
       else
          day = doy
       endif
    endif

    return

  end subroutine calc_month_day

  logical function leap_year(year, calendar)
    ! returns whether year is a leap year (with 29 days in February)

    implicit none

    ! arguments
    integer, intent(in) :: year                         ! year (starting at -4713 ,-4714 for gregorian) 
    character(len=*), intent(in), optional :: calendar  ! type of calendar: proleptic gregorian, gregorian or julian

    ! check for type of calendar
    if (gregorian(calendar)) then
       leap_year = ((mod(year,4).eq.0).and.(mod(year,100).ne.0)).or.(mod(year,400).eq.0)
    else
       leap_year = (mod(year,4).eq.0)
    endif

    return

  end function leap_year

  logical function gregorian(calendar)
    ! checks whether the argument specifies a gregorian calendar
    ! returns TRUE for gregorian and FALSE for julian calendar

    implicit none

    ! arguments
    character(len=*), intent(in), optional :: calendar  ! type of calendar: proleptic gregorian, gregorian or julian

    gregorian = .true.
    if (present(calendar)) then
       select case (trim(calendar))
       case ("julian") 
          gregorian = .false.
       case ("gregorian")
          gregorian = .true.
       case ("proleptic_gregorian")
          gregorian = .true.
       case default
          write(*,"(A,A,A)") "Unknown calendar type: ",trim(calendar),". Assuming gregorian calendar."
       end select
    endif

    return

  end function gregorian

  real(kind=8) function calc_julian_date(year, month, day, hour, minute, second, calendar)
    ! calculates the julian day from the calendar date
    ! formulas from http://en.wikipedia.org/wiki/Julian_day
    ! and http://www.tondering.dk/claus/cal/julperiod.php

    ! JD = 0 for 1 Jan 4713 BC 12:00 in the Julian calendar 
    ! JD = 0 for 24 Nov 4714 BC 12:00 in the Gregorian calendar
    ! the proleptic Gregorian calendar should be specified for Gregorian dates prior to 1582 since
    ! before this date the Julian calendar was in place
    ! 4713 BC corresponds to year -4712, so 1 BC is year 0

    implicit none

    ! arguments
    integer, intent(in) :: year                         ! year (starting at -4713 ,-4714 for gregorian) 
    integer, intent(in) :: month                        ! month 1-12
    integer, intent(in) :: day                          ! day 1-31
    integer, intent(in) :: hour                         ! hour 0-23 (utc)
    integer, intent(in) :: minute                       ! minute 0-59 (utc)
    integer, intent(in) :: second                       ! second 0-59 (utc)
    character(len=*), intent(in), optional :: calendar  ! type of calendar: proleptic gregorian, gregorian or julian

    ! local variables
    integer :: a, y, m

    ! precompute julian day parameters
    a = (14 - month) / 12
    y = year + 4800 - a
    m = month + 12*a - 3

    ! compute julian day
    if (gregorian(calendar)) then
       calc_julian_date = real(day + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045,kind=8)
    else
       calc_julian_date = real(day + (153*m + 2)/5 + 365*y + y/4 - 32083,kind=8)
    endif

    ! add fractional day
    calc_julian_date = calc_julian_date + (real(hour,kind=8) - 12.d0)/24.d0 + &
         real(minute,kind=8)/1440.d0 + real(second,kind=8)/86400.d0

    return 

  end function calc_julian_date

  subroutine calc_calendar_date( jd, year, month, day, hour, minute, second, calendar)
    ! calculates the calendar date from the julian day
    ! formulas from http://en.wikipedia.org/wiki/Julian_day
    ! and http://www.tondering.dk/claus/cal/julperiod.php

    ! JD = 0 for 1 Jan 4713 BC 12:00 in the Julian calendar 
    ! JD = 0 for 24 Nov 4714 BC 12:00 in the Gregorian calendar
    ! the proleptic Gregorian calendar should be specified for Gregorian dates prior to 1582 since
    ! before this date the Julian calendar was in place
    ! 4713 BC corresponds to year -4712, so 1 BC is year 0

    implicit none

    ! arguments
    real(kind=8), intent(in) ::  jd                      ! julian day number
    integer, intent(out) :: year                         ! year (starting at -4713 ,-4714 for gregorian)
    integer, intent(out) :: month                        ! month 1-12
    integer, intent(out) :: day                          ! day 1-31
    integer, intent(out) :: hour                         ! hour 0-23 (utc)
    integer, intent(out) :: minute                       ! minute 0-59 (utc)
    integer, intent(out) :: second                       ! second 0-59 (utc)
    character(len=*), intent(in), optional :: calendar   ! type of calendar: proleptic gregorian, gregorian or julian

    ! local variables
    integer :: a, b, c, d, e, m
    integer :: daysec    ! number of elapsed seconds within the current day

    if (gregorian(calendar)) then
       a = int(jd + 0.5d0) + 32044
       b = (4*a + 3) / 146097
       c = a - (146097*b) / 4
    else
       a = 0
       b = 0
       c = int(jd + 0.5d0) + 32082
    endif

    d = (4*c + 3) / 1461
    e = c - (1461*d) / 4
    m = (5*e + 2) / 153

    day = e - (153*m + 2) / 5 + 1
    month = m + 3 - 12 * (m/10)
    year = 100*b + d - 4800 + m/10

    ! we need to extract total integer seconds elapsed since 0 UTC first in order to avoid rounding errors
    daysec = nint(3600.d0*24.d0*(jd + 0.5d0 - int(jd + 0.5d0)))
    hour = daysec / 3600
    minute = (daysec - hour*3600) / 60
    second = daysec - hour*3600 - minute*60

    return 

  end subroutine calc_calendar_date

end module time_module
