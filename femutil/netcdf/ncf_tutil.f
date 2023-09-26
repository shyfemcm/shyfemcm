
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2019  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! convert nc files to fem files: time utilities
!
! contents :
!
! subroutine ncf_init_time(ncid,varid,nrec)	  initializes time utilities
!
! subroutine ncf_get_time(ncid,irec,atime)        returns atime for record irec
!
! subroutine ncf_print_all_time_records(ncid)     print all time records
! subroutine ncf_print_minmax_time_records(ncid)  print min/max of time records
!
! usage :
!
! call ncf_init_time(ncid,varid,nrec)		!varid is id of time records
! do irec=1,nrec
!   call ncf_get_time(ncid,irec,atime)		!atime is absolute time
!   call dts_format_abs_time(atime,dline)	!dline is string representation
! end do
!
! revision log :
!
! 16.05.2017	ggu	changed VERS_7_5_27
! 25.05.2017	ggu	changed VERS_7_5_28
! 14.11.2017	ggu	changed VERS_7_5_36
! 17.11.2017	ggu	changed VERS_7_5_37
! 05.12.2017	ggu	changed VERS_7_5_39
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.07.2018	ggu	revision control introduced
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 23.06.2022	ggu	allow for custom time correction (correct_nc_time*)
! 17.02.2023	ggu	restructured with ncf routines
!
!--------------------------------------------------------------------------

!=================================================================
	module nc_time
!=================================================================

	implicit none

	logical, save :: bcorrect = .false.	!correct time and reference
	logical, save :: bverbose = .false.	!write some debug messages
	logical, save :: bdebug = .false.	!write extra debug messages

	integer, save :: time_type = 0
	integer, save :: time_varid = 0
	integer, save :: time_nit = 0
	integer, save :: date0 = 0
	integer, save :: time0 = 0
	integer, save :: datetime0(2) = 0

	double precision, save :: atime0
	double precision, save :: time_fact

!=================================================================
	end module nc_time
!=================================================================

!*****************************************************************

        subroutine ncf_init_time(ncid,varid,nit)

! initializes time utilities - varid of time records is needed

	use nc_time
	use ncf

	implicit none

	integer ncid
	integer varid			! varid of time records
	integer nit			! total number of time records (return)

	integer var_id,att_id
	integer ifact
	integer ndims
	character*80 atext,tstring,attname
        type(var_item) :: vitem
        type(att_item) :: aitem
	type(nc_item) :: nitem

	nitem = ncf_get_nitem(ncid)
	time_varid = varid		! save for later
	var_id = varid

	attname = 'units'
	call ncf_var_inf(ncid,var_id,vitem)
	call ncf_att_id(ncid,var_id,attname,att_id)
	if( att_id == 0 ) then
	  write(6,*) 'no units attribute for time coordinate'
	  write(6,*) 'cannot parse time'
	  stop 'error stop: no units'
	end if
	call ncf_att_inf(ncid,var_id,att_id,aitem)
	atext = aitem%string

	ndims = vitem%ndims
	nit = vitem%tlen
	time_nit = nit
	if( ndims /= 1 ) stop 'error stop: more than one time dimension'

	if( bdebug ) then
	  write(6,*) 'setup_nc_time: debug'
	  write(6,*) nit
	  write(6,*) var_id
	  write(6,*) att_id
	  write(6,*) trim(atext)
	end if

	call parse_time_units(bverbose,atext &
     &				,time_type,datetime0,time_fact)

	date0 = datetime0(1)
	time0 = datetime0(2)
	call dtsini(date0,time0)
	call dts_to_abs_time(date0,time0,atime0)
	call dts_format_abs_time(atime0,tstring)

	ifact = nint(time_fact)
	if( bverbose ) then
	  write(6,*) 'setup_nc_time: ',ifact &
     &					,'  ',trim(atext) &
     &					,'  ',trim(tstring)
	end if

	call correct_nc_time_reference(ncid,bverbose)

	end

!*****************************************************************

        subroutine ncf_get_time(ncid,irec,atime)

! returns absolute time for record irec

	use nc_time
	use ncf

	implicit none

	integer ncid
	integer irec			!time record to be requested
        double precision atime		!parsed absolute time (return)

        double precision t

	call ncf_get_record(ncid,time_varid,irec,t)

	if( time_type .eq. 1 ) then
          call handle_general_time(t,atime)
	else if( time_type .eq. 2 ) then
          call handle_warf_time(t,atime)
	else
	  write(6,*) 'time_type: ',time_type
	  stop 'error stop ncf_get_time: cannot handle'
	end if

	call correct_nc_time(ncid,t,atime)

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************
! internal routines
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine parse_time_units(bverb,atext,itype,datetime0,fact)

! parses units in time attribute

	use iso8601

	implicit none

	logical bverb
	character*(*) atext
	integer itype			!type of time specification
	integer datetime0(2)		!reference time
	double precision fact		!factor to be used for time conversion

	integer ierr,off
	character*80 string

	itype = 1
	off = 1
	fact = 1.
	datetime0 = 0

	if( atext(1:10) .eq. 'days since' ) then
	  off = 12
	  fact = 86400.
	else if( atext(1:10) .eq. 'Days since' ) then
	  off = 12
	  fact = 86400.
	else if( atext(1:16) .eq. 'day as %Y%m%d.%f' ) then
	  itype = 2
	  fact = 86400.
	else if( atext(1:13) .eq. 'seconds since' ) then
	  off = 15
	  fact = 1.
	else if( atext(1:13) .eq. 'minutes since' ) then
	  off = 15
	  fact = 60.
	else if( atext(1:11) .eq. 'hours since' ) then
	  off = 13
	  fact = 3600.
	else if( atext(1:11) .eq. 'Hours since' ) then
	  off = 13
	  fact = 3600.
	else if( atext .eq. ' ' ) then	!no time coordinate
	  itype = 0
	  fact = 1.
	  write(6,*) 'no units of time coordinate found'
	  stop 'error stop parse_time_units: no unit'
	else
	  write(6,*) 'atext: ',trim(atext)
	  stop 'error stop parse_time_units: cannot parse'
	end if

	if( itype == 1 ) then
	  string = atext(off:)
	  call clean_time(string)
	  call string2date(string,datetime0,ierr)
	  if( ierr /= 0 ) then
	    write(6,*) 'error parsing time reference'
	    write(6,*) '  original string: ',trim(atext)
	    write(6,*) '   cleaned string: ',trim(string)
	    stop 'error stop parse_time_units: parsing time'
	  end if
	end if

	if( bverb ) then
	  call date2string(datetime0,string)
	  write(6,*) 'parsing date0: ',trim(atext)
	  write(6,*) 'parsed date:   ',trim(string)
	end if

	end

!*****************************************************************

	subroutine clean_time(string)

! cleans time string to make it standard

	implicit none

	character*(*) string

	integer len

	len = len_trim(string)

	if( string(len-1:len) == '.0' ) then
	  string(len-1:) = ' ' 
	  len = len - 2
	end if
	if( string(len-1:len) == ':0' ) then
	  string(len-1:) = ':00' 
	  len = len + 1
	end if

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine handle_warf_time(time,atime)

! special handling of WRF time

	use nc_time

        implicit none

        double precision time
        double precision atime

	double precision secs
        integer date

        date = time				!converts to full days
        secs = (time-date) * time_fact

	call dts_to_abs_time(date,0,atime)
	atime = atime + secs

	end

!*****************************************************************

        subroutine handle_general_time(time,atime)

! general time handling routine

	use nc_time

        implicit none

        double precision time
        double precision atime

	double precision secs

	secs = time * time_fact
	atime = atime0 + secs

        end

!*****************************************************************
!*****************************************************************
!*****************************************************************
! printing of time records
!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine ncf_print_all_time_records(ncid)

! print all time records

	use nc_time

        implicit none

        integer ncid

        integer nit,n
        double precision atime
        character*20 line

	nit = time_nit
        write(6,*) 'time records found: ',nit

        do n=1,nit
          call ncf_get_time(ncid,n,atime)
          call dts_format_abs_time(atime,line)
          write(6,*) n,atime,line
        end do

        end

!*****************************************************************

        subroutine ncf_print_minmax_time_records(ncid)

! print min/max of time records

	use nc_time

        implicit none

        integer ncid

        integer nit			!total number of time records
        double precision atime
        character*20 aline

	nit = time_nit

        if( nit == 0 ) then
          write(6,*) 'no time record found'
        else if( nit == 1 ) then
          call ncf_get_time(ncid,1,atime)
          call dts_format_abs_time(atime,aline)
          write(6,*) 'one time record found: ',atime,aline
        else
          write(6,*) 'time records found: ',nit
          call ncf_get_time(ncid,1,atime)
          call dts_format_abs_time(atime,aline)
          write(6,*) 'first time record:  ',atime,aline
          call ncf_get_time(ncid,nit,atime)
          call dts_format_abs_time(atime,aline)
          write(6,*) 'last time record:   ',atime,aline
        end if

        end

!*****************************************************************
!*****************************************************************
!*****************************************************************
! the next routines deal with bogus time references and time stamps
! in the climate files from SMHI
! time reference is often not 01.01.year but something else
! we correct it to be the first of the year
! time stamps are integers, and therefore for the hourly data the
! day is repeating 24 times with the same value
! for the monthly data we set the day to 15
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine correct_nc_time_reference(ncid,bverb)

	use nc_time

	implicit none

	integer ncid
	logical bverb

	integer nit
	integer year,month,day
	integer ifact
	character*80 tstring

	if( .not. bcorrect ) return

	call unpackdate(date0,year,month,day)

	if( month == 1 .and. day == 1 ) return

	if( month .ge. 10 ) then
	  year = year + 1
	end if
	day = 1
	month = 1
	time0 = 0

	call packdate(date0,year,month,day)
	datetime0(1) = date0
	datetime0(2) = time0

	date0 = datetime0(1)
	time0 = datetime0(2)
	call dtsini(date0,time0)
	call dts_to_abs_time(date0,time0,atime0)
	call dts_format_abs_time(atime0,tstring)

	ifact = time_fact
	if( bverb ) write(6,*) 'corrected time ref: ',ifact &
     &					,'  ',trim(tstring)

	end

!*****************************************************************

	subroutine correct_nc_time(ncid,t,atime)

	use nc_time

	implicit none

	integer ncid
	double precision t,atime

	integer nit
	integer date,time
	integer year,month,day
	integer hour,min,sec
	integer ifact
	integer, save :: last_month = 0
	integer, save :: last_t = 0
	integer, save :: last_day = 0
	integer, save :: last_hour = 0
	character*80 tstring

	if( .not. bcorrect ) return

	if( t < last_t ) then
	  last_t = t
	  last_month = 0
	  last_day = 0
	  last_hour = 0
	end if

	call dts_from_abs_time(date,time,atime)
	call unpackdate(date,year,month,day)
	call unpacktime(time,hour,min,sec)

	nit = time_nit
        !call nc_get_time_recs(ncid,nit)

	if( nit == 12 ) then		!monthly files
	  last_month = last_month + 1
	  month = last_month
	  day = 15
	else
	  if( day == last_day ) then
	    last_hour = last_hour + 1
	  else
	    last_day = day
	    last_hour = 0
	  end if
	  hour = last_hour
	end if

	!write(6,*) '+++',year,month,day,hour,t
	call packdate(date,year,month,day)
	call packtime(time,hour,min,sec)

	call dts_to_abs_time(date,time,atime)
	!call dts_format_abs_time(atime,tstring)

	end

!*****************************************************************

