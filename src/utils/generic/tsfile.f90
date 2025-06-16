
!--------------------------------------------------------------------------
!
!    Copyright (C) 2014-2015,2017-2019  Georg Umgiesser
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

! routines for reading time series
!
! revision log :
!
! 30.05.2014	ggu	changed VERS_6_1_76
! 18.06.2014	ggu	changed VERS_6_1_77
! 07.07.2014	ggu	changed VERS_6_1_79
! 18.07.2014	ggu	changed VERS_7_0_1
! 20.10.2014	ggu	integrating datetime into time series
! 30.10.2014	ggu	changed VERS_7_0_4
! 12.12.2014	ggu	changed VERS_7_0_9
! 10.02.2015	ggu	length of line set to 2048
! 26.02.2015	ggu	changed VERS_7_1_5
! 15.05.2017	ggu	new version to also read time string
! 04.11.2017	ggu	changed VERS_7_5_34
! 14.11.2017	ggu	changed VERS_7_5_36
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.04.2018	ggu	changed VERS_7_5_43
! 15.10.2018	ggu	set nvar=0 in ts_get_file_info() to avoid segfault
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.12.2019	ggu	return always time in dtime (converting from string)
! 08.12.2020	ggu	better documentation
! 17.04.2024	ggu	new routine ts_get_time_info()
!
! notes :
!
! keyword example setting date:		"#date: 20071001 0"
! time column as string:		"2007-10-01::00:00:00"
!
! possible keywords:
!
! #date: 19970101
! #info: any other info
! #vars: var1 var2 ...
! #any comment goes here
!
!*************************************************************

	function check_ts_file(file)

! checks if file is time series file

	implicit none

	logical check_ts_file
	character*(*) file	!file name

	integer nvar

	call ts_get_file_info(file,nvar)

	check_ts_file = ( nvar > 0 )

	end

!*************************************************************

	subroutine ts_get_file_info(file,nvar)

! get info on number of variables in time series file (TS)
!
! nvar is returned as the number of available data columns (except time)
! nvar <= 0 for error or no TS file

	implicit none

	character*(*) file	!file name
	integer nvar		!variables (columns) in file (except time)

	integer iunit,i,ierr
	integer datetime(2)
	double precision dtime
	real, allocatable :: f(:)
	character*80 varline

	nvar = 0
	call ts_open_file(file,nvar,datetime,varline,iunit)
	if( nvar <= 0 ) return
	if( iunit <= 0 ) return

	allocate(f(nvar))

	do i=1,3
	  call ts_read_next_record(iunit,nvar,dtime,f,datetime,ierr)
	  if( ierr /= 0 ) exit
	end do

	if( ierr > 0 ) nvar = 0	!return error only on read error, not EOF

	close(iunit)

	end

!*************************************************************

	subroutine ts_get_time_info(file,itime)

! get info on time column: relative or absolute

	use iso8601

	implicit none

	character*(*) file	!file name
	integer itime		!0: relative  1:absolute  -1:error

	integer iunit,ierr,ios
	integer is,iend,nvar
	integer datetime(2)
	double precision dtime
	character*2048 line
	character*80 stime

	integer ichafs,istod,istot

	itime = -1

	nvar = 0
	call ts_open_file(file,nvar,datetime,line,iunit)
	if( iunit <= 0 ) return
	if( nvar <= 0 ) goto 99

	do
	  read(iunit,'(a)',iostat=ios) line
	  if( ios /= 0 ) goto 99
	  if( line(1:1) /= '#' ) exit
	end do

	is = ichafs(line)
	iend = len_trim(line)
	if( iend > 2000 ) goto 99

	ios = istod(line,dtime,is)		!read time column
	if( ios == -1 ) then			!time colum may be string
	  ios = istot(line,stime,is)		!read time column as string
	  if( ios /= 1 ) goto 99
	  dtime = 0.
	  call string2date(stime,datetime,ierr)
	  if( ierr /= 0 ) goto 99
	  itime = 1
	else
	  itime = 0
	end if

   99	continue
	close(iunit)

	end

!*************************************************************

	subroutine ts_get_extra_time(file,dtime,datetime)

! gets extra time information at end of line (old format)

	use iso8601

	implicit none

	character*(*) file	!file name
	double precision dtime	!relative time of first data record
	integer datetime(2)	!extra date info on first data record

	integer nvar		!variables (columns) in file (except time)
	integer iunit
	integer ios,i,is,ierr
	character*2048 line,dummy
	character*80 varline

	integer istot,istod

	dtime = 0.
	datetime = 0

	call ts_open_file(file,nvar,datetime,varline,iunit)
	if( nvar <= 0 .or. iunit <= 0 ) return

	line = ' '
	do while( line == ' ' )
	  read(iunit,'(a)',iostat=ios) line
	  if( ios /= 0 ) return
	  if( line == ' ' ) line = ' '
	  if( line(1:1) == '#' ) line = ' '
	end do

	is = 1
	ios = istod(line,dtime,is)		!read time column
	if( ios /= 1 ) goto 1

	do i=1,nvar
	  ios = istot(line,dummy,is)
	  if( ios /= 1 ) goto 1
	end do

	ios = istot(line,dummy,is)
	if( ios == 1 ) then
	  call string2date(dummy(1:20),datetime,ierr)
	end if

    1	continue
	close(iunit)

	end

!*************************************************************

	subroutine ts_open_file(file,nvar,datetime,varline,iunit)

! opens file (on error nvar=-1 and iunit=0)

! file is only input value
! all other arguments are return values

	implicit none

	character*(*) file	!file name
	integer nvar		!variables (columns) in file (except time)
	integer datetime(2)	!info on absolute time
	character*(*) varline	!description of variables (if available)
	integer iunit		!unit of opened file

	integer ierr
	real f(nvar)
	double precision dtime

	integer ifileo

!------------------------------------------------------
! open file
!------------------------------------------------------

	nvar = 0
	datetime = 0

	iunit = ifileo(0,file,'formatted','old')

	if( iunit <= 0 ) return

	call ts_get_var_names(iunit,varline)

	call ts_read_next_record(iunit,nvar,dtime,f,datetime,ierr)

	if( ierr == 0 ) then
	  backspace(iunit)
	else
	  close(iunit)
	  nvar = -1
	  iunit = 0
	end if

!------------------------------------------------------
! end of routine
!------------------------------------------------------

	end

!*************************************************************

	subroutine ts_peek_next_record(iunit,nvar,dtime,f,datetime,ierr)

! peeks into one record of time series file
!
! see ts_read_next_record for description of arguments

	implicit none

	integer iunit
	integer nvar
	double precision dtime
	real f(nvar)
	integer datetime(2)
	integer ierr

	call ts_read_next_record(iunit,nvar,dtime,f,datetime,ierr)

	backspace(iunit)

	end

!*************************************************************

	subroutine ts_read_next_record(iunit,nvar,dtime,f,datetime,ierr)

! reads one record of time series file

! iunit must be file unit to read from
! if nvar > 0 on entry, number of variables read must be this, otherwise error
! nvar is number of variables read on return
! all other arguments are return values

	use iso8601

	implicit none

	integer iunit		!file unit (input)
	integer nvar		!variables (columns) in file (except time)
	double precision dtime	!relative time
	real f(nvar)		!data for variables
	integer datetime(2)	!info on absolute time
	integer ierr		!return code (==0 ok, <0 EOF, >0 error)

	integer i,iend,ios,is,nvar0
	double precision t
	character*2048 line
	character*80 info,stime

	logical ts_has_keyword
	integer ichafs,istod,istot,iscanf

!------------------------------------------------------
! check if file is open
!------------------------------------------------------

	ierr = 1
	if( iunit <= 0 ) return

!------------------------------------------------------
! read record
!------------------------------------------------------

	nvar0 = nvar

	!------------------------------------------------------
	! read keywords and skip empty lines
	!------------------------------------------------------

	line = ' '
	do while( line == ' ' )
	  read(iunit,'(a)',iostat=ierr) line
	  if( ierr /= 0 ) return

	  if( line == ' ' ) cycle
	  if( ts_has_keyword(line) ) then
	    call ts_parse_keyword(iunit,line,datetime,info)
	    line = ' '
	  end if
	end do

	!------------------------------------------------------
	! see if line is too long
	!------------------------------------------------------

	ierr = 2
	is = ichafs(line)
	iend = len_trim(line)
	if( iend > 2000 ) return

	!------------------------------------------------------
	! read time column
	!------------------------------------------------------

	ierr = 3
	ios = istod(line,dtime,is)		!read time column
	if( ios == -1 ) then			!time colum may be string
	  ios = istot(line,stime,is)		!read time column as string
	  if( ios /= 1 ) return
	  dtime = 0.
	  call string2date(stime,datetime,ierr)
	  call dts_to_abs_time(datetime(1),datetime(2),dtime) !absolute time
	  if( ierr /= 0 ) return
	end if

	!------------------------------------------------------
	! read rest of variables
	!------------------------------------------------------

	ierr = 4
	nvar = iscanf(line(is:),f,nvar0)	!get values, maybe only counting

	if( nvar < 0 ) nvar = -nvar-1		!read error in number -nvar
	if( nvar <= 0 ) return			!no data found
	if( nvar0 > 0 .and. nvar /= nvar0 ) return	!varying number of data

	!write(6,*) 'debug TS: ',ierr,nvar0,nvar,datetime,iunit

!------------------------------------------------------
! set error code
!------------------------------------------------------

	ierr = 0

!------------------------------------------------------
! end of routine
!------------------------------------------------------

	end

!*************************************************************
!*************************************************************
!*************************************************************

	subroutine ts_get_var_names(iunit,varline)

! returns names of variables if found

	implicit none

	integer iunit
	character*(*) varline

	integer ios,ioff
	character*80 key
	character*2048 line

	varline = ' '

	do
	  read(iunit,'(a2048)',iostat=ios) line
	  if( ios /= 0 ) exit
	  call ts_get_keyword(line,key,ioff)
	  if( ioff == 1 ) exit
	  if( key == 'vars' ) then
	    varline = line(ioff:)
	    exit
	  end if
	end do

	rewind(iunit)

	end

!*************************************************************
!*************************************************************
!*************************************************************

	function ts_has_keyword(line)

! checks if keyword on line

	implicit none

	logical ts_has_keyword
	character*(*) line

	integer i
	integer ichafs

	ts_has_keyword = .false.

	i = ichafs(line)
	if( i <= 0 ) return
	if( line(i:i) == '#' ) ts_has_keyword = .true.

	end

!*************************************************************

	subroutine ts_get_keyword(line,key,ioff)

! gets keyword from line - rest of line after ioff
!
! keyword looks like:   "#key:"
! example:		"#date: 20071001 0"

	implicit none

	character*(*) line
	character*(*) key
	integer ioff

	integer i,j

	integer ichafs

	key = ' '
	ioff = 1

	i = ichafs(line)
	if( i <= 0 ) return
	if( line(i:i) /= '#' ) return	!no keyword

	ioff = 2			!signal that we have found #

	do j=i+1,len(line)
	  if( line(j:j) == ' ' ) return	!no white space allowed
	  if( line(j:j) == ':' ) then	!end of keyword found
	    if( j-i-1 <= 0 ) return	!no keyword present
	    key = line(i+1:j-1)
	    ioff = j+1
	    return
	  end if
	end do

	return
	end

!*************************************************************

	subroutine ts_parse_keyword(iunit,line,datetime,info)

! parses keyword

	implicit none

	integer iunit
	character*(*) line
	integer datetime(2)
	character*(*) info

	character*10 key
	character*75 file
	integer ioff

	call ts_get_keyword(line,key,ioff)

	if( key == 'date' ) then		!date
	  call ts_parse_datetime(line(ioff:),datetime)
	else if( key == 'info' ) then		!info
	  call ts_parse_info(line(ioff:),info)
	else if( key == 'vars' ) then		!vars
	  call ts_parse_info(line(ioff:),info)
	else if( key == ' ' ) then		!nothing
	else
	  call filna(iunit,file)
	  write(6,*) 'not recognized keyword: ',key
	  write(6,*) 'file open at unit: ',iunit
	  write(6,*) 'file name: ',file
	  stop 'error stop ts_parse_keyword: not recognized keyword'
	end if

	end

!*************************************************************

	subroutine ts_parse_datetime(line,datetime)

! parses keyword date

	implicit none

	character*(*) line
	integer datetime(2)

	integer j
	double precision d(3)

	integer iscand

	datetime = 0

	j = iscand(line,d,2)

	if( j >= 1 ) datetime(1) = nint(d(1))
	if( j >= 2 ) datetime(2) = nint(d(2))

	if( datetime(1) > 0 .and. datetime(1) < 10000 ) then
	  datetime(1) = 10000*datetime(1) + 101
	end if

	end

!*************************************************************

	subroutine ts_parse_info(line,info)

! parses keyword info

	implicit none

	character*(*) line
	character*(*) info

	integer i
	integer ichafs

	info = ' '

	i = ichafs(line)
	if( i <= 0 ) return

	info = line(i:)

	end

!*************************************************************
!*************************************************************
!*************************************************************

	subroutine ts_test

! testing TS files

	implicit none

	integer nvar
	integer iunit
	integer datetime(2)
	integer ierr
	integer nrec,i
	double precision it
	real f(10)
	character*60 file
	character*80 varline
	character*20 line

	line = ' '
	nrec = 0
	file = 'ts_test.txt'

	call ts_get_file_info(file,nvar)
	write(6,*) 'file info: nvar = ',nvar

	call ts_open_file(file,nvar,datetime,varline,iunit)
	write(6,*) 'file open: nvar = ',nvar
	write(6,*) 'file open: iunit = ',iunit
	write(6,*) 'file open: datetime = ',datetime
	if( datetime(1) > 0 ) then
	  !call dtsini(datetime(1),datetime(2))
	end if

	do
	  call ts_read_next_record(iunit,nvar,it,f,datetime,ierr)
	  if( ierr .ne. 0 ) exit
	  if( datetime(1) > 0 ) then
	    write(6,*) 'datetime: ',datetime
	    !call dtsini(datetime(1),datetime(2))
	  end if
	  nrec = nrec + 1
	  !call dtsgf(nint(it),line)
	  !write(6,*) it,(f(i),i=1,nvar)
	  if( mod(nrec,1) .eq. 0 ) write(6,*) it,f(1),line
	end do

	write(6,*) nrec,' records read'

	end

!*************************************************************
!	program ts_test_main
!	call ts_test
!	end
!*************************************************************

