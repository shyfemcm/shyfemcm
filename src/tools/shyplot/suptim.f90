
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999,2005,2008,2010,2013-2015,2013-2015  Georg Umgiesser
!    Copyright (C) 2018-2019  Georg Umgiesser
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

!  revision log :
! 
!  12.02.1999	ggu	adapted to auto mode
!  27.05.2005	ggu	increase nrec always in oktime (even when it is the same)
!  13.11.2008	ggu	in oktime() increase irec only for new time
!  06.12.2008	ggu	in oktime() set itact to actual time
!  23.03.2010	ggu	changed v6.1.1
!  09.10.2010	ggu	in oktime() handle negative itfreq
!  15.12.2010	ggu	changed VERS_6_1_14
!  05.09.2013	ggu	new routine endtime()
!  12.09.2013	ggu	changed VERS_6_1_67
!  20.10.2014	ggu	completely restructured, old routines deleted
!  30.10.2014	ggu	changed VERS_7_0_4
!  05.11.2015	ggu	changed VERS_7_3_12
!  18.12.2018	ggu	inlude file substituted with module timlim
!  21.05.2019	ggu	changed VERS_7_5_62
! 
! ******************************************************

	module timlim

	implicit none

        integer, save :: itmin,itmax,itfreq,nrec,idto,itact
        double precision, save :: atime0,atimeact,atimemin,atimemax

	end module timlim

! ******************************************************
! ******************************************************
! ******************************************************
! 
!  itime		relative time in integer (old it)
!  dtime		relative time in double precision (in fem files)
!  atime		absolute time (seconds from 1.1.1)
! 
!  atime0	absolute time for fem time 0
! 
!  if no date/time is available then:
! 
! 		itime == dtime == atime
! 		atime0 = 0
! 
! ******************************************************
! ******************************************************
! ******************************************************

        subroutine ptime_init

!  initializes ptime

	use timlim

        implicit none

	integer ihigh
	parameter(ihigh=1000000000)
	double precision ahigh
	parameter(ahigh=4000.*365.*86400.)

!  the first 4 variables might be useless

        itmin = -ihigh
        itmax =  ihigh
	idto  = 0
	itact = itmin

!  these are still used

        itfreq = 1
	nrec = 0

!  these are the new used variables

	atime0 = 0.
	atimeact = 0.
	atimemin = -ahigh
	atimemax =  ahigh

        end

! ******************************************************

        subroutine ptime_min_max

!  sets time limits

	use timlim

        implicit none

	integer iauto
	integer itanf,itend
	double precision atanf,atend

	double precision dgetpar

	iauto = nint(dgetpar('iauto'))
	if( iauto .le. 0 ) then
	  stop 'error stop  ptime_min_max: not ready for iauto=0'
	end if

	itanf = nint(dgetpar('itanf'))
	itend = nint(dgetpar('itend'))

	if( itanf .ne. -1 ) then
	  itmin = itanf
	  call ptime_i2a(itmin,atimemin)
	end if
	if( itend .ne. -1 ) then
	  itmax = itend
	  call ptime_i2a(itmax,atimemax)
	end if

	!write(6,*) '++++++++++++++++++++++++++++++++++++++++++'
	!write(6,*) atime0
	!write(6,*) itmin,itmax,atimemin,atimemax
	!write(6,*) '++++++++++++++++++++++++++++++++++++++++++'

	itfreq = nint(dgetpar('nout'))
        if( itfreq .eq. 0 ) itfreq = 1

	atanf = dgetpar('atanf')
	atend = dgetpar('atend')

	if( atanf .gt. 0. ) atimemin = atanf
	if( atend .gt. 0. ) atimemax = atend

	write(6,*)
	write(6,*) 'Using time parameters: ',itmin,itmax,itfreq
	write(6,*) 'absolute time: ',atimemin,atimemax
	write(6,*)

        end

! ******************************************************

	subroutine ptime_info

	use timlim

	implicit none

	write(6,*) '---------- ptime_info -----------'
	write(6,*) 'itmin,itmax: ',itmin,itmax
	write(6,*) 'itact: ',itact
	write(6,*) 'itfreq: ',itfreq
	write(6,*) 'nrec: ',nrec
	write(6,*) 
	write(6,*) 'atime0: ',atime0
	write(6,*) 'atimeact: ',atimeact
	write(6,*) 'atmin,atmax: ',atimemin,atimemax
	write(6,*) '---------------------------------'

	end

! ******************************************************
! ******************************************************
! ******************************************************

	subroutine ptime_set_date_time(date,time)

	use timlim

	implicit none

	integer date,time

	atime0 = 0.

	if( date > 0 ) then
          call dtsini(date,time)
	  call dts_to_abs_time(date,time,atime0)
	end if

	end

! ******************************************************
! ******************************************************
! ******************************************************

	subroutine ptime_set_itime(it)

	use timlim

	implicit none

	integer it

	atimeact = atime0 + it

	end

! ******************************************************

	subroutine ptime_set_atime(atime)

	use timlim

	implicit none

	double precision atime

	atimeact = atime

	end

! ******************************************************

	subroutine ptime_set_dtime(dtime)

	use timlim

	implicit none

	double precision dtime

	atimeact = atime0 + dtime

	end

! ******************************************************

	subroutine ptime_get_itime(it)

	use timlim

	implicit none

	integer it

	it = nint( atimeact - atime0 )

	end

! ******************************************************

	subroutine ptime_get_atime(atime)

	use timlim

	implicit none

	double precision atime

	atime = atimeact

	end

! ******************************************************

	subroutine ptime_get_dtime(dtime)

	use timlim

	implicit none

	double precision dtime

	dtime = atimeact - atime0

	end

! ******************************************************

	subroutine ptime_i2a(it,atime)

	use timlim

	implicit none

	integer it
	double precision atime

	atime = atime0 + it

	end

! ******************************************************
! ******************************************************
! ******************************************************

        function ptime_ok()

!  is time ok?

	use timlim

        implicit none

	logical ptime_ok

	integer it
        double precision, save :: atimeold = -1
	logical, save :: bdebug = .false.
	character*20 line

	ptime_ok = .false.

	if( atimeact < atimemin ) return
	if( atimeact > atimemax ) return
	
	if( atimeact .ne. atimeold ) then !increase for new time
	  nrec = nrec + 1
	  atimeold = atimeact
	end if

	it = atimeact - atime0

	if( itfreq .gt. 0 .and. mod(nrec,itfreq) .eq. 0 ) then
	  ptime_ok = .true.
	else if( itfreq .eq. 0 ) then
	  ptime_ok = .true.
	else if( itfreq .lt. 0 .and. mod(nrec-1,-itfreq) .eq. 0 ) then
	  ptime_ok = .true.
	end if

	if( bdebug ) then
	  call dts_format_abs_time(atimeact,line)
	  write(6,*) 'testing: ',line,'  ',ptime_ok
	  write(6,*) 'ptime_ok: ',nrec,itfreq,it,atimeact
	end if

	end

! ******************************************************

        function ptime_end()

!  is time over max limit?

	use timlim

        implicit none

	logical ptime_end

	ptime_end = atimeact .gt. atimemax

	end

! ******************************************************
! ******************************************************
! ******************************************************

