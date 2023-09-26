
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-1999,2001,2003,2007,2010-2011  Georg Umgiesser
!    Copyright (C) 2015-2017,2019  Georg Umgiesser
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

! revision log :
!
! 18.11.1998	ggu	check dimensions with dimnos
! 06.04.1999	ggu	some cosmetic changes
! 03.12.2001	ggu	some extra output -> place of min/max
! 09.12.2003	ggu	check for NaN introduced
! 07.03.2007	ggu	easier call
! 23.03.2010	ggu	changed v6.1.1
! 22.11.2011	ggu	changed VERS_6_1_37
! 23.04.2015	ggu	for new version 5
! 21.05.2015	ggu	changed VERS_7_1_11
! 20.07.2015	ggu	changed VERS_7_1_81
! 27.06.2016	ggu	changed VERS_7_5_16
! 05.12.2017	ggu	changed VERS_7_5_39
! 16.02.2019	ggu	changed VERS_7_5_60
! 25.06.2021	ggu	new version to read shy lgr files
!
!**************************************************************

	program lgrinf

! reads lgr file

	implicit none

	integer nread,nin,i,it,nc
	integer mtype,nvers,ftype,iwhat,ncust
	integer nbdy,nn,nout,ierr
	integer id,ie,ies,lmax,lb
	real x,y,z,xs,ys
	integer id_special
	character*80 file
	double precision atime

	integer iapini
	integer ifem_open_file

!--------------------------------------------------------------

	ierr = 0
	nread = 0
	id_special = 4
	id_special = 1
	id_special = 0

!--------------------------------------------------------------
! open simulation
!--------------------------------------------------------------

	nin = 1
	nc = command_argument_count()
        if( nc == 0 ) stop 'no file given'
        call get_command_argument(1,file)
	open(nin,file=file,status='old',form='unformatted')

	read(nin) mtype,nvers
	write(6,*) 'mtype,nvers: ',mtype,nvers
	if( mtype /= 1617 ) stop 'error stop: mtype'
	if( nvers < 6 ) stop 'error stop: old version'
	read(nin) ftype
	write(6,*) 'ftype: ',ftype
	if( ftype /= 3 ) stop 'error stop: wrong type'

	call skip_records(nin,5,ierr)	!shy header
	if( ierr /= 0 ) goto 99
	call skip_records(nin,12,ierr)	!shy second header
	if( ierr /= 0 ) goto 99

	read(nin) ncust
	write(6,*) 'ncust: ',ncust

!--------------------------------------------------------------
! loop on data
!--------------------------------------------------------------

	do

	   read(nin,end=100) atime,nn,iwhat
	   write(6,*) atime,nn,iwhat

	   nread = nread + 1

	   do i=1,nn
	     call skip_records(nin,3,ierr)
	     if( ierr /= 0 ) goto 100
	   end do

	end do	!do while

!--------------------------------------------------------------
! end of loop on data
!--------------------------------------------------------------

  100	continue

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

	if( ierr /= 0 ) then
	  write(6,*) '*** reading problems at particle ',i
	  write(6,*)
	end if

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

   99	continue
	stop 'error stop: reading shy header'
	end

!***************************************************************

	subroutine skip_records(iu,n,ierr)

	implicit none

	integer iu,n,ierr

	integer i,ios

	ierr = 1

	do i=1,n
	  read(iu,iostat=ios)
	  if( ios /= 0 ) return
	end do

	ierr = 0

	end

!***************************************************************

