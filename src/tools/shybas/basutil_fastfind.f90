
!--------------------------------------------------------------------------
!
!    Copyright (C) 2024  Georg Umgiesser
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

! test routines for fast find routines

! revision log :
!
! 01.10.2024	ggu	written from scratch

!**************************************************************************

	subroutine test_fast_find_correctness

	use basin
	use mod_fast_find

	implicit none

	integer nsize
	integer nrun,i
	integer ief,ien
	integer nw,perc
	real xmin,ymin,xmax,ymax
	real x,y
	real rx,ry

	nsize = 10
	nrun = 10000
	nw = 1000 / 100

	call fast_find_initialize(nsize)

	call bas_get_minmax(xmin,ymin,xmax,ymax)

	do i=1,nrun
	  call random_number(rx)
	  call random_number(ry)
	  x = xmin + (xmax-xmin)*rx
	  y = ymin + (ymax-ymin)*ry
	  call fast_find_search(x,y,ief)
          call find_unique_element(x,y,ien)
	  perc = 100*i/nrun
	  if( mod(i,nw) == 0 ) write(6,*) perc,x,y,ief,ien
	  if( ief /= ien ) then
	    write(6,*) i,x,y,ief,ien
	    stop 'error stop'
	  end if
	end do

	call fast_find_finalize

	end

!**************************************************************************

	subroutine test_fast_find_speed

	use basin
	use mod_fast_find

	implicit none

	integer nsize
	integer nrun,i
	integer ief,ien
	integer nw,perc
	real time1,time2
	real xmin,ymin,xmax,ymax
	real x,y
	real rx,ry
	real, allocatable :: xv(:),yv(:)

	nsize = 30
	nrun = 10000
	nw = 1000 / 100

	allocate(xv(nrun))
	allocate(yv(nrun))

	call fast_find_initialize(nsize)

	call bas_get_minmax(xmin,ymin,xmax,ymax)

	do i=1,nrun
	  call random_number(rx)
	  call random_number(ry)
	  x = xmin + (xmax-xmin)*rx
	  y = ymin + (ymax-ymin)*ry
	  xv(i) = x
	  yv(i) = y
	end do

	call cpu_time(time1)
	do i=1,nrun
	  call fast_find_search(xv(i),yv(i),ief)
	  perc = 100*i/nrun
	end do
	call cpu_time(time2)
	write(6,*) 'cpu time fast find: ',time2-time1

	call cpu_time(time1)
	do i=1,nrun
          call find_unique_element(x,y,ien)
	  perc = 100*i/nrun
	end do
	call cpu_time(time2)
	write(6,*) 'cpu time normal find: ',time2-time1

	call fast_find_finalize

	end

!**************************************************************************

	subroutine test_fast_find_optimize

	use basin
	use mod_fast_find

	implicit none

	integer n,nsize
	integer nrun,i
	integer ief,ien
	integer nw,perc
	real time1,time2
	real xmin,ymin,xmax,ymax
	real x,y
	real rx,ry
	real, allocatable :: xv(:),yv(:)

	nrun = 1000000
	nw = 1000 / 100

	write(6,*) 'running test_fast_find_optimize: ',nrun

	allocate(xv(nrun))
	allocate(yv(nrun))

	call bas_get_minmax(xmin,ymin,xmax,ymax)

	do i=1,nrun
	  call random_number(rx)
	  call random_number(ry)
	  x = xmin + (xmax-xmin)*rx
	  y = ymin + (ymax-ymin)*ry
	  xv(i) = x
	  yv(i) = y
	end do

	do n=5,30

	nsize = n
	call fast_find_initialize(nsize)

	call cpu_time(time1)
	do i=1,nrun
	  call fast_find_search(xv(i),yv(i),ief)
	  perc = 100*i/nrun
	end do
	call cpu_time(time2)
	write(6,*) 'optimize: ',n,time2-time1
	write(66,*) n,time2-time1

	call fast_find_finalize

	end do

	end

!**************************************************************************

	subroutine test_fast_find

	implicit none

	!call test_fast_find_correctness
	!call test_fast_find_speed
	call test_fast_find_optimize

	end

!**************************************************************************

