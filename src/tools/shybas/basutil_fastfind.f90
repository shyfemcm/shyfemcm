
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

! test routines for fast_find and quad_tree routines

! revision log :
!
! 01.10.2024	ggu	written from scratch
! 04.10.2024	ggu	also adapted to quad_tree routines

!**************************************************************************

	subroutine test_fast_find_correctness

	use basin
	use mod_fast_find
	use mod_quad_tree

	implicit none

	logical bnormal
	integer nsize,iemax
	integer nrun,i
	integer ief,ien,ieq
	integer nw,perc
	real xmin,ymin,xmax,ymax
	real x,y
	real rx,ry

	bnormal = .false.
	!bnormal = .true.
	iemax = 100
	nsize = 10
	nrun = 10000
	nw = 1000 / 100

	call fast_find_initialize(nsize)
	call quad_tree_initialize(iemax)

	call bas_get_minmax(xmin,ymin,xmax,ymax)

	do i=1,nrun
	  !write(6,*) 'testing point ',i
	  call random_number(rx)
	  call random_number(ry)
	  x = xmin + (xmax-xmin)*rx
	  y = ymin + (ymax-ymin)*ry
	  call fast_find_search(x,y,ief)
	  ieq = ief
	  call quad_tree_search(x,y,ieq)
	  ien = ief
          if( bnormal ) call find_unique_element(x,y,ien)
	  perc = 100*i/nrun
	  if( mod(i,nw) == 0 ) write(6,*) perc,x,y,ien
	  if( ief /= ien .or. ieq /= ien ) then
	    write(6,*) i,x,y,ief,ien,ieq
	    stop 'error stop'
	  end if
	end do

	call fast_find_finalize
	call quad_tree_finalize

	end

!**************************************************************************

	subroutine test_fast_find_speed

	use basin
	use mod_fast_find
	use mod_quad_tree

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
	nrun = 100000
	nrun = 10000
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

	subroutine test_quad_tree_optimize

	use basin
	use mod_fast_find
	use mod_quad_tree

	implicit none

	integer n,nsize,iemax
	integer nrun,i
	integer ief,ien
	integer nw,perc
	real time1,time2
	real xmin,ymin,xmax,ymax
	real x,y
	real rx,ry
	real, allocatable :: xv(:),yv(:)

	nrun = 100000
	nw = 1000 / 100

	write(6,*) 'running test_quad_tree_optimize: ',nrun

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

	do n=50,10,-1

	iemax = n
	call quad_tree_initialize(iemax)

	call cpu_time(time1)
	do i=1,nrun
	  call quad_tree_search(xv(i),yv(i),ief)
	  perc = 100*i/nrun
	end do
	call cpu_time(time2)
	write(6,*) 'optimize: ',n,time2-time1
	write(67,*) n,time2-time1

	call quad_tree_finalize

	end do

	end

!**************************************************************************

	subroutine test_quad_tree

	use mod_quad_tree

	implicit none

	write(6,*) 'starting quad_tree_test'
	call quad_tree_initialize
	call quad_tree_plot
	call quad_tree_finalize
	write(6,*) 'finished quad_tree_test'

	end

!**************************************************************************

	subroutine test_fast_find

	implicit none

	call test_fast_find_correctness
	call test_quad_tree
	call test_quad_tree_optimize
	!call test_fast_find_optimize

	!call test_fast_find_speed
	!call test_fast_find_optimize

	end

!**************************************************************************

