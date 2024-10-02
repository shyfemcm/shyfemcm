
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

! revision log :
!
! 01.10.2024	ggu	written from scratch

!**************************************************************************

	subroutine test_fast_find

	use basin
	use mod_fast_find

	implicit none

	integer nsize
	integer nrun,i
	integer ief,ien
	real xmin,ymin,xmax,ymax
	real x,y
	real rx,ry

	nsize = 10
	nrun = 100

	call fast_find_initialize(nsize)

	call bas_get_minmax(xmin,ymin,xmax,ymax)

	call random_number(rx)
	call random_number(ry)
	x = xmin + (xmax-xmin)*rx
	y = ymin + (ymax-ymin)*ry

	do i=1,nrun
	  call fast_find_search(x,y,ief)
          call find_unique_element(x,y,ien)
	  write(6,*) x,y,ief,ien
	  if( ief /= ien ) stop 'error stop'
	end do

	end

!**************************************************************************

