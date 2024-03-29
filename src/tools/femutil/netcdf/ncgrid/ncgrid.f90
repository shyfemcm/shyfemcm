
!--------------------------------------------------------------------------
!
!    Copyright (C) 2023  Georg Umgiesser
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

! proof of concept
!
! revision log :
!
! 03.03.2023    ggu     created

!--------------------------------------------------------------------------

	program ncgrid

	implicit none

	integer nx,ny
	integer ix,iy
	integer n,nl
	double precision x,y
	double precision, allocatable :: lon(:,:),lat(:,:)
	integer, allocatable :: node(:,:)

	open(1,file='xdata.txt',status='old',form='formatted')
	read(1,*) nx,ny
	allocate(lon(nx,ny))
	read(1,*) lon
	close(1)

	open(1,file='ydata.txt',status='old',form='formatted')
	read(1,*) nx,ny
	allocate(lat(nx,ny))
	read(1,*) lat
	close(1)

	write(6,*) 'finished reading...'

	allocate(node(nx,ny))

	open(1,file='ncgrid.grd',status='unknown',form='formatted')

! write nodes

	n = 0
	do iy=1,ny
	  do ix=1,nx
	    n = n + 1
	    node(ix,iy) = n
	    x = lon(ix,iy)
	    y = lat(ix,iy)
	    write(1,1000) 1,n,0,x,y
	  end do
	end do

! write lines

	nl = 0
	do iy=1,ny
	  nl = nl + 1
	  write(1,2000) 3,nl,0,nx
	  write(1,3000) node(:,iy)
	end do
	  
	do ix=1,nx
	  nl = nl + 1
	  write(1,2000) 3,nl,0,ny
	  write(1,3000) node(ix,:)
	end do
	  
! finished

	close(1)
	write(6,*) 'data writte: ',n

	stop
 1000	format(i1,i10,i3,2f16.8)
 2000	format(i1,i10,i3,i5)
 3000	format(10i7)
	end

