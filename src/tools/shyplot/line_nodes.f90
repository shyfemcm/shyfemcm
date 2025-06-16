
!--------------------------------------------------------------------------
!
!    Copyright (C) 2009-2012,2014-2016,2018-2019  Georg Umgiesser
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

!  routines to find nodes close to line
! 
!  revision log :
! 
!  14.09.2009	ggu	routines written from scratch
!  23.03.2010	ggu	changed v6.1.1
!  01.09.2011	ggu	changed VERS_6_1_32
!  04.11.2011	ggu	changed VERS_6_1_35
!  21.06.2012	ggu	changed VERS_6_1_54
!  05.12.2014	ggu	changed VERS_7_0_8
!  12.12.2014	ggu	changed VERS_7_0_9
!  23.12.2014	ggu	changed VERS_7_0_11
!  19.01.2015	ggu	changed VERS_7_1_3
!  10.07.2015	ggu	changed VERS_7_1_50
!  17.07.2015	ggu	changed VERS_7_1_52
!  17.07.2015	ggu	changed VERS_7_1_80
!  20.07.2015	ggu	changed VERS_7_1_81
!  30.07.2015	ggu	changed VERS_7_1_83
!  17.06.2016	ggu	changed VERS_7_5_15
!  18.12.2018	ggu	changed VERS_7_5_52
!  21.05.2019	ggu	changed VERS_7_5_62
! 
! ************************************************************************
! 
!  given a line defined by 2 or more points (coordinates x,y)
!  the program returns the nodes that closest aproximate the line
! 
!  main routine to call: find_line_nodes
! 
!  needs basin information and ilinkv(1),lenkv(1),linkv(1) data structure
!  to be link with links.f sublnku.f
! 
! ****************************************************************

	subroutine find_line_nodes(nl,x,y,ndim,n,nodes)

!  find nodes close to line
! 
!  given a line defined by 2 or more points (coordinates x,y)
!  the program returns the nodes that closest aproximate the line

	implicit none

	integer nl		!total number of coordinates defining line
	real x(nl), y(nl)	!coordinates of line
	integer ndim		!dimension of nodes()
	integer n		!total number of nodes found (return)
	integer nodes(ndim)	!nodes found (return)

	integer nt,ns,i

	nt = 1

	do i=2,nl
	  nt = nt - 1	!get rid of last node - otherwise we have it twice
	  call find_segment_nodes(x(i-1),y(i-1),x(i),y(i) &
     &					,ndim,ns,nodes(nt+1))
	  nt = nt + ns	!add new nodes found
	end do

	n = nt

	end

! ****************************************************************

	subroutine find_segment_nodes(x1,y1,x2,y2,ndim,n,nodes)

!  find nodes close to one segment of line

	implicit none

	real x1,y1		!coordinates of initial point
	real x2,y2		!coordinates of final point
	integer ndim		!dimension of nodes()
	integer n		!total number of nodes found (return)
	integer nodes(ndim)	!nodes found (return)

	integer k1,k2,ka,kn

	call get_closest_node(x1,y1,k1)
	call get_closest_node(x2,y2,k2)

	call write_node_info(k1)
	call write_node_info(k2)

	ka = k1
	n = 1
	nodes(n) = ka

	do while( ka .ne. k2 )

	  call get_next_node(k1,k2,ka,kn)

	  n = n + 1
	  if( n .gt. ndim ) goto 99
	  nodes(n) = kn
	  ka = kn

	end do

	write(6,*) 'find_segment_nodes finish: ',n

	return
   99	continue
	stop 'error stop find_segment_nodes: ndim'
	end

! ****************************************************************

	subroutine get_next_node(k1,k2,ka,kn)

!  gets next node close to line segment starting from ka

	implicit none

	integer k1,k2		!start/end node of line segment
	integer ka		!last node found close to line segment
	integer kn		!next node found close to line segment (return)

	integer ndim
	parameter(ndim=100)	!must be at least maximum grade
	integer nodes(ndim)
	integer n,i,k
	real dist,distc,t,t0

	call get_nodes_around(ka,ndim,n,nodes)
	
	call get_distance_from_line(k1,k2,ka,dist,t0)
	kn = 0
	distc = 0.

	do i=1,n
	  k = nodes(i)
	  call get_distance_from_line(k1,k2,k,dist,t)
	  if( t .gt. t0 ) then				!only nodes ahead
	    if( kn .eq. 0 .or. dist .lt. distc ) then   !choose closest node
		kn = k
		distc = dist
	    end if
	  end if
	end do

	if( kn .eq. 0 ) then
	  stop 'error stop get_next_node: no node found'
	end if

	end

! ****************************************************************

	subroutine get_distance_from_line(k1,k2,kp,dist,t)

!  computes distance of kp from line given by k1,k2

	use basin

	implicit none

	integer k1,k2		!start/end node of line segment
	integer kp		!node to compute distance to line segment
	real dist		!distance of kp from line segment (return)
	real t			!position of intersection on segment (return)

!  (xs,ys) is intersection point on segment
!  (xa,ya) is vector of segment : (xa,ya) = (x2,y2) - (x1,y1)
!  (xs,ys) = (x1,y1) + t * (xa,ya)

	real x1,y1,x2,y2,xp,yp
	real xa,ya,xs,ys


	x1 = xgv(k1)
	y1 = ygv(k1)
	x2 = xgv(k2)
	y2 = ygv(k2)
	xp = xgv(kp)
	yp = ygv(kp)

	xa = x2 - x1
	ya = y2 - y1

	t = ( xa*(xp-x1) + ya*(yp-y1) ) / ( xa*xa + ya*ya )

	xs = x1 + t * xa
	ys = y1 + t * ya

	dist = sqrt( (xp-xs)**2 + (yp-ys)**2 )

	end

! ****************************************************************
! ****************************************************************
! ****************************************************************

	subroutine get_closest_node(x0,y0,kc)

!  finds closest node to coordinate (x0,y0)

	use basin

	implicit none

	real x0,y0		!coordinates of point
	integer kc		!closest node to point (return)


	integer k
	real dist,distc

	kc = 1
	distc = (xgv(1)-x0)**2 + (ygv(1)-y0)**2

	do k=2,nkn
	  dist = (xgv(k)-x0)**2 + (ygv(k)-y0)**2
	  if( dist .lt. distc ) then
	    distc = dist
	    kc = k
	  end if
	end do

	end

! ****************************************************************

	subroutine write_node_info(k)

!  writes node info for node k

	use basin

	implicit none

	integer k


	write(6,*) 'node = ',k,xgv(k),ygv(k)

	end

! ****************************************************************
! ****************************************************************
! ****************************************************************
! ****************************************************************
! ****************************************************************

	subroutine line_points_test

	use mod_geom
	use basin

	implicit none

	integer ndim
	parameter (ndim=500)

	character*80 basfil
	integer nodes(ndim)
	integer n,nb,nl
	integer nlkdi
	real x(ndim), y(ndim)

	integer iapini

	call shyfem_copyright('line_nodes - find node along line')

	call read_line(ndim,nl,x,y)

	if( iapini(1,0,0,0) .le. 0 ) stop

	call set_geom

	call find_line_nodes(nl,x,y,ndim,n,nodes)

	write(6,*) n,' nodes found in line'
	call write_line(n,nodes)
	write(6,*) 'output written to file 66 and 77'

	end

! ****************************************************************

	subroutine read_memory_basin(basin)

	implicit none

	character*(*) basin

	integer i

	basin = ' '

	open(1,file='.memory',err=99)
	read(1,*) 
	read(1,'(a)') basin 
	close(1)

	do i=len(basin),1,-1
	  if( basin(i:i) .ne. ' ' ) goto 1
	end do
    1	continue
	basin(i+1:) = '.bas'

	write(6,*) 'using basin file: ',basin(1:i+4)
	return
   99	continue
	write(6,*) 'Cannot read memory file .memory'
	stop 'error stop read_memory_basin: memory'
	end

! ****************************************************************

	subroutine write_line(n,nodes)

	use basin

	implicit none

	integer n
	integer nodes(n)

	integer i,k,kext
	real x,y

	integer ipext

! --------------------------------------------------------------
!  write grd file
! --------------------------------------------------------------

        open(66,file='line_nodes.grd',status='unknown',form='formatted')

	write(66,*)
	do i=1,n
	  k = nodes(i)
	  kext = ipext(k)
	  x = xgv(k)
	  y = ygv(k)
	  write(66,1000) 1,kext,0,x,y
	  nodes(i) = kext
	end do

	write(66,*)
	write(66,2000) 3,k,0,n
	write(66,3000) (nodes(i),i=1,n)
	write(66,*)

        close(66)

        write(6,*) 'grd file written to line_nodes.grd'

! --------------------------------------------------------------
!  write txt file
! --------------------------------------------------------------

        open(77,file='line_nodes.txt',status='unknown',form='formatted')

	do i=1,n
	  kext = nodes(i)
	  write(77,*) kext
	end do
	write(77,*) 0

        close(77)

        write(6,*) 'txt file written to line_nodes.txt'

! --------------------------------------------------------------
!  end of routine
! --------------------------------------------------------------

	return
 1000	format(i1,2i10,2f14.4)
 2000	format(i1,3i10)
 3000	format((5i10))
	end

! ****************************************************************

	subroutine read_line(ndim,nl,x,y)

	implicit none

	integer ndim,nl
	real x(ndim), y(ndim)

	integer i,k,ityp
	real xp,yp

	write(6,*) 'reading line from STDIN'

	nl = 0
    1	continue
	  !read(5,*,end=2) k,ityp,xp,yp
	  read(5,*,end=2) xp,yp,ityp
	  nl = nl + 1
	  if( nl .gt. ndim ) stop 'error stop read_line: ndim'
	  x(nl) = xp
	  y(nl) = yp
	  goto 1
    2	continue

	write(6,*) 'points read: ',nl
	do i=1,nl
	  write(6,*) x(i),y(i)
	end do

	end

! ****************************************************************

	program line_points_main
	implicit none
	call line_points_test
	end

! ****************************************************************

