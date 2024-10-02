
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

! implements a fast search method for point in element
!
! revision log :
!
! 01.10.2024	ggu	fast_find started

!**************************************************************************

!==================================================================
        module mod_fast_find
!==================================================================

        implicit none

        integer, private, save :: nsize = 0
        integer, private, save :: idmax = 0
        integer, private, save :: nbx,nby
        real, private, save :: box_size = 0
	real, private, save :: xbmin,ybmin,xbmax,ybmax
	real, private, save :: dbx,dby

	integer, save, allocatable :: n_elem_list(:)
	integer, save, allocatable :: elem_list(:,:)

!==================================================================
	contains
!==================================================================

        subroutine mod_fast_find_init

	write(6,*) 'mod_fast_find_init successfully called '

        end subroutine mod_fast_find_init

!******************************************************************

	subroutine fast_find_initialize(n)

! sets up data structure

	use basin
	use stack

	integer n

	integer nx,ny
	integer ix,iy
	integer ixmin,iymin,ixmax,iymax
	integer id
	integer ie,ii,k
	integer ne,nemax,netot,neaver
	real x(3),y(3)
	real xmin,ymin,xmax,ymax
	real xemin,yemin,xemax,yemax
	real dx,dy
	real, parameter :: eps = 0.01

!	------------------------------------------------
!	find box sizes and limits
!	------------------------------------------------

	call bas_get_minmax(xmin,ymin,xmax,ymax)

	dx = xmax - xmin
	dy = ymax - ymin
	xbmin = xmin - eps*dx
	xbmax = xmax + eps*dx
	ybmin = ymin - eps*dy
	ybmax = ymax + eps*dy
	dbx = dx
	dby = dy

	if( dx > dy ) then
	  box_size = (xbmax-xbmin)/n
	  nbx = n
	  nby = (ybmax-ybmin) / box_size
	else
	  box_size = (ybmax-ybmin)/n
	  nby = n
	  nbx = (xbmax-xbmin) / box_size
	end if

	write(6,*) 'box: ',n,nbx,nby,box_size

!	------------------------------------------------
!	initialize stack for all boxes
!	------------------------------------------------

	id = 0
	do iy=1,ny
	  do ix=1,nx
	    id = id + 1
	    call stack_init(id)
	  end do
	end do
	idmax = id

!	------------------------------------------------
!	insert elements
!	------------------------------------------------

	do ie=1,nel

          do ii=1,3
            k = nen3v(ii,ie)
            x(ii) = xgv(k)
            y(ii) = ygv(k)
          end do

          xemin = min(x(1),x(2),x(3))
          yemin = min(y(1),y(2),y(3))
          xemax = max(x(1),x(2),x(3))
          yemax = max(y(1),y(2),y(3))
	  ixmin = 1 + (xemin-xbmin)/dbx
	  ixmax = 1 + (xemax-xbmin)/dbx
	  iymin = 1 + (yemin-ybmin)/dby
	  iymax = 1 + (yemax-ybmin)/dby

	  do iy=iymin,iymax
	    do ix=ixmin,ixmax
	      id = (iy-1)*nbx + ix
	      if( id < 1 .or. id > idmax ) then
		write(6,*) 'id,idmax: ',id,idmax
		stop 'error stop fast_find_initialize: internal error (1)'
	      end if
	      call stack_push(id,ie)
	    end do
	  end do

        end do

!	------------------------------------------------
!	convert from stack to element list
!	------------------------------------------------

	nemax = 0
	netot = 0
	neaver = 0
	do id=1,idmax
	  ne = stack_fill(id)
	  nemax = max(nemax,ne)
	  netot = netot + ne
	end do
	neaver = netot / idmax

	write(6,*) 'stack ne: ',idmax,nemax,netot,neaver

	allocate(n_elem_list(idmax))
	allocate(elem_list(nemax,idmax))

	do id=1,idmax
	  call stack_get_entries(id,n,elem_list(:,id))
	  n_elem_list(id) = n
	  ne = stack_fill(id)
	  if( ne /= n ) then
	    write(6,*) 'n,ne: ',n,ne
	    stop 'error stop fast_find_initialize: internal error (w)'
	  end if
	end do

	do id=1,idmax
	  call stack_delete(id)
	end do

!	------------------------------------------------
!	end of routine
!	------------------------------------------------

	end

!******************************************************************

	subroutine fast_find_finalize

! delete data structure

	use basin

	idmax = 0
	deallocate(n_elem_list)
	deallocate(elem_list)

	end

!******************************************************************

	subroutine fast_find_search(x,y,ie)

! finds unique element for x,y - if not found returns ie==0

	use basin

	real x,y
	integer ie

	integer ix,iy,id
	integer ifound,n,iee,i,ie_ext
	integer ie_found(100)

	logical in_element

	ix = 1 + (x-xbmin)/dbx
	iy = 1 + (y-ybmin)/dby
	id = (iy-1)*nbx + ix

	if( id < 1 .or. id > idmax ) then
	  write(6,*) 'id,idmax: ',id,idmax
	  stop 'error stop fast_find_search: internal error (1)'
	end if
	
	n = n_elem_list(id)
	ifound = 0

	do i=1,n
	  iee = elem_list(i,id)
	  if( in_element(iee,x,y) ) then
	    ifound = ifound + 1
	    ie_found(ifound) = iee
	  end if
	end do

	if( ifound == 0 ) then
	  ie = 0
	else if( ifound == 1 ) then
	  ie = ie_found(1)
	else
	  ie = 0
	  ie_ext = 0
	  do i=1,ifound
	    iee = ie_found(i)
	    if( ipev(iee) > ie_ext ) then
	      ie = iee
	      ie_ext = ipev(iee)
	    end if
	  end do
	end if

	end

!******************************************************************

!==================================================================
        end module mod_fast_find
!==================================================================

