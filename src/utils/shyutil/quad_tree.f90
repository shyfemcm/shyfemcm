
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

! implements a quad tree for fast search of points in element
!
! revision log :
!
! 01.10.2024	ggu	quad_tree started
! 03.10.2024	ggu	quad_tree debugged
! 04.10.2024	ggu	quad_tree finished
! 05.10.2024	ggu	final changes

!**************************************************************************

!==================================================================
        module mod_quad_tree
!==================================================================

        implicit none

	private

        logical, save :: bdebug = .false.
        logical, save :: bverbose = .false.
        logical, save :: bplot = .false.

        type :: entry
	  integer :: nfill
	  integer, allocatable :: ielems(:)
	  real :: xbmin
	  real :: ybmin
	  real :: xbmax
	  real :: ybmax
	  integer :: nw
	  integer :: ne
	  integer :: sw
	  integer :: se
        end type entry

	type(entry), save, allocatable :: pentry(:)

        logical, save :: binit = .false.	!is initialized?
        integer, save :: idmax = 1000		!initial dimension of boxes
        integer, save :: iemax = 0		!maximum elements in one box

        integer, save :: id_divide = 0
        integer, save :: idlast = 0

	real, save, allocatable :: xemin(:)
	real, save, allocatable :: yemin(:)
	real, save, allocatable :: xemax(:)
	real, save, allocatable :: yemax(:)

	!----------------------------------------------------------
	! following are the public calls for quad_tree routines
	!----------------------------------------------------------

	public :: quad_tree_initialize	!call quad_tree_initialize[(iem)]
	public :: quad_tree_finalize	!call quad_tree_finalize
	public :: quad_tree_search	!call quad_tree_search(x,y,ie)
	public :: quad_tree_plot	!call quad_tree_plot

!==================================================================
	contains
!==================================================================

!******************************************************************

	subroutine quad_tree_initialize(iem)

! sets up data structure

	use basin
	use stack

	integer, optional, intent(in) :: iem	!max number of elements in box

	integer id
	integer ie,ii,k
	integer nfill
	integer iemax_proposed
	real xmin,ymin,xmax,ymax
	real x(3),y(3)

!	------------------------------------------------
!	start of routine
!	------------------------------------------------

	if( bverbose ) write(6,*) 'starting seting up quad_tree routine'
 
!	------------------------------------------------
!	determine iemax and if we have to re-initialize
!	------------------------------------------------

	iemax_proposed = 0
	if( present(iem) ) iemax_proposed = iem
	if( iemax_proposed < ngr ) iemax_proposed = 2*ngr

	if( binit ) then
	  if( iemax_proposed == iemax ) return		!already done
	  call quad_tree_finalize
	end if

	iemax = iemax_proposed

!	------------------------------------------------
!	initialize element boxes
!	------------------------------------------------

	allocate(xemin(nel),yemin(nel))
	allocate(xemax(nel),yemax(nel))

	do ie=1,nel
          do ii=1,3
            k = nen3v(ii,ie)
            x(ii) = xgv(k)
            y(ii) = ygv(k)
          end do
          xemin(ie) = min(x(1),x(2),x(3))
          yemin(ie) = min(y(1),y(2),y(3))
          xemax(ie) = max(x(1),x(2),x(3))
          yemax(ie) = max(y(1),y(2),y(3))
	end do

	allocate(pentry(idmax))

!	------------------------------------------------
!	initialize first box
!	------------------------------------------------

	id = 1
	idlast = 1
        call bas_get_minmax(xmin,ymin,xmax,ymax)
	call init_box(id,nel,xmin,ymin,xmax,ymax)
	pentry(id)%nfill = nel
	do ie=1,nel
	  pentry(id)%ielems(ie) = ie
	end do

	call stack_init(id_divide)
	call stack_push(id_divide,id)

!	------------------------------------------------
!	subdivide the boxes
!	------------------------------------------------

	do
	  if( .not. stack_pop(id_divide,id) ) exit
	  !write(6,*) 'dividing box ',id
	  call sub_divide(id)
	end do
	
	if( bdebug ) write(6,*) 'iemax,idmax,idlast: ',iemax,idmax,idlast
	if( bverbose ) write(6,*) 'finished seting up quad_tree routine'

	binit = .true.

!	------------------------------------------------
!	if wanted plot boxes
!	------------------------------------------------

	if( bplot ) call_quad_tree_plot

!	------------------------------------------------
!	end of routine
!	------------------------------------------------

	end

!******************************************************************

	subroutine sub_divide(id)

	use stack

	integer id

	integer nfill
	integer idnw,idne,idsw,idse
	integer ie,i
	real xmin,ymin,xmax,ymax
	real xhalf,yhalf

	nfill = pentry(id)%nfill
	if( nfill <= iemax ) return		!nothing to do

	if( idlast + 4 > idmax ) call pentry_extend

	xmin = pentry(id)%xbmin
	ymin = pentry(id)%ybmin
	xmax = pentry(id)%xbmax
	ymax = pentry(id)%ybmax

	xhalf = xmin + 0.5*(xmax-xmin)
	yhalf = ymin + 0.5*(ymax-ymin)

	idnw = idlast + 1
	idne = idlast + 2
	idsw = idlast + 3
	idse = idlast + 4
	idlast = idlast + 4

	call init_box(idnw,nfill,xmin,yhalf,xhalf,ymax)
	call init_box(idne,nfill,xhalf,yhalf,xmax,ymax)
	call init_box(idsw,nfill,xmin,ymin,xhalf,yhalf)
	call init_box(idse,nfill,xhalf,ymin,xmax,yhalf)

	do i=1,nfill
	  ie = pentry(id)%ielems(i)
	  call insert_in_box(idnw,ie)
	  call insert_in_box(idne,ie)
	  call insert_in_box(idsw,ie)
	  call insert_in_box(idse,ie)
	end do

	call stack_push(id_divide,idnw)
	call stack_push(id_divide,idne)
	call stack_push(id_divide,idsw)
	call stack_push(id_divide,idse)

	pentry(id)%nw = idnw
	pentry(id)%ne = idne
	pentry(id)%sw = idsw
	pentry(id)%se = idse

	pentry(id)%nfill = 0
	deallocate(pentry(id)%ielems)

	end

!******************************************************************

	subroutine pentry_extend

	type(entry), allocatable :: paux(:)
	
	allocate(paux(idmax))
	paux = pentry
	deallocate(pentry)
	allocate(pentry(2*idmax))
	pentry(1:idmax) = paux(1:idmax)
	idmax = 2 * idmax

	if( bdebug ) write(6,*) 'idmax extended: ',idmax

	end

!******************************************************************

	subroutine insert_in_box(idbox,ie)

	integer idbox,ie

	if( xemax(ie) < pentry(idbox)%xbmin ) return
	if( xemin(ie) > pentry(idbox)%xbmax ) return
	if( yemax(ie) < pentry(idbox)%ybmin ) return
	if( yemin(ie) > pentry(idbox)%ybmax ) return

	pentry(idbox)%nfill = pentry(idbox)%nfill + 1
	pentry(idbox)%ielems(pentry(idbox)%nfill) = ie

	end

!******************************************************************

	subroutine init_box(id,nmax,xmin,ymin,xmax,ymax)

	integer id
	integer nmax
	real xmin,ymin,xmax,ymax

	!write(6,*) 'new box ',id

	pentry(id)%nfill = 0
	allocate(pentry(id)%ielems(nmax))

	pentry(id)%xbmin = xmin
	pentry(id)%ybmin = ymin
	pentry(id)%xbmax = xmax
	pentry(id)%ybmax = ymax

	pentry(id)%nw = 0
	pentry(id)%ne = 0
	pentry(id)%sw = 0
	pentry(id)%se = 0

	end

!******************************************************************

	subroutine quad_tree_search(x,y,ie)

	use basin

	real x,y
	integer ie

	integer id,nfill
	integer i,iee,ie_ext,ifound
	integer iloop
	integer ie_found(iemax)
	real xmin,ymin,xmax,ymax
	real xhalf,yhalf

	logical in_element

	if( .not. binit ) then
	  stop 'error stop quad_tree_search: not initialized'
	end if

	id = 1			!istart from root
	iloop = 0
	ie = 0

	do
	  nfill = pentry(id)%nfill
	  !write(6,*) 'searching: ',nfill,id
	  if( nfill /= 0 ) exit

	  xmin = pentry(id)%xbmin
	  ymin = pentry(id)%ybmin
	  xmax = pentry(id)%xbmax
	  ymax = pentry(id)%ybmax

	  if( x < xmin .or. x > xmax ) return
	  if( y < ymin .or. y > ymax ) return

	  xhalf = xmin + 0.5*(xmax-xmin)
	  yhalf = ymin + 0.5*(ymax-ymin)

	  if( x > xhalf ) then
	    if( y > yhalf ) then
	      id = pentry(id)%ne
	    else
	      id = pentry(id)%se
	    end if
	  else
	    if( y > yhalf ) then
	      id = pentry(id)%nw
	    else
	      id = pentry(id)%sw
	    end if
	  end if

	  if( id == 0 ) return

	  iloop = iloop + 1
	  if( iloop > 100 ) stop 'error stop quad_tree_search: iloop too high'
	end do
	
        ifound = 0

        do i=1,nfill
          iee = pentry(id)%ielems(i)
          if( in_element(iee,x,y) ) then
            ifound = ifound + 1
	    if( ifound > iemax ) stop 'error stop quad_tree_search: iemax'
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

	subroutine quad_tree_finalize

! delete data structure

	use stack

	integer id,nfill

	if( .not. binit ) return

	call stack_delete(id_divide)

	do id=1,idlast
	  nfill = pentry(id)%nfill
	  if( nfill > 0 ) deallocate(pentry(id)%ielems)
	end do

	deallocate(pentry)

	deallocate(xemin,xemax,yemin,yemax)

	binit = .false.

	end

!******************************************************************

	subroutine quad_tree_plot

! delete data structure

	integer id,iu,nfill,n,nl
	real xmin,ymin,xmax,ymax

	if( .not. binit ) then
	  stop 'error stop quad_tree_plot: not initialized'
	end if

	iu = 67
	open(iu,file='quad_tree.grd',status='unknown',form='formatted')

	n = 0
	nl = 0

	do id=1,idlast
	  nfill = pentry(id)%nfill
	  if( nfill == 0 ) cycle
	  xmin = pentry(id)%xbmin
	  ymin = pentry(id)%ybmin
	  xmax = pentry(id)%xbmax
	  ymax = pentry(id)%ybmax
	  write(iu,1000) 1,n+1,0,xmin,ymin
	  write(iu,1000) 1,n+2,0,xmin,ymax
	  write(iu,1000) 1,n+3,0,xmax,ymax
	  write(iu,1000) 1,n+4,0,xmax,ymin
	  nl = nl + 1
	  write(iu,2000) 2,nl,0,4,n+1,n+2,n+3,n+4
	  n = n + 4
	end do

	close(iu)

	return
 1000   format(i1,2i8,2e16.8)
 2000   format(i1,7i8)
	end

!******************************************************************

!==================================================================
        end module mod_quad_tree
!==================================================================

