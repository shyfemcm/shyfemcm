
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

! tvd routines
!
! revision log :
!
! 02.08.2024	ggu	started mpi tvd
!
!******************************************************************

!==================================================================
        module shympi_tvd
!==================================================================

	logical, parameter :: btripple = .true. !handles tripple points

	integer, parameter :: nexch = 8
	integer, save :: nmax_tripple = 0
	integer, save :: itrtot = -1
	integer, save, allocatable :: ielist(:,:)
	real, save, allocatable :: buffer_tripple_in(:,:)
	real, save, allocatable :: buffer_tripple_out(:,:)
	integer, save, allocatable :: ietrp(:,:)

!==================================================================
        end module shympi_tvd
!==================================================================

	subroutine tvd_handle

	use basin
	use shympi
	use shympi_tvd
	use mod_tvd

	implicit none

	logical bsphe
	integer isphe
	integer ie,ii,j,i
	integer ia,id,inlist
	integer iafound,ianeeded,ie_ext
	integer ilist,maxlist,maxnlist
	integer, save :: nlist = 10
	real x,y
	real xmin,ymin,xmax,ymax

	real, allocatable :: raux(:,:)
	real, allocatable :: rlist(:,:)
	real, allocatable :: newlist(:,:)
	real, allocatable :: rlists(:,:,:)
	real, allocatable :: newlists(:,:,:)

	write(6,*) 'configuring tvd with mpi ',my_id

        call get_coords_ev(isphe)
        bsphe = isphe .eq. 1

	write(6,*) 'looking for elements ',isphe,nel,my_id

        do ie=1,nel
          call tvd_upwind_init(bsphe,ie)
        end do

	write(6,*) 'tvd_upwind_init finished ',my_id

	ilist = 0
	allocate(raux(nlist,nel*6))
	raux = 0.
	ia = my_id + 1

	do ie=1,nel
	  do ii=1,3
	    do j=1,3
	      if( ii == j ) cycle
	      if( ietvdup(j,ii,ie) == 0 ) then
		ilist = ilist + 1
		raux(1,ilist) = tvdupx(j,ii,ie)
		raux(2,ilist) = tvdupy(j,ii,ie)
		raux(3,ilist) = ia		!needed in this domain
		raux(4,ilist) = j
		raux(5,ilist) = ii
		raux(6,ilist) = ie
		raux(7,ilist) = ipev(ie)
	      end if
	    end do
	  end do
	end do

	call shympi_barrier

	write(6,*) 'exchanging list information ',my_id

	maxlist = shympi_max(ilist)

	write(6,*) 'array sizes: ',ilist,maxlist,my_id

	allocate(rlist(nlist,maxlist))
	allocate(rlists(nlist,maxlist,n_threads))
	rlist = 0.

	rlist(:,1:ilist) = raux(:,1:ilist)

	deallocate(raux)
	allocate(raux(nlist,maxlist*n_threads))
	raux = 0.

	inlist = 0

	call shympi_gather(nlist,rlist,rlists)

	write(6,*) 'all information collected ',my_id

	call shympi_barrier

	call bas_get_minmax(xmin,ymin,xmax,ymax)

	do ia=1,n_threads
	  id = ia - 1
	  if( id == my_id ) cycle
	  do i=1,maxlist
	    ie = nint(rlists(6,i,ia))
	    if( ie == 0 ) cycle		!end of list
	    x = rlists(1,i,ia)
	    y = rlists(2,i,ia)
	    if( x < xmin .or. x > xmax ) cycle
	    if( y < ymin .or. y > ymax ) cycle
	    call find_element(x,y,ie)
	    if( ie > 0 ) then		!element found
	      rlists(8,i,ia) = ia	!found in this domain
	      rlists(9,i,ia) = ie	
	      rlists(10,i,ia) = ipev(ie)
	      inlist = inlist + 1
	      raux(:,inlist) = rlists(:,i,ia)
	    end if
	  end do
	end do

	call shympi_barrier
	call flush(6)

	maxnlist = shympi_max(inlist)

	write(6,*) 'new list info stored ',inlist,maxnlist,my_id
	
	call shympi_barrier
	call flush(6)

	allocate(newlist(nlist,maxnlist))
	allocate(newlists(nlist,maxnlist,n_threads))
	newlist = 0.
	newlists = 0.
	newlist(:,1:inlist) = raux(:,1:inlist)

	call shympi_gather(nlist,newlist,newlists)

	write(6,*) 'all new information collected ',my_id

	do ia=1,n_threads
	  id = ia - 1
	  if( id == my_id ) cycle
	  do i=1,maxnlist
	    ie = nint(newlists(9,i,ia))
	    if( ie == 0 ) cycle		!end of list
	    iafound = nint(newlists(8,i,ia))	!found in this domain
	    if( iafound == ia ) goto 99		!cannot be
	    ianeeded = nint(newlists(3,i,ia))	!needed in this domain
	    if( ia == ianeeded ) then
	      j = nint(newlists(4,i,ia))
	      ii = nint(newlists(5,i,ia))
	      ie = nint(newlists(6,i,ia))
	      ie_ext = nint(newlists(10,i,ia))
	    end if
	  end do
	end do

	call shympi_barrier
	write(6,*) 'all new information elaborated ',my_id

	call shympi_barrier
	call flush(6)

	stop

	return
   99	continue
	write(6,*) 'error stop tvd_handle: ia=iafound'
	stop
	end

!******************************************************************

