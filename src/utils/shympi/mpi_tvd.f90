
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
! 28.08.2024	ggu	first part finding elements finished
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

	logical bsphe,bdebug,bfound
	integer isphe
	integer ie,ii,j,i,ie_new
	integer my_ia,ia,id,inlist,nfound,nchanged
	integer ie_ext,iee_ext,iee_old
	integer ilist,maxlist,maxnlist
	integer ia_found,ia_needed,ie_local,iu,iudb
	integer, save :: nlist = 10
	real x,y
	real xmin,ymin,xmax,ymax

	real, allocatable :: raux(:,:)
	real, allocatable :: rlist(:,:)
	real, allocatable :: newlist(:,:)
	real, allocatable :: rlists(:,:,:)
	real, allocatable :: newlists(:,:,:)

!----------------------------------------------------------
! starting collection of tvd information
!----------------------------------------------------------

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
	my_ia = my_id + 1

!----------------------------------------------------------
! setting up initial list with all information of this domain
!----------------------------------------------------------

	iudb = 340 + my_ia

	if( 7 > nlist ) goto 98

	do ie=1,nel
	  ie_ext = ipev(ie)
	  do ii=1,3
	    do j=1,3
	      bfound = .false.
	      if( ii == j ) cycle
	      if( ietvdup(j,ii,ie) >= 0 ) then	!we have to check all
		bfound = .true.
		ilist = ilist + 1
		raux(1,ilist) = tvdupx(j,ii,ie)
		raux(2,ilist) = tvdupy(j,ii,ie)
		raux(3,ilist) = my_ia		!needed in this domain
		raux(4,ilist) = j
		raux(5,ilist) = ii
		raux(6,ilist) = ie		!internal local element
		raux(7,ilist) = ipev(ie)	!external local element
	      end if
	    end do
	  end do
	end do

	call shympi_barrier

!----------------------------------------------------------
! find biggest list information
!----------------------------------------------------------

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

!----------------------------------------------------------
! exchange list information between domains
!----------------------------------------------------------

	call shympi_gather(nlist,rlist,rlists)

	write(6,*) 'all information collected ',my_id

	call shympi_barrier

!----------------------------------------------------------
! look for missing information needed in other domains
!----------------------------------------------------------

	inlist = 0
	if( 10 > nlist ) goto 98

	call bas_get_minmax(xmin,ymin,xmax,ymax)

	iu = 330 + my_ia

	nfound = 0

	do ia=1,n_threads
	  if( ia == my_ia ) cycle
	  do i=1,maxlist
	    ie_local = nint(rlists(6,i,ia))
	    ia_needed = nint(rlists(3,i,ia))
	    ie_ext = nint(rlists(7,i,ia))
	    ie = 0
	    iee_ext = 0
	    if( ie_local == 0 ) exit		!end of list
	    x = rlists(1,i,ia)
	    y = rlists(2,i,ia)
	    if( x < xmin .or. x > xmax ) cycle
	    if( y < ymin .or. y > ymax ) cycle
	    call find_unique_element(x,y,ie)
	    if( ie > 0 ) then			!element found
	      nfound = nfound + 1
	      iee_ext = ipev(ie)
	      rlists(8,i,ia) = my_ia		!found for this domain
	      rlists(9,i,ia) = ie		!found in this element
	      rlists(10,i,ia) = iee_ext		!external element number
	      inlist = inlist + 1
	      raux(:,inlist) = rlists(:,i,ia)
	    end if
	  end do
	end do

	call shympi_barrier
	call flush(6)

!----------------------------------------------------------
! in raux are points found in this domain needed in other domains
!----------------------------------------------------------

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

!----------------------------------------------------------
! exchanging info of found points
!----------------------------------------------------------

	write(6,*) 'all new information collected ',my_id

	if( 10 > nlist ) goto 98

	nchanged = 0

	do ia=1,n_threads
	  if( ia == my_ia ) cycle
	  do i=1,maxnlist
	    ie = nint(newlists(9,i,ia))
	    if( ie == 0 ) exit			!end of list
	    ia_found = nint(newlists(8,i,ia))	!found in this domain
	    ia_needed = nint(newlists(3,i,ia))	!needed in this domain
	    ie_ext = nint(newlists(7,i,ia))	!in this element needed
	    iee_ext = nint(newlists(10,i,ia))	!in this element found
	    if( ia_found /= ia_needed ) then
	     if( my_ia == ia_needed ) then
	      nchanged = nchanged + 1
	      j = nint(newlists(4,i,ia))
	      ii = nint(newlists(5,i,ia))
	      ie = nint(newlists(6,i,ia))
	      iee_old = ieetvdup(j,ii,ie)
	      if( iee_ext > iee_old ) then
	        ieetvdup(j,ii,ie) = iee_ext
	      end if
	     end if
	    end if
	  end do
	end do

	call shympi_barrier
	write(6,*) 'all new information elaborated ',my_id

!----------------------------------------------------------
! write out debug information
!----------------------------------------------------------

	!call write_tvd_debug(nel)

!----------------------------------------------------------
! all info exchanged
!----------------------------------------------------------

	call shympi_barrier
	call flush(6)

	!stop

	return
   98	continue
	write(6,*) 'error stop tvd_handle: nlist too small'
	stop
   99	continue
	write(6,*) 'error stop tvd_handle: ia=iafound'
	stop
	end

!******************************************************************

	subroutine write_tvd_debug(nel)

	use mod_tvd
	use shympi

	implicit none

	integer nel

	integer ie,n,j,ii,iu,ie_ext
	integer, allocatable :: iedebug(:,:)
	integer, allocatable :: ieout(:,:)
	character*80 file

	iu = 583
	file='tvd_debug.txt'

	call shympi_barrier
	call flush(6)

	write(6,*) 'collecting tvd debug information...',my_id

	allocate(iedebug(9,nel))
	allocate(ieout(9,nel_global))
	iedebug = 0
	ieout = 0

	do ie=1,nel
	  n = 0
	  do ii=1,3
	    do j=1,3
	      n = n + 1
	      iedebug(n,ie) = ieetvdup(j,ii,ie)
	    end do
	  end do
	end do

	call shympi_l2g_array(9,iedebug,ieout)

	write(6,*) 'finished collecting tvd debug information...',my_id

	if( shympi_is_master() ) then

	open(iu,file=file,status='unknown',form='formatted')
	write(iu,*) 'debug information'
	do ie=1,nel_global
	  ie_ext = ip_ext_elem(ie)
	  write(iu,'(2i8,9i7)') ie,ie_ext,ieout(:,ie)
	end do
	close(iu)
	write(6,*) 'tvd debug information written to file ',trim(file)

	end if

	call shympi_barrier
	flush(6)

	end

!******************************************************************

