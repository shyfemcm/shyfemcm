
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
! 03.09.2024	ggu	more on mpi-tvd
! 10.09.2024	ggu	data structures localized
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

        integer, private, save :: nel_tvd_mpi = 0

        integer, save :: itvd_type = 0

        real, allocatable, save :: xtvdmpiup(:,:,:)        !x-coordinate
        real, allocatable, save :: ytvdmpiup(:,:,:)        !y-coordinate
        integer, allocatable, save :: ietvdmpiup(:,:,:)    !internal elem number
        integer, allocatable, save :: ieetvdmpiup(:,:,:)   !external elem number
        integer, allocatable, save :: iatvdmpiup(:,:,:)    !domain of element

!==================================================================
        contains
!==================================================================

        subroutine mod_tvd_mpi_init(nel)

        integer nel

        if( nel == nel_tvd_mpi ) return

        if( nel_tvd_mpi > 0 ) then
          deallocate(xtvdmpiup)
          deallocate(ytvdmpiup)
          deallocate(ietvdmpiup)
          deallocate(ieetvdmpiup)
          deallocate(iatvdmpiup)
        end if

        nel_tvd_mpi = nel

        if( nel == 0 ) return

        allocate(xtvdmpiup(3,3,nel))
        allocate(ytvdmpiup(3,3,nel))
        allocate(ietvdmpiup(3,3,nel))
        allocate(ieetvdmpiup(3,3,nel))
        allocate(iatvdmpiup(3,3,nel))

	xtvdmpiup = 0.
	ytvdmpiup = 0.
	ietvdmpiup = 0
	ieetvdmpiup = 0
	iatvdmpiup = -1

        end subroutine mod_tvd_mpi_init

!==================================================================
        end module shympi_tvd
!==================================================================

	subroutine tvd_mpi_handle

	use basin
	use shympi
	use shympi_tvd
	!use mod_tvd

	implicit none

	logical bdebug,bfound
	integer ie,ii,j,i,ie_new,ix,n
	integer my_ia,ia,id,inlist,nfound,nchanged
	integer ie_ext,iee_ext,iee_old
	integer ilist,maxlist,maxnlist
	integer ia_found,ia_needed,ie_local,iu,iudb
	integer n_tvd_r,n_tvd_s
	integer, save :: nlist = 10
	real x,y
	real xmin,ymin,xmax,ymax

	real, allocatable :: raux(:,:)
	real, allocatable :: rlist(:,:)
	real, allocatable :: newlist(:,:)
	real, allocatable :: rlists(:,:,:)
	real, allocatable :: newlists(:,:,:)
	integer, allocatable :: count(:)
	integer, allocatable :: index(:)
	integer, allocatable :: n_tvd_receive(:)
	integer, allocatable :: n_tvd_send(:)
	integer, allocatable :: tvd_receive(:,:)
	integer, allocatable :: tvd_send(:,:)

!----------------------------------------------------------
! starting collection of tvd information
!----------------------------------------------------------

	write(6,*) 'configuring tvd with mpi ',my_id

	my_ia = my_id + 1
	iu = 330 + my_ia
	iudb = 340 + my_ia

	ilist = 0
	allocate(raux(nlist,nel*6))
	raux = 0.

!----------------------------------------------------------
! setting up initial list with all information of this domain
!----------------------------------------------------------

	if( 7 > nlist ) goto 98

	do ie=1,nel
	  ie_ext = ipev(ie)
	  do ii=1,3
	    do j=1,3
	      bfound = .false.
	      if( ii == j ) cycle
	      if( ietvdmpiup(j,ii,ie) >= 0 ) then	!we have to check all
		bfound = .true.
		ilist = ilist + 1
		raux(1,ilist) = xtvdmpiup(j,ii,ie)
		raux(2,ilist) = ytvdmpiup(j,ii,ie)
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
! find biggest list information and exchange
!----------------------------------------------------------

	write(6,*) 'exchanging list information ',my_id

	maxlist = shympi_max(ilist)

	write(6,*) 'array sizes: ',ilist,maxlist,my_id
	write(iudb,*) 'array sizes: ',ilist,maxlist,my_ia

	allocate(rlist(nlist,maxlist))
	allocate(rlists(nlist,maxlist,n_threads))
	rlist = 0.

	rlist(:,1:ilist) = raux(:,1:ilist)
	call shympi_gather(nlist,rlist,rlists)

	write(6,*) 'all information collected ',my_id

	call shympi_barrier

!----------------------------------------------------------
! look for missing information needed in other domains
!----------------------------------------------------------

	if( 10 > nlist ) goto 98

	call bas_get_minmax(xmin,ymin,xmax,ymax)

	nfound = 0

	deallocate(raux)
	allocate(raux(nlist,maxlist*n_threads))
	raux = 0.
	inlist = 0

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
	      rlists(9,i,ia) = ie		!internal element number
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

	write(6,*) 'new list info ',inlist,maxnlist,my_id
	write(iudb,*) 'new list info ',inlist,maxnlist,my_ia
	
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
	      iee_old = ieetvdmpiup(j,ii,ie)
	      if( iee_ext > iee_old ) then	!always take highest index
	        ieetvdmpiup(j,ii,ie) = iee_ext
	        ietvdmpiup(j,ii,ie) = ie
	        iatvdmpiup(j,ii,ie) = ia_found
	      end if
	     end if
	    end if
	  end do
	end do

	call shympi_barrier
	write(6,*) 'all new information elaborated ',nchanged,my_id

!----------------------------------------------------------
! prepare data exchange
!----------------------------------------------------------

	deallocate(raux)
	allocate(raux(nlist,9*nel))
	raux = 0.

	ilist = 0
	do ie=1,nel
	  do ii=1,3
	    do j=1,3
	        if( ii == j ) cycle
	        ia = iatvdmpiup(j,ii,ie)		!area where point found
		if( ia == my_ia ) cycle
		ilist = ilist + 1
		raux(1,ilist) = xtvdmpiup(j,ii,ie)
		raux(2,ilist) = ytvdmpiup(j,ii,ie)
		raux(3,ilist) = my_ia		!needed in this domain
		raux(4,ilist) = j
		raux(5,ilist) = ii
		raux(6,ilist) = ie		!internal local element
		raux(7,ilist) = ipev(ie)	!external local element
	        raux(8,ilist) = iatvdmpiup(j,ii,ie)	!area where point found
	        raux(9,ilist) = ietvdmpiup(j,ii,ie)	!internal remote element
	        raux(10,ilist) = ieetvdmpiup(j,ii,ie)!external remote element
	    end do
	  end do
	end do
	  
	!------------------------------------------------
	! prepare count and index
	!------------------------------------------------

	allocate(count(n_threads))
	allocate(index(n_threads))
	count = 0

	do i=1,ilist
	  ia = nint(raux(8,i))
	  count(ia) = count(ia) + 1
	end do
	if( count(my_ia) /= 0 ) stop 'error stop tvd_handle: internal (5)'

	ix = 0
	do ia=1,n_threads
	  index(ia) = ix
	  ix = ix + count(ia)
	end do
	if( ix /= ilist ) stop 'error stop tvd_handle: internal (6)'

	maxlist = shympi_max(ilist)

	write(6,*) 'new array sizes: ',ilist,maxlist,my_id
	write(iudb,*) 'new array sizes: ',ilist,maxlist,my_id
	write(iudb,*) count
	write(iudb,*) index


	deallocate(rlist)
	deallocate(rlists)
	allocate(rlist(nlist,maxlist))
	allocate(rlists(nlist,maxlist,n_threads))

	!------------------------------------------------
	! sort raux and copy to rlist
	! after this in rlist are the same entries as in raux, but sorted by ia
	!------------------------------------------------

	if( ilist > 0 ) then
	write(iudb,*) 'not sorted list: ',my_ia,ilist,maxlist
	write(iudb,*) '    i    ia     j    ii    ie   iee   iaf    ie   iee'
	do i=1,ilist
	  write(iudb,'(10i6)') i,nint(raux(3:10,i))
	end do
	end if

	rlist = 0.

	do i=1,ilist
	  ia = nint(raux(8,i))
	  ix = index(ia) + 1
	  rlist(:,ix) = raux(:,i)
	  index(ia) = ix
	end do

	if( ilist > 0 ) then
	write(iudb,*) 'sorted list: ',my_ia,ilist,maxlist
	write(iudb,*) '    i    ia     j    ii    ie   iee   iaf    ie   iee'
	do i=1,ilist
	  write(iudb,'(10i6)') i,nint(rlist(3:10,i))
	end do
	end if

	call shympi_gather(nlist,rlist,rlists)

	!------------------------------------------------
	! some checks
	!------------------------------------------------

	do ia=1,n_threads
	  do i=1,maxlist
	    ie = nint(rlists(9,i,ia))
	    if( ie == 0 ) exit			!end of list
	    ia_needed = nint(rlists(3,i,ia))
	    if( ia_needed /= ia ) then
	      write(6,*) 'ia and ia_needed are different'
	      write(6,*) ia,ia_needed
	      stop 'error stop tvd_handle: ia/=ia_needed'
	    end if
	  end do
	end do

	!------------------------------------------------
	! compute number of points to receive or send
	!------------------------------------------------

	allocate(n_tvd_receive(n_threads))
	allocate(n_tvd_send(n_threads))

	n_tvd_receive = 0
	n_tvd_send = 0

	do ia=1,n_threads
	  if( ia == my_ia ) cycle
	  do i=1,maxlist
	    ie = nint(rlists(9,i,ia))
	    if( ie == 0 ) exit			!end of list
	    ia_found = nint(rlists(8,i,ia))	!found in this domain
	    if( ia_found == my_ia ) then
	      n_tvd_send(ia) = n_tvd_send(ia) + 1
	    else
	      n_tvd_receive(ia) = n_tvd_receive(ia) + 1
	    end if
	  end do
	end do

	!if( any(n_tvd_receive/=count) ) then
	!  write(6,*) 'ia: ',my_ia
	!  write(6,*) 'receive: ',n_tvd_receive
	!  write(6,*) 'send: ',n_tvd_send
	!  write(6,*) 'count: ',count
	!  stop 'error stop tvd_handle: internal (9)'
	!end if

	n_tvd_r = maxval(n_tvd_receive)
	n_tvd_s = maxval(n_tvd_send)
	write(iudb,*) 'send/receive index start'
	write(iudb,*) n_tvd_r,n_tvd_s
	write(6,*) 'ggguuu: ',my_ia,n_tvd_r,n_tvd_s

	allocate(tvd_receive(n_tvd_r,n_threads))
	allocate(tvd_send(n_tvd_s,n_threads))
	n_tvd_receive = 0
	n_tvd_send = 0
	tvd_receive = 0
	tvd_send = 0

	do ia=1,n_threads
	  do i=1,maxlist
	    ie = nint(rlists(9,i,ia))
	    if( ie == 0 ) exit			!end of list
	    ia_found = nint(rlists(8,i,ia))	!found in this domain
	    if( ia == my_ia ) then
	      n = n_tvd_receive(ia) + 1
	      tvd_receive(n,ia) = i
	      n_tvd_receive(ia) = n
	    else
	      n = n_tvd_send(ia) + 1
	      tvd_send(n,ia) = i
	      n_tvd_send(ia) = n
	    end if
	  end do
	end do

	write(iudb,*) 'receive: ',n_tvd_r
	do ia=1,n_threads
	  n = n_tvd_receive(ia)
	  write(iudb,*) ia,n
	  if( n > 0 ) write(iudb,*) tvd_receive(1:n,ia)
	end do
	write(iudb,*) 'send: ',n_tvd_s
	do ia=1,n_threads
	  n = n_tvd_send(ia)
	  write(iudb,*) ia,n
	  if( n > 0 ) write(iudb,*) tvd_send(1:n,ia)
	end do

	write(6,*) 'all new information collected ',my_id

	call shympi_barrier

!----------------------------------------------------------
! write out debug information
!----------------------------------------------------------

	call write_tvd_debug(nel)

!----------------------------------------------------------
! all info exchanged
!----------------------------------------------------------

	call flush(iu)
	call flush(iudb)
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
	use shympi_tvd

	implicit none

	integer nel

	integer ie,n,j,ii,iu,ie_ext
	integer, allocatable :: iedebug(:,:)
	integer, allocatable :: ieout(:,:)
	integer, allocatable :: iadebug(:,:)
	integer, allocatable :: iaout(:,:)
	character*80 file

	iu = 583
	file='tvd_debug.txt'

	call shympi_barrier
	call flush(6)

	write(6,*) 'collecting tvd debug information...',my_id

	allocate(iedebug(9,nel))
	allocate(ieout(9,nel_global))
	allocate(iadebug(9,nel))
	allocate(iaout(9,nel_global))
	iedebug = 0
	iadebug = 0
	ieout = 0
	iaout = 0

	do ie=1,nel
	  n = 0
	  do ii=1,3
	    do j=1,3
	      n = n + 1
	      iedebug(n,ie) = ieetvdmpiup(j,ii,ie)
	    end do
	  end do
	  iadebug(:,ie) = reshape(ieetvdmpiup(:,:,ie),(/9/))
	  if( any( iedebug(:,ie) /= iadebug(:,ie) ) ) then
	    stop 'error stop write_tvd_debug: reshape not working'
	  end if
	  iadebug(:,ie) = reshape(iatvdmpiup(:,:,ie),(/9/))
	end do

	call shympi_l2g_array(9,iedebug,ieout)
	call shympi_l2g_array(9,iadebug,iaout)

	write(6,*) 'finished collecting tvd debug information...',my_id

	if( shympi_is_master() ) then

	open(iu,file=file,status='unknown',form='formatted')
	write(iu,*) 'debug information'
	do ie=1,nel_global
	  ie_ext = ip_ext_elem(ie)
	  write(iu,'(2i8)') ie,ie_ext
	  write(iu,'(9i7)') ieout(:,ie)
	  !write(iu,'(9i7)') iaout(:,ie)
	end do
	close(iu)
	write(6,*) 'tvd debug information written to file ',trim(file)

	end if

	call shympi_barrier
	flush(6)

	end

!******************************************************************

