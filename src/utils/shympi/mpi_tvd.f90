
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

	implicit none

	logical, parameter :: btvd = .true. !handles tripple points

	integer, save :: nmax_tvd = 0
	real, save, allocatable :: buffer_tvd_in(:,:)
	real, save, allocatable :: buffer_tvd_out(:,:)
	real, save, allocatable :: n_index_tvd(:)
	real, save, allocatable :: index_tvd(:,:)

        integer, save :: itvd_type = 0

	integer, parameter :: ind_x = 1
	integer, parameter :: ind_y = 2
	integer, parameter :: ind_ia_needed = 3
	integer, parameter :: ind_j_this = 4
	integer, parameter :: ind_ii_this = 5
	integer, parameter :: ind_ie_this = 6
	integer, parameter :: ind_iee_this = 7
	integer, parameter :: ind_ia_found = 8
	integer, parameter :: ind_ie_remote = 9
	integer, parameter :: ind_iee_remote = 10

	integer, parameter :: ind_k_this = 11
	integer, parameter :: ind_lmax_this = 12
	integer, parameter :: ind_lmax_remote = 13

	integer, parameter :: nlist = 13
	integer, parameter :: nlist_type = 1

	integer, save :: n_tvd_receive_total,n_tvd_send_total
	real, save, allocatable :: values_tvd_local(:,:)
	real, save, allocatable :: values_tvd_remote(:,:)
	integer, save :: n_my_array_list,n_my_receive_list,n_my_send_list
	integer, save, allocatable :: my_receive_list(:),my_send_list(:)
	real, save, allocatable :: my_array_receive_list(:,:)
	real, save, allocatable :: my_array_send_list(:,:)
	real, save, allocatable :: my_array_list(:,:)
	double precision, save, allocatable :: my_xi_list(:,:)

!==================================================================
        contains
!==================================================================

	subroutine jii2k(j,ii,k)

	integer j,ii,k

	integer, parameter :: mat(3,3) = reshape((/0,1,2,3,0,4,5,6,0/),(/3,3/))
	
	k = mat(j,ii)

	end subroutine

!******************************************************************

	subroutine k2jii(k,j,ii)

	integer j,ii,k

	integer, parameter :: mat(2,6) = &
     &		reshape((/2,1,3,1,1,2,3,2,1,3,2,3/),(/2,6/))
	
	j = mat(1,k)
	ii = mat(2,k)

	end subroutine

!==================================================================
        end module shympi_tvd
!==================================================================

	subroutine tvd_mpi_handle

	use basin
	use levels
	use shympi
	use shympi_tvd
	use mod_tvd

	implicit none

	logical bdebug,bfound
	integer ie,ii,j,i,ie_new,ix,n,ind,iaf,ian,k
	integer my_ia,ia,id,inlist,nfound,nchanged
	integer ie_ext,iee_ext,iee_old,ie_int
	integer ilist,maxlist,maxnlist
	integer ia_found,ia_needed,ie_local,iu,iudb
	integer n_tvd_r,n_tvd_s,n_tvd
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
	integer, allocatable :: n_tvd_index(:)
	integer, allocatable :: tvd_receive(:,:)
	integer, allocatable :: tvd_send(:,:)
	integer, allocatable :: tvd_index(:,:)
	integer, allocatable :: matrix(:,:)
	integer, allocatable :: n_tvd_receive_from(:)
	integer, allocatable :: n_tvd_send_to(:)
	integer :: n_tvd_all
	real, allocatable :: tvd_all(:,:)
	integer :: na,nr,ns
	double precision xi(3),xd,yd

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

	do ie=1,nel
	  ie_ext = ipev(ie)
	  do ii=1,3
	    do j=1,3
	      bfound = .false.
	      if( ii == j ) cycle
	      if( ietvdup(j,ii,ie) >= 0 ) then	!we have to check all
		bfound = .true.
		ilist = ilist + 1
		call jii2k(j,ii,k)
		raux(ind_x,ilist) = xtvdup(j,ii,ie)
		raux(ind_y,ilist) = ytvdup(j,ii,ie)
		raux(ind_ia_needed,ilist) = my_ia	!needed in this domain
		raux(ind_j_this,ilist) = j
		raux(ind_ii_this,ilist) = ii
		raux(ind_k_this,ilist) = k
		raux(ind_ie_this,ilist) = ie		!internal local element
		raux(ind_iee_this,ilist) = ipev(ie)	!external local element
		raux(ind_lmax_this,ilist) = ilhv(ie)	!lmax
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
	write(iudb,*) 'debug output from tvd_mpi_handle for area ',my_ia
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

	call bas_get_minmax(xmin,ymin,xmax,ymax)

	nfound = 0

	deallocate(raux)
	allocate(raux(nlist,maxlist*n_threads))
	raux = 0.
	inlist = 0

	do ia=1,n_threads
	  if( ia == my_ia ) cycle
	  do i=1,maxlist
	    ie_local = nint(rlists(ind_ie_this,i,ia))
	    ia_needed = nint(rlists(ind_ia_needed,i,ia))
	    ie_ext = nint(rlists(ind_iee_this,i,ia))
	    ie = 0
	    iee_ext = 0
	    if( ie_local == 0 ) exit		!end of list
	    x = rlists(ind_x,i,ia)
	    y = rlists(ind_y,i,ia)
	    if( x < xmin .or. x > xmax ) cycle
	    if( y < ymin .or. y > ymax ) cycle
	    call find_unique_element(x,y,ie)
	    if( ie > 0 ) then			!element found
	      nfound = nfound + 1
	      iee_ext = ipev(ie)
	      rlists(ind_ia_found,i,ia) = my_ia		!found for this domain
	      rlists(ind_ie_remote,i,ia) = ie		!internal element number
	      rlists(ind_iee_remote,i,ia) = iee_ext	!external element number
	      rlists(ind_lmax_remote,i,ia) = ilhv(ie)	!lmax
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

	nchanged = 0

	do ia=1,n_threads
	  if( ia == my_ia ) cycle
	  do i=1,maxnlist
	    ie = nint(newlists(ind_ie_remote,i,ia))
	    if( ie == 0 ) exit			!end of list
	    ia_found = nint(newlists(ind_ia_found,i,ia))!found in this domain
	    ia_needed = nint(newlists(ind_ia_needed,i,ia))!needed in this domain
	    ie_ext = nint(newlists(ind_iee_this,i,ia))	!in this element needed
	    ie_int = nint(newlists(ind_ie_remote,i,ia))	!in this element found
	    iee_ext = nint(newlists(ind_iee_remote,i,ia))!in this element found
	    if( ia_found /= ia_needed ) then
	     if( my_ia == ia_needed ) then
	      nchanged = nchanged + 1
	      j = nint(newlists(ind_j_this,i,ia))
	      ii = nint(newlists(ind_ii_this,i,ia))
	      ie = nint(newlists(ind_ie_this,i,ia))
	      iee_old = ieetvdup(j,ii,ie)
		call tvd_assert('j',j>=0.and.j<=3)
		call tvd_assert('ii',ii>=0.and.ii<=3)
		call tvd_assert('ie',ie>=0.and.ie<=nel)
	      if( iee_ext > iee_old ) then	!always take highest index
	        ieetvdup(j,ii,ie) = iee_ext
	        ietvdup(j,ii,ie) = ie_int
	        iatvdup(j,ii,ie) = ia_found
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
	        ia = iatvdup(j,ii,ie)		!area where point found
		if( ia == my_ia ) cycle
		ilist = ilist + 1
		call jii2k(j,ii,k)
		raux(ind_x,ilist) = xtvdup(j,ii,ie)
		raux(ind_x,ilist) = ytvdup(j,ii,ie)
		raux(ind_ia_needed,ilist) = my_ia	!needed in this domain
		raux(ind_j_this,ilist) = j
		raux(ind_ii_this,ilist) = ii
		raux(ind_k_this,ilist) = k
		raux(ind_ie_this,ilist) = ie		!internal local element
		raux(ind_iee_this,ilist) = ipev(ie)	!external local element
	        raux(ind_ia_found,ilist) = iatvdup(j,ii,ie)!area where found
	        raux(ind_ie_remote,ilist) = ietvdup(j,ii,ie)!internal remote
	        raux(ind_iee_remote,ilist) = ieetvdup(j,ii,ie)!external remote
	        raux(ind_lmax_this,ilist) = ilhv(ie) 	!external remote
	        raux(ind_lmax_remote,ilist) = ltvdup(j,ii,ie)!external remote
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
	  ia = nint(raux(ind_ia_found,i))
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
	write(iudb,*) 'new array sizes: ',ilist,maxlist,my_ia
	write(iudb,*) 'count: ',count
	write(iudb,*) 'index: ',index

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
	call write_array_header(iudb)
	do i=1,ilist
	  call write_array_debug(iudb,i,raux(:,i))
	end do
	end if

	rlist = 0.

	do i=1,ilist
	  ia = nint(raux(ind_ia_found,i))
	  ix = index(ia) + 1
	  rlist(:,ix) = raux(:,i)
	  index(ia) = ix
	end do

	if( ilist > 0 ) then
	write(iudb,*) 'sorted list: ',my_ia,ilist,maxlist
	call write_array_header(iudb)
	do i=1,ilist
	  call write_array_debug(iudb,i,rlist(:,i))
	end do
	end if

	call shympi_gather(nlist,rlist,rlists)

	!------------------------------------------------
	! some checks
	!------------------------------------------------

	do ia=1,n_threads
	  do i=1,maxlist
	    ie = nint(rlists(ind_iee_remote,i,ia))
	    if( ie == 0 ) exit			!end of list
	    ia_needed = nint(rlists(ind_ia_needed,i,ia))
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
	  !if( ia == my_ia ) cycle
	  do i=1,maxlist
	    ie = nint(rlists(ind_iee_remote,i,ia))
	    if( ie == 0 ) exit			!end of list
	    ia_found = nint(rlists(ind_ia_found,i,ia))	!found in this domain
	    if( ia_found == ia ) cycle
	    if( ia_found == my_ia ) then
	      n_tvd_send(ia) = n_tvd_send(ia) + 1
	    else
	      n_tvd_receive(ia) = n_tvd_receive(ia) + 1
	    end if
	  end do
	end do

	n_tvd_r = maxval(n_tvd_receive)
	n_tvd_s = maxval(n_tvd_send)
	n_tvd = max(n_tvd_r,n_tvd_s)
	write(6,*) 'n_tvd_r/s: ',n_tvd_r,n_tvd_s,n_tvd,my_ia
	write(6,*) 'n_tvd_receive: ',n_tvd_receive,my_ia
	write(6,*) 'n_tvd_send: ',n_tvd_send,my_ia

	allocate(tvd_receive(n_tvd_r,n_threads))
	allocate(tvd_send(n_tvd_s,n_threads))
	allocate(tvd_index(n_tvd,n_threads))
	allocate(tvd_all(nlist,n_tvd*n_threads))
	allocate(n_tvd_index(n_threads))
	allocate(matrix(n_threads,n_threads))
	allocate(n_tvd_send_to(n_threads))
	allocate(n_tvd_receive_from(n_threads))
	n_tvd_receive = 0
	n_tvd_send = 0
	n_tvd_index = 0
	tvd_receive = 0
	tvd_send = 0
	tvd_index = 0
	tvd_all = 0.
	matrix = 0
	n_tvd_send_to = 0
	n_tvd_receive_from = 0
	n_tvd_all = 0
	n_my_array_list = 0
	n_my_receive_list = 0
	n_my_send_list = 0

	do ia=1,n_threads
	  do i=1,maxlist
	    ie = nint(rlists(ind_ie_this,i,ia))
	    if( ie == 0 ) exit			!end of list
	    ia_needed = nint(rlists(ind_ia_needed,i,ia))!needed in this domain
	    if( ia /= ia_needed ) stop 'error stop: ia/=ia_needed'
	    ia_found = nint(rlists(ind_ia_found,i,ia))	!found in this domain
	    if( ia_found == ia_needed ) cycle
	    if( ia_found < 1 .or. ia_found > n_threads ) stop 'error stop: 6'
	    if( ia_needed < 1 .or. ia_needed > n_threads ) stop 'error stop: 7'
	    if( ia_found == my_ia ) then
	      n = n_tvd_send(ia) + 1
	      if( n > n_tvd_s ) stop 'error stop: n > n_tvd_s'
	      tvd_send(n,ia) = i
	      n_tvd_send(ia) = n
	      n = n_tvd_index(ia) + 1
	      if( n > n_tvd ) stop 'error stop: n > n_tvd'
	      tvd_index(n,ia) = i
	      n_tvd_index(ia) = n
	      n_tvd_send_to(ia_needed) =  n_tvd_send_to(ia_needed) + 1
	    else if( ia_needed == my_ia ) then
	      n = n_tvd_receive(ia) + 1
	      if( n > n_tvd_r ) stop 'error stop: n > n_tvd_r'
	      tvd_receive(n,ia) = i
	      n_tvd_receive(ia) = n
	      n = n_tvd_index(ia) + 1
	      if( n > n_tvd ) stop 'error stop: n > n_tvd'
	      tvd_index(n,ia) = i
	      n_tvd_index(ia) = n
	      n_tvd_receive_from(ia_found) =  n_tvd_receive_from(ia_found) + 1
	    end if
	    if( ia_found == my_ia .or. ia_needed == my_ia ) then
	      n_my_array_list = n_my_array_list + 1
	      if( ia_found == my_ia ) then
	        n_my_send_list = n_my_send_list + 1
	      else
	        n_my_receive_list = n_my_receive_list + 1
	      end if
	    end if
	    matrix(ia_found,ia_needed) = matrix(ia_found,ia_needed) + 1
	    if( ia_needed /= ia_found ) then
	      n_tvd_all = n_tvd_all + 1
	      tvd_all(:,n_tvd_all) = rlists(:,i,ia)
	    end if
	  end do
	end do

	if( n_my_receive_list + n_my_send_list /= n_my_array_list ) then
	  stop 'error stop: n_receive + n_send /= n_all'
	end if
	allocate(my_receive_list(n_my_receive_list))
	allocate(my_send_list(n_my_send_list))
	allocate(my_array_receive_list(nlist,n_my_receive_list))
	allocate(my_array_send_list(nlist,n_my_send_list))
	allocate(my_array_list(nlist,n_my_array_list))
	my_receive_list = 0
	my_send_list = 0
	my_array_list = 0.

	nr = 0
	ns = 0
	na = 0
	do ia=1,n_threads
	  do i=1,maxlist
	    ie = nint(rlists(ind_ie_remote,i,ia))
	    if( ie == 0 ) exit			!end of list
	    ia_needed = nint(rlists(ind_ia_needed,i,ia))!needed in this domain
	    if( ia /= ia_needed ) stop 'error stop: ia/=ia_needed'
	    ia_found = nint(rlists(ind_ia_found,i,ia))	!found in this domain
	    if( ia_found == ia_needed ) cycle
	    if( ia_needed == my_ia ) then
	      na = na + 1
	      my_array_list(:,na) = rlists(:,i,ia)
	      nr = nr + 1
	      my_receive_list(nr) = nr
	      my_array_receive_list(:,nr) = rlists(:,i,ia)
	    else if( ia_found == my_ia ) then
	      na = na + 1
	      my_array_list(:,na) = rlists(:,i,ia)
	      ns = ns + 1
	      my_send_list(ns) = ns
	      my_array_send_list(:,ns) = rlists(:,i,ia)
	    end if
	  end do
	end do

	if( na /= n_my_array_list ) stop 'error stop: na/=n_my_array_list'
	if( nr /= n_my_receive_list ) stop 'error stop: nr/=n_my_receive_list'
	if( ns /= n_my_send_list ) stop 'error stop: nr/=n_my_send_list'

	allocate(my_xi_list(3,n_my_send_list))
	my_xi_list = 0.

	do i=1,n_my_send_list
	  na = my_send_list(i)
	  xd = my_array_list(ind_x,na)
	  yd = my_array_list(ind_y,na)
	  ie = my_array_list(ind_ie_remote,na)
	  call xy2xi(ie,xd,yd,xi)
	  my_xi_list(:,i) = xi(:)
	end do

	write(iudb,*) '====================================='
	write(iudb,*) 'my list:'
	write(iudb,*) '====================================='
	write(iudb,*) 'my array list of area: ',my_ia,n_my_array_list
	call write_array_header(iudb)
	do i=1,n_my_array_list
	  call write_array_debug(iudb,i,my_array_list(:,i))
	end do
	write(iudb,*) 'my receive list of area: ',my_ia,n_my_receive_list
	call write_array_header(iudb)
	do i=1,n_my_receive_list
	  na = my_receive_list(i)
	  call write_array_debug(iudb,i,my_array_list(:,na))
	end do
	write(iudb,*) 'my send list of area: ',my_ia,n_my_send_list
	call write_array_header(iudb)
	do i=1,n_my_send_list
	  na = my_send_list(i)
	  call write_array_debug(iudb,i,my_array_list(:,na))
	end do

	write(iudb,*) '====================================='
	write(iudb,*) 'send/receive index start ',my_ia
	write(iudb,*) n_tvd_r,n_tvd_s,n_tvd

	n_tvd_receive_total = sum( n_tvd_receive )
	write(iudb,*) 'receive max: ',n_tvd_r
	write(iudb,*) 'receive dim: ',n_tvd_receive
	write(iudb,*) 'receive tot: ',n_tvd_receive_total
	do ia=1,n_threads
	  n = n_tvd_receive(ia)
	  if( n == 0 ) cycle
	  write(iudb,*) 'receiving from... n = ',n
	  call write_array_header(iudb)
	  do i=1,n
	    ind = tvd_receive(i,ia)
	    call write_array_debug(iudb,i,rlists(:,ind,ia))
	  end do
	end do
	n_tvd_send_total = sum( n_tvd_send )
	write(iudb,*) 'send max: ',n_tvd_s
	write(iudb,*) 'send dim: ',n_tvd_send
	write(iudb,*) 'send tot: ',n_tvd_send_total
	do ia=1,n_threads
	  n = n_tvd_send(ia)
	  if( n == 0 ) cycle
	  write(iudb,*) 'sending to... n = ',n
	  call write_array_header(iudb)
	  do i=1,n
	    ind = tvd_send(i,ia)
	    call write_array_debug(iudb,i,rlists(:,ind,ia))
	  end do
	end do
	do ia=1,n_threads
	  n = n_tvd_index(ia)
	  if( ia == my_ia ) then
	    if(any(tvd_receive(1:n,ia)/=tvd_index(1:n,ia))) then
		write(6,*) 'receive,index: ',n_tvd_index(ia),n_tvd_receive(ia)
		stop 'error stop: receive'
	    end if
	  else
	    if(any(tvd_send(1:n,ia)/=tvd_index(1:n,ia))) then
		write(6,*)'send,index: ', n_tvd_index(ia),n_tvd_send(ia)
		stop 'error stop: send'
	    end if
	  end if
	end do
	write(iudb,*) '====================================='
	write(iudb,*) 'matrix exchange in domain ',my_ia
	do iaf=1,n_threads
	  do ian=1,n_threads
	    n = matrix(iaf,ian)
	    if( n == 0 ) cycle
	    write(iudb,'(i5,a,i5,a,i5)') iaf,' -> ',ian,' : ',n
	  end do
	end do
	write(iudb,*) '====================================='
	write(iudb,*) 'matrix all'
	write(iudb,'(a,10i5)') ' from/to:',(i,i=1,n_threads)
	do iaf=1,n_threads
	  write(iudb,'(i5,4x,10i5)') iaf,(matrix(iaf,ian),ian=1,n_threads)
	end do
	write(iudb,*) '====================================='
	do ia=1,n_threads
	  write(iudb,9) my_ia,' receives from ',ia,' : ',n_tvd_receive_from(ia)
	end do
	do ia=1,n_threads
	  write(iudb,9) my_ia,' sends to      ',ia,' : ',n_tvd_send_to(ia)
	end do
	write(iudb,*) '====================================='

	write(6,*) 'all new information collected ',my_id

	call shympi_barrier

!----------------------------------------------------------
! write out debug information
!----------------------------------------------------------

	!call write_tvd_debug(nel)

!----------------------------------------------------------
! all info exchanged
!----------------------------------------------------------

	!call test_tvd_convert

	call flush(iu)
	call flush(iudb)
	call shympi_barrier
	call flush(6)

	!stop

	return
    9	format(i5,a,i5,a,i10)
	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine allocate_tvd_arrays(nvals)

	use shympi
	use shympi_tvd

	implicit none

	integer nvals
	integer lmax

	if( allocated(values_tvd_local) ) return

	lmax = nlv_global

	allocate( values_tvd_local(lmax,n_tvd_receive_total) )
	allocate( values_tvd_remote(lmax,n_tvd_send_total) )

	end

!******************************************************************

	subroutine tvd_prepare_remote(values)

	use basin
	use levels
	use shympi
	use shympi_tvd

	implicit none

	real values(nlvdi,nkn)

	integer i,na,ie,ii,k,l,lmax,ia
	integer my_ia
	double precision c,cacum
	double precision xi(3)

	my_ia = my_id + 1

	do i=1,n_my_send_list
	  na = my_send_list(i)
	  if( i /= na ) stop 'error stop tvd_prepare_remote: internal (0)'
	  ia = nint(my_array_send_list(ind_ia_found,na))
	  if( ia /= my_ia ) stop 'error stop tvd_prepare_remote: internal (1)'
	  ie = nint(my_array_send_list(ind_ie_remote,na))
	  if( ie > nel ) stop 'error stop tvd_prepare_remote: internal (2)'
	  !write(6,*) 'ggguuu: ',i,na,ia,my_ia
	  lmax = ilhv(ie)
	  xi(:) = my_xi_list(:,na)
	  do l=1,lmax
	    cacum = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      c = values(l,k)
	      cacum = cacum + c * xi(ii)
	    end do
	    write(6,*) 'ggguuu: ',l,na,my_ia,cacum
	    values_tvd_remote(l,na) = cacum
	  end do
	end do

	end

!******************************************************************

	subroutine tvd_exchange

	use shympi
	use shympi_tvd

	implicit none


	end

!******************************************************************

	subroutine tvd_mpi_prepare(values)

	use basin
	use levels
	use shympi
	use shympi_tvd

	implicit none

	real values(nlvdi,nkn)

	integer, save :: nvals = 1

	if( .not. bmpi ) return

	call allocate_tvd_arrays(nvals)

	call tvd_prepare_remote(values)

	call tvd_exchange

	end

!******************************************************************

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine write_tvd_debug(nel)

	use shympi
	use shympi_tvd
	use mod_tvd

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
	      iedebug(n,ie) = ieetvdup(j,ii,ie)
	    end do
	  end do
	  iadebug(:,ie) = reshape(ieetvdup(:,:,ie),(/9/))
	  if( any( iedebug(:,ie) /= iadebug(:,ie) ) ) then
	    stop 'error stop write_tvd_debug: reshape not working'
	  end if
	  iadebug(:,ie) = reshape(iatvdup(:,:,ie),(/9/))
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

	write(6,*) 'finished tvd debug...',my_id
	call shympi_barrier
	flush(6)

	end

!******************************************************************

	subroutine tvd_assert(text,bassert)

	implicit none

	character*(*) text
	logical bassert

	if( bassert ) return

	write(6,*) 'assertion violated: ',trim(text)
	stop 'error stop tvd_assert'

	end

!******************************************************************

	subroutine write_array_debug(iudb,i,array)

	use shympi_tvd

	implicit none

	integer iudb,i
	real array(nlist)

	if( nlist_type == 1 ) then
	  write(iudb,'(10i6)') i,nint(array(3:10))
	else if( nlist_type == 2 ) then
	  write(iudb,'(10i6)') i,nint(array(3:nlist))
	end if

	end

!******************************************************************

	subroutine write_array_header(iudb)

	use shympi_tvd

	implicit none

	integer iudb

	if( nlist_type == 1 ) then
	  write(iudb,*) '    i   ian     j    ii    ie   iee   iaf   ief  ieef'
	else if( nlist_type == 2 ) then
	  write(iudb,*) '    i   ian     j    ii    ie   iee   iaf   ief  ieef'
	end if

	end

!******************************************************************

	subroutine test_tvd_convert

	use shympi_tvd

	implicit none

	integer k,j,ii,m

	write(6,*) 'testing tvd convert routines...'

	do k=1,6
	  call k2jii(k,j,ii)
	  call jii2k(j,ii,m)
	  if( k /= m ) stop 'error stop text_tvd_convert'
	  write(6,*) j,ii,k,m
	end do

	end

!******************************************************************

	!program main_test_tvd
	!call test_tvd_convert
	!end

!******************************************************************

