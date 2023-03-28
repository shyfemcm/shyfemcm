
!--------------------------------------------------------------------------
!
!    Copyright (C) 2007-2008,2010-2011,2013-2015,2013-2015  Georg Umgiesser
!    Copyright (C) 2018-2019  Georg Umgiesser
!    Copyright (C) 2008,2011,2013,2016  Debora Bellafiore
!    Copyright (C) 2008  Christian Ferrarin
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

! tripple point routines
!
! revision log :
!
! 18.03.2023	ggu	adjusted horizontal diffusion for mpi
! 27.03.2023	ggu	tripple point routines copied here
! 28.03.2023	ggu	insert lmax into list
!
!******************************************************************

!==================================================================
        module shympi_tripple
!==================================================================

	integer, parameter :: nexch = 7
	integer, save :: nmax_tripple = 0
	integer, save :: itrtot = -1
	integer, save, allocatable :: ielist(:,:)
	real, save, allocatable :: buffer_tripple(:,:)

!==================================================================
        end module shympi_tripple
!==================================================================

	subroutine tripple_points_handle

	use shympi
	use shympi_tripple

	implicit none

	!call exchange_areas
	!call test_exchange_3d

	!return

	if( itrtot == -1 ) call tripple_points_init

	write(6,*) 'finished tripple_points_init: ',my_id
	call shympi_syncronize

	!if( itrtot > 0 ) stop

	call tripple_points_exchange

	!stop

	end

!******************************************************************

	subroutine tripple_points_init

! constructs list of tripple points and info on how to exchange

	use basin
	use shympi
	use shympi_tripple

	implicit none

	integer ie,itr,idn,ide,iee
	integer ii,i1,i2,ip
	integer k,k1,k2,kext1,kext2
	integer i,j,iei,iext,iint,ia,itmax
	integer, allocatable :: ies(:)

	integer, allocatable :: itrs(:)
	integer, allocatable :: iexch(:,:)
	integer, allocatable :: iexchs(:,:)
	integer, allocatable :: ips(:,:)
	integer, allocatable :: id_elem_g(:,:)

	integer ieext,ieint,ipext

	if( itrtot >= 0 ) return	!already set up

	!--------------------------------------------------
	! get general information on tripple points
	!--------------------------------------------------

	itr = 0
	do ie=1,nel_unique
	  idn = id_elem(0,ie)
	  if( idn == 3 ) then
	    write(6,1000) 'tripple point found: ',my_id,ie,id_elem(:,ie)
	    itr = itr + 1
	  end if
	end do

	call shympi_syncronize

	itrtot = shympi_sum(itr)
	write(6,*) 'total numbers of tripple points: ',my_id,itr,itrtot

	nmax_tripple = 2 * (nlv_global+1)
	allocate(buffer_tripple(nmax_tripple,itrtot))
	allocate(ielist(nexch,itrtot))
	ielist = 0

	if( itrtot == 0 ) return

	call shympi_syncronize

	!--------------------------------------------------
	! allocate arrays
	!--------------------------------------------------

	itmax = shympi_max(itr)

	allocate(ies(itmax))
	allocate(iexch(nexch,itmax))
	allocate(iexchs(nexch,n_threads))
	allocate(itrs(n_threads))
	allocate(ips(itmax,n_threads))
	iexch = 0
	ies = 0
	itr = 0
	itrs = 0

	!--------------------------------------------------
	! insert local tripple elements in ies
	!--------------------------------------------------

	do ie=1,nel_unique
	  idn = id_elem(0,ie)
	  if( idn == 3 ) then
	    itr = itr + 1
	    ies(itr) = ie
	  end if
	end do

	call shympi_gather(itr,itrs)

	!--------------------------------------------------
	! find neighbor elements
	!--------------------------------------------------

	do i=1,itr
	  ie = ies(i)
	  iexch(1,i) = ie
	  iexch(2,i) = ipev(ie)
	  iexch(3,i) = my_id
	  call find_elem_neib(ie,iext)
	  iexch(5,i) = iext
	end do

	!--------------------------------------------------
	! find id of neighbor element and remember in iexch
	!--------------------------------------------------

	call find_elem_id(nexch,itr,iexch)

	!--------------------------------------------------
	! construct pointer where to insert item in list
	!--------------------------------------------------

	ip = 0
	ips = 0
	do ia=1,n_threads
	  itr = itrs(ia)
	  do i=1,itr
	    ip = ip + 1
	    ips(i,ia) = ip	!pointer where to insert info
	  end do
	end do

	if( ip /= itrtot ) stop 'error stop tripple_points: internal (3)'

	!--------------------------------------------------
	! insert items in list (list will be identical in all domains)
	!--------------------------------------------------

	do i=1,itmax
	  call shympi_gather(iexch(:,i),iexchs)
	  do ia=1,n_threads
	    ip = ips(i,ia)
	    if( ip > 0 ) then
	      ielist(:,ip) = iexchs(:,ia)
	    end if
	  end do
	end do

	call shympi_syncronize

	!--------------------------------------------------
	! debug output
	!--------------------------------------------------

	write(6,*) 'list of tripple points:'
	do i=1,itrtot
	  write(6,'(10i8)') my_id,ielist(:,i)
	  call shympi_gather(ielist(:,i),iexchs)
	  do j=1,nexch
	    if( any( iexchs(j,:) /= iexchs(j,1) ) ) then
	      write(6,*) 'values are different: ',iexchs(j,:)
	      stop 'error stop tripple_points: internal (11)'
	    end if
	  end do
	  iext = ielist(5,i)
	  iint = ieint(iext)
	  if( iint > 0 .and. iint <= nel_unique ) then
	    if( iint /= ielist(4,i) ) then
	      write(6,'(a,10i8)') 'error: ',my_id,iext,iint,nel_unique
	      stop 'error stop tripple_points: internal (8)'
	    end if
	  end if
	end do

	call shympi_syncronize

	!--------------------------------------------------
	! end of routine
	!--------------------------------------------------

	return
 1000	format(a,6i8)
	end

!******************************************************************

	subroutine tripple_points_exchange

! exchanges information on tripple points

	use basin
	use levels
	use shympi
	use shympi_tripple

	implicit none

	integer i,itr
	integer iint,iext,id
	integer id_from,id_to
	integer nmax,lmax

	integer ieint

	do i=1,itrtot
	  write(6,'(a,10i8)') 'tr_exchange: ',my_id,ielist(:,i)
	  iint = ielist(4,i)
	  iext = ielist(5,i)
	  id_to = ielist(3,i)
	  id_from = ielist(6,i)
	  lmax = ielist(7,i)
	  itr = i
	  call send_elem_info(id_from,id_to,iint,lmax,itr)
	end do

	write(6,'(a,10i8)') 'finished tr_exchange: ',my_id

	if( .not. bmpi ) then
	  iext = 39057
	  iint = ieint(iext)
	  lmax = ilhv(iint)
	  id = 0
	  call send_elem_info(id,id,iint,lmax,0)
	end if

	write(6,'(a,10i8)') 'end of exchange: ',my_id
	call shympi_syncronize
	flush(6)
	
	end

!******************************************************************

	subroutine send_elem_info(id_from,id_to,iint,lmax,itr)

	use levels
	use evgeom
	use mod_hydro
	use shympi
	use shympi_tripple

	implicit none

	integer id_from,id_to,iint,lmax,itr

	integer iu
	integer l,n
	integer nmax
	real buffer_in(nmax_tripple)
	real buffer_out(nmax_tripple)

	integer ieext

	nmax = nmax_tripple
	n = 2*lmax + 2

	write(6,*) 'send: ',my_id,lmax,id_from,id_to,itr,n
	buffer_in(1) = lmax
	buffer_in(2) = 12. * ev(10,iint)
	buffer_in(3:lmax+2) = utlnv(1:lmax,iint)
	buffer_in(lmax+3:2*lmax+2) = utlnv(1:lmax,iint)

	if( n > nmax ) then
	  write(6,*) 'n,nmax: ',n,nmax
	  stop 'error stop send_elem_info: n>nmax'
	end if

	if( bmpi ) then
	  iu = 300 + my_id
	  write(6,*) 'send before: ',my_id,id_from,id_to
	  call shympi_receive(id_from,id_to,n,buffer_in,buffer_out)
	  write(6,*) 'send after: ',my_id
	  !buffer_tripple(1:n,itr) = buffer_out(1:n)
	  buffer_tripple(:,itr) = buffer_out(:)
	else
	  iu = 400
	  buffer_out = buffer_in
	end if

	write(iu,*) lmax,n
	write(iu,*) buffer_out(1:n)

	end

!******************************************************************
!******************************************************************
!******************************************************************
! utility routines
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine find_elem_id(nexch,itr,iexch)

! computes id and internal ie from list in iexch and inserts it there

	use levels
	use shympi

	implicit none

	integer nexch,itr
	integer iexch(nexch,itr)

	integer i,iext,ide,ie,iee,iint,lmax
	integer, allocatable :: id_elem_g(:,:)
	integer, allocatable :: ilhv_g(:)

	allocate(id_elem_g(0:3,nel_global))
	allocate(ilhv_g(nel_global))
	call shympi_l2g_array(4,id_elem,id_elem_g)
	call shympi_l2g_array(ilhv,ilhv_g)

	do i=1,itr
	  iext = iexch(5,i)
	  ide = -1
	  do ie=1,nel_global
	    if( ip_ext_elem(ie) == iext ) then
	      ide = id_elem_g(1,ie)
	      iee = ie
	      lmax = ilhv_g(ie)
	    end if
	  end do
	  if( ide == -1 ) stop 'error stop tripple_points: internal (4)'
	  do ie=1,nel_unique
	    if( iee == ip_int_elems(ie,ide+1) ) then
	      iint = ie
	    end if
	  end do
	  write(6,*) 'ide found: ',iext,iint,ide
	  iexch(4,i) = iint
	  iexch(5,i) = iext
	  iexch(6,i) = ide
	  iexch(7,i) = lmax
	end do

	end

!******************************************************************

	subroutine find_elem_neib(ie,iext)

! find neighbor element opposite to node with my_id and returns external number

	use basin
	use shympi

	implicit none

	integer ie,iext

	integer ii,i1,i2,iei
	integer k,k1,k2,kg1,kg2

	if( ie == 0 ) then
	  stop 'error stop find_elem_neib: ie == 0'
	end if

	k1 = 0
	k2 = 0
	do ii=1,3
	  k = nen3v(ii,ie)
	  if( id_node(k) == my_id ) then
	    i1 = mod(ii,3) + 1
	    i2 = mod(i1,3) + 1
	    k1 = nen3v(i1,ie)
	    k2 = nen3v(i2,ie)
	  end if
	end do

	if( k1 == 0 .or. k2 == 0 ) then
	  stop 'error stop find_elem_neib: k1 or k2 == 0'
	end if

	kg1 = ip_int_node(k1)
	kg2 = ip_int_node(k2)

	iei = 0

	do ie=1,nel_global
	  do ii=1,3
	    i1 = mod(ii,3) + 1
	    if( nen3v_global(i1,ie) == kg1 
     +			.and. nen3v_global(ii,ie) == kg2 ) then
	      if( iei /= 0 ) stop 'error stop find_elem_neib: int (1)'
	      iei = ie
	    end if
	  end do
	end do

	if( iei == 0 ) then
	  write(6,*) 'cannot find element to nodes ',kg1,kg2
	  stop 'error stop find_elem_neib: int(2)'
	end if

	iext = ip_ext_elem(iei)

	end

!******************************************************************
!******************************************************************
!******************************************************************
! old routines (not used)
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine test_exchange_3d

	use basin
	use levels
	use shympi

	implicit none

	integer ie,iee,l,lmax,idiff,ierr
	real val,diff
	real vals(nlvdi,nel)
	real valsaux(nlvdi,nel)
	real vals2d(nel)
	real vals2daux(nel)
	integer ivals2d(nel)
	integer ivals2daux(nel)

	integer ieext

	!write(6,*) 'running test_exchange_3d'

	ivals2d = 0
	ivals2daux = 0

	do ie=1,nel
	  iee = ieext(ie)
	  ivals2daux(ie) = iee
	  if( ie > nel_unique ) cycle
	  ivals2d(ie) = iee
	end do

	call shympi_exchange_2d_elem(ivals2d)

	write(550+my_id,*) nel,nel_unique,nel_local,nel_local-nel_unique
	write(540+my_id,*) nel,nel_unique,nel_local,nel_local-nel_unique

	ierr = 0
	do ie=1,nel
	  idiff = abs(ivals2d(ie)-ivals2daux(ie))
	  if( idiff > 0 ) then
	    ierr = ierr + 1
	    write(550+my_id,*) ierr,ie,ivals2d(ie),ivals2daux(ie),idiff
	  end if
	  write(540+my_id,*) ie,id_elem(:,ie)
	end do

	if( ierr > 0 ) then
	  write(6,*) 'running test_exchange_3d'
	  write(6,*) 'my_id = ',my_id,'  ierr = ',ierr
	  write(6,*) 'output in files 550+ and 540+'
	  write(6,*) '*** errors in test_exchange_3d ***'
	  stop 'error stop'
	end if

!-------------------------------------------------------------

	vals = 0.

	do ie=1,nel_unique
	  iee = ieext(ie)
	  lmax = ilhv(ie)
	  do l=1,lmax
	    val = 100*iee + l
	    vals(l,ie) = val
	  end do
	  vals2d(ie) = iee
	end do

	valsaux = vals
	vals2daux = vals2d
	call shympi_exchange_3d_elem(vals)
	call shympi_exchange_2d_elem(vals2d)

	write(670+my_id,*) nel,nel_unique,nel_local
	write(680+my_id,*) nel,nel_unique,nel_local
	write(690+my_id,*) nel,nel_unique,nel_local

	ierr = 0
	do ie=1,nel_unique
	  lmax = ilhv(ie)
	  do l=1,lmax
	    diff = abs(vals(l,ie)-valsaux(l,ie))
	    if( diff > 0 ) then
	      ierr = ierr + 1
	      write(670+my_id,*) ie,l,vals(l,ie),valsaux(l,ie),diff
	      write(680+my_id,*) ie,l,vals(l,ie),diff
	      write(690+my_id,*) ie,l,valsaux(l,ie),diff
	    end if
	  end do
	  diff = abs(vals2d(ie)-vals2daux(ie))
	  if( diff > 0 ) then
	    ierr = ierr + 1
	    write(660+my_id,*) ie,vals2d(ie),vals2daux(ie),diff
	  end if
	end do

	call shympi_check_3d_elem_r(vals,'tesssssssst')

	if( ierr > 0 ) stop 'error stop'

	end

!******************************************************************

	subroutine exchange_areas

	use mod_geom
	!use mod_internal
	use mod_hydro
	use evgeom
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

	implicit none

	integer ie,ii,iei,ia,ia0,iee
	real, parameter :: flag = -777.
	real, allocatable :: areat(:,:,:)
	real, allocatable :: areas(:,:)

	integer ieext

	allocate(areat(3,nel,n_threads))
	allocate(areas(3,nel))

	areas = flag

	do ie=1,nel
	  do ii=1,3
            iei = ieltv(ii,ie)
	    if( iei > 0 ) areas(ii,ie) = 12. * ev(10,iei)
	  end do
	end do

	call shympi_barrier
	call shympi_gather(3,areas,areat)

	ia0 = my_id + 1

	do ie=1,nel_unique
	  iee = ieext(ie)
	  do ii=1,3
	    if( areat(ii,ie,ia0) == flag ) then
	      do ia=1,n_threads
		if( ia == ia0 ) cycle
		if( areat(ii,ie,ia) /= flag ) then
		  write(6,*) ie,ii,ia0,ia,areat(ii,ie,ia)
		end if
	      end do
	    end if
	  end do
	end do

	write(6,*) 'stopping in exchange_areas'
	stop

	end

!******************************************************************
