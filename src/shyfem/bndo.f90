
!--------------------------------------------------------------------------
!
!    Copyright (C) 2001,2003-2004,2007-2012,2014-2015  Georg Umgiesser
!    Copyright (C) 2017-2019  Georg Umgiesser
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

! routines for open boundary conditions
!
! contents :
!
! subroutine bndo_init
!       sets up bndo data structure
! subroutine bndinsert(ib,area,kn)
!       inserts node kn with weight area into list (internal routine)
!
! subroutine bndo_info(iunit)
!       writes info on open boundary nodes to terminal
!
! function is_zeta_bound(k)
!	checks if node k is a zeta boundary
!
! subroutine bndo_setbc(what,nlvddi,cv,rbc,uprv,vprv)
!	sets open boundary condition for level boundaries
! subroutine bndo_impbc(what,nlvddi,cv,rbc)
!       imposes boundary conditions on open boundary
! subroutine bndo_adjbc(nlvddi,cv,uprv,vprv)
!       adjusts boundary conditions on open boundary (values on bnd already set)
!
! subroutine bndo_radiat(it,rzv)
!	imposes radiation condition for levels
!
! notes :
!
!       integer nbndo                   !total number of OB nodes
!
!	integer iopbnd(nkn)		!if >0 pointer into array irv
!
!	integer ibcnod(nrb)		!number of boundary
!	integer kbcnod(nrb)		!number of boundary node
!	integer itynod(nrb)          	!type of boundary
!	logical bexnod(nrb)          	!boundary node is external
!
!	real xynorm(2,nrb)		!normal direction for OB node
!
!	integer nopnod(nrb)		!number of internal nodes close to OB
!	integer nopnodes(ngr,nrb)	!nodes close to OB
!	real wopnodes(ngr,nrb)		!weights of nodes close to OB
!
!	iopbnd(k) = 0            no open BC
!	iopbnd(k) > 0            BC (any ibtyp)
!
! the initialization routines should be called only after the
! ... arrays kantv and ieltv have been setup
!
!--------------------------------------------------------------------------
!
! revision log :
!
! 15.01.2001	ggu	written from scratch
! 03.12.2001	ggu	LEVMX - look out for missing level of near node
! 05.12.2001	ggu	NTBC - BUG -> has not been set before
! 27.03.2003	ggu	in bndo_adjbc use ambient value (bamb)
! 13.03.2004	ggu	in bndo_adjbc only for level BC (LEVELBC)
! 05.10.2004	ggu	new routine bndo_radiat, ibtyp=31
! 31.05.2007	ggu	reset BC for flux to old type (DEBHELP)
! 23.08.2007	ggu	use iopbnd as indicator for ext/int boundary nodes
! 08.04.2008	ggu	file cleaned (new bndo_setbc, bndo_impbc)
! 17.04.2008	ggu	calls to infobnd deleted (subst by get_bnd_ipar)
! 03.09.2008	ggu	new routine bndo_info_file()
! 06.11.2008	ggu	better error handling
! 12.11.2009	ggu	new array itynod and is_zeta_bound()
! 23.03.2010	ggu	changed v6.1.1
! 17.02.2011	ggu	changed VERS_6_1_18
! 25.03.2011	ggu	bug fix in bndo_impbc() -> ibcold not initialized
! 14.04.2011	ggu	changed VERS_6_1_22
! 30.03.2012	ggu	changed VERS_6_1_51
! 18.06.2014	ggu	changed VERS_6_1_77
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 05.06.2015	ggu	changed VERS_7_1_12
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 29.09.2015	ggu	changed VERS_7_2_5
! 18.12.2015	ggu	changed VERS_7_3_17
! 05.12.2017	ggu	changed VERS_7_5_39
! 03.04.2018	ggu	changed VERS_7_5_43
! 19.04.2018	ggu	changed VERS_7_5_45
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 26.04.2022	ggu	implementing OB in more than one domain
! 03.05.2022	ggu	lots of debug code integrated
! 22.03.2023	ggu	relax error condition on one node open boundary
! 01.04.2023	ggu	new array bexnod, handle fix boundary condition (bfix)
! 02.04.2023	ggu	merged files subbndo.f and mod_bndo.f
! 21.04.2023	ggu	call bexnod() only for existing boundaries
! 29.10.2024	ggu	check if boundary node is unique
! 12.06.2025	ggu	extra info on error
!
!--------------------------------------------------------------------------

!==================================================================
	module mod_bndo
!==================================================================

	implicit none

	integer, private, save  :: nrb_bndo = 0
	integer, private, save  :: ngr_bndo = 0

	integer, save :: nbndo = 0	!total number of OB nodes
	integer, save :: ndebug = 0	!unit number for debug messages

	integer, save :: kbcdim = 0	!to be deleted later
	integer, save :: kopdim = 0	!to be deleted later

	integer, allocatable, save :: nopnod(:)	!number of nodes close to OB
	integer, allocatable, save :: ibcnod(:)	!number of boundary
	integer, allocatable, save :: kbcnod(:)	!number of boundary node
	integer, allocatable, save :: itynod(:)	!type of boundary
	logical, allocatable, save :: bexnod(:)	!is boundary external (1 or 2)
	integer, allocatable, save :: nopnodes(:,:)	!nodes close to OB
	real, allocatable, save :: xynorm(:,:)		!normal direction
	real, allocatable, save :: wopnodes(:,:)	!weights

	real, parameter :: flag_bndo = -999.

!==================================================================
	contains
!==================================================================

	subroutine mod_bndo_init(ngr,nrb)

	integer ngr,nrb

	integer nlk,naux

        if( ngr == ngr_bndo .and. nrb == nrb_bndo ) return

        if( ngr_bndo > 0 ) then
          deallocate(nopnod)
          deallocate(ibcnod)
          deallocate(kbcnod)
          deallocate(itynod)
          deallocate(bexnod)
          deallocate(xynorm)
          deallocate(nopnodes)
          deallocate(wopnodes)
        end if

        ngr_bndo = ngr
        nrb_bndo = nrb

	kbcdim = nrb		!to be deleted later
	kopdim = ngr		!to be deleted later

	if( ngr == 0 ) return

	naux = max(1,nrb)

        allocate(nopnod(naux))
        allocate(ibcnod(naux))
        allocate(kbcnod(naux))
        allocate(itynod(naux))
        allocate(bexnod(naux))
        allocate(xynorm(2,naux))
        allocate(nopnodes(ngr,naux))
        allocate(wopnodes(ngr,naux))

	end subroutine mod_bndo_init

	subroutine mod_bndo_info

	write(6,*) 'mod_bndo_info ================='
	write(6,*) 'ngr_bndo: ',ngr_bndo
	write(6,*) 'nrb_bndo: ',nrb_bndo
	write(6,*) 'nbndo: ',nbndo
	write(6,*) 'mod_bndo_info end ================='

	end subroutine mod_bndo_info

!==================================================================
	end module mod_bndo
!==================================================================


!
!***********************************************************************

	subroutine bndo_init

! sets up bndo data structure

	use mod_bndo
	use mod_bound_geom
	use mod_geom
	use basin
	use shympi

	implicit none

	logical bexternal
	logical berror,bdebug
	integer k,nodes,itype
	integer i,ibc
	integer inext,ilast,knext,klast
	integer ie,n,ngood,ie_mpi,id
	integer ii,iii,ib,in,kn,nb,j,ip
	integer nbc
	integer iunit,kint,kext
	real area
	real dx,dy

	integer nkbnds,itybnd,kbnds,ipext,ipint,nbnds
	real areaele

!----------------------------------------------------------
! set up array iopbnd
!----------------------------------------------------------

	iopbnd(:) = 0

	nbndo = 0
        ndebug = 0              !unit number for debug (in module)
	bdebug = ( ndebug > 0 )

	nbc = nbnds()
	if( bdebug ) write(ndebug,*) 'boundaries: ',nbc,bmpi

	do ibc = 1,nbc
	  nodes = nkbnds(ibc)
	  itype = itybnd(ibc)
	  if( bdebug ) write(ndebug,*) 'boundary: ',ibc,nodes,itype

	  ngood = 0
	  bexternal = ( itype .ge. 1 .and. itype .le. 2                & 
      &      		    .or. itype .ge. 31 .and. itype .le. 39 )

	  do i=1,nodes
	    k = kbnds(ibc,i)
	    if( k <= 0 ) then
	      if( bmpi ) cycle			!node can be in other domain
	      write(6,*) ibc,i,k
	      stop 'error stop bndo_init: no such node'
	    end if
	    ngood = ngood + 1
	    nbndo = nbndo + 1
	    if( nbndo .gt. kbcdim ) goto 99
	    if( iopbnd(k) > 0 ) then
	      ip = iopbnd(k)
	      write(6,*) 'boundary node already inserted'
	      write(6,*) 'boundary: ',ibc
	      write(6,*) 'node already in boundary: ',ibcnod(ip)
	      write(6,*) 'kint,kext: ',k,ipext( k)
	      write(6,*) 'boundary node is in two boundaries... not possible'
	      stop 'error stop bndo_init: not unique boundary node'
	    end if
	    iopbnd(k) = nbndo
	    ibcnod(nbndo) = ibc
	    kbcnod(nbndo) = k
	    itynod(nbndo) = itype
	    bexnod(nbndo) = bexternal		!flag this node as external
	    if( bdebug ) write(ndebug,*) ibc,i,k,nbndo
            !if( iopbnd(k) .ne. i ) then
	  end do

	  if( ngood == 0 ) then
	    !boundary not in domain
	  else if( ngood == nodes ) then
	    !boundary fully in domain
	  else if( itype == 2 ) then	!boundary only partially in domain
	    write(6,*) 'ngood,nodes: ',ngood,nodes
	    write(6,*) 'boundary is only partially in domain'
	    write(6,*) 'cannot handle flux OB in different domains yet'
	    stop 'error stop bndo_init: boundary between domains'
	  else			!zeta boundary - should be able to handle
	    !
	  end if
	end do

!----------------------------------------------------------
! set up normal direction
!----------------------------------------------------------

	do i=1,nbndo
	  k = kbcnod(i)
	  ibc = ibcnod(i)

	  itype = itybnd(ibc)
	  bexternal = bexnod(i)

	  if( .not. bexternal ) cycle

	  knext = kantv(1,k)
	  klast = kantv(2,k)

	  inext = 0
	  ilast = 0
	  if( knext > 0 ) inext = iopbnd(knext)
	  if( klast > 0 ) ilast = iopbnd(klast)

!	  -------------------------------
!	  internal consistency check
!	  -------------------------------

	  if( iopbnd(k) .ne. i ) then
	    stop 'error stop bndo_init: internal error (0)'
	  end if
	  if( inext .gt. 0 ) then
	   if( kbcnod(inext) .ne. knext ) then
	    stop 'error stop bndo_init: internal error (1)'
	   end if
	  end if
	  if( ilast .gt. 0 ) then
	   if( kbcnod(ilast) .ne. klast ) then
	    stop 'error stop bndo_init: internal error (2)'
	   end if
	  end if

!	  -------------------------------
!	  adjacent boundary nodes must be of same boundary
!	  -------------------------------

	  if( inext .gt. 0 ) then
	    if( ibcnod(inext) .ne. ibc ) goto 98
	  end if
	  if( ilast .gt. 0 ) then
	    if( ibcnod(ilast) .ne. ibc ) goto 98
	  end if

!	  -------------------------------
!	  get normal direction (might be wrong in case of domain border)
!	  -------------------------------

	  if( inext .gt. 0 .and. ilast .gt. 0 ) then	!inner node in OB
	    dx = xgv(knext) - xgv(klast)
	    dy = ygv(knext) - ygv(klast)
	  else if( inext .gt. 0 ) then			!first node in OB
	    dx = xgv(knext) - xgv(k)
	    dy = ygv(knext) - ygv(k)
	  else if( ilast .gt. 0 ) then			!last node in OB
	    dx = xgv(k) - xgv(klast)
	    dy = ygv(k) - ygv(klast)
	  else
	    id = id_node(k)
	    if( id == my_id ) then
	      write(6,*) 'One node open boundary not permitted'
	      write(6,*) 'node in boundary: ',i
	      write(6,*) 'node number (internal): ',k
	      write(6,*) 'node number (external): ',ipext(k)
	      write(6,*) 'domain: ',my_id
	      write(6,*) 'boundary: ',ibc
	      write(6,*) 'knext,klast: ',knext,klast
	      stop 'error stop bndo'
	    else
	      dx = 0.
	      dy = 0.
	    end if
	  end if

	  xynorm(1,i) = -dy			!x-component
	  xynorm(2,i) = dx			!y-component
	end do
  
!----------------------------------------------------------
! set up weights and node list
!----------------------------------------------------------

	do i=1,nbndo
	  nopnod(i) = 0
	end do

	do ie_mpi=1,nel

	  ie = ip_sort_elem(ie_mpi)
	  n = 3
	  area = areaele(ie)

	  do ii=1,n
	    k = nen3v(ii,ie)
	    ib = iopbnd(k)
	    if( ib < 1 ) cycle
	    bexternal = bexnod(ib)
	    if( bexternal ) then		!insert inner nodes
	      do iii=1,n
		kn = nen3v(iii,ie)
	        in = iopbnd(kn)
		if( in .le. 0 ) then		!only inner nodes
		  call bndinsert(ib,area,kn)
		end if
	      end do
	    end if
	  end do

	end do

!----------------------------------------------------------
! scale weights to unit
!----------------------------------------------------------

	iunit = 730 + my_id
	kext = 6651
	kext = 0
	kint = ipint(kext)

	berror = .false.

	do i=1,nbndo

	  bexternal = bexnod(i)
	  if( .not. bexternal ) cycle

	  nb = nopnod(i)
	  if( nb .le. 0 ) then
	    k = kbcnod(i)
	    ibc = ibcnod(i)
	    !write(6,*) i,k,ipext(k)
	    write(6,*) '*** No inner nodes for boundary node'
	    write(6,*) '      boundary = ',ibc,'   node = ',ipext(k)
	    berror = .true.
	  end if

	  area = 0.
	  do j=1,nb
	    area = area + wopnodes(j,i)
	  end do
	  do j=1,nb
	    if( area .gt. 0. ) then
	      wopnodes(j,i) = wopnodes(j,i) / area
	    end if
	  end do
	  
	  bdebug = ( kint == kbcnod(i) )
	  if( bdebug ) then
	    iunit = 730 + my_id
	    write(iunit,*) '--------- bndo_init ------------'
	    write(iunit,*) i,kint,kext,nb,id_node(k)
	    write(iunit,*) (wopnodes(j,i),j=1,nb)
	    write(iunit,*) '--------- end bndo_init ------------'
	  end if
	end do

	if( berror ) stop 'error stop bndo'

	write(6,*) 'finished setting up bndo_init, nbndo = ',nbndo

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	return
   98	continue
	write(6,*) 'different boundary : ',ibc
	stop 'error stop bndo'
   99	continue
	write(6,*) 'dimension error kbcdim : ',kbcdim
	stop 'error stop bndo'
	end

!***********************************************************************

	subroutine bndinsert(ib,area,kn)

! inserts node kn with weight area into list (internal routine)

	use mod_bndo

	implicit none

	integer ib		!nodal index
	real area		!weight
	integer kn		!node number to insert

	integer nb,j,k

	nb = nopnod(ib)

	do j=1,nb
	  k = nopnodes(j,ib)
	  if( k .eq. kn ) then		!already there -> add weight
	    wopnodes(j,ib) = wopnodes(j,ib) + area
	    return			!done
	  end if
	end do

	nb = nb + 1
	if( nb .gt. kopdim ) then
	  write(6,*) 'Too much inner neighbors: ',nb
	  write(6,*) 'Please raise kopdim in subbndo.h'
	  stop 'error stop bndo'
	end if

	nopnodes(nb,ib) = kn
	wopnodes(nb,ib) = area
	nopnod(ib) = nb

	end

!***********************************************************************

	subroutine bndo_info_file(file)

! writes bndo info to file

	use shympi

	implicit none

	character*(*) file

	if( file .eq. ' ' ) return

        if(shympi_is_master()) then
	  open(1,file=file,status='unknown',form='formatted')
	  call bndo_info(1)
	  close(1)
	end if

	end

!***********************************************************************

	subroutine bndo_info(iunit)

! writes info on open boundary nodes to terminal

	use mod_bndo
	use mod_bound_geom
	use levels

	implicit none

        integer iunit

	integer i,k,ibc,nb,j,iu
	integer itybnd,ipext

        iu = iunit
        if( iu .le. 0 ) iu = 6

	write(iu,*) '--------------------------------'
	write(iu,*) 'Information on open boundary nodes'
	write(iu,*) 'Total number of open boundary nodes: ',nbndo

	do i=1,nbndo

	  k = kbcnod(i)
	  ibc = ibcnod(i)
	  nb = nopnod(i)

	  if( iopbnd(k) .ne. i ) then
	    write(6,*) i,iopbnd(k)
	    stop 'error stop bndo_info: internal error (11a)'
	  end if

	  write(iu,*) '-------------------------------- bndo_info'
	  write(iu,*) i,k,ipext(k),ibc,itybnd(ibc),ilhkv(k)
	  write(iu,*) nb
	  write(iu,*) (ipext(nopnodes(j,i)),j=1,nb)
	  write(iu,*) (wopnodes(j,i),j=1,nb)
	  write(iu,*) (ilhkv(nopnodes(j,i)),j=1,nb)
	  write(iu,*) (xynorm(j,i),j=1,2)
	end do

	write(iu,*) '--------------------------------'

	end

!***********************************************************************

	function is_zeta_bound(k)

! checks if node k is a zeta boundary

	use mod_bndo
	use mod_bound_geom

	implicit none

	logical is_zeta_bound
	integer k

	integer ip

	is_zeta_bound = .false.
	ip = iopbnd(k)

	if( ip .gt. 0 ) then
	  if( itynod(ip) .eq. 1 ) then
	    is_zeta_bound = .true.
	  end if
	end if
	
	end

!***********************************************************************
!***********************************************************************
!***********************************************************************

        subroutine bndo_setbc(what,nlvddi,cv,rbc,uprv,vprv)

! sets open boundary condition for level boundaries
!
! simply calls bndo_impbc() and bndo_adjbc()

	use basin
	use shympi

        implicit none

        character*(*) what      !conz/temp/salt or else
        integer nlvddi
        real cv(nlvddi,nkn)
        real rbc(nlvddi,nkn)	!boundary condition (3D)
	real uprv(nlvddi,nkn)
	real vprv(nlvddi,nkn)

	logical bdebug
	integer k,iunit,kext
	integer ipint

	iunit = 730 + my_id
	kext = 6651
	kext = 0
	k = ipint(kext)
	bdebug = ( k > 0 )

!----------------------------------------------------------
! simply imposes whatever is in rbc
!----------------------------------------------------------

	if( bdebug ) then
	write(iunit,*) '------------ bndo_setbc 1 -------'
	write(iunit,*) cv(:,k)
	write(iunit,*) rbc(:,k)
	end if
	
        call bndo_impbc(what,nlvddi,cv,rbc)

!----------------------------------------------------------
! adjusts for ambient value, no gradient or outgoing flow
!----------------------------------------------------------

	if( bdebug ) then
	write(iunit,*) '------------ bndo_setbc 2 -------'
	write(iunit,*) cv(:,k)
	write(iunit,*) rbc(:,k)
	end if
	
	call bndo_adjbc(what,nlvddi,cv,uprv,vprv)

	if( bdebug ) then
	write(iunit,*) '------------ bndo_setbc 3 -------'
	write(iunit,*) cv(:,k)
	write(iunit,*) rbc(:,k)
	end if
	
!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	end

!***********************************************************************

        subroutine bndo_impbc(what,nlvddi,cv,rbc)

! imposes boundary conditions on open boundary

	use mod_bndo
	use mod_bound_geom
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        character*(*) what      !conz/temp/salt or else
        integer nlvddi
        real cv(nlvddi,nkn)
        real rbc(nlvddi,nkn)	!boundary condition (3D)

        logical blevel,bfix
        logical bdebug
        integer i,j,k,l
        integer ibc,ibcold
        integer nb,lmax,lmin
        integer ibtyp
        real value,rb,flag

        integer ifemopa

        bdebug = .true.
        bdebug = .false.

	ibcold = 0

	flag = flag_bndo

        do i=1,nbndo

          k = kbcnod(i)
          nb = nopnod(i)
          ibc = ibcnod(i)

	  if( ibc .ne. ibcold ) then
	    call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	    ibcold = ibc
	  end if

	  blevel = ibtyp .eq. 1
	  bfix = ibtyp .eq. 5
          lmax = ilhkv(k)
          lmin = jlhkv(k)

          if( iopbnd(k) .ne. i ) then
	    write(6,*) i,iopbnd(k)
	    stop 'error stop bndo_impbc: internal error (11b)'
	  end if

	  if( blevel .or. bfix ) then
            do l=lmin,lmax
	      rb = rbc(l,k)
              if( rb /= flag ) cv(l,k) = rb
            end do
	  end if

        end do

        end

!***********************************************************************

	subroutine bndo_adjbc(what,nlvddi,cv,uprv,vprv)

! adjusts boundary conditions on open boundary (values on bnd already set)
!
! adjusts for ambient value, no gradient or outgoing flow
!
! this is only done on level boundaries ( ibtyp == 1 )

	use mod_bndo
	use mod_bound_geom
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

	implicit none

        character*(*) what	!conz/temp/salt or else
	integer nlvddi
	real cv(nlvddi,nkn)
	real uprv(nlvddi,nkn)
	real vprv(nlvddi,nkn)

	logical bgrad0
	logical blevel
	logical bdebug
	logical bdggu
	logical bout,bamb
	logical binside
	integer i,j,k,l
	integer ibtyp,igrad0
	integer ibc,ibcold
	integer nb,nlev,flev,ko
        integer ntbc,nlevko
	real dx,dy
	real scal,bc,weight,tweight
	real value,flag
	character*20 aline

	integer kint,kext,iunit
	integer ipext,ipint

	integer ifemopa

	bdebug = .true.
	bdebug = .false.

	ibcold = 0

	flag = flag_bndo

	bgrad0 = .false.
	blevel = .false.

	iunit = 730 + my_id
	kext = 6651
	kext = 0
	kint = ipint(kext)

	if( bdebug ) then
	  if( ndebug .eq. 0 ) then
	    ndebug = ifemopa('bndo_adjbc (91)','.bndo','form','unknown')
            call bndo_info(ndebug)
	  end if
	  call get_act_timeline(aline)
	  write(ndebug,*) 'bndo_adjbc ........... ',what,aline
	end if

	do i=1,nbndo

	  k = kbcnod(i)
	  nb = nopnod(i)
	  ibc = ibcnod(i)
	  bdebug = ( k == kint )

	  if( iopbnd(k) .ne. i ) then
	    write(6,*) i,iopbnd(k)
	    stop 'error stop bndo_adjbc: internal error (11c)'
	  end if

	  if( ibc .ne. ibcold ) then
	    call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	    call get_bnd_ipar(ibc,'igrad0',igrad0)
	    bgrad0 = igrad0 .gt. 0
	    blevel = ibtyp .eq. 1
	    ibcold = ibc
	  end if

          if( blevel ) then       !LEVELBC !DEBHELP

	  dx = xynorm(1,i)
	  dy = xynorm(2,i)
	  nlev = ilhkv(k)
	  flev = jlhkv(k)

	  do l=flev,nlev
	    scal = dx * uprv(l,k) + dy * vprv(l,k)
	    bout = scal .le. 0.				!outgoing flow
	    bamb = cv(l,k) == flag			!make ambient value
	    binside = bgrad0 .or. bout .or. bamb
	    if( binside ) then				!take from inside
	      bc = 0.
              ntbc = 0
              tweight = 0.
	      do j=1,nb
		ko = nopnodes(j,i)
                nlevko = ilhkv(ko)
                if( l .le. nlevko ) then        !LEVMX  !only if level exists
                  ntbc = ntbc + 1
		  weight = wopnodes(j,i)
                  tweight = tweight + weight
		  bc = bc + weight * cv(l,ko)
                end if
	      end do
	    else				!impose boundary value
              ntbc = nb                         !NTBC - BUG -> has not been set
	      tweight = 1.			!prob useless
	      bc = cv(l,k)
	    end if

            if( ntbc .gt. 0 ) then              !LEVMX  !at least 1 node found
	      value = bc / tweight
	      if( cv(l,k) .le. -5555. ) then	!differential value
		value = value - 10000. - cv(l,k)
	      end if
	      cv(l,k) = value
            else if( l .gt. 1 ) then            !take from above
              cv(l,k) = cv(l-1,k)
            else
	      write(6,*) i,k,ipext(k),nb,ntbc,ibc,binside
	      write(6,*) nlev,ntbc
              stop 'error stop bndo_adjbc: internal error (5)'
            end if

	  end do

          end if

	end do

	end

!***********************************************************************

	subroutine bndo_radiat(it,rzv)

! imposes radiation condition for levels

	use mod_bndo
	use mod_bound_geom
	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer it
	real rzv(nkn)

	logical bdebug
	integer i,j,k,l
	integer ibc,ibcold
	integer nb,nlev,ko
        integer ntbc,nlevko
	integer ibtyp
	real bc,weight,tweight

	integer ipext

	integer ifemopa

	bdebug = .true.
	bdebug = .false.
	ibcold = 0

	if( bdebug ) then
	  if( ndebug .eq. 0 ) then
	    ndebug = ifemopa('bndo_adjbc (91)','.bndo','form','unknown')
            call bndo_info(ndebug)
	  end if
	  write(ndebug,*) 'bndo_adjbc ........... ','radiat',it
	end if

	do i=1,nbndo

	  k = kbcnod(i)
	  nb = nopnod(i)
	  ibc = ibcnod(i)

	  if( ibc .ne. ibcold ) then
	    call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	    ibcold = ibc
	  end if
          !write(78,*) i,k,nb,ibc,ibtyp

	  if( iopbnd(k) .ne. i ) then
	    write(6,*) i,iopbnd(k)
	    stop 'error stop bndo_radiat: internal error (11d)'
	  end if

          if( ibtyp .eq. 31 ) then       !radiation condition only

	    bc = 0.
            ntbc = 0
            tweight = 0.
	    do j=1,nb
		ko = nopnodes(j,i)
                ntbc = ntbc + 1
		weight = wopnodes(j,i)
                tweight = tweight + weight
		bc = bc + weight * znv(ko)
	    end do

            if( ntbc .gt. 0 ) then         !LEVMX  !at least 1 node found
	      rzv(k) = bc / tweight
            else
              stop 'error stop bndo_radiat: internal error (5)'
            end if

          end if

	end do

	end

!***********************************************************************

