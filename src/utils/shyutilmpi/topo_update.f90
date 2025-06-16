
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003,2008-2009,2011,2015-2019  Georg Umgiesser
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

! topological set up routines
!
! contents :
!
! subroutine setnod
!			 sets array inodv
!
! revision log :
!
! 01.08.2003	ggu	created from sublnk.f
! 13.08.2003	ggu	in update_geom do not call setweg and setnod
! 06.11.2008	ggu	better error handling
! 06.04.2009	ggu	nlidim -> nlkdim
! 23.03.2010	ggu	changed v6.1.1
! 02.12.2011	ggu	print_bound_nodes() for writing out boundary nodes
! 09.12.2011	ggu	changed VERS_6_1_38
! 30.03.2012	ggu	changed VERS_6_1_51
! 12.09.2013	ggu	changed VERS_6_1_67
! 18.06.2014	ggu	changed VERS_6_1_77
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 04.05.2015	ggu	remove equivalence - use winkv as local array
! 20.05.2015	ggu	new call to mklenkii to set up lenkiiv
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 10.10.2015	ggu	changed VERS_7_3_2
! 16.12.2015	ggu	changed VERS_7_3_16
! 28.04.2016	ggu	changed VERS_7_5_9
! 26.09.2017	ggu	changed VERS_7_5_32
! 05.12.2017	ggu	changed VERS_7_5_39
! 24.01.2018	ggu	changed VERS_7_5_41
! 19.04.2018	ggu	changed VERS_7_5_45
! 16.02.2019	ggu	changed VERS_7_5_60
! 20.05.2020	ggu	new way to compute link structure (still experimental)
! 28.05.2020	ggu	some more changes in constructing link structure
! 13.04.2022	ggu	new call to make_links (ibound)
! 07.11.2024	ggu	updated call to update_ielt()
!
!*****************************************************************

	subroutine update_geom

! updates geometrical array (ieltv)

	use mod_geom
	use mod_geom_dynamic
	use basin

        implicit none

! local
        integer n

!-------------------------------------------------------------
! update ieltv
!-------------------------------------------------------------

        call update_ielt(nkn,nel,inodv,ieltv,nen3v)

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!*****************************************************************

	subroutine exchange_ieltv

! exchanges ieltv structure

	use mod_geom
	use mod_geom_dynamic
	use basin
	use shympi

	implicit none

        integer ie,ii
	integer iaux(3,nel)
	integer iiaux(nel)
	integer i,ia,ic,nc

!-------------------------------------------------------------
! start exchanging
!-------------------------------------------------------------

	call shympi_comment('exchanging ieltv')
	iiaux(:) = ieltv(1,:)
	call shympi_exchange_2d_elem(iiaux)
	iaux(1,:) = iiaux(:)
	iiaux(:) = ieltv(2,:)
	call shympi_exchange_2d_elem(iiaux)
	iaux(2,:) = iiaux(:)
	iiaux(:) = ieltv(3,:)
	call shympi_exchange_2d_elem(iiaux)
	iaux(3,:) = iiaux(:)
	call shympi_barrier

!-------------------------------------------------------------
! print info
!-------------------------------------------------------------

        write(my_unit,*) 'printing ghost elems: ' // 'ieltv'
        write(my_unit,*) 'n_ghost_areas = ',n_ghost_areas,my_id

	do ia=1,n_ghost_areas
          ic = ghost_areas(1,ia)
          nc = ghost_areas(4,ia)
          write(my_unit,*) 'elems outer: ',ic,nc
          do i=1,nc
            ie = ghost_elems_out(i,ia)
            write(my_unit,*) ie,ipev(ie),(ieltv(ii,ie),ii=1,3)
          end do
          nc = ghost_areas(5,ia)
          write(my_unit,*) 'elems inner: ',ic,nc
          do i=1,nc
            ie = ghost_elems_in(i,ia)
            write(my_unit,*) ie,ipev(ie),(ieltv(ii,ie),ii=1,3)
          end do
        end do

	do ie=1,nel
	  do ii=1,3
	    if( ieltv(ii,ie) /= iaux(ii,ie) ) then
	      write(my_unit,*) 'ieltv: ',ie,ieltv(ii,ie),iaux(ii,ie)
	    end if
	  end do
	end do

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!****************************************************************

        subroutine setnod

! sets (dynamic) array inodv
!
! inodv 
!	 0: internal node  
!	>0: open boundary node
!       -1: boundary node  
!	-2: out of system
!
! if open boundary node, inodv(k) is number of boundary (ggu 15.11.2001)

	use mod_geom_dynamic
	use evgeom
	use basin
	use shympi

        implicit none

	double precision, parameter :: winmax = 359.99
        integer ie,ii,k,n
        integer ibc,ibtyp
	integer nbc,ndry
        double precision winkv(nkn)

        integer ipext
	integer itybnd,nkbnds,kbnds,nbnds

! initialize array to hold angles

	ndry = 0
	winkv = 0.

! sum angles

        do ie=1,nel
          if(iwegv(ie).eq.0) then !element is in system
            do ii=1,3
              k=nen3v(ii,ie)
              winkv(k)=winkv(k)+ev(10+ii,ie)
            end do
	  else
	    ndry = ndry + 1
          end if
	end do

	!if( ndry > 0 ) then
	!  write(6,*) 'dry elements: ',ndry,' / ',nel
	!end if

        !call shympi_comment('shympi_elem: exchange winkv')
        call shympi_exchange_and_sum_2d_nodes(winkv)

! set up inodv

        do k=1,nkn
          if(winkv(k).gt.winmax) then     !internal node
            inodv(k)=0
          else if(winkv(k).eq.0.) then  !out of system
            inodv(k)=-2
          else                          !boundary node
            inodv(k)=-1
          end if
        end do

! now mark open boundary nodes

	nbc = nbnds()

	do ibc=1,nbc
          ibtyp = itybnd(ibc)
	  n = nkbnds(ibc)
          if(ibtyp.ge.3) then       !$$ibtyp3	!$$ibtyp4
            do ii=1,n
              k=kbnds(ibc,ii)
	      if( k <= 0 ) cycle
              if(inodv(k).eq.-1) inodv(k)=ibc
            end do
          else if(ibtyp.gt.0) then
            do ii=1,n
              k=kbnds(ibc,ii)
	      if( k <= 0 ) cycle
              if(inodv(k).ne.-1) goto 99
              inodv(k)=ibc
            end do
          end if
        end do

	!call shympi_comment('exchanging inodv')
	call shympi_exchange_2d_node(inodv)
	!call shympi_barrier

        return
   99   continue
        write(6,*) 'error for open boundary node'
	write(6,*) ibc,n,ii,k,ipext(k),inodv(k)
        if( inodv(k) .eq. 0 ) then
          write(6,*) 'the node is an inner node'
        else if( inodv(k) .eq. -2 ) then
          write(6,*) 'the node is not in the system (dry)'
        else
          write(6,*) 'the node has already been flagged as an'
          write(6,*) 'open boundary node (listed twice?)'
        end if
        write(6,*) 'external node number : ',ipext(k)
        write(6,*) 'boundary flag : ',inodv(k)
        stop 'error stop setnod: open boundary node'
        end

!****************************************************************
!****************************************************************
!****************************************************************

	function is_internal_node(k)

	use mod_geom_dynamic

	implicit none

	logical is_internal_node
	integer k

	is_internal_node = inodv(k) .eq. 0

	end

!****************************************************************

	function is_boundary_node(k)

	use mod_geom_dynamic

	implicit none

	logical is_boundary_node
	integer k

	is_boundary_node = inodv(k) .ne. 0 .and. inodv(k) .ne. -2

	end

!****************************************************************

	function is_open_boundary_node(k)

	use mod_geom_dynamic

	implicit none

	logical is_open_boundary_node
	integer k

	is_open_boundary_node = inodv(k) .gt. 0

	end

!****************************************************************

	function is_dry_node(k)

	use mod_geom_dynamic

	implicit none

	logical is_dry_node
	integer k

	is_dry_node = inodv(k) .eq. -2

	end

!*****************************************************************

