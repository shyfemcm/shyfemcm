
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

! mpi check (debug) routines
!
! revision log :
!
! 08.10.2024	ggu	copied from mpi tvd

!==================================================================
        module mod_shympi_check
!==================================================================

	implicit none

! iu_debug writes debug output to fort.500
! iuunf_debug writes debug output to fort.501

	private

	logical, save :: bcheck = .true.	!allows for writing to fort.500
	integer, save :: iu_debug = 0		!writes fort.500 if > 0
	integer, save :: iuunf_debug = 501	!writes fort.501 if > 0
	integer, save :: ifreq_debug = 1	!frequency of debug in fort.500

	logical, save :: bdebout = .false.	!allows for writing to fort.500
	logical, save :: bform = .false.	!allows for writing to fort.500
	logical, save :: bunform = .false.	!allows for writing to fort.500
	logical, save :: bout = .false.		!allows for writing to fort.500
	logical, save :: b2d = .false.		!allows for writing to fort.500

	integer, save :: nsize_check = 0	!nkn or nel
	integer, save :: nlvdi_check = 0	!nlv
	logical, save :: belem_check = .false.	!element or node
	real, save, allocatable :: mpi_check_aux(:,:)
	logical, save, allocatable :: blevel(:)

	! the following are internal routines - do not use directly

	!public :: mpi_check_initialize
	!public :: mpi_check_accum
	!public :: mpi_check_set
	!public :: mpi_check_finalize
	!public :: mpi_check_all

	! the following are public routines

	public :: mpi_check_node	!call mpi_check_node(what,is,array)
	public :: mpi_check_elem	!call mpi_check_elem(what,is,array)

        INTERFACE mpi_check_node
        MODULE PROCEDURE   mpi_check_node_2d &
     &                   , mpi_check_node_3d
        END INTERFACE

        INTERFACE mpi_check_elem
        MODULE PROCEDURE   mpi_check_elem_2d &
     &                   , mpi_check_elem_3d
        END INTERFACE

!==================================================================
        contains
!==================================================================

	subroutine mpi_check_initialize(dtime,what,isact,belem,b2dim)

	use basin
	use levels
	use shympi

	implicit none

	double precision dtime
	character*(*) what
	integer isact
	logical belem
	logical, optional :: b2dim

	integer l,ld,lmax,lmax_form,lmax_unf,nsize,nsize_global
	integer, save :: icall = 0
	character*80 what_aux

	if( .not. bcheck ) return

	bform = iu_debug > 0
	bunform = iuunf_debug > 0
	bout = bform .or. bunform

	b2d = .false.
	if( present(b2dim) ) b2d = b2dim

	if( .not. bout ) return

	!---------------------------------------------------------
	! allocate arrays
	!---------------------------------------------------------

	nsize = nkn
	nsize_global = nkn_global
	if( belem ) then
	  nsize = nel
	  nsize_global = nel_global
	end if
	if( allocated(mpi_check_aux) ) deallocate(mpi_check_aux)
	if( allocated(blevel) ) deallocate(blevel)
	nlvdi_check = nlvdi
	if( b2d ) nlvdi_check = 1
	allocate( mpi_check_aux(nlvdi_check,nsize) )
	allocate( blevel(nlv_global) )
	mpi_check_aux = 0.
	blevel = .false.

	nsize_check = nsize
	belem_check = belem

	icall = icall + 1

	!---------------------------------------------------------
	! decide if output has to be written
	!---------------------------------------------------------

	if( ifreq_debug < 0 ) then
	  bdebout = .false.
	else if( ifreq_debug > 0 ) then
	  bdebout = .false.
	  if( mod(icall,ifreq_debug) == 0 ) bdebout = .true.
	else
	  bdebout = .true.
	end if

	if( .not. bdebout ) return

	!---------------------------------------------------------
	! choose levels to be written
	!---------------------------------------------------------

	lmax_unf = nlv_global
	lmax_form = 0
	ld = 1
	if( lmax_unf > 1 ) ld = lmax_unf/2

	lmax = 0
	blevel = .false.
	do l=1,lmax_unf,ld
	  blevel(l) = .true.
	  lmax = lmax + 1
	end do
	lmax_form = lmax

	if( b2d ) then
	  lmax_unf = 1
	  lmax_form = 1
	  blevel = .false.
	  blevel(1) = .true.
	end if

	!---------------------------------------------------------
	! open file and write header
	!---------------------------------------------------------

	if( shympi_is_master() ) then
	  write(iu_debug,'(a,f14.2,3i8,l2,2a)') &
     &				'mpi_check: ',dtime &
     &				,isact,nsize_global,lmax_form &
     &				,belem,'  ',trim(what)
	  if( bform ) then
	    if( icall == 0 ) write(iu_debug,*) 'written by mpi_check_initialize'
	    write(iu_debug,'(a,f14.2,3i8,l2,2a)') &
     &				'time = ',dtime &
     &				,isact,nsize_global,lmax_form &
     &				,belem,'  ',trim(what)
	    flush(iu_debug)
	  end if
	  if( bunform ) then
	    what_aux(1:80) = what
	    write(iuunf_debug) &
     &				dtime &
     &				,isact,nsize_global,lmax_unf &
     &				,belem,what_aux
	    flush(iuunf_debug)
	  end if
	end if

	!---------------------------------------------------------
	! end of routine
	!---------------------------------------------------------

	end

!******************************************************************

	subroutine mpi_check_accum(i,l,conu)

	use shympi

	implicit none

	integer i,l
	real conu

	if( .not. bcheck ) return
	if( .not. bout ) return

	mpi_check_aux(l,i) = conu

	end

!******************************************************************

	subroutine mpi_check_set(array)

	use basin
	use levels
	use shympi

	implicit none

	real array(nlvdi_check,nsize_check)

	if( .not. bcheck ) return
	if( .not. bout ) return

	mpi_check_aux(:,:) = array(:,:)

	end

!******************************************************************

	subroutine mpi_check_finalize

	use basin
	use levels
	use shympi

	implicit none

	logical belem
	integer l,lmax,i,i_ext
	integer nsize,nsize_global
	real, allocatable :: aux(:),auxg(:)
	integer, save :: icall = 0
	integer, save :: nrec = 0
	integer, save :: nrec_form = 0
	integer, save :: nrec_unf = 0

	!write(6,*) 'bform:',bform,bunform,bout,bcheck,bdebout

	if( .not. bcheck ) return
	if( .not. bdebout ) return
	if( .not. bout ) return

	!bout = .false.		!if false checks how may records are written
	!bout = .true.		!if false checks how may records are written
	icall = icall + 1

	if( .not. bdebout ) return

	belem = belem_check
	nsize = nkn
	nsize_global = nkn_global
	if( belem ) then
	  nsize = nel
	  nsize_global = nel_global
	end if
	allocate(aux(nsize))
	allocate(auxg(nsize_global))

	lmax = nlv_global
	if( b2d ) lmax = 1

	!write(6,*) 'bform:',bform,bunform,bout

	do l=1,lmax
	  auxg = 0.
	  aux = 0.
	  if( l <= nlv ) aux(:) = mpi_check_aux(l,:)
	  call shympi_l2g_array(aux,auxg)
	  if( shympi_is_master() ) then
	    nrec_unf = nrec_unf + 1
	    if( bform .and. blevel(l) ) then
	      nrec_form = nrec_form + 1
	      write(iu_debug,*) 'level = ',l,nsize_global,nrec_form,icall
	    end if
	    if( bunform ) write(iuunf_debug) l,nsize_global,nrec_unf,icall
	    if( .not. bout ) cycle
	    if( bform ) then
	    do i=1,nsize_global
	      if( belem ) then
	        i_ext = ip_ext_elem(i)
	      else
	        i_ext = ip_ext_node(i)
	      end if
	      if( bform .and. blevel(l) ) then
	        write(iu_debug,*) i,i_ext,auxg(i)
	      end if
	      !if( bunform ) write(iuunf_debug) i,i_ext,auxg(i)
	    end do
	    end if
	    if( bunform ) then
	      if( belem ) then
	        write(iuunf_debug) ip_ext_elem
	      else
	        write(iuunf_debug) ip_ext_node
	      end if
	      write(iuunf_debug) auxg
	    end if
	  end if
	end do

	if( shympi_is_master() ) then
	  if( bform ) flush(iu_debug)
	  if( bunform ) flush(iuunf_debug)
	end if

	end

!******************************************************************

	subroutine mpi_check_all(dtime,what,isact,belem,array,bdim)

	use basin
	use levels
	use shympi

	implicit none

	double precision dtime
	character*(*) what
	integer isact
	logical belem
	real array(nlvdi,nsize_check)
	!real array(:,:)
	logical, optional :: bdim

	logical bbb

	bbb = .false.
	if( present(bdim) ) bbb = bdim

	call mpi_check_initialize(dtime,what,isact,belem,bbb)
	call mpi_check_set(array)
	call mpi_check_finalize

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine mpi_check_elem_3d(what,isact,array)

	implicit none

	character*(*) what
	integer isact
	real array(:,:)

	logical, save :: belem = .true.
	logical, save :: b2d = .false.
	double precision dtime

	call get_act_dtime(dtime)
	call mpi_check_all(dtime,what,isact,belem,array,b2d)

	end

!******************************************************************

	subroutine mpi_check_node_3d(what,isact,array)

	implicit none

	character*(*) what
	integer isact
	real array(:,:)

	logical, save :: belem = .false.
	logical, save :: b2d = .false.
	double precision dtime

	call get_act_dtime(dtime)
	call mpi_check_all(dtime,what,isact,belem,array,b2d)

	end

!******************************************************************

	subroutine mpi_check_elem_2d(what,isact,array)

	implicit none

	character*(*) what
	integer isact
	real array(:)

	logical, save :: belem = .true.
	logical, save :: b2d = .true.
	double precision dtime

	call get_act_dtime(dtime)
	call mpi_check_all(dtime,what,isact,belem,array,b2d)

	end

!******************************************************************

	subroutine mpi_check_node_2d(what,isact,array)

	implicit none

	character*(*) what
	integer isact
	real array(:)

	logical, save :: belem = .false.
	logical, save :: b2d = .true.
	double precision dtime

	call get_act_dtime(dtime)
	call mpi_check_all(dtime,what,isact,belem,array,b2d)

	end

!==================================================================
        end module mod_shympi_check
!==================================================================

