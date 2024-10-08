
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
        module shympi_check
!==================================================================

	implicit none

! iu_debug writes debug output to fort.500
! iuunf_debug writes debug output to fort.501

	logical, save :: bcheck = .false.	!allows for writing to fort.500
	logical, save :: bdebout = .false.	!allows for writing to fort.500
	integer, save :: iu_debug = 500		!writes fort.500 if > 0
	integer, save :: iuunf_debug = 501	!writes fort.501 if > 0
	integer, save :: ifreq_debug = 1	!frequency of debug in fort.500
	real, save, allocatable :: mpi_check_aux(:,:)
	logical, save, allocatable :: blevel(:)

!==================================================================
        contains
!==================================================================


!==================================================================
        end module shympi_check
!==================================================================

	subroutine mpi_check_initialize(dtime,what,isact)

	use basin
	use levels
	use shympi
	use shympi_check

	implicit none

	double precision dtime
	character*(*) what
	integer isact
	integer l,lmax,lmax_form,lmax_unf
	integer itime
	integer, save :: icall = 0
	character*80 what_aux

	if( .not. bcheck ) return
	if( iu_debug <= 0 ) return

	if( icall == 0 ) then
	  allocate( mpi_check_aux(nlvdi,nel) )
	  allocate( blevel(nlv_global) )
	end if
	mpi_check_aux = 0.

	icall = icall + 1

	if( ifreq_debug < 0 ) then
	  bdebout = .false.
	else if( ifreq_debug > 0 ) then
	  bdebout = .false.
	  if( mod(icall,ifreq_debug) == 0 ) bdebout = .true.
	else
	  bdebout = .true.
	end if

	lmax_unf = 0
	lmax_form = 0
	if( bdebout ) then	!we write the output
	  lmax_unf = nlv_global
	  lmax = 0
	  blevel = .false.
	  do l=1,lmax_unf,lmax_unf/2
	    blevel(l) = .true.
	    lmax = lmax + 1
	  end do
	  lmax_form = lmax
	end if

	if( shympi_is_master() ) then
	  if( icall == 0 ) write(iu_debug,*) 'written by check_debug_initialize'
	  write(iu_debug,*) 'time = ',dtime,isact,lmax_form,'  ',trim(what)
	  flush(iu_debug)
	  what_aux(1:80) = what
	  itime = nint(dtime)
	  write(iuunf_debug) dtime,isact,lmax_unf,what_aux
	  flush(iuunf_debug)
	end if

	end

!******************************************************************

	subroutine mpi_check_accum(ie,l,conu)

	use shympi
	use shympi_check

	implicit none

	integer ie,l
	real conu

	if( .not. bcheck ) return
	if( iu_debug <= 0 ) return

	mpi_check_aux(l,ie) = conu

	end

!******************************************************************

	subroutine mpi_check_set(conu)

	use basin
	use levels
	use shympi
	use shympi_check

	implicit none

	real conu(nlvdi,nel)

	if( .not. bcheck ) return
	if( iu_debug <= 0 ) return

	mpi_check_aux(:,:) = conu(:,:)

	end

!******************************************************************

	subroutine check_check_finalize

	use basin
	use levels
	use shympi
	use shympi_check

	implicit none

	logical bout
	integer l,lmax,ie,ie_ext
	real, allocatable :: aux(:),auxg(:)
	integer, save :: icall = 0
	integer, save :: nrec = 0
	integer, save :: nrec_form = 0
	integer, save :: nrec_unf = 0

	if( .not. bcheck ) return
	if( .not. bdebout ) return
	if( iu_debug <= 0 ) return

	bout = .false.		!if false checks how may records are written
	bout = .true.		!if false checks how may records are written
	icall = icall + 1

	if( .not. bdebout ) return

	allocate(aux(nel))
	allocate(auxg(nel_global))

	lmax = nlv_global

	do l=1,lmax
	  auxg = 0.
	  aux = 0.
	  if( l <= nlv ) aux(:) = mpi_check_aux(l,:)
	  call shympi_l2g_array(aux,auxg)
	  if( shympi_is_master() ) then
	    nrec_unf = nrec_unf + 1
	    if( blevel(l) ) then
	      nrec_form = nrec_form + 1
	      write(iu_debug,*) 'level = ',l,nel_global,nrec_form,icall
	    end if
	    write(iuunf_debug) l,nel_global,nrec_unf,icall
	    if( .not. bout ) cycle
	    do ie=1,nel_global
	      ie_ext = ip_ext_elem(ie)
	      if( blevel(l) ) then
	        write(iu_debug,*) ie,ie_ext,auxg(ie)
	      end if
	      write(iuunf_debug) ie,ie_ext,auxg(ie)
	    end do
	  end if
	end do

	if( shympi_is_master() ) then
	  flush(iu_debug)
	  flush(iuunf_debug)
	end if

	end

!******************************************************************

