
!--------------------------------------------------------------------------
!
!    Copyright (C) 2004,2010-2011,2015,2019  Georg Umgiesser
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

! routines dealing with histogram
!
! contents :
!
! subroutine histo_init(nbin,bin0,dbin,rbin)
! subroutine histo_insert(value)
! subroutine histo_final(ic)
! 
! revision log :
!
! 16.03.2004	ggu	routines written from scratch
! 23.03.2010	ggu	changed v6.1.1
! 31.05.2011	ggu	changed VERS_6_1_23
! 09.01.2015	ggu	changed VERS_7_0_12
! 16.02.2019	ggu	changed VERS_7_5_60
! 28.09.2023	ggu	include substituted with module
! 22.09.2024	ggu	revisited - allocate arrays
!
!****************************************************************

!================================================================
	module mod_histo
!================================================================

	implicit none

	private

        integer, save :: ncbin = 0
        integer, allocatable, save :: icount(:)
        real, allocatable, save :: abin(:)

        INTERFACE histo_init
        MODULE PROCEDURE histo_init_auto, histo_init_bins
        END INTERFACE

	public :: histo_init,histo_insert,histo_get_bins,histo_final

!================================================================
	contains
!================================================================

	subroutine histo_alloc(nbin)

	integer nbin

	if( ncbin /= 0 ) deallocate(icount,abin)
	allocate(icount(ncbin+1),abin(ncbin))
	ncbin = nbin
	icount = 0

	end

!****************************************************************

        subroutine histo_init_auto(nbin,bin0,dbin)

! sets up icount and abin

        implicit none

        integer nbin            !total number of bins
        real bin0               !first bin (limit)
        real dbin               !regular bin size (0 => use rbin)

        integer i

	call histo_alloc(nbin)

        do i=1,nbin
          abin(i) = bin0 + (i-1) * dbin
        end do

        end

!****************************************************************

        subroutine histo_init_bins(nbin,rbin)

! sets up icount and abin

        implicit none

        integer nbin            !total number of bins
        real rbin(nbin)         !bin size limits (upper)

	call histo_alloc(nbin)

	abin = rbin

        end

!****************************************************************

        subroutine histo_insert(value)

        implicit none

        real value

        integer i

        do i=1,ncbin
          if( value .le. abin(i) ) then
            icount(i) = icount(i) + 1
            return
          end if
        end do

        i = ncbin+1
        icount(i) = icount(i) + 1

        end

!****************************************************************

        subroutine histo_get_bins(ac)

        implicit none

        real ac(ncbin)

	ac = abin

        end

!****************************************************************

        subroutine histo_final(ic)

        implicit none

        integer ic(ncbin+1)

	ic = icount

        end

!================================================================
	end module mod_histo
!================================================================

