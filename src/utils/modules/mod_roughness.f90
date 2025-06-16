
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2016,2019  Georg Umgiesser
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

! revision log :
!
! 10.07.2015	ggu	changed VERS_7_1_50
! 28.04.2016	ggu	changed VERS_7_5_9
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 22.09.2020    ggu     correct warnings for PGI compiler

!**************************************************************************

!============================================================
	module mod_roughness
!============================================================

	implicit none

	integer, private, save :: nkn_roughness = 0

	real, allocatable, save :: z0bk(:)	!bottom roughness on nodes
	real, allocatable, save :: z0sk(:)	!surface roughness on nodes

	double precision, parameter :: z0bmin = 1.e-4
	double precision, parameter :: z0smin = 0.02

	real, parameter :: z0bini = 0.03*0.03
	real, parameter :: z0sini = 0.02

!============================================================
	contains
!============================================================

	subroutine mod_roughness_init(nkn)

	integer nkn

	if( nkn == nkn_roughness ) return

	if( nkn_roughness > 0 ) then
	  deallocate(z0bk)
	  deallocate(z0sk)
	end if

	nkn_roughness = nkn

	if( nkn == 0 ) return

	allocate(z0bk(nkn))
	allocate(z0sk(nkn))

	z0bk = z0bini
	z0sk = z0sini

	end subroutine mod_roughness_init

!============================================================
	end module mod_roughness
!============================================================

