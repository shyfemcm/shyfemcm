
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998,2010,2014,2019  Georg Umgiesser
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

! dummy subroutines
!
! revision log :
!
! 17.06.1998	ggu	nrdskp introduced (needed to skip section)
! 23.03.2010	ggu	changed v6.1.1
! 05.12.2014	ggu	changed VERS_7_0_8
! 16.02.2019	ggu	changed VERS_7_5_60
!
!***********************************************************

	subroutine colrd
	implicit none
	call nrdskp
	end

!***********************************************************

	subroutine legrd
	implicit none
	call nrdskp
	end

!***********************************************************
 
