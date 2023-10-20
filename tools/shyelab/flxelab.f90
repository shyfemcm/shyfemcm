
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-1999,2001,2003,2007-2008,2010  Georg Umgiesser
!    Copyright (C) 2015,2019  Georg Umgiesser
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
! 18.11.1998	ggu	check dimensions with dimnos
! 06.04.1999	ggu	some cosmetic changes
! 03.12.2001	ggu	some extra output -> place of min/max
! 09.12.2003	ggu	check for NaN introduced
! 07.03.2007	ggu	easier call
! 08.11.2008	ggu	do not compute min/max in non-existing layers
! 07.12.2010	ggu	write statistics on depth distribution (depth_stats)
! 06.05.2015	ggu	noselab started
! 05.06.2015	ggu	many more features added
! 05.10.2015	ggu	flxelab started
! 12.10.2015	ggu	changed VERS_7_3_3
! 16.02.2019	ggu	changed VERS_7_5_60
!
!**************************************************************

	program flxelab_main

	call flxelab

	end

!**************************************************************

