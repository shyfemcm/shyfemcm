
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

!--------------------------------------------------------------------------
!
! trace point routines
!
! revision log :
!
! 30.08.1995    ggu     $$AUST - austausch coefficient introduced
!
!--------------------------------------------------------------------------

!================================================================
        module mod_trace_point
!================================================================

	implicit none

        logical, parameter :: btrace = .false.		!enables trace point
        logical, save      :: btrace_do = .true.	!user trace point
        logical, parameter :: btrace_master = .false.	!only master

!================================================================
        contains
!================================================================

	subroutine trace_point(text)

	use shympi

	character*(*) text

	logical bwrite

	if( .not. btrace ) return
	if( .not. btrace_do ) return

	if( btrace_master ) then
	  bwrite = ( my_id == 0 )
	else
	  bwrite = .true.
	end if
	  
	if( bwrite ) then
	  write(6,*) trim(text),my_id
	  flush(6)
	end if

	call shympi_barrier

	end subroutine trace_point

!****************************************************************

	subroutine set_trace_point(bt)

	logical bt

	btrace_do = bt

	end subroutine set_trace_point

!================================================================
        end module mod_trace_point
!================================================================

