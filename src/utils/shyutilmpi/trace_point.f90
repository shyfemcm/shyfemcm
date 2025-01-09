
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
! 15.03.2024    ggu     written from scratch
! 10.04.2024    ggu     some documentation
!
!--------------------------------------------------------------------------

!================================================================
        module mod_trace_point
!================================================================

	implicit none

        logical, parameter :: btrace_enable = .true.	!enables trace point
        logical, parameter :: btrace_master = .true.	!only master writes

        logical, save      :: btrace_do = .true.	!user trace point
        logical, save      :: btrace = btrace_enable	!final decsion

	integer, parameter :: iut = 999			!unit for writing tp

!================================================================
        contains
!================================================================

	subroutine trace_point(text)

	use shympi

	character*(*) text

	logical bwrite

	if( .not. btrace ) return

	if( btrace_master ) then
	  bwrite = ( my_id == 0 )
	else
	  bwrite = .true.
	end if
	  
	if( .not. btrace_master ) call shympi_barrier

	if( bwrite ) then
	  write(6,*) 'trace_point: ',trim(text),my_id
	  write(iut,*) 'trace_point: ',trim(text),my_id
	  flush(6)
	  flush(iut)
	end if

	if( .not. btrace_master ) call shympi_barrier

	end subroutine trace_point

!****************************************************************

	subroutine set_trace_point(bt)

	use shympi

	logical bt

	btrace_do = bt

	btrace = btrace_do .and. btrace_enable
	if( btrace_master .and. my_id /= 0 ) btrace = .false.

	end subroutine set_trace_point

!================================================================
        end module mod_trace_point
!================================================================

