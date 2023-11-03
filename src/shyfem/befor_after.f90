
--------------------------------------------------------------------------
!
!    Copyright (C) 1994-1995,1997-2020  Georg Umgiesser
!    Copyright (C) 2008,2010,2014,2018  Christian Ferrarin
!    Copyright (C) 2016  William McKiver
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

! utility routines for modules (befor, after)
!
! contents :
!
!       subroutine do_befor
!       subroutine do_after
!
! revision log :
!
! 28.09.2023    gml     newly created from subnsh.f (iostr.f90)
! 29.10.2023    ggu     recovered log

!----------------------------------------------------------------------
        module befor_after
!----------------------------------------------------------------------

!----------------------------------------------------------------------
        contains
!----------------------------------------------------------------------

	subroutine do_befor

! to do in time loop before time step

	implicit none

	include 'modules.h'

	call modules(M_BEFOR)

	call tide_vuf
        call tideforce       !tidal potential !ccf
	call adjust_chezy

end subroutine do_befor

!********************************************************************

	subroutine do_after

! to do in time loop after time step

	implicit none

	include 'modules.h'

	double precision dtime

	call modules(M_AFTER)

	call get_act_dtime(dtime)

!	call wrouta
	call wrousa
!	call wrexta
	!call wrflxa
	!call wrvola(dtime)
	call wrboxa

        call resid
        call rmsvel

        call rst_write_restart

!        call tsmed
	call ts_shell

!	call wrnetcdf		!output in netcdf format - not supported

	call custom(dtime)

end subroutine do_after

!----------------------------------------------------------------------
        end module befor_after
!----------------------------------------------------------------------
