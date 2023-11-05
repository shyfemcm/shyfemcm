
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2013-2014,2018-2019  Georg Umgiesser
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
! 23.03.2010	ggu	changed v6.1.1
! 10.05.2013	ggu	changed VERS_6_1_64
! 28.01.2014	ggu	changed VERS_6_1_71
! 30.10.2014	ggu	changed VERS_7_0_4
! 18.12.2018	ggu	changed VERS_7_5_52
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!***************************************************************
!
! look in subcst for initialization
!
! cstinit	M_INIT		start of hp/ht
! cstcheck	M_CHECK		after read of STR file
! cstsetup	M_SETUP		after setup of basic arrays (ev,...)
!
! look in subnsh for read, write and others
!
! prilog	M_PRINT		after everything has been setup (before loop)
! pritst	M_TEST		only for debug
!
! dobefor	M_BEFOR		beginning of each time loop
! doafter	M_AFTER		end of each time loop
!
! nlsh2d	M_READ		read in STR file
!
!***************************************************************
!
! modules still to transform:
!
! wrouta
! resid
! rmsvel
! adjust_chezy (bottom friction)
! prwnds
! prarea (checy)
! prclos
! proxy, prlgr
! prbnds
! pripar
! biocos
! locous, locspc
!
!***************************************************************

	subroutine modules(mode)

! handles module framework

	implicit none

	integer mode		!mode of call

!------------------------------------------------------------
! modules
!------------------------------------------------------------

	call mod_ext(mode)		!extra nodes
!	call mod_ets(mode)		!extra time series nodes
!	call mod_box(mode)		!boxes
	call mod_flx(mode)		!fluxes through sections
!	call mod_vol(mode)		!volumes in areas

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	end

!***************************************************************

