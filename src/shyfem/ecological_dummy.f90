
!--------------------------------------------------------------------------
!
!    Copyright (C) 2008,2011-2012,2014,2016,2018-2019  Georg Umgiesser
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

! ecological dummy module
!
! revision log :
!
! 29.04.2008	ggu	bfm model integrated in main branch
! 18.02.2011	ggu	general framework for ecological model
! 17.02.2012	ggu	changed VERS_6_1_45
! 24.02.2012	ggu	changed VERS_6_1_46
! 19.12.2014	ggu	changed VERS_7_0_10
! 09.09.2016	ggu	changed VERS_7_5_17
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
!
!**************************************************************

        subroutine ecological_module

! general interface to ecological module

        implicit none

	integer ibfm,ibio
	real getpar

	integer, save :: icall = 0

	if( icall .eq. -1 ) return

        ibfm = nint(getpar('ibfm'))
	ibio = nint(getpar('ibio'))

	ibfm = 0
	if( ibfm .gt. 0 ) then
	  write(6,*) 'BFM module has not been linked'
	  write(6,*) 'ibfm = ',ibfm
	  write(6,*) 'You must enable this feature in Rules.make'
	  stop 'error stop ecological_module: ibfm'
	end if

	if( ibio .gt. 0 ) then
	  write(6,*) 'No ecological module has been linked'
	  write(6,*) 'ibio = ',ibio
	  write(6,*) 'You must enable this feature in Rules.make'
	  stop 'error stop ecological_module: ibio'
	end if

	icall = -1

        end

!**************************************************************

        subroutine write_restart_eco(iunit)
        implicit none
	integer iunit
	end
        subroutine skip_restart_eco(iunit)
        implicit none
	integer iunit
	end
        subroutine read_restart_eco(iunit)
        implicit none
	integer iunit
	call skip_restart_eco(iunit)
	end

!**************************************************************

