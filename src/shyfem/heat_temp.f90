
!--------------------------------------------------------------------------
!
!    Copyright (C) 2004,2006,2009-2011,2014-2016,2019  Georg Umgiesser
!    Copyright (C) 2016  Christian Ferrarin
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

! heat flux module (temperature module)
!
! contents :
!
! subroutine heat2t(dt,dh,qs,qrad,albedo,ts,tsnew)
!               computes new sea temperature
!
! revision log :
!
! 16.08.2004	ggu	heat2t copied from subqfxt.f
! 23.03.2006	ggu	changed time step to real
! 11.11.2009	ggu	new routine make_albedo(), pass albedo to heat2t
! 23.03.2010	ggu	changed v6.1.1
! 01.06.2011	ggu	use constant albedo
! 05.11.2014	ggu	changed VERS_7_0_5
! 05.12.2014	ggu	changed VERS_7_0_8
! 26.02.2015	ggu	changed VERS_7_1_5
! 04.05.2016	ccf	do not pass albedo into heat2t
! 25.05.2016	ggu	changed VERS_7_5_10
! 16.02.2019	ggu	changed VERS_7_5_60
!
!*****************************************************************************

        subroutine heat2t(dt,dh,qs,qrad,ts,tsnew)

! computes new sea temperature
!
! radiation is positive if into the water

	use heat_const

        implicit none

        real dt                 !time step
        real dh                 !layer depth
        real qs                 !solar radiation corrected
        real qrad               !other radiation
        real ts                 !old temperature
        real tsnew              !new temperature

        real ct
        real qseff

!--------------------------------------------------
! general constants
!--------------------------------------------------

        ct = cpw*rhow*dh          !heat capacity

!--------------------------------------------------
! new temperature
!--------------------------------------------------

        tsnew = ts + (qs+qrad)*dt/ct

!--------------------------------------------------
! end of routine
!--------------------------------------------------

        end

!*****************************************************************************

        subroutine make_albedo(temp,albedo)

        implicit none

        real temp       !water temperature
        real albedo     !albedo

        real albed0,albed4
        save albed0,albed4

        double precision dgetpar

        integer icall
        save icall
        data icall /0/

        if( icall .eq. 0 ) then
          albed0 = dgetpar('albedo')
          albed4 = dgetpar('albed4')
          icall = 1
        end if

        albedo = albed0

        if( temp .lt. 4. ) then
          albedo = albed4
          !albedo = 0.5
          !albedo = 0.75
        end if

        end

!*****************************************************************************

