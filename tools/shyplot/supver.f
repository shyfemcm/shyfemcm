
!--------------------------------------------------------------------------
!
!    Copyright (C) 1991,1998,2010,2014-2019  Georg Umgiesser
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

!  Vertical velocities
! 
!  revision log :
! 
!  27.08.1991	ggu	(from scratch)
!  14.08.1998	ggu	w = 0 at open boundary nodes
!  20.08.1998	ggu	some documentation
!  23.03.2010	ggu	changed v6.1.1
!  26.03.2010	ggu	HACK for bwater -> compute always
!  22.04.2010	ggu	changed VERS_6_1_5
!  23.12.2014	ggu	changed VERS_7_0_11
!  19.01.2015	ggu	changed VERS_7_1_3
!  10.07.2015	ggu	changed VERS_7_1_50
!  17.07.2015	ggu	changed VERS_7_1_52
!  17.07.2015	ggu	changed VERS_7_1_80
!  20.07.2015	ggu	changed VERS_7_1_81
!  25.05.2016	ggu	changed VERS_7_5_10
!  12.01.2017	ggu	changed VERS_7_5_21
!  18.12.2018	ggu	changed VERS_7_5_52
!  13.03.2019	ggu	changed VERS_7_5_61
!  21.05.2019	ggu	changed VERS_7_5_62
!  21.10.2021	ggu	compute and keep in wauxv average on mid layer
! 
! ******************************************************************

        subroutine make_vertical_velocity

!  computes vertical velocities
! 
!  from sp256w in new3di.F
! 
!  velocities are computed on S/T points (top and bottom of layer)
!  bottom velocity of the whole column is assumed to be 0
!  -> maybe change this
! 
!  computes volume differences and from these computes vertical
!  velocities at every time step so that the computed velocities
!  satisfy the continuity equation for every single time step
! 
!  wlnv is computed horizontally at a node and vertically
!  it is at the center of the layer -> there are nlv velocities
!  computed
! 
!  b,c are 1/m, (phi is dimensionless)
!  aj is m**2
!  utlnv... is m**2/s
!  dvol is in m**3/s
!  vv is m**2 (area)
! 
!  wlnv is first used to accumulate volume difference -> dvol
!  at the end it receives the vertical velocity
! 
!  wlnv (dvol)   aux array for volume difference
!  vv            aux array for area

	use mod_hydro_plot
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin

        implicit none

	logical byaron
        integer k,ie,ii,kk,l
        integer ilevel
	integer inwater
        real aj,wbot,wtop,ff

!  initialize

	byaron = .true.
	byaron = .false.

        do k=1,nkn
          do l=0,nlv
            wauxv(l,k)=0.
            wlnv(l,k) = 0.
          end do
        end do

!  compute difference of velocities for each layer
! 
!  f(ii) > 0 ==> flux into node ii

	inwater = 0

        do ie=1,nel
         !if( bwater(ie) ) then           !FIXME	!not working
	  inwater = inwater + 1
          aj=4.*ev(10,ie)               !area of triangle / 3
          ilevel = ilhv(ie)
          do l=1,ilevel
            do ii=1,3
                kk=nen3v(ii,ie)
                ff = utlnv(l,ie)*ev(ii+3,ie) + vtlnv(l,ie)*ev(ii+6,ie)
	if( byaron .and. kk .eq. 2088 .and. l .eq. 8 ) then
		ff = ff + 10./(6.*3.*aj)
	end if
                wlnv(l,kk) = wlnv(l,kk) + 3. * aj * ff
                wauxv(l,kk)=wauxv(l,kk)+aj
            end do
          end do
         !end if
        end do

	!write(6,*) '******** inwater = ',inwater

!  from vel difference get absolute velocity (w_bottom = 0)
!        -> wlnv(nlv,k) is already in place !
!        -> wlnv(nlv,k) = 0 + wlnv(nlv,k)
!  w of bottom of last layer must be 0 ! -> shift everything up
!  wlnv(nlv,k) is always 0
! 
!  dividing dvol(m**3/s) by area (wauxv) gives vertical velocity

        do k=1,nkn
          wbot = 0.
          do l=nlv,1,-1
            wtop = wlnv(l,k)
            wlnv(l,k) = wbot
            wbot = wbot + wtop
          end do
          wlnv(0,k) = wbot
        end do

        do k=1,nkn
          do l=1,nlv
            if( wauxv(l,k) .gt. 0. ) then
              wlnv(l-1,k) = wlnv(l-1,k) / wauxv(l,k)
            end if
          end do
        end do

! average to mid layer

        do k=1,nkn
          do l=1,nlv
              wauxv(l-1,k) = 0.5 * ( wlnv(l-1,k) + wlnv(l,k) )
          end do
          wauxv(nlv,k) = 0.
        end do

!  set w to zero at open boundary nodes (new 14.08.1998)
! 
!  FIXME -> only for ibtyp = 1,2 !!!!

!         do k=1,nkn
!           if( inodv(k) .gt. 0 ) then    !open boundary node
!             do l=0,nlv
!                wlnv(l,k) = 0.
!             end do
!           end if
!         end do

        return
        end

! ******************************************************************

