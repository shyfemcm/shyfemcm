
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

! computes areas and volumes of layers

! revision log :
!
! 10.10.2015	ggu	changed VERS_7_3_2
! 25.05.2016	ggu	changed VERS_7_5_10
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 09.11.2024	ggu	copied, computes areas and volumes of layers

!**************************************************************************

	subroutine bas_layers(slayers)

	use basin
	use mod_depth
	use evgeom

	implicit none

	character*(*) slayers

	integer ie,l,lmax
	double precision hlev,depth,ldown,lup
	double precision areae,area,volume,area_tot,volume_tot
	double precision, allocatable :: levels(:)
	real, save :: hflag = -999.

!-------------------------------------------------------
! make depth values
!-------------------------------------------------------

        call makehev(hev)
        call makehkv_minmax(hkv,-1)     !use minimum depth

!-----------------------------------------------------------------
! prepare levels
!-----------------------------------------------------------------

	lmax = 0
	call make_levels(slayers,lmax,levels)
	allocate(levels(0:lmax))
	levels = 0.
	call make_levels(slayers,lmax,levels(1:lmax))

!-----------------------------------------------------------------
! compute volumes and areas
!-----------------------------------------------------------------

	write(6,*) 'areas and volumes of layers:'
	write(6,'(a)') ' layer            area          volume'

	area_tot = 0.
	volume_tot = 0.

	do l=1,lmax
	  area = 0.
	  volume = 0.
	  lup = levels(l-1)
	  ldown = levels(l)
	  do ie=1,nel
	    depth = hev(ie)
	    areae = 12.*ev(10,ie)
	    if( depth <= lup ) cycle
	    hlev = ldown - lup
	    if( depth < ldown ) hlev = depth - lup
	    area = area + areae
	    volume = volume + hlev*areae
	  end do
	  write(6,1000) l,area,volume
	  volume_tot = volume_tot + volume
	  if( l == 1 ) area_tot = area
	end do

	write(6,1100) ' total',area_tot,volume_tot

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
 1000	format(i6,2e16.8)
 1100	format(a6,2e16.8)
	end

!*******************************************************************

	subroutine make_levels(string,lmax,levels)

	implicit none

	character*(*) string
	integer lmax
	double precision levels(lmax)

	integer n
	integer iscand

	n = iscand(string,levels,lmax)
	if( lmax == 0 ) lmax = n

	end

!*******************************************************************

