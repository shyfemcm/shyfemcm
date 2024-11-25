
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
	use levels
	use mod_depth
	use evgeom

	implicit none

	character*(*) slayers

	integer ie,l,lmax,llmax
	double precision thick,depth,hup
	double precision areae,area,volume,area_tot,volume_tot
	double precision, allocatable :: layers(:)
	real, save :: hflag = -999.

!-----------------------------------------------------------------
! prepare levels
!-----------------------------------------------------------------

	write(6,*) 'rrr: ',trim(slayers)
	lmax = 0
	call parse_levels(slayers,lmax,layers)
	allocate(layers(0:lmax))
	layers = 0.
	call parse_levels(slayers,lmax,layers(1:lmax))

	call create_levels(lmax,layers(1:lmax))

!-----------------------------------------------------------------
! compute volumes and areas
!-----------------------------------------------------------------

	write(6,*) 'areas and volumes of layers:'
	write(6,'(a)') ' layer     depth            area          volume'

	area_tot = 0.
	volume_tot = 0.

	do l=1,nlv
	  area = 0.
	  volume = 0.
	  hup = 0.
	  if( l > 1 ) hup = hlv(l-1)
	  do ie=1,nel
	    llmax = ilhv(ie)
	    if( l > llmax ) cycle
	    depth = hev(ie)
	    thick = hldv(l)
	    if( l == llmax ) thick = depth - hup
	    areae = 12.*ev(10,ie)
	    area = area + areae
	    volume = volume + thick*areae
	  end do
	  write(6,1000) l,layers(l),area,volume
	  volume_tot = volume_tot + volume
	  if( l == 1 ) area_tot = area
	end do

	write(6,1100) ' total          ',area_tot,volume_tot

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
 1000	format(i6,f10.2,2e16.8)
 1100	format(a16,2e16.8)
	end

!*******************************************************************

	subroutine parse_levels(string,lmax,layers)

	implicit none

	character*(*) string
	integer lmax
	double precision layers(lmax)

	integer n
	integer iscand

	write(6,*) 'hhh: ',trim(string)
	n = iscand(string,layers,lmax)
	write(6,*) 'hhh nnn: ',n,lmax
	if( lmax == 0 ) lmax = n

	end

!*******************************************************************

	subroutine create_levels(lmax,layers)

	use basin
	use levels
	use mod_depth

	implicit none

	integer lmax
	double precision layers(lmax)

	integer ie,l
	real depth
	real hup,hdown

        call makehev(hev)
	call levels_init(nkn,nel,lmax)
	hlv = 0.
	hldv = 0.
	hlv = layers(1:lmax)

	hup = 0.
	do l=1,lmax
	  hdown = hlv(l)
	  hldv(l) = hdown - hup
	  hup = hdown
	end do

	nlv = 0
	do ie=1,nel
	  depth = hev(ie)
	  do l=1,lmax
	    if( depth <= layers(l) ) exit
	  end do
	  if( l > lmax ) goto 99
	  ilhv(ie) = l
	  nlv = max(nlv,l)
	end do

	return
   99	continue
	write(6,*) 'depth = ',depth
	write(6,*) 'lmax,layers: ',lmax,layers
	write(6,*) 'not enough layers for handling depth of basin'
	write(6,*) 'add one more layer >= maximum depth of basin'
	stop 'error stop create_levels: not enough levels'
	end

!*******************************************************************

