
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017,2019  Georg Umgiesser
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
! 05.12.2017	ggu	changed VERS_7_5_39
! 14.02.2019	ggu	changed VERS_7_5_56

	program zinit

! prepares initial z condition

	use basin

        include 'param.h'

!
	integer iapini
!
	include 'simul.h'

	real, allocatable :: hv(:)
!
!----------------------------------------------------------------
! k1,k2		nodes that define x/y of node 1/2
! z1,z2		initial z values at nodes 1/2
! mode		0: nothing  1: linear  2: cosine
!
!	data k1,k2 / 5,221 /
!	data k1,k2 / 5,167 /
!	data k1,k2 / 5,317 /
!	data z1,z2 / -.30,+.30 /
!	data z1,z2 / -2.30,+2.30 /
!	data z1,z2 / -1.30,+1.30 /
!	data mode /2/
!	data z1,z2 / .10,0. /
!	data mode /1/
!	data z1,z2 / -0.0475,+0.0475 /
!	data mode /1/
!
! bas092 -----------------------------
	data k1,k2 / 5,221 /
	data z1,z2 / -.30,+.30 /
	data mode /2/
! venlag61 -----------------------------
!	data k1,k2 / 4289,4340 /
!	data z1,z2 / -.40,+.40 /
!	data mode /2/
!----------------------------------------------------------------
!
	if(iapini(1,0,0,0).eq.0) stop

	allocate(hv(nkn))
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
	pi=acos(-1.)
!
	kk1 = ipext(k1)
	kk2 = ipext(k2)
	if( kk1 .le. 0 ) stop 'error stop: no node'
	if( kk2 .le. 0 ) stop 'error stop: no node'

	x1=xgv(kk1)
	y1=ygv(kk1)
	x2=xgv(kk2)
	y2=ygv(kk2)
	rl = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)
	zm = (z1+z2)/2.
	za = (z1-z2)/2.
!
	write(6,*) '------------ (kext,kint,z,x,y) -------------'
	write(6,1000) 'first  node : ',k1,kk1,z1,x1,y1
	write(6,1000) 'second node : ',k2,kk2,z2,x2,y2
 1000	format(1x,a,2i7,3f12.3)
!
	do i=1,nkn
	  xi=xgv(i)
	  yi=ygv(i)
!	  r varies linearily from 0 to 1 (r is distance from (x1,y1))
	  r = ( (xi-x1)*(x2-x1) + (yi-y1)*(y2-y1) ) / rl
	  if(mode.eq.1) then
!	    hv is linear -- hv(r=0)=z1 , hv(r=1)=z2
	    hv(i) = z1 + (z2-z1)*r
	  else if(mode.eq.2) then
!	    hv is cosinus -- hv(r=0)=z1 , hv(r=1)=z2 , hv(r=1/2)=(z1+z2)/2
	    hv(i) = zm  +  za * cos( r*pi )
	  else
!	    hv is average 
	    hv(i)=(z1+z2)*0.5
	  end if
	end do

	write(6,*) (hv(i),i=1,nkn)

	open(55,file='new.ini',status='unknown',form='unformatted')
	write(55) (hv(i),i=1,nkn)
	close(55)

	write(6,*) 'data written to file new.ini'

	end
