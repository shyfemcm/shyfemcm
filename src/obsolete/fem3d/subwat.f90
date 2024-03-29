
!--------------------------------------------------------------------------
!
!    Copyright (C) 1992,1994-1995,1997-1999,2001,2003  Georg Umgiesser
!    Copyright (C) 2009-2010,2014-2016,2019  Georg Umgiesser
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

! water volume, surface and flux routines
!
! contents :
!
! subroutine volnod(k,vol,dz)	inputs volume vol into finite volume (node) k
! subroutine zrise(dz)		inputs water distributed over total surface
! subroutine surel(ielem,dz)	inputs water distributed over surface
! subroutine connod(k,dz,con,coe)	changes concentration in node k
! subroutine conele(ielem,dz,con,coe)	changes concentration in element ielem
!
! subroutine volz3(k,dvol)                      inputs water volume
!
! subroutine volco0(k,lmax,nlvddi,s,area,vol,vol0,svol)
!       computes volume and total mass of concentration in column of node k
! subroutine volno0(k,lmax,nlvddi,s,dvol,dcon)
!       inputs concentration in finite volume (node) k
!
! revision log :
!
! 05.08.1992	ggu	$$ibtyp3 - implementation of ibtyp=3 (volnod) (connod)
! 19.01.1994	ggu	$$conz - implementation of concentration (vol/connod)
! 20.01.1994	ggu	$$lumpc - evaluate conz at node (vol/connod)
! 24.03.1994	ggu	from volnod (surel,conele)
! 07.04.1995	ggu	copied from volno3 (volz3)
! 06.08.1997	ggu	use zenv for water level (volz3)
! 04.12.1997	ggu	concentration not adjusted anymore (volnod,surel)
! 04.12.1997	ggu	concentration adjusted in own routine (connod,conele)
! 30.04.1998	ggu	routines from subflx.f merged into this file
! 20.06.1998	ggu	two dimensional common blocks regolarized (nen3v,...)
! 21.08.1998	ggu	xv eliminated
! 22.10.1999	ggu	file cleaned, added 3D routines
! 16.01.2001	ggu	new concentration unique for node (connod)
! 20.11.2001	ggu	input concentration with ambient value ($AMB0)
! 25.03.2003	ggu	new routine zrise
! 28.04.2009	ggu	links re-structured 
! 23.03.2010	ggu	changed v6.1.1
! 18.06.2014	ggu	changed VERS_6_1_77
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_2
! 19.01.2015	ggu	changed VERS_7_1_3
! 05.05.2015	ggu	changed VERS_7_1_10
! 05.06.2015	ggu	changed VERS_7_1_12
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 16.12.2015	ggu	changed VERS_7_3_16
! 28.04.2016	ggu	changed VERS_7_5_9
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 21.05.2019	ggu	changed VERS_7_5_62
! 
!*****************************************************************

	subroutine volnod(k,vol,dz)

! inputs volume vol into finite volume (node) k by changing water level z
!
! k	node where to input
! vol	volume to input
! dz	achieved water level change (return)

	use mod_geom
	use mod_hydro
	use evgeom
	use basin

	implicit none

! arguments
	integer k
	real vol,dz
! local
	real area
	integer ie,i,ii,nl
	integer ibase
	integer elems(maxlnk)

	integer ithis

	call get_elems_around(k,maxlnk,nl,elems)

	area=0.
	do i=1,nl
	  ie=elems(i)
	  if(ie.le.0) stop 'error stop volnod: internal error'
	  area=area+ev(10,ie)
	end do
	area=area*4.

	dz=vol/area

	do i=1,nl
	  ie=elems(i)
	  ii = ithis(k,ie)
	  zenv(ii,ie)=zenv(ii,ie)+dz
	  zeov(ii,ie)=zeov(ii,ie)+dz		!FIXME ?????
	end do

	end

!*****************************************************************

	subroutine zrise(dz)

! inputs water distributed over total surface
!
! dz		rise of water level to achieve

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

! arguments
	real dz
! local
	integer ie,ii,k

	do ie=1,nel
	  do ii=1,3
	    zenv(ii,ie)=zenv(ii,ie)+dz
	    zeov(ii,ie)=zeov(ii,ie)+dz
	  end do
	end do

	do k=1,nkn
	  znv(k) = znv(k) + dz
	  zov(k) = zov(k) + dz
	end do

	end

!*****************************************************************

	subroutine surel(ielem,dz)

! inputs water distributed over surface to element ie
!
! ielem		element to input water, 0: all elements
! dz		rise of water level to achieve

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

! arguments
	integer ielem
	real dz
! local
	integer ie,ii
	integer ie1,ie2

	if(ielem.le.0) then
	  ie1=1
	  ie2=nel
	else
	  ie1=ielem
	  ie2=ielem
	end if

	do ie=ie1,ie2
	  do ii=1,3
	    zenv(ii,ie)=zenv(ii,ie)+dz
	    zeov(ii,ie)=zeov(ii,ie)+dz
	  end do
	end do

	end

!*****************************************************************

	subroutine connod(k,dz,con,coe)

! changes concentration according to a water level rise dz with conz. con
! in node k on variable coe
!
! k	node where to input
! dz	water level change 
! con	concentraion of water injected
! coe	variable to change

	use mod_geom
	use mod_hydro
	use evgeom
	use basin

	implicit none

! arguments
	integer k
	real con,dz
	real coe(3,1)
! local
	integer ie,i,ii,nl
	integer ibase
	integer elems(maxlnk)
	real depth
	real area,vol,dvol,voltot,cnew
	real massold,massnew
	real dvoltot

	integer ithis

	call get_elems_around(k,maxlnk,nl,elems)

	voltot = 0.
	dvoltot = 0.
	cnew = 0.
	massold = 0.
	massnew = 0.

	do i=1,nl
	  ie=elems(i)
	  if(ie.le.0) stop 'error stop connod: internal error'
	  ii = ithis(k,ie)
	  area = 4. * ev(10,ie)
	  depth=hm3v(ii,ie)+zenv(ii,ie)-dz		!$$lumpc
	  vol = area * depth
	  dvol = area * dz
	  dvoltot = dvoltot + dvol
	  voltot = voltot + vol + dvol
	  massold = massold + coe(ii,ie)*vol
	  cnew = cnew + coe(ii,ie)*vol + con*dvol
	end do

	do i=1,nl
	  ie=elems(i)
	  ii = ithis(k,ie)
	  coe(ii,ie) = cnew / voltot
	end do

	massnew = cnew
!	write(6,*) 'ggu0: ',dvoltot
!	write(6,*) 'ggu1: ',massold,massnew,massnew-massold

	massnew = 0.
	do i=1,nl
	  ie=elems(i)
	  ii = ithis(k,ie)
	  area = 4. * ev(10,ie)
	  depth=hm3v(ii,ie)+zenv(ii,ie)		!$$lumpc
	  vol = area * depth
	  massnew = massnew + coe(ii,ie)*vol
	end do
!	write(6,*) 'ggu2: ',massold,massnew,massnew-massold

	end

!*****************************************************************

	subroutine conele(ielem,dz,con,coe)

! changes concentration according to water level rise dz with conz. con
! in element ielem on variable coe
!
! ielem		element to input water, 0: all elements
! dz		rise of water level
! con		conzentration of injected water
! coe		variable to change

	use mod_hydro
	use basin

	implicit none

! arguments
	integer ielem
	real dz,con
	real coe(3,1)
! local
	real depth
	integer ie,ii
	integer ie1,ie2

	if(ielem.le.0) then
	  ie1=1
	  ie2=nel
	else
	  ie1=ielem
	  ie2=ielem
	end if

	do ie=ie1,ie2
	  do ii=1,3
	    depth=hm3v(ii,ie)+zeov(ii,ie)
	    coe(ii,ie)=(depth*coe(ii,ie)+dz*con)/(depth+dz)
	  end do
	end do

	end

!*****************************************************************
!
	subroutine volz3(k,dvol)
!
! inputs water volume into finite volume
! ( 3d version )
!
! k	node
! dvol	water volume change (for whole time step)

	use mod_hydro
	use evgeom
	use basin

	implicit none

! arguments
	integer k
	real dvol
! local
        integer ie,ii
        real area,zz

        area=0.

        do ie=1,nel
          do ii=1,3
            if(nen3v(ii,ie).eq.k) then
                area=area+ev(10,ie)
            end if
          end do
        end do

        area=4.*area
	zz = dvol/area
	znv(k) = znv(k) + zz

        do ie=1,nel
          do ii=1,3
            if(nen3v(ii,ie).eq.k) then
		zenv(ii,ie) = zenv(ii,ie) + zz
            end if
          end do
        end do

	end

!***************************************************************

	subroutine volco0(k,lmax,nlvddi,s,area,vol,vol0,svol)

! computes volume and total mass of concentration in column of node k
! + volume of upper layer
! new version that computes only up to layer lmax (lmax > 0)

	use levels

	implicit none

! arguments
	integer k		!node defining column			(in)
	integer lmax		!maximum level to which compute		(in)
	integer nlvddi		!vertical dimension			(in)
	real s(nlvddi,1)	!variable (temperature, salinity,...)	(in)
	real area		!area of column 			(out)
	real vol		!total volume of column			(out)
	real vol0		!volume of first layer			(out)
	real svol		!total mass of s in column		(out)
! local
	integer l,ilevel,nlev
	integer mode
	real volume
! functions
	real volnode,areanode

	nlev = lmax
	if( lmax .le. 0 ) nlev = nlvddi		!all levels
	ilevel = min(nlev,ilhkv(k))

	mode = +1	!new time level
	vol=0.
	svol=0.

	do l=1,ilevel
	  volume = volnode(l,k,mode)
	  vol = vol + volume
	  svol = svol + volume * s(l,k)
	end do

	vol0 = volnode(1,k,mode)
	area = areanode(1,k)

	end

!******************************************************

	subroutine volno0(k,lmax,nlvddi,s,dvol,dcon)

! inputs concentration in finite volume (node) k
! ( 3d version ) -> only up to layer lmax

	use levels

	implicit none

! arguments
	integer k		!node defining volume
	integer lmax		!maximum level to which introduce
	integer nlvddi		!vertical dimension
	real s(nlvddi,1)	!variable (temperature, salinity,...)
	real dvol		!change in water volume (whole time step)
	real dcon		!concentration of new volume dvol
! local
	logical debug
	integer l,nlev
	real area,vol,vol0,svol
	real alpha,beta,gamma

	debug=.false.

        if( dcon .le. -990. ) return !input with ambient concentration $AMB0

	call volco0(k,lmax,nlvddi,s,area,vol,vol0,svol)

	nlev = ilhkv(k)
	if( lmax .gt. 0 .and. nlev .gt. lmax ) nlev = lmax

	alpha=dvol/(vol+dvol)
	beta=vol0/(vol0+dvol)
	gamma=dvol*svol/(vol*(vol0+dvol))
	gamma=dvol*dcon+alpha*(svol-dcon*vol)
	gamma=gamma/(vol0+dvol)

	do l=1,nlev
	  if(l.eq.1) then
	    s(l,k)=beta*s(l,k)+alpha*beta*(dcon-s(l,k))+gamma
	  else
	    s(l,k)=s(l,k)+alpha*(dcon-s(l,k))
	  end if
	  if(debug) write(6,*) k,l,s(l,k)
	end do

	end

!******************************************************

