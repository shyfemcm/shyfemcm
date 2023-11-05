
!--------------------------------------------------------------------------
!
!    Copyright (C) 1993,1998,2010,2015,2019  Georg Umgiesser
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

! tidal constants for the north sea model
!
! contents :
!
! subroutine roger(rzv,dt,iweich)	tides for nwe-shelf
! subroutine s6768(it,zov,nkn)		special output to files 67/68
!
! revision log :
!
! 27.10.1993	ggu	written
! 02.11.1993	ggu	$$TEST -> xgv/ygv/only 2 calls !!
! 05.11.1993	ggu	compo introduced to select components
! 18.06.1998	ggu	more documentation
! 23.03.2010	ggu	changed v6.1.1
! 19.01.2015	ggu	changed VERS_7_1_3
! 05.05.2015	ggu	changed VERS_7_1_10
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
!
!******************************************************************
!
	subroutine roger(rzv,irv,irb,dt,iweich)
!
! implements tides for north west shelf
!
! rzv		z-boundary condition for nodes
! irv		pointer into rzv
! irb		total number of boundary nodes (entries in irv)
! dt		time step in seconds
! iweich	switch:  
!			iweich=0 : first call, only set up   
!			iweich>0 : compute new zeta
!
! ndim		number of boundary nodes to be created
! ndimda	number of values read from file 
!
! compo		1 if component is to be used, 0 if not to be used
!
! file contains 8 tidal constituents
! for every partial tide (K1, M2...) two sets of values are given
! which are the amplification of the cosine and sin part of the signal
! for every data set 182 points are given
! 69 points for north boundary (west to east)
! 82 points for west boundary (north to south)
! 32 points for south boundary (west to east)
! ==> total of 183 points
! the two corner points (NW,SW) are counted twice, 
! so there is a total of 181 points
! the data set contains 182 points, because the NW points is given twice
! therefore: points 1 == point 70 in data set
! corner points:  NW -> 1 and 70 , SW -> 151
! end points: NE -> 69 , SE -> 182
!
! FEM: every point has been converted into one quadratic cell
! this cell has been divided into two triangles
! therefore the points of the border are 3 more
! ==> 181 + 3 = 184
!
	use basin

	implicit none
!
! arguments
	integer iweich
	integer irb,irv(1)
	real dt
	real rzv(1)
! parameters
	integer ndim,lanz,ndimda
	character*80 tiden
	real rad
	parameter (ndim=184,ndimda=182,lanz=8)
	parameter (tiden='/d4/georg/hamsom/daten/rogtides.data')
	parameter (rad=3.14159/180.)
! common (just for debug)	!$$TEST
	include 'param.h'
! local
	real hsing(lanz,ndim),hcosg(lanz,ndim)
	real speed(lanz)
	real phasin(lanz),gfact(lanz)
	real cosinc(lanz),sininc(lanz),costom(lanz),sintom(lanz)
	integer itpoin(ndim)
	integer compo(lanz)
	character*2 tides(lanz)
	integer l,n,i,j,kint,kaux
	integer ncall
	real tom,park,zet
! functions
	integer ipint,ipext
! save & data
	save ncall
	save hsing,hcosg,speed,phasin,gfact
	save cosinc,sininc,costom,sintom
	save itpoin,tides,compo
	data ncall	/ 0 /
	data phasin	/ &
     &			 275.9223, 258.7572, 311.9180,  51.2723 &
     &			,328.4926, 311.3276,   0.0000, 282.9712 &
     &			/
	data gfact	/ &
     &			   1.1702,   1.1702,   1.0000,   1.1052 &
     &			,  0.9665,   0.9665,   1.0000,   1.2902 &
     &			/
!                         q1  o1  p1  k1  n2  m2  s2  k2
!	data compo	/ 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 /
	data compo	/ 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 /	!m2
!	data compo	/ 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 /	!k1
!	data compo	/ 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 /	!o1
!
	ncall=ncall+1
!
! first call, set up
!
	if(iweich.eq.0) then
!
! read in tidal data (ndimda data values for each component)
!
	write(6,*)
	write(6,*) 'reading from file :'
	write(6,*) tiden
	open(1,file=tiden,status='old',form='formatted')
	do l=1,lanz
	  read(1,'(a2,10x,f12.7)') tides(l),speed(l)
	  read(1,'(8f10.6)') (hcosg(l,n),n=1,ndimda)
	  read(1,'(8f10.6)') (hsing(l,n),n=1,ndimda)
	  write(6,*) tides(l),'  ',speed(l),compo(l)
	end do
	close(1)
!
! set value for south west corner and last point (east) on south border
!
	do l=1,lanz
	  hcosg(l,ndim-1)=hcosg(l,ndimda)	!east extension of last point
	  hsing(l,ndim-1)=hsing(l,ndimda)
	  hcosg(l,ndim)=hcosg(l,151)		!SW corner from point to north
	  hsing(l,ndim)=hsing(l,151)
	end do
!
! set increment and initial phase
!
	do l=1,lanz
	  tom=speed(l)*rad*dt/3600.
	  cosinc(l)=cos(tom)
	  sininc(l)=sin(tom)
	  costom(l)=cos(phasin(l)*rad)
	  sintom(l)=sin(phasin(l)*rad)
	end do
!
! note that point 1 and point 70 are the same !!! (redundant)
! point 183 is extension of 182 (see befor)
!
	do n=1,ndim
	  if(n.le.69) then		!north (west to east)
	    i=n+1			!nw-corner not included
	    j=1
	  else if(n.le.151) then	!west (north to south)
	    i=1				!now nw-corner included
	    j=n-69
	  else if(n.le.183) then	!south (+extension) (west to east)
	    i=n-151+1
	    j=83			!sw-corner not included
	  else				!sw-corner (node number 184)
	    i=1
	    j=83
	  end if
	  itpoin(n)=ipint(100*i+j)
	end do
!
! convert to pointer into irv
!
	do n=1,ndim
	  kint=itpoin(n)
	  j=0
	  do i=1,irb
	    if(irv(i).eq.kint) j=i
	  end do
	  itpoin(n)=j
	end do
!
!ccccccccccccccccccccc
!
	write(6,*) 'boundary conditions (roger) : ',irb
	do n=1,ndim
	  kaux=itpoin(n)
	  if(kaux.gt.0) then
	    kint=irv(kaux)
	    write(6,*) kaux,ipext(kint),xgv(kint),ygv(kint)
	  end if
	end do
!
	else	!if(iweich.ne.0) then
!
! normal call - set zeta
!
!	if(ncall.gt.1) return	!$$TEST
!
	do l=1,lanz
	  park=costom(l)
	  costom(l)=costom(l)*cosinc(l)-sintom(l)*sininc(l)
	  sintom(l)=sintom(l)*cosinc(l)+  park   *sininc(l)
	end do
!
	end if
!
	do n=1,ndim
	 if(itpoin(n).gt.0) then
	  zet=0.
	  do l=1,lanz
	    if(compo(l).ne.0) then
	      zet=zet+gfact(l)* &
     &		(costom(l)*hcosg(l,n)+sintom(l)*hsing(l,n))
	    end if
	  end do
	  rzv(irv(itpoin(n)))=zet
!	  write(6,*) itpoin(n),ipext(itpoin(n)),zet
	 end if
	end do
!
	return
	end
!
!***********************************************************
!
	subroutine s6768(it,zov,nkn)
!
! special output
!
	implicit none
!
	integer ndim
	parameter (ndim=200)
!
	integer it,nkn
	real zov(1)
!
	logical bfirst
	integer irvt1,irvt2,i,ir
	integer imin,imax
	real zmin,zmax
	integer irv1(ndim),irv2(ndim)
!
	integer ipint,ipext
!
	save bfirst,irvt1,irvt2,irv1,irv2
	data bfirst /.true./
!
	if(bfirst) then
		bfirst=.false.
!
		ir=0
		do i=7003,303,-100
		  ir=ir+1
		  irv1(ir)=i
		end do
		write(6,*) 's6768 : ',ir
		do i=304,381,1
		  ir=ir+1
		  irv1(ir)=i
		end do
		write(6,*) 's6768 : ',ir
		do i=481,3381,100
		  ir=ir+1
		  irv1(ir)=i
		end do
		write(6,*) 's6768 : ',ir
		irvt1=ir
!
		ir=0
		do i=7005,505,-100
		  ir=ir+1
		  irv2(ir)=i
		end do
		write(6,*) 's6768 : ',ir
		do i=506,579,1
		  ir=ir+1
		  irv2(ir)=i
		end do
		write(6,*) 's6768 : ',ir
		do i=679,3379,100
		  ir=ir+1
		  irv2(ir)=i
		end do
		write(6,*) 's6768 : ',ir
		irvt2=ir
!
		do i=1,irvt1
		  irv1(i)=ipint(irv1(i))
		end do
		do i=1,irvt2
		  irv2(i)=ipint(irv2(i))
		end do
	else
		zmin=zov(1)
		zmax=zmin
		imin=1
		imax=imin
		do i=2,nkn
		  if(zov(i).gt.zmax) then
			imax=i
			zmax=zov(i)
		  end if
		  if(zov(i).lt.zmin) then
			imin=i
			zmin=zov(i)
		  end if
		end do
		write(6,*) it,ipext(imin),ipext(imax),zmin,zmax

		write(67,*) it,1067,irvt1
		write(67,*) (zov(irv1(i)),i=1,irvt1)
		write(68,*) it,1068,irvt2
		write(68,*) (zov(irv2(i)),i=1,irvt2)
	end if
!
	return
	end
