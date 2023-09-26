
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2019  Georg Umgiesser
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
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!*********************************************************************
!
	subroutine setlst(ip,rkey,n,rflag)
!
	implicit none
!
! arguments
	integer ip(1),n
	real rkey(1),rflag
! common
	integer nmax,nins
	real rlast
	common /lstloc/ nmax,nins,rlast
! local
	integer i
!
	nmax=n
	rlast=rflag
	nins=0
!
	do i=1,n
	  ip(i)=0
	  rkey(i)=rflag
	end do
!
	return
	end
!
!*********************************************************************
!
	subroutine inslst(ip,rkey,ipact,rkact)
!
	implicit none
!
! arguments
	integer ip(*),ipact
	real rkey(*),rkact
! common
	integer nmax,nins
	real rlast
	common /lstloc/ nmax,nins,rlast
! local
	integer i
!
	if(rkact.le.rlast) return
!
	nins=nins+1
!
	do i=nmax-1,1,-1
		if(rkact.le.rkey(i)) goto 1
		rkey(i+1)=rkey(i)
		ip(i+1)=ip(i)
	end do
    1	continue
	rkey(i+1)=rkact
	ip(i+1)=ipact
!
	rlast=rkey(nmax)
!
	return
	end
!
	
!
!*********************************************************************
!
	subroutine maxlst(n)
!
	implicit none
!
! arguments
	integer n
! common
	integer nmax,nins
	real rlast
	common /lstloc/ nmax,nins,rlast
!
	n=nins
	if(nins.gt.nmax) n=nmax
!
	return
	end
