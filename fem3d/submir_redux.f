
!--------------------------------------------------------------------------
!
!    Copyright (C) 1997,2010,2012,2015,2019  Georg Umgiesser
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

! linear equation solvers
!
! contents :
!
!      SUBROUTINE LEQT1B (A,N,NLC,NUC,IA,B,M,IB,IJOB,XL,IER)
!      SUBROUTINE LEQT2B (A,N,NLC,NUC,IA,B,M,IB,IJOB,U,IU,XL,IER)
!      SUBROUTINE LEQ1PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,IER)
!      SUBROUTINE LEQ2PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,WK,IER)
!      SUBROUTINE LUDAPB (A,N,NC,IA,UL,IU,D1,D2,IER)
!      SUBROUTINE LUELPB (UL,B,N,NC,IA,X)
!      SUBROUTINE LUREPB (A,N,NC,IA,UL,IU,B,X,IDGT,RES,IER)
!      SUBROUTINE DEQT1B (A,N,NLC,NUC,IA,B,M,IB,IJOB,XL,IER)
!      SUBROUTINE DEQT2B (A,N,NLC,NUC,IA,B,M,IB,IJOB,U,IU,XL,IER)
!      SUBROUTINE DEQ1PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,IER)
!      SUBROUTINE DEQ2PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,WK,IER)
!      SUBROUTINE DUDAPB (A,N,NC,IA,UL,IU,D1,D2,IER)
!      SUBROUTINE DUELPB (UL,B,N,NC,IA,X)
!      SUBROUTINE DUREPB (A,N,NC,IA,UL,IU,B,X,IDGT,RES,IER)
!      SUBROUTINE UERTST (IER,NAME)
!      SUBROUTINE UGETIO (IOPT,NIN,NOUT)
!      SUBROUTINE GELB(R,A,M,N,MUD,MLD,EPS,IER)
!      SUBROUTINE DGELB(R,A,M,N,MUD,MLD,EPS,IER)
!      SUBROUTINE MCHB(R,A,M,N,MUD,IOP,EPS,IER)
!      SUBROUTINE DMCHB(R,A,M,N,MUD,IOP,EPS,IER)
!      subroutine loctst(i,j,n,m)
!      function locmy(i,j,n,m)
!      function locimm(i,j,n,m)
!      function loccer(i,j,n,m)
!      function locssp(i,j,n,m)
!      function locsps(i,j,n,m)
!
! revision log :
!
! 03.04.1997	ggu	general - compiler warnings for gfortran
! 24.05.1997	ggu	general - compiler warnings -> call to dp routines
! 23.03.2010	ggu	changed v6.1.1
! 30.03.2012	ggu	changed VERS_6_1_51
! 01.06.2012	ggu	changed VERS_6_1_53
! 18.12.2015	ggu	changed VERS_7_3_17
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 09.11.2021	ggu	reduced to esential routines
!
!*************************************************************************
!
        function locssp(i,j,n,m)
!
! access ssp routines
!
! (i,j)   position of element in square matrix (row,column)
! n       dimension of square matrix
! m       band width of square matrix
! locssp  position of element in band matrix
!
        implicit none
!
	integer locssp
        integer i,j,n,m
!
	if(i-j.gt.m.or.j-i.gt.m) then
	  locssp = 0
	else if(n.le.m) then	!this is for a full matrix
	  locssp = n*(i-1) + j
        else if(i.lt.m) then
          locssp = 2*m*i - m + j - m*(m+1)/2 + (m-i)*(m-i+1)/2
        else if(i.gt.n-m+1) then
          locssp = 2*m*i - m + j - m*(m+1)/2 &
     &                - (i-(n-m+1))*(i-(n-m))/2
        else
          locssp = 2*m*i - m + j - m*(m+1)/2
        end if
!
        return
        end
!
!*************************************************************************
!
        function locsps(i,j,n,m)
!
! access ssp routines (symmetric compressed storage mode)
!
! only main and upper diagonals - if an element in the lower
! ...diagonals is referenced, 0 is returned
!
! (i,j)   position of element in square matrix (row,column)
! n       dimension of square matrix
! m       band width of square matrix
! locsps  position of element in band matrix
!
! original formula : locsps = (1+m)*(i-1)+abs(j-i)+1
!
        implicit none
!
	integer locsps
        integer i,j,n,m
!
	if(i.gt.j.or.j-i.gt.m) then
	  locsps = 0
        else if(i.gt.n-m+1) then
          locsps = m*(i-1)+j &
     &                - (i-(n-m+1))*(i-(n-m))/2
        else
          locsps = m*(i-1)+j
        end if
!
        return
        end

!*************************************************************************

