
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

!-----------------------------------------------------------
! integer variables
!-----------------------------------------------------------

!	ifunit			unit of flux file
!	itfold,itfnew,itfact	time for old/new/actual level
!	itperiod		periodic heat flux condition 
!				(0 if none, default)
!				(31536000 for one year - 365 days)
!	irhumid			relative humidity given 
!				(0 for wet bulb given)
!				(1 for relative humidity, default)

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2019	ggu	changed VERS_7_5_60

!-----------------------------------------------------------


        integer ifunit,itfold,itfnew,itfact,itperiod,irhumid
        common /qflxi/ ifunit,itfold,itfnew,itfact,itperiod,irhumid

!-----------------------------------------------------------
! values of old time level
!-----------------------------------------------------------

        real qsold,taold,tbold,uwold,ccold &
     &			,urold,pold,eold,rold,qold
        common /qflxro/ qsold,taold,tbold,uwold,ccold &
     &			,urold,pold,eold,rold,qold

!-----------------------------------------------------------
! values of new time level
!-----------------------------------------------------------

        real qsnew,tanew,tbnew,uwnew,ccnew &
     &			,urnew,pnew,enew,rnew,qnew
        common /qflxrn/ qsnew,tanew,tbnew,uwnew,ccnew &
     &			,urnew,pnew,enew,rnew,qnew

!-----------------------------------------------------------
! values of actual time level
!-----------------------------------------------------------

        real qsact,taact,tbact,uwact,ccact &
     &			,uract,pact,eact,ract,qact
        common /qflxra/ qsact,taact,tbact,uwact,ccact &
     &			,uract,pact,eact,ract,qact

!-----------------------------------------------------------
! save common blocks
!-----------------------------------------------------------

        save /qflxi/, /qflxro/, /qflxrn/, /qflxra/

!-----------------------------------------------------------
! end of header
!-----------------------------------------------------------

