
!--------------------------------------------------------------------------
!
!    Copyright (C) 2000,2010,2019  Georg Umgiesser
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

! routines for pre processing routines
!
! contents :
!
! subroutine sp190(nen,nel,nform,mbw)           determine bandwidth mbw
!
! revision log :
!
! 15.05.2000	ggu	subroutine sp191 removed
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2019	ggu	changed VERS_7_5_60
!
!**********************************************************

	subroutine sp190(nen,nel,nform,mbw)

! determine bandwidth mbw
!
! nen           node index
! nel           number of elements
! nform         number of nodes per element
! mbw           bandwidth (return value)

	dimension nen(1)

	mh=0

	do ie=1,nel
	do i=1,nform
	k=nen(nform*(ie-1)+i)
	if(k.gt.0) then
		do ii=i+1,nform
		kk=nen(nform*(ie-1)+ii)
		if(kk.gt.0) then
			mm=iabs(kk-k)
			if(mm.gt.mh) mh=mm
		end if
		end do
	end if
	end do
	end do

	mbw=mh

	end

!**********************************************************

