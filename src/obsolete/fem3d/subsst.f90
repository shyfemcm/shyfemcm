
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

!**************************************************************************

! time scaling routines
!
! contents :
!
! function tconvs(itscal)		shell for tconv
! function tconv(scale1,scale2)		returns factor to scale time value
!
!*************************************************************
!
        function tconvs(itscal)
!
! shell for tconv
!
! itscal	switch to find out in what time unit
!		...has to be converted; it is assumed, that
!		...the original time unit is seconds
! tconv		factor time values have to be multiplied with
!		...(on error 1. is returned)
!
! possible values for itscal are :
!			0: sec   1: min    2:hour
!			3: day   4: month  5:year
!
        implicit none
!
	real tconvs
        integer itscal
        real tconv
!
	character*3 unit
!
	if(itscal.eq.0) then
		unit='sec'
	else if(itscal.eq.1) then
		unit='min'
	else if(itscal.eq.2) then
		unit='hou'
	else if(itscal.eq.3) then
		unit='day'
	else if(itscal.eq.4) then
		unit='mon'
	else if(itscal.eq.5) then
		unit='yea'
	end if
!
	tconvs=tconv('sec',unit)
!
	return
	end
!
!*************************************************************
!
        function tconv(scale1,scale2)
!
! returns factor to scale time value to an other unit
!
! scale1	unit of time
! scale2	unit to be converted in
! tconv		factor time values have to be multiplied with
!		...(on error 1. is returned)
!
! possible values for scale1/2 are :
!		sec,min,hour,day,month,year
!		sek,min,stunde,tag,monat,jahr
!		sec,min,ora,giorno,mese,anno
! abbreviations are possible (m=min,s=sec)
!
        implicit none
!
	real tconv
	character*(*) scale1,scale2
!
	character*2 scale
        integer i,iscale(2)
        real rscale(6)
!
	data rscale /1.,60.,3600.,86400.,2592000.,933120000./
!
	do i=1,2
!
	if(i.eq.1) then
		scale=scale1(1:2)
	else
		scale=scale2(1:2)
	end if
!
	call uplow(scale,'low')
!
	if(scale(1:1).eq.'s') then
		if(scale(2:2).eq.'t') then
			iscale(i)=3
		else
			iscale(i)=1
		end if
	else if(scale(1:1).eq.'m') then
		if(scale(2:2).eq.'e'.or.scale(2:2).eq.'o') then
			iscale(i)=5
		else
			iscale(i)=2
		end if
	else if(scale(1:1).eq.'y'.or.scale(1:1).eq.'j' &
     &			.or.scale(1:1).eq.'a') then
		iscale(i)=6
	else if(scale(1:1).eq.'d'.or.scale(1:1).eq.'t' &
     &			.or.scale(1:1).eq.'g') then
		iscale(i)=4
	else if(scale(1:1).eq.'h'.or.scale(1:1).eq.'o') then
		iscale(i)=3
	else
		iscale(i)=0
	end if
!
	end do
!
	if(iscale(1).ne.0.and.iscale(2).ne.0 &
     &				.and.iscale(1).ne.iscale(2)) then
		tconv=rscale(iscale(1))/rscale(iscale(2))
	else
		tconv=1.
	end if
!
	return
	end
