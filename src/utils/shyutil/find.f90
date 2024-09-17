
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-1999,2004-2005,2009,2011-2012  Georg Umgiesser
!    Copyright (C) 2014-2017,2019  Georg Umgiesser
!    Copyright (C) 2011  Marco Bajo
!    Copyright (C) 2011  Christian Ferrarin
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

! routines for interpolation onto regular grid
!
! contents :
!
! subroutine find_close_elem(ieold,xp,yp,ielem)
!		finds element for point (xp,yp) starting from ieold
! subroutine find_elem_from_old(ieold,xp,yp,ielem)
!		finds element for point (xp,yp) starting from ieold
! subroutine find_element(xp,yp,ielem)
!		finds element for point (xp,yp)
! function in_element(ie,xp,yp)
!		checks if point (xp,yp) is in element ie
!
! revision log :
!
! 18.11.1998	ggu	routine commented
! 18.11.1998	ggu	routine setgeo introduced
! 19.11.1998	ggu	routines a2char, prchar added
! 19.10.1999	ggu	routine mkmask added from subutl
! 25.11.2004	ggu	new routines femintp and elemintp for interpolation
! 14.03.2005	ggu	new routines for interpolation in element
! 11.03.2009	ggu	new helper routine getgeoflag()
! 12.06.2009	ggu	passing to double precision, intrid, bug bug_f_64bit
! 26.01.2011	ggu&mbj	handling extrapolation in am2av()
! 27.01.2011	ggu&ccf	bug fix in find_elem_from_old() BUG_27.01.2011
! 31.03.2011	ggu	new routine elemmask()
! 24.11.2011	ggu	new routine find_close_elem()
! 20.06.2012	ggu	new routine get_scal_elem()
! 07.10.2012	ggu	new routine av2fm()
! 10.10.2012	ggu	new routine fm2am2d() and fm2am3d()
! 26.10.2012	ggu	bug fix: do not access not existing storage
! 30.05.2014	ggu	in av2amk() do not interpolate for flag values
! 07.07.2014	ggu	new routine intp_reg()
! 25.09.2015	ggu	new routines intp_reg_nodes(), intp_reg_elems()
! 05.05.2016	ggu	file restructured (module)
! 14.05.2016	ggu	allow for extension of grid -> bregextend
! 23.06.2016	ggu	allow for eps in computing box
! 23.09.2016	ggu	allow for eps in computing box and reg intp
! 23.04.2017	ggu	new routine intp_reg_single_nodes()
! 25.05.2017	ggu	changed VERS_7_5_28
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 18.10.2019	ggu	cleaned contents
! 28.08.2024	ggu	new routines find_unique_element(), find_elements()
!
! notes :
!
! pxareg,pyareg         coordinates of lower left point of matrix
! pxdreg,pydreg         grid size of matrix
! pzlreg                value of z for land points
!
!******************************************************
!******************************************************
!******************************************************

	subroutine find_close_elem(ieold,xp,yp,ielem)

! finds element for point (xp,yp) starting from ieold
!
! uses data structure ev and ieltv
! honors material boundary

	use mod_geom
	use evgeom
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ieold
	real xp,yp
	integer ielem	!element number on return

	logical binit,bdebug
	integer ie,ii,iside,lmax,loop
	real xi,ximin
	double precision a(3),b(3),c(3)

	logical in_element

	bdebug = ieold == 73
	bdebug = .false.
	lmax = 10

	call is_init_ev(binit)

!-------------------------------------------------------------
! check if old element is given -> if not test all elements
!-------------------------------------------------------------

	if( bdebug ) write(6,*) 'ggu_xi (1) ',ieold
	if( ieold .le. 0 .or. ieold .gt. nel ) then
	  call find_element(xp,yp,ielem)
	  return
	end if

	if( bdebug ) write(6,*) 'ggu_xi (2) ',ieold
	if( .not. binit ) then
	  call find_elem_from_old(ieold,xp,yp,ielem)
	  return
	end if

!-------------------------------------------------------------
! start from old element
!-------------------------------------------------------------

	if( bdebug ) write(6,*) 'ggu_xi (3) ',ieold

	loop = 0
	ie = ieold
	do while( ie .gt. 0 )
	  iside = 0
	  ximin = 1.e+30
	  call xi_abc(ie,a,b,c)
	  do ii=1,3
	    xi = a(ii) + b(ii)*xp + c(ii)*yp
	    if( bdebug ) write(6,*) 'ggu_xiii ',ie,ii,xi
	    if( xi .lt. ximin ) then
	      ximin = xi
	      iside = ii
	    end if
	  end do
	  if( ximin .ge. 0. ) then
	    ielem = ie
	    return
	  end if
	  if( iside .le. 0 .or. iside .gt. 3 ) then
	    if( bdebug ) write(6,*) '******** ',iside,ie,ximin,xi
	    ie = 0
	  else
	    ie = ieltv(iside,ie)
	    if( bdebug ) write(6,*) 'ggu_xiii iterate',ie,iside,ximin
	  end if
	  loop = loop + 1
	  if( loop .gt. lmax ) ie = 0
	end do

	if( bdebug ) write(6,*) 'ggu_xi final',ie
	ielem = 0

	end

!******************************************************

	subroutine find_elem_from_old(ieold,xp,yp,ielem)

! finds element for point (xp,yp) starting from ieold
!
! finds elements also if behind material barrier

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	logical bdebug
	integer ieold
	real xp,yp
	integer ielem	!element number on return

	logical in_element
	integer iem,iep,iloop

	bdebug = ieold == 73
	bdebug = .false.

	if( bdebug ) write(6,*) 'starting find_elem_from_old'

!-------------------------------------------------------------
! check if old element is given -> if not test all elements
!-------------------------------------------------------------

	if( ieold .le. 0 .or. ieold .gt. nel ) then
	  call find_element(xp,yp,ielem)
	  return
	end if

!-------------------------------------------------------------
! check if in old element
!-------------------------------------------------------------

	if( in_element(ieold,xp,yp) ) then
	  ielem = ieold
	  return
	end if

!-------------------------------------------------------------
! start from old element going upwards and downwards
!-------------------------------------------------------------

	iloop = 0
	iem = ieold-1
	if( iem .lt. 1 ) iem = nel		!BUG_27.01.2011
	iep = ieold+1
	if( iep .gt. nel ) iep = 1		!BUG_27.01.2011

	do while( iem .ne. ieold .and. iep .ne. ieold )
	  iloop = iloop + 1
	  if( bdebug ) write(6,*) 'looking for elem ',iem,iep,iloop
	  if( in_element(iem,xp,yp) ) then
	    if( bdebug ) write(6,*) 'iem found: ',iem
	    ielem = iem
	    return
	  end if
	  iem = iem - 1
	  if( iem .lt. 1 ) iem = nel

	  if( in_element(iep,xp,yp) ) then
	    if( bdebug ) write(6,*) 'iep found: ',iep
	    ielem = iep
	    return
	  end if
	  iep = iep + 1
	  if( iep .gt. nel ) iep = 1
	end do

!-------------------------------------------------------------
! no element found
!-------------------------------------------------------------

	if( bdebug ) write(6,*) 'no element found for ',ieold

	ielem = 0

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!******************************************************

	subroutine find_element(xp,yp,ielem)

! finds element for point (xp,yp)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real xp,yp
	integer ielem	!element number on return


	integer ie
	logical in_element

	do ie=1,nel
	  if( in_element(ie,xp,yp) ) then
		  ielem = ie
		  return
	  end if
	end do

	ielem = 0

	end

!******************************************************

	subroutine find_unique_element(xp,yp,ielem)

! finds element for point (xp,yp)
! if there are more than one element it uses the one 
! with the highest external number

	use basin, only : nkn,nel,ngr,mbw,ipev

	implicit none

	real xp,yp
	integer ielem	!element number on return

	integer n
	integer, parameter :: ndim = 20
	integer ielems(ndim)
	integer i,imax,nmax,next

	ielem = 0

	n = ndim
	call find_elements(xp,yp,n,ielems)

	if( n == 0 ) return
	if( n == 1 ) then
	  ielem = ielems(1)
	  return
	end if

	nmax = 0
	imax = 0

	do i=1,n
	  next = ipev(ielems(i))
	  if( next <= 0 ) stop 'error stop find_unique_element: internal (1)'
	  if( next > nmax ) then
	    nmax = next
	    imax = i
	  end if
	end do

	if( imax == 0 ) stop 'error stop find_unique_element: internal (2)'
	ielem = ielems(imax)

	end

!******************************************************

	subroutine find_elements(xp,yp,n,ielems)

! finds all element for point (xp,yp)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real xp,yp
	integer n		!dim enter, elems found on return
	integer ielems(n)	!element numbers on return

	integer ie,ndim
	logical in_element

	ndim = n
	n = 0

	do ie=1,nel
	  if( in_element(ie,xp,yp) ) then
	    n = n + 1
	    if( n > ndim ) stop 'error stop find_elements: dimension'
	    ielems(n) = ie
	  end if
	end do

	end

!******************************************************

	function in_element(ie,xp,yp)

! checks if point (xp,yp) is in element ie

	use basin

	implicit none

	logical in_element
	integer ie
	real xp,yp

	integer ii,k,in
	real xmin,ymin,xmax,ymax
	real x(3),y(3)

	integer intri

	in_element = .false.

	do ii=1,3
	    k = nen3v(ii,ie)
	    x(ii) = xgv(k)
	    y(ii) = ygv(k)
	end do

	xmin = min(x(1),x(2),x(3))
	ymin = min(y(1),y(2),y(3))
	xmax = max(x(1),x(2),x(3))
	ymax = max(y(1),y(2),y(3))

	if( xp .ge. xmin .and. xp .le. xmax ) then
	  if( yp .ge. ymin .and. yp .le. ymax ) then
		in = intri(x,y,xp,yp)
		if( in .gt. 0 ) in_element = .true.
	  end if
	end if

	end

!******************************************************

