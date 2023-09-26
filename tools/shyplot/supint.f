
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999,2002,2004-2005,2008-2011,2013-2016  Georg Umgiesser
!    Copyright (C) 2018-2019  Georg Umgiesser
!    Copyright (C) 2012  Christian Ferrarin
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

!  interactive routines for plotsim
! 
!  revision log :
! 
!  12.02.1999	ggu	adapted to auto mode
!  29.01.2002	ggu	new routine getisec()
!  17.03.2004	ggu	new routine okvar()
!  02.03.2005	ggu	new routines set_flag and get_flag
!  17.09.2008	ggu	comments for level = -1
!  06.12.2008	ggu	in extlev set not-existing values to flag
!  14.09.2009	ggu	new way to determine if section plot in getisec()
!  23.03.2010	ggu	changed v6.1.1
!  18.08.2011	ggu	make vsect bigger
!  31.08.2011	ggu	new plotting eos
!  01.09.2011	ggu	changed VERS_6_1_32
!  23.02.2012	ccf	allow plotting also for last layer
!  13.06.2013	ggu	scans varnam to decide what to plot
!  19.06.2013	ggu	changed VERS_6_1_66
!  05.09.2013	ggu	handle variable choice better
!  12.09.2013	ggu	changed VERS_6_1_67
!  28.01.2014	ggu	changed VERS_6_1_71
!  21.10.2014	ggu	changed VERS_7_0_3
!  05.12.2014	ggu	changed VERS_7_0_8
!  19.01.2015	ggu	changed VERS_7_1_2
!  19.01.2015	ggu	changed VERS_7_1_3
!  05.05.2015	ggu	changed VERS_7_1_10
!  17.07.2015	ggu	changed VERS_7_1_80
!  20.07.2015	ggu	changed VERS_7_1_81
!  18.12.2015	ggu	changed VERS_7_3_17
!  25.05.2016	ggu	changed VERS_7_5_10
!  18.12.2018	ggu	changed VERS_7_5_52
!  21.05.2019	ggu	changed VERS_7_5_62
! 
! **********************************************************
! **********************************************************
! **********************************************************
! **********************************************************

!==================================================================
        module mod_plot
!==================================================================

!  initializes actual level
! 
!  -1	bottom
!   0	integrated
!  >0	level

	implicit none

	integer, save :: level3 = 0
	integer, save :: ivsect = 0
	integer, save :: ivar3 = 0
	real, save :: flagco = -999.

	integer, save :: icall = 0

!==================================================================
        end module mod_plot
!==================================================================

	subroutine init_plot

	use mod_plot
	use para

	implicit none

	integer ivar
	character*80 name,vsect
	real getpar

	if( icall > 0 ) return

	if( para_has_name('level') ) then
	  level3 = nint(getpar('level'))
	end if

	if( para_has_name('vsect') ) then
	  call getfnm('vsect',vsect)
	  ivsect = 0
	  if( vsect .ne. ' ' ) ivsect = 1
	end if

	ivar = 0
	name = ' '
	if( para_has_name('ivar') ) then
	  ivar = nint(getpar('ivar'))
	end if
	if( para_has_name('varnam') ) then
	  call getfnm('varnam',name)
	end if
	if( ivar > 0 .and. name /= ' ' ) then
	  write(6,*) 'You cannot give both ivar and varnam'
	  stop 'error stop init_plot: non compatible variables'
	end if
	if( name .ne. ' ' ) call string2ivar(name,ivar)
	ivar3 = ivar

	icall = 1

	end

! **********************************************************

	subroutine setlev( level )

!  set actual level

	use mod_plot

	implicit none

	integer level

	level3 = level

	end

! **********************************************************

	function getlev()

!  get actual level

	use mod_plot

	implicit none

	integer getlev

	getlev = level3

	end

! **********************************************************
! **********************************************************
! **********************************************************

        function getisec()

!  is it a vertical section

	use mod_plot

        implicit none

        integer getisec

	getisec = ivsect

        end

! **********************************************************
! **********************************************************
! **********************************************************

	subroutine setvar(ivar)

!  set actual variable

	use mod_plot

	implicit none

	integer ivar

	ivar3 = ivar

	end

! **********************************************************

	function getvar()

!  get actual variable

	use mod_plot

	implicit none

	integer getvar

	getvar = ivar3

	end

! **********************************************************

	function okvar(ivar)

!  shall we plot this variable ?

	use mod_plot

	implicit none

	logical okvar
        integer ivar

	okvar = ivar3 .eq. ivar .or. ivar3 .le. 0

	end

! **********************************************************

       subroutine checkvar(ivar)

!  checks what variable has to be plotted
!  returns in ivar the variable to be plotted

       implicit none

       integer ivar

       integer ivar3
       integer getvar

       if( ivar .gt. 0 ) then  !ivar given - must be equal to ivar3
         ivar3 = getvar()
         if( ivar3 .gt. 0 .and. ivar3 .ne. ivar ) goto 99
         call setvar(ivar)
       else
         ivar = getvar()
       end if

       if( ivar .le. 0 ) then
         write(6,*) 'Do not know what to plot: ivar = ',ivar
         stop 'error stop checkvar: no ivar given'
       end if

       return
   99  continue
       write(6,*) 'ivar3 = ',ivar3
       write(6,*) 'ivar  = ',ivar
       stop 'error stop checkvar: different values of ivar3 and ivar'
       end

! **********************************************************
! **********************************************************
! **********************************************************

	subroutine extnlev(level,nlvddi,nkn,p3,p2)

!  extract level from 3d array (nodes)

	use levels

	implicit none

	integer level		!level to extract
	integer nlvddi		!vertical dimension of p3
	integer nkn		!number of nodes
	real p3(nlvddi,nkn)	!3d array
	real p2(nkn)		!2d array filled on return

	call extlev(level,nlvddi,nkn,ilhkv,p3,p2)

	end

! **********************************************************

	subroutine extelev(level,nlvddi,nel,p3,p2)

!  extract level from 3d array (elements)

	use levels

	implicit none

	integer level		!level to extract
	integer nlvddi		!vertical dimension of p3
	integer nel		!number of elements
	real p3(nlvddi,nel)	!3d array
	real p2(nel)		!2d array filled on return

	call extlev(level,nlvddi,nel,ilhv,p3,p2)

	end

! **********************************************************

	subroutine extlev(level,nlvddi,n,ilv,p3,p2)

!  extract level from 3d array

	implicit none

	integer level		!level to extract
	integer nlvddi		!vertical dimension of p3
	integer n		!number values
        integer ilv(n)
	real p3(nlvddi,n)	!3d array
	real p2(n)		!2d array filled on return

	integer i,lmax,lact
        real flag

	if( level .gt. nlvddi ) then
	  write(6,*) 'level, nlvddi : ',level,nlvddi
	  stop 'error stop extlev: level'
	end if

        call get_flag(flag)

        if( level .lt. -1 ) then
	  stop 'error stop extlev: internal error (1)'
	end if

	if( level .eq. 0 ) then
	  call intlev(nlvddi,n,ilv,p3,p2)		!integrate
	else
	  do i=1,n
	    lmax = ilv(i)
	    lact = level
	    if( lact .eq. -1 ) lact = lmax
	    p2(i) = flag
            if( lact .le. ilv(i) ) p2(i) = p3(lact,i)
	  end do
	end if

	end

! **********************************************************

	subroutine intlev(nlvddi,n,ilv,p3,p2)

!  integrate over water column

	implicit none

	integer nlvddi		!vertical dimension of p3
	integer n		!number of nodes
        integer ilv(n)
	real p3(nlvddi,n)	!3d array
	real p2(n)		!2d array filled on return

	integer i,l,lmax
	real value

	lmax = 0

	do i=1,n
	  lmax = ilv(i)
	  if( lmax .eq. 1 ) then	!2d
	    p2(i) = p3(1,i)
	  else				!primitive method of averaging
	    if( lmax .gt. nlvddi ) goto 99
	    if( lmax .le. 0 ) goto 99
	    value = 0.
	    do l=1,lmax
	      value = value + p3(l,i)
	    end do
	    p2(i) = value / lmax
	  end if
	end do

	return
   99	continue
	write(6,*) 'lmax,nlvddi : ',lmax,nlvddi
	stop 'error stop intlev : error in lmax'
	end

! **********************************************************
! **********************************************************
! **********************************************************

	subroutine set_flag(flag)

	use mod_plot

	implicit none

	real flag

	flagco = flag

	end

	subroutine get_flag(flag)

	use mod_plot

	implicit none

	real flag

	flag = flagco

	end

! **********************************************************

