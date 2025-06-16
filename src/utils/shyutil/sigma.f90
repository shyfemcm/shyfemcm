
!--------------------------------------------------------------------------
!
!    Copyright (C) 2011-2017,2019  Georg Umgiesser
!    Copyright (C) 2012  Debora Bellafiore
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

! sigma utilities for output files
!
! revision log :
!
! 07.11.2011	ggu	layer thickness for hybrid coordinates
! 14.11.2011	ggu	new sigma routines copied to this file
! 22.11.2011	ggu	changed VERS_6_1_37
! 02.12.2011	ggu	bug fix in init_sigma_info() for nlv == 1
! 09.12.2011	ggu	changed VERS_6_1_38
! 16.12.2011	ggu	check for non-initialized data structure (blockdata)
! 19.12.2011	ggu	bug fix in init_sigma_info(): call set_sigma_info()
! 24.01.2012	ggu	changed VERS_6_1_41
! 27.01.2012	dbf&ggu	changes to get_layer_thickness()
! 27.01.2012	dbf&ggu	new routine compute_sigma_info()
! 30.03.2012	ggu	changed VERS_6_1_51
! 25.01.2013	ggu	changed VERS_6_1_62
! 03.05.2013	ggu	changed VERS_6_1_63
! 17.05.2013	ggu	layer_thickness for elem and node, general routine
! 17.05.2013	ggu	new routine get_bottom_of_layer()
! 13.06.2013	ggu	changed VERS_6_1_65
! 05.09.2013	ggu	new call interface to get_layer_thickness()
! 12.09.2013	ggu	changed VERS_6_1_67
! 15.05.2014	ggu	changed VERS_6_1_75
! 25.06.2014	ggu	error stop if computed layer thickness is <= 0
! 07.07.2014	ggu	changed VERS_6_1_79
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 15.02.2015	ggu	in get_layer_thickness() handle last layer correctly
! 31.07.2015	ggu	changed VERS_7_1_84
! 18.12.2015	ggu	changed VERS_7_3_17
! 19.02.2016	ggu	changed VERS_7_5_2
! 01.05.2016	ggu	changes in get_layer_thickness(): exit from loop
! 14.05.2016	ggu	substitute blockdata/common with module
! 25.05.2016	ggu	changed VERS_7_5_10
! 07.10.2017	ggu	substitute get_bottom_of_layer with get_depth_of_layer
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 29.04.2022	ggu	in compute_sigma_info() check for nlv<1
! 09.05.2023    lrp     introduce top layer index variable
! 05.06.2023    lrp     introduce z-star
! 09.03.2025	ggu	avoid out of bounds access (hlvaux)
!
! notes :
!
! this file is used also in femplot
!
!	get_sigma_info (supdep.f,suplin.f)
!	get_layer_thickness (supdep.f,suplin.f)
!	init_sigma_info (supout.f)
!
!******************************************************************
!******************************************************************
!******************************************************************

!==================================================================
        module sigma
!==================================================================

	implicit none

	integer, save :: nlv_com    = -1
	integer, save :: nsigma_com = -1
	real, save ::    hsigma_com = 10000.

!==================================================================
        end module sigma
!==================================================================

!******************************************************************

	subroutine check_sigma_initialized

	use sigma

	implicit none

	if( nlv_com .le. 0 ) then
	  write(6,*) 'nlv_com: ',nlv_com
	  stop 'error stop check_sigma_initialized: not initialized'
	end if

	end

!******************************************************************

	subroutine get_sigma_info(nlv,nsigma,hsigma)

	use sigma

	implicit none

	integer nlv
	integer nsigma
	real hsigma

	call check_sigma_initialized

	nlv    = nlv_com
	nsigma = nsigma_com
	hsigma = hsigma_com

	end

!******************************************************************

	subroutine set_sigma_info(nlv,nsigma,hsigma)

	use sigma

	implicit none

	integer nlv
	integer nsigma
	real hsigma

	nlv_com    = nlv
	nsigma_com = nsigma
	hsigma_com = hsigma

	end

!******************************************************************

	subroutine init_sigma_info(nlv,hlv)

	implicit none

	integer nlv
	real hlv(1)

	integer nsigma
	real hsigma

	call compute_sigma_info(nlv,hlv,nsigma,hsigma)
	call set_sigma_info(nlv,nsigma,hsigma)

	end

!******************************************************************
!******************************************************************
!******************************************************************
! next routines can be used without using routines above (common)
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine compute_sigma_info(nlv,hlv,nsigma,hsigma)

	implicit none

	integer nlv		!total number of layers
	real hlv(nlv)		!layer structure
	integer nsigma		!total number of sigma layers (return)
	real hsigma		!closing depth of hybrid layers (return)

	integer l

!---------------------------------------------------------
! scan depth structure
!---------------------------------------------------------

	if( nlv < 1 ) then
	  write(6,*) 'nlv = ',nlv
	  stop 'error stop compute_sigma_info: nlv<1'
	end if

	hsigma = 10000.
        l = 2                           !HACK for nlv == 1
        if( nlv .eq. 1 ) goto 1

	do l=2,nlv
	  if( hlv(l) .gt. hlv(l-1) ) goto 1
	end do

!---------------------------------------------------------
! only sigma layers found
!---------------------------------------------------------

	if( hlv(nlv) .ne. -1 ) then
          write(6,*) 'nlv,hlv(nlv): ',nlv,hlv(nlv)
	  write(6,*) (hlv(l),l=1,nlv)
	  stop 'error stop compute_sigma_info: internal error (1)'
	end if
	nsigma = nlv
	return

!---------------------------------------------------------
! zeta or hybrid levels found
!
! this algorithm cannot handle hybrid levels with only 2 sigma layers
!---------------------------------------------------------

    1	continue
	if( l .eq. 2 ) then	!only zeta levels
	  nsigma = 0
	else			!hybrid levels
	  nsigma = l
	  hsigma = hlv(l)
	end if

!---------------------------------------------------------
! end of routine
!---------------------------------------------------------

	end

!******************************************************************

	subroutine get_layer_thickness(lmax,lmin,nsigma,nadapt,     &
      &				       hsigma,hadapt,z,h,hlv,hdl)

! returns layer thickness - works also for lmax higher than actual layers
!
! works also for lmax higher than actual layers
! in this case the last values for hl are 0

	implicit none

	integer lmax		!index of bottom layer
	integer lmin		!index of surface layer
	integer nsigma		!total number of sigma layers
        integer nadapt          !number of z-surface-adaptive layers
	real hsigma		!closing depth of hybrid layers
	real hadapt		!closing depth of z-surface-adaptive layers
	real z			!water level
	real h			!total depth
	real hlv(lmax)		!layer structure
	real hdl(lmax)		!layer thickness computed (return)

	logical bdebug,berror
	integer ii,l
	real zmed
	real htot,hsig,hzad,htop,hbot,den
	real hlvaux(0:lmax)

	bdebug = .true.
	bdebug = .false.
	berror = .false.

	zmed = z
	htot = h
	hdl = 0.

!---------------------------------------------------------
! compute level structure of sigma levels
!---------------------------------------------------------

	hsig = min(htot,hsigma) + zmed

	hbot = 0.
	do l=1,nsigma
	  htop = hbot
	  hbot = hlv(l)
	  if( l .eq. nsigma ) hbot = -1.
	  hdl(l) = -hsig * (hbot-htop)
	end do

!---------------------------------------------------------
! compute level structure of z-surface-adaptive levels
!---------------------------------------------------------

	hlvaux(0) = 0.
	hlvaux(1:lmax) = hlv

        hzad = hadapt + zmed
	den  = hadapt - hlvaux(lmin-1)		!zstar def
	if (nadapt+lmin-1.eq.lmax) then 
	  hzad = htot + zmed
          den  = htot - hlvaux(lmin-1)	
	end if

        hbot = 0.	
        do l=lmin,lmin+nadapt-1
          htop = hbot
          hbot = hlvaux(l)
	  if( l .eq. lmax ) hbot = htot 
          hdl(l) = hzad * (hbot-htop)/den
        end do

	if( bdebug ) write(6,*) l,hsig,hzad,lmax,nsigma,nadapt
	if( bdebug ) write(6,*) 'hdl: ',hdl

!---------------------------------------------------------
! compute level structure of zeta and/or hybrid levels
!---------------------------------------------------------

	if( lmax .gt. (nsigma+nadapt) ) then	!also zeta coordinates
	  if( lmax .eq. lmin ) then		!just one layer
	    hdl(lmin) = htot + zmed
	  else
	    if( nsigma .ne. 0 ) hbot = hsigma
	    if( nadapt .ne. 0 ) hbot = hadapt
	    if( nsigma .eq. 0 .and.  nadapt .eq. 0 ) hbot = -zmed
	    if( bdebug ) write(6,*) nsigma,lmax,zmed,hbot
	    do l=nsigma+nadapt+lmin,lmax
	      if( hbot == htot ) exit	!no more layers
	      htop = hbot
	      hbot = hlvaux(l)
	      if( bdebug ) write(6,*) l,htop,hbot,htot
	      if( hbot .gt. htot ) hbot = htot	!last layer
	      hdl(l) = hbot - htop
	      if( hdl(l) .le. 0. ) berror = .true.
	    end do
	    if( htot > hbot ) hdl(lmax) = hdl(lmax) + htot - hbot
	  end if
	end if
	if( bdebug ) write(6,*) 'hdl: ',hdl

	!if( berror ) goto 99

!---------------------------------------------------------
! end of routine
!---------------------------------------------------------

	return
   99	continue
	write(6,*) 'error computing layer thickness'
	write(6,*) 'lmax,nsigma: ',lmax,nsigma
	write(6,*) 'hsigma,z,h: ',hsigma,z,h
	write(6,*) 'hlv: '
	write(6,*) (hlv(l),l=1,lmax)
	write(6,*) 'hd: '
	write(6,*) (hdl(l),l=1,lmax)
	stop 'error stop get_layer_thickness: 0 thickness'
	end

!******************************************************************

	subroutine get_depth_of_layer(bcenter,lmax,z,hl,hz)

! computes depth of layer (center if bcenter == .true., else bottom)

	implicit none

	logical bcenter	!compute depth at center of layer (else bottom)
	integer lmax	!total number of layers
	real z		!water level
	real hl(lmax)	!layer thickness (from get_layer_thickness)
	real hz(lmax)	!depth at layer depth/center (return)

	integer l
	real htop,hbot

	htop = -z

	do l=1,lmax
	  hbot = htop + hl(l)
	  hz(l) = hbot
	  if( bcenter ) hz(l) = 0.5*(htop+hbot)
	  htop = hbot
	end do

	end

!******************************************************************

	subroutine adjust_layer_index(nel,nlv,hev,hlv,ilhv)

	integer nel,nlv
	real hev(nel)
	real hlv(nlv)
	integer ilhv(nel)

	integer ie,l
	integer nlvaux,nsigma
	real hsigma,z,h
	real hdl(nlv)

	z = 0.
	call get_sigma_info(nlvaux,nsigma,hsigma)

	do ie=1,nel
	  h = hev(ie)
	  call get_layer_thickness(nlv,1,nsigma,0,hsigma,0.,z,h,hlv,hdl)
	  do l=1,nlv
	    if( hdl(l) == 0. ) exit
	  end do
	  ilhv(ie) = l - 1
	end do

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine compute_iztype(iztype)

! computes type of vertical coordinates

	implicit none

	integer iztype		!type of vertical coordinates (return)

	integer nlv,nsigma
	real hsigma

	call get_sigma_info(nlv,nsigma,hsigma)

	if( nsigma .eq. 0 ) then		! z-coordinates
	  iztype = 1
	else if( nlv .eq. nsigma ) then		! sigma-coordinates
	  iztype = 2
	else					! hybrid-coordinates
	  iztype = 3
	end if

	end
	
!******************************************************************

	subroutine sigma_test

	integer, parameter :: ndim = 5
	integer lmax,nsigma
	real hsigma,z,h
	real hlv(ndim)
	real hdl(ndim)

	hlv = (/2,4,6,8,10/)

	nsigma = 0
	hsigma = 10000.
	z = 0.
	h = 4.2

	lmax = 5
	call get_layer_thickness(lmax,1,nsigma,0,hsigma,0.,z,h,hlv,hdl)
	write(6,*) lmax,hdl

	lmax = 2
	call get_layer_thickness(lmax,1,nsigma,0,hsigma,0.,z,h,hlv,hdl)
	write(6,*) lmax,hdl

	end

!******************************************************************
!	program sigma_main
!	call sigma_test
!	end
!******************************************************************

