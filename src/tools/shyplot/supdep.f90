
!--------------------------------------------------------------------------
!
!    Copyright (C) 2000,2008-2011,2013-2019  Georg Umgiesser
!    Copyright (C) 2012  Debora Bellafiore
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

!  routines for averaging depth
! 
!  contents :
! 
!  subroutine mkht(hetv,href)		makes hetv (elem depth of actual layer)
!  subroutine mkht3(nlvddi,het3v,href)	makes het3v (3D depth structure)
!  function hlthick(l,lmax,hl)		layer thickness
! 
!  revision log :
! 
!  26.05.2000	ggu	routines written from scratch
!  17.09.2008	ggu	routine mkht changed for layer = -1
!  13.10.2009	ggu	new routine mkht3
!  13.10.2009	ggu	new routine mkht3
!  23.03.2010	ggu	changed v6.1.1
!  17.12.2010	ggu	substituted hv with hkv, new routine hlthick()
!  30.03.2011	ggu	new routine mkareafvl()
!  14.04.2011	ggu	changed VERS_6_1_22
!  14.11.2011	ggu	use get_layer_thickness() for layer structure
!  22.11.2011	ggu	changed VERS_6_1_37
!  23.02.2012	ccf	bug fix in mkht (get_layer_thicknes, get_sigma_info)
!  16.03.2012	dbf	bug fix in mkht3 (get_layer_thicknes, get_sigma_info)
!  10.06.2013	ggu	bug fix in mkht,mkht3 (get_layer_thicknes_e)
!  05.09.2013	ggu	adapt to new get_layer_thickness()
!  12.09.2013	ggu	changed VERS_6_1_67
!  23.12.2014	ggu	changed VERS_7_0_11
!  19.01.2015	ggu	changed VERS_7_1_2
!  19.01.2015	ggu	changed VERS_7_1_3
!  05.05.2015	ggu	changed VERS_7_1_10
!  10.07.2015	ggu	changed VERS_7_1_50
!  17.07.2015	ggu	changed VERS_7_1_80
!  20.07.2015	ggu	changed VERS_7_1_81
!  25.05.2016	ggu	changed VERS_7_5_10
!  27.05.2016	ggu	mkhkv,mkhev deleted
!  12.01.2017	ggu	changed VERS_7_5_21
!  25.10.2018	ggu	changed VERS_7_5_51
!  18.12.2018	ggu	changed VERS_7_5_52
!  13.03.2019	ggu	changed VERS_7_5_61
!  21.05.2019	ggu	changed VERS_7_5_62
!  26.07.2023   lrp     introduce zstar in shyplot
!  28.09.2023   lrp     bug fix for zstar (forgotten parameter)
! 
! ******************************************************************

	subroutine mkht(hetv,href)

!  makes hetv (elementwise depth of actual layer)
! 
!  uses level to decide what to do

	use mod_depth
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real hetv(nel)
	real href

	logical bdebug
	integer ie,ii
	integer level,lmax,lmin
	real zeta,zmin
	integer nsigma,nadapt,nlvaux
	real hsigma,hadapt
        real hl(nlv)

	integer getlev
	real hlthick

! -------------------------------------------------------------------
!  initialization
! -------------------------------------------------------------------

        bdebug = .true.
        bdebug = .false.

	level = getlev()
        call get_sigma_info(nlvaux,nsigma,hsigma)

! -------------------------------------------------------------------
!  handle different kind of levels
! -------------------------------------------------------------------

	do ie=1,nel
	  lmax = ilhv(ie)
	  lmin = 1
	  call compute_levels_on_element(ie,zenv,zeta)
	  zmin = minval(zenv(:,ie))   !min: z-adapt coords works with zmin
	  call compute_zadapt_info(zmin,hlv,nsigma,lmax,lmin, &
     & 				   nadapt,hadapt)
	  call get_layer_thickness(lmax,lmin,nsigma,nadapt, &
     &				   hsigma,hadapt,zeta,hev(ie),hlv,hl)
	  hetv(ie) = hlthick(level,lmax,hl)
	end do

! -------------------------------------------------------------------
!  debug output
! -------------------------------------------------------------------

	if( bdebug ) then
	  ie = 1644
	  write(6,*) 'debugging mkht...'
	  write(6,*) ilhv(ie),hetv(ie),hev(ie)
	  do ii=1,ilhv(ie)
	    write(6,*) hlv(ii)
	  end do
	end if

! -------------------------------------------------------------------
!  end of routine
! -------------------------------------------------------------------

	end

! ******************************************************************

	subroutine mkht3(nlvddi,het3v,href)

!  makes het3v (3D depth structure)

	use mod_depth
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

!  arguments
	integer nlvddi
	real het3v(nlvddi,nel)
	real href
!  local
	logical bdebug
	integer ie,ii
	integer l,lmax,lmin
	integer nlvaux,nsigma,nadapt
	real hsigma,hadapt
	real zeta,zmin

! -------------------------------------------------------------------
!  initialization
! -------------------------------------------------------------------

        bdebug = .true.
        bdebug = .false.

        call get_sigma_info(nlvaux,nsigma,hsigma)

! -------------------------------------------------------------------
!  compute layer thickness
! -------------------------------------------------------------------

	do ie=1,nel
	  lmax = ilhv(ie)
	  lmin = 1
	  call compute_levels_on_element(ie,zenv,zeta)
          zmin = minval(zenv(:,ie))   !min: z-adapt coords works with zmin
          call compute_zadapt_info(zmin,hlv,nsigma,lmax,lmin, &
     &                             nadapt,hadapt)
	  call get_layer_thickness(lmax,lmin,nsigma,nadapt, &
     &			hsigma,hadapt,zeta,hev(ie),hlv,het3v(1,ie))
!	  call get_layer_thickness_e(ie,lmax,bzeta,nsigma,hsigma
!     +				,het3v(1,ie))
	end do

! -------------------------------------------------------------------
!  end of routine
! -------------------------------------------------------------------

	end

! ******************************************************************

	function hlthick(l,lmax,hl)

!  computes thickness of layer l
! 
!  works also for sigma layers

	implicit none

	real hlthick		! layer thickness (return)
	integer l		! layer to compute thickness
	integer lmax		! maximum layers
	real hl(lmax)		! level thickness

	integer ll
	real hm

	if( l .eq. 0 ) then			! total layer
	  hm = 0.
	  do ll=1,lmax
	    hm = hm + hl(ll)
	  end do
	  hlthick = hm
	else if( l .ge. 1 .and. l .le. lmax ) then
	  hlthick = hl(l)
	else if( l .eq. -1 ) then		! bottom layer
	  hlthick = hl(lmax)
	else
	  hlthick = 0.
	end if

	end

! ******************************************************************

	subroutine mkareafvl

!  makes area of finite volume

	use mod_hydro_plot
	use evgeom
	use basin

	implicit none

	integer ie,ii,k
	real afvl

	do k=1,nkn
	  arfvlv(k) = 0.
	end do

	do ie=1,nel
	  afvl = 4. * ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    arfvlv(k) = arfvlv(k) + afvl
	  end do
	end do

	end

! ******************************************************************

