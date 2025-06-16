
!--------------------------------------------------------------------------
!
!    Copyright (C) 2013,2015-2019  Georg Umgiesser
!    Copyright (C) 2016  Christian Ferrarin
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

! utility routines for flux computations
!
! contents :
!
! subroutine flxscs(n,kflux,iflux,az,fluxes)	flux through sections
! subroutine flxsec(n,kflux,iflux,az,fluxes)	flux through section
!
! subroutine flxini				initializes flux routines
! subroutine flx_init(kfluxm,kflux,nsect,iflux)	sets up array iflux
! subroutine flxinf(m,kflux,iflux)		sets up one info structure
! function igtnsc(k1,k2)			gets number of internal section
!
! revision log :
!
! 09.05.2013	ggu	separated from subflxa.f
! 14.05.2013	ggu	deleted error check between 2d and 3d computation
! 13.06.2013	ggu	changed VERS_6_1_65
! 19.01.2015	ggu	changed VERS_7_1_2
! 19.01.2015	ggu	changed VERS_7_1_3
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 16.12.2015	ggu	changed VERS_7_3_16
! 18.12.2015	ggu	changed VERS_7_3_17
! 15.04.2016	ggu	changed VERS_7_5_8
! 26.10.2016	ccf	bug fix in flxsec
! 12.01.2017	ggu	changed VERS_7_5_21
! 30.03.2017	ggu	changed accumulator to time step dt, not number of calls
! 04.02.2018	ggu	new routines with accumulator in double
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 27.05.2022	ggu	prepared for mpi use, fluxes now double
! 28.10.2023	ggu	bug fix for node at boundary with ibtyp==3
!
! notes :
!
! These routines can also be used internally to compute the flux
! over various sections. The following calling sequence must be respected:
!
! call flx_init(kfluxm,kflux,nsect,iflux)		initializes iflux
!
! call flxscs(kfluxm,kflux,iflux,az,fluxes) computes fluxes 
!
! Initialization can be done anytime.
!
!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************


!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************


!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine fluxes_init(nlvddi,nsect,nlayers,tr,masst)

! initializes nr and masst

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	real tr
	real masst(0:nlvddi,3,nsect)

	integer i,l,lmax

        tr = 0.
	masst = 0.

	end

!******************************************************************

	subroutine fluxes_accum(nlvddi,nsect,nlayers,dt,tr,masst,fluxes)

! accumulates fluxes into masst

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	real dt
	real tr
	real masst(0:nlvddi,3,nsect)
	real fluxes(0:nlvddi,3,nsect)

        tr = tr + dt
	masst = masst + fluxes * dt

	end

!******************************************************************

	subroutine fluxes_aver(nlvddi,nsect,nlayers,tr,masst,fluxes)

! averages masst and puts result into fluxes

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	real tr
	real masst(0:nlvddi,3,nsect)
	real fluxes(0:nlvddi,3,nsect)

	if( tr == 0. ) return
        fluxes = masst / tr

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine fluxes_init_d(nlvddi,nsect,nlayers,tr,masst)

! initializes nr and masst

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	double precision tr
	double precision masst(0:nlvddi,3,nsect)

	integer i,l,lmax

        tr = 0.
	masst = 0.

	end

!******************************************************************

	subroutine fluxes_accum_d(nlvddi,nsect,nlayers,dt,tr,masst,fluxes)

! accumulates fluxes into masst

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	real dt
	double precision tr
	double precision masst(0:nlvddi,3,nsect)
	double precision fluxes(0:nlvddi,3,nsect)

        tr = tr + dt
	masst = masst + fluxes * dt

	end

!******************************************************************

	subroutine fluxes_aver_d(nlvddi,nsect,nlayers,tr,masst,fluxes)

! averages masst and puts result into fluxes

	implicit none

	integer nlvddi,nsect
	integer nlayers(nsect)
	double precision tr
	double precision masst(0:nlvddi,3,nsect)
	double precision fluxes(0:nlvddi,3,nsect)

	if( tr == 0. ) return
        fluxes = masst / tr

	end

!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine flxscs(kfluxm,kflux,iflux,az,fluxes,is,scalar)

! computes flux through all sections and returns them in fluxes
!
! flux are divided into total, positive and negative

	use levels, only : nlvdi,nlv

	implicit none

	integer kfluxm
	integer kflux(kfluxm)
	integer iflux(3,kfluxm)
	real az
	double precision fluxes(0:nlvdi,3,*)	!computed fluxes (return)
	integer is			!type of scalar (0=mass)
	real scalar(nlvdi,*)

	integer nnode,ifirst,ilast,ntotal
	integer ns
	logical nextline

	nnode = 0
	ns = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
	  !write(66,*) 'section ',ns,ntotal,is
	  call flxsec(ntotal,kflux(ifirst),iflux(1,ifirst),az &
     &				,fluxes(0,1,ns),is,scalar)
	end do

	end

!******************************************************************

	subroutine flxsec(n,kflux,iflux,az,fluxes,is,scalar)

! computes flux through one section and returns it in fluxes

	use levels, only : nlvdi,nlv

	implicit none

	integer n
	integer kflux(n)
	integer iflux(3,n)
	real az
	double precision fluxes(0:nlvdi,3)	!computed fluxes (return)
	integer is				!type of scalar (0=mass)
	real scalar(nlvdi,*)

	integer i,k,l,lkmax
	integer istype,iafter,ibefor
	real port,port2d,sport
	real flux(nlvdi)

	fluxes = 0.

	do i=1,n
		k = kflux(i)
		if( k <= 0 ) cycle
		istype = iflux(1,i)
		ibefor = iflux(2,i)
		iafter = iflux(3,i)
		if( istype == 0 ) cycle		!not in domain

		!port2d = 0.
		!call flx2d(k,ibefor,iafter,istype,az,port2d)

		flux = 0.
		call flx3d(k,ibefor,iafter,istype,az,lkmax,flux)

		do l=1,lkmax
		  port = flux(l)
		  sport = port
		  if( is .gt. 0 ) sport = port * scalar(l,k)	!not mass

		  fluxes(l,1) = fluxes(l,1) + sport
		  if( port .gt. 0. ) then
		    fluxes(l,2) = fluxes(l,2) + sport
		  else
		    fluxes(l,3) = fluxes(l,3) - sport
		  end if
		end do
	end do

! compute vertical integrated fluxes

        fluxes(0,1) = sum(fluxes(1:nlvdi,1))
        fluxes(0,2) = sum(fluxes(1:nlvdi,2))
        fluxes(0,3) = sum(fluxes(1:nlvdi,3))

	end
	  
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine flx_init(kfluxm,kflux,nsect,iflux)

! does basic checks and sets up array iflux

	implicit none

        integer kfluxm		!total number of nodes in kflux
        integer kflux(kfluxm)	!nodes in sections
	integer nsect		!number of section (return)
	integer iflux(3,kfluxm)	!internal array for flux computation (return)

	integer ifirst,ilast,nnode,ntotal,ns

	integer klineck
	logical nextline

!----------------------------------------------------------
! check nodes for compatibility
!----------------------------------------------------------

	nsect = klineck(kfluxm,kflux)

	if( nsect .lt. 0 ) then
	  write(6,*) 'errors setting up fluxes ',nsect
	  write(6,*) 'error happens in klineck '
	  stop 'error stop : flx_init'
	end if

!----------------------------------------------------------
! now set info structure for sections
!----------------------------------------------------------

	ns = 0
	nnode = 0
	iflux = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
!	  write(6,*) kfluxm,nnode,ifirst,ilast,ntotal
	  call flxinf(ns,ntotal,kflux(ifirst),iflux(1,ifirst))
	  call flxtst(ns,ntotal,kflux(ifirst),iflux(1,ifirst))
	end do

	!stop 'flxtst output'

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	end

!******************************************************************

	subroutine flxtst(ns,n,kflux,iflux)

! tests structure iflux(3,1) for one section

	use shympi

	implicit none

	integer ns
	integer n
	integer kflux(n)
	integer iflux(3,n)

	integer i,k,nmax,nt,id
	character*80 string
	character*80 format
	integer ids(n)

	nmax = min(80,n)
	nmax = n
	nt = 0

	do i=1,nmax
	  k = kflux(i)
	  if( k > 0 ) then
	    nt = nt + 1
	    id = id_node(k)
	    !string(i:i) = '1'
	  else
	    id = -1
	    !string(i:i) = '0'
	  end if
	  ids(i) = id
	end do

	write(format,'(a,i4,a)') '(a,',n+3,'i4)'
	!write(6,*) 'format: ',trim(format)
	call shympi_syncronize
	if( nt > 0 ) then
	  write(6,format) 'flxtst: ',my_id,ns,n,ids
	end if
	call shympi_syncronize

	end

!******************************************************************

	subroutine flxinf(ns,n,kflux,iflux)

! sets up info structure iflux(3,1) for one section

	implicit none

	integer ns
	integer n
	integer kflux(n)
	integer iflux(3,n)

	integer i,k
	integer nt
	integer ktype
	integer kafter,kbefor
	integer ngood
	logical berror
	character*15 what

	integer igtnsc,flxtype

	ngood = 0

	do i=1,n
	  k = kflux(i)
	  if( k == 0 ) stop 'error stop flxinf: (1)'
	  if( k < 0 ) cycle
	  ngood = ngood + 1
	  if( n == 1 ) then	!one node section ... must be ibtyp == 3
	    ktype = 1		!fake inner node
	  else
	    ktype = flxtype(k)
	  end if

	  iflux(1,i) = ktype

	  kbefor = 0
	  if( i .ne. 1 ) kbefor = kflux(i-1)
	  kafter = 0
	  if( i .ne. n ) kafter = kflux(i+1)

	  iflux(2,i) = igtnsc(k,kbefor)
	  iflux(3,i) = igtnsc(k,kafter)
	end do

	nt = n
	berror = .false.
        if( ngood == nt ) then
	  what = 'sect full'
        else if( ngood == 0 ) then
	  what = 'sect empty'
        else
	  what = 'sect partial'
        end if
        !write(6,*) what,ns,nt,ngood

	end

!******************************************************************

	function igtnsc(k1,k2)

! gets number of internal section in link index of k1

	use mod_geom

	implicit none

	integer igtnsc
	integer k1,k2
	integer elems(maxlnk)

	integer k,i,n,ie

	integer knext,kbhnd

	call get_elems_around(k1,maxlnk,n,elems)

!	deal with no section given

	igtnsc = 0
	if( k2 .le. 0 ) return

!	section in the middle

	do i=1,n
	  igtnsc = i
	  ie = elems(i)
	  k = knext(k1,ie)
	  if( k .eq. k2 ) return
	end do

!	in case we are on the boundary

	i = n
	igtnsc = igtnsc + 1
	ie = elems(i)
	k = kbhnd(k1,ie)
	if( k .eq. k2 ) return

!	no node found

	write(6,*) k1,k2
	write(6,*) k1,n
	write(6,*) (elems(i),i=1,n)
	call get_elems_around(k2,maxlnk,n,elems)
	write(6,*) k2,n
	write(6,*) (elems(i),i=1,n)
	stop 'error stop igtnsc: internal error (2)'
	end
	      
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************

	subroutine get_nlayers(kfluxm,kflux,nsect,nlayers,nlmax)

! computes maximum numer of layers for sections

	use levels

	implicit none

	integer kfluxm
	integer kflux(kfluxm)
	integer nsect
	integer nlayers(nsect)	!total number of layers for sections (return)
	integer nlmax		!maximum layers for all sections (return)

	integer ns
	integer nnode,ifirst,ilast
	integer i,k,l,lmax

	logical nextline

	ns = 0
	nnode = 0
	nlmax = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  !write(6,*) 'get_nlayers: ',ns,ifirst,ilast,ilast-ifirst+1
	  lmax = 0
	  do i=ifirst,ilast
	    k = kflux(i)
	    if( k <= 0 ) cycle
	    l = ilhkv(k)
	    lmax = max(lmax,l)
	  end do
	  nlayers(ns) = lmax
	  nlmax = max(nlmax,lmax)
	end do

	if( ns /= nsect ) then
	  write(6,*) 'ns,nsect: ',ns,nsect
	  stop 'error stop get_nlayers: ns /= nsect'
	end if

	end

!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************

