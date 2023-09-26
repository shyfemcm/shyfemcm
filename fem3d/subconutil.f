
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-1999,2001-2012,2014-2019  Georg Umgiesser
!    Copyright (C) 2008  Christian Ferrarin
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

! routines for concentration (utilities) (old subcon1.f)
!
! contents :
!
! subroutine conini(nlvddi,c,cref,cstrat)		sets initial conditions
!
! subroutine conmima(nlvddi,c,cmin,cmax)                 computes min/max
! subroutine conmimas(nlvddi,c,cmin,cmax)                computes scalar min/max
!
! revision log :
!
! 19.08.1998	ggu	call to conzfi changed
! 20.08.1998	ggu	makew removed (routine used is sp256w)
! 24.08.1998	ggu	levdbg used for debug
! 26.08.1998	ggu	subroutine convol, tstvol transferred to newchk
! 26.08.1998	ggu	all subroutines re-written more generally
! 26.01.1999	ggu	can be used also with 2D routines
! 16.11.2001	ggu	subroutine conmima and diffstab
! 05.12.2001	ggu	new routines diffstab,diffstab1,difflimit
! 11.10.2002	ggu	commented diffset
! 09.09.2003	ggu	new routine con3bnd
! 13.03.2004	ggu	new routines set_c_bound, distribute_vertically
! 13.03.2004	ggu	exec routine con3bnd() only for level BC (LEVELBC)
! 14.03.2004	ggu	new routines open_b_flux
! 05.01.2005	ggu	routine to write 2d nos file into subnosa.f
! 07.01.2005	ggu	routine diffwrite deleted
! 14.01.2005	ggu	new file for diffusion routines (copied to subdif.f)
! 23.03.2006	ggu	changed time step to real
! 31.05.2007	ggu	reset BC of flux type to old way (DEBHELP)
! 07.04.2008	ggu	deleted set_c_bound
! 08.04.2008	ggu	cleaned, deleted distribute_vertically, open_b_flux
! 09.10.2008	ggu&ccf	call to confop changed -> nlv
! 23.03.2010	ggu	changed v6.1.1
! 14.07.2011	ggu	changed VERS_6_1_27
! 01.06.2012	ggu	changed VERS_6_1_53
! 20.01.2014	ggu	new writing format for nos files in confop, confil
! 28.01.2014	ggu	changed VERS_6_1_71
! 26.11.2014	ggu	changed VERS_7_0_7
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 26.02.2015	ggu	changed VERS_7_1_5
! 05.05.2015	ggu	changed VERS_7_1_10
! 05.06.2015	ggu	changed VERS_7_1_12
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 29.09.2015	ggu	changed VERS_7_2_5
! 18.12.2015	ggu	changed VERS_7_3_17
! 28.04.2016	ggu	changed VERS_7_5_9
! 07.06.2016	ggu	changed VERS_7_5_12
! 03.11.2017	ggu	new routines to write shy files scalar_output_*()
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.04.2018	ggu	changed VERS_7_5_43
! 06.07.2018	ggu	changed VERS_7_5_48
! 16.02.2019	ggu	changed VERS_7_5_60
! 20.03.2022	ggu	started discommissioning file
! 21.03.2022	ggu	only some subroutines save to this new file
! 22.03.2022	ggu	cmed routines transfered from subres.f to here
! 05.04.2022	ggu	initialize only existing layers
! 27.10.2022	ggu	conmima also working for 2d arrays
!
!*****************************************************************

	subroutine conini(nlvddi,c,cref,cstrat,hdko)

! sets initial conditions (with stratification)

	use basin, only : nkn,nel,ngr,mbw
	use levels

	implicit none

	integer nlvddi		!vertical dimension of c
	real c(nlvddi,nkn)	!variable to initialize
	real cref		!reference value
	real cstrat		!stratification [conc/km]
	real hdko(nlvddi,nkn)	!layer thickness
! local
	integer k,l,lmax
	real depth,hlayer

	c = 0.

	do k=1,nkn
	  depth=0.
	  lmax = ilhkv(k)
	  do l=1,lmax
	    hlayer = 0.5 * hdko(l,k)
	    depth = depth + hlayer
	    c(l,k) = cref + cstrat*depth/1000.
	    depth = depth + hlayer
	  end do
	end do

	end

!*************************************************************
!*************************************************************
!*************************************************************

        subroutine conmima(nlvddi,c,cmin,cmax)

! computes min/max for scalar field

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

! arguments
	integer nlvddi		!vertical dimension of c
	real c(nlvddi,nkn)	!tracer (conz,salt,temp,...)
        real cmin,cmax
! local
	integer k,l,lmax
	real cc
        logical debug
        integer kcmin,lcmin,kcmax,lcmax

        debug = .false.
        cmin = c(1,1)
        cmax = c(1,1)

	do k=1,nkn
	  lmax=min(nlvddi,ilhkv(k))	!works also for 2d arrays (nlvddi==1)
	  do l=1,lmax
	    cc = c(l,k)
            if( debug ) then
              if( cc .lt. cmin ) then
                    kcmin = k
                    lcmin = l
              end if
              if( cc .gt. cmax ) then
                    kcmax = k
                    lcmax = l
              end if
            end if
            cmin = min(cmin,cc)
            cmax = max(cmax,cc)
	  end do
	end do

        if( debug ) then
          write(6,*) 'conmima min: ',kcmin,lcmin,cmin
          write(6,*) 'conmima max: ',kcmax,lcmax,cmax
        end if

        end

!********************************************************************
!********************************************************************
!********************************************************************

	subroutine cmed_accum_2d(cvec,cmed,cmin,cmax)

! accumulates scalar values

	use basin, only : nkn,nel,ngr,mbw

	implicit none

! parameter

	real cvec(nkn)			!array with concentration
	double precision cmed(nkn)	!average
	real cmin(nkn)			!minimum
	real cmax(nkn)			!maximum

! local
	logical bdebug
	integer nout,id
	integer nvar,nr
	integer idtc,itmc,itc,it
	integer i,k
	real c
	double precision rr

!-------------------------------------------------------------
! accumulate results
!-------------------------------------------------------------

	do k=1,nkn
	    c = cvec(k)
	    cmed(k) = cmed(k) + c
	    if( c .lt. cmin(k) ) cmin(k) = c
	    if( c .gt. cmax(k) ) cmax(k) = c
	end do

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!********************************************************************

	subroutine cmed_reset_2d(nr,cmed,cmin,cmax)

! resets scalar values

	use basin

	implicit none

	integer nr
	double precision cmed(nkn)	!average
	real cmin(nkn)			!minimum
	real cmax(nkn)			!maximum

	real, parameter :: high = 1.e+30

	nr = 0
	cmed = 0.
	cmin = high
	cmax = -high

	end

!********************************************************************

	subroutine cmed_accum(cvec,cmed,cmin,cmax)

! accumulates scalar values

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

! parameter

	real cvec(nlvdi,nkn)			!array with concentration
	double precision cmed(nlvdi,nkn)	!average
	real cmin(nlvdi,nkn)			!minimum
	real cmax(nlvdi,nkn)			!maximum

! local
	logical bdebug
	integer nout,id
	integer nvar,nr
	integer idtc,itmc,itc,it
	integer i,k,l,lmax
	real c
	double precision rr

	bdebug = .true.
	bdebug = .false.

	if( bdebug ) write(6,*) it,nout,id,nvar,nr,idtc,itmc,itc,nlv

!-------------------------------------------------------------
! accumulate results
!-------------------------------------------------------------

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    c = cvec(l,k)
	    cmed(l,k) = cmed(l,k) + c
	    if( c .lt. cmin(l,k) ) cmin(l,k) = c
	    if( c .gt. cmax(l,k) ) cmax(l,k) = c
	  end do
	end do

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!********************************************************************

	subroutine cmed_reset(nr,cmed,cmin,cmax)

! resets scalar values

	use levels, only : nlvdi,nlv
	use basin

	implicit none

	integer nr
	double precision cmed(nlvdi,nkn)	!average
	real cmin(nlvdi,nkn)			!minimum
	real cmax(nlvdi,nkn)			!maximum

	real, parameter :: high = 1.e+30

	nr = 0
	cmed = 0.
	cmin = high
	cmax = -high

	end

!********************************************************************

	subroutine ts_shell

	use mod_ts
	use levels, only : nlvdi,nlv
	use basin

	implicit none

! local
	integer idtc,itmc,itsmed
	integer id,nvar,idc,nr
	double precision dtime
	double precision rr
! function
	real getpar
	logical has_output_d,is_over_output_d,next_output_d
! save
	double precision, save, allocatable :: tacu(:,:)
	double precision, save, allocatable :: sacu(:,:)
	real, save, allocatable :: tmin(:,:)
	real, save, allocatable :: tmax(:,:)
	real, save, allocatable :: smin(:,:)
	real, save, allocatable :: smax(:,:)
	real, save, allocatable :: raux(:,:)

	integer, save :: icall = 0
	double precision, save :: da_out(4) = 0.

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then

	  itsmed=nint(getpar('itsmed'))
	  if( itsmed .le. 0 ) then
	    icall = -1
	    return
	  end if

          call init_output_d('itmcon','idtcon',da_out)
	  call increase_output_d(da_out)
          if( has_output_d(da_out) ) then
            nvar = 2*3
            call shyfem_init_scalar_file('tsmed',nvar,.false.,id)
            da_out(4) = id
          end if

	  allocate(tacu(nlvdi,nkn),tmin(nlvdi,nkn),tmax(nlvdi,nkn))
	  allocate(sacu(nlvdi,nkn),smin(nlvdi,nkn),smax(nlvdi,nkn))
	  allocate(raux(nlvdi,nkn))

	  call cmed_reset(nr,tacu,tmin,tmax)
	  call cmed_reset(nr,sacu,smin,smax)

	  icall = 1
	end if

	if( .not. is_over_output_d(da_out) ) return

	nr = nr + 1
	call cmed_accum(tempv,tacu,tmin,tmax)
	call cmed_accum(saltv,sacu,smin,smax)

        if( .not. next_output_d(da_out) ) return

        id = nint(da_out(4))
	call get_act_dtime(dtime)
	rr=1./nr

	idc = 160
	raux = tacu * rr
	call shy_write_scalar_record(id,dtime,idc+1,nlvdi,raux)
	call shy_write_scalar_record(id,dtime,idc+2,nlvdi,tmin)
	call shy_write_scalar_record(id,dtime,idc+3,nlvdi,tmax)

	idc = 170
	raux = sacu * rr
	call shy_write_scalar_record(id,dtime,idc+1,nlvdi,raux)
	call shy_write_scalar_record(id,dtime,idc+2,nlvdi,smin)
	call shy_write_scalar_record(id,dtime,idc+3,nlvdi,smax)

	call cmed_reset(nr,tacu,tmin,tmax)
	call cmed_reset(nr,sacu,smin,smax)

	end

!********************************************************************

