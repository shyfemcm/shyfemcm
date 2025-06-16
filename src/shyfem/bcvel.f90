
!--------------------------------------------------------------------------
!
!    Copyright (C) 2008-2010,2012-2019  Georg Umgiesser
!    Copyright (C) 2008  Debora Bellafiore
!    Copyright (C) 2014  Christian Ferrarin
!    Copyright (C) 2016  Ivan Federico
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

! fix or nudge velocities at open boundary
!
! contents :
!
! subroutine bclfix_ini		initialization of bclfix routines
! subroutine bclfix		fixes velocities on open boundaries
!
! revision log :
!
! 15.09.2008	ggu	written from scratch
! 03.11.2008	ggu&dbf	nudging implemented
! 12.11.2008	ggu	handle sigma coordinates
! 06.12.2008	ggu	read nbfix from STR
! 19.01.2009	ggu	no error stop in initializing when nbfix=0
! 23.03.2009	ggu	tramp from start of simulation
! 23.03.2010	ggu	changed v6.1.1
! 16.12.2010	ggu	bsigma renamed to bosigma
! 24.01.2012	ggu	changed VERS_6_1_41
! 27.01.2012	ggu	changed VERS_6_1_43
! 26.06.2012	ggu	changed VERS_6_1_55
! 25.01.2013	ggu	changed VERS_6_1_62
! 18.06.2014	ggu	changed VERS_6_1_77
! 29.10.2014	ccf	rewritten for 7_0_3, vel file for each boundary
! 05.11.2014	ggu	changed VERS_7_0_5
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 26.02.2015	ggu	changed VERS_7_1_5
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 30.07.2015	ggu	changed VERS_7_1_83
! 23.09.2015	ggu	changed VERS_7_2_4
! 13.07.2016	ivn	bug fix setting up ielfix
! 09.09.2016	ggu	changed VERS_7_5_17
! 31.03.2017	ggu	changed VERS_7_5_24
! 05.12.2017	ggu	changed VERS_7_5_39
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 03.06.2022	ggu	prepared for mpi (not yet ready)
! 09.05.2023    lrp     introduce top layer index variable
! 14.10.2024	ggu	fixed INTEL_BUG
! 13.11.2024	ggu	marked old code with INTEL_BUG_OLD
!
!*****************************************************************

        subroutine bclfix_ini

! initialization of bclfix routines

	use mod_bclfix
	use mod_internal
	use basin
	use shympi

        implicit none 

	real tnudge	!relaxation time for nudging [s]

        integer ie,l,i,ii,k,n,nn,nf,ie_mpi
	integer ibc,nodes
	integer nbc
	integer iflag(nkn)
       
	integer nkbnds,kbnds,nbnds
	integer ieext
	logical bdebug

	bdebug = .true.
	bdebug = .false.

!------------------------------------------------------------------
! initialize arrays
!------------------------------------------------------------------

	iuvfix = 0
	tnudgev = 0.
	ielfix = 0

!------------------------------------------------------------------
! loop over boundaries
!------------------------------------------------------------------

	nbc = nbnds()

	do ibc = 1, nbc

          call get_bnd_par(ibc,'tnudge',tnudge)

	  if (tnudge .lt. 0. ) cycle

	  nodes = nkbnds(ibc)

	  !------------------------------------------------------------------
	  ! flag boundary nodes
	  !------------------------------------------------------------------

	  iflag = 0

	  do i=1,nodes
	    k = kbnds(ibc,i)
	    if( k > 0 ) iflag(k) = 1
	  end do

	  !------------------------------------------------------------------
	  ! set up ielfix and tnudgev/iuvfix
	  !------------------------------------------------------------------

	  do ie_mpi=1,nel

	    ie = ip_sort_elem(ie_mpi)

	    n = 0
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( iflag(k) .ne. 0 ) then
	        n = n + 1
	        ielfix(n,ie) = k
	      end if
	    end do

	    if( n .gt. 0 ) then			!nudging or fixing
	      ielfix(0,ie) = n        		!total number of nodes for ele
  	      if( tnudge == 0. ) then
		iuvfix(ie) = 1
  	      else if( tnudge > 0. ) then
		tnudgev(ie) = tnudge
	      else
		stop 'error stop bclfix_ini: internal error (1)'
	      end if
	    end if

	  end do

!------------------------------------------------------------------
! end loop over boundaries
!------------------------------------------------------------------

	end do

!------------------------------------------------------------------
! debug output
!------------------------------------------------------------------

	if( .not. bdebug ) return

	write(6,*) 'bclfix_ini has been initialized'

	nf = 0
	nn = 0
	do ie=1,nel
	  if( ielfix(0,ie) .ne. 0 ) then
	    if( iuvfix(ie) .ne. 0 ) then
	     write(6,*) 'Fix ie = ',ie,ieext(ie),'nn = ',ielfix(0,ie)
	     nf = nf + 1
	    else
	     write(6,*) 'Nudge ie = ',ie,ieext(ie),'nn = ',ielfix(0,ie)
	     nn = nn + 1
	    end if
	  end if
	end do

	write(6,*) '-------------'
	write(6,*) 'bclfix_ini has found ',nf,' elements to be fixed'
	write(6,*) 'bclfix_ini has found ',nn,' elements to be nudged'
	write(6,*) '-------------'

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

        end

!*****************************************************************

        subroutine bclfix

! fix or nudge  velocities on open boundaries

	use mod_bclfix
	use mod_internal
	use mod_geom_dynamic
	use mod_layer_thickness
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

        implicit none 

	double precision tnudge	!relaxation time for nudging [s]
	double precision tramp	!time for smooth init
	double precision h,alpha,uexpl,vexpl		!INTEL_BUG
        double precision u(nlvdi),v(nlvdi)		!INTEL_BUG
	double precision dfact				!INTEL_BUG
	!real tnudge	!relaxation time for nudging [s]!INTEL_BUG_OLD
	!real tramp	!time for smooth init		!INTEL_BUG_OLD
	!real h,alpha,uexpl,vexpl			!INTEL_BUG_OLD
        !real u(nlvdi),v(nlvdi)				!INTEL_BUG_OLD
	!real dfact					!INTEL_BUG_OLD

        integer ie,l,i,k,ii,n,ie_mpi
	integer lmax,lmin
	integer nbc
        
	integer nintp,nvar
	real cdef(2)
        double precision dtime0,dtime,ddtime
        integer, save, allocatable :: idvel(:)
        character*10 what

	integer nbnds

	integer, save :: icall = 0

        if( icall .eq. -1 ) return

!------------------------------------------------------------------
! initialize open boundary conditions
!------------------------------------------------------------------

	if( icall .eq. 0 ) then

          icall = -1
	  do ie = 1,nel
	    if ( ielfix(0,ie) .ne. 0 ) icall = 1
	  end do

          if( icall .eq. -1 ) return

          nbc = nbnds()
          allocate(idvel(nbc))
	  idvel = 0

	  call get_first_dtime(dtime0)
          nintp   = 2
          nvar    = 2
          cdef    = 0.
          what    = 'vel'
          call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv &
     &                      ,cdef,idvel)

	end if

!------------------------------------------------------------------
! read boundary velocities from file and store in u/vbound
!------------------------------------------------------------------

	call get_act_dtime(dtime)

        call bnds_read_new(what,idvel,dtime)
        call bnds_trans_new(what,idvel,dtime,1,nkn,nlv,nlvdi,ubound)
        call bnds_trans_new(what,idvel,dtime,2,nkn,nlv,nlvdi,vbound)

!------------------------------------------------------------------
! simulate smooth initial discharge with tramp
!------------------------------------------------------------------

        tramp = 43200.          !smooth init
        tramp = 0.              !smooth init

        alpha = 1.
        if( tramp .gt. 0. ) then
          call get_passed_dtime(ddtime)
          alpha = ddtime/tramp
          if( alpha .gt. 1. ) alpha = 1.
        end if

!------------------------------------------------------------------
! nudge or fix velocities in elements
!------------------------------------------------------------------

	do ie_mpi=1,nel
	  ie = ip_sort_elem(ie_mpi)
	  n = ielfix(0,ie)
	  if( n .gt. 0 ) then

	    lmax = ilhv(ie)
            lmin = jlhv(ie)

	    do l=1,lmax
	      u(l) = 0.
	      v(l) = 0.
	    end do

	    do i=1,n
	      k = ielfix(i,ie)
	      do l=lmin,lmax
	        u(l) = u(l) + ubound(l,k)
	        v(l) = v(l) + vbound(l,k)
	      end do
	    end do

	    do l=lmin,lmax
	      u(l) = u(l) / n
	      v(l) = v(l) / n
	    end do

	    tnudge = tnudgev(ie)

	    if (iuvfix(ie) .eq. 1 ) then	!fix velocities
              do l=lmin,lmax
                h = hdeov(l,ie)
                ulnv(l,ie) = u(l)*alpha
                vlnv(l,ie) = v(l)*alpha
                utlnv(l,ie) = ulnv(l,ie)*h
                vtlnv(l,ie) = vlnv(l,ie)*h
              end do
	    else if( tnudge > 0 ) then		!nudge velocities
	      do l=lmin,lmax
                h = hdeov(l,ie)
		dfact = alpha*h/tnudge		!INTEL_BUG
                uexpl = (ulnv(l,ie)-u(l))*dfact
                vexpl = (vlnv(l,ie)-v(l))*dfact
	        fxv(l,ie) = fxv(l,ie) + uexpl
	        fyv(l,ie) = fyv(l,ie) + vexpl
	      end do
	    end if
	  end if
	end do

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

	return
   99	continue
	stop 'error stop bclfix: internal error (1)'
        end

!*****************************************************************

