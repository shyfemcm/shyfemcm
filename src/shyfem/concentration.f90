
!--------------------------------------------------------------------------
!
!    Copyright (C) 1994,1996,1998-2019  Georg Umgiesser
!    Copyright (C) 2009  Andrea Cucco
!    Copyright (C) 2012,2014,2016  Christian Ferrarin
!    Copyright (C) 2013  Debora Bellafiore
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

! routines for concentration
!
! contents :
!
!-------------------------------------------------------------
!
! subroutine scal_adv(what,ivar
!     +                          ,scal,ids
!     +                          ,rkpar,wsink
!     +                          ,difhv,difv,difmol)
!		shell for scalar (for parallel version)

! subroutine scal_adv_nudge(what,ivar
!     +				,scal,ids
!     +				,rkpar,wsink
!     +                          ,difhv,difv,difmol
!     +				,sobs,robs,rtauv)
!		shell for scalar with nudging (for parallel version)
!
! subroutine scal_adv_fact(what,ivar,fact
!     +                          ,scal,ids
!     +                          ,rkpar,wsink,wsinkv,rload,load
!     +                          ,difhv,difv,difmol)
!		shell for scalar (for parallel version)
!		special version for cohesive sediments with factor
!
!-------------------------------------------------------------
!
! subroutine scal3sh(what,cnv,nlvddi,rcv,cobs,robs,rtauv,rkpar
!					,wsink,wsinkv,rload,load
!     +                                 ,difhv,difv,difmol)
!		shell for scalar T/D
!
! subroutine conz3d(cn1,co1
!     +                  ,ddt
!     +                  ,rkpar,difhv,difv
!     +                  ,difmol,cbound
!     +                  ,itvd,gradxv,gradyv
!     +                  ,cobs,robs,rtauv
!     +                  ,wsink,wsinkv
!     +                  ,rload,load
!     +                  ,azpar,adpar,aapar
!     +                  ,istot,isact
!     +                  ,nlvddi,nlv)
!
! subroutine conzstab(cn1,co1
!     +                  ,ddt
!     +                  ,robs,rtauv,wsink,wsinkv
!     +                  ,rkpar,difhv,difv
!     +                  ,difmol,azpar
!     +                  ,adpar,aapar
!     +                  ,sindex
!     +                  ,istot,isact
!     +                  ,nlvddi,nlv)
!
!-------------------------------------------------------------
!
! subroutine massconc(mode,cn,nlvddi,mass)
!		computes total mass of conc
!
! subroutine check_scal_bounds(cnv,cmin,cmax,eps,bstop)
!		checks if scalar is out of bounds
!
! subroutine assert_min_max_property(cnv,cov,sbconz,rmin,rmax,eps)
!		checks min/max property
!
! subroutine stb_histo(it,nlvddi,nkn,ilhkv,cwrite)
!		writes histogram info about stability index
!
!-------------------------------------------------------------
!
! notes:
!
!	dispersion:
!	  scal_adv(...,cnv,ids,...)
!
!	scal_adv(...,ids,...)
!	  bnds_trans_new(...,ids,r3v,...)	(transfers BC to matrix)
!         call scal3sh(...,cnv,r3v,...)
!	
!	scal3sh(...,cnv,r3v,...)
!         call make_scal_flux(...,cnv,r3v,sbconz,...)
!	  do
!           call conz3d(...,cnv,sbconz,...)
!           call assert_min_max_property(...,cnv,sbconz,...)
!           call bndo_setbc(what,nlvddi,cnv,rcv,uprv,vprv)
!	  end do
!
!-------------------------------------------------------------
!
! revision log :
!
! 09.01.1994	ggu	(from scratch)
! 19.01.1994	ggu	$$flux - flux conserving property
! 20.01.1994	ggu	$$iclin - iclin not used to compute volume
! 20.01.1994	ggu	$$lumpc - evaluate conz. nodewise
! 03.02.1994	ggu	$$itot0 - exception for itot=0 or 3
! 04.02.1994	ggu	$$fact3 - factor 3 missing in transport
! 04.02.1994	ggu	$$azpar - azpar used to compute transport
! 04.02.1994	ggu	$$condry - comute conz also in dry areas
! 07.02.1994	ggu	$$istot - istot for fractional time step
! 01.06.1994	ggu	restructured for 3-d model
! 18.07.1994	ggu	$$htop - use htop instead of htopo for mass cons.
! 09.04.1996	ggu	$$rvadj adjust rv in certain areas
! 14.08.1998	ggu	rkpar/rvpar -> chpar/cvpar
! 14.08.1998	ggu	use ilhkv to scan vertical levels on node
! 14.08.1998	ggu	$$LEV0 - bug fix : vertical level 0 used
! 19.08.1998	ggu	call to conzfi changed
! 26.08.1998	ggu	cleaned up conzsh
! 26.08.1998	ggu	conz uses zeov,zenv for water level
! 28.10.1999	ggu	names changed
! 07.03.2000	ggu	constant vertical eddy coefficient subst. with difv
! 20.06.2000	ggu	pass difmol to conz3d and use it
! 05.12.2001	ggu	variable horizontal diffusion, limit on dif.coef.
! 11.10.2002	ggu	file cleaned, t/shdif are set equal
! 11.10.2002	ggu	con3sh removed, conzstab better commented
! 14.10.2002	ggu	rstot re-introduced as rstol
! 09.09.2003	ggu	call to scal3sh changed -> new arg nlvddi
! 10.03.2004	ggu	call conwrite() to write stability param to nos file
! 13.03.2004	ggu	new boundary conditions through flux (cbound)
! 15.10.2004	ggu	boundary conditions back to old
! 02.12.2004	ggu	return also sindex in conzstab
! 17.01.2005	ggu	new routines with difhv
! 17.01.2005	ggu	get_stability and get_stab_index in this file
! 03.03.2005	ggu	new 3d boundary arrays implemented
! 16.08.2005	ggu	TVD algorithm implemented (gradxv,gradyv,grad_tvd,btvd)
! 04.11.2005	ggu	TVD changes from andrea integrated
! 07.11.2005	ggu	parameter itvd introduced for TVD
! 07.11.2005	ggu	sinking velocity wsink introduced in call to scal3sh
! 11.11.2005	ggu	bug fix in grad_tvd (ggx/ggy in layer loop now)
! 11.11.2005	ggu	new routine grad_2d()
! 16.02.2006	ggu	set w to zero at surface and bottom (WZERO)
! 23.03.2006	ggu	changed time step to real
! 08.08.2007	ggu	new parameter istot_max
! 23.08.2007	ggu	test for boundary nodes using routines in testbndo.h
! 18.09.2007	ggu	new subroutine check_scal
! 01.10.2007	ggu	Hack for ssurface -> set to 0 or -999 (temp)
! 17.03.2008	ggu	new open boundary routines introduced
! 08.04.2008	ggu	treatment of boundaries changed
! 22.04.2008	ggu	parallelization: scal_adv, scal_bnd
! 22.04.2008	ggu	local saux, sbflux, no explh, cl{c|m|p}e
! 22.04.2008	ggu	new routine scal_adv_fact for cohesive sediments
! 23.04.2008	ggu	call to bnds_set_def() changed
! 28.04.2008	ggu	rstol deleted
! 28.04.2008	ggu	new routines for stability, s/getistot deleted
! 28.04.2008	ggu	conz3sh into own file
! 24.06.2008	ggu	rstol re-introduced
! 08.11.2008	ggu	BUGFIX in conz3d (vertical velocity)
! 11.11.2008	ggu	conzstab cleaned
! 19.11.2008	ggu	changes in advect_stability() - incomplete
! 06.12.2008	ggu	in conzstab changed wprv, new routine write_elem_info()
! 27.01.2009	aac	bugs in TVD scheme fixed
! 24.03.2009	ggu	more bugs in TVD scheme fixed
! 31.03.2009	ggu	TVD algorithm tested and cleaned
! 20.04.2009	ggu	test for parallel execution (parallel_test)
! 13.10.2009	ggu	write_elem_info() substituted with check_elem()
! 12.11.2009	ggu	make_scal_flux() into internal time loop for stability
! 12.11.2009	ggu	in conz_stab loop over all k that are not z-boundaries
! 16.02.2010	ggu	use wdiff also in stab, use point sources in stab
! 16.02.2010	ggu	min/max property, sbconz is passed to conz3d
! 19.02.2010	ggu	restructured, stab routines into newstab.f
! 10.03.2010	ggu	in assert_min_max_property() check all nodes (also BC)
! 11.03.2010	ggu	in assert_min_max_property() do not check ibtyp=1
! 12.03.2010	ggu	in assert_min_max_property() limit error messages
! 22.03.2010	ggu	bug fix for evaporation (distr. sources) BUG_2010_01
! 15.12.2010	ggu	new routine vertical_flux_ie() for vertical tvd
! 26.01.2011	ggu	nudging implemented (scal_adv_nudge, cobs, robs)
! 16.02.2011	ggu	pass robs to info_stability()
! 23.03.2011	ggu	new parameter itvdv
! 25.03.2011	ggu	error check for aapar and itvdv
! 14.04.2011	ggu	changed VERS_6_1_22
! 31.05.2011	ggu	changed VERS_6_1_23
! 01.06.2011	ggu	wsink for stability integrated
! 12.07.2011	ggu	run over nlv, not nlvddi, vertical_flux() for lmax>1
! 15.07.2011	ggu	call vertical_flux() anyway (BUG)
! 30.03.2012	ggu	changed VERS_6_1_51
! 21.06.2012	ggu&ccf	variable vertical sinking velocity integrated
! 26.06.2012	ggu	changed VERS_6_1_55
! 25.10.2012	ggu	changed VERS_6_1_59
! 13.06.2013	ggu	changed VERS_6_1_65
! 12.09.2013	ggu	changed VERS_6_1_67
! 03.12.2013	ggu&dbf	bug fix for horizontal diffusion
! 15.05.2014	ggu	write min/max error only for levdbg >= 3
! 18.06.2014	ggu	changed VERS_6_1_77
! 27.06.2014	ggu	changed VERS_6_1_78
! 07.07.2014	ggu	changed VERS_6_1_79
! 10.07.2014	ggu	only new file format allowed
! 18.07.2014	ggu	changed VERS_7_0_1
! 13.10.2014	ggu	changed VERS_7_0_2
! 20.10.2014	ggu	accept ids from calling routines
! 22.10.2014	ccf	load in call to scal3sh
! 05.11.2014	ggu	changed VERS_7_0_5
! 12.12.2014	ggu	changed VERS_7_0_9
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 26.02.2015	ggu	changed VERS_7_1_5
! 20.05.2015	ggu	accumulate over nodes (for parallel version)
! 05.06.2015	ggu	changed VERS_7_1_12
! 10.07.2015	ggu	changed VERS_7_1_50
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_52
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 30.07.2015	ggu	changed VERS_7_1_83
! 23.09.2015	ggu	changed VERS_7_2_4
! 30.09.2015	ggu	routine cleaned, no reals in conz3d
! 23.10.2015	ggu	changed VERS_7_3_9
! 26.10.2015	ggu	critical omp sections introduced (eliminated data race)
! 26.10.2015	ggu	mass check only for levdbg > 2
! 09.11.2015	ggu	changed VERS_7_3_13
! 20.11.2015	ggu	changed VERS_7_3_15
! 16.12.2015	ggu	changed VERS_7_3_16
! 18.12.2015	ggu	changed VERS_7_3_17
! 01.04.2016	ggu	most big arrays moved from stack to allocatable
! 20.10.2016	ccf	pass rtauv for differential nudging
! 12.01.2017	ggu	changed VERS_7_5_21
! 05.12.2017	ggu	changed VERS_7_5_39
! 03.02.2018	ggu	sindex did not use rstol for stability
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.04.2018	ggu	changed VERS_7_5_43
! 03.04.2018	ggu	changed VERS_7_5_44
! 19.04.2018	ggu	changed VERS_7_5_45
! 23.04.2018	ggu	exchange mpi inside loop for istot>1
! 11.05.2018	ggu	compute only unique nodes (needed for zeta layers)
! 30.05.2018	ggu	better debug output in conzstab (idtstb,itmstb)
! 01.06.2018	ggu	stability of scalar revised - aa > 0 possible again
! 01.06.2018	ggu	implicit nudging (relaxation) (ANT)
! 06.07.2018	ggu	changed VERS_7_5_48
! 11.10.2018	ggu	caux substituted with load,cobs,rtauv (bug inout)
! 14.02.2019	ggu	check for negative scalar
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 21.05.2019	ggu	changed VERS_7_5_62
! 31.05.2021	ggu	possibly write stability index to inf file
! 15.02.2022	ggu	new routine limit_scalar() implemented
! 19.02.2022	ggu	write nodes where limit is exceeded
! 06.04.2022	ggu	adapted to regular assembling over elems (ie_mpi)
! 07.04.2022	ggu	debug code (kdebug)
! 03.05.2022	ggu	exchanging twice around bndo_setbc() -> improve
! 09.05.2023    lrp     introduce top layer index variable
! 31.05.2023    ggu     in conzstab, use ie_mpi, run over nkn_unique (bug fix)
! 27.09.2024    ggu     btvddebug introduced
!
!*********************************************************************

	subroutine scal_adv(what,ivar &
     &				,scal,ids &
     &				,rkpar,wsink &
     &                          ,difhv,difv,difmol)

! shell for scalar (for parallel version)

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

        character*(*) what
	integer ivar
        real scal(nlvdi,nkn)
        integer ids(*)
        real rkpar
	real wsink
        real difhv(nlvdi,nel)
	real difv(0:nlvdi,nkn)
        real difmol

        real, allocatable :: r3v(:,:)
        real, allocatable :: load(:,:)
        real, allocatable :: cobs(:,:)
        real, allocatable :: rtauv(:,:)
        real, allocatable :: wsinkv(:,:)

	double precision dtime
	real robs,rload
        integer iwhat,ichanm
	character*10 whatvar,whataux

	allocate(r3v(nlvdi,nkn))
	allocate(load(nlvdi,nkn))
	allocate(cobs(nlvdi,nkn))
	allocate(rtauv(nlvdi,nkn))
	allocate(wsinkv(0:nlvdi,nkn))

	robs = 0.
	rload = 0.
	r3v = 0.
	load = 0.
	cobs = 0.
	rtauv = 0.
	wsinkv = 1.

!--------------------------------------------------------------
! make identifier for variable
!--------------------------------------------------------------

!$OMP CRITICAL
	whatvar = what
	if( ivar .ne. 0 ) then
          write(whataux,'(i2)') ivar
          whatvar = what // '_' // whataux
	end if
        iwhat = ichanm(whatvar)
!$OMP END CRITICAL

!--------------------------------------------------------------
! transfer boundary conditions of var ivar to 3d matrix r3v
!--------------------------------------------------------------

	call get_act_dtime(dtime)

	call bnds_trans_new(whatvar(1:iwhat) &
     &			,ids,dtime,ivar,nkn,nlv,nlvdi,r3v)

!--------------------------------------------------------------
! do advection and diffusion
!--------------------------------------------------------------

        call scal3sh(whatvar(1:iwhat) &
     &				,scal,nlvdi &
     &                          ,r3v,cobs,robs,rtauv &
     &				,rkpar,wsink,wsinkv,rload,load &
     &                          ,difhv,difv,difmol)

	deallocate(r3v,load,cobs,rtauv,wsinkv)

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*********************************************************************

	subroutine scal_adv_nudge(what,ivar &
     &				,scal,ids &
     &				,rkpar,wsink &
     &                          ,difhv,difv,difmol &
     &				,sobs,robs,rtauv)

! shell for scalar with nudging (for parallel version)

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

        character*(*) what
	integer ivar
        real scal(nlvdi,nkn)
        integer ids(*)
        real rkpar
	real wsink
        real difhv(nlvdi,nel)
	real difv(0:nlvdi,nkn)
        real difmol
	real sobs(nlvdi,nkn)		!observations
	real robs
	real rtauv(nlvdi,nkn)		!varible relaxation coefficient

        real, allocatable :: r3v(:,:)
        real, allocatable :: wsinkv(:,:)
	real, allocatable :: load(:,:)

	integer ierr,l,k,lmax
	real eps
	real rload
	double precision dtime
        integer iwhat,ichanm
	character*10 whatvar,whataux

	allocate(r3v(nlvdi,nkn),load(nlvdi,nkn))
	allocate(wsinkv(0:nlvdi,nkn))

	rload = 0.
	r3v = 0.
	wsinkv = 1.
	load = 0.

!--------------------------------------------------------------
! make identifier for variable
!--------------------------------------------------------------

	whatvar = what
	if( ivar .ne. 0 ) then
          write(whataux,'(i2)') ivar
          whatvar = what // whataux
	end if
        iwhat = ichanm(whatvar)

!--------------------------------------------------------------
! transfer boundary conditions of var ivar to 3d matrix r3v
!--------------------------------------------------------------

	call get_act_dtime(dtime)

	call bnds_trans_new(whatvar(1:iwhat) &
     &			,ids,dtime,ivar,nkn,nlv,nlvdi,r3v)

!--------------------------------------------------------------
! do advection and diffusion
!--------------------------------------------------------------

        call scal3sh(whatvar(1:iwhat) &
     &				,scal,nlvdi &
     &                          ,r3v,sobs,robs,rtauv &
     &				,rkpar,wsink,wsinkv,rload,load &
     &                          ,difhv,difv,difmol)

	deallocate(r3v,load)
	deallocate(wsinkv)

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*********************************************************************

	subroutine scal_adv_fact(what,ivar,fact &
     &				,scal,ids &
     &				,rkpar,wsink,wsinkv,rload,load &
     &                          ,difhv,difv,difmol)

! shell for scalar (for parallel version)
!
! special version with factor for BC, variable sinking velocity and loads

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

        character*(*) what
	integer ivar
	real fact			!factor for boundary condition
        real scal(nlvdi,nkn)
        integer ids(*)
        real rkpar
	real wsink
	real wsinkv(0:nlvdi,nkn)
	real rload			!load factor (1 for load given)
	real load(nlvdi,nkn)		!load [kg/s]
        real difhv(nlvdi,nel)
	real difv(0:nlvdi,nkn)
        real difmol

        real, allocatable :: r3v(:,:)
        real, allocatable :: cobs(:,:)
        real, allocatable :: rtauv(:,:)

	double precision dtime
	real robs
        integer iwhat,ichanm
	character*20 whatvar,whataux

	allocate(r3v(nlvdi,nkn))
	allocate(cobs(nlvdi,nkn))
	allocate(rtauv(nlvdi,nkn))

	robs = 0.
	r3v = 0.
	cobs = 0.
	rtauv = 0.

!--------------------------------------------------------------
! make identifier for variable
!--------------------------------------------------------------

	whatvar = what
	if( ivar .ne. 0 ) then
          write(whataux,'(i2)') ivar
          whatvar = what // whataux
	end if
        iwhat = ichanm(whatvar)

!--------------------------------------------------------------
! transfer boundary conditions of var ivar to 3d matrix r3v
!--------------------------------------------------------------

	call get_act_dtime(dtime)

	call bnds_trans_new(whatvar(1:iwhat) &
     &			,ids,dtime,ivar,nkn,nlv,nlvdi,r3v)

!--------------------------------------------------------------
! multiply boundary condition with factor
!--------------------------------------------------------------

	if( fact .ne. 1. ) then
	  call mult_scal_bc(r3v,fact)
	end if

!--------------------------------------------------------------
! do advection and diffusion
!--------------------------------------------------------------

        call scal3sh(whatvar(1:iwhat) &
     &				,scal,nlvdi &
     &                          ,r3v,cobs,robs,rtauv &
     &				,rkpar,wsink,wsinkv,rload,load &
     &                          ,difhv,difv,difmol)

	deallocate(r3v,cobs,rtauv)

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*********************************************************************

	subroutine scal_bnd(what,t,bnd3)

! sets boundary conditions for scalar - not used anymore - to be deleted

	implicit none

        character*(*) what
	real t
	real bnd3(1,1)

	stop 'error stop: call to scal_bnd not allowed'

	end

!*********************************************************************
!*********************************************************************
!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine scal3sh(what,cnv,nlvddi,rcv,cobs,robs,rtauv,rkpar &
     &					,wsink,wsinkv,rload,load &
     &					,difhv,difv,difmol)

! shell for scalar T/D

	use mod_hydro_print
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
	use shympi
	use shympi_debug
	use mkonst

	implicit none

! arguments
        character*(*) what
        real cnv(nlvddi,nkn)
	integer nlvddi		!vertical dimension
        real rcv(nlvddi,nkn)	!boundary condition (value of scalar)
	real cobs(nlvddi,nkn)	!observations (for nudging)
	real robs		!use nudging
	real rtauv(nlvddi,nkn)	!observations (for nudging)
        real rkpar
	real wsink
	real wsinkv(0:nlvddi,nkn)
	real rload
	real load(nlvddi,nkn)		!load [kg/s]
        real difhv(nlvddi,nel)
	real difv(0:nlvddi,nkn)
        real difmol
! parameters
	integer istot_max
	!parameter ( istot_max = 100 )
	!parameter ( istot_max = 200 )
	!parameter ( istot_max = 300 )
	parameter ( istot_max = 1000 )
! common
! local
        real, allocatable :: saux(:,:)		!aux array
        real, allocatable :: sbflux(:,:)	!flux boundary conditions
        real, allocatable :: sbconz(:,:)	!conz boundary conditions
	real, allocatable :: gradxv(:,:)	!gradient in x for tvd
	real, allocatable :: gradyv(:,:)	!gradient in y for tvd

	logical btvd,btvd1,btvd2,btvddebug
	integer isact
	integer istot
	integer itvd
	integer itvdv
	integer iuinfo
	integer iunit,k
	integer levdbg
        real dt,dtstep
	real eps
        real sindex
	real mass,massold,massdiff
	real azpar,adpar,aapar
	real ssurface
	real wsinkl				!local sinking
	double precision dtime,dtime_act
	character*20 aline
! function
	real getpar
	integer ipint

!-------------------------------------------------------------
! start of routine
!-------------------------------------------------------------

	iunit = 888 + my_id
	iunit = 0

!-------------------------------------------------------------
! initialization
!-------------------------------------------------------------

	call getaz(azpar)
	adpar=getpar('adpar')
	aapar=getpar('aapar')
	itvd=nint(getpar('itvd'))	!horizontal tvd scheme
	itvdv=nint(getpar('itvdv'))	!vertical tvd scheme
	levdbg = nint(getpar('levdbg'))
	btvd = itvd .gt. 0
	btvd1 = itvd .eq. 1
	btvd2 = itvd .eq. 2
	btvddebug = .true.
	btvddebug = btvddebug .and. btvd2

	allocate(saux(nlvddi,nkn))
	allocate(sbflux(nlvddi,nkn),sbconz(nlvddi,nkn))
	allocate(gradxv(nlvddi,nkn),gradyv(nlvddi,nkn))

	saux = 0.
	sbflux = 0.
	sbconz = 0.
	gradxv = 0.
	gradyv = 0.

	wsinkl = wsink
	!wsinkl = wsink * 10.

!$OMP CRITICAL
	if(shympi_is_master()) then
          call getinfo(iuinfo)  !unit number of info file
	end if
!$OMP END CRITICAL

	eps = 1.e-5
	eps = 1.e-4
	eps = 1.e-2

!-------------------------------------------------------------
! check stability criterion -> set istot
!-------------------------------------------------------------

	call get_timestep(dt)
	call get_act_dtime(dtime)
	call get_act_timeline(aline)

	saux = 0.
	call scalar_stability(dt,robs,rtauv,wsinkl,wsinkv,rkpar, &
     &					sindex,istot,saux)

!$OMP CRITICAL
	if(shympi_is_master()) then
          !write(iuinfo,*) 'stability_',what,': ',aline,sindex,istot
	end if
!$OMP END CRITICAL

        if( istot .gt. istot_max ) then
	    call scalar_info_stability(dt,robs,rtauv,wsinkl,wsinkv &
     &				,rkpar,sindex,istot,saux)
            write(6,*) 'istot  = ',istot,'   sindex = ',sindex
            stop 'error stop scal3sh: istot index too high'
        end if

!-------------------------------------------------------------
! set up TVD scheme
!-------------------------------------------------------------

	!call tvd_init(itvd)	!is called in shyfem

!-------------------------------------------------------------
! set up flux boundary conditions (temporary) -> put in sbflux
!-------------------------------------------------------------

	ssurface = 0.	!FIXME - HACK
	if( what .eq. 'temp' ) ssurface = -999.

! make_scal_flux has been moved inside time loop of isact (see below)
! to increase stability when treating outflow flux boundary
! this is needed only for discharge < 0 in order to use always the
! ambient tracer concentration

	!write(iunit,*) 'istot = ',istot,dtime

!-------------------------------------------------------------
! transport and diffusion
!-------------------------------------------------------------

	call massconc(-1,cnv,nlvddi,massold)

	do isact=1,istot

	  dtstep = -((istot-isact)*dt)/istot
	  dtime_act = dtime + dtstep		!why is dtstep negative?

	  call make_scal_flux(what,isact,rcv,cnv,sbflux,sbconz,ssurface)
	  !call check_scal_flux(what,cnv,sbconz)

	  if( what /= 'temp' ) then
	    where( cnv < 1.e-15 ) cnv = 0.
	    if( any(cnv < 0.) ) then
	      write(6,*) 'negative scalar: ',what,isact &
     &		,count(cnv<0.),minval(cnv)
	    end if
	  end if

	  if( btvd1 ) call tvd_grad_3d(cnv,gradxv,gradyv,saux,nlvddi)
	  if( btvddebug ) call tvd_debug_initialize(dtime_act,what,isact)

          !call conz3d_orig( &
          call conz3d_omp(       &
     &           cnv &
     &          ,saux &
     &          ,dt &
     &          ,rkpar,difhv,difv,difmol &
     &          ,sbconz &
     &		,itvd,itvdv,gradxv,gradyv &
     &		,cobs,robs,rtauv &
     &		,wsinkl,wsinkv &
     &		,rload,load &
     &		,azpar,adpar,aapar &
     &          ,istot,isact &
     &          ,nlvddi,nlv &
     &               )

	  call assert_min_max_property(cnv,saux,sbconz,gradxv,gradyv,eps)

	  call limit_scalar(what,dtstep,cnv)

          call shympi_exchange_3d_node(cnv)
          call bndo_setbc(what,nlvddi,cnv,rcv,uprv,vprv)
          call shympi_exchange_3d_node(cnv)

	  if( btvddebug ) call tvd_debug_finalize
	end do

        !if( shympi_is_parallel() .and. istot > 1 ) then
        !  write(6,*) 'cannot handle istot>1 with mpi yet'
        !  stop 'error stop scal3sh: istot>1'
        !end if
        !call shympi_comment('exchanging scalar: '//trim(what))
        !call shympi_exchange_3d_node(cnv)
        !call shympi_barrier

!-------------------------------------------------------------
! check total mass
!-------------------------------------------------------------

	if( levdbg > 2 ) then
	  call massconc(+1,cnv,nlvddi,mass)
	  massdiff = mass - massold
!$OMP CRITICAL
          if(shympi_is_master())then
	    write(iuinfo,1000) 'scal3sh_'//trim(what)//': ',aline &
     &                          ,mass,massold,massdiff
	  end if
!$OMP END CRITICAL
	end if

	deallocate(saux)
	deallocate(sbflux,sbconz)
	deallocate(gradxv,gradyv)

	if( iunit > 0 ) flush(iunit)

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

        return
 !1000   format(a,a,a,2i10,3d13.5)
 1000   format(a,a,3d13.5)
	end

!**************************************************************

        subroutine conz3d_orig(cn1,co1 &
     &			,ddt &
     &                  ,rkpar,difhv,difv &
     &			,difmol,cbound &
     &		 	,itvd,itvdv,gradxv,gradyv &
     &			,cobs,robs,rtauv &
     &			,wsink,wsinkv &
     &			,rload,load &
     &			,azpar,adpar,aapar &
     &			,istot,isact &
     &			,nlvddi,nlev)
!
! computes concentration
!
! cn     new concentration
! co     old concentration              !not used !FIXME
! caux   aux vector
! clow	 lower diagonal of vertical system
! chig	 upper diagonal of vertical system
! ddt    time step
! rkpar  horizontal turbulent diffusivity
! difhv  horizontal turbulent diffusivity (variable between elements)
! difv   vertical turbulent diffusivity
! difmol vertical molecular diffusivity
! cbound boundary condition (mass flux) [kg/s] -> now concentration [kg/m**3]
! itvd	 type of horizontal transport algorithm used
! itvdv	 type of vertical transport algorithm used
! gradxv,gradyv  gradient vectors for TVD algorithm
! cobs	 observations for nudging
! robs	 use observations for nuding (real)
! rtauv	 relaxation time for nuding (real)
! wsink	 factor for settling velocity
! wsinkv variable settling velocity [m/s]
! rload	 factor for loading
! load   load (source or sink) [kg/s]
! azpar  time weighting parameter
! adpar  time weighting parameter for vertical diffusion (ad)
! aapar  time weighting parameter for vertical advection (aa)
! istot	 total inter time steps
! isact	 actual inter time step
! nlvddi	 dimension in z direction
! nlv	 actual needed levels
!
! solution of purely diffusional part :
!
! dC/dt = a*laplace(C)    with    c(x,0+)=delta(x)
!
! C(x,t) =  (4*pi*a*t)**(-n/2) * exp( -|x|**2/(4*a*t) )
!
! for n-dimensions and
!
! C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )
!
! for 1 dimension
!
! the solution is normalized, i.e.  int(C(x,t)dx) = 1 over the whole area
!
! DPGGU -> introduced double precision to stabilize solution

	use mod_bound_geom
	use mod_geom
	use mod_depth
	use mod_layer_thickness
	use mod_diff_aux
	use mod_bound_dynamic
	use mod_area
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin
	use shympi

	implicit none
!
! arguments
	integer nlvddi,nlev
        real cn1(nlvddi,nkn),co1(nlvddi,nkn)		!DPGGU
        real difv(0:nlvddi,nkn)
        real difhv(nlvddi,nel)
	real difmol
        real cbound(nlvddi,nkn)
	integer itvd
	integer itvdv
	real gradxv(nlvddi,nkn)
	real gradyv(nlvddi,nkn)
	real cobs(nlvddi,nkn)
	real robs
	real rtauv(nlvddi,nkn)
	real wsink
	real wsinkv(0:nlvddi,nkn)
	real rload
        real load(nlvddi,nkn)                      !ccf_load
        real ddt,rkpar
        real azpar,adpar,aapar			!$$azpar
	integer istot,isact
! local
	logical bdebug,bdebug1,btvdv
	integer k,ie,ii,l,iii,ll,ibase,ntot,ie_mpi
	integer lstart
	integer ilevel,jlevel
	integer itot,isum	!$$flux
	logical berror
	integer kn(3)
        integer ip(3,3)
        integer n,i,ipp
	integer icount,icc,iunit
	integer elems(maxlnk)
        double precision mflux,qflux,cconz
	double precision loading
	double precision us,vs
	double precision az,azt
	double precision aa,aat,ad,adt
	double precision an,ant
	double precision aj,rk3,rv,aj4,aj12
	double precision hmed,hmbot,hmtop
	double precision hmotop,hmobot,hmntop,hmnbot
	double precision rvptop,rvpbot
	double precision dt,w,aux
	double precision b(3),c(3),f(3)
	double precision rso,rsn,rsot,rsnt,rstot
	double precision hn,ho

	double precision cn(nlvddi,nkn)		!DPGGU	!FIXME
	double precision co(nlvddi,nkn)
	double precision cdiag(nlvddi,nkn)
	double precision clow(nlvddi,nkn)
	double precision chigh(nlvddi,nkn)

	double precision cexpl
	double precision cbm,ccm
	double precision fw(3),fd(3)
	double precision fl(3),fnudge(3)
	double precision flux_tot,flux_tot1,flux_top,flux_bot
        double precision wdiff(3),waux
! local (new)
	double precision clc(nlvddi,3), clm(nlvddi,3), clp(nlvddi,3)
	double precision cle(nlvddi,3)

	double precision cclc(nlvddi,3,nel)
	double precision cclm(nlvddi,3,nel)
	double precision cclp(nlvddi,3,nel)
	double precision ccle(nlvddi,3,nel)

	double precision cl(0:nlvddi+1,3)
	double precision wl(0:nlvddi+1,3)
	double precision vflux(0:nlvddi+1,3)
	double precision cob(0:nlvddi+1,3)
	double precision rtau(0:nlvddi+1,3)
	double precision finu(0:nlvddi+1,3)

	double precision hdv(0:nlvddi+1)
	double precision haver(0:nlvddi+1)
	double precision hnew(0:nlvddi+1,3)
	double precision hold(0:nlvddi+1,3)
	double precision htnew(0:nlvddi+1,3)
	double precision htold(0:nlvddi+1,3)
	double precision present(0:nlvddi+1)

	double precision cauxn(nlvddi)	!FIXME
	double precision cauxd(nlvddi)
	double precision cauxh(nlvddi)
	double precision cauxl(nlvddi)
! tvd
	logical btvd
	integer ic,kc,id,kdebug,ippp
	integer ies
	integer iext
	double precision fls(3)
        double precision wws

! functions
	integer ipint,ieint
	integer ipext,ieext
	integer ithis

        if(nlv.ne.nlev) stop 'error stop conz3d_orig: nlv/=nlev'

!----------------------------------------------------------------
! initialize variables and parameters
!----------------------------------------------------------------

	kdebug = ipint(3371)
	kdebug = -1
	iunit = 888 + my_id
	!write(iunit,*) kdebug,nkn_inner,nkn

        bdebug1 = .true.
        bdebug1 = .false.
	bdebug=.false.
	berror=.false.

	btvd = itvd .gt. 0
	btvdv = itvdv .gt. 0

	az=azpar		!$$azpar
	azt=1.-az
	ad=adpar
	adt=1.-ad
	aa=aapar
	aat=1.-aa
	an=0.			!nudging parameter (1 for implicit) $ANT
	ant=1.-an

	rstot=istot		!$$istot
	rso=(isact-1)/rstot
	rsn=(isact)/rstot
	rsot=1.-rso
	rsnt=1.-rsn

	dt=ddt/rstot

	if( btvdv .and. aa .ne. 0. ) then
	  write(6,*) 'aapar = ',aapar,'  itvdv = ',itvdv
	  write(6,*) 'Cannot use implicit vertical advection'
	  write(6,*) 'together with vertical TVD scheme.'
	  write(6,*) 'Please set either aapar = 0 (explicit) or'
	  write(6,*) 'itvdv = 0 (no vertical TVD) in the STR file.'
	  stop 'error stop conz3d: vertical tvd scheme'
	end if

!	----------------------------------------------------------------
!	global arrays for accumulation of implicit terms
!	----------------------------------------------------------------

	do k=1,nkn
          do l=1,nlv
	    co1(l,k)=cn1(l,k)	!COLD
	    co(l,k)=cn1(l,k)	!DPGGU
            cn(l,k)=0.
            cdiag(l,k)=0.
            clow(l,k)=0.
            chigh(l,k)=0.
          end do
	end do

	if( kdebug > 0 ) then
	  write(iunit,*) 'init --------',ipext(kdebug),ilhkv(kdebug)
	  write(iunit,*) co1(1,kdebug)
	  write(iunit,*) co(1,kdebug)
	end if
	
!	----------------------------------------------------------------
!	aux elements inside element
!	----------------------------------------------------------------

!	these are aux arrays (bigger than needed) to avoid checking for
!	what layer we are in -> we never get out of bounds

        do l=0,nlv+1
	  hdv(l) = 0.		!layer thickness
          haver(l) = 0.
	  present(l) = 0.	!1. if layer is present
	  do ii=1,3
	    hnew(l,ii) = 0.	!as hreal but with zeta_new
	    hold(l,ii) = 0.	!as hreal but with zeta_old
	    cl(l,ii) = 0.	!concentration in layer
	    wl(l,ii) = 0.	!vertical velocity
	    vflux(l,ii) = 0.	!vertical velocity
	  end do
	end do

!	these are the local arrays for accumulation of implicit terms
!	(maybe we do not need them, but just to be sure...)
!	after accumulation we copy them onto the global arrays

        do l=1,nlv
	  do ii=1,3
	    cle(l,ii) = 0.
	    clc(l,ii) = 0.
	    clm(l,ii) = 0.
	    clp(l,ii) = 0.
	  end do
	end do

!	----------------------------------------------------------------
!	define vertical velocities
!	----------------------------------------------------------------

!----------------------------------------------------------------
! loop over elements
!----------------------------------------------------------------

        do ie_mpi=1,nel

	ie = ip_sort_elem(ie_mpi)

	do ii=1,3
          k=nen3v(ii,ie)
	  kn(ii)=k
	  b(ii)=ev(ii+3,ie)
	  c(ii)=ev(ii+6,ie)
	end do

	aj=ev(10,ie)    !area of triangle / 12
	aj4=4.*aj
	aj12=12.*aj
        ilevel=ilhv(ie)
        jlevel=jlhv(ie)

!	----------------------------------------------------------------
!	set up vectors for use in assembling contributions
!	----------------------------------------------------------------

        do l=jlevel,ilevel
	  hdv(l) = hdeov(l,ie)		!use old time step -> FIXME
          !haver(l) = 0.5 * ( hdeov(l,ie) + hdenv(l,ie) )
          haver(l) = rso*hdenv(l,ie) + rsot*hdeov(l,ie)
	  present(l) = 1.
	  do ii=1,3
	    k=kn(ii)
	    hn = hdknv(l,k)		! there are never more layers in ie
	    ho = hdkov(l,k)		! ... than in k
            htold(l,ii) = ho
            htnew(l,ii) = hn
	    hold(l,ii) = rso * hn + rsot * ho
	    hnew(l,ii) = rsn * hn + rsnt * ho
	    cl(l,ii) = co(l,k)
	    cob(l,ii) = cobs(l,k)	!observations
	    rtau(l,ii) = rtauv(l,k)	!observations
	    wl(l,ii) = wlnv(l,k) - wsink * wsinkv(l,k)
	  end do
	end do

	do l=ilevel+1,nlv
	  present(l) = 0.
	end do

!	----------------------------------------------------------------
!	set vertical velocities in surface and bottom layer
!	----------------------------------------------------------------

!	we do not set wl(0,ii) because otherwise we loose concentration
!	through surface
!
!	we set wl(ilevel,ii) to 0 because we are on the bottom
!	and there should be no contribution from this element
!	to the vertical velocity

	do ii=1,3
	  wl(ilevel,ii) = 0.
	end do

!	----------------------------------------------------------------
!	compute vertical fluxes (w/o vertical TVD scheme)
!	----------------------------------------------------------------

	wws = 0.	!sinking velocity alread in wl
	call vertical_flux_ie(btvdv,ie,ilevel,jlevel, &
     &			      dt,wws,cl,wl,hold,vflux)

!----------------------------------------------------------------
! loop over levels
!----------------------------------------------------------------

        do l=jlevel,ilevel

        us=az*utlnv(l,ie)+azt*utlov(l,ie)             !$$azpar
        vs=az*vtlnv(l,ie)+azt*vtlov(l,ie)

        rk3 = 3. * rkpar * difhv(l,ie)

	cbm=0.
	ccm=0.
	itot=0
	isum=0
	do ii=1,3
	  k=kn(ii)
	  f(ii)=us*b(ii)+vs*c(ii)	!$$azpar
	  if(f(ii).lt.0.) then	!flux out of node
	    itot=itot+1
	    isum=isum+ii
	  end if
	  cbm=cbm+b(ii)*cl(l,ii)
	  ccm=ccm+c(ii)*cl(l,ii)

!	  ----------------------------------------------------------------
!	  initialization to be sure we are in a clean state
!	  ----------------------------------------------------------------

	  fw(ii) = 0.
	  cle(l,ii) = 0.
	  clc(l,ii) = 0.
	  clm(l,ii) = 0.
	  clp(l,ii) = 0.

!	  ----------------------------------------------------------------
!	  contributions from horizontal diffusion
!	  ----------------------------------------------------------------

          waux = 0.
          do iii=1,3
            waux = waux + wdifhv(iii,ii,ie) * cl(l,iii)
          end do
          wdiff(ii) = waux

!	  ----------------------------------------------------------------
!	  contributions from vertical diffusion
!	  ----------------------------------------------------------------

!	  in fd(ii) is explicit contribution
!	  the sign is for the term on the left side, therefore
!	  fd(ii) must be subtracted from the right side
!
!	  maybe we should use real layer thickness, or even the
!	  time dependent layer thickness

	  rvptop = difv(l-1,k) + difmol
	  rvpbot = difv(l,k) + difmol
	  !hmtop = 2. * rvptop * present(l-1) / (hdv(l-1)+hdv(l))
	  !hmbot = 2. * rvpbot * present(l+1) / (hdv(l)+hdv(l+1))
	  hmotop =2.*rvptop*present(l-1)/(hold(l-1,ii)+hold(l,ii))
	  hmobot =2.*rvpbot*present(l+1)/(hold(l,ii)+hold(l+1,ii))
	  hmntop =2.*rvptop*present(l-1)/(hnew(l-1,ii)+hnew(l,ii))
	  hmnbot =2.*rvpbot*present(l+1)/(hnew(l,ii)+hnew(l+1,ii))

	  fd(ii) = adt * (  &
     &			(cl(l,ii)-cl(l+1,ii))*hmobot - &
     &			(cl(l-1,ii)-cl(l,ii))*hmotop &
     &			  )

	  clc(l,ii) = clc(l,ii) + ad * ( hmntop + hmnbot )
	  clm(l,ii) = clm(l,ii) - ad * ( hmntop )
	  clp(l,ii) = clp(l,ii) - ad * ( hmnbot )

!	  ----------------------------------------------------------------
!	  contributions from vertical advection
!	  ----------------------------------------------------------------

!	  in fw(ii) is explicit contribution
!	  the sign is for the term on the left side, therefore
!	  fw(ii) must be subtracted from the right side
!	  fw is positive if flux out of node
!
!	  if we are in last layer, w(l,ii) is zero
!	  if we are in first layer, w(l-1,ii) is zero (see above)

	  w = wl(l-1,ii)		!top of layer
	  if( l .eq. jlevel ) w = 0.	!surface -> no transport (WZERO)
	  if( w .ge. 0. ) then
	    fw(ii) = aat*w*cl(l,ii)
	    flux_top = w*cl(l,ii)
	    clc(l,ii) = clc(l,ii) + aa*w
	  else
	    fw(ii) = aat*w*cl(l-1,ii)
	    flux_top = w*cl(l-1,ii)
	    clm(l,ii) = clm(l,ii) + aa*w
	  end if

	  w = wl(l,ii)			!bottom of layer
	  if( l .eq. ilevel ) w = 0.	!bottom -> handle flux elsewhere (WZERO)
	  if( w .gt. 0. ) then
	    fw(ii) = fw(ii) - aat*w*cl(l+1,ii)
	    flux_bot = w*cl(l+1,ii)
	    clp(l,ii) = clp(l,ii) - aa*w
	  else
	    fw(ii) = fw(ii) - aat*w*cl(l,ii)
	    flux_bot = w*cl(l,ii)
	    clc(l,ii) = clc(l,ii) - aa*w
	  end if

	  flux_tot1 = aat * ( flux_top - flux_bot )
	  flux_tot = aat * ( vflux(l-1,ii) - vflux(l,ii) )

	  fw(ii) = flux_tot
	end do

!	----------------------------------------------------------------
!	contributions from horizontal advection (only explicit)
!	----------------------------------------------------------------
!
!	f(ii) > 0 ==> flux into node ii
!	itot=1 -> flux out of one node
!		compute flux with concentration of this node
!	itot=2 -> flux into one node
!		for flux use conz. of the other two nodes and
!		minus the sum of these nodes for the flux of this node

	if(itot.eq.1) then	!$$flux
	  fl(1)=f(1)*cl(l,isum)
	  fl(2)=f(2)*cl(l,isum)
	  fl(3)=f(3)*cl(l,isum)
	else if(itot.eq.2) then
	  isum=6-isum
	  fl(1)=f(1)*cl(l,1)
	  fl(2)=f(2)*cl(l,2)
	  fl(3)=f(3)*cl(l,3)
	  fl(isum) = 0.
	  fl(isum) = -(fl(1)+fl(2)+fl(3))
	  isum=6-isum		!reset to original value
	else			!exception	$$itot0
	  fl(1)=0.
	  fl(2)=0.
	  fl(3)=0.
	end if

!	----------------------------------------------------------------
!	horizontal TVD scheme start
!	----------------------------------------------------------------

        if( btvd ) then
	  iext = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( is_external_boundary(k) ) iext = iext + 1
	  end do

          if( iext .eq. 0 ) then
	    call tvd_fluxes(ie,l,itot,isum,dt,cl,co1,gradxv,gradyv,f,fl)
	  end if
	end if

!	----------------------------------------------------------------
!	horizontal TVD scheme finish
!	----------------------------------------------------------------

!	----------------------------------------------------------------
!	contributions from nudging
!	----------------------------------------------------------------

	do ii=1,3
	  fnudge(ii) = robs * rtau(l,ii) * ( cob(l,ii) - ant * cl(l,ii) )
	  finu(l,ii) = an * robs * rtau(l,ii)	!implicit contribution
	end do

!	----------------------------------------------------------------
!	sum explicit contributions
!	----------------------------------------------------------------

	do ii=1,3
          hmed = haver(l)                    !new ggu   !HACK
	  cexpl = aj4 * ( hold(l,ii)*cl(l,ii) &
     &				+ dt *  (  &
     &					    hold(l,ii)*fnudge(ii) &
     &					  + 3.*fl(ii)  &
     &					  - fw(ii) &
     &					  - rk3*hmed*wdiff(ii) &
     &					  - fd(ii) &
     &					) &
     &		         )
	  cle(l,ii) = cle(l,ii) + cexpl
	  k = nen3v(ii,ie)
	  if( k == kdebug ) then
	    write(iunit,*) 'cexpl'
	    write(iunit,*) hold(l,ii)*cl(l,ii)
	    write(iunit,*) hold(l,ii)*fnudge(ii)
	    write(iunit,*) fl(ii)
	    write(iunit,*) 'fw: ',fw(ii)
	    write(iunit,*) fd(ii),rk3*hmed*wdiff(ii)
	  end if
	  !k=kn(ii)
	  !cn(l,k) = cn(l,k) + cexpl
	end do

	end do		! loop over l

!----------------------------------------------------------------
! end of loop over l
!----------------------------------------------------------------

!----------------------------------------------------------------
! set up implicit contributions
!----------------------------------------------------------------

! cdiag is diagonal of tri-diagonal system
! chigh is high (right) part of tri-diagonal system
! clow is low (left) part of tri-diagonal system
!
! clp -> bottom
! clm -> top

	do ii=1,3
	  clm(1,ii) = 0.
	  clp(ilevel,ii) = 0.
	end do

        !do l=1,ilevel
	!  do ii=1,3
	!    k=kn(ii)
	!    cn(l,k)    = cn(l,k)    +            cle(l,ii)
	!    clow(l,k)  = clow(l,k)  + aj4 * dt * clm(l,ii)
	!    chigh(l,k) = chigh(l,k) + aj4 * dt * clp(l,ii)
	!    cdiag(l,k) = cdiag(l,k) + aj4 * dt * clc(l,ii)
	!    cdiag(l,k) = cdiag(l,k) + aj4 * hnew(l,ii)
	!  end do
	!end do

	do ii=1,3
          do l=jlevel,ilevel
	    ccle(l,ii,ie) =            cle(l,ii)
	    cclm(l,ii,ie) = aj4 * dt * clm(l,ii)
	    cclp(l,ii,ie) = aj4 * dt * clp(l,ii)
	    cclc(l,ii,ie) = aj4 * ( dt * clc(l,ii)  &
     &				+ (1.+dt*finu(l,ii)) * hnew(l,ii) )
	  end do
          do l=ilevel+1,nlv
	    ccle(l,ii,ie) = 0.
	    cclm(l,ii,ie) = 0.
	    cclp(l,ii,ie) = 0.
	    cclc(l,ii,ie) = 0.
	  end do
	end do

	end do		! loop over ie

!----------------------------------------------------------------
! end of loop over elements
!----------------------------------------------------------------

! in cdiag, chigh, clow is matrix (implicit part)
! if explicit calculation, chigh=clow=0 and in cdiag is volume of node [m**3]
! in cnv is mass of node [kg]
! for explicit treatment, cnv/cdiag gives new concentration [kg/m**3]

!----------------------------------------------------------------
! accumulate contributions on each node
!----------------------------------------------------------------

! here we could make just one loop over nkn
! this would include this loop, the boundary conditions
! and then  the solution of the vertical system
! in this case we would only need the following arrays:
! cn(l),clow(l),chigh(l),cdiag(l) (one dimensional arrays over the vertical)

	do ie_mpi=1,nel
	  ie = ip_sort_elem(ie_mpi)
	  ilevel = ilhv(ie)
          jlevel = jlhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    do l=jlevel,ilevel
	      cn(l,k)    = cn(l,k)    + ccle(l,ii,ie)
	      clow(l,k)  = clow(l,k)  + cclm(l,ii,ie)
	      chigh(l,k) = chigh(l,k) + cclp(l,ii,ie)
	      cdiag(l,k) = cdiag(l,k) + cclc(l,ii,ie)
	    end do
	  end do
	end do
	  
	!do k=1,nkn
	!  call get_elems_around(k,maxlnk,n,elems)
	!  ilevel = ilhkv(k)
	!  do i=1,n
	!    ie = elems(i)
	!    ii = ithis(k,ie)
	!    if( ii == 0 .or. nen3v(ii,ie) /= k ) then
	!      stop 'error stop: cannot find ii...'
	!    end if
	!    do l=1,ilevel
	!      cn(l,k)    = cn(l,k)    + ccle(l,ii,ie)
	!      clow(l,k)  = clow(l,k)  + cclm(l,ii,ie)
	!      chigh(l,k) = chigh(l,k) + cclp(l,ii,ie)
	!      cdiag(l,k) = cdiag(l,k) + cclc(l,ii,ie)
	!    end do
	!  end do
	!end do

        !call shympi_comment('shympi_elem: exchange scalar')
	if( shympi_partition_on_elements() ) then
          call shympi_exchange_and_sum_3d_nodes(cn)
          call shympi_exchange_and_sum_3d_nodes(cdiag)
          call shympi_exchange_and_sum_3d_nodes(clow)
          call shympi_exchange_and_sum_3d_nodes(chigh)
	end if

!----------------------------------------------------------------
! integrate boundary conditions
!----------------------------------------------------------------

       ntot = nkn
       if( shympi_partition_on_nodes() ) ntot = nkn_unique

! in case of negative flux (qflux<0) must check if node is OBC (BUG_2010_01)

	do k=1,ntot
	  ilevel = ilhkv(k)
	  jlevel = jlhkv(k)
	  do l=jlevel,ilevel
            !mflux = cbound(l,k)		!mass flux has been passed
	    cconz = cbound(l,k)		!concentration has been passed
	    qflux = mfluxv(l,k)
	    if( qflux .lt. 0. .and. is_boundary(k) ) cconz = cn1(l,k)
	    mflux = qflux * cconz

            cn(l,k) = cn(l,k) + dt * mflux	!explicit treatment

	    loading = rload*load(l,k)
            if( loading .eq. 0. ) then
	      !nothing
	    else if( loading .gt. 0. ) then    		!treat explicit
              cn(l,k) = cn(l,k) + dt * loading
            else !if( loading .lt. 0. ) then		!treat quasi implicit
	      if( cn1(l,k) > 0. ) then
                cdiag(l,k) = cdiag(l,k) - dt * loading/cn1(l,k)
	      end if
            end if
	  end do
	end do

!----------------------------------------------------------------
! compute concentration for each node (solve system)
!----------------------------------------------------------------

	if( ( aa .eq. 0. .and. ad .eq. 0. ) .or. ( nlv .eq. 1 ) ) then

	  if( nlv .gt. 1 ) then
	    write(6,*) 'conz: computing explicitly ',nlv
	  end if

	  do k=1,ntot
	   ilevel = ilhkv(k)
	   jlevel = jlhkv(k)
	   do l=jlevel,ilevel
	    if(cdiag(l,k).ne.0.) then
	      cn(l,k)=cn(l,k)/cdiag(l,k)
	    end if
	   end do
	  end do

	else

	do k=1,ntot
	  ilevel = ilhkv(k)
	  jlevel = jlhkv(k)
	  aux=1./cdiag(jlevel,k)
	  chigh(jlevel,k)=chigh(jlevel,k)*aux
	  cn(jlevel,k)=cn(jlevel,k)*aux
	  do l=jlevel+1,ilevel
	    aux=1./(cdiag(l,k)-clow(l,k)*chigh(l-1,k))
	    chigh(l,k)=chigh(l,k)*aux
	    cn(l,k)=(cn(l,k)-clow(l,k)*cn(l-1,k))*aux
	  end do
	  lstart = ilevel-1
	  do l=lstart,jlevel,-1	!$$LEV0 bug 14.08.1998 -> ran to 0
	    cn(l,k)=cn(l,k)-cn(l+1,k)*chigh(l,k)
	  end do
	end do

	end if

	do k=1,ntot		!DPGGU
          do l=1,nlv
	    cn1(l,k)=cn(l,k)
	  end do
	end do

	if( kdebug > 0 ) then
	  write(iunit,*) 'end --------',ipext(kdebug),ilhkv(kdebug)
	  write(iunit,*) cn(1,kdebug)
	end if
	
!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!*****************************************************************

        subroutine conzstab( &
     &			ddt &
     &			,robs,rtauv,wsink,wsinkv &
     &                  ,rkpar,difhv,difv &
     &			,difmol,azpar &
     &			,adpar,aapar &
     &                  ,sindex &
     &			,istot,isact &
     &			,nlvddi,nlev)

! checks stability
!
! cn     new concentration
! co     old concentration
! caux   aux vector
! clow	 lower diagonal of vertical system
! chig	 upper diagonal of vertical system
! ddt    time step
! robs	 factor for nudging
! rtauv	 variable relaxation coefficient
! rkpar  horizontal turbulent diffusivity
! difhv  horizontal turbulent diffusivity (variable between elements)
! difv   vertical turbulent diffusivity
! difmol vertical molecular diffusivity
! azpar  time weighting parameter
! adpar  time weighting parameter for vertical diffusion (ad)
! aapar  time weighting parameter for vertical advection (aa)
! sindex stability index
! istot	 total inter time steps
! isact	 actual inter time step
! nlvddi	 dimension in z direction
! nlv	 actual needed levels
!
! solution of purely diffusional part :
!
! dC/dt = a*laplace(C)    with    c(x,0+)=delta(x)
!
! C(x,t) =  (4*pi*a*t)**(-n/2) * exp( -|x|**2/(4*a*t) )
!
! for n-dimensions and
!
! C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )
!
! for 1 dimension
!
! the solution is normalized, i.e.  int(C(x,t)dx) = 1 over the whole area
!
! DPGGU -> introduced double precision to stabilize solution

	use mod_bound_geom
	use mod_depth
	use mod_layer_thickness
	use mod_diff_aux
	use mod_bound_dynamic
	use mod_area
	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin
	use shympi
	use mkonst

	implicit none
!
! arguments
	integer nlvddi,nlev
        !real cn1(nlvddi,nkn),co1(nlvddi,nkn)		!DPGGU
        real difv(0:nlvddi,nkn)
        real difhv(nlvddi,nel)
	real difmol
        real ddt,rkpar,azpar,adpar,aapar			!$$azpar
	real robs,wsink
        real rtauv(nlvddi,nkn)
	real wsinkv(0:nlvddi,nkn)
	integer istot,isact
! common
! local
	logical bdebug,bdebug1
	integer k,ie,ii,l,iii,id,ie_mpi
	integer lstart
	integer ilevel,jlevel
	logical berror
	integer kn(3)
        real sindex,rstol,raux
	double precision us,vs
	double precision az,azt
	double precision aa,aat,ad,adt
	double precision aj,rk3,rv,aj4
	double precision hmed,hmbot,hmtop
	double precision rvptop,rvpbot
	double precision dt,w,aux
	double precision aux1,aux2,aux3,aux4,aux5
	double precision b(3),c(3),f(3)
	double precision rso,rsn,rsot,rsnt,rstot
	double precision hn,ho
        double precision wdiff(3)
	double precision dtime

	character*20 aline

	double precision difabs,difrel,volold,volnew,flxin,flxtot,diff
	double precision stabind,stabadv,stabdiff,stabvert,stabpoint
	double precision voltot

!------------------------------------------------------------
! big arrays
!------------------------------------------------------------
	double precision, allocatable :: cn(:,:)
	double precision, allocatable :: co(:,:)
	double precision, allocatable :: cdiag(:,:)
	double precision, allocatable :: clow(:,:)
	double precision, allocatable :: chigh(:,:)
        real, allocatable :: cwrite(:,:)
        real, allocatable :: c2write(:)
        real, allocatable :: saux(:,:)
!------------------------------------------------------------
! end of big arrays
!------------------------------------------------------------

	double precision cexpl
	double precision cadv,cviadv,cvoadv,chdiff,cvdiff
	double precision ciadv,coadv,chadv
	double precision fwin(3),fwout(3),fd(3)
	double precision fl(3)
! local (new)
	double precision wl(0:nlvddi+1,3)
!
	double precision hdv(0:nlvddi+1)
	double precision haver(0:nlvddi+1)
	double precision hnew(0:nlvddi+1,3)
	double precision hold(0:nlvddi+1,3)
	double precision present(0:nlvddi+1)

        integer kstab
	real dtorig

	integer, save :: icall = 0
	integer, save :: iuinfo = 0
	double precision, save :: da_out(4) = 0

! functions
	logical is_zeta_bound,openmp_in_parallel,openmp_is_master
	logical has_output_d,next_output_d
	real getpar

        if(nlv.ne.nlev) stop 'error stop conzstab: nlv/=nlev'

!-----------------------------------------------------------------
! allocation
!-----------------------------------------------------------------

	allocate(cn(nlvddi,nkn),co(nlvddi,nkn),cdiag(nlvddi,nkn))
	allocate(clow(nlvddi,nkn),chigh(nlvddi,nkn))
	allocate(cwrite(nlvddi,nkn),saux(nlvddi,nkn))
	allocate(c2write(nkn))

!-----------------------------------------------------------------
! global initialization
!-----------------------------------------------------------------

	if( icall == 0 ) then
	 if( openmp_is_master() ) then
          call init_output_d('itmstb','idtstb',da_out)
	  id = 0
          if( has_output_d(da_out) ) then
            call shyfem_init_scalar_file('stb',1,.true.,id)	!1 var, 2d
            da_out(4) = id
          end if
	  icall = 1
	  if( id > 0 ) call info_output_d('conz_stab',da_out)
	 end if
	 call getinfo(iuinfo)  !unit number of info file
	end if

!-----------------------------------------------------------------
! initialization
!-----------------------------------------------------------------

        bdebug1 = .true.
        bdebug1 = .false.
	bdebug=.false.
	berror=.false.

        if( bdebug1 ) then
                write(6,*) 'debug parameters in conz3d'
		write(6,*) ddt,rkpar,difmol,azpar,adpar,aapar
                write(6,*) istot,isact,nlvddi,nlv
                write(6,*) nkn,nel
        end if

	az=azpar		!$$azpar
	azt=1.-az
	ad=adpar
	adt=1.-ad
	aa=aapar
	aat=1.-aa

	!if( aa .ne. 0. .and. nlv .gt. 1 ) then
	!  write(6,*) 'aapar = ',aapar
	!  write(6,*) 'Cannot use implicit vertical advection.'
	!  write(6,*) 'This might be resolved in a future version.'
	!  write(6,*) 'Please set aapar = 0 in the STR file.'
	!  stop 'error stop conzstab: implicit vertical advection'
	!end if

!	-----------------------------------------------------------------
!	 fractional time step
!	-----------------------------------------------------------------

	rstot=istot		!$$istot
	rso=(isact-1)/rstot
	rsn=(isact)/rstot
	rsot=1.-rso
	rsnt=1.-rsn

	dt=ddt/rstot

!	-----------------------------------------------------------------
!	 initialize global arrays for accumulation of implicit terms
!	-----------------------------------------------------------------

	do k=1,nkn
          do l=1,nlv
	    !co(l,k)=cn1(l,k)	!DPGGU	!not used for stability
            cn(l,k)=0.          !Malta
            co(l,k)=0.
	    if( mfluxv(l,k) .gt. 0. ) co(l,k) = mfluxv(l,k)	!point sources
            cdiag(l,k)=0.
            clow(l,k)=0.
            chigh(l,k)=0.
            cwrite(l,k)=0.
	    saux(l,k)=0.
          end do
	end do

!	-----------------------------------------------------------------
!	these are aux arrays (bigger than needed) to avoid checking for
!	what layer we are in -> we never get out of bounds
!	-----------------------------------------------------------------

        do l=0,nlv+1
	  hdv(l) = 0.		!layer thickness
          haver(l) = 0.
	  present(l) = 0.	!1. if layer is present
	  do ii=1,3
	    hnew(l,ii) = 0.	!as hreal but with zeta_new
	    hold(l,ii) = 0.	!as hreal but with zeta_old
	    !cl(l,ii) = 0.	!concentration in layer
	    wl(l,ii) = 0.	!vertical velocity
	  end do
	end do

!-----------------------------------------------------------------
! loop over elements
!-----------------------------------------------------------------

        do ie_mpi=1,nel

	ie = ip_sort_elem(ie_mpi)

	do ii=1,3
          k=nen3v(ii,ie)
	  kn(ii)=k
	  b(ii)=ev(ii+3,ie)
	  c(ii)=ev(ii+6,ie)
	end do

	aj=ev(10,ie)    !area of triangle / 12
	aj4=4.*aj
        ilevel=ilhv(ie)
        jlevel=jlhv(ie)

! set up vectors for use in assembling contributions

        do l=jlevel,ilevel
	  hdv(l) = hdeov(l,ie)		!use old time step -> FIXME
          haver(l) = 0.5 * ( hdeov(l,ie) + hdenv(l,ie) )
	  present(l) = 1.
	  do ii=1,3
	    k=kn(ii)
	    hn = hdknv(l,k)		! there are never more layers in ie
	    ho = hdkov(l,k)		! ... than in k
	    hold(l,ii) = rso * hn + rsot * ho
	    hnew(l,ii) = rsn * hn + rsnt * ho
	    !cl(l,ii) = co(l,k)
	    wl(l,ii) = wlnv(l,k) - wsink * wsinkv(l,k)
	  end do
	end do

	do l=ilevel+1,nlv
	  present(l) = 0.
	end do

!	set vertical velocities in surface and bottom layer
!
!	we do not set wl(0,ii) because otherwise we loose concentration
!	through surface
!
!	we set wl(ilevel,ii) to 0 because we are on the bottom
!	and there should be no contribution from this element
!	to the vertical velocity

	do ii=1,3
	  wl(ilevel,ii) = 0.
	end do

!-----------------------------------------------------------------
! loop over levels
!-----------------------------------------------------------------

        do l=jlevel,ilevel

        us=az*utlnv(l,ie)+azt*utlov(l,ie)             !$$azpar
        vs=az*vtlnv(l,ie)+azt*vtlov(l,ie)

        rk3 = 3. * rkpar * difhv(l,ie)

	do ii=1,3

	  k=kn(ii)

!	  --------------------------------------------------------
!	  new weights for horizontal diffusion
!	  --------------------------------------------------------

          wdiff(ii) = wdifhv(ii,ii,ie)

!	  --------------------------------------------------------
!	  contributions from horizontal advection
!	  --------------------------------------------------------

	  f(ii)=us*b(ii)+vs*c(ii)	!$$azpar - positive if flux into node

!	  --------------------------------------------------------
!	  contributions from vertical advection
!	  --------------------------------------------------------
!
!	  in fw(ii) is explicit contribution
!	  now using fwin and fwout
!	  fw is positive if out of layer
!
!	  if we are in last layer, wl(l,ii) is zero
!	  if we are in first layer, wl(l-1,ii) is zero (see above)
!	  wl already contains sinking velocity

	  fwin(ii) = 0.
	  fwout(ii) = 0.

	  w = wl(l-1,ii)		!top of layer
	  if( w .gt. 0. ) then          !out
	    fwout(ii) = fwout(ii) + aat*w
	  else
	    fwin(ii) = fwin(ii) - aat*w
	  end if

	  w = wl(l,ii)			!bottom of layer
	  if( w .gt. 0. ) then
	    fwin(ii) = fwin(ii) + aat*w
	  else				!out
	    fwout(ii) = fwout(ii) - aat*w
	  end if

!	  --------------------------------------------------------
!	  contributions from vertical diffusion
!	  --------------------------------------------------------
!
!	  in fd(ii) is explicit contribution
!	  the sign is for the term on the left side, therefore
!	  fd(ii) must be subtracted from the right side
!
!	  maybe we should use real layer thickness, or even the
!	  time dependent layer thickness

	  rvptop = difv(l-1,k) + difmol
	  rvpbot = difv(l,k) + difmol
	  hmtop = 2. * rvptop * present(l-1) / (hdv(l-1)+hdv(l))
	  hmbot = 2. * rvpbot * present(l+1) / (hdv(l)+hdv(l+1))

          fd(ii) = adt * ( hmtop + hmbot )

	end do

!	--------------------------------------------------------
!	sum explicit contributions
!	--------------------------------------------------------
!
!	chigh contains explicit advection (horizontal and vertical)
!	clow contains explicit horizontal diffusion
!	cn contains flux due to explicit vertical diffusion
!	co contains flux due to point sources and nudging

	coadv = 0.
	ciadv = 0.
	do ii=1,3
	  k=kn(ii)
          hmed = hold(l,ii)                      !new ggu   !HACK
          chdiff = dt * aj4 * rk3 * hmed * wdiff(ii)	!bug fix 12.2.2010
	  cvdiff = dt * aj4 * fd(ii)
          chadv = dt * aj4 * 3. * f(ii)
	  cviadv = dt * aj4 * fwin(ii)
	  cvoadv = dt * aj4 * fwout(ii)
          if( chadv < 0. ) then             !flux out of node
	    coadv = cvoadv - chadv
	  else
	    ciadv = cviadv + chadv
          end if
	  cadv = max(coadv,ciadv)		!here we use the absolute value
	  chigh(l,k) = chigh(l,k) + cadv
	  clow(l,k) = clow(l,k) + chdiff
          cn(l,k) = cn(l,k) + cvdiff
          co(l,k) = co(l,k) + dt * aj4 * hmed * robs * rtauv(l,k) !nudging
	end do

	end do		! loop over l

!	--------------------------------------------------------
!	sum volumes
!	--------------------------------------------------------
!
!	cdiag contains volume of finite node

        do l=jlevel,ilevel
	  do ii=1,3
	    k=kn(ii)
            hmed = min(hold(l,ii),hnew(l,ii))
	    cdiag(l,k) = cdiag(l,k) + aj4 * hmed
	  end do
	end do

	end do		! loop over ie

!-----------------------------------------------------------------
! compute stability
!
! cdiag		volume of cell
! chigh		flux due to advection
! clow		flux due to horizontal diffusion
! cn		flux due to vertical diffusion (explicit)
! co		flux due to point sources and nudging
!-----------------------------------------------------------------

        stabind = 0.		!total max stability index
        stabadv = 0.		!advective max stability index
        stabdiff = 0.		!diffusive max stability index
        stabvert = 0.		!vertical max stability index
        stabpoint = 0.		!point source max stability index
        kstab = 0		!node with highest stabind

	do k=1,nkn_unique
	  bdebug1 = k .eq. -1
	  ilevel = ilhkv(k)
	  jlevel = jlhkv(k)
          if( is_zeta_bound(k) ) cycle
	  do l=jlevel,ilevel
            voltot = cdiag(l,k)
            flxtot = chigh(l,k) + clow(l,k) + cn(l,k) + co(l,k)
	    if( bdebug1 ) write(99,*) k,l,voltot,flxtot
            if( voltot .gt. 0. ) then
                  aux1 = flxtot / voltot
                  if( aux1 .gt. stabind ) kstab = k
                  stabind = max(stabind,aux1)
		  cwrite(l,k) = aux1		!save for write
                  aux2 = chigh(l,k) / voltot
                  stabadv = max(stabadv,aux2)
		  saux(l,k) = aux2		!for adv. stab.
                  aux3 = clow(l,k) / voltot
                  stabdiff = max(stabdiff,aux3)
                  aux4 = cn(l,k) / voltot
                  stabvert = max(stabvert,aux4)
                  aux5 = co(l,k) / voltot
                  stabpoint = max(stabpoint,aux5)
	          if( bdebug1 ) write(99,*) aux1,aux2,aux3,aux4,aux5
            end if
	  end do
	end do

	raux = stabind			!FIXME - SHYFEM_FIXME
        stabind = shympi_max(raux)

!-----------------------------------------------------------------
! in stabind is stability index (advection and diffusion)
! in cdiag is the volume of the node
! in cwrite is the value of the stability index for each node
!
! istot  is saved and returned from subroutine (number of iterations)
! sindex is saved and returned from subroutine (stability index)
!-----------------------------------------------------------------

	rstol = getpar('rstol')
        istot = 1 + stabind / rstol
        sindex = stabind / rstol

	!write(iuinfo,*) 'stability_scalar: ',istot,sindex,stabind

	!if( .not. openmp_in_parallel() ) then
	if( openmp_is_master() ) then
          if( next_output_d(da_out) ) then
            id = nint(da_out(4))
	    call get_act_dtime(dtime)
	    call get_act_timeline(aline)
	    write(6,*) 'writing STB at ',aline,' (more info in fort.197)'
	    forall(k=1:nkn) c2write(k) = maxval(cwrite(:,k))
            call shy_write_scalar_record2d(id,dtime,75,c2write)
	    write(197,1010) aline,kstab &
     &			,stabind,stabadv,stabdiff,stabvert,stabpoint
	  end if
	end if

!-----------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------

	deallocate(cn,co,cdiag)
	deallocate(clow,chigh)
	deallocate(cwrite,saux)

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
 1010	format(a20,i10,5f10.3)
	end

!*****************************************************************
!*****************************************************************
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine massconc(mode,cn,nlvddi,mass)

! computes total mass of conc

	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

	implicit none

! arguments
	integer mode
	integer nlvddi
	real cn(nlvddi,nkn)
	real mass
! local
	integer k,l,lmax,lmin
        double precision vol
	double precision sum,masstot
	real volnode

!SHYMPI_ELEM - should be total nodes to use - FIXME shympi

        masstot = 0.

        do k=1,nkn
	  lmax = ilhkv(k)
	  lmin = jlhkv(k)
          sum = 0.
          do l=lmin,lmax
            vol = volnode(l,k,mode)
            sum = sum + cn(l,k) * vol
          end do
          masstot = masstot + sum
        end do

	mass = masstot

        mass = shympi_sum(mass)
        !call shympi_comment('massconc: shympi_sum(masstot)')
	!shympi on elems FIXME

	end

!**************************************************************

        subroutine check_scal_bounds(cnv,cmin,cmax,eps,bstop)

! checks if scalar is out of bounds

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        real cnv(nlvdi,nkn)
        real cmin,cmax
	real eps
	logical bstop		!stop simulation if true


	logical berror
        integer k,l,lmax,kext
        real cc,cn,cx,cd

        integer ipext

	berror = .false.
	cd = cmax - cmin
	if( cd .le. 0. ) return

	cn = cmin - eps*cd
	cx = cmax + eps*cd

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cc = cnv(l,k)
            if( cc .lt. cn .or. cc .gt. cx ) then
                kext = ipext(k)
                write(6,*) 'scalar out of bounds: ',k,kext,l,lmax,cc
		berror = .true.
            end if
          end do
        end do

	if( bstop .and. berror ) then
	  stop 'error stop check_scal_bounds: out of bounds'
	end if

        end

!*****************************************************************

	subroutine assert_min_max_property(cnv,cov,sbconz,rmin,rmax,eps)

! checks min/max property

	use mod_bound_geom
	use mod_bound_dynamic
	use levels
	use basin

	implicit none

        real cnv(nlvdi,nkn)			!new concentration
        real cov(nlvdi,nkn)			!old concentration
        real sbconz(nlvdi,nkn)		!conz boundary conditions
	real rmin(nlvdi,nkn)		!aux arrray to contain min
	real rmax(nlvdi,nkn)		!aux arrray to contain max
	real eps

	logical bwrite,bstop
	integer k,ie,l,ii,lmax,lmin,ierr
	integer levdbg
	real amin,amax,c,qflux,dmax
	real drmax,diff
	real dt
	character*20 aline

	integer ipext
	logical is_zeta_bound
	real getpar

	bwrite = .true.		! write every violation
	bwrite = .false.		! write every violation
	bstop = .false.		! stop after error

	levdbg = nint(getpar('levdbg'))

!---------------------------------------------------------------
! vertical contribution (normally implicit -> whole column)
!---------------------------------------------------------------

	do k=1,nkn
	  lmax = ilhkv(k)
	  lmin = jlhkv(k)
	  amin = +1.e+30
	  amax = -1.e+30
	  do l=lmin,lmax
	    amin = min(amin,cov(l,k))
	    amax = max(amax,cov(l,k))
	  end do
	  do l=lmin,lmax
	    rmin(l,k) = amin
	    rmax(l,k) = amax
	  end do
	end do

!---------------------------------------------------------------
! point sources
!---------------------------------------------------------------

	do k=1,nkn
	  lmax = ilhkv(k)
	  lmin = jlhkv(k)
	  do l=lmin,lmax
	    qflux = mfluxv(l,k)
	    if( qflux .gt. 0. ) then
	      c = sbconz(l,k)
	      rmin(l,k) = min(rmin(l,k),c)
	      rmax(l,k) = max(rmax(l,k),c)
	    end if
	  end do
	end do

!---------------------------------------------------------------
! horizontal contribution
!---------------------------------------------------------------

	do ie=1,nel
	  lmax = ilhv(ie)
	  lmin = jlhv(ie)
	  do l=lmin,lmax
	    amin = +1.e+30
	    amax = -1.e+30
	    do ii=1,3
	      k = nen3v(ii,ie)
	      amin = min(amin,cov(l,k))
	      amax = max(amax,cov(l,k))
	    end do
	    do ii=1,3
	      k = nen3v(ii,ie)
	      rmin(l,k) = min(rmin(l,k),amin)
	      rmax(l,k) = max(rmax(l,k),amax)
	    end do
	  end do
	end do

!---------------------------------------------------------------
! check with new concentration
!---------------------------------------------------------------

	ierr = 0
	dmax = 0.
	drmax = 0.

	call get_act_timeline(aline)

	do k=1,nkn
	 !if( .not. is_external_boundary(k) ) then	!might be relaxed
	 if( .not. is_zeta_bound(k) ) then	!might be relaxed
	  lmax = ilhkv(k)
	  lmin = jlhkv(k)
	  do l=lmin,lmax
	    c = cnv(l,k)
	    !rm1 = rmin(l,k)
	    !rm2 = rmax(l,k)
	    if( c .lt. rmin(l,k) .or. c .gt. rmax(l,k) ) then
	      amin = rmin(l,k)
	      amax = rmax(l,k)
	      diff = max(rmin(l,k)-c,c-rmax(l,k))
	      dmax = max(dmax,diff)
	      drmax = max(drmax,diff*(amax-amin))
	      if( bwrite ) then
	        write(6,*) 'min/max property violated: '
	        write(6,*) '   ',l,k,ipext(k)
	        write(6,*) '   ',c,amin,amax
	      end if
	      ierr = ierr + 1
	    end if
	  end do
	 end if
	end do

!---------------------------------------------------------------
! check error condition
!---------------------------------------------------------------

	if( ierr .gt. 0 ) then
	  if( drmax .gt. eps .and. dmax .gt. eps ) then
	    if( levdbg .ge. 3 ) then
	      write(6,*) 'min/max error: ',aline,ierr,dmax,drmax
	    end if
	  end if
	  !write(94,*) 'min/max error violated: ',aline,ierr,dmax,drmax
	  if( bstop ) then
	    stop 'error stop assert_min_max_property: violation'
	  end if
	end if

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine stb_histo(it,nlvddi,nkn,ilhkv,cwrite)

! writes histogram info about stability index

	use mod_histo

        implicit none

        integer it
        integer nlvddi,nkn
        integer ilhkv(nkn)
        real cwrite(nlvddi,nkn)

        integer k,l,lmax

        integer, parameter :: nbin = 11
        integer ic(nbin+1)
        real, save :: bins(nbin) = (/1.,2.,5.,10.,15.,20.,30.,40.,50.,75.,100./)

        call histo_init(nbin,bins)

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            call histo_insert(cwrite(l,k))
          end do
        end do

        call histo_final(ic)

        write(98,*) it,ic

        end

!*****************************************************************

	subroutine limit_scalar(what,dtstep,cnv)

	use levels
	use basin

	implicit none

	character*(*) what
	real dtstep
        real cnv(nlvdi,nkn)			!new concentration

	integer, save :: icall = 0
	logical, save :: btlimit0,btlimit1
	logical, save :: bslimit0,bslimit1
	logical, save :: bclimit0,bclimit1
	real, save :: tlimit0,slimit0,climit0
	real, save :: tlimit1,slimit1,climit1
	real, parameter :: flag = -999.
	integer, parameter :: iu = 0		!set this to iu/=0 for write
	!integer, parameter :: iu = 777		!set this to iu/=0 for write
	logical, save :: bwrite = .false.
	logical, allocatable, save :: mask(:,:)
	logical blimit0,blimit1
	real limit0,limit1
	integer ic
	real thresh
	real getpar

	if( icall == 0 ) then
	  icall = 1

	  bwrite = ( iu > 0 )

	  allocate(mask(nlvdi,nkn))

	  btlimit0 = .false.
	  btlimit1 = .false.
	  tlimit0 = getpar('tlimit0')
	  tlimit1 = getpar('tlimit1')

	  bslimit0 = .false.
	  bslimit1 = .false.
	  slimit0 = getpar('slimit0')
	  slimit1 = getpar('slimit1')

	  bclimit0 = .false.
	  bclimit1 = .false.
	  climit0 = getpar('climit0')
	  climit1 = getpar('climit1')

	  if( tlimit0 /= flag ) btlimit0 = .true.
	  if( tlimit1 /= flag ) btlimit1 = .true.
	  if( btlimit0 .and. btlimit1 .and. tlimit0 > tlimit1 ) goto 99
	  
	  if( slimit0 /= flag ) bslimit0 = .true.
	  if( slimit1 /= flag ) bslimit1 = .true.
	  if( bslimit0 .and. bslimit1 .and. slimit0 > slimit1 ) goto 99
	  
	  if( climit0 /= flag ) bclimit0 = .true.
	  if( climit1 /= flag ) bclimit1 = .true.
	  if( bclimit0 .and. bclimit1 .and. climit0 > climit1 ) goto 99
	  
	  write(6,*) 'limiting scalars has been set up'
	  if( btlimit0 ) write(6,*) 'limiting min temp: ',tlimit0
	  if( btlimit1 ) write(6,*) 'limiting max temp: ',tlimit1
	  if( bslimit0 ) write(6,*) 'limiting min salt: ',slimit0
	  if( bslimit1 ) write(6,*) 'limiting max salt: ',slimit1
	  if( bclimit0 ) write(6,*) 'limiting min conz: ',climit0
	  if( bclimit1 ) write(6,*) 'limiting max conz: ',climit1
	end if

	if( what == 'temp' ) then
	  !write(6,*) 'limiting temp: ',tlimit0,tlimit1
	  blimit0 = btlimit0
	  blimit1 = btlimit1
	  limit0 = tlimit0
	  limit1 = tlimit1
	else if( what == 'salt' ) then
	  !write(6,*) 'limiting salt: ',slimit0,slimit1
	  blimit0 = bslimit0
	  blimit1 = bslimit1
	  limit0 = slimit0
	  limit1 = slimit1
	else if( what == 'conz' ) then
	  !write(6,*) 'limiting conz: ',climit0,climit1
	  blimit0 = bclimit0
	  blimit1 = bclimit1
	  limit0 = climit0
	  limit1 = climit1
	else
	  blimit0 = .false.
	  blimit1 = .false.
	end if

	if( blimit0 ) then
	  if( bwrite ) then
	    mask = .false.
	    where( cnv < limit0 ) mask = .true.
	    call write_out_of_limit(iu,dtstep,what,' < ',limit0,mask,cnv)
	  end if
	  ic = 0
	  !ic=count(cnv<limit0)
	  if( ic > 0 ) write(6,*) 'ggu0 before ',what,ic
	  where( cnv < limit0 ) cnv = limit0
	  !ic=count(cnv<limit0)
	  if( ic > 0 )write(6,*) 'ggu0 after ',what,ic
	end if
	if( blimit1 ) then
	  if( bwrite ) then
	    mask = .false.
	    where( cnv > limit1 ) mask = .true.
	    call write_out_of_limit(iu,dtstep,what,' > ',limit1,mask,cnv)
	  end if
	  ic = 0
	  !ic=count(cnv>limit1)
	  if( ic > 0 ) write(6,*) 'ggu1 before ',what,ic
	  where( cnv > limit1 ) cnv = limit1
	  !ic=count(cnv>limit1)
	  if( ic > 0 ) write(6,*) 'ggu1 after ',what,ic
	end if
	
	return
   99	continue
	write(6,*) 'error setting limiter...'
	write(6,*) 'either both limiters are set or are left as flag'
	write(6,*) 'tlimit: ',tlimit0,tlimit1
	write(6,*) 'slimit: ',slimit0,slimit1
	write(6,*) 'tlimit: ',climit0,climit1
	end

!*****************************************************************

	subroutine write_out_of_limit(iu,dtstep,what,gtlt,limit,mask,cnv)

	use levels
	use basin

	implicit none

	integer iu
	real dtstep
	character*(*) what,gtlt
	real limit
	logical mask(nlvdi,nkn)
	real cnv(nlvdi,nkn)

	integer k,l
	double precision dtime
	character*20 aline

	if( count( mask ) == 0 ) return

	call get_act_dtime(dtime)
	call get_act_timeline(aline)

	dtime = dtime + dtstep
	write(iu,*) 'time ',dtime,' (',aline,') ' &
     &			,trim(what),trim(gtlt),limit
	write(iu,*) '   internal node       layer        value'

	do k=1,nkn
	  do l=1,nlvdi
	    if( mask(l,k) ) then
	      write(iu,*) '    ',k,l,cnv(l,k)
	    end if    
	  end do
	end do

	end

!*****************************************************************

