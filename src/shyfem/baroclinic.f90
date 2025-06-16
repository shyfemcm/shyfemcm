
!--------------------------------------------------------------------------
!
!    Copyright (C) 1994-1995,1998,2001-2005,2007-2020  Georg Umgiesser
!    Copyright (C) 2007,2009,2012  Debora Bellafiore
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

! dd: newbcl.f,v 1.37 2010-03-08 17:46:45 georg Exp $
!
! baroclinic routines
!
! contents :
!
! subroutine barocl(mode)		amministrates the baroclinic time step
! subroutine rhoset_shell(iterwant)	sets rho iterating to real solution
! subroutine rhoset(resid)		computes rhov
! subroutine convectivecorr             convective adjustment
! subroutine getts(l,k,t,s)             accessor routine to get T/S
!
! revision log :
!
! 09.01.1994	ggu	(from scratch)
! 30.08.1995	ggu	$$AUST - austausch coefficient introduced
! 11.10.1995	ggu	$$BCLBND - boundary condition for barocliic runs
! 19.08.1998	ggu	call to barcfi changed
! 20.08.1998	ggu	can initialize S/T from file
! 24.08.1998	ggu	levdbg used for debug
! 26.08.1998	ggu	init, bnd and file routines substituted with con..
! 30.01.2001	ggu	eliminated compile directives
! 05.12.2001	ggu	horizontal diffusion variable, limit diffusion coef.
! 05.12.2001	ggu	compute istot, more debug info
! 11.10.2002	ggu	diffset introduced, shpar = thpar
! 10.08.2003	ggu	qfluxr eliminated (now in subn11.f)
! 10.08.2003	ggu	rhov and bpresv are initialized here
! 04.03.2004	ggu	in init for T/S pass number of vars (inicfil)
! 15.03.2004	ggu	general clean-up, bclint() deleted, new scal3sh
! 17.01.2005	ggu	new difhv implemented
! 15.03.2005	ggu	new diagnostic routines implemented (diagnostic)
! 15.03.2005	ggu	new 3d boundary conditions implemented
! 05.04.2005	ggu	some changes in routine diagnostic
! 07.11.2005	ggu	sinking velocity wsink introduced in call to scal3sh
! 08.06.2007	ggu&dbf	restructured for new baroclinic version
! 04.10.2007	ggu	bug fix -> call qflux3d with dt real
! 17.03.2008	ggu	new open boundary routines introduced
! 08.04.2008	ggu	treatment of boundaries slightly changed
! 22.04.2008	ggu	advection parallelized, no saux1v...
! 23.04.2008	ggu	call to bnds_set_def() changed
! 12.06.2008	ggu	s/tdifhv deleted
! 09.10.2008	ggu	new call to confop
! 12.11.2008	ggu	new initialization, check_layers, initial nos file
! 13.01.2009	ggu&dbf	changes in reading file in ts_next_record()
! 13.10.2009	ggu	in rhoset bug computing pres
! 13.11.2009	ggu	only initialize T/S if no restart, new rhoset_shell
! 19.01.2010	ggu	different call to has_restart() 
! 23.03.2010	ggu	changed v6.1.1
! 16.12.2010	ggu	sigma layers introduced (maybe not finished)
! 26.01.2011	ggu	read in obs for t/s (tobsv,sobsv)
! 28.01.2011	ggu	parameters changed in call to ts_nudge()
! 17.02.2011	ggu	changed VERS_6_1_18
! 04.03.2011	ggu	better error message for rhoset_shell
! 23.03.2011	ggu	changed VERS_6_1_21
! 31.03.2011	ggu	only write temp/salt if computed
! 14.04.2011	ggu	changed VERS_6_1_22
! 31.05.2011	ggu	changed VERS_6_1_23
! 15.07.2011	ggu	changed VERS_6_1_28
! 04.11.2011	ggu	adapted for hybrid coordinates
! 07.11.2011	ggu	hybrid changed to resemble code in newexpl.f
! 11.11.2011	ggu	restructured ts_next_record() and diagnostic()
! 22.11.2011	ggu	bug fix in ts_file_open() -> bhashl
! 02.12.2011	ggu	adapt ts_file_open() for barotropic version (ihashl)
! 09.12.2011	ggu	changed VERS_6_1_38
! 27.01.2012	dbf&ggu	changes for hybrid in ts_file_open,ts_next_record
! 10.02.2012	ggu	bug in call to ts_next_record (called with nlvddi)
! 23.02.2012	ccf	do noy check depth structure
! 09.03.2012	dbf	bug fix in ts_next_record: ilhkv was real
! 16.03.2012	ggu	changed VERS_6_1_48
! 21.03.2012	ggu	changed VERS_6_1_50
! 30.03.2012	ggu	changed VERS_6_1_51
! 21.06.2012	ggu	changed VERS_6_1_54
! 08.10.2012	ggu	changed VERS_6_1_58
! 31.10.2012	ggu	open and next_record transfered to subtsuvfile.f
! 05.11.2012	ggu	changed VERS_6_1_60
! 05.09.2013	ggu	limit salinity to [0,...]
! 12.09.2013	ggu	changed VERS_6_1_67
! 25.10.2013	ggu	changed VERS_6_1_68
! 25.03.2014	ggu	new offline
! 27.06.2014	ggu	changed VERS_6_1_78
! 07.07.2014	ggu	changed VERS_6_1_79
! 10.07.2014	ggu	only new file format allowed
! 18.07.2014	ggu	changed VERS_7_0_1
! 20.10.2014	ggu	pass ids to scal_adv()
! 05.11.2014	ggu	changed VERS_7_0_5
! 26.11.2014	ggu	changed VERS_7_0_7
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 10.02.2015	ggu	call to bnds_read_new() introduced
! 26.02.2015	ggu	changed VERS_7_1_5
! 21.05.2015	ggu	changed VERS_7_1_11
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 30.07.2015	ggu	changed VERS_7_1_83
! 18.09.2015	ggu	changed VERS_7_2_3
! 23.09.2015	ggu	changed VERS_7_2_4
! 29.09.2015	ggu	changed VERS_7_2_5
! 15.10.2015	ggu	added new calls for shy file format
! 22.10.2015	ggu	changed VERS_7_3_7
! 23.10.2015	ggu	changed VERS_7_3_9
! 26.10.2015	ggu	bug fix for parallel code (what was not set)
! 05.11.2015	ggu	changed VERS_7_3_12
! 09.11.2015	ggu	changed VERS_7_3_13
! 28.04.2016	ggu	changed VERS_7_5_9
! 07.06.2016	ggu	changed VERS_7_5_12
! 10.06.2016	ggu	not used routines deleted
! 14.06.2016	ggu	open and write of file in own subroutine
! 27.06.2016	ggu	bug fix: irho was not saved
! 11.10.2016	ggu	changed VERS_7_5_20
! 12.01.2017	ggu	changed VERS_7_5_21
! 04.11.2017	ggu	changed VERS_7_5_34
! 03.04.2018	ggu	changed VERS_7_5_43
! 19.04.2018	ggu	changed VERS_7_5_45
! 05.10.2018	ggu	new diagnostic routine ts_dia()
! 16.10.2018	ggu	changed VERS_7_5_50
! 18.12.2018	ggu	changed VERS_7_5_52
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 21.05.2019	ggu	changed VERS_7_5_62
! 14.02.2020	ggu	nudging enhanced with reading of tau values
! 05.03.2020	ggu	finished new nudging routines
! 27.03.2020	ggu	cleaned new nudging routines
! 17.09.2020    ggu     renamed sigma to sigma_stp
! 24.03.2021    ggu     more diagnostic output in ts_dia()
! 31.05.2021    ggu     changes in ts_dia(), debug section
! 09.04.2022    ggu     new iterwant in rhoset_shell for mpi
! 15.10.2022	ggu	bug fix: rhov was recomputed after restart
! 15.10.2022	ggu	bpresv deleted
! 02.04.2023	ggu	min/max of T/S computed correctly for mpi
! 09.05.2023    lrp     introduce top layer index variable
! 03.12.2024    ggu     new info_output framework
! 05.02.2025    ggu     new info_format for T/S
! 24.04.2025    ggu     new value for ibarcl: ibercl == 5
!
! notes :
!
! initialization of T/S values
!
! if restart 
!	T/S not initialized (are initialized with restart values)
! else if bdiag
!	initialized with ts_diag values
! else
!	initialized with ts_init values
! end if
!
! if bobs 
!	read ts_obs values (for nudging)
! end if
!
! normal call
!
! if bdiag
!	read ts_diag values into ts
! else if bobs
!	read ts_obs values
! end if
!
! important: bobs does not imply initialization of ts with ts_obs values
!		only if ts_init is given ts is initialized
!
!*****************************************************************

	subroutine barocl(mode)

! amministrates the baroclinic time step
!
! mode : =0 initialize  >0 normal call
!
	use mod_layer_thickness
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use mod_hydro
	use levels
	use basin
	use shympi
	use pkonst
	use mkonst
	use mod_trace_point
	use mod_info_output

	implicit none
!
! arguments
	integer mode
! common

! local
	logical debug
	logical bgdebug
        logical bstop
	logical binitial_nos
	logical boff
	integer levdbg
	integer ie,ks
	integer idtext,itmext
	integer imin,imax
	integer nintp,nvar
	integer nbc
	integer id
	real cdef(1)
	real xmin,xmax
	real salref,temref,sstrat,tstrat
	real shpar,thpar
        real s
	real dt
        real gamma,gammax
	real mass
	real wsink
	real robs
	real array(4)
	real, allocatable :: rho_aux1(:,:)
	real, allocatable :: rho_aux2(:,:)
	double precision dtime0,dtime
	integer isact,l,k,lmax
	integer kspec
	integer icrst
	integer ftype
	real stot,ttot,smin,smax,tmin,tmax,rmin,rmax
	real saver,taver,vtot
	double precision v1,v2,mm
	character*80 file
	character*20 aline
	character*4 what
	character*80 format
! functions
!	real sigma
	real getpar
	double precision dq
	integer iround
	integer nbnds
	logical rst_use_restart
	logical has_output_d,next_output_d

	integer tid
	!integer openmp_get_thread_num
	
	double precision theatold,theatnew
	double precision theatconv1,theatconv2,theatqfl1,theatqfl2
! save
	logical, save :: badvect,bobs,bbarcl,bdiag
        integer, save :: iuinfo = 0
	integer, save :: ibarcl
        integer, save :: itemp,isalt,irho
	real, save :: difmol
        double precision, save :: da_out(4)
	integer, save, allocatable :: idtemp(:),idsalt(:)
	integer, save :: icall = 0

!----------------------------------------------------------
! parameter setup and check
!----------------------------------------------------------

	if(icall.eq.-1) return

        call is_offline(2,boff)
        !if( boff ) write(6,*) 'TS reading from offline...'
        if( boff ) return

	levdbg = nint(getpar('levdbg'))
	debug = levdbg .ge. 3
        bgdebug = .false.
	binitial_nos = .true.

!----------------------------------------------------------
! initialization
!----------------------------------------------------------

	if(icall.eq.0) then	!first time

		call trace_point('initalizing barocl')
		call shympi_bdebug('initalizing barocl')

		ibarcl=iround(getpar('ibarcl'))
		if(ibarcl.le.0) icall = -1
		if(ibarcl.gt.5) then
		  write(6,*) 'Value of ibarcl not allowed: ',ibarcl
	          stop 'error stop barocl: ibarcl'
		end if
		if(icall.eq.-1) return

		bdiag = ibarcl .eq. 2			!diagnostic run
		badvect = .not. bdiag			!must advect T/S
		bobs = ibarcl .eq. 4 .or. ibarcl .eq. 5	!nudging
		bbarcl = ibarcl .ne. 3			!baroclinic run
		if( ibarcl .eq. 5 ) bbarcl = .false.

		salref=getpar('salref')
		temref=getpar('temref')
		sstrat=getpar('sstrat')
		tstrat=getpar('tstrat')
		difmol=getpar('difmol')
                itemp=iround(getpar('itemp'))
                isalt=iround(getpar('isalt'))
                irho=iround(getpar('irho'))

!		--------------------------------------------
!		initialize saltv,tempv
!		--------------------------------------------

		call get_first_dtime(dtime0)
		call trace_point('initalizing S/T')

		if( .not. rst_use_restart(3) ) then   !no restart of T/S values
		  saltv = 0.
		  tempv = 0.
		  call conini(nlvdi,saltv,salref,sstrat,hdkov)
		  call conini(nlvdi,tempv,temref,tstrat,hdkov)

		  if( bdiag ) then
		    call ts_diag(dtime0,nlvdi,nlv,nkn,tempv,saltv)
		  else
		    call ts_init(dtime0,nlvdi,nlv,nkn,tempv,saltv)
		  end if
		end if

		if( bobs ) then
	  	  call ts_nudge(dtime0,nlvdi,nlv,nkn,tobsv,sobsv &
     &						,ttauv,stauv)
		end if

!		--------------------------------------------
!		initialize observations and relaxation times
!		--------------------------------------------

		tobsv = 0.
		sobsv = 0.
		ttauv = 0.
		stauv = 0.

!		--------------------------------------------
!		initialize open boundary conditions
!		--------------------------------------------

		call trace_point('initalizing OBC of S/T')

		nbc = nbnds()
		allocate(idtemp(nbc))
		allocate(idsalt(nbc))
		idtemp = 0
		idsalt = 0

                nintp = 2
                nvar = 1
                cdef(1) = 0.
		what = 'temp'
		call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv &
     &					,cdef,idtemp)
		what = 'salt'
		call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv &
     &					,cdef,idsalt)

!		--------------------------------------------
!		initialize rhov
!		--------------------------------------------

		call trace_point('initalizing rhov')

!		rhov depends on pressure and viceversa
!		-> we iterate to the real solution

		if( .not. rst_use_restart(3) ) then   !no restart of rho values
		  rhov = 0.		!rhov is rho^prime => 0
		  !call ts_dia('init before rhoset_shell')
		  call rhoset_shell(3)
		  !call ts_dia('init after rhoset_shell')
		end if

!		--------------------------------------------
!		initialize output files
!		--------------------------------------------

		call trace_point('initalizing T/S output')

		call bcl_open_output(da_out,itemp,isalt,irho)
		call trace_point('writing T/S output')
		call bcl_write_output(dtime0,da_out,itemp,isalt,irho)

		if( shympi_is_master() ) then
                  call getinfo(iuinfo)
		end if

		call shympi_bdebug('finished initalizing barocl')
	end if

	icall=icall+1

	if(mode.eq.0) return

	call trace_point('barocl normal call')

!----------------------------------------------------------
! normal call
!----------------------------------------------------------

	!call ts_dia('begin normal call')

	call get_act_dtime(dtime)
	call get_act_timeline(aline)

	wsink = 0.
	robs = 0.
	if( bobs ) robs = 1.

	shpar=getpar('shpar')   !out of initialization because changed
	thpar=getpar('thpar')

	if( bdiag ) then
	  call ts_diag(dtime,nlvdi,nlv,nkn,tempv,saltv)
	else if( bobs ) then
	  call ts_nudge(dtime,nlvdi,nlv,nkn,tobsv,sobsv,ttauv,stauv)
	end if

!----------------------------------------------------------
! salt and temperature transport and diffusion
!----------------------------------------------------------

	if( badvect ) then

	  call openmp_get_thread_num(tid)
	  !write(6,*) 'number of thread of temp: ',tid

          if( itemp .gt. 0 ) then
		what = 'temp'
	        call bnds_read_new(what,idtemp,dtime)
	  end if
          if( isalt .gt. 0 ) then
		what = 'salt'
	        call bnds_read_new(what,idsalt,dtime)
	  end if
	  
	  !call ts_dia('before T/D')

!$OMP TASK PRIVATE(what,dtime) FIRSTPRIVATE(thpar,wsink,robs,itemp) &
!$OMP&     SHARED(idtemp,tempv,difhv,difv,difmol,tobsv,ttauv) &
!$OMP&     DEFAULT(NONE) &
!$OMP&     IF(itemp > 0)

          if( itemp .gt. 0 ) then
		what = 'temp'
                call scal_adv_nudge(what,0 &
     &                          ,tempv,idtemp &
     &                          ,thpar,wsink &
     &                          ,difhv,difv,difmol &
     &				,tobsv,robs,ttauv)
	  end if

!$OMP END TASK

!	  call openmp_get_thread_num(tid)
!	  !write(6,*) 'number of thread of salt: ',tid

!$OMP TASK PRIVATE(what,dtime) FIRSTPRIVATE(shpar,wsink,robs,isalt) &
!$OMP&     SHARED(idsalt,saltv,difhv,difv,difmol,sobsv,stauv) &
!$OMP&     DEFAULT(NONE) &
!$OMP&     IF(isalt > 0)

          if( isalt .gt. 0 ) then
		what = 'salt'
                call scal_adv_nudge(what,0 &
     &                          ,saltv,idsalt &
     &                          ,shpar,wsink &
     &                          ,difhv,difv,difmol &
     &				,sobsv,robs,stauv)
          end if

!$OMP END TASK
!$OMP TASKWAIT

	  !call ts_dia('after T/D')

	end if

	if( binfo ) then
	  format = '(a,e18.9,3f8.3)'
          if( itemp .gt. 0 ) then
  	    !call tsmass(tempv,+1,nlvdi,ttot) 
	    call scalar_mass(+1,tempv,ttot,vtot)
      	    call conmima(nlvdi,tempv,tmin,tmax)
	    ttot = shympi_sum(ttot)
	    vtot = shympi_sum(vtot)
	    tmin = shympi_min(tmin)
	    tmax = shympi_max(tmax)
	    taver = ttot / vtot
	    array = (/ttot,tmin,tmax,taver/)
!$OMP CRITICAL
	    call info_output(' temp','none',4,array,.true.,format)
	    !if( iuinfo > 0 ) then
  	    !  write(iuinfo,*) 'temp: ',aline,ttot,tmin,tmax
	    !end if
!$OMP END CRITICAL
	  end if
          if( isalt .gt. 0 ) then
  	    !call tsmass(saltv,+1,nlvdi,stot) 
	    call scalar_mass(+1,saltv,stot,vtot)
       	    call conmima(nlvdi,saltv,smin,smax)
	    stot = shympi_sum(stot)
	    vtot = shympi_sum(vtot)
	    smin = shympi_min(smin)
	    smax = shympi_max(smax)
	    saver = stot / vtot
	    array = (/stot,smin,smax,saver/)
!$OMP CRITICAL
	    call info_output(' salt','none',4,array,.true.,format)
	    !if( iuinfo > 0 ) then
  	    !  write(iuinfo,*) 'salt: ',aline,stot,smin,smax
	    !end if
!$OMP END CRITICAL
	  end if
	end if

!----------------------------------------------------------
! heat flux through surface	!ccf --> moved to meteo_force
!----------------------------------------------------------

	call compute_heat_flux

!----------------------------------------------------------
! compute rhov
!----------------------------------------------------------

	call ts_dia('normal before rhoset_shell')
	allocate(rho_aux1(nlvdi,nkn),rho_aux2(nlvdi,nkn))
	rho_aux1 = rhov
	call rhoset1
	rho_aux2 = rhov
	rhov = rho_aux1
	call rhoset_shell(1)
	!call ts_dia('normal after rhoset_shell')

!----------------------------------------------------------
! compute min/max
!----------------------------------------------------------

	call stmima(saltv,nkn,nlvdi,ilhkv,smin,smax)
	call stmima(tempv,nkn,nlvdi,ilhkv,tmin,tmax)
	call stmima(rhov,nkn,nlvdi,ilhkv,rmin,rmax)

!----------------------------------------------------------
! write results to file
!----------------------------------------------------------

	call bcl_write_output(dtime,da_out,itemp,isalt,irho)

!----------------------------------------------------------
! debug
!----------------------------------------------------------

	ks = 30280
	ks = 29988
	ks = 0

	if( ks > 0 ) then
	  call check_node(ks)
	  call check_elems_around_node(ks)
	  call check_nodes_around_node(ks)
	end if

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	end

!********************************************************
!********************************************************
!********************************************************

	subroutine rhoset_shell(iterwant)

! sets rho iterating to real solution

	use mod_ts
	use shympi

	implicit none

	integer iterwant	!want these iterations, if 0 do as needed

	logical biter
	integer itermax,iter,iu
	real eps,resid,resid_old

	itermax = 10
	eps = 1.e-7

	biter = .true.
	iter = 0
	resid = 0.
	resid_old = 0.

	do while( biter )
	  resid_old = resid
          call rhoset(resid)
	  iter = iter + 1
	  if( iterwant > 0 ) then 		!impose number of iterations
	    if( iter .lt. iterwant ) cycle	!must do more
	    if( iter .ge. iterwant ) exit	!finished
	  end if
	  if( resid .lt. eps ) biter = .false.
	  if( abs(resid-resid_old) .lt. eps ) biter = .false.
	  if( iter .gt. itermax ) biter = .false.
	end do

	if( iter .gt. itermax ) then
	  write(6,*) '*** warning: max iterations in rhoset_shell '
	  write(6,*) '    strange values of T/S have been found'
	  write(6,*) '    resid,iter: ',resid,iter
	  call tsrho_check
	end if

	!iu = 666+my_id
	!write(iu,*) my_id,iter
	!flush(iu)

	call shympi_exchange_3d_node(rhov)

	end

!********************************************************

	subroutine rhoset(resid)

! computes rhov
!
! 1 bar = 100 kPascal ==> factor 1.e-5
! pres = rho0*g*(zeta-z) + presbc
! with presbc = int_{z}^{zeta}(g*rho_prime)dz
! and rho_prime = rho - rho_0 = sigma - sigma_0
!
! in rhov()   is rho_prime (=sigma_prime)
!
! presbc and rhov() are given at node and layer center

	use mod_layer_thickness
	use mod_ts
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use pkonst

	implicit none

	real resid
! common

! local
	logical bdebug,debug,bsigma
	integer k,l,lmax,lmin
	integer nresid,nsigma
	real sigma0,rho0,pres,hsigma
	real depth,hlayer,hh
	real rhop,presbt,presbc,dpresc
	real salt
	double precision dresid
! functions
	real sigma_stp

	rho0 = rowass
	sigma0 = rho0 - 1000.

	debug=.false.
	bdebug=.false.

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	if(debug) write(6,*) sigma0,rowass,rho0

	nresid = 0
	dresid = 0.

	do k=1,nkn
	  depth = 0.
	  presbc = 0.
	  lmin = jlhkv(k)
	  lmax = ilhkv(k)

	  do l=lmin,lmax
	    bsigma = l .le. nsigma

	    hlayer = hdkov(l,k)
	    if( .not. bsigma ) hlayer = hldv(l)

	    hh = 0.5 * hlayer
	    depth = depth + hh
	    rhop = rhov(l,k)			!rho^prime

	    dpresc = rhop * grav * hh		!differential bc. pres.
	    presbc = presbc + dpresc            !baroclinic pres. (mid-layer)
	    presbt = rho0 * grav * depth	!barotropic pressure

	    pres = 1.e-5 * ( presbt + presbc )	!pressure in bars (BUG)
	
	    salt = max(0.,saltv(l,k))
	    rhop = sigma_stp(salt,tempv(l,k),pres) - sigma0
	    call set_rhomud(k,l,rhop)

	    nresid = nresid + 1
	    dresid = dresid + (rhov(l,k)-rhop)**2

	    rhov(l,k) = rhop

	    depth = depth + hh
	    presbc = presbc + dpresc		!baroclinic pres. (bottom-lay.)
	  end do
	end do

	resid = dresid/nresid

	return
	end

!********************************************************

	subroutine rhoset1

! computes rhov
!
! 1 bar = 100 kPascal ==> factor 1.e-5
! pres = rho0*g*(zeta-z) + presbc
! with presbc = int_{z}^{zeta}(g*rho_prime)dz
! and rho_prime = rho - rho_0 = sigma - sigma_0
!
! in rhov()   is rho_prime (=sigma_prime)
!
! presbc and rhov() are given at node and layer center

	use mod_layer_thickness
	use mod_ts
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use pkonst

	implicit none

! common

! local
	logical bdebug,debug,bsigma
	integer k,l,lmax
	integer nresid,nsigma
	integer iter,iter_max
	real sigma0,rho0,pres,hsigma
	real depth,hlayer,hh
	real rhop,presbt,presbc,dpresc
	real salt
	real resid
	real eps
	double precision dresid
! functions
	real sigma_stp

	iter_max = 10
	eps = 1.e-7

	rho0 = rowass
	sigma0 = rho0 - 1000.

	debug=.false.
	bdebug=.false.

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	if(debug) write(6,*) sigma0,rowass,rho0

	do k=1,nkn
	 iter = 0
	 lmax = ilhkv(k)
	 do
	  iter = iter + 1
	  depth = 0.
	  presbc = 0.
	  dresid = 0.
	  do l=1,lmax
	    bsigma = l .le. nsigma

	    hlayer = hdkov(l,k)
	    if( .not. bsigma ) hlayer = hldv(l)

	    hh = 0.5 * hlayer
	    depth = depth + hh
	    rhop = rhov(l,k)			!rho^prime

	    dpresc = rhop * grav * hh		!differential bc. pres.
	    presbc = presbc + dpresc            !baroclinic pres. (mid-layer)
	    presbt = rho0 * grav * depth	!barotropic pressure

	    pres = 1.e-5 * ( presbt + presbc )	!pressure in bars (BUG)
	
	    salt = max(0.,saltv(l,k))
	    rhop = sigma_stp(salt,tempv(l,k),pres) - sigma0
	    call set_rhomud(k,l,rhop)

	    dresid = dresid + (rhov(l,k)-rhop)**2

	    rhov(l,k) = rhop

	    depth = depth + hh
	    presbc = presbc + dpresc		!baroclinic pres. (bottom-lay.)
	  end do
	  resid = dresid/lmax
	  if( resid < eps ) exit
	  if( iter > iter_max ) goto 99
	 end do
	end do

	return
   99	continue
	write(6,*) 'error iterating rho'
	write(6,*) 'k,lmax: ',k,lmax
	write(6,*) 'resid,eps: ',resid,eps
	write(6,*) '#  layer     salt    temp    rho'
	do l=1,lmax
	  write(6,*) l,saltv(l,k),tempv(l,k),rhov(l,k)
	end do
	stop 'error stop rhoset1: too many iterations'
	end

!*******************************************************************	

	subroutine tsrho_check

! checks values of t/s/rho

	use mod_ts
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real smin,smax,tmin,tmax,rmin,rmax
	character*30 text

	text = '*** tsrho_check'

	call stmima(saltv,nkn,nlvdi,ilhkv,smin,smax)
	call stmima(tempv,nkn,nlvdi,ilhkv,tmin,tmax)
	call stmima(rhov,nkn,nlvdi,ilhkv,rmin,rmax)

	write(6,*) 'S   min/max: ',smin,smax
	write(6,*) 'T   min/max: ',tmin,tmax
	write(6,*) 'Rho min/max: ',rmin,rmax

	write(6,*) 'checking for Nans...'
        call check2Dr(nlvdi,nlv,nkn,tempv,-30.,+70.,text,'tempv')
        call check2Dr(nlvdi,nlv,nkn,saltv, -1.,+70.,text,'saltv')
        call check2Dr(nlvdi,nlv,nkn,rhov,-2000.,+2000.,text,'rhov')

	end

!*******************************************************************	
!*******************************************************************	
!*******************************************************************	

	subroutine ts_diag(dtime,nlvddi,nlv,nkn,tempv,saltv)

	implicit none

	double precision dtime
	integer nlvddi
	integer nlv
	integer nkn
	real tempv(nlvddi,nkn)
	real saltv(nlvddi,nkn)

	character*80 tempf,saltf
	character*80 string
	integer, save :: idtemp,idsalt
	integer, save :: icall = 0

	call getfnm('tempobs',tempf)
	call getfnm('saltobs',saltf)

	if( icall .eq. 0 ) then
	  string = 'temp diag'
	  call ts_open(string,tempf,dtime,nkn,nlv,idtemp)
	  string = 'salt diag'
	  call ts_open(string,saltf,dtime,nkn,nlv,idsalt)
	  icall = 1
	end if

        call ts_next_record(dtime,idtemp,nlvddi,nkn,nlv,tempv)
        call ts_next_record(dtime,idsalt,nlvddi,nkn,nlv,saltv)

	end

!*******************************************************************	

	subroutine ts_nudge(dtime,nlvddi,nlv,nkn,tobsv,sobsv &
     &					,ttauv,stauv)

	implicit none

	double precision dtime
	integer nlvddi
	integer nlv
	integer nkn
	real tobsv(nlvddi,nkn)
	real sobsv(nlvddi,nkn)
	real ttauv(nlvddi,nkn)
	real stauv(nlvddi,nkn)

	character*80 tempf,saltf,ttauf,stauf
	character*80 string
	real ttaup,staup
	logical, save :: btnudge,bsnudge
	integer, save :: idtemp,idsalt
	integer, save :: idttau,idstau
	integer, save :: icall = 0

	real getpar

	call getfnm('tempobs',tempf)
	call getfnm('saltobs',saltf)
	call getfnm('temptau',ttauf)
	call getfnm('salttau',stauf)
	ttaup = getpar('temptaup')
	staup = getpar('salttaup')

	if( icall .eq. 0 ) then
	  string = 'temp tau'
	  write(6,'(a)') 'ts_nudge: opening file for '//trim(string)
	  call ts_nudge_get_tau(string,ttauf,ttaup,dtime,nkn,nlv,idttau)
	  btnudge = idttau > 0 .or. ttaup > 0.
	  if( btnudge ) then
	    if( ttaup > 0. ) ttauv = 1./ttaup
	    string = 'temp nudge'
	    write(6,'(a)') 'ts_nudge: opening file for '//trim(string)
	    call ts_open(string,tempf,dtime,nkn,nlv,idtemp)
	  end if

	  string = 'salt tau'
	  write(6,'(a)') 'ts_nudge: opening file for '//trim(string)
	  call ts_nudge_get_tau(string,stauf,staup,dtime,nkn,nlv,idstau)
	  bsnudge = idstau > 0 .or. staup > 0.
	  if( bsnudge ) then
	    if( staup > 0. ) stauv = 1./staup
	    string = 'salt nudge'
	    write(6,'(a)') 'ts_nudge: opening file for '//trim(string)
	    call ts_open(string,saltf,dtime,nkn,nlv,idsalt)
	  end if

	  if( btnudge .or. bsnudge ) then
	    write(6,*) 'nudging has been initialized'
	    if( btnudge ) write(6,*) '  nudging for temperature is active'
	    if( bsnudge ) write(6,*) '  nudging for salinity is active'
	  else
	    write(6,*) 'nudging requested but no files found'
	    write(6,*) 'files to be set:'
	    write(6,*) 'tempobs and saltobs for observations'
	    write(6,*) 'temptau and salttau for relaxation time scale'
	    write(6,*) 'in alternative parameters to be set:'
	    write(6,*) 'temptaup and salttaup for relaxation time scale'
	    stop 'error stop ts_nudge: no files found for nudging'
	  end if

	  icall = 1
	end if

	if( btnudge ) then
	  if( idttau > 0 ) then
            call ts_next_record(dtime,idttau,nlvddi,nkn,nlv,ttauv)
	    where( ttauv > 0. ) ttauv = 1./ttauv
	  end if
          call ts_next_record(dtime,idtemp,nlvddi,nkn,nlv,tobsv)
	end if

	if( bsnudge ) then
	  if( idstau > 0 ) then
            call ts_next_record(dtime,idstau,nlvddi,nkn,nlv,stauv)
	    where( stauv > 0. ) stauv = 1./stauv
	  end if
          call ts_next_record(dtime,idsalt,nlvddi,nkn,nlv,sobsv)
	end if

	end

!*******************************************************************	

	subroutine ts_nudge_get_tau(string,file,tau,dtime,nkn,nlv,id)

! returns id - if >0 read from file, if =0 tau has time scale

	implicit none

	character*(*) string,file
	real tau
	double precision dtime
	integer nkn,nlv
	integer id

	if( file /= ' ' ) then
	  !write(6,'(a)') 'ts_nudge: opening file for '//trim(string)
	  call ts_open(string,file,dtime,nkn,nlv,id)
	else if( tau < 0. ) then	!no tau given
	  id = -1
	  !goto 99
	else				!tau specified - do not read file
	  id = 0
	end if

	return
   99	continue
	write(6,*) '*** error preparing for nudging: '//trim(string)
	write(6,*) 'no nudging time scale given...'
	write(6,*) 'please provide filename in temptau and salttau'
	write(6,*) 'or set nudging time scale using parameters'
	write(6,*) 'temptaup and salttaup'
	stop 'error stop ts_nudge_get_tau: error getting tau'
	end

!*******************************************************************	

	subroutine ts_open(string,file,dtime,nkn,nlv,id)

	implicit none

	character*(*) string,file
	double precision dtime
	integer nkn,nlv
	integer id

	logical bexist

	call ts_file_exists(file,bexist)
	if( .not. bexist ) goto 99
	call ts_file_open(file,dtime,nkn,nlv,id)
	if( id <= 0 ) goto 99
	call ts_file_descrp(id,string)

	return
   99	continue
	write(6,*) '*** error opening file for '//trim(string)
	if( file == ' ' ) then
	  write(6,*) 'no file given...'
	else if( .not. bexist ) then
	  write(6,*) 'file does not exisit: ',trim(file)
	else
	  write(6,*) 'error opening file: ',trim(file)
	end if
	stop 'error stop ts_open: error opening file'
	end

!*******************************************************************	

	subroutine ts_init(dtime,nlvddi,nlv,nkn,tempv,saltv)

! initialization of T/S from file

	use shympi

	implicit none

	double precision dtime
        integer nlvddi
        integer nlv
        integer nkn
        real tempv(nlvddi,nkn)
        real saltv(nlvddi,nkn)

        character*80 tempf,saltf
        character*80 string
        integer id

	call getfnm('tempin',tempf)
	call getfnm('saltin',saltf)

	if( tempf .ne. ' ' ) then
          write(6,'(a)') 'initializing temperature from file ' &
     &                        //trim(tempf)
	  string = 'temp init'
	  call ts_open(string,tempf,dtime,nkn,nlv,id)
          call ts_next_record(dtime,id,nlvddi,nkn,nlv,tempv)
	  call ts_file_close(id)
	  call shympi_exchange_3d_node(tempv)
          write(6,*) 'temperature initialized from file ',trim(tempf)
	end if

	if( saltf .ne. ' ' ) then
          write(6,'(a)') 'initializing salinity from file ' &
     &                        //trim(saltf)
	  string = 'salt init'
	  call ts_open(string,saltf,dtime,nkn,nlv,id)
          call ts_next_record(dtime,id,nlvddi,nkn,nlv,saltv)
	  call ts_file_close(id)
	  call shympi_exchange_3d_node(saltv)
          write(6,*) 'salinity initialized from file ',trim(saltf)
	end if

	end

!*******************************************************************	
!*******************************************************************	
!*******************************************************************	

	subroutine bcl_open_output(da_out,itemp,isalt,irho)

! opens output of T/S

	use levels
	use mod_trace_point

	implicit none

	double precision da_out(4)
	integer itemp,isalt,irho

	logical b2d
	integer nvar,id
	logical has_output_d

	b2d = .false.
	nvar = 0
	if( itemp .gt. 0 ) nvar = nvar + 1
	if( isalt .gt. 0 ) nvar = nvar + 1
	if( irho  .gt. 0 ) nvar = nvar + 1

	call init_output_d('itmcon','idtcon',da_out)

	if( has_output_d(da_out) ) then
	  call trace_point('initalizing scalar file of T/S')
	  call shyfem_init_scalar_file('ts',nvar,b2d,id)
	  call trace_point('finished initalizing scalar file of T/S')
	  da_out(4) = id
	end if

	end

!*******************************************************************	

	subroutine bcl_write_output(dtime,da_out,itemp,isalt,irho)

! writes output of T/S

	use levels
	use mod_ts

	implicit none

	double precision dtime
	double precision da_out(4)
	integer itemp,isalt,irho

	integer id
	logical next_output_d
	
	if( next_output_d(da_out) ) then
	  id = nint(da_out(4))
	  if( isalt .gt. 0 ) then
	    call shy_write_scalar_record(id,dtime,11,nlvdi,saltv)
	  end if
	  if( itemp .gt. 0 ) then
	    call shy_write_scalar_record(id,dtime,12,nlvdi,tempv)
	  end if
	  if( irho  .gt. 0 ) then
	    call shy_write_scalar_record(id,dtime,13,nlvdi,rhov)
	  end if
	  call shy_sync(id)
	end if

	end

!*******************************************************************	
!*******************************************************************	
!*******************************************************************	

	subroutine getts(l,k,t,s)

! accessor routine to get T/S

	use mod_ts

        implicit none

        integer k,l
        real t,s

        t = tempv(l,k)
        s = saltv(l,k)

        end

!******************************************************************
!*******************************************************************	
!*******************************************************************	

	subroutine ts_dia(string)

	use mod_ts
	use levels
	use basin

	implicit none

	character*(*) string

	logical berr
	integer ke,k,l,lmax
	real tmin,tmax,smin,smax
	real t,s
	character*20 aline
	logical, save :: bwrite = .false.
	logical, save :: bstop = .false.
	real, parameter :: tstop = 100.
	real, parameter :: twrite = 50.

	integer ipext

	tmin = minval(tempv)
	tmax = maxval(tempv)
	smin = minval(saltv)
	smax = maxval(saltv)

	if( tmin < -tstop ) bstop = .true.
	if( tmax >  tstop ) bstop = .true.
	if( smin < -tstop ) bstop = .true.
	if( smax >  tstop ) bstop = .true.

	call get_act_timeline(aline)

	if( bstop ) then
	  write(6,*) 'ts_dia: ',trim(string),'  ',aline
	  write(6,*) 'saltv (min/max): ',smin,smax
	  write(6,*) 'tempv (min/max): ',tmin,tmax
	  write(6,*) 'for list of nodes see file fort.166'
	  bwrite = .true.
	end if

	if( bwrite ) then
	  write(166,*) 'ts_dia: ',trim(string),'  ',aline
	  write(166,*) 'saltv (min/max): ',smin,smax
	  write(166,*) 'tempv (min/max): ',tmin,tmax
	  write(166,*) 'list of nodes with strange values:'
	  write(166,*) '      layer     kintern     kextern' // &
     &				'            t                s'
	  do k=1,nkn
	    lmax = ilhkv(k)
	    ke = ipext(k)
	    do l=1,lmax
	      s = saltv(l,k)
	      t = tempv(l,k)
	      berr = .false.
	      if( t < -twrite .or. t > twrite ) berr = .true.
	      if( s < -twrite .or. s > twrite ) berr = .true.
	      if( berr ) then
		write(166,*) l,k,ke,t,s
	      end if
	    end do
	  end do
	end if

	if( bstop ) stop 'error stop ts_dia'

	end

!*******************************************************************	

