
!--------------------------------------------------------------------------
!
!    Copyright (C) 2007-2008,2010-2011,2013-2015,2013-2015  Georg Umgiesser
!    Copyright (C) 2018-2019  Georg Umgiesser
!    Copyright (C) 2008,2011,2013,2016  Debora Bellafiore
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

! explicit term routines
!
! contents :
!
! subroutine set_explicit
! subroutine viscous_stability(ahpar,ahstab)	computes stability for viscosity
! subroutine set_diff_horizontal
! subroutine set_advective
! subroutine set_semi_lagrange
! subroutine set_barocl
! subroutine set_barocl_new
! subroutine set_barocl_new1
!
! revision log :
!
! 01.05.2007	ggu	new file -> all explicit terms here
! 28.09.2007	ggu	semi-lagrangian part introduced
! 16.04.2008	ggu	bugfix in set_barocl (real do indices!!)
! 14.07.2008	ggu&ccf	ahpar is real in set_diff_horizontal_new
! 03.11.2008	ggu&dbf	nudging implemented (call to bclnudge)
! 09.11.2008	ggu	set_barocl_new (cleaned version of set_barocl)
! 19.11.2008	ggu	new set_diff_horizontal_new1(), viscous_stability()
! 19.02.2010	ggu	in viscous_stability() for dt=1
! 26.02.2010	ggu	new call to momentum_viscous_stability()
! 26.02.2010	ggu	set_advective() cleaned up
! 26.02.2010	ggu	new momentum_advective_stability()
! 08.03.2010	ggu	run only down to avail layers (bug fix)
! 23.03.2010	ggu	changed v6.1.1
! 14.04.2010	ggu	changed v6.1.4
! 16.12.2010	ggu	barocl preconditioned for sigma layers, but not finshed
! 27.01.2011	ggu	changed VERS_6_1_17
! 20.05.2011	ggu	compute statistics of stability, no stab in dry elemes
! 31.05.2011	ggu	changed VERS_6_1_23
! 14.07.2011	ggu	changed VERS_6_1_27
! 25.08.2011	dbf&ggu	baroclinic gradient for sigma level integrated
! 01.09.2011	ggu	changed VERS_6_1_32
! 18.10.2011	ggu	changed VERS_6_1_33
! 25.10.2011	dbf&ggu	bug fix in set_barocl_new_interface (psigma)
! 04.11.2011	ggu	adapted for hybrid coordinates
! 10.05.2013	dbf&ggu	new routines for vertical advection (bvertadv)
! 10.05.2013	dbf&ggu	new routines for non-hydro
! 25.05.2013	ggu	new version for vertical advection (bvertadv)
! 13.06.2013	ggu	changed VERS_6_1_65
! 12.09.2013	ggu	changed VERS_6_1_67
! 13.09.2013	dbf&ggu	new sigma layer adjustment integrated
! 25.10.2013	ggu	changed VERS_6_1_68
! 05.12.2013	ggu	changed VERS_6_1_70
! 28.01.2014	ggu	changed VERS_6_1_71
! 10.04.2014	ggu	use rlin and rdistv to determin advective contribution
! 05.05.2014	ggu	changed VERS_6_1_74
! 18.06.2014	ggu	changed VERS_6_1_77
! 07.07.2014	ggu	changed VERS_6_1_79
! 05.11.2014	ggu	changed VERS_7_0_5
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 17.04.2015	ggu	only one routine set_diff_horizontal()
! 10.07.2015	ggu	changed VERS_7_1_50
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 18.09.2015	ggu	use momentx/yv to store advective terms, not aux arrays
! 23.09.2015	ggu	changed VERS_7_2_4
! 25.09.2015	ggu	new call to set_nudging()
! 12.10.2015	ggu	changed VERS_7_3_3
! 14.06.2016	dbf	diff. vertical momentum advection schemes implemented
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 26.05.2020	ggu	use rdistv now on elements and ruseterm
! 30.03.2022	ggu	compiler bug with PGI (PGI_ggguuu) - no solution
! 04.04.2022	ggu	exchange momentx/yv arrays
! 08.04.2022	ggu	ie_mpi introduced computing advective terms
! 09.04.2022	ggu	ie_mpi also in baroclinic section, some debug code
! 15.10.2022	ggu	some new debug code, bpresxv,bpresyv local
! 21.10.2022	ggu	allocate big array saux that was on stack
! 18.03.2023	ggu	adjusted horizontal diffusion for mpi
! 27.03.2023	ggu	tripple point routines in new file submpi_tripple.f
! 29.03.2023	ggu	handle tripple points, new horizontal diffusion routine
! 29.03.2023	ggu	exchange rindex between domains
! 09.05.2023    lrp     introduce top layer index variable
! 24.05.2023    ggu     momentum_viscous_stability(): must run over nel_unique
! 05.06.2023    lrp     introduce z-star
! 13.04.2024    ggu     use tripple_point_set_values, also in hydro_stability
! 13.10.2024    ggu     bug fix for INTEL_BUG
! 01.12.2024    ggu     in momentum_advective_stability() do not reduce
! 24.04.2025    ggu     handle new value for ibarcl == 5
!
! notes :
!
! sign of explicit term is computed for left hand side
!
!******************************************************************

	subroutine set_explicit

	use mod_internal
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none
        
	integer ie,l
        
        logical bbarcl,blin
        integer ilin,itlin,ibarcl
	double precision dtime
	real rlin
        real getpar
	logical bnohyd

!-------------------------------------------
! parameters
!-------------------------------------------

        ilin = nint(getpar('ilin'))
        rlin = getpar('rlin')
        itlin = nint(getpar('itlin'))
        ibarcl = nint(getpar('ibarcl'))
        bbarcl = ibarcl .gt. 0 .and. ibarcl .ne. 3 .and. ibarcl .ne. 5
	blin = ilin /= 0
	call nonhydro_get_flag(bnohyd)
	call get_act_dtime(dtime)
	
	!write(6,*) 'in explicit: ',blin,bbarcl,bnohyd

!-------------------------------------------
! initialization
!-------------------------------------------

	fxv = 0.
	fyv = 0.

!-------------------------------------------
! fix or nudge boundary velocities
!-------------------------------------------

	call bclfix

!-------------------------------------------
! horizontal diffusion
!-------------------------------------------

	call set_diff_horizontal

!-------------------------------------------
! advective (non-linear) terms
!-------------------------------------------

        if( .not. blin ) then
          if( itlin .eq. 0 ) then
	    call set_advective(rlin)
	  else if( itlin .eq. 1 ) then
	    call set_semi_lagrange
	  else
	    write(6,*) 'itlin = ',itlin
	    stop 'error stop set_explicit: no such option'
	  end if
	end if

!-------------------------------------------
! baroclinic contribution
!-------------------------------------------

	if( bbarcl ) call set_barocl_new_interface

!-------------------------------------------
! non-hydrostatic contribution (experimental)
!-------------------------------------------

	if( bnohyd ) call nonhydro_set_explicit

!-------------------------------------------
! nudging of water levels and velocities
!-------------------------------------------

	call set_nudging

!-------------------------------------------
! end of routine
!-------------------------------------------

	end

!******************************************************************

        subroutine momentum_viscous_stability_old(ahpar,rindex,dstab)

! computes stability for viscosity
!
! stability is computed for dt == 1

        use mod_geom
        use mod_internal
        use mod_diff_visc_fric
        use evgeom
        use levels
        use basin
	use shympi

        implicit none

        real ahpar
        real rindex
        real dstab(nlvdi,nel)

        logical bdebug
        integer ie,ii,iei,l,lmax,lmin,lmaxi
        real u,v,ui,vi
        real anu,ax,ay
        real area,areai
        real dt
        real a,ai,amax,afact,r

        rindex = 0.
        if( ahpar .le. 0 ) return

        amax = 0.

        do ie=1,nel

          lmax = ilhv(ie)
          lmin = jlhv(ie)
          area = 12. * ev(10,ie)
          r = rdistv(ie)

          do l=lmin,lmax

            a = 0.
            do ii=1,3
              iei = ieltv(ii,ie)
              if( iei .le. 0 ) iei = ie

              lmaxi = ilhv(iei)
              !if( l > lmaxi ) cycle
              areai = 12. * ev(10,iei)

              anu = ahpar * difhv(l,ie)
              ai = 2. * anu / ( area + areai )
              a = a + ai
            end do

            a = a * r
            amax = max(amax,a)
            dstab(l,ie) = a

          end do

        end do

        !rindex = shympi_max(amax)
	rindex = amax

        end

!******************************************************************

	subroutine momentum_viscous_stability(ahpar,rindex,dstab)

! computes stability for viscosity
!
! stability is computed for dt == 1

	use mod_geom
	use mod_internal
	use mod_diff_visc_fric
	use evgeom
	use levels
	use basin
	use shympi
	use shympi_tripple

	implicit none

	real ahpar
	real rindex
	real dstab(nlvdi,nel)

	integer ie,ii,iei,l,lmax,lmin,lmaxi
	integer itr
	integer, save :: icall = 0
	real u,v,ui,vi
	real anu,ax,ay
	real area,areai
	real dt
	real a,ai,amax,afact,r
	real al(nlvdi)

	!real, allocatable :: dstab2d(:),dstab2dg(:),rdistvg(:)
	!real, allocatable :: raux(:),rauxg(:)
	!integer, allocatable :: iext(:),iextg(:)

	rindex = 0.
	if( ahpar .le. 0 ) return

	amax = 0.
	dstab = 0.

	!allocate(dstab2d(nel),dstab2dg(nel_global))
	!allocate(iext(nel),iextg(nel_global))
	!allocate(raux(nel),rauxg(nel_global))
	!allocate(rdistvg(nel_global))

! in the next call we initialize lmaxi and areai, dstab is a dummy

	if( icall == 0 ) then		!just do it the first time
	  call tripple_points_init	!just in case
	  call tripple_point_set_values(nlvdi,evdim,nel,ev,dstab,dstab)
	  call tripple_points_exchange
	end if

	do ie=1,nel_unique

	  lmax = ilhv(ie)
          lmin = jlhv(ie)
	  area = 12. * ev(10,ie)
	  r = rdistv(ie)
	  al(1:lmax) = 0.

	  do ii=1,3

            iei = ieltv(ii,ie)
            if( iei .le. 0 ) iei = ie

	    lmaxi = ilhv(iei)
            areai = 12. * ev(10,iei)
	    itr = ietrp(ii,ie)
	    if( itr > 0 ) then
	      call tripple_point_get_la(itr,lmaxi,areai)
	    end if

	    do l=lmin,lmax
	      !if( l > lmaxi ) exit
	      anu = ahpar * difhv(l,ie)
              ai = 2. * anu / ( area + areai )
              al(l) = al(l) + ai
	    end do

          end do

	  al(1:lmax) = r * al(1:lmax)
	  dstab(1:lmax,ie) = al(1:lmax)
	  a = maxval(al(1:lmax))
	  !dstab2d(ie) = a
	  !raux(ie) = maxval(difhv(:,ie))
	  amax = max(amax,a)

	end do

	!rindex = shympi_max(amax)
	rindex = amax

! debug section

	icall = icall + 1

	end

!******************************************************************

	subroutine set_diff_horizontal

	use mod_geom
	use mod_internal
	use mod_diff_visc_fric
	use mod_hydro
	use evgeom
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi
	use shympi_tripple

	implicit none

	integer ie,ii,iei,l,lmax,lmaxi
	integer iext,ies,iu,ineib,nl,itr,ieneib
	integer noslip
	real u,v,ui,vi
	real anu,ahpar,ax,ay
	real area,areai
	real dt
	real a,ai,afact
	real rdist!,rcomp,ruseterm
	double precision dtime
	logical bnoslip,bdebug
	real, parameter :: axymax = 1.e+2

	real uiv(nlv_global)
	real viv(nlv_global)

	integer ieext,ieint
	real getpar

	call get_timestep(dt)
	call get_act_dtime(dtime)

        ahpar = getpar('ahpar')
	if( ahpar .le. 0 ) return

	nl = nlv_global

	call tripple_points_init	!just in case
	call tripple_point_set_values(nlvdi,evdim,nel,ev,utlov,vtlov)
	call tripple_points_exchange

        noslip = nint(getpar('noslip'))
	bnoslip = noslip .ne. 0

	do ie=1,nel_unique

          !rcomp = rcomputev(ie)           !use terms (custom elements)
          !ruseterm = min(rcomp,rdist)     !use terms (both)

	  lmax = ilhv(ie)
	  area = 12. * ev(10,ie)
          rdist = rdistv(ie)              !use terms (distance from OB)

	  do ii=1,3
            afact = 1.
            iei = ieltv(ii,ie)
            if( bnoslip .and. iei .eq. 0 ) afact = -1.
            if( iei .le. 0 ) iei = ie
	    lmaxi = ilhv(iei)
            areai = 12. * ev(10,iei)

	    itr = ietrp(ii,ie)
	    if( itr > 0 ) then
	      !write(6,*) 'handling tripple point ',itr
	      call tripple_point_get_values(itr,nl,lmaxi,areai,uiv,viv)
	    end if

	    do l=1,lmax
	      if( l > lmaxi ) exit		!INTEL_BUG

	      u  = utlov(l,ie)
	      v  = vtlov(l,ie)

	      ui = afact * utlov(l,iei)
	      vi = afact * vtlov(l,iei)
	      if( itr > 0 ) then
		ui = afact * uiv(l)
		vi = afact * viv(l)
	      end if

	      anu = ahpar * difhv(l,ie)
              ai = 2. * anu / ( area + areai )

	      ax = rdist * ai * ( ui - u )
	      ay = rdist * ai * ( vi - v )

	      ! we insert with minus sign because f is on left side

	      fxv(l,ie) = fxv(l,ie) - ax
	      fyv(l,ie) = fyv(l,ie) - ay
	    end do
          end do

	end do

	call shympi_exchange_3d_elem(fxv)		!ggu_diff
	call shympi_exchange_3d_elem(fyv)

	end

!******************************************************************

	subroutine set_diff_horizontal_old

! old routine - cannot handle tripple points

	use mod_geom
	use mod_internal
	use mod_diff_visc_fric
	use mod_hydro
	use evgeom
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

	implicit none

	integer ie,ii,iei,l,lmax,lmaxi
	integer iext,ies,iu
	integer noslip
	real u,v,ui,vi
	real anu,ahpar,ax,ay
	real area,areai
	real dt
	real a,ai,amax,afact
	real rdist!,rcomp,ruseterm
	logical bnoslip,bdebug

	integer ieext,ieint
	real getpar

	call get_timestep(dt)

	call tripple_points_init	!just in case
	call tripple_point_set_values(nlvdi,evdim,nel,ev,utlov,vtlov)
	call tripple_points_exchange

        ahpar = getpar('ahpar')
	if( ahpar .le. 0 ) return

        noslip = nint(getpar('noslip'))
	bnoslip = noslip .ne. 0

	amax = 0.

	do ie=1,nel_unique
	!do ie=1,nel

          rdist = rdistv(ie)              !use terms (distance from OB)
          !rcomp = rcomputev(ie)           !use terms (custom elements)
          !ruseterm = min(rcomp,rdist)     !use terms (both)

	  lmax = ilhv(ie)
	  area = 12. * ev(10,ie)

	  do l=1,lmax
	    u  = utlov(l,ie)
	    v  = vtlov(l,ie)

	    a = 0.
	    do ii=1,3
              iei = ieltv(ii,ie)
              afact = 1.
              if( bnoslip .and. iei .eq. 0 ) afact = -1.
              if( iei .le. 0 ) iei = ie

	      lmaxi = ilhv(iei)
	      !if( l > lmaxi ) cycle
              areai = 12. * ev(10,iei)

	      anu = ahpar * difhv(l,ie)
              ai = 2. * anu / ( area + areai )
              a = a + ai

	      ui = afact * utlov(l,iei)
	      vi = afact * vtlov(l,iei)

	      ax = rdist * ai * ( ui - u )
	      ay = rdist * ai * ( vi - v )

	      fxv(l,ie) = fxv(l,ie) - ax	!minus because f is on left side
	      fyv(l,ie) = fyv(l,ie) - ay
	    end do
	    amax = max(amax,a)
          end do

	end do

	call shympi_exchange_3d_elem(fxv)		!ggu_diff
	call shympi_exchange_3d_elem(fyv)

	end

!******************************************************************

        subroutine set_momentum_flux

! sets arrays momentx/yv

	use mod_layer_thickness
	use mod_hydro_print
	use mod_hydro
	use mod_internal
	use evgeom
	use levels
	use basin
	use shympi

        implicit none

        integer ii,ie,k,l,lmax,lmin,ie_mpi
        real b,c
        real ut,vt
        real uc,vc
        real up,vp
        real um,vm
        real f,h
	real xadv,yadv
	real area,vol

	real, allocatable :: saux(:,:)

!---------------------------------------------------------------
! initialization
!---------------------------------------------------------------

	allocate(saux(nlvdi,nkn))

	saux = 0.
	momentxv = 0.
	momentyv = 0.

!---------------------------------------------------------------
! accumulate momentum that flows into nodes (weighted by flux)
!---------------------------------------------------------------

	do ie_mpi=1,nel
	  ie = ip_sort_elem(ie_mpi)
	  lmax = ilhv(ie)
	  lmin = jlhv(ie)
	  do l=lmin,lmax
            h = hdeov(l,ie)
	    ut = utlov(l,ie)
	    vt = vtlov(l,ie)
            do ii=1,3
                k = nen3v(ii,ie)
                b = ev(3+ii,ie)
                c = ev(6+ii,ie)
                f = ut * b + vt * c	! f>0 => flux into node
                if( f .gt. 0. ) then
		  saux(l,k) = saux(l,k) + f
		  momentxv(l,k) = momentxv(l,k) + f * ut
		  momentyv(l,k) = momentyv(l,k) + f * vt
                end if
	    end do
          end do
	end do

!---------------------------------------------------------------
! compute average momentum for every node
!---------------------------------------------------------------

	do k=1,nkn
	  lmax = ilhkv(k)
	  lmin = jlhkv(k)
	  do l=lmin,lmax
            h = hdkov(l,k)
	    if( saux(l,k) .gt. 0 ) then		!flux into node
	      momentxv(l,k) = momentxv(l,k) / saux(l,k)
	      momentyv(l,k) = momentyv(l,k) / saux(l,k)
	    else				!only flux out of node
	      momentxv(l,k) = uprv(l,k) * h
	      momentyv(l,k) = vprv(l,k) * h
	    end if
	  end do
	end do

!---------------------------------------------------------------
! exchange arrays
!---------------------------------------------------------------

	call shympi_exchange_3d_node(momentxv)
	call shympi_exchange_3d_node(momentyv)

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	end

!******************************************************************

        subroutine set_advective(rlin)

	use mod_internal
	use mod_layer_thickness
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin
	use shympi

        implicit none

	real rlin		!strength of advection terms - normally 1

	real zxadv,zyadv
	real wtop,wbot

        integer ihwadv  	! vertical advection for momentum
        integer ii,ie,k,l,lmax,lmin,ie_mpi
        real b,c
        real ut,vt
        real uc,vc
        real up,vp
        real um,vm
        real f,h
	real xadv,yadv
	real area,vol
        real getpar
        real wlay,dzbb,dz,dztt,ubot,utop,vbot,vtop
	real rdist,rcomp,ruseterm

	!write(6,*) 'set_advective called...'

!---------------------------------------------------------------
! initialization
!---------------------------------------------------------------

        ihwadv = nint(getpar('ihwadv'))

!---------------------------------------------------------------
! accumulate momentum that flows into nodes (weighted by flux)
!---------------------------------------------------------------

        call set_momentum_flux	!sets aux arrays momentx/yv

!---------------------------------------------------------------
! compute advective contribution
!---------------------------------------------------------------

	uc = 0.
	vc = 0.

	do ie_mpi=1,nel
	  ie = ip_sort_elem(ie_mpi)
          wtop = 0.0
	  lmax = ilhv(ie)
	  lmin = jlhv(ie)

          rdist = rdistv(ie)              !use terms (distance from OB)
          rcomp = rcomputev(ie)           !use terms (custom elements)
          ruseterm = min(rcomp,rdist)     !use terms (both)

	  do l=lmin,lmax

!	    ---------------------------------------------------------------
!	    horizontal advection
!	    ---------------------------------------------------------------

	    area = 12. * ev(10,ie)
            h = hdeov(l,ie)
	    vol = area * h
  	    ut = utlov(l,ie)
  	    vt = vtlov(l,ie)
	    !this throws a floating point exception with PGI (PGI_ggguuu)
	    !write(6,*) 'PGI_ggguuu adv a ',ie,l,lmax,h
	    !write(6,*) 'PGI_ggguuu adv h = ',h
	    !write(6,*) 'PGI_ggguuu adv b ',ut,vt,uc,vc
	    !write(6,*) 'PGI_ggguuu adv c ',ut*h,vt*h
	    !write(6,*) 'PGI_ggguuu adv d ',ut/h,vt/h
	    !flush(6)	!PGI_ggguuu
            uc = ut / h
            vc = vt / h
	    xadv = 0.
	    yadv = 0.
	    wbot = 0.
            do ii=1,3
                k = nen3v(ii,ie)
                b = ev(3+ii,ie)
                c = ev(6+ii,ie)
		wbot = wbot + wlov(l,k)
                up = momentxv(l,k) / h		!NEW
                vp = momentyv(l,k) / h
                f = ut * b + vt * c
                if( f .lt. 0. ) then	!flux out of node => into element
                  xadv = xadv + f * ( up - uc )
                  yadv = yadv + f * ( vp - vc )
                end if
            end do
	    
	    zxadv = 0.
	    zyadv = 0. 
	   
!	    ---------------------------------------------------------------
!	    vertical advection
!	    ---------------------------------------------------------------

	    if( ihwadv > 0 ) then	!compute vertical momentum advection
	      wbot = wbot / 3.
	      if( l .eq. lmax ) wbot = 0.

              if(ihwadv == 1) then	!use upwind scheme
  	        if(wtop.gt.0.) then
    	          zxadv = wtop * ulov(l,ie)
	          zyadv = wtop * vlov(l,ie)
                else
	          zxadv = wtop * ulov(l-1,ie)
	          zyadv = wtop * vlov(l-1,ie)
                end if

	        if(wbot.gt.0.) then
	          zxadv = zxadv - wbot * ulov(l+1,ie)
	          zyadv = zyadv - wbot * vlov(l+1,ie)
                else
	          zxadv = zxadv - wbot * ulov(l,ie)
                  zyadv = zyadv - wbot * vlov(l,ie)
                end if
              else if(ihwadv == 2) then	!use centered scheme
                dz = hdeov(l,ie)
                if (l .eq. 1) then
                  dzbb = hdeov(l+1,ie)
                  utop = 0.0
                  ubot = (ulov(l,ie)*dz+ulov(l+1,ie)*dzbb)/(dz+dzbb)
                  vtop = 0.0
                  vbot = (vlov(l,ie)*dz+vlov(l+1,ie)*dzbb)/(dz+dzbb)
                else if (l .eq. lmax) then
                  dztt = hdeov(l-1,ie)
                  utop = (ulov(l-1,ie)*dztt+ulov(l,ie)*dz)/(dztt+dz)
                  ubot = 0.0
                  vtop = (vlov(l-1,ie)*dztt+vlov(l,ie)*dz)/(dztt+dz)
                  vbot = 0.0
                else
                  dztt = hdeov(l-1,ie)
                  dzbb = hdeov(l+1,ie)
                  utop = (ulov(l-1,ie)*dztt+ulov(l,ie)*dz)/(dztt+dz)
                  ubot = (ulov(l,ie)*dz+ulov(l+1,ie)*dzbb)/(dz+dzbb)
                  vtop = (vlov(l-1,ie)*dztt+vlov(l,ie)*dz)/(dztt+dz)
                  vbot = (vlov(l,ie)*dz+vlov(l+1,ie)*dzbb)/(dz+dzbb)
                end if
                wlay = (wtop + wbot)/2.0
                zxadv = zxadv + wlay * (utop - ubot)
                zyadv = zyadv + wlay * (vtop - vbot)
	      else
		write(6,*) 'ihwadv = ',ihwadv
		stop 'error stop set_advective: ihwadv not supported'
              end if
	      wtop = wbot
	    end if

!	    ---------------------------------------------------------------
!	    total contribution
!	    ---------------------------------------------------------------

	    f = rlin * ruseterm
	    fxv(l,ie) = fxv(l,ie) + f * (xadv + zxadv)
	    fyv(l,ie) = fyv(l,ie) + f * (yadv + zyadv)
	  end do
	end do

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

        end

!******************************************************************

	subroutine momentum_advective_stability(rlin,rindex,astab)

! computes courant number of advective terms in momentum equation
!
! stability is computed for dt == 1

	use mod_internal
	use mod_geom_dynamic
	use mod_layer_thickness
	use mod_hydro
	use evgeom
	use levels
	use basin
	use shympi

	implicit none

	real rlin		   !factor for advective terms - normally 1
	real rindex		   !stability index (return)
	real astab(nlvdi,nel)      !stability matrix (return)

	integer ie,l,ii,k,lmax,lmin,iweg
	real cc,cmax
	real ut,vt
	real area,h,vol
	real b,c,f,ftot,r

	cmax = 0.
	!call compute_stability_stats(-1,cc)

	do ie=1,nel
	  area = 12. * ev(10,ie)
	  lmax = ilhv(ie)
	  lmin = jlhv(ie)
	  iweg = iwegv(ie)
	  do l=lmin,lmax

            h = hdenv(l,ie)
	    vol = area * h

	    r = rdistv(ie)
  	    ut = utlnv(l,ie)
  	    vt = vtlnv(l,ie)

	    ftot = 0.
            do ii=1,3
                k = nen3v(ii,ie)
                b = ev(3+ii,ie)
                c = ev(6+ii,ie)
                f = ut * b + vt * c
                if( f .lt. 0. ) ftot = ftot - f
            end do

	    cc = rlin*r*area*ftot/vol
	    if( iweg .gt. 0 ) cc = 0.	! dry element
	    astab(l,ie) = cc
	    cmax = max(cmax,cc)
	    !call compute_stability_stats(0,cc)

	  end do
	end do

	!rindex = shympi_max(cmax)
	rindex = cmax
	!call compute_stability_stats(1,cc)

	end

!******************************************************************

	subroutine compute_stability_stats(what,cc)

! computes histogram of stability of elements

	implicit none

	integer what
	real cc

	integer ndim,dbin
	parameter( ndim = 10 , dbin = 10 )
	real eps
	parameter( eps = 1.e-5 )

	integer i,idt,it
	character*20 aline
	integer, save :: bin(0:ndim)

	if( what .lt. 0 ) then		!initialize
	  do i=0,ndim
	    bin(i) = 0
	  end do
	else if( what .eq. 0 ) then	!accumulate
	  if( cc .gt. eps ) then
	    idt = 1. / cc
	  else
	    idt = 9999999
	  end if
	  i = idt / dbin
	  i = max(i,0)
	  i = min(i,ndim)
	  bin(i) = bin(i) + 1
	else				!write out
	  call get_act_timeline(aline)
	  write(97,1000) aline,(bin(i),i=0,ndim)
	end if
	
	return
 1000	format(a20,11i6)
	end

!******************************************************************

	subroutine set_semi_lagrange

	use mod_internal
	use basin, only : nkn,nel,ngr,mbw

        implicit none
         
	integer ie,l
	real xadv,yadv,dt
        real uadv(nel),vadv(nel)

	call get_timestep(dt)

        call back_trace(uadv,vadv)

	l = 1			!only for one layer

	do ie=1,nel
          xadv = uadv(ie) / dt
          yadv = vadv(ie) / dt

	  fxv(l,ie) = fxv(l,ie) + xadv
	  fyv(l,ie) = fyv(l,ie) + yadv
	end do

	end

!******************************************************************

        subroutine set_barocl

	use mod_internal
	use mod_layer_thickness
	use mod_ts
	use evgeom
	use levels
	use basin
	use pkonst

        implicit none
         
        integer k,l,ie,ii			!was BUG
        real dt
        real rrho0
        real salref,temref,sstrat,tstrat
        real hlayer
        real hhi

        real xbcl,ybcl
        integer lmax

        real rhop,presbt,presbcx,presbcy,dprescx,dprescy,br,cr!dbf
        real b,c
	real, allocatable :: bpresxv(:,:), bpresyv(:,:)

        call get_timestep(dt)

	allocate(bpresxv(nlv,nel))
	allocate(bpresyv(nlv,nel))
	bpresxv = 0.
	bpresyv = 0.

        do ie=1,nel
            presbcx = 0.
            presbcy = 0.
            lmax=ilhv(ie)
            !print*,lmax,' lmax'
            do l=1,lmax            
                hlayer = 0.5 * hdeov(l,ie)
                
		br = 0.
		cr = 0.                 
                do ii=1,3                 
                        k = nen3v(ii,ie)
                        rhop = rhov(l,k) ! rho^prime for each node of element 
                        !print*,'rhov ', l,k,rhov(l,k)
                        b = ev(3+ii,ie)!gradient in x della funz di forma
                        c = ev(6+ii,ie)!gradient in y della funz di forma
                        br = br + (b*rhop) 
                        cr = cr + (c*rhop)
                end do
                presbcx = presbcx + br*hlayer
                presbcy = presbcy + cr*hlayer
                bpresxv(l,ie) = presbcx
                bpresyv(l,ie) = presbcy
                presbcx = presbcx + br*hlayer
                presbcy = presbcy + cr*hlayer
           end do
        end do
        
        rrho0=1./rowass
        !print*,'rowass ', rowass, grav,hhi 
        do ie=1,nel
            lmax=ilhv(ie)
                do l=1,lmax
                     hhi = hdeov(l,ie)
                     !print*, 'hhi ',hhi,bpresxv(l,ie),l,ie
                     xbcl =  rrho0*grav*hhi*bpresxv(l,ie)
                     ybcl =  rrho0*grav*hhi*bpresyv(l,ie)
                     fxv(l,ie) = fxv(l,ie) +xbcl
                     fyv(l,ie) = fyv(l,ie) +ybcl
                     !write(6,*)'fxv ',fxv,ie
                enddo
        enddo

        end

!**********************************************************************

        subroutine set_barocl_new

! computes baroclinic contribution centered on layers
!
! cannot use this for sigma levels

	use mod_internal
	use mod_layer_thickness
	use mod_ts
	use evgeom
	use levels
	use basin
	use pkonst

        implicit none
         
	logical bsigma
        integer k,l,ie,ii,lmax,lmin
        double precision hlayer,hhi
        double precision xbcl,ybcl
        double precision raux,rhop,presbcx,presbcy
        double precision b,c,br,cr

        raux=grav/rowass
	call get_bsigma(bsigma)

	if( bsigma ) then
	  stop 'error stop set_barocl_new: cannot use with sigma levels'
	end if

        do ie=1,nel
          presbcx = 0.
          presbcy = 0.
	  lmin = ilmv(ie)
          lmax = ilhv(ie)
          do l=1,lmax
            hhi = hdeov(l,ie)
            hhi = hldv(l)
            hlayer = 0.5 * hhi
                
	    br = 0.
	    cr = 0.                 
            do ii=1,3                 
              k = nen3v(ii,ie)
              rhop = rhov(l,k)		!rho^prime for each node of element 
              b = ev(3+ii,ie)		!gradient in x
              c = ev(6+ii,ie)		!gradient in y
              br = br + b * rhop
              cr = cr + c * rhop
            end do

            presbcx = presbcx + br * hlayer
            presbcy = presbcy + cr * hlayer

            xbcl =  raux * hhi * presbcx
            ybcl =  raux * hhi * presbcy
            fxv(l,ie) = fxv(l,ie) + xbcl
            fyv(l,ie) = fyv(l,ie) + ybcl

            presbcx = presbcx + br * hlayer
            presbcy = presbcy + cr * hlayer
          end do
        end do
        
        end

!**********************************************************************

        subroutine set_barocl_old

! do not use this routine !

	use mod_internal
	use mod_layer_thickness
	use mod_ts
	use evgeom
	use levels
	use basin
	use pkonst

        implicit none
         
        integer k,l,ie,ii,lmax,lmin
        double precision hlayer,hhi
        double precision xbcl,ybcl
        double precision raux,rhop,presbcx,presbcy
        double precision b,c,br,cr

	double precision px(0:nlvdi)
	double precision py(0:nlvdi)

	stop 'error stop set_barocl_old: do not use this routine'

        raux=grav/rowass

        do ie=1,nel
          presbcx = 0.
          presbcy = 0.
	  lmin = jlhv(ie)
          lmax = ilhv(ie)

	  px(0) = presbcx
	  py(0) = presbcy
          do l=lmin,lmax
            hhi = hdeov(l,ie)
            hhi = hldv(l)
            hlayer = 0.5 * hhi
            hlayer = hhi
                
	    br = 0.
	    cr = 0.                 
            do ii=1,3                 
              k = nen3v(ii,ie)
              rhop = rhov(l,k)		!rho^prime for each node of element 
              b = ev(3+ii,ie)		!gradient in x
              c = ev(6+ii,ie)		!gradient in y
              br = br + b * rhop
              cr = cr + c * rhop
            end do

            presbcx = presbcx + br * hlayer
            presbcy = presbcy + cr * hlayer
	    px(l) = presbcx
	    py(l) = presbcy

          end do

          do l=1,lmax
	    presbcx = 0.5*(px(l) + px(l-1))
	    presbcy = 0.5*(py(l) + py(l-1))
            xbcl =  raux * hhi * presbcx
            ybcl =  raux * hhi * presbcy
            fxv(l,ie) = fxv(l,ie) + xbcl
            fyv(l,ie) = fyv(l,ie) + ybcl
	  end do
        end do
        
        end

!**********************************************************************

        subroutine set_barocl_new_interface

! computes baroclinic contribution centered on interfaces
!
! this routine works with Z and sigma layers

	use mod_internal
	use mod_layer_thickness
	use mod_ts
	use mod_hydro
	use evgeom
	use levels
	use basin
	use shympi
	use shympi_debug
	use pkonst

        implicit none
         
!---------- DEB SIG
	real hkk
	!real hkko(0:nlvdi,nkn)	!depth of interface at node
	!real hkkom(0:nlvdi,nkn)	!average depth of layer at node
	real, allocatable :: hkko(:,:)	!depth of interface at node
	real, allocatable :: hkkom(:,:)	!average depth of layer at node
	!real, allocatable :: aux2d(:)	!temporary
	!real, allocatable :: aux3d(:,:)	!temporary
	real hele
	real helei
	real alpha,aux,bn,cn,bt,ct
	real h,hd,hu,rd,ru
	real brl,crl

	integer lkmax,ld,lu
	integer laux,ll,ls,nn,nb !DEB
	integer llup(3),lldown(3)
!---------- DEB SIG

	logical bdebug
	logical bsigma,bzadapt,badapt,bmoveinterface,bsigadjust
        integer k,l,ie,ii,lmax,lmin,nsigma,ie_mpi
	real hsigma,hdep,htot
        double precision hlayer,hint,hhk,hh,hhup,htint
	double precision dzdx,dzdy,zk
        double precision xbcl,ybcl
        double precision raux,rhop,presbcx,presbcy
        double precision b,c,br,cr,brup,crup,brint,crint
	double precision rhoup,psigma
	double precision b3,c3
	double precision rdist,rcomp,ruseterm
	double precision dtime

	integer ipext,ieext
	integer nzadapt
	integer nadapt(4)
	real hadapt(4)

	bsigadjust = .false.		!regular sigma coordinates
	bsigadjust = .true.		!interpolate on horizontal surfaces

	call get_act_dtime(dtime)

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0
	call get_nzadapt_info(nzadapt)
	bzadapt = nzadapt .gt. 1
	bmoveinterface = bsigma .or. bzadapt	!interfaces not aligned to geopotentials

        raux=grav/rowass
	psigma = 0.
	
	allocate(hkko(0:nlvdi,nkn))
	allocate(hkkom(0:nlvdi,nkn))
	!allocate(aux3d(nlvdi,nkn))
	!allocate(aux2d(nkn))

	if( bmoveinterface .and. bsigadjust ) then	!-------------- DEB SIG
	  do k=1,nkn
	    lmax=ilhkv(k)
	    hkko(0,k)=-zov(k)	!depth of interface on node
	    hkkom(0,k)=-zov(k)	!depth of mid layer on node (0 not used)
	    hkk=0.
	    hkk=-zov(k)		!ggu
	    do l=1,lmax
	      hkk=hkk+hdkov(l,k)
	      hkko(l,k)=hkk
	      hkkom(l,k)=(hkko(l,k)+hkko(l-1,k))/2.
            end do
	  end do
	end if
	 
        do ie_mpi=1,nel
	  ie = ip_sort_elem(ie_mpi)
          call get_zadapt_info(ie,nadapt,hadapt)
          rdist = rdistv(ie)              !use terms (distance from OB)
          rcomp = rcomputev(ie)           !use terms (custom elements)
          ruseterm = min(rcomp,rdist)     !use terms (both)

          presbcx = 0.
          presbcy = 0.
	  lmin = jlhv(ie)
          lmax = ilhv(ie)
	  brup=0.
	  crup=0.
	  hhup=0.
          do l=lmin,lmax		!loop over layers to set up interface l-1
	    bsigma = l .le. nsigma
	    badapt = l .le. (nadapt(4)+lmin-1)
	    bmoveinterface = bsigma .or. badapt

	    htint = 0.				!depth of layer top interface
	    if( l .gt. 1 ) htint = hlv(l-1)

            hlayer = hdeov(l,ie)		!layer thickness
	    if( (.not. bsigma) .and. (.not. badapt)) hlayer = hldv(l)

            hh = 0.5 * hlayer
	    hint = hh + hhup			!interface thickness
                
	    if( bmoveinterface .and. bsigadjust ) then	!-------------- DEB SIG
	      hele = 0.
	      helei = 0.
	      do ii=1,3
                k = nen3v(ii,ie)
	        hele=hele+hkko(l,k)+hkko(l-1,k)	!depth of mid layer in element
	        helei=helei+hkko(l-1,k)		!depth of interface in element
	      end do

	      hele=hele/6.			!depth of mid layer l
	      helei=helei/3.			!depth of interface l-1

	      do ii=1,3
                k = nen3v(ii,ie)   
		lkmax = ilhkv(k)
	        if(helei.lt.hkko(l-1,k))then	!look upwards
		  do ll=l-1,1,-1
	            if(helei.gt.hkko(ll-1,k)) exit
		  end do
		  if( ll .le. 0 ) ll = 1
                else if(helei.gt.hkko(l,k))then	!look downwards
		  do ll=l+1,lkmax
	            if(helei.lt.hkko(ll,k)) exit
		  end do
		  if( ll .gt. lkmax ) ll = lkmax
		else				!inside layer l
		  ll = l
	        end if
		!interface l-1 is inside layer ll
		if( helei.lt.hkkom(ll,k) ) then	!find part of layer (up or down)
		  llup(ii) = ll-1
		  if( ll .eq. 1 ) llup(ii) = 1
		  lldown(ii) = ll
		else
		  llup(ii) = ll
		  lldown(ii) = ll+1
		  if( ll .eq. lkmax ) lldown(ii) = lkmax
		end if
		!do final check just to be sure (may be commented)
		if( lldown(ii) .eq. 1 ) then
		  if( helei .gt. hkkom(1,k) ) goto 99
		else if( llup(ii) .eq. lkmax ) then
		  if( helei .lt. hkkom(lkmax,k) ) goto 99
		else
		  if( helei .gt. hkkom(lldown(ii),k) ) goto 99
		  if( helei .lt. hkkom(llup(ii),k) ) goto 99
		end if
	      end do
	    end if

	    nn = 0 
	    nb = 0
	    brl = 0.
	    crl = 0.                 
	    br = 0.
	    cr = 0.                 
	    dzdx = 0.
	    dzdy = 0.
	    psigma = 0.

            do ii=1,3
	      k = nen3v(ii,ie)
	      rhop = rhov(l,k)		!rho^prime for each node of element
	      rhoup = rhop
	      if( l.gt.1) rhoup = rhov(l-1,k)
	      lkmax = ilhkv(k)
              b = ev(3+ii,ie)		!gradient in x
              c = ev(6+ii,ie)		!gradient in y

	      if( l .eq. nsigma .or. l .eq. nadapt(4) ) then !last mov layer
	        brl = brl + b * rhop
	        crl = crl + c * rhop
	      end if

	      if( bmoveinterface .and. bsigadjust ) then 
		lu = llup(ii)
		ld = lldown(ii)
		if( ld .eq. 1 ) then		!above surface
		  rhop = rhov(1,k)
		else if( lu .ge. lkmax ) then	!below bottom
		  nb = nb + 1
		  nn = nn + ii
		  rhop = rhov(lkmax,k)
		else				!do interpolation
		  !hu = hkko(lu,k)
		  !hd = hkko(ld,k)
		  hu = hkkom(lu,k) !DEB
		  hd = hkkom(ld,k) !DEB
		  ru = rhov(lu,k)
		  rd = rhov(ld,k)
		  h = helei
		  alpha = (h-hu)/(hd-hu)
		  rhop = alpha*rd + (1.-alpha)*ru
		end if
	      end if

              br = br + b * rhop
              cr = cr + c * rhop

              if (bmoveinterface) then
	       if( bsigadjust ) then
		psigma = 0.
	       else
                psigma = psigma + (rhoup-rhop)/hint
		htot = hm3v(ii,ie) 
                if (bsigma) then 
		   hsigma = htot !min(htot,hsigma)	  !FIXME lrp: this should be the min(htot,hsigma)
		   hdep = hsigma + zov(k)
                   hhk = -htint * hdep
                   zk = -(hhk) !- zov(k))                 !FIXME lrp: add free-surface
		else        
		  if (l.le.(nadapt(ii)+lmin-1)) then 
 		    hdep = hadapt(ii) + zov(k)
		    if ( nadapt(ii).eq.lmax ) then 
		      hdep = htot + zov(k)   		  !bottom layer
		      hadapt(ii) = htot
	            end if  
                    hhk = htint/hadapt(ii) * hdep 
		    zk = 0.
		    if (l.gt.1) zk = -(hhk - zov(k))      !transform depth in z		    
		  else
		    zk = htint
		  end if
		end if 
                dzdx = dzdx + b * zk
                dzdy = dzdy + c * zk
	       end if
              end if
            end do

	    if( bmoveinterface .and. bsigadjust ) then 
              if(nb.eq.2)then
	        brint = brup
	        crint = crup
	      elseif(nb.eq.1)then
	        b3 = ev(3+nn,ie)
	        c3 = ev(6+nn,ie)
	        aux=1./(c3*c3+b3*b3)
	        bn = aux*(brup*b3+crup*c3)*b3
	        cn = aux*(brup*b3+crup*c3)*c3
	        bt = br - aux*(br*b3+cr*c3)*b3
	        ct = cr - aux*(br*b3+cr*c3)*c3
	        brint = bn + bt
	        crint = cn + ct
              else  
	        brint = br
	        crint = cr
	      end if
	    else			  !zeta layer
              if( l .eq. 1 ) then         !surface layer ... treat differently
                brint = br
                crint = cr
              else
                brint = 0.5*(br+brup)
                crint = 0.5*(cr+crup)
              end if
	    end if

	    brup=br
	    crup=cr
	    if( l .eq. nsigma .or. l .eq. nadapt(4) ) then
	      brup=brl
	      crup=crl
	    end if
            hhup=hh
	    psigma = psigma / 3.

            presbcx = presbcx + hint * ( brint - dzdx * psigma )
	    presbcy = presbcy + hint * ( crint - dzdy * psigma )

            xbcl =  ruseterm * raux * hlayer * presbcx
            ybcl =  ruseterm * raux * hlayer * presbcy

            fxv(l,ie) = fxv(l,ie) + xbcl 
            fyv(l,ie) = fyv(l,ie) + ybcl
          end do
        end do
        
	deallocate(hkko)
	deallocate(hkkom)

	return
   99	continue
	write(6,*) ie,k,ii
	write(6,*) llup(ii),lldown(ii)
	write(6,*) helei,hkkom(lldown(ii),k),hkkom(llup(ii),k)
	stop 'error stop set_barocl_new_interface: internal error'
        end

!**********************************************************************
!**********************************************************************
!**********************************************************************

