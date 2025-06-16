
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-2000,2003-2006,2008-2020  Georg Umgiesser
!    Copyright (C) 2007-2008,2014  Christian Ferrarin
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

! routines for initialization
!
! contents :
!
! subroutine inicfil(name,var,nvar)
!				initializes nodal value variable from file
! subroutine inic2fil(name,var,nvar)
!				initializes nodal value variable from file (2D)
!
! revision log :
!
! 20.08.1998	ggu	initialization routines copied from subn11/new11
! 20.08.1998	ggu	new routine inicfil to init nodal values from file
! 06.11.1998	ggu	computation of hkv/hev taken out of sp211 -> cstset
! 22.10.1999	ggu	igrv, ngrv eliminated (subst by ilinkv, lenkv)
! 07.03.2000	ggu	arrays for vertical turbulent diffusion coefficient
! 07.03.2000	ggu	useless parts commented
! 20.06.2000	ggu	useless parts deleted
! 21.06.2000	ggu	dzreg introduced
! 08.08.2000	ggu	hlvmin is now percentage of last layer thickness
! 25.03.2003	ggu	voldist adjusted for 3D, new routine surinl
! 07.08.2003	ggu	deleted sp212
! 09.08.2003	ggu	cleaned up (sp211 is completely different)
! 14.08.2003	ggu	error in check_ilevels (integer instead of real)
! 14.08.2003	ggu	error in check_ilevels (must test .lt.)
! 14.08.2003	ggu	new routine set_depth -> transferred from cstcheck
! 14.08.2003	ggu	renamed sp211 into init_vertical and init_others
! 14.08.2003	ggu	new routines init_z, init_const_z, init_file_z
! 14.08.2003	ggu	new routines init_uvt
! 02.09.2003	ggu	routine inicfil restructured for lmax=0
! 03.09.2003	ggu	some more checks for check_ilevels
! 04.03.2004	ggu	change in inicfil() for more variables
! 05.03.2004	ggu	LMAX - changed bound from 0 to 1
! 02.09.2004	ggu	rdrst transfered to subrst.f
! 15.03.2005	ggu	set_austausch() and check_austausch() eliminated
! 05.04.2005	ggu	minimum levels (for baroclinic terms) -> set_min_levels
! 28.11.2005	ggu	in set_depth use makehkv to compute hkv
! 16.02.2006	ggu	describe file format in inicfil()
! 15.11.2006	ggu	new parameters to construct streched vert. coordinates
! 27.08.2007	ccf	variable coriolis from spherical coordinates (isphe)
! 17.03.2008	ggu	better error message if missing levels
! 07.04.2008	ggu	deleted surinl, volinl, voldist
! 12.11.2008	ggu	set_last_layer() rewritten, new adjust_k_depth()
! 06.12.2008	ggu&ccf	small bug fix in set_last_layer()
! 21.01.2009	ggu	cleaned up, new read_scalar()
! 27.01.2009	ggu	mzreg and hlvmax deleted (not used)
! 24.03.2009	ggu	bug fix: in set_last_layer() do not adjust for 2D
! 21.04.2009	ggu	new routine inic2fil()
! 23.03.2010	ggu	changed v6.1.1
! 09.04.2010	ggu	changed v6.1.3
! 28.09.2010	ggu	bug fix in init_coriolis() for isphe=1
! 08.10.2010	ggu	bug fix in init_coriolis() -> ym not set for isphe=1
! 16.12.2010	ggu	big restructering for sigma levels (look for bsigma)
! 27.01.2011	ggu	changed VERS_6_1_17
! 21.02.2011	ggu	error check for dzreg, nsigma and levels
! 01.03.2011	ggu	changed VERS_6_1_20
! 23.03.2011	ggu	new routine adjust_spherical(), error check isphe
! 14.04.2011	ggu	changed VERS_6_1_22
! 24.08.2011	ggu	eliminated hbot for sigma due to run time error
! 25.10.2011	ggu	hlhv eliminated
! 04.11.2011	ggu	adapted for hybrid coordinates
! 10.11.2011	ggu	changed VERS_6_1_36
! 11.11.2011	ggu	bug fix in adjust_levels: for zeta levels set nlv
! 22.11.2011	ggu	changed VERS_6_1_37
! 12.12.2011	ggu	changed VERS_6_1_39
! 29.08.2012	ggu	changed VERS_6_1_56
! 10.05.2013	ggu	changed VERS_6_1_64
! 13.06.2013	ggu	changed VERS_6_1_65
! 23.08.2013	ggu	renamed adjust_k_depth() to make_hkv() -> to subdep
! 23.08.2013	ggu	depth routines to subdep.f
! 05.09.2013	ggu	in adjust_levels() allow for nlv==1
! 12.09.2013	ggu	changed VERS_6_1_67
! 25.10.2013	ggu	changed VERS_6_1_68
! 31.10.2014	ccf	initi_z0 for zos and zob
! 05.11.2014	ggu	changed VERS_7_0_5
! 05.12.2014	ggu	changed VERS_7_0_8
! 12.12.2014	ggu	changed VERS_7_0_9
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_2
! 19.01.2015	ggu	changed VERS_7_1_3
! 30.04.2015	ggu	changed VERS_7_1_9
! 21.05.2015	ggu	changed VERS_7_1_11
! 25.05.2015	ggu	file cleaned and prepared for module
! 05.06.2015	ggu	changed VERS_7_1_12
! 10.07.2015	ggu	changed VERS_7_1_50
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 10.10.2015	ggu	changed VERS_7_3_2
! 05.11.2015	ggu	can now initialize z,u,v from file
! 28.04.2016	ggu	changed VERS_7_5_9
! 30.05.2016	ggu	changes in set_last_layer(), possible bug fix ilytyp==1
! 28.06.2016	ggu	coriolis computation changed -> yc=(ymax-ymin)/2.
! 05.10.2016	ggu	changed VERS_7_5_19
! 05.12.2017	ggu	changed VERS_7_5_39
! 24.01.2018	ggu	changed VERS_7_5_41
! 03.04.2018	ggu	changed VERS_7_5_43
! 19.04.2018	ggu	changed VERS_7_5_45
! 11.05.2018	ggu	new exchange_vertical() to compute global nlv,hlv
! 06.07.2018	ggu	changed VERS_7_5_48
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 12.02.2020	ggu	better error messages in set_last_layer()
! 02.06.2021	ggu	call levels_reinit() changed to levels_hlv_reinit()
! 20.07.2021	ggu	test if file has been opened for velocities
! 21.10.2022	ggu	in init_vertical() bug fix - update sigma_info (GGUBS)
! 28.04.2023	ggu	possible nkn=nel bug flagged with GGU_NKN_NEL
! 03.03.2024	ggu	bug fix in init_file_uv() - reset np
! 10.05.2024	ggu	set spherical in bas file
!
! notes :
!
! for information on hybrid levels see adjust_levels()
!
!**************************************************************

	subroutine init_vertical

! set up time independent vertical vectors

	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

	implicit none

	integer nlv_est,nlv_read,nlv_final
	integer nlv_e,nlv_k
	integer nlvaux,nsigma
	real hmax,hsigma
	real, allocatable :: hlv_aux(:)

	write(6,*) 'setting up vertical structure'

!------------------------------------------------------------------
! sanity check
!------------------------------------------------------------------

	call get_hmax_global(hmax)
	write(6,*) 'maximum depth: ',hmax

	nlv_est = nlv
	call estimate_nlv(nlv_est,hmax)
	write(6,*) 'nlv,nlv_est,nlvdi: ',nlv,nlv_est,nlvdi

	call check_nlv

	if( nlv > 0 ) then
	  allocate(hlv_aux(nlv))
	  hlv_aux(1:nlv) = hlv(1:nlv)
	  call levels_hlv_init(0)
	end if
	call levels_init(nkn,nel,nlv_est)
	if( nlv > 0 ) then
	  hlv(1:nlv) = hlv_aux(1:nlv)
	  deallocate(hlv_aux)
	end if

!------------------------------------------------------------------
! levels read in from $levels section
!------------------------------------------------------------------

	call adjust_levels(hmax)	!sets hlv, hldv, nlv, sigma_info, etc.

!------------------------------------------------------------------
! set up layer vectors
!------------------------------------------------------------------

	call set_ilhv		!sets nlv, ilhv (elemental)
	call set_last_layer	!adjusts nlv, ilhv, hm3v
	call set_ilhkv		!sets ilhkv (nodal)
	call set_ilmkv		!sets ilmkv (nodal)
	call exchange_levels	!copies from other domains and sets nlv
	call set_ilmv		!sets ilmv (elemental)

!------------------------------------------------------------------
! compute final nlv
!------------------------------------------------------------------

	nlv_e = maxval(ilhv)
	nlv_k = maxval(ilhkv)
	nlv = max(nlv_e,nlv_k)

!------------------------------------------------------------------
! check data structure
!------------------------------------------------------------------

	call get_sigma_info(nlvaux,nsigma,hsigma) !to change nlv info !GGUBS
	call set_sigma_info(nlv,nsigma,hsigma)
	
	nlv_final = nlv
	call levels_hlv_reinit(nlv_final)

	call check_vertical
	call shympi_set_hlv(nlv,hlv)

	write(6,*) 'init_vertical: nlvdi = ',nlvdi,'  nlv = ',nlv

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

	end

!**************************************************************

	subroutine exchange_vertical(nlv,hlv)

! exchanges vertical structure between domains
!
! to be deleted...
! should be really deleted... is never called FIXME

	use shympi

	implicit none

	integer nlv
	real hlv(nlv)

	integer ia
	real, parameter :: flag = -1.35472E+10
	real, allocatable :: hlvs(:,:)

        nlv_global = shympi_max(nlv)

	allocate(hlv_global(nlv_global))
	allocate(hlvs(nlv_global,n_threads))

	hlv_global = flag
	hlv_global(1:nlv) = hlv
	call shympi_gather(hlv_global,hlvs)

	do ia=1,n_threads
	  if( hlvs(nlv_global,ia) /= flag ) exit
	end do
	if( ia > n_threads ) then
	  write(6,*) 'error setting global nlv'
	  write(6,*) nlv_global,nlv
	  write(6,*) hlvs
	  stop 'error stop exchange_vertical: global nlv'
	end if

	hlv_global = hlvs(:,ia)

	!write(6,*) 'exchange_vertical: global nlv set: ',nlv,nlv_global
	!write(6,*) hlv_global

	end

!**************************************************************

	subroutine init_others

! set up various arrays (coriolis, eddy)

	implicit none

!------------------------------------------------------------------
! set others
!------------------------------------------------------------------

	call init_coriolis
	call set_eddy

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

	end

!*****************************************************************

	subroutine set_eddy

! sets vertical eddy coefficient

	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,l
	real vistur,diftur

	real getpar

!------------------------------------------------------------------
! get parameters
!------------------------------------------------------------------

	vistur=getpar('vistur')
	diftur=getpar('diftur')

!------------------------------------------------------------------
! set eddy coefficient
!------------------------------------------------------------------

	visv=vistur
	difv=diftur

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

	end

!*****************************************************************

	subroutine check_eddy

! checks vertical eddy coefficient

	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,l
	real v,d

!------------------------------------------------------------------
! check eddy coefficient
!------------------------------------------------------------------

	do k=1,nkn
	  do l=0,nlv
	    v = visv(l,k)
	    d = difv(l,k)
	    if( v .lt. 0. .or. v .gt. 1.e+5 ) goto 99
	    if( d .lt. 0. .or. d .gt. 1.e+5 ) goto 99
	  end do
	end do

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

	return
   99	continue
	write(6,*) k,l,v,d
	stop 'error stop check_eddy: error in values'
	end

!*****************************************************************

	subroutine check_vertical

! checks arrays containing vertical structure

	use levels
	use shympi

	implicit none

	call check_nlv
	call check_hlv
	call check_levels
	call check_ilevels

        if(bmpi_debug) then
	  call shympi_check_2d_node(ilhkv,'ilhkv')
	  call shympi_check_2d_node(ilmkv,'ilmkv')
	  call shympi_check_2d_elem(ilhv,'ilhv')
	  call shympi_check_2d_elem(ilmv,'ilmv')
        end if

	end

!*****************************************************************

	subroutine check_hlv

	use levels

	implicit none

	integer l

	write(6,*) 'check_hlv: ',nlv,nlvdi
	write(6,'(5g14.6)') (hlv(l),l=1,nlv)

	end

!*****************************************************************

	subroutine check_nlv

! checks nlv and associated parameters

	use levels, only : nlvdi,nlv

	implicit none

	write(6,*) 'check_nlv : ',nlvdi,nlv

	if(nlv.gt.nlvdi) stop 'error stop check_nlv: level dimension'

	end

!*****************************************************************

	subroutine estimate_nlv(nlv_est,hmax)

! estimates maximum value for nlv

	use basin

	implicit none

	integer nlv_est		!nlv_read on entry, estimate on return
	real hmax

	integer ie,ii
	integer nsigma,nreg
	real hsigma,dzreg

	real getpar

	dzreg = getpar('dzreg')
	call get_sigma(nsigma,hsigma)

	nreg = 0
	if( dzreg > 0 ) nreg = hmax/dzreg

	!nlv_est = nlv_est + nsigma + nreg + 1
	nlv_est = nlv_est + nsigma + nreg
	nlv_est = max(nlv_est,1)

	end

!*****************************************************************

	subroutine adjust_levels(hmax)

! adjusts levels read in from $levels section
!
! creates hlv if not set
! from hlv creates hldv
! needs hm3v to compute hmax
!
! strategy:
!
! set hsigma < hmax to ask for hybrid levels
! set nsigma > 1 to ask for sigma coordinates
! set dzreg > 0. to ask for regular zeta levels
!
! only nsigma				only sigma levels
! only dzreg				only regular zeta levels
! nsigma, dzreg and hsigma		hybrid levels
!
! only sigma levels			only sigma levels
! only zeta levels			only zeta levels
! sigma levels, dzreg and hsigma	hybrid levels
! zeta levels, nsigma and hsigma	hybrid levels
! sigma and zeta levels and hsigma	hybrid levels

	use levels
	use basin

	implicit none

	real hmax

	logical bsigma,bhybrid,bzeta
	integer l,ie,ii,nsigma,second
	real dzreg,hl,fact,hsigma
	real hbot,htop

	real getpar

	write(6,*) 'adjust layer structure'

!--------------------------------------------------------------
! create hlv values
!--------------------------------------------------------------

	second = 2
	dzreg = getpar('dzreg')
	call get_sigma(nsigma,hsigma)
	write(6,*) 'nlv,nsigma,hsigma: ',nlv,nsigma,hsigma

	bsigma = nsigma .gt. 0
	bhybrid = hsigma .lt. hmax
	bzeta = dzreg .gt. 0.

	if( nsigma .eq. 1 ) goto 92

	if( nlv .le. 0 ) then		! hlv not set -> must set

	  if( bsigma ) then		!sigma layers
	    if( nsigma > nlvdi ) goto 86
	    call make_sigma_levels(nsigma,hlv)
	    nlv = nsigma
	  end if

	  if( bhybrid .or. bzeta ) then	!zeta levels
	    if( bhybrid .and. .not. bzeta ) goto 99
	    if( bhybrid .and. .not. bsigma ) goto 90
	    nlv = nlvdi
	    if( bhybrid ) then
	      call make_zeta_levels(nsigma,hsigma,dzreg,nlv,hlv)
	    else if( .not. bsigma ) then	!no sigma layers given
	      call make_zeta_levels(0,0.,dzreg,nlv,hlv)
	    else
	      goto 89
	    end if
	  else if( .not. bsigma ) then		!just one layer
	    nlv = 1
	    hlv(1) = hsigma
	  end if

	else				!level section present

	  if( nlv .eq. 1 ) then		!just one level (barotropic)
	    if( bsigma .or. bhybrid ) goto 98
	  else if( hlv(second) .gt. hlv(1) ) then	!zeta layers given
	    if( bzeta ) goto 97
	    if( bsigma ) then		!put sigma layers on top
	      if( .not. bhybrid ) goto 88
	      if( hlv(1) .ne. hsigma ) goto 96
	      if( nsigma + nlv - 1 .gt. nlvdi ) goto 95
	      do l=1,nlv
	        hlv(nsigma+l-1) = hlv(l)
	      end do
	      call make_sigma_levels(nsigma,hlv)
	      hlv(nsigma) = hsigma
	      nlv = nlv + nsigma - 1
	    end if
	  else if( hlv(nlv) .eq. -1. ) then	!only sigma layers given
	    if( bsigma ) goto 91
	    nsigma = nlv
	    if( bhybrid ) then		!put zeta levels below
	      if( .not. bzeta ) goto 94
	      nlv = nlvdi
	      call make_zeta_levels(nsigma,hsigma,dzreg,nlv,hlv)
	    end if
	  else				!both sigma and zeta levels given
	    if( bsigma ) goto 91
	    if( bzeta ) goto 97
	    if( .not. bhybrid ) goto 88
	    do l=2,nlv
	      if( hlv(l) .gt. hlv(l-1) ) exit
	    end do
	    if( l > nlv ) goto 93	!not found
	    if( hlv(l-1) .eq. -1. ) then
	      nsigma = l-1
	      hlv(l-1) = hsigma
	    else if( hlv(l) .eq. hsigma ) then
	      nsigma = l
	    else
	      goto 93
	    end if
	    if( hlv(nsigma) .eq. hlv(nsigma+1) ) then	!double entry
	      do l=nsigma+1,nlv
		hlv(l-1) = hlv(l)
	      end do
	      nlv = nlv - 1
	    else if( hlv(nsigma) .ge. hlv(nsigma+1) ) then
	      goto 93
	    end if
	  end if

	end if

	call set_sigma(nsigma,hsigma)		!uses putpar
	call set_sigma_info(nlv,nsigma,hsigma)	!sets internal data structure

!--------------------------------------------------------------
! create hldv values
!--------------------------------------------------------------

	hldv(1)=hlv(1)
	do l=2,nlv
	  htop = hlv(l-1)
	  hbot = hlv(l)
	  if( l .eq. nsigma ) hbot = -1.
	  hldv(l) = hbot - htop
	end do

	write(6,*) 'adjust_levels: '
	write(6,*) 'nlv,nsigma,hsigma: ',nlv,nsigma,hsigma
	write(6,'(5g14.6)') 'hlv:  ',(hlv(l),l=1,nlv)
	write(6,'(5g14.6)') 'hldv: ',(hldv(l),l=1,nlv)

!--------------------------------------------------------------
! check hlv and hldv values
!--------------------------------------------------------------

	hbot = hlv(nlv)
	if( hbot /= -1. .and. hmax > hbot ) goto 87

	call check_levels

	write(6,*) 'finished adjusting layer structure ',nlv

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	return
   86	continue
	write(6,*) 'nsigma,nlvdi: ',nsigma,nlvdi
	stop 'error stop adjust_levels: dimension too small'
   87	continue
	write(6,*) 'nlv,hlv(nlv),hmax: ',nlv,hlv(nlv),hmax
	stop 'error stop adjust_levels: not enough layers'
   88	continue
	write(6,*) 'hsigma: ',hsigma
	write(6,*) 'for hybrid levels hsigma must be set'
	stop 'error stop adjust_levels: hsigma'
   89	continue
	write(6,*) 'nsigma,dzreg: ',nsigma,dzreg
	write(6,*) 'with pure sigma layers cannot give dzreg'
	stop 'error stop adjust_levels: bsigma and dzreg'
   90	continue
	write(6,*) 'nsigma,hsigma: ',nsigma,hsigma
	write(6,*) 'with hsigma set nsigma must be > 0'
	stop 'error stop adjust_levels: nsigma and hsigma'
   91	continue
	write(6,*) 'nlv,nsigma: ',nlv,nsigma
	write(6,*) 'You cannot give both sigma levels and nsigma'
	stop 'error stop adjust_levels: nsigma'
   92	continue
	write(6,*) 'nsigma must be > 1'
	stop 'error stop adjust_levels: nsigma = 1'
   93	continue
	write(6,*) 'error in level structure'
	write(6,*) (hlv(l),l=1,nlv)
	stop 'error stop adjust_levels: hlv'
   94	continue
	write(6,*) 'for hybrid levels dzreg must be > 0'
	stop 'error stop adjust_levels: dzreg'
   95	continue
	write(6,*) 'nlv,nsigma,nlvdi: ',nlv,nsigma,nlvdi
	write(6,*) 'not enough space to add zeta levels'
	stop 'error stop adjust_levels: nlvdi'
   96	continue
	write(6,*) 'hlv(1),hsigma: ',hlv(1),hsigma
	write(6,*) 'for hybrid levels first zeta level must be hsigma'
	stop 'error stop adjust_levels: hsigma'
   97	continue
	write(6,*) 'nlv,dzreg: ',nlv,dzreg
	write(6,*) 'You cannot give both zeta levels and dzreg'
	stop 'error stop adjust_levels: dzreg'
   98	continue
	write(6,*) 'only one level given in level section'
	write(6,*) hlv(1)
	write(6,*) 'cannot have sigma levels...'
	stop 'error stop adjust_levels: nlv = 1'
   99	continue
	write(6,*) 'for hybrid levels dzreg must be > 0'
	stop 'error stop adjust_levels: dzreg = 0'
	end

!*****************************************************************

	subroutine check_levels

! checks arrays hlv and hldv

	use levels

	implicit none

! local
	logical bstop,bsigma
	integer l,nsigma,levmin
	real h,hd,fact,hsigma
	real hbot,htop

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	bstop = .false.

!--------------------------------------------------------------
! check hlv values
!--------------------------------------------------------------

	if( nsigma .gt. 0 ) then
	  h = hlv(nsigma)
	  if( h .ne. -1. .and. h .ne. hsigma ) then
	    write(6,*) h,hsigma
	    stop 'error stop check_levels: hsigma'
	  end if
	end if

	hbot = -hlv(1)
	do l=2,nsigma
	  htop = hbot
	  hbot = -hlv(l)
	  if( l .eq. nsigma ) hbot = 1.
	  if( hbot .le. htop ) then
	    write(6,*) 'Error in level values for level : ',l
	    write(6,*) '   hlv(l-1),hlv(l) :',hlv(l-1),hlv(l)
	    bstop = .true.
	  end if
	end do

	levmin = nsigma + 1
	if( levmin .eq. 1 ) levmin = 2
	do l=levmin,nlv
	  if( hlv(l) .le. hlv(l-1) ) then
	    write(6,*) 'Error in level values for level : ',l
	    write(6,*) '   hlv(l-1),hlv(l) :',hlv(l-1),hlv(l)
	    bstop = .true.
	  end if
	end do

	if( bstop ) stop 'error stop check_levels: level error'

!--------------------------------------------------------------
! check hldv values
!--------------------------------------------------------------

	hbot = hlv(1)
	do l=2,nlv
	  htop = hbot
	  hbot = hlv(l)
	  if( l .eq. nsigma ) hbot = -1.
	  hd=hbot-htop
	  if( hd .ne. hldv(l) ) then
	    write(6,*) 'Error in dlevel values for level : ',l
	    write(6,*) '   hd,hldv(l) :',hd,hldv(l)
	    bstop = .true.
	  end if
	  hbot = hlv(l)
	end do
	if( bstop ) stop 'error stop check_levels: dlevel error'

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*****************************************************************

	subroutine exchange_levels

! exchanges level info with other domains - sets nlv

	use levels
	use shympi

	implicit none

	!call shympi_comment('exchanging ilhkv, ilmkv')
	call shympi_exchange_2d_node(ilhkv)
	call shympi_exchange_2d_node(ilmkv)
	!call shympi_barrier

	end

!*****************************************************************

	subroutine set_ilhv

! sets nlv and ilhv - only needs hm3v and hlv, hev is temporary array

	use levels
	use basin

	implicit none

! local
	logical bsigma
	integer ie,ii,l,lmax,nsigma
	real hsigma

	real h,hmax,hm
	real hev(nel)		!local

	lmax=0
	hmax = 0.

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	do ie=1,nel
	  hm = 0.
	  do ii=1,3
	    hm = hm + hm3v(ii,ie)
	  end do
	  hm = hm / 3.
	  hev(ie) = hm
	  hmax = max(hmax,hm)
	end do

	do ie=1,nel

	  h=hev(ie)

	  if( bsigma .and. h .le. hsigma ) then	!only sigma levels
	    l = nsigma
	  else
	    do l=nsigma+1,nlv
	      if(hlv(l).ge.h) exit
	    end do
	    if( l .gt. nlv ) goto 99
	  end if

	  ilhv(ie)=l
	  lmax = max(lmax,l)

	end do

	nlv = lmax

	write(6,*) 'finished setting ilhv and nlv'
	write(6,*) 'nsigma,hsigma: ',nsigma,hsigma
	write(6,*) 'nlv,lmax,hmax: ',nlv,lmax,hmax
	write(6,'(5g14.6)') (hlv(l),l=1,nlv)

	return
   99	continue
	write(6,*) ie,l,nlv,h,hlv(nlv)
	write(6,*) 'maximum basin depth: ',hmax
	write(6,*) 'maximum layer depth: ',hlv(nlv)
	write(6,'(5g14.6)') (hlv(l),l=1,nlv)
	stop 'error stop set_ilhv: not enough layers'
	end

!*****************************************************************

	subroutine set_ilhkv

! set ilhkv array - only needs ilhv

	use levels
	use basin
	use shympi

	implicit none

	integer ie,ii,k,l

	do k=1,nkn
	  ilhkv(k)=0
	end do

	do ie=1,nel
	  l=ilhv(ie)
	  do ii=1,3
	    k=nen3v(ii,ie)
	    if(l.gt.ilhkv(k)) ilhkv(k)=l
	  end do
	end do

	if( shympi_partition_on_elements() ) then
          !call shympi_comment('shympi_elem: exchange ilhkv - max')
          call shympi_exchange_2d_nodes_max(ilhkv)
	else
          call shympi_exchange_2d_node(ilhkv)
	end if

	end

!*****************************************************************

	subroutine set_ilmkv

! set minimum number of levels for node

	use levels
	use basin
	use shympi

	implicit none

	integer ie,ii,k,l
	integer lmin,lmax

	ilmkv = huge(1)

	do ie=1,nel
	  l=ilhv(ie)
	  do ii=1,3
	    k=nen3v(ii,ie)
	    if(l.lt.ilmkv(k)) ilmkv(k)=l
	  end do
	end do

	if( shympi_partition_on_elements() ) then
          !call shympi_comment('shympi_elem: exchange ilmkv - min')
          call shympi_exchange_2d_nodes_min(ilmkv)
	else
          call shympi_exchange_2d_node(ilmkv)
	end if

	do k=1,nkn
	  lmin = ilmkv(k)
	  lmax = ilhkv(k)
	  if( lmin .gt. lmax .or. lmin .le. 0 ) then
	    stop 'error stop set_min_levels: lmin'
	  end if
	end do

	end

!*****************************************************************

	subroutine set_ilmv

! set minimum number of levels for elems

	use levels
	use basin

	implicit none

	integer ie,ii,k,lmin

	do ie=1,nel
	  lmin = huge(1)
	  do ii=1,3
	    k=nen3v(ii,ie)
	    if(lmin.gt.ilmkv(k)) lmin = ilmkv(k)
	  end do
	  ilmv(ie) = lmin
	end do

	end

!*****************************************************************

	subroutine check_ilevels

! checks arrays ilhv and ilhkv

	use levels
	use basin

	implicit none

	logical bsigma,bspure
	integer nsigma
	integer ie,ii,k,lmax,lk
	real hmax,hsigma

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	hmax = 0.
	do ie=1,nel
	  lmax = ilhv(ie)
	  if( lmax .le. 0 ) goto 99
	  do ii=1,3
	    hmax = max(hmax,hm3v(ii,ie))
	  end do
	end do

	bspure = bsigma .and. hmax .le. hsigma	!pure sigma coordinates

	do ie=1,nel
	  lmax=ilhv(ie)
	  do ii=1,3
	    k=nen3v(ii,ie)
	    lk = ilhkv(k)
	    if( lk .le. 0 ) goto 98
	    if( lk .lt. lmax ) goto 98
	  end do
	  if( bspure .and. lmax .ne. nsigma ) goto 96
	  if( bsigma .and. lmax .lt. nsigma ) goto 96
	end do

	do k=1,nkn
	  lmax=ilhkv(k)
	  if( lmax .le. 0 ) goto 97
	  if( bspure .and. lmax .ne. nsigma ) goto 96
	  if( bsigma .and. lmax .lt. nsigma ) goto 96
	end do

	return
   96	continue
	write(6,*) ie,k,lmax,nsigma
	stop 'error stop check_ilevels: error in vertical structure (4)'
   97	continue
	write(6,*) k,lmax
	stop 'error stop check_ilevels: error in vertical structure (3)'
   98	continue
	write(6,*) ie,lmax,k,lk
	stop 'error stop check_ilevels: error in vertical structure (2)'
   99	continue
	write(6,*) ie,lmax
	stop 'error stop check_ilevels: error in vertical structure (1)'
	end

!*****************************************************************

	subroutine set_last_layer

! sets last layer thickness
!
! adjusts nlv, hm3v, ilhv

	use levels
	use basin

	implicit none

	logical bwrite
	logical b2d,bsigma,binsigma
	logical bdepth,blayer
	integer ie,l,ii
	integer ihtot,lmax
	integer ilytyp
        integer ic1,ic2,ic3
	integer nsigma
	real h,hold,hnew,hlast
	real hlvmin
	real hmin,hmax,hsigma
	real hm
	!double precision hm

	integer ieext
	real getpar

	bwrite = .false.

!------------------------------------------------------------
! see if 2D application
!------------------------------------------------------------

	b2d = nlv .eq. 1			!2D application

!------------------------------------------------------------
! check if sigma levels
!------------------------------------------------------------

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

!	if( bsigma ) then		!sigma layers
!	  write(6,*) 'set_last_layer (nlv used) : ',nlv
!	  write(6,*) 'sigma layers detected'
!	  return
!	end if

!------------------------------------------------------------
! set hlvmin
!------------------------------------------------------------

!	ilytyp: 
!	  0=no adjustment  
!	  1=adjust to full layers (change depth)
!         2=adjust to full layers only if h<hlvmin (change depth)
!         3=add to last layer if h<hlvmin (keep depth but change layer)

	hlvmin = getpar('hlvmin')		!min percentage of last layer
	ilytyp = nint(getpar('ilytyp'))

	if( hlvmin .gt. 1. .or. hlvmin .lt. 0. ) goto 98

!------------------------------------------------------------
! adjust last layer thickness
!------------------------------------------------------------

	ic1 = 0
	ic2 = 0
	ic3 = 0

	ihtot = 0
	lmax = 0

	do ie=1,nel

	  hm = sum(hm3v(:,ie))/3.

	  l = ilhv(ie)
	  binsigma = l .le. nsigma
	  hold = hm				!original depth of element
	  hnew = hold
	  h = hold
	  if( l .gt. 1 ) h = h - hlv(l-1)	!actual last layer thickness
	  if( binsigma ) h = -hold*hldv(l)	!sigma layer
	  hlast = hldv(l)			!regular last layer thickness
	  hmin = hlvmin * hlast

	  if( l .gt. 1 .and. h .le. 0. ) goto 99

	  bdepth = .false.			!adjust depth? (implies blayer)
	  blayer = .false.			!adjust layer?

	  if( l .gt. 1 ) then			!only for more than 1 layer
	    if( ilytyp .eq. 0 ) then
!		no adjustment
	    else if( ilytyp .eq. 1 ) then
	      if( h < hlast ) then		!last layer not full layer
		bdepth = .true.
	      end if
	    else if( ilytyp .eq. 2 ) then
	      if( h < hmin ) then		!last layer too small
		bdepth = .true.
	      end if
	    else if( ilytyp .eq. 3 ) then
	      if( h < hmin ) then		!add to layer above
		blayer = .true.
	      end if
	    else
	      write(6,*) 'ilytyp = ',ilytyp
	      stop 'error stop set_last_layer: internal error (9)'
	    end if
	  end if

	  if( b2d .or. binsigma ) then		!dont for 2D or in sigma layer
	    bdepth = .false.
	    blayer = .false.
	  end if

	  if( blayer .or. bdepth ) l = l -1

	  if( bdepth ) then		!depth and layer have changed
	    h = hldv(l)
	    hnew = hlv(l)
	    hm3v(:,ie) = hnew
	  else if( blayer ) then	!only layer has changed
	    h = hldv(l) + h
	  end if

	  ilhv(ie) = l

	  if( blayer .or. bdepth ) then		!write to terminal
	    if( bwrite ) then
	      if( ihtot .eq. 0 ) then
		  write(6,*) 'set_last_layer: Adjustment of depth values'
	      end if
	      write(6,'(2i8,3f12.2)') ie,ieext(ie),l,hnew,hold,h
	    end if
	    ihtot = ihtot + 1
	  end if

	  hlast = h			!finally set last layer thickness
	  lmax = max(lmax,l)

	  if( hlv(l) .ne. hnew ) ic1 = ic1 + 1
	  if( hldv(l) .ne. hlast ) ic2 = ic2 + 1
	  if( l .eq. 1 ) ic3 = ic3 + 1
	end do

!------------------------------------------------------------
! adjust nlv and write final statistics
!------------------------------------------------------------

	nlv=lmax

	write(6,*) 'set_last_layer (nlv used) : ',nlv
	write(6,*) 'Total number of elements with adjusted depth: ' &
     &			,ihtot
	write(6,*) 'Incomplete depth:     ',ic1
	write(6,*) 'Differing last layer: ',ic2
	write(6,*) 'One layer elements:   ',ic3

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	return
   98	continue
	write(6,*) '*** error setting up layer structure...'
	write(6,*) 'hlvmin is given in fraction of last layer'
	write(6,*) 'therefore:  0 <= hlvmin <= 1'
	write(6,*) 'hlvmin = ',hlvmin
	stop 'error stop set_last_layer: hlvmin'
   99	continue
	write(6,*) '*** error setting up layer structure...'
	write(6,*) 'ie,l,hold,h: ',ie,l,hold,h
	write(6,*) 'nsigma,hsigma: ',nsigma,hsigma
	write(6,*) 'hlv: ',nlv,(hlv(l),l=1,nlv)
	if( nsigma > 0 .and. hold <= 0. ) then
	  write(6,*) 'total depth: ',hold
	  write(6,*) 'cannot yet handle salt marshes with sigma layers'
	end if
	stop 'error stop set_last_layer: layer depth <=0 (internal error)'
	end

!*****************************************************************

	subroutine init_coriolis

! sets coriolis parameter

	use mod_internal
	use basin
        use coordinates
	use shympi
	use pkonst
	use mkonst

	implicit none

	real omega2	!double frequency of earth rotation
	parameter ( omega2 = 2.0 * 0.729E-4 )
	real rearth	!radius of earth
	parameter ( rearth = 6371000. )

	logical bgeo
	integer k,ie,ii,ie_mpi,iu
	integer icor
	integer isphe
	real yc,ym,y,ymin,ymax,dlat
	real aux1,aux2,rad
	real getpar
	real, dimension(nkn)	:: yaux

! if coordinates are cartesian (isphe=0) then icor determines 
! how Coriolis is used:
!
! 0:	no Coriolis (default)
! 1:	beta-plane (linear varying Coriolis parameter)
! 2:	f-plane (constant Coriolis parameter)
!
! please note that in case of icor=1 or 2 also the parameter dlat (the average 
! latitude of the basin) has to be set
!
! with spherical coordinates (isphe=1) the default is to always use
! Coriolis with latitude read from basin. If you really do not want to 
! use Coriolis, then please set icor = -1. The parameter dlat is not needed.

! with cartesian coordinates (isphe=0) the default is to use a constant 
! latitude (dlat) for Coriolis (icor > 0). If you want a spatially varying 
! Coriolis parameter you have to convert the cartesian coordinates to 
! spherical setting the basin projection (iproj > 0)

	icor=nint(getpar('icor'))	!flag how to use Coriolis
	call get_coords_ev(isphe)

	rad = pi / 180.
	dlat = dcor			! average latitude
	fcorv = 0.

	if( icor < 0 ) return		! no coriolis

	bgeo = ( isphe .eq. 1 .or. iproj .ne. 0 )  !use geographical coords

	if( bgeo ) then
	  yaux = ygeov
	else
	  yaux = ygv
	end if
       
! next is handled differently - shympi FIXME - might break compatibility

	ymin = shympi_min(yaux)
	ymax = shympi_max(yaux)
	yc = (ymax+ymin)/2.
	!yc   = sum(yaux)/nkn

	if( bgeo ) dlat = yc		! get directly from basin

	!iu = 670 + my_id
	!write(iu,*) my_id,icor,bgeo,isphe
	!write(iu,'(i10,3f16.10)') my_id,ymin,ymax,yc

	aux1 = 0.
	aux2 = 0.

	if(icor.eq.1) then	!beta plane
		aux1 = omega2 * sin(dlat*rad)
		aux2 = omega2 * cos(dlat*rad) / rearth
	else if(icor.eq.2) then	!f plane
		aux1 = omega2 * sin(dlat*rad)
		aux2 = 0.
	end if

	fcor = aux1		! coriolis value of average latitude

	write(6,*) 'pi, dlat     : ',pi,dlat
	write(6,*) 'icor, fcor   : ',icor,fcor
	write(6,*) 'f_0, beta    : ',aux1,aux2
	write(6,*) 'yc,ymin,ymax : ',yc,ymin,ymax

	do ie_mpi=1,nel
          ie = ip_sort_elem(ie_mpi)

	  ym=0.
	  do ii=1,3
	    ym=ym+yaux(nen3v(ii,ie))
	  end do
	  ym=ym/3.
	  if( bgeo ) then			!spherical
	    !fcorv(ie) = omega2*cos(ym*rad)	!BUG
	    fcorv(ie) = omega2*sin(ym*rad)
	  else if( isphe .eq. 0 ) then		!cartesian
  	    fcorv(ie)=aux1+aux2*(ym-yc)
	  else
	    write(6,*) 'isphe = ',isphe
	    if( isphe .eq. -1 ) write(6,*) '...not initialized'
	    stop 'error stop init_coriolis: value for isphe not allowed'
	  end if
	end do

	if( bextra_exchange ) then
	  call shympi_exchange_2d_elem(fcorv)
	end if

	end

!*****************************************************************

	subroutine check_coriolis

! checks coriolis parameter

	use mod_internal
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ie
	real f,fmax

	fmax = 2.0 * 0.729E-4

	do ie=1,nel
	  f = fcorv(ie)
	  if( fmax - abs(f) .lt. 0. ) then
	    write(6,*) ie,f,fmax
	    stop 'error stop check_coriolis: f too big'
	  end if
	end do

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine inicfil(name,var,nvar)

! initializes nodal value variable from file

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) name		!name of variable
	real var(nlvdi,nkn,1)		!variable to set
        integer nvar

	integer it
	character*80 file

	call getfnm(name,file)
	if( file .eq. ' ' ) return	!nothing to initialize

	write(6,*) 'this call is not supported anymore...'
	write(6,*) name
	write(6,*) file
	write(6,*) 'please use tracer_file_init()'
	stop 'error stop inicfil: internal error - unsupported call'

	end

!*****************************************************************

	subroutine inic2fil(name,var,nvar)

! initializes nodal value variable from file (2D version)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) name		!name of variable
	real var(nkn,1)			!variable to set
        integer nvar

	integer it
	integer nlvdi
	character*80 file

	call getfnm(name,file)
	if( file .eq. ' ' ) return	!nothing to initialize

	write(6,*) 'this call is not supported anymore...'
	write(6,*) name
	write(6,*) file
	write(6,*) 'please use tracer_file_init()'
	stop 'error stop inic2fil: internal error - unsupported call'

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine init_z(zconst)

! initializes water levels (znv and zenv)

	implicit none

	real zconst		!constant z value to impose

	character*80 name
	logical rst_use_restart

!--------------------------------------------------------
! see if we already have hydro restart data
!--------------------------------------------------------

	if( rst_use_restart(1) ) return

!--------------------------------------------------------
! get name of file
!--------------------------------------------------------

        call getfnm('zinit',name)

!--------------------------------------------------------
! initialize from file or with constant
!--------------------------------------------------------

	if(name.ne.' ') then
	  call init_file_z(name)
	else
	  call init_const_z(zconst)
	end if

!--------------------------------------------------------
! transfer nodal values to element values
!--------------------------------------------------------

	call setzev

!--------------------------------------------------------
! end of routine
!--------------------------------------------------------

	end

!*****************************************************************

	subroutine init_const_z(zconst)

! initializes water level with constant

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real zconst		!constant z value to impose

	integer k,ie,ii

	znv = zconst
	zenv = zconst

        write(6,*) 'Water levels initialized with constant z = ',zconst

	end

!*******************************************************************

	subroutine init_file_z(name)

! initializes water level from file

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw
	use intp_fem_file

	implicit none

	character*(*) name	!file name

	integer nb,k
	integer nvar,nintp,np,lmax,ibc
	integer idzeta
	double precision dtime
	integer nodes(1)
	real zconst(1)
	character*10 what

	call get_first_dtime(dtime)
        nodes = 0
        nvar = 1
        nintp = 2
        np = nkn
        lmax = 0
        ibc = 0                         !no lateral boundary
        what = 'zeta init'
        zconst = 0.

	write(6,'(a)') 'Initializing water levels...'
        call iff_init(dtime,name,nvar,np,lmax,nintp &
     &                          ,nodes,zconst,idzeta)
        call iff_set_description(idzeta,ibc,what)

        lmax = 1
        call iff_read_and_interpolate(idzeta,dtime)
        call iff_time_interpolate(idzeta,dtime,1,np,lmax,znv)

	call iff_forget_file(idzeta)

	write(6,*) 'Initial water levels read from file : '
	write(6,*) trim(name)

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine init_uvt

! initializes transport in levels

	use mod_hydro

	implicit none

	character*80 name
	logical rst_use_restart

!--------------------------------------------------------
! see if we already have hydro restart data
!--------------------------------------------------------

	if( rst_use_restart(1) ) return

!--------------------------------------------------------
! get name of file
!--------------------------------------------------------

        call getfnm('uvinit',name)

!--------------------------------------------------------
! initialize from file or with constant
!--------------------------------------------------------

	if(name.ne.' ') then
	  call init_file_uv(name)
	  call vtot
	else
	  utlnv = 0.
	  vtlnv = 0.
	  call ttov
	end if

!--------------------------------------------------------
! end of routine
!--------------------------------------------------------

	end

!*******************************************************************

	subroutine init_file_uv(name)

! initializes water level from file

	use mod_hydro
	use mod_hydro_vel
	use mod_hydro_print
	use basin, only : nkn,nel,ngr,mbw
	use intp_fem_file
	use levels

	implicit none

	character*(*) name	!file name

	integer nb,k
	integer nvar,nintp,np,lmax,ibc
	integer ntype,iformat
	integer idvel
	double precision dtime
	integer nodes(1)
	real uvconst(2)
	character*10 what

	integer fem_file_regular

	call get_first_dtime(dtime)
        nodes = 0
        nvar = 2
        nintp = 2
        np = nkn		!velocities must be on node - relax later
        lmax = nlvdi
        ibc = 0                 !no lateral boundary
        what = 'uv init'
        uvconst = 0.

	write(6,*) 'Initializing velocities...'

	call fem_file_test_formatted(name,np,nvar,ntype,iformat)
	if( nvar == 0 ) then
	  write(6,*) 'error opening or reading file'
	  write(6,*) 'file name ',trim(name)
	  stop 'error stop init_file_uv: open or read'
	else if( nvar /= 2 ) then
	  write(6,*) 'expecting 2 variables, found ',nvar
	  write(6,*) 'read error from file ',trim(name)
	  stop 'error stop init_file_uv: nvar'
	end if
	if( fem_file_regular(ntype) > 0 ) then
	  !np = nel
	  write(6,*) 'np = ',np
	else if( np /= nkn .and. np /= nel ) then
	  write(6,*) 'unexpected value for np found ',np
	  write(6,*) 'possible values: ',nkn,nel
	  write(6,*) 'read error from file ',trim(name)
	  stop 'error stop init_file_uv: np'
	end if

        np = nkn		!velocities must be on node - relax later

        call iff_init(dtime,name,nvar,np,lmax,nintp &
     &                          ,nodes,uvconst,idvel)
        call iff_set_description(idvel,ibc,what)

        call iff_read_and_interpolate(idvel,dtime)
	if( np == nkn ) then	!possible GGU_NKN_NEL bug
          call iff_time_interpolate(idvel,dtime,1,np,lmax,uprv)
          call iff_time_interpolate(idvel,dtime,2,np,lmax,vprv)
	  call prtouv
	else
          call iff_time_interpolate(idvel,dtime,1,np,lmax,ulnv)
          call iff_time_interpolate(idvel,dtime,2,np,lmax,vlnv)
	end if

	call make_prvel

	call iff_forget_file(idvel)

	write(6,*) 'Initial velocities read from file : '
	write(6,*) trim(name)

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine init_z0

! initializes surface z0sk(k) and bottom z0bn(k) roughness 

	use mod_roughness

	implicit none

	z0sk = z0sini
	z0bk = z0bini

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine set_spherical

	implicit none

	integer isphe
	real getpar

	isphe = nint(getpar('isphe'))
	if( isphe == -1 ) return	!no info given

	call bas_set_spherical(isphe)
	call set_coords_ev(isphe)

	end

!*******************************************************************

	subroutine adjust_spherical

	implicit none

	integer isphe
	real rsphe

	call get_coords_ev(isphe)
	rsphe = isphe
	call putpar('isphe',rsphe)

	end

!*******************************************************************

        subroutine print_spherical

	use shympi

        implicit none

        integer isphe

	if( .not. shympi_output() ) return

	call get_coords_ev(isphe)

        write(6,*) 'setting for coordinates: isphe = ',isphe
        if( isphe .eq. 0 ) then
          write(6,*) 'using cartesian coordinates'
        else
          write(6,*) 'using lat/lon coordinates'
        end if

        end

!*******************************************************************

