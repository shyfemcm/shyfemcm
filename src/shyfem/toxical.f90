
!--------------------------------------------------------------------------
!
!    Copyright (C) 2006,2008,2010,2012,2014-2015,2014-2015  Georg Umgiesser
!    Copyright (C) 2018-2019  Georg Umgiesser
!    Copyright (C) 2006  Francesca De Pascalis
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

! atoxi3d - Toxical routines from ARPAV - shell
!
! contents :
!
! revision log :
!
! 15.02.2006	ggu&fdp	new routine atoxi3d for ARPAV (from bio3d)
! 23.03.2006	ggu	ntot eliminated
! 23.03.2006	ggu	changed time step to real
! 17.04.2008	ggu	new open boundary conditions 
! 22.04.2008	ggu	advection parallelized
! 23.04.2008	ggu	call to bnds_set_def() changed
! 09.10.2008	ggu	new call to confop
! 23.03.2010	ggu	changed v6.1.1
! 30.03.2012	ggu	changed VERS_6_1_51
! 21.10.2014	ggu	converted to new boundary treatment
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 26.02.2015	ggu	changed VERS_7_1_5
! 17.07.2015	ggu	changed VERS_7_1_53
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 30.07.2015	ggu	changed VERS_7_1_83
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.03.2022	ggu	upgraded to da_out
! 22.03.2022	ggu	upgraded to new cmed
!
! notes :
!
! State variables used: (Wasp)
!
! todo :
!
! - * FIXME -> number of levels nlvdim, and not 1 (done)
! - wind speed and air temp is not available -> introduce (wmeteo)
!
! this routine has to be completely revised
!
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!
!***********************************************

	subroutine atoxi3d

! toxi module ARPAV

	use mod_diff_visc_fric
	use levels
	use basin
	use mkonst

	implicit none

	integer, parameter :: nlvdim = 1
	integer, parameter :: nkndim = 1

	integer nstate
	parameter( nstate = 1 )

	real e(nlvdim,nkndim,nstate)	!state vector
	real eb(nlvdim,nkndim,nstate)	!boundary vector of state vectors
	save e,eb

        real tstot(nstate)              !for mass test

        character*10 what

	integer k,i,l,lmax
	integer itoxi
        integer id
	integer nintp,ivar
	integer nbc
	real t,s
	real u,v
	real dt		!time step in seconds

	integer, save, allocatable :: idtoxi(:)

	real eaux(nstate)
	real einit(nstate)
	real ebound(nstate)
	save einit,ebound

	integer iunit
	integer j
	real rlux,rluxaux,itot,fday
	real dtt,dttday
	real area,vol
	real oxysat
	real getpar
	integer iround
	integer ieint,ipint
	integer nvar
	double precision dtime0,dtime

        integer mode
        real ai,lsurf

	logical bcheck
	logical bsedim
	integer ie,ii
	integer kspec
	integer idc,is
	real d
	real cbod,nh3,krear,sod
	real vel
	real windspeed,tempair
	real tday,t0,tsec
	real stp
        real mass
	real wsink

	integer nbnds

	double precision, save :: da_out(4)

	logical has_output_d,next_output_d

	integer, save :: icall = 0
	integer iespecial,inspecial
	save iespecial,inspecial
	real rkpar,difmol
	save rkpar,difmol

!------------------------------------------------------------------
!	initial and boundary conditions  [mg/l]			??
!
! 	 data einit /0.0, 0., 0.0, 0.0, 0.,   0.,0.,0.0,0.0/
 	 data einit /0.0/
 	 data ebound /0.0/
!
!------------------------------------------------------------------

        what = 'toxi'

!-------------------------------------------------------------------
! initialization
!-------------------------------------------------------------------

	if( icall .le. -1 ) return

	if( icall .eq. 0 ) then
	  itoxi = iround(getpar('itoxi'))
	  if( itoxi .le. 0 ) icall = -1
	  if( icall .le. -1 ) return
	  icall = 1

	  stop 'error stop atoxi3d: not adapted yet for new framework'

!         --------------------------------------------------
!	  initialize state variables with einit
!         --------------------------------------------------

	  do k=1,nkn		!loop on nodes
            lmax = ilhkv(k)
            do l=1,lmax
	      do is=1,nstate
	        e(l,k,is) = einit(is)
              end do
	    end do
          end do

!         --------------------------------------------------
!	  initialize state variables from external file
!         --------------------------------------------------

          call inicfil('toxi',e,nstate)

!         --------------------------------------------------
!	  set boundary conditions for all state variables
!         --------------------------------------------------

          nbc = nbnds()
          allocate(idtoxi(nbc))
          idtoxi = 0

	  call get_first_dtime(dtime0)
	  nintp = 2
	  nvar = nstate
          call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv &
     &                          ,ebound,idtoxi)

	  !call bnds_init(what,tox3dn,nintp,nstate,nb3dim,toxiarr,ebound)
	  !call bnds_set_def(what,nb3dim,toxiarr) !nvar != 1
	  !call bnds_print('init of '//what,nb3dim,toxiarr)

!         --------------------------------------------------
!	  initialize eco model
!         --------------------------------------------------

	  call atoxi_ini

!         --------------------------------------------------
!	  parameters for transport/diffusion resolution
!         --------------------------------------------------

          rkpar=getpar('chpar')
          difmol=getpar('difmol')

!         --------------------------------------------------
!	  initialize output 
!         --------------------------------------------------

          call init_output_d('itmcon','idtcon',da_out)
	  !call increase_output_d(da_out)
          if( has_output_d(da_out) ) then
            nvar = nstate
            call shyfem_init_scalar_file('tox',nvar,.false.,id)
            da_out(4) = id
          end if

	  write(6,*) 'toxi model initialized...'

	end if

!-------------------------------------------------------------------
! normal call
!-------------------------------------------------------------------

	kspec = -100
	!kspec = 930
	bcheck = .true.
	bcheck = .false.
	wsink = 0.

!	-------------------------------------------------------------------
!	time management
!	-------------------------------------------------------------------

	t0 = 0.
	call get_act_dtime(dtime)
	call get_timestep(dt)
	tsec = dtime

!	-------------------------------------------------------------------
!	loop on nodes for biological reactor
!	-------------------------------------------------------------------

        mode = +1               !new time level for volume and depth

	do k=1,nkn		!loop on nodes

          lmax = ilhkv(k)
          !call getmeteo(k,tempair,windspeed)    !meteo FIXME
          !call wmeteo(tempair,windspeed)      !meteo FIXME
          rlux = 1.

          do l=1,lmax
            call dvanode(l,k,mode,d,vol,area)   !gets depth, volume and area
            call getts(l,k,t,s)                 !gets temp and salt
            call getuv(l,k,u,v)                 !gets velocities u/v
            vel = sqrt(u*u+v*v)

            id = 1000*k+l

	    do is=1,nstate
	      eaux(is) = e(l,k,is)
	    end do

	    !call atoxi(id,tsec,dt,d,t,eaux)

	    do is=1,nstate
	      e(l,k,is) = eaux(is)
	    end do
          end do

	end do

!	-------------------------------------------------------------------
!	advection and diffusion
!	-------------------------------------------------------------------

	if( bcheck ) call check_toxi('before advection',e)

	call bnds_read_new(what,idtoxi,dtime)

!$OMP PARALLEL PRIVATE(i)
!$OMP DO SCHEDULE(DYNAMIC)

	do is=1,nstate

          call scal_adv(what,is &
     &                          ,e(1,1,is),idtoxi &
     &                          ,rkpar,wsink &
     &                          ,difhv,difv,difmol)

          call tsmass (e(1,1,is),1,nlvdim,tstot(is)) !mass control

	end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL

	if( bcheck ) call check_toxi('after advection',e)
	!write(86,*) dtime,tstot(1)

!	-------------------------------------------------------------------
!	write of results (file BIO)
!	-------------------------------------------------------------------

        if( next_output_d(da_out) ) then
          id = nint(da_out(4))
          idc = 120
	  do is=1,nstate
	    idc = idc + 1
            call shy_write_scalar_record(id,dtime,idc,nlvdi,e(:,:,is))
	  end do
        end if

!	call toxi_av_shell(e)		!aver/min/max of state vars

	if( bcheck ) call check_toxi('at end',e)

!	-------------------------------------------------------------------
!	debug output
!	-------------------------------------------------------------------

!        k = 100
!        l = 1
!        call getts(l,k,t,s)
!        call writeet(95,it,k,l,e,t,s,nlvdim,nkndim,nstate)

!	-------------------------------------------------------------------
!	end of routine
!	-------------------------------------------------------------------

	end

!*************************************************************

	subroutine writeet(iunit,it,k,l,e,t,s,nlvdim,nknddi,nstate)

! formatted write for debug

	implicit none

	integer iunit
	integer it,k,l
	integer nlvdim,nknddi,nstate
	real e(nlvdim,nknddi,nstate)
	real t
	real s

	integer i

	write(iunit,'(i10,11f12.4)') it, &
     &			(e(l,k,i),i=1,nstate), &
     &			t,s

	end

!*************************************************************

	subroutine toxi_av_shell(e)

! computes and writes average/min/max of bio variables
!
! id = 260
!
! e(1) average	== 261
! e(1) min	== 262
! e(1) max	== 263
! e(2) average	== 264
! ...

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

! parameter

	integer, parameter :: nstate = 1

	real e(nlvdi,nkn,nstate)	!state vector

! local
	integer itsmed
	integer id,nvar,idc,is
	double precision dtime
	double precision rr
! function
	real getpar
	logical has_output_d,is_over_output_d,next_output_d
! save
        double precision, save :: da_out(4)
        double precision, save, allocatable :: bioacu(:,:,:)
        real, save, allocatable :: biomin(:,:,:)
        real, save, allocatable :: biomax(:,:,:)
        real, save, allocatable :: raux(:,:)

	integer, save :: icall = 0
	integer, save :: nr = 0

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
            nvar = nstate*3
            call shyfem_init_scalar_file('bav',nvar,.false.,id)
            da_out(4) = id
          end if

          allocate(biomin(nlvdi,nkn,nstate),biomax(nlvdi,nkn,nstate))
	  allocate(bioacu(nlvdi,nkn,nstate))
          allocate(raux(nlvdi,nkn))

	  do is=1,nstate
            call cmed_reset(nr,bioacu(:,:,is) &
     &			,biomin(:,:,is),biomax(:,:,is))
	  end do

	  icall = 1
	end if

        if( .not. is_over_output_d(da_out) ) return

        nr = nr + 1
	do is=1,nstate
          call cmed_accum(e(:,:,is),bioacu(:,:,is) &
     &			,biomin(:,:,is),biomax(:,:,is))
	end do

        if( .not. next_output_d(da_out) ) return

        id = nint(da_out(4))
        call get_act_dtime(dtime)
        rr=1./nr

        idc = 260
	do is=1,nstate
	  idc = idc + 1
          raux = bioacu(:,:,is) * rr
          call shy_write_scalar_record(id,dtime,idc,nlvdi,raux)
          raux = biomin(:,:,is)
          call shy_write_scalar_record(id,dtime,idc,nlvdi,raux)
          raux = biomax(:,:,is)
          call shy_write_scalar_record(id,dtime,idc,nlvdi,raux)
	end do

	do is=1,nstate
          call cmed_reset(nr,bioacu(:,:,is) &
     &			,biomin(:,:,is),biomax(:,:,is))
	end do

	end

!*************************************************************

	subroutine check_toxi(title,e)

! checks bio vars

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer, parameter :: nlvdim = 1
	integer, parameter :: nkndim = 1

	integer nstate
	parameter( nstate = 1 )

	character*(*) title
	real e(nlvdim,nkndim,nstate)	!state vector


        character*20 text
	integer i

	write(6,*) 'check_toxi: ',title

        text = '*** bio check e     '
	do i=1,nstate
          write(text(18:19),'(i2)') i
          call check2Dr(nlvdim,nlv,nkn,e(1,1,i),0.,1.e+20,text,title)
	end do

	end

!*************************************************************

