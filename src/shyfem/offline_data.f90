
!--------------------------------------------------------------------------
!
!    Copyright (C) 2013-2015,2017-2020  Georg Umgiesser
!    Copyright (C) 2015  Christian Ferrarin
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

! routines for offline data handling
!
! revision log :
!
! 13.06.2013	ggu	new routines written from scratch
! 17.06.2013	ggu	eliminated compiler warnings
! 25.03.2014	ggu	new offline (for T/S)
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 01.04.2015	ggu	changed VERS_7_1_7
! 06.05.2015	ccf	write offline to .off file
! 06.05.2015	ccf	read offline from offlin file in section name
! 21.05.2015	ggu	changed VERS_7_1_11
! 10.07.2015	ggu	changed VERS_7_1_50
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 05.11.2015	ggu	revisited and checked
! 09.11.2015	ggu	changed VERS_7_3_13
! 16.11.2015	ggu	changed VERS_7_3_14
! 29.03.2017	ggu	bug fix - input file opened on unit 1
! 05.12.2017	ggu	changed VERS_7_5_39
! 12.11.2018	ggu	linear arrays introduced
! 18.12.2018	ggu	changed VERS_7_5_52
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 04.07.2019	ggu	solved problem for vertical velocity (1 index off)
! 17.02.2020	ggu	femtime eliminated
! 28.04.2020	ggu	routines dealing with records in new file mod_offline.f
! 20.03.2022	ggu	upgraded to suboutputd.f
! 15.03.2024	ggu	time to iitime; prepared for turbulence (bturb)
! 10.04.2024	ggu	finished offline
!
! notes :
!
!	uses still time as integer -> pass to double
!
!****************************************************************

	subroutine handle_offline(mode)

! handles offline version

!-----------------------------------------------------
!
! parameters:
!
!	mode - input parameter
!
!	mode = 1	write to file
!	mode = 2	read from file
!
!	idtoff - parameter set in STR file
!
!	idtoff = 0	nothing (no offline routines called)
!	idtoff > 0	write offline data file with time step idtoff
!	idtoff < 0	reads offline data from file
!	idtoff = -1	uses offline hydro results
!	idtoff = -2	uses offline T/S results
!	idtoff = -4	uses offline turbulence results
!
! combinations are possible: -3,-7
! everything is used by idtoff = -7
!
!-----------------------------------------------------

	use levels, only : nlvdi,nlv,ilhv,ilhkv
	use basin, only : nkn,nel,ngr,mbw
	use mod_offline
	use shympi

	implicit none

	integer mode

	logical bwrite
	integer itstart,it,itime_first
	integer ierr,ig,iu
	double precision dtime
	real dt
	character*60 name
	double precision dtoff,tmoff,toff
	integer, save :: iwhat = 0
        integer ifemop, ifileo
	real getpar

	if( icall .lt. 0 ) return

	bwrite = shympi_is_master()

!-------------------------------------------------------------
! initialize
!-------------------------------------------------------------

	call get_act_dtime(dtime)
	it = nint(dtime)

	if( bfirst ) then
	  ioffline = 0

	  ! this only works with small time values (integer)
          call convert_date_d('itmoff',tmoff)
          call convert_time_d('idtoff',dtoff)
	  call adjust_itmidt_d(tmoff,dtoff,toff)
	  itmoff = nint(tmoff)
	  idtoff = nint(dtoff)
	  itoff = nint(toff)
	  itoff = itmoff + idtoff

	  if( bwrite ) then
	    write(6,*) 'offline init:',itmoff,idtoff,it,itoff
	  end if

  	  if( idtoff .eq. 0 ) iwhat = 0		!nothing
	  if( idtoff .gt. 0 ) iwhat = 1		!write
	  if( idtoff .lt. 0 ) iwhat = 2		!read

          !if( iwhat == 1 .and. .not. has_output_d(da_out) ) iwhat = 0

	  if( iwhat .le. 0 ) icall = -1
	  if( idtoff .eq. 0 ) icall = -1
	  if( icall .lt. 0 ) return

	  ioff_mode = iwhat

	  !if( shympi_is_parallel() ) then
	  !  stop 'error stop offline: not ready for mpi'
	  !end if

	  call mod_offline_init(nkn,nel,nlvdi)
	  call off_init_vertical(nkn,nel,ilhv,ilhkv)
	  call off_reset

	  if( iwhat .eq. 1 ) then		! writing offline
            iu = ifemop('.off','unform','new')
            if( iu .le. 0 ) then
              write(6,*) 'iu = ',iu
              stop 'error stop handle_offline: cannot open output file'
            end if
	    iuoff = iu
	    if( bwrite ) write(6,*) 'Start writing offline file'
	  else					! reading offline
	    !if( shympi_is_parallel() ) then
	    !  stop 'error stop offline: reading not ready for mpi'
	    !end if
            call getfnm('offlin',name)
	    if( name == ' ' ) then
              write(6,*) '*** No offline file given'
              write(6,*) '*** please specify offlin in section $name'
              stop 'error stop handle_offline: cannot open input file'
	    end if
            iu = ifileo(0,name,'unformatted','old')
            if( iu .le. 0 ) then
              write(6,*) '*** Cannot find offline file: '
              write(6,*) trim(name)
              stop 'error stop handle_offline: cannot open input file'
            end if
	    iuoff = iu
	    if( bwrite ) then
              write(6,*) '--------------------------------------'
              write(6,*) 'reading offline from file: ',trim(name)
              write(6,*) '--------------------------------------'
	    end if
	  end if

	  ioffline = -idtoff

	  if( ioffline >= 4 ) then
	    write(6,*) 'warning offline'
	    write(6,*) 'requested reading of turbulence'
	    write(6,*) 'ioffline = ',ioffline
	    write(6,*) 'cannot yet do turbulence from offline'
	    !stop 'error stop handle_offline: no turbulence'
	    ioffline = ioff_max
	  end if
	  
	  bfirst = .false.
	end if

	if( it < itmoff ) return

!-------------------------------------------------------------
! do different modes (iwhat is 1 for write and 2 for read)
!-------------------------------------------------------------

	if( mode .ne. iwhat ) return

	if( mode .eq. 1 ) then			! writing

!	  -------------------------------------------------------------
!	  accumulate and write data
!	  -------------------------------------------------------------

	  call get_timestep(dt)
	  call off_accum(dt)

	  if( icall .eq. 0 ) then	!write first record
	    call off_aver
	    call off_write(iuoff,itmoff)
	    call off_reset
	    icall = 1
	  end if

	  if( it .lt. itoff ) return

	  call off_aver
	  call off_write(iuoff,it)
	  call off_reset
	  itoff = itoff + idtoff

	else if( mode .eq. 2 ) then		! reading

!	  -------------------------------------------------------------
!	  read data and put into hydro structures
!	  -------------------------------------------------------------

	  if( icall .eq. 0 ) then
	    do ig=1,nintp
	      call off_read(iuoff,ig,ierr)
	      if( ierr .ne. 0 ) goto 97
	    end do
	    call can_do_offline
	    itime_first = iitime(1)		!first time available in file
	    if( it .lt. itime_first ) goto 99
	    call get_timestep(dt)
	    if( it .eq. itmoff ) then
	      itstart = it
	    else
	      itstart = max(it-nint(dt),itmoff)
	    end if
	    call off_intp_all(iuoff,itstart)
	    icall = 1
	  end if

	  call off_intp_all(iuoff,it)

	  !call off_check(1)
	  !call off_check(2)

	else

!	  -------------------------------------------------------------
!	  error in mode
!	  -------------------------------------------------------------

	  write(6,*) 'mode = ',mode,'  iwhat = ',iwhat
	  stop 'error stop handle_offline: value for mode not allowed'

	end if

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	return
   97	continue
	write(6,*) iitime
	write(6,*) nintp,ig
	stop 'error stop handle_offline: read error at start'
   99	continue
	write(6,*) it,iitime
	stop 'error stop handle_offline: no time available'
	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine off_intp_all(iu,it)

	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use mod_diff_visc_fric
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use mod_offline
	use mod_trace_point

	implicit none

	integer iu
	integer it

	logical boff,bhydro,bts,bturb
	integer ierr
	integer ip,i,itnext
	integer ilhkw(nkn)

	integer, save :: ieof = 0

!	---------------------------------------------------------
!	initialize
!	---------------------------------------------------------

	ip = 2
	if( nintp .eq. 4 ) ip = 3

	call is_offline(1,bhydro)		!hydro
	call is_offline(2,bts)			!T/S
	call is_offline(4,bturb)		!turbulence

	!if( btrace ) write(6,*) 'trace offline: ',bhydro,bts,bturb

!	---------------------------------------------------------
!	find new records for time
!	---------------------------------------------------------

	do while( ieof .eq. 0 .and. it .gt. iitime(ip) )
	  call off_peek_next_record(iu,itnext,ieof)
	  if( ieof .ne. 0 ) exit
	  call off_copy
	  call off_read(iu,nintp,ierr)
	end do

	if( it .gt. iitime(nintp) ) goto 99

	!write(67,*) it,(iitime(i),i=1,nintp)
	!write(6,*) it,bhydro,bts,iwhat

!	---------------------------------------------------------
!	pre processing
!	---------------------------------------------------------

	if( bhydro ) then
	  call copy_uvz
	  call copy_depth
	end if

!	---------------------------------------------------------
!	interpolation
!	---------------------------------------------------------

	ilhkw = ilhkv + 1	!one more vertical value for wn

	if( bhydro ) then
	  !call trace_point('offline interpolating hydro')
	  call off_intp(nintp,it,iitime,nlvdi,nel,ilhv,nel,ut,utlnv)
	  call off_intp(nintp,it,iitime,nlvdi,nel,ilhv,nel,vt,vtlnv)
	  call off_intp(nintp,it,iitime,1,3*nel,ilhv,3*nel,ze,zenv)
	  call off_intp(nintp,it,iitime,nlvdi+1,nkn,ilhkw,nkn,wn,wlnv)
	  call off_intp(nintp,it,iitime,1,nkn,ilhkv,nkn,zn,znv)
	end if

	if( bts ) then
	  !call trace_point('offline interpolating T/S')
	  call off_intp(nintp,it,iitime,nlvdi,nkn,ilhkv,nkn,sn,saltv)
	  call off_intp(nintp,it,iitime,nlvdi,nkn,ilhkv,nkn,tn,tempv)
	end if

	if( bturb ) then
	  !call trace_point('offline interpolating turbulence')
	  call off_intp(nintp,it,iitime,nlvdi+1,nkn,ilhkw,nkn,vd,visv)
	  call off_intp(nintp,it,iitime,nlvdi+1,nkn,ilhkw,nkn,dd,difv)
	end if

!	---------------------------------------------------------
!	post processing
!	---------------------------------------------------------

	if( bhydro ) then
	  call make_new_depth
	  call uvint
          call ttov
          call make_prvel
	end if

	if( bts ) then
	  call rhoset_shell(1)
	end if

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	return
   99	continue
	write(6,*) 'time to interpolate: it = ',it
	write(6,*) 'time values available in time(): '
	write(6,*) (iitime(i),i=1,nintp)
	stop 'error stop off_intp_all: no such time'
	end

!****************************************************************

	subroutine off_intp(nintp,it,time,nlvddi,ndim,il,n,dval,rval)

	implicit none

	integer nintp
	integer it
	integer time(nintp)
	integer nlvddi,ndim
	integer il(ndim)
	integer n
	double precision dval(nlvddi,ndim,nintp)
	real rval(nlvddi,ndim)

	integer l,lmax,i,j
	real x(4),y(4),t

	real intp_neville

	if( nintp .lt. 2 .or. nintp .gt. 4 ) then
	  write(6,*) 'nintp = ',nintp
	  stop 'error stop off_intp: nintp not possible'
	end if

	t = it
	do j=1,nintp
	  x(j) = time(j)
	end do

	do i=1,n
	  lmax = 1
	  if( nlvddi .gt. 1 ) lmax = il(i)
	  do l=1,lmax
	    do j=1,nintp
	      y(j) = dval(l,i,j)
	    end do
	    rval(l,i) = intp_neville(nintp,x,y,t)
	  end do
	end do

	end

!****************************************************************

	subroutine off_intp4(it)

! not used

	use mod_hydro_vel
	use mod_hydro
	!use levels
	!use basin, only : nkn,nel,ngr,mbw
	use mod_offline

	implicit none

	integer it

	integer ie,ii,k,l,lmax,i,nintpol
	real x(4),y(4),t

	real intp_neville

	nintpol = 4

	t = it
	do i=1,nintpol
	  x(i) = iitime(i)
	end do
	
	do ie=1,nel_off
	  lmax = ile(ie)
	  do l=1,lmax
	    do i=1,nintpol
	      y(i) = ut(l,ie,i)
	    end do
	    utlnv(l,ie) = intp_neville(nintpol,x,y,t)
	    do i=1,nintpol
	      y(i) = vt(l,ie,i)
	    end do
	    vtlnv(l,ie) = intp_neville(nintpol,x,y,t)
	  end do
	  do ii=1,3
	    do i=1,nintpol
	      y(i) = ze(ii,ie,i)
	    end do
	    zenv(ii,ie) = intp_neville(nintpol,x,y,t)
	  end do
	end do

	do k=1,nkn_off
	  lmax = ilk(k)
	  do l=0,lmax
	    do i=1,nintpol
	      y(i) = wn(l,k,i)
	    end do
	    wlnv(l,k) = intp_neville(nintpol,x,y,t)
	  end do
	  do i=1,nintpol
	    y(i) = zn(k,i)
	  end do
	  znv(k) = intp_neville(nintpol,x,y,t)
	end do

	end

!****************************************************************

	subroutine off_intp2(it)

! not used

	use mod_hydro_vel
	use mod_hydro
	!use levels
	!use basin, only : nkn,nel,ngr,mbw
	use mod_offline

	implicit none

	integer it

	integer ie,ii,k,l,lmax
	integer it1,it2
	double precision rr,rt

	it1 = iitime(1)
	it2 = iitime(2)

	rr = 0.
	if( it2 .gt. it1 ) rr = float(it-it1)/float(it2-it1)
	rt = 1. - rr
	
	do ie=1,nel_off
	  lmax = ile(ie)
	  do l=1,lmax
	    utlnv(l,ie) = rt*ut(l,ie,1) + rr*ut(l,ie,2)
	    vtlnv(l,ie) = rt*vt(l,ie,1) + rr*vt(l,ie,2)
	  end do
	  do ii=1,3
	    zenv(ii,ie) = rt*ze(ii,ie,1) + rr*ze(ii,ie,2)
	  end do
	end do

	do k=1,nkn_off
	  lmax = ilk(k)
	  do l=0,lmax
	    wlnv(l,k) = rt*wn(l,k,1) + rr*wn(l,k,2)
	  end do
	  znv(k) = rt*zn(k,1) + rr*zn(k,2)
	end do

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine off_copy

	!use levels
	!use basin, only : nkn,nel,ngr,mbw
	use mod_offline

	implicit none

	integer ie,ii,k,l,lmax
	integer ito,ifrom

	do ito=1,nintp-1

	  ifrom = ito + 1

	  iitime(ito) = iitime(ifrom)

	  ut(:,:,ito) = ut(:,:,ifrom)
	  vt(:,:,ito) = vt(:,:,ifrom)
	  ze(:,:,ito) = ze(:,:,ifrom)

	  zn(:,ito) = zn(:,ifrom)
	  wn(:,:,ito) = wn(:,:,ifrom)
	  sn(:,:,ito) = sn(:,:,ifrom)
	  tn(:,:,ito) = tn(:,:,ifrom)
	  vd(:,:,ito) = vd(:,:,ifrom)
	  dd(:,:,ito) = dd(:,:,ifrom)
	  
	end do

	end

!****************************************************************
	
	subroutine off_reset

	!use levels
	!use basin, only : nkn,nel,ngr,mbw
	use mod_offline

	implicit none

	integer ie,ii,k,l,lmax

	dtr = 0.

	ut(:,:,1) = 0.
	vt(:,:,1) = 0.
	ze(:,:,1) = 0.
	zn(:,1) = 0.
	wn(:,:,1) = 0.
	sn(:,:,1) = 0.
	tn(:,:,1) = 0.
	vd(:,:,1) = 0.
	dd(:,:,1) = 0.

	end

!****************************************************************
	
	subroutine off_accum(dt)

! this subroutine should accumulate, but now it only takes a snapshot

	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use mod_diff_visc_fric
	!use levels
	!use basin, only : nkn,nel,ngr,mbw
	use mod_offline

	implicit none

	real dt

	integer ie,ii,k,l,lmax
	double precision dtt

	dtt = dt
	dtr = dtr + dtt
	dtr = 1.

	do ie=1,nel_off
	  lmax = ile(ie)
	  do l=1,lmax
	    !ut(l,ie,1) = ut(l,ie,1) + utlnv(l,ie) * dtt
	    !vt(l,ie,1) = vt(l,ie,1) + vtlnv(l,ie) * dtt
	    ut(l,ie,1) = utlnv(l,ie)
	    vt(l,ie,1) = vtlnv(l,ie)
	  end do
	  do ii=1,3
	    !ze(ii,ie,1) = ze(ii,ie,1) + zenv(ii,ie) * dtt
	    ze(ii,ie,1) = zenv(ii,ie)
	  end do
	end do

	do k=1,nkn_off
	  lmax = ilk(k)
	  do l=1,lmax
	    !wn(l,k,1) = wn(l,k,1) + wlnv(l,k) * dtt
	    wn(l,k,1) = wlnv(l,k)
	    !sn(l,k,1) = sn(l,k,1) + saltv(l,k) * dtt
	    sn(l,k,1) = saltv(l,k)
	    !tn(l,k,1) = tn(l,k,1) + tempv(l,k) * dtt
	    tn(l,k,1) = tempv(l,k)
	    vd(l,k,1) = visv(l,k)
	    dd(l,k,1) = difv(l,k)
	  end do
	  !zn(k,1) = zn(k,1) + znv(k) * dtt
	  zn(k,1) = znv(k)
	  !wn(0,k,1) = wn(0,k,1) + wlnv(0,k) * dtt
	  wn(0,k,1) = wlnv(0,k)
	  vd(0,k,1) = visv(0,k)
	  dd(0,k,1) = difv(0,k)
	end do

	end

!****************************************************************
	
	subroutine off_aver

	!use levels
	!use basin, only : nkn,nel,ngr,mbw
	use mod_offline

	implicit none

	integer ie,ii,k,l,lmax
	double precision rr

	rr = 0.
	if( dtr .gt. 0. ) rr = 1. / dtr

	do ie=1,nel_off
	  lmax = ile(ie)
	  do l=1,lmax
	    ut(l,ie,1) = ut(l,ie,1) * rr
	    vt(l,ie,1) = vt(l,ie,1) * rr
	  end do
	  do ii=1,3
	    ze(ii,ie,1) = ze(ii,ie,1) * rr
	  end do
	end do

	do k=1,nkn_off
	  lmax = ilk(k)
	  do l=1,lmax
	    wn(l,k,1) = wn(l,k,1) * rr
	    sn(l,k,1) = sn(l,k,1) * rr
	    tn(l,k,1) = tn(l,k,1) * rr
	    vd(l,k,1) = vd(l,k,1) * rr
	    dd(l,k,1) = dd(l,k,1) * rr
	  end do
	  zn(k,1) = zn(k,1) * rr
	  wn(0,k,1) = wn(0,k,1) * rr
	  vd(0,k,1) = vd(0,k,1) * rr
	  dd(0,k,1) = dd(0,k,1) * rr
	end do

	end

!****************************************************************
!****************************************************************
!****************************************************************

        subroutine off_print_debug_var

        use mod_offline

        implicit none

	logical, parameter :: bprintoffvar = .false.
	logical is_time_for_offline

	if( .not. bprintoffvar ) return

	if( is_time_for_offline() ) then
	  write(6,*) 'off_print_debug_var'
	  call print_debug_var
	end if

        end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine off_check(ig)

	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use mod_offline

	implicit none

	integer ig

	integer ie,ii,k,l,lmax
	integer ierr
	real utmax,umax,zmax,wmax,smax,tmax

	ierr = 0
	utmax = 10000.
	umax = 10.
	zmax = 10.
	wmax = 10.
	smax = 100.
	tmax = 100.

	do ie=1,nel_off
	  lmax = ile(ie)
	  do l=1,lmax
	    call off_check_val('ut',ie,l,real(ut(l,ie,ig)),utmax,ierr)
	    call off_check_val('vt',ie,l,real(vt(l,ie,ig)),utmax,ierr)
	    call off_check_val('utlnv',ie,l,utlnv(l,ie),utmax,ierr)
	    call off_check_val('vtlnv',ie,l,vtlnv(l,ie),utmax,ierr)
	    call off_check_val('utlov',ie,l,utlov(l,ie),utmax,ierr)
	    call off_check_val('vtlov',ie,l,vtlov(l,ie),utmax,ierr)
	    call off_check_val('ulnv',ie,l,ulnv(l,ie),umax,ierr)
	    call off_check_val('vlnv',ie,l,vlnv(l,ie),umax,ierr)
	    call off_check_val('ulov',ie,l,ulov(l,ie),umax,ierr)
	    call off_check_val('vlov',ie,l,vlov(l,ie),umax,ierr)
	  end do
	  do ii=1,3
	    call off_check_val('ze',ie,ii,real(ze(ii,ie,ig)),zmax,ierr)
	  end do
	end do

	do k=1,nkn_off
	  lmax = ilk(k)
	  do l=1,lmax-1
	    call off_check_val('wn',k,l,real(wn(l,k,ig)),wmax,ierr)
	  end do
	  do l=1,lmax
	    call off_check_val('sn',k,l,real(sn(l,k,ig)),smax,ierr)
	    call off_check_val('tn',k,l,real(tn(l,k,ig)),tmax,ierr)
	  end do
	  call off_check_val('zn',k,0,real(zn(k,ig)),zmax,ierr)
	end do

	if( ierr .gt. 0 ) then
	  write(6,*) 'errors checking variables read from file'
	  write(6,*) iitime(ig),ig,ierr
	  stop 'error stop off_check: out of range'
	else
	  !write(6,*) 'finished offline error check... ok... ',it
	end if

	end

!****************************************************************

	subroutine off_check_val(what,iek,l,val,vmax,ierr)

	implicit none

	character*(*) what
	integer iek,l
	real val,vmax
	integer ierr

	if( abs(val) .gt. vmax ) then
	  write(6,*) what,iek,l,val
	  ierr = ierr + 1
	end if

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine off_write(iu,it)

        use mod_offline
        use shympi

	implicit none

	integer iu,it

	if( shympi_is_master() ) then
          write(6,*) 'writing offline record for time ',it
	end if

	call off_write_record(iu,it)

	end 

!****************************************************************

        subroutine off_read(iu,ig,ierr)

        use mod_offline
        use shympi

        implicit none

        integer iu,ig
        integer ierr

	integer it

        call off_read_record(iu,ig,it,ierr)

	if( shympi_is_master() ) then
          write(6,*) 'reading offline record for time ',it
	end if

	end

!****************************************************************

