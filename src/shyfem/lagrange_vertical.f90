
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2019  Georg Umgiesser
!    Copyright (C) 2015,2017  Christian Ferrarin
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

! deal with vertical velocities
!
! revision log :
!
! 01.04.2015	ggu	changed VERS_7_1_7
! 23.04.2015	ggu	internal coodinates finished
! 06.05.2015	ccf	included settling velocity for lagrangian
! 08.05.2015	ggu&ccf	bug fix in track_xi_next_element
! 14.05.2015	ccf	bug fix in track_xi
! 21.05.2015	ggu	changed VERS_7_1_11
! 05.06.2015	ggu	changed VERS_7_1_12
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 16.11.2015	ggu	changed VERS_7_3_14
! 15.02.2016	ggu	handle particles on vertical wall gracefully
! 19.02.2016	ggu	changed VERS_7_5_3
! 22.02.2016	ggu	changed VERS_7_5_4
! 09.05.2017	ggu	changed VERS_7_5_26
! 26.05.2017	ccf	handle particles on open boundary
! 25.10.2018	ggu	changed VERS_7_5_51
! 16.02.2019	ggu	changed VERS_7_5_60
! 22.05.2023	ggu	get_layer_thickness() was missing an argument
! 05.06.2023    lrp     introduce z-star
! 28.09.2023    lrp     bug fix for zstar (wrong declaration nadapt/hadapt)

!
!******************************************************

	subroutine getzvel(ie,z0,l0,is_in,is_out,a_in,a_out,w)

! returns vertical velocity to be used in lagrangian model

	implicit none

	integer ie	!element number
	real z0		!depth of particle
	integer l0	!layer in which particle is in [1-nlv]
	integer is_in	!side from which particle comes [1-3]
	integer is_out	!side to which particle goes [1-3]
	real a_in	!relative postion of side is_in [0-1]
	real a_out	!relative postion of side is_out [0-1]
	real w		!computed vertical velocity for ie and l0 (return)

	integer ii
	real win,wout

!	-----------------------------------------------
!	velocity at entering and exiting point
!	-----------------------------------------------

	if( ie .le. 0 ) then
	  w = 0.
	  return
	end if

	ii = mod(is_in,3) + 1
	write(6,*) 'gguuut: ',ie,l0,ii,a_in
	call getzvel_point(ie,l0,ii,a_in,win)

	ii = mod(is_out,3) + 1
	call getzvel_point(ie,l0,ii,a_out,wout)

!	-----------------------------------------------
!	average velocity in element
!	-----------------------------------------------

	w = 0.5 * ( win + wout )

!	-----------------------------------------------
!	end of routine
!	-----------------------------------------------

	end

!******************************************************

	subroutine getzvel_point(ie,l0,is,a,w)

! returns vertical velocity at point given by ie,l0,is,a

	use mod_hydro_vel
	use basin

	implicit none

	integer ie	!element number
	integer l0	!layer in which particle is in
	integer is	!side where particle is [1-3]
	real a		!relative postion of side
	real w		!computed vertical velocity for ie and l0 (return)

	integer ii,k1,k2
	real wo1,wo2,wn1,wn2
	real w1,w2

!	-----------------------------------------------
!	velocity at first node
!	-----------------------------------------------

	ii = is

	k1 = nen3v(ii,ie)
	wo1 = wlov(l0,k1) + wlov(l0-1,k1)
	wn1 = wlnv(l0,k1) + wlnv(l0-1,k1)
	w1 = 0.25 * ( wo1 + wn1 )

!	-----------------------------------------------
!	velocity at second node
!	-----------------------------------------------

	ii = mod(ii,3) + 1

	k2 = nen3v(ii,ie)
	wo2 = wlov(l0,k2) + wlov(l0-1,k2)
	wn2 = wlnv(l0,k2) + wlnv(l0-1,k2)
	w2 = 0.25 * ( wo2 + wn2 )

!	-----------------------------------------------
!	interpolate velocity
!	-----------------------------------------------

	w = (1.-a)*w1 + a*w2

!	-----------------------------------------------
!	end of routine
!	-----------------------------------------------

	end

!******************************************************

	subroutine lagr_layer_thickness(ie,lmax,hl,htot,htotz)

! computes layer thickness for element ie

	use mod_depth
	use mod_hydro
	use levels

	implicit none

	integer ie		!element number
	integer lmax		!max number of layers (dimension in, actual out)
	real hl(lmax)		!layer thickness (return)
	real htot,htotz		!total depth without and with zeta (return)

	integer nlev,nsigma,ii,lmax_act
	real hsigma
	real z,h
	integer nadapt(4)
	real hadapt(4)

        !call compute_sigma_info(nlev,hlv,nsigma,hsigma)
	call get_sigma_info(nlev,nsigma,hsigma)
        call get_zadapt_info(ie,nadapt,hadapt) 
	lmax_act = ilhv(ie)
	if( lmax_act > lmax ) goto 99
	if( lmax_act > nlev ) goto 99
	lmax = lmax_act

	h = hev(ie)
	z = 0.
	do ii=1,3
	  z = z + zenv(ii,ie) + zeov(ii,ie)
	end do
	z = z / 6.

	call get_layer_thickness(lmax,1,nsigma,nadapt(4), &
     &				 hsigma,hadapt(4),z,h,hlv,hl)
	htot = h
	htotz = h + z

	return
   99	continue
	write(6,*) 'nlev,lmax,lmax_act: ',nlev,lmax,lmax_act
	stop 'error stop lagr_layer_thickness: incompatible layers'
	end

!*******************************************************************

        subroutine vertpos(zn0,deltat,layd,w,zn1,ztime,addl)

	implicit none

	real zn0		!initial relative vertical position in layer
	real deltat		!time spent in element(horizontal)
	real w			!average vertical velocity
	real zn1		!final relative vertical position in layer
	real ztime		!time spent in element (vertical)
        integer addl		!relative movement to next layer [-1,0,+1]
        real layd		!layer thickness

! return are zn1,ztime,addl

! lstd		max vertical distance to be traveled
! cpstd		computed distance to be traveled

	real cpdst,lstd

	lstd = 0
        if(w.gt.0) lstd=layd*(1-zn0)
        if(w.lt.0) lstd=-layd*zn0

        cpdst=deltat*w
        zn1=zn0+(cpdst/layd)
        
	if( w == 0 ) then
	  ztime = 2*deltat
	else
          ztime=lstd/w
	end if

        if(ztime.le.deltat)then
          if(w.gt.0)addl=1
          if(w.lt.0)addl=-1
        else
          addl=0
        end if
	
	if(zn1.gt.1)zn1=0
        if(zn1.lt.0)zn1=1

        end

!************************************************************

	subroutine getalfa(side,d,near,far,re)

	implicit none

        integer side		!number of side of element [1-3]
	real d			!relative distance from closest vertex
	integer near		!number of closest vertex [1-3]
	integer far		!to be deleted...
	real re			!relative distance to first point of side

	integer ind

        ind=mod(side,3)+1
	
	if(ind.eq.near)then
	  re=d
	else
	  re=1-d
	end if
	
	end 

!************************************************************
!************************************************************
!************************************************************
!
! z = 0		top of layer
! z = 1		bottom of layer
!
!************************************************************

	subroutine track_xi(id,iel,lb,sv,xi,z,time)

	use mod_lagrange

	implicit none

	integer id			!id of particle
	integer iel			!element (in and out)
	integer lb			!level (in and out)
	double precision sv		!sinking velocity
	double precision xi(3)		!natural coordinates in triangle
	double precision z		!rel vert pos: 0=top, 1=bottom
	real time

	logical bdebug
	logical bsurf,bbott
	integer iflux,lmax,ii
	integer ieorig,lborig
	double precision s,alpha
	double precision dist,dh,dv,ds
	double precision t,th,tv,tt
	double precision vel,w
	double precision dz,hd
	double precision xis(3)
	double precision xie(3)

	logical track_xi_on_material_boundary
	logical track_xi_has_layer

	blgrdebug = id == 172
	blgrdebug = id == 1
	blgrdebug = id == 0
	!blgrdebug = id == 9939
	!blgrdebug = .true.
	bdebug = .false.
	bdebug = .true.
	bdebug = blgrdebug
	!bdebug = id == 9939

	ieorig = iel
	lborig = lb

	!write(6,*) 'tracking particle: ',id,time
	call track_xi_check('start track_xi',id,xi)

	!-----------------------------------------------
	! get advection information
	!-----------------------------------------------

	call track_xi_get_vertical(id,iel,lb,lmax,hd,w)
	w = w - sv				!total vertical velocity

	!-----------------------------------------------
	! handle particles on surface or on bottom
	!-----------------------------------------------

	bsurf = lb .eq. 1 .and. z .eq. 0.d0	!particle on surface
	bbott = lb .eq. lmax .and. z .eq. 1.d0	!particle on bottom

	if( bsurf .and. w > 0.d0 ) w = 0.d0
	if( bbott .and. w < 0.d0 ) w = 0.d0
	if( blgrsurf .or. blgr2d ) w = 0.d0	!advection only in surface layer

	if ( w > 0.d0 .and. z == 0.d0 ) then
	  lb = lb - 1
	  z  = 1.d0
	end if
	if ( w < 0.d0 .and. z == 1.d0 ) then
	  lb = lb + 1
	  z  = 0.d0
	end if

	!-----------------------------------------------
	! gets flux and vel information for element and layer
	!-----------------------------------------------

	call track_xi_get_flux(iel,lb,iflux,alpha,vel)

	if( vel < 0.d0 ) then
	  write(6,*) 'vel < 0'
	  write(6,*) vel,iel,lb,lmax
	  write(6,*) (flux3d(lb,ii,iel),ii=1,3)
	end if

	!-----------------------------------------------
	! determine time to leave element (th)
	!-----------------------------------------------

	th = 2.*time

	if( bdebug ) then
	  write(6,*) '======================================='
	  write(6,*) 'track_xi start debugging: ',id
	  write(6,*) '======================================='
	  write(6,*) iel,lb,lmax,iflux
	  write(6,*) alpha,vel,z
	  write(6,*) xi
	  !call track_xi_info_element(iel)
	end if

	if( vel > 0. ) then
	  if( iflux > 0 ) then		!out of flux node
	    call xit_start_end(iflux,alpha,xi,xis,xie,s)
	  else if( iflux < 0 ) then	!into flux node
	    call xit_start_end(-iflux,alpha,xi,xie,xis,s)
	    s = 1.d0 - s
	  else
	    time = 0.
	    return
	    !goto 99
	  end if
	  call xi_dist(iel,xis,xie,dist)
	  dh = dist*(1.d0 - s)		!distance to travel in element
	  th = dh / vel			!time to arrive at edge of element
	else
	  s = 1.
	  dist = 0.
	  dh = 0.
	  xis = xi
	  xie = xi
	end if

	!-----------------------------------------------
	! determine time to leave layer (tv)
	!-----------------------------------------------

	tv = 2.*time

	if( w /= 0. ) then
	  if( w > 0. ) then
	    dv = hd*z
	    tv = dv / w
	  else
	    dv = hd*(1.-z)
	    tv = dv / (-w)
	  end if
	else
	  dv = 0.
	end if

	!-----------------------------------------------
	! what happens first? (tt is total time available)
	!-----------------------------------------------

	tt = time
	t = tt
	if( th > 0. ) t = min(t,th)
	if( tv > 0. ) t = min(t,tv)
	!t = min(th,tv,tt)

	if( bdebug ) then
	  write(6,*) 'times for advection'
	  write(6,*) th,tv,tt
	  write(6,*) dist,s,dh
	  write(6,*) xis
	  write(6,*) xie
	end if

	!-----------------------------------------------
	! if landing on material boundary - artificially slow down particle
	!CCF IN THIS CASE WE DO NOT ALLOW BEACHING BY ADVECTION
	!-----------------------------------------------

	if( track_xi_on_material_boundary(iel,xie) ) then
	  th = 2.*min(tv,tt)	!longer than other times
	  t = min(th,tv,tt)
	  if( bdebug ) then
	    write(6,*) 'material boundary: ',id,iel,lb
	    write(6,*) th,tv,tt
	  end if
	end if

	!-----------------------------------------------
	! handle horizontal advection
	!-----------------------------------------------

	if( th == 0. ) then		!body not moving or laying on side
          if (vel > 0.) then
            xi = xie
            call track_xi_next_element(iel,xi)    !for particles laying on side
	    if( .not. track_xi_has_layer(iel,lb) ) then	!just do vert adv
	      iel = ieorig
	      xi = xie
            end if
          end if
	else if( th > t ) then		!body remains in element
	  ds = (1.-s)*t/th
	  s = s + ds
	  if( s > 1. ) goto 97
	  xi = (1.-s)*xis + s*xie
	  call track_xi_check('after advection track_xi 1',id,xi)
	else				!body reaches border of element
	  xi = xie
	  call track_xi_check('before advection track_xi 2',id,xi)
	  call track_xi_next_element(iel,xi)
	  call track_xi_check('after advection track_xi 2',id,xi)
	  if( .not. track_xi_has_layer(iel,lb) ) then	!just do vert adv
	    iel = ieorig
	    xi = xie
	    if( th == 0. ) then
	      t = min(tv,tt)	!just advect vertically
	      if( bdebug ) then
	        write(6,*) 'waiting for vertical: ',id
	        write(6,*) 'waiting... ',t,tt
	        write(6,*) 'waiting... ',iel,w
	        write(6,*) 'waiting... ',lb,z
	      end if
	    end if
	  end if
	end if

	call track_xi_check('after advection track_xi',id,xi)

	if( bdebug ) then
	  write(6,*) 'after horizontal advection'
	  write(6,*) ieorig,iel
	  write(6,*) lborig,lb,lmax
	  write(6,*) th,t,iel
	  write(6,*) s
	  write(6,*) xi
	end if

	!-----------------------------------------------
	! handle vertical advection
	!-----------------------------------------------

	if( tv == 0. ) then		!body not moving
	  !
	else if( tv > t ) then		!body remains in layer
	  dz = w*t/hd
	  z = z - dz
	else
	  if( w > 0. ) then
	    if( lb > 1 ) then
	      lb = lb - 1
	      z = 1.
	    else
	      z = 0.
	    end if
	  else
	    if( lb < lmax ) then
	      lb = lb + 1
	      z = 0.
	    else
	      z = 1.
	    end if
	  end if
	end if

	if( bdebug ) then
	  write(6,*) 'after vertical advection'
	  write(6,*) ieorig,iel
	  write(6,*) lborig,lb,lmax
	  write(6,*) w,z
	end if

	!-----------------------------------------------
	! compute remaining time and check for error
	!-----------------------------------------------

	if( iel .ne. ieorig .and. iel > 0 ) then
	  call track_xi_adjust_layer(id,ieorig,iel,lb,z)
	end if

	time = tt - t

	if( bdebug ) then
	  write(6,*) 'final time: ',time
	  write(6,*) iel,lb,z
	  write(6,*) 'track_xi end debugging'
	end if

	if( z < 0. .or. z > 1. ) goto 98

	!-----------------------------------------------
	! end of routine
	!-----------------------------------------------

	return 
   97	continue
	write(6,*) 's,t,th: ',s,t,th
	write(6,*) 'ds,sorig: ',ds,s-ds
	stop 'error stop track_xi: internal error (3)'
   98	continue
	write(6,*) 'id,iel ',id,iel
	write(6,*) 'tv,t ',tv,t
	write(6,*) 'dz,w,hd: ',dz,w,hd
	write(6,*) 'z,lb,lmax ',z,lb,lmax
	stop 'error stop track_xi: internal error (2)'
   99	continue
	write(6,*) 'id,iel ',id,iel
	write(6,*) 'vel,iflux ',vel,iflux
	write(6,*) 'dz,w,hd: ',dz,w,hd
	write(6,*) 'z,lb,lmax ',z,lb,lmax
	stop 'error stop track_xi: internal error (1)'
	end

!************************************************************

	subroutine track_xi_next_element(ie,xi)

! copies internal coordinates to new element - avoid falling on vertex

	use mod_lagrange
	use mod_geom

	implicit none

	integer ie
	double precision xi(3)

	logical bdebug
	integer ii,in,it
	integer ia1,ia2,ia3
	integer ib1,ib2,ib3
	integer ieb
	double precision r,eps
	double precision xiaux(3)

	eps = 1.e-5
	bdebug = blgrdebug

	!---------------------------------------
	! check xi - how many zeros
	!---------------------------------------

	in = 0
	it = 0
	do ii=1,3
	  if( xi(ii) == 0. ) then
	    in = in + 1
	    it = it + ii
	  end if
	end do

	if( bdebug ) then
	  write(6,*) 'track_xi_next_element: ',in,it
	end if

	!---------------------------------------
	! assign pointers - ia1 indicates side of particle
	!---------------------------------------

	if( in == 1 ) then	!normal mode - particle on side
	  ia1 = it
	  ia2 = mod(ia1,3) + 1
	  ia3 = mod(ia2,3) + 1
	else if( in == 2 ) then	!particle on vertex - must move
	  if( bdebug ) write(6,*) 'in==2 ',xi
	  r = 0.
	  do while( r == 0. )
	    call random_number(r)
	    r = eps*(r-0.5)
	  end do
	  it = 6 - it
	  if( r < 0. ) then
	    ia1 = mod(it+1,3) + 1
	    ia2 = mod(ia1,3) + 1
	    ia3 = mod(ia2,3) + 1
	    xi(ia2) = 1. + r
	    xi(ia3) = -r
	  else
	    ia1 = mod(it,3) + 1
	    ia2 = mod(ia1,3) + 1
	    ia3 = mod(ia2,3) + 1
	    xi(ia2) = r
	    xi(ia3) = 1. - r
	  end if
	else
	  write(6,*) 'xi not on side or impossible'
	  write(6,*) ie,in,it
	  write(6,*) xi
	  stop 'error stop track_xi_next_element: erroneous xi'
	end if

	!---------------------------------------
	! sanity checks (may be removed)
	!---------------------------------------

	if( xi(ia1) /= 0. ) goto 99
	if( abs(xi(ia2)+xi(ia3)-1.) > eps ) goto 99

	!---------------------------------------
	! look for neighboring element
	!---------------------------------------

	ieb = ieltv(ia1,ie)
	if( ieb < 1 ) then	!no element or open boundary
	  ie = -ie
	  return
	end if

	if( bdebug ) then
	  write(6,*) ia1,ia2,ia3,ieb
	end if

	!---------------------------------------
	! look for side of particle in new element
	!---------------------------------------

	do ii=1,3
	  if( ieltv(ii,ieb) == ie ) exit
	end do
	if( ii > 3 ) goto 99

	!---------------------------------------
	! copy xi to new elements
	!---------------------------------------

	ib1 = ii
	ib2 = mod(ib1,3) + 1
	ib3 = mod(ib2,3) + 1

	xiaux = xi

	xi(ib1) = 0.
	xi(ib2) = xiaux(ia3)
	xi(ib3) = xiaux(ia2)

	ie = ieb

	!---------------------------------------
	! end of routine
	!---------------------------------------

	return
   99	continue
	write(6,*) ie,in,it,ia1,ieb
	write(6,*) xi
	write(6,*) (ieltv(ii,ieb),ii=1,3)
	stop 'error stop track_xi_next_element: erroneous ieltv'
	end

!************************************************************

	function track_xi_on_material_boundary(ie,xi)

! checks if particle is on material boundary

	use mod_geom

	implicit none

	logical track_xi_on_material_boundary
	integer ie
	double precision xi(3)

	integer ii

	track_xi_on_material_boundary = .false.

	do ii=1,3
	  if( xi(ii) == 0. .and. ieltv(ii,ie) == 0 ) then
	    track_xi_on_material_boundary = .true.
	  end if
	end do

	end

!************************************************************

	function track_xi_has_layer(iel,lb)

	use levels

	implicit none

	logical track_xi_has_layer
	integer iel,lb

	track_xi_has_layer = ( ilhv(iel) >= lb )

	end

!************************************************************

	subroutine track_xi_adjust_layer(id,ieorig,iel,lb,z)

! adjusts layer when passing from one element to the next
!
! is only temporary - must also adjust w

	use levels

	implicit none

	integer id
	integer ieorig
	integer iel
	integer lb
	double precision z

	integer lb2
	double precision z2

	integer lmax

	lmax = ilhv(iel)

	!call track_new_vert_pos(id,ieorig,iel,lb,lb2,z,z2)
	!lb = lb2
	!z = z2

	if( lb > lmax ) then
	  lb = lmax
	  z = 1.
	  write(6,*) 'warning: wrong layer: ',iel,lb,lmax
	end if

	end

!************************************************************
!************************************************************
!************************************************************

	subroutine track_new_vert_pos(id,ie1,ie2,lb1,lb2,z1,z2)

	use levels

	implicit none

	integer id
	integer ie1,ie2
	integer lb1,lb2
	double precision z1,z2

	logical bdebug
	integer ierror
	integer nlvaux,nsigma
	integer lmax1,lmax2,l
	real hsigma
	real hl1(nlv)
	real hl2(nlv)
	real hp,hp1,hp2,hd
	real htot1,htot2,htotz1,htotz2

	real, parameter :: eps = 1.e-5

	call get_sigma_info(nlvaux,nsigma,hsigma)

	if( hsigma < 10000 ) then
	  stop 'error stop track_layer_thickness: hsigma'
	end if
	if( nsigma > 0 ) then
	  lb2 = lb1
	  z2 = z1
	  return
	end if

	bdebug = id == 341

	lmax1 = nlv
	call lagr_layer_thickness(ie1,lmax1,hl1,htot1,htotz1)
	if( lb1 > lmax1 ) goto 99
	if( z1 < 0. .or. z1 > 1. ) goto 99

	hd = 0.
	do l=1,lb1-1
	  hd = hd + hl1(l)
	end do
	hp1 = hd + hl1(lb1)*z1

	lmax2 = nlv
	call lagr_layer_thickness(ie2,lmax2,hl2,htot2,htotz2)

	hp2 = hp1 * (htotz2/htotz1)

	if( hp2 < 0 .or. hp2 > htot2 ) goto 99

	hd = 0.
	do l=1,lmax2
	  if( hd + hl2(l) >= hp2 ) exit
	  hd = hd + hl2(l)
	end do

	lb2 = l
	if( lb2 > lmax2 ) then	!HACK
	  z2 = -1.
	else
	  z2 = (hp2-hd)/hl2(l)
	end if

	if( lb2 > lmax2 .and. abs(z1-1.) < eps ) then
	  lb2 = lmax2
	  z2 = 0.
	else if( z2 > 1. .and. abs(z1-1.) < eps ) then
	  z2 = 1.
	end if
	if( z2 > 1. .and. abs(z2-1.) < eps ) z2 = 1.

	if( z2 < 0. .or. z2 > 1. ) goto 99
	if( lb2 > lmax2 ) goto 99
	!if( bdebug ) goto 99

	return
   99	continue
	ierror = 0
	write(6,*) 'debug track_new_vert_pos: ',id
	write(6,*) ie1,lb1,lmax1,z1
	write(6,*) ie2,lb2,lmax2,z2
	write(6,*) (hl1(l),l=1,lmax1)
	write(6,*) (hl2(l),l=1,lmax2)
	write(6,*) hd
	write(6,*) hp1,htot1,htotz1
	write(6,*) hp2,htot2,htotz2
	if( lb1 > lmax1 .or. lb2 > lmax2 ) ierror = 1
	if( hp1 < 0. .or. hp1 > htotz1 ) ierror = 2
	if( hp2 < 0. .or. hp2 > htotz2 ) ierror = 2
	if( z1 < 0. .or. z1 > 1. ) ierror = 3
	if( z2 < 0. .or. z2 > 1. ) ierror = 3
	write(6,*) 'internal error: ',ierror
	stop 'error stop track_new_vert_pos: internal error'
	end

!************************************************************
!************************************************************
!************************************************************

	subroutine track_xi_get_flux(iel,lb,iflux,alpha,vel)

! gets flux and vel information for element and layer

	use mod_lagrange
	use mod_geom
	use mod_hydro_vel

	implicit none

	integer iel			!element number
	integer lb			!layer
	integer iflux			!node of flux line
	double precision alpha		!fraction of opsosite side
	double precision vel		!velocity in element (always positive)

	logical bdebug
	integer ii,in,io,nn,no,inext,imax
	real az,azt,u,v
	double precision flux(3)
	double precision fp,fm,fmmax,fact

	bdebug = .false.
	bdebug = blgrdebug

	az = azlgr
	azt = 1. - az

	flux = flux3d(lb,:,iel)		!fluxes over side into element

	do ii=1,3
	  if( ieltv(ii,iel) .eq. 0 ) then	!closed boundary - set to 0
	    if( flux(ii) .ne. 0. ) then
	      write(6,*) 'material boundary: ',iel,ii,flux(ii)
	      flux(ii) = 0.
	    end if
	  end if
	end do

	u = azt*ulov(lb,iel) + az*ulnv(lb,iel)
	v = azt*vlov(lb,iel) + az*vlnv(lb,iel)
	vel = sqrt(u*u+v*v)

	if( bdebug ) then
	  write(6,*) 'track_xi_get_flux: ',iel,lb
	  write(6,*) flux
	  write(6,*) vel
	end if

	in = 0
	io = 0
	nn = 0
	no = 0
	fp = 0.
	fm = 0.
	do ii=1,3
	  if( flux(ii) > 0. ) then	!flux into element
	    in = in + 1
	    nn = nn + ii
	    fp = fp + flux(ii)
	  else if( flux(ii) < 0. ) then
	    io = io + 1
	    no = no + ii
	    fm = fm - flux(ii)
	  end if
	end do

! fact is introduced to compense for flux if fp /= fm
! always scale to outgoing flux fm

	if( in == 1 .and. io > 0 ) then
	  iflux = -nn			!into this node
	  inext = mod(nn,3) + 1
	  fact = fm/fp
	  alpha = 1 + flux(inext)/fm
	else if( in == 2 .and. io > 0 ) then
	  iflux = 6 - nn		!out from this node
	  inext = mod(iflux,3) + 1
	  fact = fm/fp
	  alpha = 1 - flux(inext)/fp
	else if( io == 0 ) then		!only ingoing - do not advect
	  iflux = 0
	  alpha = 0.
	  vel = 0.
	  fact = 0.
	else if( in == 0 .and. io > 0 ) then	!only outgoing - special
	  imax = 0
	  fmmax = 0.
	  do ii=1,3
	    if( flux(ii) < fmmax ) then
	      fmmax = flux(ii)
	      imax = ii
	    end if
	  end do
	  iflux = imax
	  inext = mod(iflux,3) + 1
	  alpha = 1 - flux(inext)/fm
	  fact = 1.
	else				!just to be sure...
	  iflux = 0
	  alpha = 0.
	  vel = 0.
	  fact = 0.
	end if

	vel = vel * fact

	if( bdebug ) then
	  write(6,*) in,nn
	  write(6,*) iflux,inext
	  write(6,*) fp,fm,alpha
	  write(6,*) 'track_xi_get_flux end'
	end if

	if( bback ) iflux = -iflux

	alpha = dmin1(alpha,1.d0)
	alpha = dmax1(alpha,0.d0)

	end

!************************************************************

	subroutine track_xi_get_vertical(id,iel,lb,lmax,hd,w)

	use mod_hydro_vel
	use levels, only : nlvdi,nlv
	use basin

	implicit none

	integer id			!id of particle
	integer iel			!element number
	integer lb			!layer
	integer lmax			!maximum layers in element (return)
	double precision hd		!layer thickness (return)
	double precision w		!vertical velocity (return)

	integer ii,k
	real wo,wn
	real hl(nlv)
	real htot,htotz

	lmax = nlv
	call lagr_layer_thickness(iel,lmax,hl,htot,htotz)

	if( lb > lmax ) goto 99

	hd = hl(lb)

	w = 0.
	do ii=1,3
	  k = nen3v(ii,iel)
	  wo = wlov(lb,k) + wlov(lb-1,k)
	  wn = wlnv(lb,k) + wlnv(lb-1,k)
	  w = w + wo + wn
	end do

	w = w / 12.

	return
   99	continue
	write(6,*) 'id,iel,l,lmax: ',id,iel,lb,lmax
	stop 'error stop track_xi_get_vertical: no such layer'
	end

!************************************************************

	subroutine track_xi_check(text,id,xi)

	implicit none

	character*(*) text
	integer id
	double precision xi(3)

	logical berror
	integer ii
	double precision tot,eps

	berror = .false.
	eps = 1.d-5

	tot = 0.d0
	do ii=1,3
	  tot = tot + xi(ii)
	  if( xi(ii) > 1.d0 ) berror = .true.
	  if( xi(ii) < 0.d0 ) berror = .true.
	end do
	if( abs(tot-1.d0) > eps ) berror = .true.

	if( berror ) then
	  write(6,*) text
	  write(6,*) id
	  write(6,*) xi
	  stop 'error stop track_xi_check: wrong xi'
	end if

	end

!************************************************************

	subroutine track_xi_info_element(ie)

! prints information on element

	use mod_geom
	use mod_hydro_vel
	use levels, only : nlvdi,nlv
	use basin

	implicit none

	integer ie

	integer ii,k,lmax,l
	real hl(nlv)
	real htot,htotz

	write(6,*) 'info on element ',ie

	do ii=1,3
	  k = nen3v(ii,ie)
	  write(6,*) ii,k,xgv(k),ygv(k)
	end do

	write(6,*) (ieltv(ii,ie),ii=1,3)

	lmax = nlv
	call lagr_layer_thickness(ie,lmax,hl,htot,htotz)
	write(6,*) lmax,htot,htotz
	write(6,*) (hl(l),l=1,lmax)

	write(6,*) (ulov(l,ie),l=1,lmax)
	write(6,*) (vlov(l,ie),l=1,lmax)
	write(6,*) (ulnv(l,ie),l=1,lmax)
	write(6,*) (vlnv(l,ie),l=1,lmax)

	write(6,*) 'info on element end'

	end

!************************************************************

