c
c $Id: subn11.f,v 1.69 2010-03-22 15:29:31 georg Exp $
c
c boundary condition routines
c
c contents :
c
c subroutine sp111(mode)			set up boundary/initial cond.
c subroutine tilt				tilting of boundary surface
c
c subroutine initilt(ibc)       		finds tilting node in node list
c subroutine iniflux(kranf,krend)		initializes flux boundary
c
c subroutine mvalue(ibc,it,rmu,rmv)		computes momentum value
c subroutine bvalue(ibc,it,rwv)			computes z value for boundary
c function cvalue(ibc,it,what)			computes conz value for bound.
c subroutine rain( t , dt )			simulates rain increasing z
c
c revision log :
c
c revised 05.08.92 by ggu   $$ibtyp3 - implementation of ibtyp=3
c revised 31.08.92 by ggu   $$rqv1   - rqv initialized in sp159b, not here
c revised 05.09.92 by ggu   $$close1 - rqv initialized here
c revised 27.10.93 by ggu   $$roger  - ibtyp=70 (nwe-shelf)
c revised 05.11.93 by ggu   $$roger  - call to roger has been commented
c revised 11.01.94 by ggu   $$restart - restart is reading all variables
c revised 20.01.94 by ggu   $$zeov - initialize zeov
c revised 20.01.94 by ggu   $$conz - impl. of conz with bnd(12,.)
c revised 31.01.94 by ggu   $$conzbc - impl. of bc for conz with rcv
c revised 24.03.94 by ggu   $$surel - impl. of distributed source/sink
c revised 02.04.96 by ggu   $$exxqq - interpolation for more files ok
c revised 25.06.97 by ggu   complete restructured -> new subroutines
c revised 18.09.97 by ggu   $$FLUX3 - special treatment of type 3 boundary
c revised 23.09.97 by ggu   concentration boundary as file implemented
c revised 03.12.97 by ggu   $$TS - temp/salt implemented
c 20.05.1998	ggu	rain from file implemented
c 22.05.1998	ggu	local variable t introduced
c 22.05.1998	ggu	corrected bug for rain (surel called twice)
c 28.05.1998    ggu     new ruv, rvv (momentum input) (incomplete)
c 20.06.1998    ggu     more on momentum input (mvalue)
c 29.06.1998    ggu     !$$momin1 - bug fix for momentum input
c 13.07.1998    ggu     !initialize ruv, rvv by node
c 14.07.1998    ggu     finally momentum input finished -> new exxqq
c 20.08.1998    ggu     iextpo finally eliminated
c 20.08.1998    ggu     2d dependent routines copied to subini.f
c 21.08.1998    ggu     xv eliminated
c 24.08.1998	ggu	BC for concentration is bnd(20,..)
c 24.08.1998	ggu	BC for maximum input level is bnd(12,..) -> levmax
c 05.11.1998	ggu	slightly restructured restart
c 18.02.1999	ggu	allow for boundary type 0
c 20.04.1999	ggu	converted to stress instead of wind (tauxnv...)
c 01.12.1999	ggu	handle negative boundary type
c 07.05.2001	ggu	introduced variable zfact
c 24.10.2001	ggu	ignore boundary with type 0
c 10.08.2003	ggu	better commented, new routines meteo_init, meteo_force
c 14.08.2003	ggu	initialization of z and uv made explicit
c 10.03.2004	ggu	RQVDT - value in rqv is now discharge [m**3/s]
c 03.09.2004	ggu	restart taken out to ht
c 11.10.2004	ggu	new ibtyp=31, multiply with zfact also for sin (ZFACT)
c 11.03.2005	ggu	new boundary routines b3dvalue, c3dvalue (3D bounds)
c 17.02.2006	ggu	new routine get_bflux()
c 23.03.2006    ggu     changed time step to real
c 18.09.2007    ggu     new set_mass_flux, bug fix in dist_3d
c 25.09.2007    ggu     routines deleted: [mbc]value, testbc
c 02.10.2007    ggu     bug fix in make_scal_flux: surface flux only in layer 1
c 08.10.2007    ggu     deleted commented lines
c 17.03.2008    ggu     name of some variables changed, new scalar values
c 07.04.2008    ggu     deleted c2dvalue, set_scal_bc
c 07.04.2008    ggu     differential input introduced (-5555)
c 09.04.2008    ggu     only level boundary array, re-arranged
c 10.04.2008    ggu&ccf new call to init_z0b()
c 17.04.2008    ggu     lmin introduced in set_mass_flux (negative levmax)
c 17.04.2008    ggu     evaporation introduced, rain adjusted
c 18.04.2008    ggu     rain simplified, bugfix
c 22.04.2008    ggu     in make_scal_flux do not alter mfluxv (parallel code)
c 03.06.2008    ggu     levmin introduced
c 24.03.2009    ggu     bug fix for rain; rain0d; new rain2distributed()
c 27.03.2009    ggu     new routine adjust_mass_flux() for dry nodes
c 31.03.2009    ggu     in make_scal_flux() bug fix for negative fluxes
c 02.04.2009    ggu     use get_bnd_(i)par() for special routines
c 03.04.2009    ggu     set intpol depending on ibtyp if intpol==0
c 20.04.2009    ggu     new routine z_tilt (tilting around a fixed value)
c 27.04.2009    ggu     can use ktilt with ztilt
c 19.01.2010    ggu     in make_scal_flux() return also sconz at bound
c 26.02.2010    ggu     bug fix in make_scal_flux() - sconz was out of l-loop
c 01.03.2010    deb     tramp introduced
c 10.03.2010    ggu     new routine check_scal_flux()
c 22.03.2010    ggu     in make_scal_flux() forgotten to initialize sconz
c
c***************************************************************

	subroutine sp111(mode)

c set up boundary and initial conditions
c
c mode		1 : first call, initialize b.c.
c		2 : read in b.c.

	implicit none

        include 'param.h'

	integer mode

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it

	real eps1,eps2,pi,flag,high,higi
	common /mkonst/ eps1,eps2,pi,flag,high,higi

	real rrv(1)
	common /rrv/rrv
	integer irv(1)
	common /irv/irv
	real rzv(1), rqv(1)
	common /rzv/rzv, /rqv/rqv
	real rqpsv(1), rqdsv(1)
	common /rqpsv/rqpsv, /rqdsv/rqdsv
	real mfluxv(nlvdim,1)
	common /mfluxv/mfluxv

	integer nlvdi,nlv
	common /level/ nlvdi,nlv

        real ruv(1)
        real rvv(1)
        common /ruv/ruv
        common /rvv/rvv

        real bnd3(nb3dim,0:nbcdim)     !boundary array for water level
        real rwv(nb3dim)
        save bnd3

	logical bimpose
	integer kranf,krend,k,kn
	integer ibc,ibtyp
        integer nk,i,kk,kindex,iv
        integer nbdim,nsize
        integer iunrad,ktilt
	integer ip,l
	real rw,const,aux
	real dt
	real conz,temp,salt
	real conzdf,tempdf,saltdf
c	real dz
	real rmu,rmv
	real getpar,rwint
	real conz3,temp3,salt3
	real tramp,alpha
	integer iround
        integer nkbnds,kbnds,itybnd

	integer ipext,kbndind

	tramp = 0.
	!tramp = 86400. !DEB

	call get_timestep(dt)

	if( mode .eq. 1 ) goto 1
	if( mode .eq. 2 ) goto 2
	stop 'error stop: internal error sp111 (1)'

c---------------------------------------------------------------
c first call %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c---------------------------------------------------------------

    1	continue

c	-----------------------------------------------------
c       initialize meteo
c	-----------------------------------------------------

	call meteo_init

c	-----------------------------------------------------
c       initialize tilted and flux boundaries
c	-----------------------------------------------------

	do ibc=1,nbc
	  ibtyp=itybnd(ibc)
          nk = nkbnds(ibc)

          bimpose = ibtyp .ge. 1 .and. ibtyp .le. 2

          if( bimpose .and. nk .le. 1 ) goto 95 !$$ibtyp3

	  if( bimpose ) then    		!$$FLUX3 - not for ibtyp=3
	    call initilt(ibc)           	!find nodes for tilting
	    call iniflux(ibc)   		!set up rlv,rhv,rrv,ierv
	  end if
	end do

c	-----------------------------------------------------
c       initialization of aux array
c	-----------------------------------------------------

c maybe use directly bnds_init ?? FIXME

	do ibc=0,nbc
	  do i=1,nb3dim
	    bnd3(i,ibc) = 0.
	  end do
	end do

c	-----------------------------------------------------
c       determine constant for z initialization
c	-----------------------------------------------------

	const=getpar('const')	!constant initial z value

	do ibc=1,nbc
	  ibtyp=itybnd(ibc)

	  if(const.eq.flag.and.ibtyp.eq.1) then
	        call get_bnd_ipar(ibc,'nbdim',nbdim)
	        nk = nkbnds(ibc)
		nsize = 0
		if( nbdim .gt. 0 ) nsize = nk
	        call b3dvalue(ibc,itanf,nsize,nb3dim,bnd3(1,ibc),rwv)
		const=rwv(1)	!prepare constant z value for start
	  end if

	  if(ibtyp.eq.70) then	!nwe-shelf	!$$roger - special b.c.
c	    call roger(rzv,dt,0)
c	    call roger(rzv,irv,nrb,dt,0)
	  end if
	end do

	if(const.eq.flag) const=0.
	call putpar('const',const)

c	-----------------------------------------------------
c       initialize variables or restart
c	...the variables that have to be set are zenv, utlnv, vtlnv
c	-----------------------------------------------------

	call init_z(const)		!initializes zenv
	call init_uvt			!initializes utlnv, vtlnv
	call init_z0b			!initializes bottom roughness

c	-----------------------------------------------------
c       finish
c	-----------------------------------------------------

c next only for radiation condition
c
c        write(78,*) 46728645,1
c        write(78,*) 50,0
c        write(78,*) 0,50,-0.01,0.01
c        write(79,*) 46728645,1
c        write(79,*) 50,0
c        write(79,*) 0,50,-0.01,0.01

	return

c---------------------------------------------------------------
c normal call for boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%
c---------------------------------------------------------------

    2	continue

c	-----------------------------------------------------
c	initialize node vectors with boundary conditions
c	-----------------------------------------------------

	do k=1,nkn
          rzv(k)=flag
          rqv(k)=0.	!$$rqv1 !$$close1	[m**3/s]
          rqpsv(k)=0.	!fluxes - point sources [m**3/s]
          rqdsv(k)=0.	!fluxes - distributed sources through surface [m**3/s]
          ruv(k)=0.	!momentum input
          rvv(k)=0.
	end do

c	-----------------------------------------------------
c	loop over boundaries
c	-----------------------------------------------------

        call bndo_radiat(it,rzv)

	do ibc=1,nbc

          call get_bnd_ipar(ibc,'ibtyp',ibtyp)
          call get_bnd_ipar(ibc,'nbdim',nbdim)

          nk = nkbnds(ibc)   !total number of nodes of this boundary

	  if( ibtyp .le. 3 ) then
	   if( nbdim .eq. 0 ) then
	      nsize = 0
	      call b3dvalue(ibc,it,nsize,nb3dim,bnd3(1,ibc),rwv)
              rw = rwv(1)
	      call dist_horizontal(1,rwv,nk,rw)
	    else
	      nsize = nk
	      call b3dvalue(ibc,it,nsize,nb3dim,bnd3(1,ibc),rwv)
              call aver_horizontal(1,rwv,nk,rw)
	    end if
	    call set_bnd_par(ibc,'zval',rw)
	  else if( ibtyp .eq. 31 ) then
	    nsize = 0
	    call b3dvalue(ibc,it,nsize,nb3dim,bnd3(1,ibc),rwv)
	    call dist_horizontal(1,rwv,nk,rwv(1))
	  else if( ibtyp .eq. 4 ) then
            stop 'error stop: momentum input is broken...'
	  end if

	  alpha = 1.
	  if( tramp .gt. 0. .and. it-itanf .le. tramp ) then
	     alpha = (it-itanf) / tramp
	  end if

	  do i=1,nk

             kn = kbnds(ibc,i)
	     rw = rwv(i)

	     if(ibtyp.eq.1) then		!z boundary
               rzv(kn)=rw
	     else if(ibtyp.eq.2) then		!q boundary
	       kindex = kbndind(ibc,i)
               rqpsv(kn)=alpha*rw*rrv(kindex)	!BUGFIX 21-08-2002, RQVDT
             else if(ibtyp.eq.3) then		!$$ibtyp3 - volume inlet
               rqpsv(kn) = alpha*rw
             else if(ibtyp.eq.4) then		!momentum input
	       ruv(kn) = rmu
	       rvv(kn) = rmv
             else if(ibtyp.eq.31) then		!zero gradient for z
               ! already done...
	       call get_bnd_ipar(ibc,'ktilt',ktilt)
               rzv(ktilt)=rw
               if( ibc .eq. 1 ) then
                 iunrad = 78
               else
                 iunrad = 79
               end if
               if( i .eq. 1 ) write(iunrad,*) it,nk,ktilt,rw
               write(iunrad,*) rzv(kn)
             else if(ibtyp.eq.32) then		!for malta 
c	       nothing
             else if(ibtyp.eq.0) then		!switched off
c	       nothing
             else if(ibtyp.lt.0) then		!closed
c	       nothing
	     else
c               kranf,krend not available...
c	       call zspeci(ibtyp,kranf,krend,rw)	!for radiation...
	       write(6,*) 'boundary = ',ibc,'   type = ',ibtyp
	       stop 'error stop sp111: Unknown boundary type'
	     end if

	  end do

	end do

c	-----------------------------------------------------
c	tilting
c	-----------------------------------------------------

	call tilt
	call z_tilt

c	-----------------------------------------------------
c	meteo forcing					!$$surel
c	-----------------------------------------------------

	call meteo_force

c	-----------------------------------------------------
c	set mass flux -> fills mfluxv and integrates to rqv
c	-----------------------------------------------------

	call set_mass_flux

c	-----------------------------------------------------
c	testing
c	-----------------------------------------------------

c	call tsbnds

c -----------------------------------------------------------
c end of routine
c -----------------------------------------------------------


	return
   95	continue
	write(6,*) 'One node boundary not allowed'
	write(6,*) 'Boundary :',ibc
	write(6,*) 'type     :',ibtyp
	stop 'error stop : sp111'
	end

c**************************************************************

	subroutine z_tilt

c tilting of boundary surface

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it

	integer irv(1)
	common /irv/irv
	real rzv(1)
	common /rzv/rzv
        real xgv(1),ygv(1)
        common /xgv/xgv,/ygv/ygv

	integer ibc,ibtyp,ktilt
	integer k,kranf,krend,kn1,kn2
	real dx,dy,ztilt,z
	double precision rltot,rltot1,rltot2,rl

	integer itybnd

	do ibc=1,nbc
          ibtyp = itybnd(ibc)
          call kanfend(ibc,kranf,krend)
	  call get_bnd_par(ibc,'ztilt',ztilt)
	  call get_bnd_ipar(ibc,'ktilt',ktilt)
	  if( ztilt .ne. 0. .and. ibtyp .eq. 1 ) then
	    rltot = 0.
	    rltot1 = 0.
	    do k=kranf+1,krend
		kn2=irv(k)
		kn1=irv(k-1)
		dx=xgv(kn2)-xgv(kn1)
		dy=ygv(kn2)-ygv(kn1)
	        rltot = rltot + sqrt(dx*dx+dy*dy)
		if( k .le. ktilt ) rltot1 = rltot
	    end do
	    rltot2 = rltot - rltot1
	    rl = 0.
            rzv(irv(kranf)) = rzv(irv(kranf)) - ztilt
	    do k=kranf+1,krend
		kn2=irv(k)
		kn1=irv(k-1)
		dx=xgv(kn2)-xgv(kn1)
		dy=ygv(kn2)-ygv(kn1)
	        rl = rl + sqrt(dx*dx+dy*dy)
		if( ktilt .le. 0 ) then			!no fixed node
	          z = - ztilt + (rl/rltot) * 2. * ztilt
		else
		  if( rl .le. rltot1 ) then
	            z = - ztilt + (rl/rltot1) * ztilt
		  else
	            z = ((rl-rltot1)/rltot2) * ztilt
		  end if
		end if
                rzv(kn2) = rzv(kn2) + z
	    end do
	  end if
	end do

	end

c**************************************************************

	subroutine tilt

c tilting of boundary surface

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it

	real grav,fcor,dcor,dirn,rowass,roluft
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
	real eps1,eps2,pi,flag,high,higi
	common /mkonst/ eps1,eps2,pi,flag,high,higi

	real rrv(1)
	common /rrv/rrv
	integer irv(1)
	common /irv/irv
	real rzv(1)
	common /rzv/rzv
	real rhv(1)
	common /rhv/rhv

        real znv(1)
        common /znv/znv
        real up0v(1), vp0v(1)
        common /up0v/up0v, /vp0v/vp0v

        real xgv(1),ygv(1)
        common /xgv/xgv,/ygv/ygv

        real tauxnv(1),tauynv(1)
        common /tauxnv/tauxnv,/tauynv/tauynv
        real ppv(1)
        common /ppv/ppv

	integer ibc,ibtyp,kranf,krend,ktilt,k,kn1,kn2
	real roinv,f,ginv,dx,dy,taux,tauy
	real u,v,z,h,p1,p2,b,hh,ztilt
	real getpar
	integer iround,itybnd

	roinv=1./rowass
	f=fcor
	ginv=1./grav

	do ibc=1,nbc
         ibtyp = itybnd(ibc)
         call kanfend(ibc,kranf,krend)
	 call get_bnd_ipar(ibc,'ktilt',ktilt)
	 call get_bnd_par(ibc,'ztilt',ztilt)

	 if( ztilt .ne. 0 ) then
		!nothing
	 else if(ktilt.gt.0.and.ibtyp.eq.1) then
	   do k=ktilt+1,krend,1
		kn2=irv(k)
		kn1=irv(k-1)
		dx=xgv(kn2)-xgv(kn1)
		dy=ygv(kn2)-ygv(kn1)
		taux=(tauxnv(kn1)+tauxnv(kn2))*.5
		tauy=(tauynv(kn1)+tauynv(kn2))*.5
		u=(up0v(kn1)+up0v(kn2))*.5
		v=(vp0v(kn1)+vp0v(kn2))*.5
		z=(znv(kn1)+znv(kn2))*.5
		h=(rhv(k-1)+rhv(k))*.5
		p1=ppv(kn1)
		p2=ppv(kn2)
		b=1./(h+z)
		hh=-roinv*(p2-p1)+b*(taux*dx+tauy*dy)+f*(v*dx-u*dy)
		rzv(kn2)=rzv(kn1)+ginv*hh
	   end do

	   do k=ktilt-1,kranf,-1
		kn2=irv(k+1)
		kn1=irv(k)
		dx=xgv(kn2)-xgv(kn1)
		dy=ygv(kn2)-ygv(kn1)
		taux=(tauxnv(kn1)+tauxnv(kn2))*.5
		tauy=(tauynv(kn1)+tauynv(kn2))*.5
		u=(up0v(kn1)+up0v(kn2))*.5
		v=(vp0v(kn1)+vp0v(kn2))*.5
		z=(znv(kn1)+znv(kn2))*.5
		h=(rhv(k+1)+rhv(k))*.5
		p1=ppv(kn1)
		p2=ppv(kn2)
		b=1./(h+z)
		hh=-roinv*(p2-p1)+b*(taux*dx+tauy*dy)+f*(v*dx-u*dy)
		rzv(kn1)=rzv(kn2)-ginv*hh
	   end do
	 end if
	end do

	end

c**************************************************************

	subroutine initilt(ibc)

c finds tilting node in boundary node list

	implicit none

	integer ibc

	logical berr
	integer kranf,krend
	integer ktilt,i
        integer kb

	integer ipext,iround,kbnd

	call get_bnd_ipar(ibc,'ktilt',ktilt)
	if(ktilt.le.0) return

        call kanfend(ibc,kranf,krend)

	berr = .true.
	do i=kranf,krend
           kb = kbnd(i)
	   if( kb .eq. ktilt ) then
		call set_bnd_ipar(ibc,'ktilt',i)
		berr = .false.
	   end if
	end do

	if( berr ) then
	  write(6,*) 'Node number for tilting not in boundary node list'
	  write(6,*) 'ktilt :',ipext(ktilt)
	  stop 'error stop : initilt'
	end if

	end

c******************************************************************

	subroutine iniflux(ibc)

c initializes flux boundary

	implicit none

        integer ibc

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v
	integer irv(1)
	common /irv/irv
	integer ierv(2,1)
	common /ierv/ierv
	real rrv(1)
	common /rrv/rrv
	real rlv(1), rhv(1)
	common /rlv/rlv, /rhv/rhv
	real hm3v(3,1)
	common /hm3v/hm3v
	real xgv(1),ygv(1)
	common /xgv/xgv,/ygv/ygv

	integer kranf,krend
	integer ie,i,k1,k2,kk1,kk2,ii1,ii2
	real fm,dx,dy,rl,h1,h2,fl
	integer ipext

        call kanfend(ibc,kranf,krend)

	if( krend-kranf .le. 0 ) return

	fm=0
	do i=kranf,krend
	   rrv(i)=0.
	   rhv(i)=0.
	end do

	do i=kranf,krend-1

	  k1=irv(i)
	  k2=irv(i+1)

	  do ie=1,nel
	   do ii1=1,3
	      ii2=mod(ii1,3)+1
	      kk1=nen3v(ii1,ie)
	      kk2=nen3v(ii2,ie)
	      if(k1.eq.kk1.and.k2.eq.kk2) goto 33
	   end do
	  end do

   33	  continue

	  if( k1 .ne. kk1 .or. k2 .ne. kk2 ) then
	    write(6,*) 'Cannot locate boundary nodes in element index'
	    write(6,*) 'node 1,node 2 :',ipext(k1),ipext(k2)
	    write(6,*) '(Are you sure that boundary nodes are given'
	    write(6,*) '   in anti-clockwise sense ?)'
	    stop 'error stop : iniflux'
	  end if

	  dx=xgv(k1)-xgv(k2)
	  dy=ygv(k1)-ygv(k2)
	  rl=sqrt(dx*dx+dy*dy)
	  h1=hm3v(ii1,ie)
	  h2=hm3v(ii2,ie)
	  fl=rl*(h1+h2)/2.
	  fm=fm+fl

	  rrv(i)=rrv(i)+rl*(2.*h1+h2)/6.
	  rrv(i+1)=rrv(i+1)+rl*(h1+2.*h2)/6.

	  rlv(i)=rl

	  rhv(i)=rhv(i)+h1
	  rhv(i+1)=rhv(i+1)+h2

	  ierv(1,i)=ie
	  ierv(2,i)=mod(ii2,3)+1

	end do

	do i=kranf,krend
	   rrv(i)=rrv(i)/fm
	end do
	do i=kranf+1,krend-1
	   rhv(i)=rhv(i)/2.
	end do

	end

c*********************************************************************

	subroutine rain0d( it )

c simulates rain (0D) increasing the water level

	implicit none

	integer it

	integer ndim
	parameter (ndim=100)
	real barray(ndim)
	save barray

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real metrain(1)
	common /metrain/metrain

	character*(80)	file
	integer k
	real rw
	real t,zdist

	real getpar

	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. -1 ) return

c---------------------------------------------------------------
c initialization
c---------------------------------------------------------------

	if( icall .eq. 0 ) then		!initialize
	   
	   zdist = getpar('zdist')
	   call getfnm('rain',file)

	   if( zdist .eq. 0. .and. file .eq. ' ' ) icall = -1
	   if( icall .eq. -1 ) return

	   call exffild(file,ndim,barray,zdist)

	   icall = 1

	end if

c---------------------------------------------------------------
c normal call
c---------------------------------------------------------------

	t = it
	call exfintp(barray,t,rw)

	do k=1,nkn
	  metrain(k) = rw
	end do

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c*******************************************************************

	subroutine rain2distributed

c adjustes distributed source with rain

	implicit none

	integer it

	real zconv
	parameter( zconv = 1.e-3 / 86400. )

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real metrain(1)
	common /metrain/metrain
	real rqdsv(1)
	common /rqdsv/rqdsv

	integer k

c---------------------------------------------------------------
c rain is in mm/day -> convert to m/s (water level change)
c---------------------------------------------------------------

	do k=1,nkn
	  rqdsv(k) = rqdsv(k) + metrain(k) * zconv
	end do

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c*******************************************************************

	subroutine convert_distributed

c converts distributed source from [m/s] to [m**3/s]

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real rqpsv(1), rqdsv(1)
	common /rqpsv/rqpsv, /rqdsv/rqdsv
	integer nen3v(3,1)
	common /nen3v/nen3v
	include 'ev.h'
	real v1v(1)
	common /v1v/v1v

	integer k,ie,ii
	real area3

	do k=1,nkn
	  v1v(k) = 0.
	end do

	do ie=1,nel
	  area3 = 4. * ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    v1v(k) = v1v(k) + area3
	  end do
	end do

	do k=1,nkn
	  rqdsv(k) = rqdsv(k) * v1v(k)
	end do

	end

c*******************************************************************

	subroutine evap_init

c initializes evaporation mass flux

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real evapv(1)
	common /evapv/evapv

	integer k

	do k=1,nkn
	  evapv(k) = 0.	
	end do

	end

c*******************************************************************

	subroutine evap_set

c adds evaporation mass flux to distributed source

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real evapv(1)
	common /evapv/evapv
	real rqdsv(1)
	common /rqdsv/rqdsv

	integer k,ievap
	real getpar

	ievap = nint(getpar('ievap'))
	if( ievap .le. 0 ) return

	do k=1,nkn
	  rqdsv(k) = rqdsv(k) - evapv(k)
	end do

	end

c*******************************************************************

	subroutine meteo_init

c initializes meteo variables

	integer imreg
	real getpar

	imreg = nint(getpar('imreg'))

	if( imreg .eq. 1 ) then
	  call meteo_regular
	else
	  call windad
	  call qflux_init
	end if

	call evap_init

	end

c*******************************************************************

	subroutine meteo_force

c update meteo variables and admin rain/evaporation

	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it

	integer imreg
	real getpar

	imreg = nint(getpar('imreg'))

	if( imreg .eq. 1 ) then
	  call meteo_regular
	else
	  call windad			!wind
	  call qflux_read(it)		!heat flux
	  call rain0d(it)			!rain
	end if

	call rain2distributed		!copy rain to distributed source
	call evap_set			!add evaporation
	call convert_distributed	!convert from [m/s] to [m**3/s]

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine b3dvalue(ibc,it,nsize,ndim,array,rw)

c computes z value for boundary ibc

	implicit none

	integer ibc		!boundary number
	integer it		!time for which value is needed
        integer nsize		!total number of values needed
        integer ndim		!dimension for following arrays
        real array(1)		!array with information on boundary
        real rw(1)		!array with results

	integer iqual,intpol,nbunit,ibtyp
	integer iaux
	integer nvar,ndata
	integer i
	real rit,zfact
	real rw0
        character*80 file

	logical debug

	debug = .false.
	!if( ibc .eq. 13 ) debug = .true.

	if( debug ) write(91,*) 'entering... b3dvalue (before): ',ibc

	call get_bnd_par(ibc,'zfact',zfact)
	call get_bnd_ipar(ibc,'iqual',iqual)

	nvar = 1
	ndata = nvar*max(nsize,1)

        call exfunit(array,nbunit)

	if( nbunit .eq. 0 ) then
		call get_boundary_file(ibc,'zeta',file)
		call get_bnd_ipar(ibc,'ibtyp',ibtyp)
		call get_bnd_ipar(ibc,'intpol',intpol)
		if( debug ) then
		  write(6,*) 'b3dvalue: ',ibc,intpol,nsize,ndim,file
		end if
		if( intpol .le. 0 ) then
		  intpol = 2
		  if( ibtyp .eq. 1 ) intpol = 4
		  write(6,*) 'intpol set: ',ibc,ibtyp,intpol
		end if
                call exffil(file,intpol,nvar,nsize,ndim,array)
	        call exfunit(array,nbunit)	!need unit
		write(6,*) 'intpol final set: ',ibc,ibtyp,intpol
		call set_bnd_ipar(ibc,'intpol',intpol)
		if( debug ) then
		  write(91,*) 'initializing z boundary ',ibc
		  write(91,*) 'file :',file
                  call exfinfo(91,array)
		end if
	end if

	!write(6,*) 'entering... b3dvalue (after): ',ibc

	rit = it

	if( iqual .ne. 0 ) then		!iqual is broken
		stop 'error stop b3dvalue: iqual is broken'
	else if( nbunit .gt. 0 ) then	!read boundary files
		call exfintp(array,rit,rw)
                do i=1,ndata
		  rw(i) = rw(i) * zfact
                end do
		if( debug ) then
		  !write(86,*) it,ibc,(rw(i),i=1,nvar)
		  !write(88,*) it,ibc,(rw(i),i=1,nvar)
		  write(91,*) 'z boundary ',ibc
                  call exfinfo(91,array)
		end if
	else				!constant b.c.
		call get_oscil(ibc,rit,rw0)
		do i=1,ndata
		  rw(i) = rw0 * zfact
		end do
	end if

	if( debug ) write(91,*) 'z boundary ',ibc,nvar,(rw(i),i=1,nvar)

	return
	end

c**********************************************************************

	subroutine set_mass_flux

c sets up (water) mass flux array mfluxv (3d) and rqv (vertically integrated)

	implicit none

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        integer ilhkv(1)
        common /ilhkv/ilhkv
	real rqv(1)
	common /rqv/rqv
	real rqpsv(1), rqdsv(1)			![m**3/s]
	common /rqpsv/rqpsv, /rqdsv/rqdsv
	real mfluxv(nlvdim,1)			![m**3/s]
	common /mfluxv/mfluxv

	logical debug
	integer i,k,l,lmin,lmax,nk,ibc,mode
	integer ibtyp,levmax,levmin
	real flux,vol,voltot,fluxtot,fluxnode
	real vols(nkndim)

	integer nkbnds,kbnds
	real volnode		!function to compute volume of node

c------------------------------------------------------------------
c initialize arrays and parameter
c------------------------------------------------------------------

	mode = -1		!old time step
	debug = .true.
	debug = .false.

	do k=1,nkn
	  do l=1,nlv
	    mfluxv(l,k) = 0.
	  end do
	end do

c------------------------------------------------------------------
c loop over boundaries for point sources -> use volumes as weight
c------------------------------------------------------------------

	do ibc=1,nbc

          nk = nkbnds(ibc)
	  call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	  call get_bnd_ipar(ibc,'levmin',levmin)
	  call get_bnd_ipar(ibc,'levmax',levmax)

	  if( ibtyp .lt. 2 .or. ibtyp .gt. 3 ) nk = 0		!skip

	  if(debug) then
	    write(6,*) 'computing mass flux: ',ibc,ibtyp,nk,levmax
	  end if

	  do i=1,nk
            k = kbnds(ibc,i)
	    lmax = ilhkv(k)
	    if( levmax .gt. 0 ) lmax = min(lmax,levmax)
	    lmin = 1
	    if( levmin .gt. 0 ) lmin = max(lmin,levmin)

	    if( lmin .gt. lmax ) goto 98

	    voltot = 0.
	    do l=lmin,lmax
	      vol = volnode(l,k,mode)
	      vols(l) = vol
	      voltot = voltot + vol
	    end do

	    if( voltot .le. 0. ) goto 99

	    flux = rqpsv(k)
	    if(debug) write(6,*) '   ',k,lmin,lmax,flux,voltot
	    do l=lmin,lmax
	      mfluxv(l,k) = flux * vols(l) / voltot
	    end do
	  end do

	end do

c------------------------------------------------------------------
c add distributed sources
c------------------------------------------------------------------

	do k=1,nkn
	  mfluxv(1,k) = mfluxv(1,k) + rqdsv(k)	!rain, evaporation
	  !lmax = ilhkv(k)
	  !mfluxv(lmax,k) = gwf		!here distributed ground water flow
	end do

c------------------------------------------------------------------
c compute total flux for check and integrate flux into rqv
c------------------------------------------------------------------

	fluxtot = 0.
	do k=1,nkn
	  fluxnode = 0.
	  lmax = ilhkv(k)
	  do l=1,lmax
	    flux = mfluxv(l,k)
	    !if( debug .and. flux .gt. 0. ) write(6,*) '  flux: ',k,l,flux
	    fluxnode = fluxnode + flux
	  end do
	  rqv(k) = fluxnode
	  fluxtot = fluxtot + fluxnode
	end do

	if( debug ) write(6,*) '  total flux: ',fluxtot

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	return
   98	continue
	write(6,*) 'lmin > lmax'
   99	continue
	write(6,*) 'ibc = ',ibc
	write(6,*) 'i = ',i
	write(6,*) 'k = ',k
	write(6,*) 'ilhkv(k) = ',ilhkv(k)
	write(6,*) 'levmin = ',levmin
	write(6,*) 'levmax = ',levmax
	write(6,*) 'lmin = ',lmin
	write(6,*) 'lmax = ',lmax
	stop 'error stop set_mass_flux: voltot = 0'
	end

c**********************************************************************

	subroutine adjust_mass_flux

c adjusts mass flux for dry nodes

	implicit none

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        integer inodv(1)
        common /inodv/inodv
	real rqv(1)
	common /rqv/rqv
	real rqpsv(1), rqdsv(1)
	common /rqpsv/rqpsv, /rqdsv/rqdsv
	real mfluxv(nlvdim,1)
	common /mfluxv/mfluxv

	integer k,l

        logical iskout
        iskout(k) = inodv(k).eq.-2

	do k=1,nkn
	  if( iskout(k) ) then
	    do l=1,nlv
	      mfluxv(l,k) = 0.
	    end do
	    rqv(k) = 0.
	    rqdsv(k) = 0.
	  end if
	end do

	end

c**********************************************************************

	subroutine make_scal_flux(what,r3v,scal,sflux,sconz,ssurf)

c computes scalar flux from fluxes and concentrations

	implicit none

	include 'param.h'

	character*(*) what
	real r3v(nlvdim,1)	!concentration for boundary condition
	real scal(nlvdim,1)	!concentration of scalar
	real sflux(nlvdim,1)	!mass flux for each finite volume (return)
	real sconz(nlvdim,1)	!concentration for each finite volume (return)
	real ssurf		!value of scalar for surface flux

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real eps1,eps2,pi,flag,high,hihi
        common /mkonst/ eps1,eps2,pi,flag,high,hihi
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        integer ilhkv(1)
        common /ilhkv/ilhkv
	real rqpsv(1), rqdsv(1)
	common /rqpsv/rqpsv, /rqdsv/rqdsv
	real mfluxv(nlvdim,1)
	common /mfluxv/mfluxv

	integer k,l,lmax,ks
	real flux,conz
	real surf_flux
	real getpar

	ks = 2827
	ks = 2757
	ks = 2831
	ks = -1
	!ks = nint(getpar('kref'))	!not working - here global, but local

	do k=1,nkn
	  lmax = ilhkv(k)
	  surf_flux = rqdsv(k)
	  do l=1,lmax
	    sflux(l,k) = 0.
	    sconz(l,k) = 0.
	    flux = mfluxv(l,k)
	    if( l .eq. 1 ) flux = flux - surf_flux	!without surface flux
	    conz = r3v(l,k)
	    if( flux .ne. 0. .or. conz .ne. flag ) then
	      if( flux .ne. 0. .and. conz .eq. flag ) goto 99
	      if( conz .le. -990. ) conz = scal(l,k)	!ambient value
	      if( conz .le. -5555. ) conz = scal(l,ks) - 10000. - conz !diff
	      if( flux .lt. 0. ) conz = scal(l,k)	!ambient value (BUG)
	      sflux(l,k) = flux * conz
	      sconz(l,k) = conz			!bug fix - was out of loop
	    end if
	  end do
	  conz = ssurf
	  if( ssurf .le. -990 ) conz = scal(1,k)
	  if( ssurf .le. -5555 ) conz = scal(1,k) - 10000. - ssurf !diff
	  sflux(1,k) = sflux(1,k) + surf_flux * conz
	  ! next should be sconz(1,k) = conz if surf_flux is eliminated
	  sconz(1,k) = 0.
	  if( mfluxv(1,k) .ne. 0 ) sconz(1,k) = sflux(1,k) / mfluxv(1,k)
	end do

	return
   99	continue
	write(6,*) 'what: ',what
	write(6,*) 'k,l: ',k,l
	write(6,*) 'flux,conz: ',flux,conz
	stop 'error stop make_scal_flux: boundary condition mismatch'
	end

c**********************************************************************

	subroutine check_scal_flux(what,scal,sconz)

c checks scalar flux

	implicit none

	include 'param.h'

	character*(*) what
	real scal(nlvdim,1)	!concentration of scalar
	real sconz(nlvdim,1)	!concentration for each finite volume

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real eps1,eps2,pi,flag,high,hihi
        common /mkonst/ eps1,eps2,pi,flag,high,hihi
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        integer ilhkv(1)
        common /ilhkv/ilhkv
	real mfluxv(nlvdim,1)
	common /mfluxv/mfluxv

	integer k,l,lmax,ks
	real cconz,qflux,mflux

	write(46,*) 'check_scal_flux ',what,it

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cconz = sconz(l,k)         !concentration has been passed
            qflux = mfluxv(l,k)
            if( qflux .lt. 0. ) cconz = scal(l,k)
            mflux = qflux * cconz
	    if( qflux .ne. 0 ) then
	      write(46,1000) k,l,mflux,qflux,cconz,scal(l,k)
	    end if
          end do
        end do

	return
 1000	format(2i10,4f10.4)
	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine init_scal_bc(r3v)

c initializes array for scalar boundary condition

	implicit none

	include 'param.h'

	real r3v(nlvdim,nkndim)

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /level/ nlvdi,nlv
	real eps1,eps2,pi,flag,high,higi
	common /mkonst/ eps1,eps2,pi,flag,high,higi

	integer k,l

	do k=1,nkn
	  do l=1,nlv
	    r3v(l,k) = flag
	  end do
	end do

	end

c*******************************************************************

	subroutine mult_scal_bc(r3v,value)

c multiplies array for scalar boundary condition with value

	implicit none

	include 'param.h'

	real r3v(nlvdim,nkndim)
	real value

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /level/ nlvdi,nlv
	real eps1,eps2,pi,flag,high,higi
	common /mkonst/ eps1,eps2,pi,flag,high,higi

	integer k,l

	do k=1,nkn
	  do l=1,nlv
	    if( r3v(l,k) .ne. flag ) r3v(l,k) = r3v(l,k) * value
	  end do
	end do

	end


c*******************************************************************

	subroutine dist_3d(nlvdim,r3v,kn,nbdim,values)

	implicit none

	integer nlvdim
	real r3v(nlvdim,1)
	integer kn
	integer nbdim
	real values(1)

	integer l,lmax

	if( nbdim .eq. 0 ) then
	  lmax = 1
	else
	  lmax = min(nbdim,nlvdim)
	end if

	do l=1,lmax
	  r3v(l,kn) = values(l)
	end do
	  
	do l=lmax+1,nlvdim
	  !r3v(l,kn) = values(nbdim)	!BUGFIX
	  r3v(l,kn) = r3v(lmax,kn)
	end do
	  
	end

c**********************************************************************

	subroutine dist_horizontal(nlvdim,r3v,n,value)

	implicit none

	integer nlvdim
	real r3v(nlvdim,1)
	integer n
	real value

	integer k

	do k=1,n
	  r3v(1,k) = value
	end do
	  
	end

c**********************************************************************

        subroutine aver_horizontal(nlvdim,r3v,n,value)

        implicit none

        integer nlvdim
        real r3v(nlvdim,1)
        integer n
        real value

        integer k

	value = 0.
        do k=1,n
          value = value + r3v(1,k)
        end do
	value = value / n

	end

c**********************************************************************

	subroutine print_scal_bc(r3v)

c prints non-flag entries of array for scalar boundary condition

	implicit none

	include 'param.h'

	real r3v(nlvdim,nkndim)

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /level/ nlvdi,nlv
	real eps1,eps2,pi,flag,high,higi
	common /mkonst/ eps1,eps2,pi,flag,high,higi

	integer k,l
	real value

	do k=1,nkn
	  do l=1,nlv
	    value = r3v(l,k)
	    if( value .ne. flag ) then
		write(6,*) 'print_scal_bc: ',k,l,value
	    end if
	  end do
	end do

	end

c**********************************************************************

	subroutine get_bflux(k,flux)

c returns boundary flux of node k

	implicit none

	integer k	!node
	real flux	!flux in node k (return)

	real rqv(1)
	common /rqv/rqv

	flux = rqv(k)

	end

c**********************************************************************
