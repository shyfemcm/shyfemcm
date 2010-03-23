c
c $Id: lagrange_cont.f,v 1.5 2009-09-14 08:20:57 georg Exp $
c
c simulates continuous release over open boundaries
c
c revision log :
c
c 12.12.2007    ggu	written from scratch
c 12.06.2008    ggu	initialize also z
c 28.08.2009    ggu	new call to find_elems_to_segment (before line_elems)
c
c*******************************************************************

	subroutine lagr_continuous_release_shell

c manages continuous release of particles

	implicit none

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	real pps
	integer itmin,itmax

	pps = 0.5
	pps = 0.
	itmin = 86400
	itmax = 3*86400

	if( pps .le. 0. ) return

	if( it .ge. itmin .and. it .le. itmax ) then
	  call lagr_continuous_release(pps)
	end if

	end

c*******************************************************************

	subroutine lagr_continuous_release(pps)

c continuous release

	implicit none

	real pps	!particles per second

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer k1,k2
	integer ibc,nk,i,ibtyp,np
	real totdist,dxy,part,dt

	integer nkbnds,kbnds
	real dist_node

	call get_timestep(dt)

	part = pps*dt
	if( part .le. 0. ) return

	call compress_particles

	do ibc=1,nbc

	  nk = nkbnds(ibc)
	  call get_bnd_ipar(ibc,'ibtyp',ibtyp)

	  if( ibtyp .eq. 1 ) then	!only for level boundaries

	    totdist = 0.
	    do i=2,nk
	      k1 = kbnds(ibc,i-1)
	      k2 = kbnds(ibc,i)
	      dxy = dist_node(k1,k2)
	      totdist = totdist + dxy
	    end do

	    do i=2,nk
	      k1 = kbnds(ibc,i-1)
	      k2 = kbnds(ibc,i)
	      dxy = dist_node(k1,k2)
	      np = nint(part*dxy/totdist)
	      call create_parts(np,k1,k2)
	    end do

	  end if

	end do

	end

c*******************************************************************

	subroutine create_parts(np,k1,k2)

	implicit none

	integer np,k1,k2

	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv

	integer i,ie1,ie2
	real x1,y1,x2,y2,dx,dy
	real rl,x,y

	real ggrand

	if( np .le. 0 ) return

	x1 = xgv(k1)
	y1 = ygv(k1)
	x2 = xgv(k2)
	y2 = ygv(k2)

	dx = x1 - x2
	dy = y1 - y2

	call find_elems_to_segment(k1,k2,ie1,ie2)
	if( ie1 .eq. 0 .or. ie2 .ne. 0 ) then
	  write(6,*) k1,k2,ie1,ie2
	  stop 'error stop create_parts: error in boundary'
	end if
	  
	write(6,*) 'create_parts: ',np,k1,k2,ie1

	do i=1,np
	  rl = ggrand(77)
	  x = x1 + rl*dx
	  y = y1 + rl*dy
	  call insert_particle(ie1,x,y)
	end do

	end

c*******************************************************************
