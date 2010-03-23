c
c $Id: bastreat.f,v 1.14 2010-02-26 17:35:06 georg Exp $
c
c revision log :
c
c 06.04.1999    ggu     completely restructured
c 04.06.1999    ggu     new statistics are computed
c 08.09.2003    ggu     mode 5 -> write depth values from elements
c 23.09.2004    ggu     interpolq() changed for bathy interpolation
c 02.10.2004    ggu     interpole() for exponential interpolation
c 01.11.2004    ggu     whole program simplyfied
c 06.12.2008    ggu     smoothing introduced
c 06.04.2005    ggu     read param.h
c 29.05.2009    ggu     does only depth limiting and smoothing
c 20.11.2009    ggu     possibility to smooth only on specific areas
c
c****************************************************************

        program bastreat

c performs modifications on basin
c
c takes care of lat/lon coordinates

	implicit none

	include 'param.h'

	integer ndim

        character*80 descrp
        common /descrp/ descrp

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        real xgv(nkndim), ygv(nkndim)
        real hm3v(3,neldim)
        common /xgv/xgv, /ygv/ygv
        common /hm3v/hm3v

        integer nen3v(3,neldim)
        integer ipev(neldim), ipv(nkndim)
        integer iarv(neldim)
        common /nen3v/nen3v
        common /ipev/ipev, /ipv/ipv
        common /iarv/iarv

        real hev(neldim)
        common /hev/hev
        real hkv(nkndim)
        common /hkv/hkv

        real raux(neldim)
        integer iaux(neldim)
        integer ipaux(nkndim)

	include 'evmain.h'

        character*40 bfile,gfile,nfile
        character*60 line
	integer node,nit
	integer mode,np,n,niter,i
        integer ner,nco,nknh,nelh,nli
	integer nlidim,nlndim
	integer ike,idepth
	real hmin,hmax
	real f(5)
	logical bstop

	integer iscanf

	hmin = -99999.
	hmax = 99999.

c-----------------------------------------------------------------
c what to do
c-----------------------------------------------------------------

        write(6,*)
        write(6,*) 'I need the name of the grid file '
        write(6,*)
	write(6,*) 'Enter file name: '
	read(5,'(a)') gfile
        if( gfile .eq. ' ' ) stop
	write(6,*) 'grid is read from file : ', gfile
        write(6,*)

        write(6,*)
        write(6,*) 'Limiting values for depth (<CR> for no limit):'
        write(6,*)
	write(6,*) 'Enter hmin,hmax: '
	read(5,'(a)') line
	n = iscanf(line,f,2)
	if( n .lt. 0 ) goto 97
	if( n .gt. 0 ) hmin = f(1)
	if( n .gt. 1 ) hmax = f(2)
	write(6,*) 'hmin,hmax :',hmin,hmax
        write(6,*)

        write(6,*)
        write(6,*) 'Number of iterations for smoothing:'
        write(6,*)
	write(6,*) 'Enter niter (<CR> for no smoothing): '
	read(5,'(i10)') niter
	write(6,*) 'niter :',niter
        write(6,*)

	if( niter .gt. 0 ) then

        write(6,*)
        write(6,*) 'Enter parameters for smoothing:'
        write(6,*)
        write(6,*) 'Either alpha or a1,h1,a2,h2'
        write(6,*)
	write(6,*) 'Enter parameters: '
	read(5,'(a)') line
	n = iscanf(line,f,4)
	if( n .ne. 1 .and. n .ne. 4 ) goto 96
	if( n .eq. 1 ) then
	  f(3) = f(1)
	  f(2) = 0.
	  f(4) = 10000.
	end if
	write(6,*) 'parameters used :',(f(i),i=1,4)
        write(6,*)

	end if

c-----------------------------------------------------------------
c read in bathymetry file
c-----------------------------------------------------------------

	!np = ndim
	!call readbat(bfile,np,xp,yp,dp)

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        ner = 6
        bstop = .false.

        nlidim = 0
        nlndim = 0
        call rdgrd(
     +                   gfile
     +                  ,bstop
     +                  ,nco,nkn,nel,nli
     +                  ,nkndim,neldim,nlidim,nlndim
     +                  ,ipv,ipev,iaux
     +                  ,iaux,iarv,iaux
     +                  ,hkv,hev,raux
     +                  ,xgv,ygv
     +                  ,nen3v
     +                  ,iaux,iaux
     +                  )

        if( bstop ) stop 'error stop rdgrd'

c        call rdgrd(gfile,ner,bstop,nco,nkn,nknh,nel,nelh,nli
c     +                  ,nkndim,neldim,nliread
c     +                  ,ipv,ipev,ianv,iarv,nen3v,xgv,ygv,hev,hkv)
c
c        call ex2in(nkn,nel,ipv,ipaux,nen3v,ner,bstop)

        call ex2in(nkn,3*nel,nlidim,ipv,ipaux,nen3v,iaux,bstop)
        if( bstop ) stop 'error stop ex2in'

c-----------------------------------------------------------------
c handling depth and coordinates
c-----------------------------------------------------------------

	call check_coords                       !sets lat/lon flag
	call set_depth(nknh,nelh)

	ike = 1
	if( nknh .gt. 0 ) ike = 2
	if( nknh .gt. 0 .and. nknh .ne. nkn ) goto 99
	if( nelh .gt. 0 .and. nelh .ne. nel ) goto 99
	if( nknh .eq. 0 .and. nelh .eq. 0 ) goto 99
	if( nknh .gt. 0 .and. nelh .gt. 0 ) goto 99

c-----------------------------------------------------------------
c general info
c-----------------------------------------------------------------

        write(6,*)
        write(6,*) ' nkn  = ',nkn, '  nel  = ',nel
        write(6,*) ' nknh = ',nknh,'  nelh = ',nelh
        write(6,*)

c-----------------------------------------------------------------
c node_test
c-----------------------------------------------------------------

	call node_test
	call sp110a

c-----------------------------------------------------------------
c smooth
c-----------------------------------------------------------------

	call limit_depth(ike,hmin,hmax)

	if( niter .gt. 0 ) then
	  call smooth_bathy(ike,niter,f)
	end if

	call copy_depth(ike)	!copyt to nodes/elements

c-----------------------------------------------------------------
c write
c-----------------------------------------------------------------

        nfile = 'bastreat.grd'
        open(1,file=nfile,status='unknown',form='formatted')
        call wrgrd(1,ike)
        close(1)
        write(6,*) 'file has been written to ',nfile

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	stop
   96	continue
	write(6,*) n,(f(i),i=1,n)
	write(6,*) 'there must be either 1 or 4 parameters'
	stop 'error stop bastreat: error in smoothing parameters'
   97	continue
	write(6,*) line
	stop 'error stop bastreat: read error'
   99	continue
	write(6,*) 'nelh,nknh: ',nelh,nknh
	stop 'error stop bastreat: error in parameters'
	end

c*******************************************************************

	subroutine copy_depth(ike)

c copies depth values from elems/nodes to nodes/elems

	implicit none

	integer ike

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,neldim)
        common /nen3v/nen3v

        real hev(neldim)
        common /hev/hev
        real hkv(nkndim)
        common /hkv/hkv

	integer k,ie,ii
	real depth
	integer ic(nkndim)

	if( ike .eq. 1 ) then

	  do k=1,nkn
	    ic(k) = 0
	    hkv(k) = 0.
	  end do

	  do ie=1,nel
	    do ii=1,3
	      k = nen3v(ii,ie)
	      hkv(k) = hkv(k) + hev(ie)
	    end do
	  end do

	  do k=1,nkn
	    hkv(k) = hkv(k) / ic(k)
	  end do

	else

	  do ie=1,nel
	    depth = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      depth = depth + hkv(k)
	    end do
	    hev(ie) = depth / 3.
	  end do

	end if

	end

c*******************************************************************

	subroutine set_depth(nknh,nelh)

c handles depth values

	implicit none

	integer idepth
	integer nknh,nelh

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real hev(neldim)
        common /hev/hev
        real hkv(nkndim)
        common /hkv/hkv

	integer k,ie

	nknh = 0
	nelh = 0

	do k=1,nkn
	  if( hkv(k) .gt. -990 ) nknh = nknh + 1
	end do
	do ie=1,nel
	  if( hev(ie) .gt. -990 ) nelh = nelh + 1
	end do

	end

c*******************************************************************

        function areatr(ie)

c determination of area of element
c
c ie            number of element (internal)
c areatr        element area (return value)

	real areatr
	integer ie

	include 'param.h'

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	real aj
	integer ii,i1,i2,k1,k2

        aj=0.
        do ii=1,3
          i1=mod(ii,3)+1
          i2=mod(i1,3)+1
          k1=nen3v(i1,ie)
          k2=nen3v(i2,ie)
          aj=aj+xgv(k1)*ygv(k2)-xgv(k2)*ygv(k1)
        end do

        areatr = aj / 2.

        end

c*******************************************************************

	subroutine wrgrd(iunit,ike)

c writes grd file from bas

	implicit none

	integer iunit
	integer ike

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

        integer nen3v(3,neldim)
        integer ipev(neldim), ipv(nkndim)
        integer iarv(neldim)
        common /nen3v/nen3v
        common /ipev/ipev, /ipv/ipv
        common /iarv/iarv

        real hev(neldim)
        common /hev/hev
        real hkv(nkndim)
        common /hkv/hkv

	integer k,ie,ii

	do k=1,nkn
	  if( ike .eq. 1 ) then
	    write(iunit,1000) 1,ipv(k),0,xgv(k),ygv(k)
	  else
	    write(iunit,1000) 1,ipv(k),0,xgv(k),ygv(k),hkv(k)
	  end if
	end do

	write(iunit,*)

	do ie=1,nel
	  if( ike .eq. 1 ) then
	    write(iunit,1100) 2,ipev(ie),iarv(ie)
     +		,3,(ipv(nen3v(ii,ie)),ii=1,3),hev(ie)
	  else
	    write(iunit,1100) 2,ipev(ie),iarv(ie)
     +		,3,(ipv(nen3v(ii,ie)),ii=1,3)
	  end if
	end do

	return
 1000	format(i1,2i10,3e16.8)
 1100	format(i1,2i10,i4,3i10,e16.8)
	end

c*******************************************************************

	subroutine node_test

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	integer ie,ii,k

	write(6,*) 'node_testing ... ',nel,nkn
	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( k .le. 0 ) write(6,*) ie,ii,k
	    iii = mod(ii,3) + 1
	    k1 = nen3v(iii,ie)
	    if( k .eq. k1 ) write(6,*) ie,(nen3v(iii,ie),iii=1,3)
	  end do
	end do
	write(6,*) 'end of node_testing ... '

	end

c*******************************************************************

	subroutine triab(x,y,area,x0,y0)

c computes area and center point of triangle

	implicit none

	real x(3)
	real y(3)
	real area,x0,y0

	area = 0.5 * ( (x(2)-x(1))*(y(3)-y(1)) 
     +			- (x(3)-x(1))*(y(2)-y(1)) )

	x0 = (x(1)+x(2)+x(3))/3.
	y0 = (y(1)+y(2)+y(3))/3.

	end

c*******************************************************************

	subroutine prepare_on_element(nt,xt,yt,at,ht)

	implicit none

	integer nt
	real xt(1)
	real yt(1)
	real at(1)
	real ht(1)

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        real hm3v(3,neldim)
        common /hm3v/hm3v
	real hev(neldim)
        common /hev/hev
	real hkv(nkndim)
        common /hkv/hkv

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	integer ie,ii,k
	real area,x0,y0
	real xaux(3),yaux(3)

	nt = nel

	do ie=1,nel

	  do ii=1,3
	    k = nen3v(ii,ie)
	    xaux(ii) = xgv(k)
	    yaux(ii) = ygv(k)
	  end do
	  call triab(xaux,yaux,area,x0,y0)

	  xt(ie) = x0
	  yt(ie) = y0
	  at(ie) = area
	  ht(ie) = hev(ie)

	end do

	end

c*******************************************************************

	subroutine prepare_on_node(nt,xt,yt,at,ht)

	implicit none

	integer nt
	real xt(1)
	real yt(1)
	real at(1)
	real ht(1)

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        real hm3v(3,neldim)
        common /hm3v/hm3v
	real hev(neldim)
        common /hev/hev
	real hkv(nkndim)
        common /hkv/hkv

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	integer ie,ii,k
	real area,x0,y0
	real xaux(3),yaux(3)

	nt = nkn

	do k=1,nkn
	  at(k) = 0.
	  ht(k) = 0.
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    xaux(ii) = xgv(k)
	    yaux(ii) = ygv(k)
	  end do
	  call triab(xaux,yaux,area,x0,y0)
	  at(k) = at(k) + area
	  ht(k) = ht(k) + 1.
	end do

	do k=1,nkn
	  xt(k) = xgv(k)
	  yt(k) = ygv(k)
	  at(k) = at(k) / (3.*ht(k))
	  ht(k) = hkv(k)
	end do

	end

c*******************************************************************

        subroutine check_coords

c checks if coordinates are lat/lon

        implicit none

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

        integer k,isphe
        real xmin,xmax,ymin,ymax

        xmin = xgv(1)
        xmax = xgv(1)
        ymin = ygv(1)
        ymax = ygv(1)

        do k=1,nkn
          xmin = min(xmin,xgv(k))
          xmax = max(xmax,xgv(k))
          ymin = min(ymin,ygv(k))
          ymax = max(ymax,ygv(k))
        end do

        isphe = 1
        if( xmin .lt. -180. .or. xmax .gt. 360. ) isphe = 0
        if( ymin .lt. -180. .or. ymax .gt. 180. ) isphe = 0
        !call set_dist(isphe)

        call set_coords_ev(isphe)
        write(6,*) 'setting for coordinates: ',isphe
        if( isphe .ne. 0 ) write(6,*) 'using lat/lon coordinates'

        end

c*******************************************************************

	subroutine limit_depth(ike,hmin,hmax)

	implicit none

	integer ike
	real hmin,hmax

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real hev(neldim)
        common /hev/hev
	real hkv(nkndim)
        common /hkv/hkv

	integer ie,k

	if( ike .eq. 1 ) then		!elementwise

	  do ie=1,nel
	    if( hev(ie) .gt. hmax ) hev(ie) = hmax
	    if( hev(ie) .lt. hmin ) hev(ie) = hmin
	  end do

	else				!nodewise

	  do k=1,nkn
	    if( hkv(k) .gt. hmax ) hkv(k) = hmax
	    if( hkv(k) .lt. hmin ) hkv(k) = hmin
	  end do

	end if

	end

c*******************************************************************

	subroutine smooth_bathy(ike,niter,f)

c smoothes depth values

	implicit none

	integer ike
	integer niter
	real f(4)

	if( ike .eq. 1 ) then
	  call smooth_bathy_elem(niter,f)
	else
	  stop 'error stop: cannot yet smooth on nodes'
	end if

	end

c*******************************************************************

	subroutine smooth_bathy_elem(niter,f)

c smoothes depth values

	implicit none

	integer niter
	real f(4)

	include 'param.h'
	include 'evmain.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
	real hev(neldim)
        common /hev/hev
	real hkv(nkndim)
        common /hkv/hkv

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	integer ie,ii,k,i,nok
	real x,y,d,h,hold,hnew,ao
	real hmin
	real alpha,beta
        real h1,h2,a1,a2
	integer iaux,inum,ityp
	real xt(3), yt(3)
	real v1v(nkndim)
	integer ihev(neldim)
	logical inconvex

c--------------------------------------------------------------
c set parameters
c--------------------------------------------------------------

	a1 = f(1)
	h1 = f(2)
	a2 = f(3)
	h2 = f(4)

	write(6,*) 'smoothing bathymetry: ',niter
	write(6,*) '  params: ',(f(i),i=1,4)

c--------------------------------------------------------------
c iterate over smoother
c--------------------------------------------------------------

	do i=1,niter

c	  -----------------------------------------------
c	  compute values at nodes (averages of element)
c	  -----------------------------------------------

	  do k=1,nkn
	    hkv(k) = 0.
	    v1v(k) = 0.
	  end do

	  do ie=1,nel
	    ao = ev(10,ie)
	    do ii=1,3
	      k = nen3v(ii,ie)
	      hkv(k) = hkv(k) + ao * hev(ie)
	      v1v(k) = v1v(k) + ao
	    end do
	  end do

	  do k=1,nkn
	    hkv(k) = hkv(k) / v1v(k)
	  end do

c	  -----------------------------------------------
c	  average to element
c	  -----------------------------------------------

	  nok = 0

	  do ie=1,nel
	    h = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      h = h + hkv(k)
	    end do
	    hnew = h / 3.
	    hold = hev(ie)

            beta = (hold-h1)/(h2-h1)
            beta = max(beta,0.)
            beta = min(beta,1.)

            alpha = a1 + beta * (a2-a1)

	    !call coords_ok(ie,alpha)	!customize to smooth on specific areas
	    if( alpha .gt. 0. ) nok = nok + 1
            !write(6,*) ie,hold,alpha

	    hev(ie) = (1.-alpha) * hold + alpha * hnew
	  end do

c	  -----------------------------------------------
c	  write to terminal
c	  -----------------------------------------------

	  write(6,*) 'pass ',i,' of ',niter,'  elements smoothed ',nok

	end do

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

	subroutine coords_ok(ie,alpha)

	implicit none

	integer ie
	real alpha

	real x,y,x1,y1,x2,y2

	x1 = 103242.
	y1 = 62770.
	x2 = 104133.
	y2 = 64304.

	x1 = 0.
	y1 = 0.
	x2 = 10000.
	y2 = 10000.

	call baric(ie,x,y)

	if( x .ge. x1 .and. x .le. x2 ) then
	  if( y .ge. y1 .and. y .le. y2 ) then
	    return				!ok - keep alpha
	  end if
	end if

	alpha = 0.0

	end

c*******************************************************************
