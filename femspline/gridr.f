c
c smoothing program
c
c contents :
c
c revision log :
c 
c 01.02.2002    ggu     handle periodic line
c 06.08.2002    ggu     read also depth value
c 23.09.2004    ggu     adapted for malta, bug fix
c 19.10.2005    ggu     documentation and description
c
c description :
c
c This program smothes and reduces the number of points on a line
c
c The program ask for the GRD file that contains the lines
c It also asks for the smoothing parameter (sigma) and reduction
c	parameter (reduct)
c The units of these values must be in the units of the line coordinates.
c 
c Smoothing:
c
c The program uses a gaussian curve as the smoothing operator, with
c its standard deviation given by sigma. Therefore, if (x,y) coordinates
c of the lines are in meters and sigma=100, the lengthscale for the 
c smoothing is about 100 meters, which means that points closer than 100
c meters to the point will have much more effect than points further away
c than 100 m.
c
c Reduction:
c
c The program uses the parameter reduct to decide if a point should
c be kept or not. If reduct=100, than all points closer than 100 m
c to the reference point are discarded. At the end a line is left with
c a distance between the points of minimum 100 m.
c
c Adaption of resolution:
c
c If the points of the line have no depth values, then the 
c smoothing and reduction is working as described above. If, however,
c depth values do exists, then they are used to change the resolution
c of the line points. For a value of 2, twice the points in the reduction
c algorithm are retained. In other words, the value for reduct is changed
c locally to a value of 100/2=50, and only points closer as 50 are
c eliminated. Line points without depth receive an implicit value of 1.
c The value -1 indicates that this point should never be eliminated 
c nor should it be smoothed.
c
c Please note that the change in resolution is applicable only for the
c reduction algorithm, and that the smoothing algorithm actually does
c only use -1 for not smoothing.
c
c***********************************************************************

	program gridr

c smoothing program
c
c same as gridf but smooth on length of line, not points

	implicit none

	integer nkndim,neldim,nlidim,nlndim,ndim
	parameter( nkndim = 50000 )
	parameter( neldim = 10000 )
	parameter( nlidim = 1000 )
	parameter( nlndim = 50000 )
	parameter( ndim = 20000 )

	integer ner
	integer nco,nkn,nel,nli
	logical bstop
	character*80 file

	integer ipv(nkndim)
	integer ipev(neldim)
	integer iaev(neldim)
	integer ianv(nkndim)
	integer nen3v(3,neldim)

	real xgv(nkndim)
	real ygv(nkndim)
	real hkv(nkndim)
	real hev(neldim)
	real hlv(nlidim)

	integer iplv(nlidim)
	integer ialrv(nlidim)
	integer ipntlv(0:nlidim)
	integer inodlv(nlndim)

	logical bperiod
	integer nl,nt
	real xt(ndim)
	real yt(ndim)
	real ht(ndim)

	integer l
	integer nline
	integer nnode
	real sigma
	real reduct

c-----------------------------------
c sigma		smoothing parameter (size of gaussian kernel)
c reduct	reduction of points in line (in meters)
c-----------------------------------
	sigma = 2.0
	sigma = 1.5
	sigma = 500.
	sigma = 200.
	reduct = 4.0
	reduct = 200.
c----------------------------------- hakata
	sigma = 200.
	reduct = 200.
c----------------------------------- circle
	sigma = 0.0
	reduct = 0.1
c----------------------------------- lido
	sigma = 0.
	reduct = 100.
c----------------------------------- curonian lagoon
	sigma = 200.
	reduct = 200.
c----------------------------------- brasile
	sigma = 0.005
	reduct = 0.005
c----------------------------------- malta
	sigma = 0.1
	reduct = 0.2
	sigma = 0.1 * 245.
	reduct = 0.2 * 245.
c----------------------------------- hue
	sigma = 30.
	reduct = 300.
c-----------------------------------

c------------------------------------------------------
c read in parameters and file name
c------------------------------------------------------

	file = ' '
	write(6,*) 'Enter file name with line: '
	read(5,'(a)') file
	if( file .eq. ' ' ) stop 'no file given'

	write(6,*) 'Enter sigma and reduct: '
	read(5,*) sigma,reduct

	write(6,*) 'file   :',file
	write(6,*) 'sigma  :',sigma
	write(6,*) 'reduct :',reduct

c------------------------------------------------------

	nline = 0
	nnode = 0

	bstop = .false.

        call rdgrd(
     +                   file
     +                  ,bstop
     +                  ,nco,nkn,nel,nli
     +                  ,nkndim,neldim,nlidim,nlndim
     +                  ,ipv,ipev,iplv
     +                  ,ianv,iaev,ialrv
     +                  ,hkv,hev,hlv
     +                  ,xgv,ygv
     +                  ,nen3v
     +                  ,ipntlv,inodlv
     +                  )

        if( bstop ) stop 'error stop rdgrd'

	write(6,*) 'comments : ',nco
	write(6,*) 'nodes    : ',nkn
	write(6,*) 'elements : ',nel
	write(6,*) 'lines    : ',nli

	call insnod(nkn,ipv)
c	call prline(nli,iplv,ialrv,ipntlv,inodlv)
c	call prlixy(nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv)

	open(99,file='smooth.grd',status='unknown',form='formatted')
	open(98,file='reduce.grd',status='unknown',form='formatted')

	do l=1,nli
	  nl = ndim
	  call extrli(l,nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv,hkv
     +				,xt,yt,ht,nl,nt)
	  call mkperiod(xt,yt,nl,bperiod)
	  call intpdep(ht,nl,bperiod)
	  nline = iplv(l)
	  call smooth(sigma,xt,yt,ht,nl,bperiod)
	  call wrline(99,nline,nnode,nl,xt,yt,ht,nt,bperiod)
	  call reduce(reduct,xt,yt,ht,nl)
	  call wrline(98,nline,nnode,nl,xt,yt,ht,nt,bperiod)
	end do

	write(6,*) 'routine finished...'

	close(99)
	close(98)

	write(6,*) 'files smooth.grd and reduce.grd written'

	end

c**********************************************************

	subroutine extrli(l,nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv,hkv
     +		,xt,yt,ht,nl,nt)

	implicit none

	integer l
	integer nli
	integer iplv(1)
	integer ialrv(1)
	integer ipntlv(0:1)
	integer inodlv(1)
	real xgv(1)
	real ygv(1)
	real hkv(1)
	integer nl
	real xt(1)
	real yt(1)
	real ht(1)
	integer nt

	integer nvert,i,ibase
	integer node,ier,k

	integer retriv

	if( l .gt. nli ) stop 'error stop extrli : l > nli'

	nvert = ipntlv(l) - ipntlv(l-1)
	ibase = ipntlv(l-1)

	if( nvert .gt. nl ) goto 98

	write(6,*) 'extracting line ',l,iplv(l),nvert,nt

	do i=1,nvert
	    node = inodlv(ibase+i)
	    ier = retriv(node,k)
	    if( ier .lt. 0 ) goto 99
	    xt(i) = xgv(k)
	    yt(i) = ygv(k)
	    ht(i) = hkv(k)
	end do

	nl = nvert
	nt = ialrv(l)

	return
   98	continue
	write(6,*) 'nvert = ',nvert
	write(6,*) 'nl = ',nl
	stop 'error stop extrli : dimension nl'
   99	continue
	stop 'error stop extrli: hash routines'
	end

c**********************************************************

	subroutine mkperiod(xt,yt,nl,bperiod)

c decides if line is periodic or not

	implicit none

	real xt(1)
	real yt(1)
	integer nl
	logical bperiod

	bperiod = .false.

	if( xt(1) .eq. xt(nl) .and. yt(1) .eq. yt(nl) ) then
	  bperiod = .true.
	  nl = nl - 1
	end if

	end

c**********************************************************

	subroutine prlixy(nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv)

	implicit none

	integer nli
	integer iplv(1)
	integer ialrv(1)
	integer ipntlv(0:1)
	integer inodlv(1)
	real xgv(1)
	real ygv(1)

	integer l,nvert,i,ibase
	integer node,ier,k

	integer retriv

	do l=1,nli
	  nvert = ipntlv(l) - ipntlv(l-1)
	  ibase = ipntlv(l-1)
	  write(6,*) 'line : ',l,iplv(l),ialrv(l),nvert
	  do i=1,nvert
	    node = inodlv(ibase+i)
	    ier = retriv(node,k)
	    if( ier .lt. 0 ) goto 99
	    write(6,*) i,node,xgv(k),ygv(k)
	  end do
	end do

	return
   99	continue
	stop 'error stop prlixy: hash routines'
	end

c**********************************************************

	subroutine prline(nli,iplv,ialrv,ipntlv,inodlv)

	implicit none

	integer nli
	integer iplv(1)
	integer ialrv(1)
	integer ipntlv(0:1)
	integer inodlv(1)

	integer l,nvert,i,ibase

	do l=1,nli
	  nvert = ipntlv(l) - ipntlv(l-1)
	  ibase = ipntlv(l-1)
	  write(6,*) 'line : ',l,iplv(l),ialrv(l),nvert
	  write(6,*) (inodlv(ibase+i),i=1,nvert)
	end do

	end

c**********************************************************

	subroutine insnod(nkn,ipv)

c inserts nodes into hash table

	implicit none

	integer nkn
	integer ipv(1)

	integer k,ier
	integer insert

	call init

	do k=1,nkn
	  ier = insert(ipv(k),k)
	  if( ier .lt. 0 ) goto 99
	end do

	return
   99	continue
	stop 'error stop insnod: hash routines'
	end

c********************************************************

	subroutine distxy(ndim,nl,xt,yt,ht,dxy)

c sets up distance between nodes

c dxy(i) -> distance to next node (i,i+1)
c dxy(i-1) -> distance to previous node (i,i-1)

	implicit none

	integer ndim
	integer nl
	real xt(1)
	real yt(1)
	real ht(1)
	real dxy(-ndim:2*ndim)

	integer i
	real xo,yo,xn,yn
	real dx,dy,dist
	real ho,hn,h

	xn = xt(nl)
	yn = yt(nl)
	hn = ht(nl)

	do i=1,nl
	  xo = xn
	  yo = yn
	  ho = hn
	  xn = xt(i)
	  yn = yt(i)
	  hn = ht(i)

	  dx = xn - xo
	  dy = yn - yo
	  dist = sqrt( dx*dx + dy*dy )

	  dxy(i-1) = dist
	  dxy(i-1+nl) = dist
	  dxy(i-1-nl) = dist
	end do

	dxy(2*nl) = dxy(0)

	end
	  
c********************************************************

	subroutine smooth(sigma,xt,yt,ht,nl,bperiod)

c smoothing

	real sigma
	real xt(1)
	real yt(1)
	real ht(1)
	integer nl
	logical bperiod

	integer ndim
	parameter(ndim=20000)

	integer i
	integer ngk
	real dist,dx,dy,distot
	real gk(-ndim:ndim)
	real raux(-ndim:2*ndim)
	real dxy(-ndim:2*ndim)

	if( nl .gt. ndim ) stop 'error stop smooth: ndim'

        if( sigma .le. 0. ) return

c set up dxy

	call distxy(ndim,nl,xt,yt,ht,dxy)

	distot = 0.
	do i=1,nl-1
	  distot = distot + dxy(i)
	end do

	dx = xt(1) - xt(nl)
	dy = yt(1) - yt(nl)
	dist = sqrt ( dx**2 + dy**2 )
	write(6,*) 'distance : ',nl,dist,distot,dist/distot,bperiod

c create smoothing kernel

c	ngk = ndim
c	call gkernel(sigma,ndim,ngk,gk)
c	write(6,*) 'Gaussian kernel : ',sigma,ngk

c smooth x and y

	call grsmooth(ndim,nl,sigma,xt,raux,dxy,ht,bperiod)
	call grsmooth(ndim,nl,sigma,yt,raux,dxy,ht,bperiod)

	end

c********************************************************

	subroutine gkernel(sigma,ndim,ngk,gk)

c gaussian smoothing

	real sigma
	integer ndim
	integer ngk
	real gk(-ndim:ndim)

	integer i
	real pi,eps,a,b,x,val

        pi = 4.*atan(1.)
        eps = 1.e-7

        a = 1. / ( sqrt(2.*pi) * sigma )
        b = 1. / ( 2 * sigma * sigma )

	do i=0,ndim
          x = i
          val = a * exp( -x**2 * b )
	  gk(i) = val
	  gk(-i) = val
	  if( val .lt. eps ) goto 1
        end do
    1	continue
	ngk = i-1

c	do i=-ngk,ngk
c	  write(6,*) i,gk(i)
c	end do

	end

c********************************************************

	subroutine grsmooth(ndim,nl,sigma,rt,raux,dxy,ht,bperiod)

c gaussian smoothing

	integer ndim
	integer nl
	real sigma
	real rt(1)
	real raux(-ndim:2*ndim)
	real dxy(-ndim:2*ndim)
	real ht(1)
	logical bperiod

	integer i,j
	real val
	real rkern,tkern
	real rintv

c set up circular array

	do i=1,nl
	  val = rt(i)
	  raux(i) = val
	  raux(i+nl) = val
	  raux(i-nl) = val
	end do

c set up gaussian kernel parameters

        pi = 4.*atan(1.)
        eps = 1.e-7

        a = 1. / ( sqrt(2.*pi) * sigma )
        b = - 1. / ( 2 * sigma * sigma )

c convolution

	do i=1,nl

	  rintv = 0.5 * ( dxy(i-1) + dxy(i) )
	  rkern = a
	  weight = rkern * rintv
	  tkern = weight
	  val = raux(i) * weight
	  jmin = i - nl
	  jmax = i + nl
	  if( .not. bperiod ) then
	    jmin = max(2,jmin)
	    jmax = min(nl-1,jmax)
	  end if
c	  write(6,*) i,i,weight,rkern,tkern,val

	  rkern = 1.
	  dist = 0.
	  j = i
	  do while( rkern .gt. eps .and. j .lt. jmax )
	    j = j + 1
	    dist = dist + dxy(j-1)
	    rintv = 0.5 * ( dxy(j-1) + dxy(j) )
	    rkern = a * exp( b * dist * dist )
	    weight = rkern * rintv
	    tkern = tkern + weight
	    val = val + raux(j) * weight
c	    write(6,*) i,j,weight,rkern,tkern,val,dist
	  end do

c	  write(6,*) i,j-i,nl,dist,rkern

	  rkern = 1.
	  dist = 0.
	  j = i
	  do while( rkern .gt. eps .and. j .gt. jmin )
	    j = j - 1
	    dist = dist + dxy(j)
	    rintv = 0.5 * ( dxy(j-1) + dxy(j) )
	    rkern = a * exp( b * dist * dist )
	    weight = rkern * rintv
	    tkern = tkern + weight
	    val = val + raux(j) * weight
c	    write(6,*) i,j,weight,rkern,tkern,val,dist
	  end do

c	  write(6,*) i,i-j,nl,dist,rkern

c	  write(6,*) '*** ',i,val,rt(i),tkern

	  if( tkern .gt. 0. .and. ht(i) .ge. 0. ) then
	    rt(i) = val / tkern
	  else
	    rt(i) = raux(i)
	  end if

	end do

	end

c********************************************************

	subroutine gsmooth(ndim,nl,rt,raux,ngk,gk)

c gaussian smoothing

	integer ndim
	integer nl
	real rt(1)
	real raux(-ndim:2*ndim)
	integer ngk
	real gk(-ndim:ndim)

	integer i,j
	real val

c set up circular array

	do i=1,nl
	  val = rt(i)
	  raux(i) = val
	  raux(i+nl) = val
	  raux(i-nl) = val
	end do

c convolution

	do i=1,nl
	  val = 0.
	  do j=-ngk,ngk
	    val = val + raux(i+j) * gk(j)
	  end do
	  rt(i) = val
	end do

	end

c********************************************************

	subroutine wrline(iunit,nline,nnode,nl,xt,yt,ht,nt,bperiod)

c write line

	integer iunit
	integer nline
	integer nnode
	integer nl
	real xt(1)
	real yt(1)
	real ht(1)
	integer nt
	logical bperiod

	integer i
	integer ntot

	if( nl .le. 0 ) return

	do i=1,nl
	  write(iunit,'(i1,2i10,3e14.6)') 1,nnode+i,nt,xt(i),yt(i),ht(i)
	end do
	  
c	nline = nline + 1

	ntot = nl
	if( bperiod ) ntot = ntot + 1

	write(iunit,'(i1,3i10)') 3,nline,nt,ntot
	write(iunit,*) (nnode+i,i=1,nl)
	if( bperiod ) then
	  write(iunit,*) nnode+1
	end if

	nnode = nnode + nl

	end

c********************************************************

	  subroutine reduce(reduct,xt,yt,ht,nl)

c reduces points in line

	implicit none

	real reduct
	real xt(1)
	real yt(1)
	real ht(1)
	integer nl

	integer i,nnew
	real rr,rtot
	real dist,dx,dy
	real h

c very simplicistic approach

	rr = reduct
	rtot = 1.

	if( nl .lt. 3 * rr ) then
	  rr = nl / 3.
	end if

	nnew = 0
	rtot = 0.
	if( ht(1) .lt. 0. .or. ht(1) .ge. 0. ) then	!bug fix 23.09.2004
	    nnew = nnew + 1
	    xt(nnew) = xt(1)
	    yt(nnew) = yt(1)
	    ht(nnew) = ht(1)
	end if
	do i=2,nl
	  dx = xt(i) - xt(i-1)
	  dy = yt(i) - yt(i-1)
	  h = ht(i)
	  dist = sqrt( dx*dx + dy*dy )
	  !write(6,*) dist,h,dist/h,rtot,reduct
	  if( h .gt. 0. ) dist = dist / h
	  rtot = rtot + dist
	  if( rtot .gt. reduct .or. h .lt. 0. ) then
	    nnew = nnew + 1
	    xt(nnew) = xt(i)
	    yt(nnew) = yt(i)
	    ht(nnew) = ht(i)
	    rtot = 0.
	  end if
	end do

c	if( nnew .lt. 3 ) stop 'error stop reduce: line too short...'
	if( nnew .lt. 3 ) then
	  write(6,*) 'line too short -> eliminated'
	  nnew = 0
	end if

	write(6,*) 'reduce: new nodes = ',nnew
	nl = nnew

	end

c********************************************************

	  subroutine intpdep(ht,nl,bperiod)

c interpolates depth values for points in line
c
c negative values are not interpolated and are left alone

	implicit none

	real ht(1)
	integer nl
	logical bperiod

	integer i
	integer ifirst,ilast,inext,nonzero
	integer nval
	real value,vfirst,vlast
	real aux
	logical bdebug

	bdebug = .false.

c------------------------------------------------------------
c look for values greater than 0
c------------------------------------------------------------

	nonzero = 0
	ifirst = 0
	ilast = 0

	do i=1,nl
	  if( ht(i) .gt. 0. ) then
	    nonzero = nonzero + 1
	    if( nonzero .eq. 1 ) ifirst = i
	    ilast = i
	  end if
	end do

c------------------------------------------------------------
c if no value or 1 value found -> set everything constant
c------------------------------------------------------------

	if( nonzero .le. 1 ) then
	  if( nonzero .eq. 0 ) then
	    value = 1.
	  else
	    value = ht(ifirst)
	  end if
	  do i=1,nl
	    if( ht(i) .ge. 0. .or. ht(i) .lt. -990. ) ht(i) = value
	  end do
	  return
	end if

c------------------------------------------------------------
c interpolate extreme values
c------------------------------------------------------------

	vfirst = ht(ifirst)
	vlast = ht(ilast)
	nval = nl - ilast + ifirst
	aux = (vfirst-vlast)/nval
	inext = 0

	if( bdebug ) write(6,*) bperiod
	if( bdebug ) write(6,*) ifirst,ilast,nval,vfirst,vlast

	i = 1
	do while( ht(i) .lt. 0 )
	  i = i + 1
	end do
	if( bperiod ) then
	  inext = nl - ilast + i
	  value = vlast + inext * aux
	else
	  value = ht(ifirst)
	end if
	if( bdebug ) write(6,*) i,inext,value,ht(i)
	ht(i) = value
	ifirst = i

	i = nl
	do while( ht(i) .lt. 0 )
	  i = i - 1
	end do
	if( bperiod ) then
	  inext = i - ilast
	  value = vlast + inext * aux
	else
	  value = ht(ilast)
	end if
	if( bdebug ) write(6,*) i,inext,value,ht(i)
	ht(i) = value
	ilast = i

c------------------------------------------------------------
c interpolate central values
c------------------------------------------------------------

	do while( ifirst .lt. ilast )
	  do i=ifirst+1,nl
	    if( ht(i) .gt. 0. ) goto 1
	  end do
	  stop 'error stop intpdep: impossible branch (1)'
    1	  continue
	  inext = i
	  nval = inext - ifirst
	  if( bdebug ) write(6,*) ifirst,inext,ilast,ht(ifirst),ht(inext)
	  do i=ifirst+1,inext-1
	    value = ht(ifirst) + (i-ifirst) * (ht(inext)-ht(ifirst))/nval
	    if( ht(i) .ge. 0. ) ht(i) = value
	    if( bdebug ) write(6,*) i,value
	  end do
	  ifirst = inext
	end do

c------------------------------------------------------------
c do some sanity check
c------------------------------------------------------------

	do i=1,nl
	  if( ht(i) .eq. 0 ) then
	    write(6,*) i,nl,ht(i)
	    stop 'error stop intpdep: values not interpolated'
	  end if
	  if( bdebug ) write(6,*) i,ht(i)
	end do

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	end

c********************************************************

