c
c $Id: basinf.f,v 1.20 2010-03-22 15:29:31 georg Exp $
c
c revision log :
c
c 06.04.1999	ggu	completely restructured
c 04.06.1999	ggu	new statistics are computed
c 28.11.2005	ggu	new call to makehkv
c 31.05.2007	ggu	added area and volume frequency curve
c 24.08.2007	ggu	added new routine write_grd_from_bas
c 06.04.2009    ggu     read param.h
c 12.06.2009    ggu     areatr in double precision - new algorithm
c 01.03.2010    ggu     new routine basqual() to compute grid quality
c 22.03.2010    ggu     write external element number in basqual()
c
c****************************************************************

        program basinf

c writes information on basin about nodes and elements

	implicit none

	include 'param.h'

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

	real haux(nkndim)
	real hkv(nkndim)
	common /hkv/hkv
	real hev(neldim)
	common /hev/hev
	include 'evmain.h'

	logical bnode,belem
	integer iapini

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        if( iapini(1,nkndim,neldim,0) .le. 0 ) stop

	call set_ev

c-----------------------------------------------------------------
c general info
c-----------------------------------------------------------------

        write(6,*)
        write(6,*) ' nkn = ',nkn,'  nel = ',nel
        write(6,*) ' mbw = ',mbw,'  ngr = ',ngr
        write(6,*)
        write(6,*) ' dcor = ',dcor,'  dirn = ',dirn
        write(6,*)

c-----------------------------------------------------------------
c specific info
c-----------------------------------------------------------------

	call basstat
	!call basqual		!grid quality
        call makehev(hev)
        call makehkv(hkv,haux)
	call freqdep

c-----------------------------------------------------------------
c loop until no node and element is given
c-----------------------------------------------------------------

	bnode = .true.
	belem = .true.

	do while( bnode .or. belem )

	   call nodeinfo(bnode)
	   if( .not. bnode .and. .not. belem ) goto 1

	   call eleminfo(belem)
	   if( .not. bnode .and. .not. belem ) goto 1

	end do

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

    1	continue

	!call write_grd_from_bas	!write grd from bas

	end

c*******************************************************************

	subroutine nodeinfo(bnode)

c info on node number

	implicit none

	logical bnode

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

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

	real hkv(nkndim)
	common /hkv/hkv

	integer ie,ii,in
	integer kext,kint
	integer ipext,ipint
	logical bloop

	bnode = .false.
	bloop = .true.

c look for node and give info

	do while( bloop )

        write(6,*) 'Enter node number : '
        read(5,'(i10)') kext
	if( kext .gt. 0 ) bnode = .true.
	kint = ipint(kext)

	if( kext .le. 0 ) then
	   bloop = .false.
        else if(kint.le.0) then
           write(6,*) ' no node number : ',kext
           if(kext.le.nkn) then
              write(6,*) '(intern : ',kext
     +                          ,' extern : ',ipext(kext),')'
           end if
           write(6,*)
	else
           write(6,*) ' extern : ',kext,' intern : ',kint
           write(6,*) '(intern : ',kext,' extern : ',ipext(kext),')'
           write(6,*) ' (x,y)  : ',xgv(kint),ygv(kint)
           write(6,*) ' depth  : ',hkv(kint)
           write(6,2200)

           do ie=1,nel
              in=0
              do ii=1,3
                 if(nen3v(ii,ie).eq.kint) in=ii
              end do
              if(in.gt.0) then
                 write(6,2000)   ipev(ie)
     +                          ,(ipext(nen3v(ii,ie)),ii=1,3)
     +                          ,(hm3v(ii,ie),ii=1,3)
     +                          ,iarv(ie),hm3v(in,ie)
              end if
	   end do
           write(6,*)
	end if

	end do

	return
 2200   format(/1x,'element',8x,'nodes',14x,'depth in element'
     +          ,3x,'  area',4x,'depth of node')
 2000   format(1x,i6,2x,3i6,2x,3f8.2,2x,i5,5x,f8.2)
	end

c*****************************************************************

	subroutine eleminfo(belem)

c info on element number

	implicit none

	logical belem

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

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

	integer ie,ii,k
	integer eext,eint
	integer ipext,ieext,ieint
	logical bloop
	real areatr

	belem = .false.
	bloop = .true.

	do while( bloop )

        write(6,*) 'Enter element number : '
        read(5,'(i10)') eext
        if(eext.gt.0) belem = .true.
	eint = ieint(eext)

	if( eext .le. 0 ) then
	  bloop = .false.
        else if(eint.le.0) then
           write(6,*) ' no element number : ',eext
           if(eext.le.nel) then
              write(6,*) '(intern : ',eext
     +                     ,' extern : ',ieext(eext),')'
           end if
           write(6,*)
	else
           write(6,*) ' extern : ',eext,' intern : ',eint
           write(6,*) '(intern : ',eext,' extern : ',ieext(eext),')'
           write(6,*)

	   ie = eint
           do ii=1,3
              k=nen3v(ii,ie)
              write(6,*) ' (x,y) : ',xgv(k),ygv(k)
	   end do

           write(6,3200)
           write(6,3000) ieext(ie)
     +                  ,(ipext(nen3v(ii,ie)),ii=1,3)
     +                  ,(hm3v(ii,ie),ii=1,3)
     +                  ,areatr(ie)
     +                  ,iarv(ie)
           write(6,*)
	end if

	end do

        return
 3200   format(/1x,'element',8x,'nodes',14x,'depth in element'
     +          ,3x,'   area',3x,'  area code')
 3000           format(1x,i6,2x,3i6,2x,3f8.2,2x,e10.2,2x,i5)
        end

c*****************************************************************

	subroutine basstat

c writes statistics on basin

	implicit none

	include 'param.h'

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

	integer ie,ii,k
	integer imin,imax
	real area,amin,amax,atot
        real vtot
	real aptot,vptot
	real x(3),y(3)
	real xmin,xmax,ymin,ymax
	real dxmax,dymax
	real h
        real xx,yy
        real dist,distmin
        integer i,k1,k2

	real areatr

c-----------------------------------------------------------------
c area code
c-----------------------------------------------------------------

	imin = iarv(1)
	imax = imin

	do ie=1,nel
	  imin = min(imin,iarv(ie))
	  imax = max(imax,iarv(ie))
	end do

	write(6,*) 'Area code min/max:      ',imin,imax

c-----------------------------------------------------------------
c node numbers
c-----------------------------------------------------------------

	imin = ipv(1)
	imax = imin

	do k=1,nkn
	  imin = min(imin,ipv(k))
	  imax = max(imax,ipv(k))
	end do

	write(6,*) 'Node number min/max:    ',imin,imax

c-----------------------------------------------------------------
c element numbers
c-----------------------------------------------------------------

	imin = ipev(1)
	imax = imin

	do ie=1,nel
	  imin = min(imin,ipev(ie))
	  imax = max(imax,ipev(ie))
	end do

	write(6,*) 'Element number min/max: ',imin,imax

c-----------------------------------------------------------------
c area
c-----------------------------------------------------------------

	amin = areatr(1)
	amax = amin
	atot = 0.
        vtot = 0.
	aptot = 0.
        vptot = 0.

	do ie=1,nel
	  area = areatr(ie)
	  atot = atot + area
	  amin = min(amin,area)
	  amax = max(amax,area)
          h = 0.
          do ii=1,3
            h = h + hm3v(ii,ie)
          end do
          vtot = vtot + area * h / 3.
	  if( h .gt. 0. ) then		!only positive depths
	    aptot = aptot + area
            vptot = vptot + area * h / 3.
	  end if
	end do

	write(6,*) 'Area min/max:           ',amin,amax
	write(6,*) 'Total area:             ',atot,aptot
	write(6,*) 'Total volume:           ',vtot,vptot

c-----------------------------------------------------------------
c coordinates
c-----------------------------------------------------------------

	xmin = xgv(1)
	ymin = ygv(1)
	xmax = xgv(1)
	ymax = ygv(1)

	do k=1,nkn
	  xmin = min(xmin,xgv(k))
	  ymin = min(ymin,ygv(k))
	  xmax = max(xmax,xgv(k))
	  ymax = max(ymax,ygv(k))
	end do

	write(6,*) 'X-Coordinates min/max:  ',xmin,xmax
	write(6,*) 'Y-Coordinates min/max:  ',ymin,ymax

c-----------------------------------------------------------------
c size of elements
c-----------------------------------------------------------------

	dxmax = 0.
	dymax = 0.

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    x(ii) = xgv(k)
	    y(ii) = ygv(k)
	  end do
	  xmin = min(x(1),x(2),x(3))
	  xmax = max(x(1),x(2),x(3))
	  ymin = min(y(1),y(2),y(3))
	  ymax = max(y(1),y(2),y(3))
	  dxmax = max(dxmax,xmax-xmin)
	  dymax = max(dymax,ymax-ymin)
	end do

	write(6,*) 'Element dxmax/dymax:    ',dxmax,dymax

c-----------------------------------------------------------------
c depth
c-----------------------------------------------------------------

	amin = 999999.
	amax = -amin

	do ie=1,nel
	  h = 0
	  do ii=1,3
	    h = h + hm3v(ii,ie)
	  end do
	  h = h / 3.
	  amin = min(amin,h)
	  amax = max(amax,h)
	end do

	write(6,*) 'Depth min/max:          ',amin,amax
	write(6,*) 'Depth average:          ',vtot/atot,vptot/aptot

c-----------------------------------------------------------------
c minimum distance of nodes
c-----------------------------------------------------------------

        distmin = (xmax-xmin)**2 + (ymax-ymin)**2
        distmin = 2*distmin
        k1 = 0
        k2 = 0

        do k=1,nkn
          xx = xgv(k)
          yy = ygv(k)
          do i=k+1,nkn
            dist = (xx-xgv(i))**2 + (yy-ygv(i))**2
            if( dist .lt. distmin ) then
              k1 = k
              k2 = i
              distmin = dist
            end if
          end do
        end do

        distmin = sqrt(distmin)
	write(6,*) 'min node distance:      ',distmin,ipv(k1),ipv(k2)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	write(6,*)

	end

c*******************************************************************

	subroutine freqdep

c writes frequency distribution of depth

	implicit none

	include 'param.h'

	integer ndim
	parameter (ndim=10000)

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

	integer ie,i
	integer imax,ih
	real area,vol
	real hmin,hmax,dh,h,fr
	real fa,fv

	double precision freqa(0:ndim)
	double precision freqv(0:ndim)

	real areatr

c-----------------------------------------------------------------
c area code
c-----------------------------------------------------------------

	hmin = -1.
	hmax = 40.
	dh = 0.10
	dh = 0.01
	imax = (hmax-hmin)/dh
	if( imax .gt. ndim ) then
	  write(6,*) imax,ndim,hmin,hmax,dh
	  stop 'error stop freqdep: ndim'
	end if

	write(6,*) 'computing depth frequency curve ',imax

	do i=0,imax
	  freqa(i) = 0.
	  freqv(i) = 0.
	end do

	do ie=1,nel
	  h = hev(ie)
	  area = areatr(ie)
	  vol = area * h
	  ih = (h - hmin) / dh
	  !write(6,*) ie,area,h,ih
	  if( ih .lt. 0 .or. ih .gt. imax ) then
	    if( ih .lt. 0 ) then
		ih = 0
	    else if( ih .gt. imax ) then
		ih = imax
	    else		!never
	        write(6,*) ih,imax,hmin,hmax,h
	        stop 'error stop freqdep: index out of range'
	    end if
	  end if
	  do i=ih,imax
	    freqa(i) = freqa(i) + area
	    freqv(i) = freqv(i) + vol
	  end do
	end do

	area = freqa(imax)	!total area
	vol  = freqv(imax)	!total volume
	write(6,*) 'total area/vol for freq curve: ',area,vol

	do i=0,imax
	  freqa(i) = freqa(i) / area
	  freqv(i) = freqv(i) / vol
	end do

	do i=0,imax
	  h = hmin + i*dh
	  fa = freqa(i) * 100.
	  fv = freqv(i) * 100.
	  !write(6,*) i,h,fa
	  write(66,*) h,fa
	  write(67,*) h,fv
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

	integer ii,i1,i2,k1,k2
	double precision f,x(3),y(3)

        do ii=1,3
          k=nen3v(ii,ie)
	  x(ii) = xgv(k)
	  y(ii) = ygv(k)
        end do

	f = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))

        areatr = f / 2.D0

        end

c*******************************************************************

	subroutine write_grd_from_bas

c writes grd file extracting info from bas file

	implicit none

	include 'param.h'

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

	integer k,ie,ii,ia
	real x,y,h

	open(8,file='bas.grd',status='unknown',form='formatted')

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  !write(8,2000) 1,k,0,x,y
	  write(8,2000) 1,ipv(k),0,x,y
	end do

	do ie=1,nel
	  h = hm3v(1,ie)
	  ia = iarv(ie)
	  !write(8,2100) 2,ie,0,3,(nen3v(ii,ie),ii=1,3),h
	  write(8,2100) 2,ipev(ie),ia,3,(ipv(nen3v(ii,ie)),ii=1,3),h
	end do
	  
	close(8)

	return
 2000	format(i1,i10,i5,2e14.6)
 2100	format(i1,i10,2i5,3i10,e14.6)
	end

c*******************************************************************

	subroutine basqual

c writes statistics on grid quality

	implicit none

	include 'param.h'

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
	include 'evmain.h'

	real areav(nkndim)

	integer ie,ii,k
	integer imin,imax
	real area,amin,amax,atot
        real vtot
	real aptot,vptot
	real x(3),y(3)
	real xmin,xmax,ymin,ymax
	real dxmax,dymax
	real h
        real xx,yy
        real dist,distmin
	real fmax,fmin,a,f
        integer i,k1,k2
	integer iemin,iemax

	integer ieext
	real areatr

c-----------------------------------------------------------------
c initialization
c-----------------------------------------------------------------

	do k=1,nkn
	  areav(k) = 0.
	end do

c-----------------------------------------------------------------
c compute fraction
c-----------------------------------------------------------------

	do ie=1,nel
	  area = 4.*ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    areav(k) = areav(k) + area
	  end do
	end do

	fmax = 0
	fmin = 10.
	do ie=1,nel
	  area = 12.*ev(10,ie)
	  amax = 0.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    a = areav(k)
	    amax = max(amax,a)
	  end do
	  f = amax/area
	  if( f .gt. 10. ) then
	    write(6,*) 'bad quality of element: ',ie,ieext(ie),f
	  end if
	  if( f .gt. fmax ) then
	    fmax = f
	    iemax = ie
	  end if
	  if( f .lt. fmin ) then
	    fmin = f
	    iemin = ie
	  end if
	end do

	write(6,*) 'Grid quality: (internal/external element number)'
	write(6,*) '   minimum: ',iemin,ieext(iemin),fmin
	write(6,*) '   maximum: ',iemax,ieext(iemax),fmax
	write(6,*)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************
