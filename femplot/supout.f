c
c $Id: supout.f,v 1.15 2010-02-26 15:29:19 georg Exp $
c
c routines for reading data files
c
c revision log :
c
c 31.10.2003  ggu     new routines velclose(), resetsim()
c 31.10.2003  ggu     new routines for handling wind
c 22.09.2004  ggu     new routines for handling ous file
c 22.09.2004  ggu     use OUS file for 3D, try to use alsways 3D (level=0)
c 05.10.2004  ggu     adjustments -> use always 3D data structure
c 14.03.2007  ggu     new routines for wave plotting
c 17.09.2008  ggu     new routine level_e2k to compute ilhkv from ilhv
c 09.10.2009  ggu     read also pressure from wind file
c 13.10.2009  ggu     set nlv once file is read
c 23.02.2010  ggu     change in reading wind file
c
c**********************************************************
c**********************************************************
c**********************************************************

        subroutine velopen

c opens velocity file (2d and 3d)

        implicit none

        logical is2d

        if( is2d() ) then
          call outopen
        else
          call ousopen
        end if

        end

c**********************************************************

        function velnext(it)

c gets next velocity field (2d and 3d)

        implicit none

        logical velnext
        integer it

        logical outnext, ousnext, is2d

        if( is2d() ) then
          velnext = outnext(it)
        else
          velnext = ousnext(it)
        end if

        end

c**********************************************************

        subroutine velclose

c closes velocity file

        implicit none

        logical is2d

        if( is2d() ) then
          call outclose
        else
          call ousclose
        end if

        end

c******************************************************
c******************************************************
c******************************************************

	function is2d()

	implicit none

	logical is2d
	integer getlev

c -1	use 2D model
c  0	use 3D model, barotropic currents
c >0	use 3D, level

	is2d = getlev() .lt. 0

	end

c******************************************************

        subroutine resetsim

c resets mask data structure

	implicit none

	logical bwater(1)
	common /bwater/bwater
	logical bkwater(1)
	common /bkwater/bkwater

        call initmask(bwater)
	call nodemask(bwater,bkwater)

        end

c******************************************************

	subroutine prepsim

c prepares simulation for use

	implicit none

	logical bwater(1)
	common /bwater/bwater
	logical bkwater(1)
	common /bkwater/bkwater
        integer ilhv(1)
        common /ilhv/ilhv
	real znv(1)
	common /znv/znv

	integer level
	real href,hzmin

	logical is2d
	integer getlev
	real getpar

c---------------------------------------------------
c copy data structures -> we use only 3D
c---------------------------------------------------

	if( is2d() ) then
	  !call from2to3
	else
	  !call from3to2
	end if

c---------------------------------------------------
c set up mask of water points
c---------------------------------------------------

        href = getpar('href')
        hzmin = getpar('hzmin')
	level = getlev()

        call initmask(bwater)			!true for all elements
        call drymask(bwater,znv,href,hzmin)	!false if znv/zenv not equal
        call levelmask(bwater,ilhv,level)	!element has this level

	call nodemask(bwater,bkwater)		!copy element to nodal mask

c---------------------------------------------------
c end of routine
c---------------------------------------------------

	end

c******************************************************

	subroutine from3to2

c sets up 2D data structures from 3D
c
c -> do not use anymore !!!!! (outdated) , may be deleted
c
c	nlv,nlvdi    		nlv=1 for 2d
c	hlv			hlv(1) = 10000 for 2d
c
c	hev			is set after basin read, newly set by 3d
c	ilhv			ilhv(ie) = 1 for 2d
c	znv,utlnv,vtlnv 	-> set zenv, usnv, vsnv

	implicit none

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real znv(1)
        common /znv/znv
        real zenv(3,1)
        common /zenv/zenv
        integer nen3v(3,1)
        common /nen3v/nen3v
        real xv(3,1)
        common /xv/xv
        integer ilhv(1)
        common /ilhv/ilhv
        real utlnv(nlvdim,1)
        common /utlnv/utlnv
        real vtlnv(nlvdim,1)
        common /vtlnv/vtlnv
        real usnv(1), vsnv(1)
        common /usnv/usnv, /vsnv/vsnv

	integer k,ie,ii,l,lmax
	real utot,vtot

c we do not use xv anymore

	do k=1,nkn
	  xv(1,k) = 0.
	  xv(2,k) = 0.
	  xv(3,k) = znv(k)
	end do

c we already have read zenv from file ous

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    zenv(ii,ie) = znv(k)
	  end do
	end do

c interpolation is now done in plo3vel

	do ie=1,nel
	  utot = 0.
	  vtot = 0.
	  lmax = ilhv(ie)
	  do l=1,lmax
	    utot = utot + utlnv(l,ie)
	    vtot = vtot + vtlnv(l,ie)
	  end do
	  usnv(ie) = utot
	  vsnv(ie) = vtot
	end do

	end

c******************************************************

	subroutine from2to3

c sets up 3D data structures from 2D
c
c outdated -> do not call - arrays are set in outnext()
c
c	xv,zenv,usnv,vsnv -> set znv

	implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real znv(1)
        common /znv/znv
        real xv(3,1)
        common /xv/xv
        integer ilhv(1)
        common /ilhv/ilhv

	integer k,ie

	do k=1,nkn
	  znv(k) = xv(3,k)
	end do

	do ie=1,nel
	  ilhv(ie) = 1
	end do

	end

c******************************************************
c******************************************************
c******************************************************

        subroutine windini

        implicit none

        integer nunit,iform
        common /winwin/ nunit,iform
	save /winwin/

	integer icall
	save icall
	data icall /0/

	if( icall .ne. 0 ) return

	icall = 1

	nunit = 0
        iform = 0

        end

c******************************************************

        subroutine windclose

        implicit none

        integer nunit,iform
        common /winwin/ nunit,iform

	if( nunit .gt. 0 ) close(nunit)

        end

c******************************************************

	subroutine windopen

	implicit none

        integer nunit,iform
        common /winwin/ nunit,iform

	integer ideffi

	call windini

	nunit = ideffi('datdir','runnam','.win','unform','old')
	if( nunit .le. 0 ) then
		stop 'error stop windopen: cannot open WIN file'
	end if

        call timeset(0,0,0)

        end

c******************************************************

	function windnext(it)

	implicit none

	logical windnext
	integer it

        integer nunit,iform
        common /winwin/ nunit,iform

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer ierr,n

        real uv(1), vv(1)
        real pres(1)
        common /uv/uv, /vv/vv
        common /pres/pres

        n = nkn
        call rdwin(nunit,it,n,uv,vv,pres,ierr)

	if( ierr .gt. 0 ) then
		stop 'error stop windnext: error reading data record'
	else if( ierr .lt. 0 ) then
		windnext = .false.
	else
		windnext = .true.
	end if

        end

c******************************************************
c******************************************************
c******************************************************

        subroutine waveini

        implicit none

        integer nunit,iform
        common /wavwav/ nunit,iform
	save /wavwav/

	integer icall
	save icall
	data icall /0/

	if( icall .ne. 0 ) return

	icall = 1

	nunit = 0
        iform = 0	! 0 = wave height   1 = wave period

        end

c******************************************************

        subroutine waveclose

        implicit none

        integer nunit,iform
        common /wavwav/ nunit,iform

	if( nunit .gt. 0 ) close(nunit)

        end

c******************************************************

	subroutine waveopen

	implicit none

	character*80 descrp
        common /descrp/ descrp
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv
        real hlv(1), hev(1)
        common /hlv/hlv, /hev/hev
        integer ilhkv(1)
        common /ilhkv/ilhkv

        integer nunit,iform
        common /wavwav/ nunit,iform

	integer nknaux,nelaux,nlvaux
	integer nvers,nvar

	integer ideffi

	call waveini

	nunit = ideffi('datdir','runnam','.wav','unform','old')
	if( nunit .le. 0 ) then
		stop 'error stop waveopen: cannot open WAV file'
	end if

	nvers = 3
	call rhnos(nunit,nvers
     +                          ,nkn,nel,nlv
     +                          ,nknaux,nelaux,nlvaux,nvar
     +                          ,ilhkv,hlv,hev
     +                          ,descrp
     +                          )

	if( nknaux .ne. nkn ) goto 99
	if( nelaux .ne. nel ) goto 99
	if( nlvaux .ne. nlv ) goto 99
	if( nvar .ne. 3 ) goto 99

        call timeset(0,0,0)

	return
   99	continue
	write(6,*) nkn,nel,nlv
	write(6,*) nknaux,nelaux,nlvaux
	write(6,*) nvar
	stop 'error stop waveopen: parameter mismatch'
        end

c******************************************************

	function wavenext(it)

	implicit none

	logical wavenext
	integer it

        integer nunit,iform
        common /wavwav/ nunit,iform

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer ilhkv(1)
        common /ilhkv/ilhkv

        integer ierr,nlvdim,ivar

        real uv(1), vv(1)
        real v1v(1), v2v(1), v3v(1)
        common /uv/uv, /vv/vv
        common /v1v/v1v
        common /v2v/v2v
        common /v3v/v3v

	wavenext = .false.
	nlvdim = 1
	call rdnos(nunit,it,ivar,nlvdim,ilhkv,v1v,ierr)
	if( ierr .gt. 0 ) goto 99
	if( ierr .lt. 0 ) return
	if( ivar .ne. 31 ) goto 97
	call rdnos(nunit,it,ivar,nlvdim,ilhkv,v2v,ierr)
	if( ierr .gt. 0 ) goto 99
	if( ierr .lt. 0 ) goto 98
	if( ivar .ne. 32 ) goto 97
	call rdnos(nunit,it,ivar,nlvdim,ilhkv,v3v,ierr)
	if( ierr .gt. 0 ) goto 99
	if( ierr .lt. 0 ) goto 98
	if( ivar .ne. 33 ) goto 97

	wavenext = .true.

	if( iform .eq. 0 ) then
	  call polar2xy(nkn,v1v,v3v,uv,vv)
	else
	  call polar2xy(nkn,v2v,v3v,uv,vv)
	end if

	return
   99	continue
	stop 'error stop wavenext: error reading data'
   98	continue
	stop 'error stop wavenext: not enough records'
   97	continue
	stop 'error stop wavenext: records not in order'
        end

c******************************************************

	subroutine polar2xy(n,speed,dir,uv,vv)

	implicit none

	integer n
	real speed(1), dir(1)
	real uv(1), vv(1)

	integer i
	real rad,a

	rad = atan(1.) / 45.

	do i=1,n
	  a = dir(i)
          a = 90. - a + 180.
          do while( a .lt. 0. )
            a = a + 360.
          end do
          a = mod(a,360.)

	  uv(i) = speed(i) * cos( rad * a )
	  vv(i) = speed(i) * sin( rad * a )
	end do

	end

c******************************************************
c******************************************************
c******************************************************

	subroutine outini

	implicit none

	integer nunit,nvers
	common /outout/ nunit,nvers
	save /outout/

	integer icall
	save icall
	data icall /0/

	if( icall .ne. 0 ) return

	icall = 1

	nunit = 0
	nvers = 0

	end

c******************************************************

	subroutine outclose

	implicit none

	integer nunit,nvers
	common /outout/ nunit,nvers

	if( nunit .gt. 0 ) close(nunit)

	end

c******************************************************

	subroutine outopen

	implicit none

	integer nunit,nvers
	common /outout/ nunit,nvers

	character*80 descrp
        common /descrp/ descrp
	integer itanf,itend,idt,idtout
	integer ierr
	real href,hzoff

	character*80 file

	integer rfout,ideffi

	call outini

	nunit = ideffi('datdir','runnam','.out','unform','old')
	if( nunit .le. 0 ) then
		stop 'error stop outopen: cannot open OUT file'
	end if

        write(6,*) 'file OUT opened'
        inquire(nunit,name=file)
        write(6,*) 'Reading file ...'
        write(6,*) file

	ierr=rfout(nunit,nvers,itanf,itend,idt,idtout,href,hzoff,descrp)
	if( ierr .ne. 0 ) then
		stop 'error stop outopen: error reading header'
	end if

        call putpar('href',href)
        call putpar('hzoff',hzoff)

        write(6,*)
        write(6,*)   descrp
        write(6,*)
        write(6,*) ' itanf =',itanf,'  itend =',itend
        write(6,*) '   idt =',idt,  ' idtout =',idtout
        write(6,*)
        write(6,*) ' nvers =',nvers
        write(6,*) '  href =',href, '  hzoff =',hzoff
        write(6,*)

	call timeset(itanf,itend,idtout)

	end

c******************************************************

	function outnext(it)

	implicit none

	logical outnext
	integer it

	integer nunit,nvers
	common /outout/ nunit,nvers

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xv(3,1)
        real zenv(3,1)
        real usnv(1), vsnv(1)
        common /xv/xv
        common /zenv/zenv
        common /usnv/usnv, /vsnv/vsnv

        real znv(1)
        common /znv/znv
        integer ilhv(1)
        common /ilhv/ilhv
        integer ilhkv(1)
        common /ilhkv/ilhkv

	integer ierr,nk,ne
        integer k,ie

	integer rdout

	nk = nkn
	ne = nel

	ierr=rdout(nunit,nvers,it,nk,ne,xv,zenv,usnv,vsnv)

	if( ierr .gt. 0 ) then
		!stop 'error stop outnext: error reading data record'
		write(6,*) '*** outnext: error reading data record'
		outnext = .false.
	else if( nk .ne. nkn .or. ne .ne. nel ) then
		stop 'error stop outnext: internal error (1)'   !???
		nkn = nk
		nel = ne
	else if( ierr .lt. 0 ) then
		outnext = .false.
	else
		outnext = .true.
	end if

	do k=1,nkn
	  znv(k) = xv(3,k)
	end do

	do ie=1,nel
	  ilhv(ie) = 1
	end do

	do k=1,nkn
	  ilhkv(k) = 1
	end do

	end

c******************************************************
c******************************************************
c******************************************************

	subroutine ousini

c initializes internal data structure for OUS file

	implicit none

	integer nunit
	common /ousous/ nunit
	save /ousous/

	integer icall
	save icall
	data icall /0/

	if( icall .ne. 0 ) return

	icall = 1

	nunit = 0

	end

c******************************************************

	subroutine ousclose

c closes OUS file

	implicit none

	integer nunit
	common /ousous/ nunit

	if( nunit .gt. 0 ) close(nunit)
	nunit = 0

	end

c******************************************************

        subroutine ousinfo(nvers,nkn,nel,nlv)

c returns info on parameters

        implicit none

	integer nunit
	common /ousous/ nunit

        integer nvers,nkn,nel,nlv

        call getous(nunit,nvers,nkn,nel,nlv)

        end

c******************************************************

	subroutine ousopen

c opens OUS file and reads header

	implicit none

	integer nunit
	common /ousous/ nunit

	character*80 descrp
        common /descrp/ descrp
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        real hlv(1), hev(1)
        common /hlv/hlv, /hev/hev
        integer ilhv(1)
        common /ilhv/ilhv

	character*80 file

	integer nvers
	integer nknaux,nelaux,nlvaux,nvar
	integer ierr
        integer i
        real href,hzmin
	integer ideffi

c initialize routines

	call ousini

c open file

	nunit = ideffi('datdir','runnam','.ous','unform','old')
	if( nunit .le. 0 ) goto 96

        write(6,*) 'file OUS opened'
        inquire(nunit,name=file)
        write(6,*) 'Reading file ...'
        write(6,*) file

c read first header

	nvers = 1
	call rfous(nunit,nvers,nknaux,nelaux,nlvaux
     +                  ,href,hzmin,descrp,ierr)

	if( ierr .ne. 0 ) goto 98

        write(6,*)
        write(6,*)   descrp
        write(6,*)
        write(6,*) ' nvers = ',nvers
        write(6,*) '   nkn = ',nknaux, '   nel = ',nelaux
        write(6,*) '   nlv = ',nlvaux
        write(6,*) '  href = ',href
        write(6,*) ' hzmin = ',hzmin
        write(6,*)

	if( nkn .ne. nknaux ) goto 99
	if( nel .ne. nelaux ) goto 99
	if( nlvdi .lt. nlvaux ) goto 99

	nlv = nlvaux

c read second header

	call rsous(nunit,ilhv,hlv,hev,ierr)

	if( ierr .ne. 0 ) goto 97

	call level_e2k		!computes ilhkv

c initialize time

	call timeset(0,0,0)

c        write(66,*) nvers,nknaux,nelaux,nlvaux
c        write(66,*) href,hzmin
c        write(66,*) descrp
c        write(66,*) (hlv(i),i=1,nlvaux)
c        write(66,*) (ilhv(i),i=1,nelaux)
c        write(66,*) (hev(i),i=1,nelaux)

c end

	return
   96	continue
	stop 'error stop ousopen: cannot open OUS file'
   97	continue
	stop 'error stop ousopen: error reading second header'
   98	continue
	stop 'error stop ousopen: error reading first header'
   99	continue
	write(6,*) 'error in parameters :'
	write(6,*) 'nkn : ',nkn,nknaux
	write(6,*) 'nel : ',nel,nelaux
	write(6,*) 'nlv : ',nlvdi,nlvaux
	stop 'error stop ousopen'
	end

c******************************************************

	function ousnext(it)

c reads next OUS record - is true if a record has been read, false if EOF

	implicit none

	logical ousnext		!true if record read, flase if EOF
	integer it		!time of record

        include 'param.h'

	integer nunit
	common /ousous/ nunit
        integer nlvdi,nlv
        common /level/ nlvdi,nlv
        integer ilhv(1)
        common /ilhv/ilhv
        real znv(1)
        common /znv/znv
        real zenv(3,1)
        common /zenv/zenv
        real utlnv(nlvdim,1)
        common /utlnv/utlnv
        real vtlnv(nlvdim,1)
        common /vtlnv/vtlnv

	integer ierr

	if( nlvdim .ne. nlvdi ) stop 'error stop ousnext: nlvdim'

	call rdous(nunit,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)

c set return value

	if( ierr .gt. 0 ) then
		!stop 'error stop ousnext: error reading data record'
		write(6,*) '*** ousnext: error reading data record'
		ousnext = .false.
	else if( ierr .lt. 0 ) then
		ousnext = .false.
	else
		ousnext = .true.
	end if

c end

	end

c******************************************************
c******************************************************
c******************************************************
c 3d-version of concentration
c******************************************************
c******************************************************
c******************************************************

	subroutine nosini

c initializes internal data structure for NOS file

	implicit none

	integer nunit
	common /nosnos/ nunit
	save /nosnos/

	integer icall
	save icall
	data icall /0/

	if( icall .ne. 0 ) return

	icall = 1

	nunit = 0

	end

c******************************************************

	subroutine nosclose

c closes NOS file

	implicit none

	integer nunit
	common /nosnos/ nunit

	if( nunit .gt. 0 ) close(nunit)
	nunit = 0

	end

c******************************************************

	subroutine nosopen(type)

c opens NOS file and reads header

	implicit none

	character*(*) type

	integer nunit
	common /nosnos/ nunit

	character*80 descrp
        common /descrp/ descrp
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        real hlv(1), hev(1)
        common /hlv/hlv, /hev/hev
        integer ilhkv(1)
        common /ilhkv/ilhkv

	character*80 file

	integer nvers
	integer nknaux,nelaux,nlvaux,nvar
	integer ierr
	integer ideffi

c initialize routines

	call nosini

c open file

	nunit = ideffi('datdir','runnam',type,'unform','old')
	if( nunit .le. 0 ) then
		stop 'error stop conopen: cannot open NOS file'
        else
                write(6,*) 'File opened :'
                inquire(nunit,name=file)
                write(6,*) file
                write(6,*) 'Reading file ...'
	end if

c read first header

	nvers = 3
	call rfnos(nunit,nvers,nknaux,nelaux,nlvaux,nvar,descrp,ierr)

	if( ierr .ne. 0 ) then
		stop 'error stop nosopen: error reading first header'
	end if

        write(6,*)
        write(6,*)   descrp
        write(6,*)
        write(6,*) ' nvers = ', nvers
        write(6,*) '   nkn = ',nknaux, '   nel = ',nelaux
        write(6,*) '   nlv = ',nlvaux, '  nvar = ',  nvar
        write(6,*)

	if( nkn .ne. nknaux ) goto 99
	if( nelaux .ne. 0 .and. nel .ne. nelaux ) goto 99
	if( nlvdi .lt. nlvaux ) goto 99

	nlv = nlvaux

c read second header

	call rsnos(nunit,ilhkv,hlv,hev,ierr)

	if( ierr .ne. 0 ) then
		stop 'error stop nosopen: error reading second header'
	end if

c initialize time

	call timeset(0,0,0)

c end

	return
   99	continue
	write(6,*) 'error in parameters :'
	write(6,*) 'nkn : ',nkn,nknaux
	write(6,*) 'nel : ',nel,nelaux
	write(6,*) 'nlv : ',nlvdi,nlvaux
	stop 'error stop nosopen'
	end

c******************************************************

	function nosnext(it,ivar,nlvdim,array)

c reads next NOS record - is true if a record has been read, false if EOF

	implicit none

	logical nosnext		!true if record read, flase if EOF
	integer it		!time of record
	integer ivar		!type of variable
	integer nlvdim		!dimension of vertical coordinate
	real array(nlvdim,1)	!values for variable

	integer nunit
	common /nosnos/ nunit
        integer nlvdi,nlv
        common /level/ nlvdi,nlv
        integer ilhkv(1)
        common /ilhkv/ilhkv

	integer ierr

	if( nlvdim .ne. nlvdi ) stop 'error stop nosnext: nlvdim'

	call rdnos(nunit,it,ivar,nlvdim,ilhkv,array,ierr)

c set return value

	if( ierr .gt. 0 ) then
		!stop 'error stop nosnext: error reading data record'
		write(6,*) '*** nosnext: error reading data record'
		nosnext = .false.
	else if( ierr .lt. 0 ) then
		nosnext = .false.
	else
		nosnext = .true.
	end if

c end

	end

c******************************************************

	subroutine level_e2k

c computes max level at nodes from elements

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer nen3v(3,1)
	common /nen3v/nen3v
        integer ilhv(1)
        common /ilhv/ilhv
        integer ilhkv(1)
        common /ilhkv/ilhkv

	integer ie,ii,k,lmax

	do k=1,nkn
	  ilhkv(k) = 0
	end do

	do ie=1,nel
	  lmax = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ilhkv(k) = max(ilhkv(k),lmax)
	  end do
	end do

	end

c******************************************************
