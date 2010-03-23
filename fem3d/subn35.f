c
c $Id: subn35.f,v 1.17 2009-02-04 15:26:54 georg Exp $
c
c parameter changing area routines
c
c contents :
c
c subroutine bofric(uv,vv,hl)	computes bottom friction
c function cdf(h,z0)		computes cd from h and z0
c
c subroutine rdarea		reads area section (chezy) from STR file
c subroutine ckarea		checks values for chezy parameters
c subroutine prarea		prints chezy values to log file
c subroutine tsarea		prints test message to terminal
c subroutine inarea		initializes chezy values
c
c revision log :
c
c revised on 31.08.88 by ggu  (writes real chezy on czv)
c revised on 29.11.88 by ggu  (new chezy, iarv array)
c revised on 12.04.90 by ggu  (href)
c revised on 03.06.90 by ggu  (austausch)
c revised on 26.06.97 by ggu  (implicit none, useless parts deleted)
c 25.05.1998	ggu	documentation started
c 21.08.1998	ggu	xv eliminated
c 25.05.1999	ggu	new routine bofric
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 09.08.2003    ggu     bofric now returns array with friction
c 10.08.2003    ggu     completely restructured, counter from 0 to nczdum
c 04.09.2003    ggu     bug fix for missing return in get_chezy_values
c 11.01.2005    ggu     ausv eliminated (was not used anymore)
c 02.04.2007    ggu     in check -> warning only for cz=0 and Chezy/Strickler
c 10.12.2008    ggu     re-organized, deleted sp135r(), use bottom_friction()
c 29.01.2009    ggu     ausdef eliminated (chezy(5,.) is not used anymore)
c
c***********************************************************
c
c-------------------------------------------------------------------
c
c DOCS  START   P_friction
c
c DOCS  FRICTION		Bottom friction
c
c The friction term in the momentum equations can be written as
c $Ru$ and $Rv$ where $R$ is the variable friction coefficient and
c $u,v$ are the velocities in $x,y$ direction respectively.
c The form of $R$ can be specified in various ways. The value of 
c |ireib| is choosing between the formulations. In the parameter
c input file a value $\lambda$ is specified that is used in 
c the formulas below.
c
c |ireib|	Type of friction used (default 0):
c		\begin{description}
c		\item[0] No friction used
c		\item[1] $R=\lambda$ is constant
c		\item[2] $\lambda$ is the Strickler coefficient.
c			 In this formulation $R$ is written as
c			 $R = \frac{g}{C^2} \frac{\vert u \vert}{H}$
c			 with $C=k_s H^{1/6}$ and $\lambda=k_s$ is
c			 the Strickler coefficient. In the above
c			 formula $g$ is the gravitational acceleration,
c			 $\vert u \vert$ the modulus of the current velocity
c			 and $H$ the total water depth.
c		\item[3] $\lambda$ is the Chezy coefficient.
c			 In this formulation $R$ is written as
c			 $R = \frac{g}{C^2} \frac{\vert u \vert}{H}$
c			 and $\lambda=C$ is the Chezy coefficient.
c		\item[4] $R=\lambda/H$ with $H$ the total water depth
c		\item[5] $R=\lambda\frac{\vert u \vert}{H}$
c		\end{description}
c |czdef|	The default value for the friction parameter $\lambda$.
c		Depending on the value of |ireib| the coefficient $\lambda$
c		is describing linear friction, constant drag coefficient
c		or a Chezy or Strickler
c		form of friction (default 0).
c |iczv|	Normally $R$ is evaluated at every time step (|iczv| = 1).
c		If for some reason this behavior is not desirable,
c		|iczv| = 0 evaluates the value of $R$ only before the
c		first time step, keeping it constant for the
c		rest of the simulation. (default 1)
c
c The value of $\lambda$ may be specified for the whole basin through
c the value of |czdef|. For more control over the friction parameter
c it can be also specified in section |area| where the friction
c parameter depending on the type of the element may be varied. Please
c see the paragraph on section |area| for more information.
c
c DOCS  END
c
c-------------------------------------------------------------------
c
c***********************************************************

	subroutine bottom_friction

c computes bottom friction

	implicit none

	include 'param.h'

	real drittl
	parameter(drittl=1./3.)

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real grav,fcor,dcor,dirn,rowass,roluft
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        integer nen3v(3,1)
        common /nen3v/nen3v
        integer ilhv(1)
        common /ilhv/ilhv
        real czv(1)
        common /czv/czv
        real rfricv(1)
        common /rfricv/rfricv
        real z0bk(1)
        common /z0bk/z0bk

        real utlov(nlvdim,1),vtlov(nlvdim,1)
        common /utlov/utlov, /vtlov/vtlov
        real hdeov(nlvdim,1)
        common /hdeov/hdeov

	integer ie,ii,k,lmax
	integer ireib
	real hzg
	real hzoff
	real uso,vso,uv
	real rlamb
	real rfric,raux,rr,ss

	real getpar,cdf

c-------------------------------------------------------------------
c get variables
c-------------------------------------------------------------------

	rlamb = getpar('czdef')			! not used
	hzoff = getpar('hzoff')
	ireib = nint(getpar('ireib'))

c-------------------------------------------------------------------
c loop on elements
c-------------------------------------------------------------------

	do ie=1,nel

c         ----------------------------------------------------------
c	  get transport in layer
c         ----------------------------------------------------------

	  lmax = ilhv(ie)

          uso = utlov(lmax,ie)
          vso = vtlov(lmax,ie)
	  uv = sqrt(uso*uso+vso*vso)

c         ----------------------------------------------------------
c	  set total depth
c         ----------------------------------------------------------

	  hzg = hdeov(lmax,ie)
          if( hzg .lt. hzoff ) hzg = hzoff

c         ----------------------------------------------------------
c	  get friction parameter
c         ----------------------------------------------------------

	  rfric = czv(ie)

c         ----------------------------------------------------------
c	  compute friction
c         ----------------------------------------------------------

	  if(ireib.eq.0) then
		rr = 0.
          else if(ireib.eq.1) then
                rr = rfric
	  else if(ireib.eq.2) then		! Strickler
		raux = grav/((rfric**2)*(hzg**drittl))
		rr = raux*uv/(hzg*hzg)
	  else if(ireib.eq.3) then		! Chezy
		raux = grav/(rfric**2)
		rr = raux*uv/(hzg*hzg)
          else if(ireib.eq.4) then
                rr = rfric/hzg
	  else if(ireib.eq.5) then		! constant drag coefficient
		rr = rfric*uv/(hzg*hzg)
	  else if(ireib.eq.6) then		! rfric is z0
                raux = cdf(hzg,rfric)
		rr = raux*uv/(hzg*hzg)
          else if(ireib.eq.7) then		! mixed Strickler / drag
                if( rfric .ge. 1. ) then
		  raux = grav/((rfric**2)*(hzg**drittl))
                else
		  raux = rfric
                end if
		rr = raux*uv/(hzg*hzg)
          else if(ireib.eq.8) then		! use z0 computed by sedtrans
                ss = 0.
                do ii=1,3
                  k = nen3v(ii,ie)
                  ss = ss + z0bk(k)
                end do
                ss = ss / 3.
                raux = cdf(hzg,ss)
		rr = raux*uv/(hzg*hzg)
	  else
		write(6,*) 'unknown friction : ',ireib
		stop 'error stop bottom_friction'
	  end if

	  rfricv(ie) = rr

	end do

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	end

c***********************************************************

        function cdf(h,z0)

c computes cd from h and z0

        implicit none

        real cdf
        real h,z0

        real kappa,cds

        kappa = 0.4

        cds = kappa / log( (z0+0.5*h) / z0 )

        cdf = cds*cds

        end

c***********************************************************
c***********************************************************
c***********************************************************
c***********************************************************
c***********************************************************

	subroutine n_chezy_values(nareas)

	implicit none

	integer nareas

	include 'chezy.h'

	nareas = nczdum

	end

c***********************************************************

	subroutine get_chezy_values(iar,valin,valout)

	implicit none

	include 'chezy.h'

	integer iar
	real valin,valout

	if( iar .gt. nczdum ) goto 99

	valin = czdum(1,iar)
	valout = czdum(2,iar)

	return
   99	continue
	write(6,*) 'iar,nczdum: ',iar,nczdum
	stop 'error stop get_chezy_values'
	end

c***********************************************************

	subroutine set_chezy_values(iar,valin,valout)

	implicit none

	include 'chezy.h'

	integer iar
	real valin,valout

	if( iar .gt. nczdum ) goto 99

	czdum(1,iar) = valin
	czdum(2,iar) = valout

	return
   99	continue
	write(6,*) 'iar,nczdum: ',iar,nczdum
	stop 'error stop set_chezy_values'
	end

c***********************************************************
c***********************************************************
c***********************************************************
c***********************************************************
c***********************************************************

	subroutine set_chezy

c initializes chezy arrays

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	include 'chezy.h'

	integer iarv(1)
	common /iarv/iarv
	real czv(1)
	common /czv/czv

	integer ie,iar

	do ie=1,nel
	    iar=iarv(ie)
	    czv(ie)=czdum(6,iar)
	end do

	end

c***********************************************************

	subroutine init_chezy

c initializes chezy arrays

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	include 'chezy.h'

	integer iarv(1)
	common /iarv/iarv
	real czv(1)
	common /czv/czv

	logical bdebug
	integer i

	bdebug = .true.
	bdebug = .false.

	do i=0,nczdum
	  czdum(6,i)=czdum(1,i)
	end do

	if( bdebug ) call print_chezy

	call set_chezy

	end

c***********************************************************

	subroutine adjust_chezy

c adjusts chezy arrays

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	include 'chezy.h'

	integer iarv(1)
	common /iarv/iarv
	real czv(1)
	common /czv/czv
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
        real up0v(1), vp0v(1)
        common /up0v/up0v, /vp0v/vp0v

	logical bdebug
	integer i,k1,k2
	integer iczv
	real dx,dy,scal

	real getpar

	bdebug = .true.
	bdebug = .false.

	iczv=nint(getpar('iczv'))
	if( iczv .eq. 0 ) return	!chezy is not adjusted

	do i=0,nczdum
	    if(czdum(2,i).eq.0.) then
		czdum(6,i)=czdum(1,i)
	    else
		k1=nint(czdum(3,i))
		k2=nint(czdum(4,i))
		dx=xgv(k2)-xgv(k1)
		dy=ygv(k2)-ygv(k1)
		scal=dx*up0v(k1)+dy*vp0v(k1)
		if(scal.ge.0.) then
			czdum(6,i)=czdum(1,i)
		else
			czdum(6,i)=czdum(2,i)
		end if
	    end if
	end do

	if( bdebug ) call print_chezy

	call set_chezy

	end

c***********************************************************

	subroutine print_chezy

c prints chezy arrays

	implicit none

	include 'chezy.h'

	integer i
	integer iunit

	iunit = 6

	write(iunit,*) 'Values for chezy (czv) :'
	do i=0,nczdum
	  write(iunit,*) i,czdum(6,i)
	end do

	end

c***********************************************************

	subroutine check_chezy

c checks chezy arrays

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	include 'chezy.h'

	integer iarv(1)
	common /iarv/iarv
	real czv(1)
	common /czv/czv

	integer ie,iar
	integer i,j,k
	real cz

	do i=0,nczdum
	  cz = czdum(1,i)
	  if( cz .lt. 0. .or. cz .gt. 1.e+10 ) goto 99
	  cz = czdum(2,i)
	  if( cz .lt. 0. .or. cz .gt. 1.e+10 ) goto 99
	  if( cz .gt. 0. ) then
	    k = nint(czdum(3,i))
	    if( k .lt. 0. .or. k .gt. nkn ) goto 99
	    k = nint(czdum(4,i))
	    if( k .lt. 0. .or. k .gt. nkn ) goto 99
	  end if
	  cz = czdum(6,i)
	  if( cz .lt. 0. .or. cz .gt. 1.e+10 ) goto 99
	end do

	do ie=1,nel
	    iar=iarv(ie)
	    if( iar .gt. nczdum ) goto 98
	    cz=czdum(6,iar)
	    if( cz .lt. 0. .or. cz .gt. 1.e+10 ) goto 98
	end do

	return
   98	continue
	write(6,*) 'ie,iar,nczdum,cz: ',ie,iar,nczdum,cz
	write(6,*) (czdum(j,iar),j=1,6)
	stop 'error stop check_chezy: error in values (1)'
   99	continue
	write(6,*) 'i,iar,nczdum: ',i,i-1,nczdum
	write(6,*) (czdum(j,i),j=1,6)
	stop 'error stop check_chezy: error in values (2)'
	end

c***********************************************************
c***********************************************************
c***********************************************************
c***********************************************************
c***********************************************************

	subroutine rdarea

c reads area section (chezy) from STR file

	implicit none

	include 'chezy.h'

	character*80 line
	integer ianz,iar,i
	real f(20)

	integer nrdlin,iscan,iround

c DOCS	START	S_area
c
c DOCS	AREA	description of area code (Chezy)

	call inarea		! just to be sure

	do while( nrdlin(line) .ne. 0 )
		ianz = iscan(line,1,f)
		if( ianz .gt. 0 ) then
			iar = iround(f(1))
                        if(iar.lt.0) goto 88
			if( ianz .gt. 7 ) goto 86
                        if(iar.gt.nczdum) nczdum=iar
                        if(iar.gt.nardim) goto 75
			do i=2,ianz
			  czdum(i-1,iar) = f(i)
			end do
		else if( ianz .lt. 0 ) then
			goto 98
		end if
	end do

	return
   75   continue
        write(6,*) 'Dimension error for nardim'
        write(6,*) '   iar :',   iar
        write(6,*) 'nardim :',nardim
        stop 'error stop : rdarea'
   86   continue
        write(6,*) 'Too much data on line'
        write(6,*) line
        stop 'error stop : rdarea'
   88   continue
        write(6,*) 'Negative area code = ',iar,' not allowed'
        write(6,*) line
        stop 'error stop : rdarea'
   98   continue
        write(6,*) 'Read error in line :'
	write(6,*) line
        stop 'error stop : rdarea'
	end

c***********************************************************

	subroutine ckarea

c checks values for chezy parameters

	implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	include 'chezy.h'

        integer iarv(1)
        common /iarv/iarv

	integer i,knode,knodeh,ireib
	logical bstop,bpos
	real czdef

	integer iround,ipint
	real getpar

	bstop = .false.

c get default values

        ireib=nint(getpar('ireib'))
	!if( ireib .le. 0 ) return

	bpos = ireib .gt. 1 .and. ireib .ne. 5	!must be strictly positive

        czdef=getpar('czdef')

        if(czdum(1,0).gt.0.) czdef=czdum(1,0)

c compute maximum value of area code

	nczdum = 0
	do i=1,nel
          if(iarv(i).gt.nczdum) nczdum=iarv(i)
        end do

        if(nczdum.gt.nardim) then
                write(6,*) 'section AREA : dimension nardim  ',nardim
                write(6,*) '               maximum area code ',nczdum
                nczdum=nardim
                bstop=.true.
        end if

c check read in values

        do i=0,nczdum

         if(czdum(1,i).eq.-1.) czdum(1,i)=czdef

         if(czdum(1,i).lt.0.) then
                write(6,*) 'Friction value cannot be negative:'
                write(6,*) 'area = ',i,'  chezy = ',czdum(1,i)
                bstop=.true.
	 end if

	 if( bpos .and. czdum(1,i).eq.0.) then
                write(6,*) 'Friction value must be positive:'
                write(6,*) 'area = ',i,'  chezy = ',czdum(1,i)
                bstop=.true.
	 end if

         if(czdum(2,i).eq.-1.) czdum(2,i)=0.

         if(czdum(3,i).eq.-1. .or. czdum(3,i).eq.0.) then
           czdum(3,i)=0.
         else
           knodeh=iround(czdum(3,i))
           knode=ipint(knodeh)          !$$EXTINW
           if(knode.le.0) then
                write(6,*) 'section AREA : node not found ',knodeh
                bstop=.true.
           end if
           czdum(3,i)=knode
         end if

         if(czdum(4,i).eq.-1. .or. czdum(4,i).eq.0.) then
           czdum(4,i)=0.
         else
           knodeh=iround(czdum(4,i))
           knode=ipint(knodeh)          !$$EXTINW
           if(knode.le.0) then
                write(6,*) 'section AREA : node not found ',knodeh
                bstop=.true.
           end if
           czdum(4,i)=knode
         end if

         czdum(6,i)=0.

        end do

	if( bstop ) stop 'error stop ckarea'

	end

c***********************************************************

	subroutine prarea

c prints chezy values to log file

	implicit none

	include 'chezy.h'

	integer ianf,i
	integer ipext,iround

        ianf=0
        if(czdum(1,0).eq.0) ianf=1
        write(6,*)
        write(6,1007)

        do i=ianf,nczdum
            if(czdum(2,i).ne.0.) then			!with two chezy
                write(6,1008) i,czdum(1,i),czdum(2,i)
     +                          ,ipext(iround(czdum(3,i)))
     +                          ,ipext(iround(czdum(4,i)))
            else					!just one chezy
                write(6,1008) i,czdum(1,i)
            end if
        end do

	return
 1007   format(' area,cz1,cz2,k1,k2 : ')
 1008   format(i5,2f7.1,2i7,e12.4)
	end

c***********************************************************

	subroutine tsarea

c prints test message to terminal

	implicit none

	include 'chezy.h'

	integer j,i

        write(6,*) '/chezy/'
        write(6,*) nczdum
        do j=0,nczdum
            write(6,'(1x,6e12.4)') (czdum(i,j),i=1,6)
        end do

	end

c***********************************************************

	subroutine inarea

c initializes chezy values

	implicit none

	include 'chezy.h'

	integer i,j

	nczdum = 0
        do i=0,nardim
         do j=1,6
           czdum(j,i) = -1.
         end do
	end do

	end

c***********************************************************
