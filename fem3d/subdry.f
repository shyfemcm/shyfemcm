c
c $Id: subdry.f,v 1.14 2008-11-20 10:51:34 georg Exp $
c
c routines handling flooding and drying
c
c contents :
c
c subroutine setweg(iweich,iw,hzmin)	sets array iwegv
c subroutine setuvd			sets velocities in dry areas
c subroutine setzev			sets array zenv from znv
c subroutine setznv			sets array znv from zenv
c subroutine zuniq(zv,av)		makes z values unique
c
c revision log :
c
c revised 01.07.92 by ggu   $$lump  - lumping of matrix
c revised 05.08.92 by ggu   $$ibtyp3 - implementation of ibtyp=3
c 27.03.1998	ggu	eliminated /bnd/, /irv/
c 27.04.1998	ggu	$$NKNEL - do not call nknel in tstlnk
c 08.05.1998	ggu	new routine pntfla -> absolute element index
c 20.05.1998	ggu	hard coded unit 88 substituted (use ifileo)
c 14.07.1998	ggu	$$ibtyp4 - boundary type 4 integrated
c 21.08.1998	ggu	file sublnk splitted into lnk/dry
c 21.08.1998	ggu	routine setczg removed
c 21.08.1998    ggu     xv eliminated
c 21.11.2001    ggu     extra bdebug in setweg
c 21.11.2001    ggu     more debug information in setuvd
c 13.08.2003    ggu     new routine setznv
c 03.09.2004    ggu     setznv: do not stop if znv is not unique (restart)
c 22.09.2004    ggu     debug in setweg()
c 23.03.2006    ggu     changed time step to real
c
c*****************************************************************
c
        subroutine setweg(iweich,iw)
c
c sets array iwegv
c
c iweich 	flag
c	-1: initialize  
c	 0: first call  
c	 1: only take away (hzmin)
c	 2: only take away (hzoff)    
c	 3: only add
c iw		if on return iw != 0 iwegv has changed
c
c 1 is used to check, if iteration has to be repeated
c 2 is just to switch off an element in time, so no iteration
c	has to be repeated
c generally it is true that : 0 <= hzmin <= zmin <= hzoff <= hzon
c
c iwegv   0:all nodes wet   >0:number of nodes dry -> out of system
c
c revised 12.01.94 by ggu   $$hzon  - use new variable hzon
c
        implicit none
c
c arguments
        integer iweich,iw
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer iwegv(1),nen3v(3,1)
        common /iwegv/iwegv, /nen3v/nen3v
        real hm3v(3,1)
        common /hm3v/hm3v
        real zenv(3,1)
        common /zenv/zenv
	real znv(1)
	common /znv/znv
	include 'ev.h'
c local
        integer ie,ii,iwh,iweg,k
        integer iespec
        real hzg,hzmin,hzoff,hzon,volmin,aomega,zmin
c functions
        real getpar
c	integer ieint,iround,ipint
c aux
	logical debug,bnodry
        logical bdebug
c	real rz,rh,rx,rxmin
c	integer k,iex
c	integer iz(30),ih(30),izh(30),ix(30),ixx(30)
        integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it
c
	debug=.false.
	bnodry=.true.		!drying is not allowd (for debug)
	bnodry=.false.
        bdebug=.true.
        bdebug=.false.

        iespec = 4574
        iespec = -1
c
	if( iweich .ge. 0 ) then
          hzmin=getpar('hzmin')
          hzoff=getpar('hzoff')
          hzon=getpar('hzon')			!$$hzon
          volmin=getpar('volmin')
	end if
c
        iwh=0
c
        if(iweich.eq.-1) then                   !initialize iwegv
          do ie=1,nel
            iwegv(ie)=0
          end do
        else if(iweich.eq.0) then               !first call
          do ie=1,nel
            aomega=ev(10,ie)
            zmin=hzmin+volmin/(4.*aomega)
c            zmin=hzoff+volmin/(4.*aomega)   !$$new 24.07.92
            if(zmin.gt.hzoff) zmin=hzmin+0.5*(hzmin+hzoff) !just in case...
c
            iweg=0
            do ii=1,3
              hzg = zenv(ii,ie) + hm3v(ii,ie)
              if(hzg.lt.hzoff) iweg=iweg+1
              if(hzg.lt.zmin) zenv(ii,ie)=zmin-hm3v(ii,ie)
              if(bdebug.and.hzg.lt.hzoff) then
                      write(6,*) 'setweg 0: ',ie,iweg,hzg
                      write(6,*) znv(nen3v(ii,ie)),hm3v(ii,ie)
              end if
            end do
            iwegv(ie)=iweg
          end do
        else if(iweich.eq.1) then               !only take away (hzmin)
          do ie=1,nel
            debug = ie .eq. iespec
            iweg=0
            do ii=1,3
              hzg = znv(nen3v(ii,ie)) + hm3v(ii,ie)
              if(hzg.lt.hzmin) iweg=iweg+1
              if(debug .or. bdebug.and.hzg.lt.hzmin) then
                      write(6,*) 'setweg 1: ',ie,iweg,hzg
                      write(6,*) znv(nen3v(ii,ie)),hm3v(ii,ie)
              end if
            end do
            if(iweg.gt.0) then
              if(iwegv(ie).eq.0) iwh=iwh+1
	      if( bnodry ) goto 77
              iwegv(ie)=iweg
            end if
          end do
        else if(iweich.eq.2) then               !only take away (hzoff)
          do ie=1,nel
            debug = ie .eq. iespec
            iweg=0
            do ii=1,3
              hzg = znv(nen3v(ii,ie)) + hm3v(ii,ie)
              if(hzg.lt.hzoff) iweg=iweg+1
              if(debug .or. bdebug.and.hzg.lt.hzoff) then
                      write(6,*) 'setweg 2: ',ie,iweg,hzg
                      write(6,*) znv(nen3v(ii,ie)),hm3v(ii,ie)
              end if
            end do
            if(iweg.gt.0) then
              if(iwegv(ie).eq.0) iwh=iwh+1
	      if( bnodry ) goto 77
              iwegv(ie)=iweg
            end if
          end do
        else if(iweich.eq.3) then               !only add
          do ie=1,nel
            debug = ie .eq. iespec
            iweg=0
            do ii=1,3
              hzg = znv(nen3v(ii,ie)) + hm3v(ii,ie)
              if(hzg.lt.hzon) iweg=iweg+1 !$$hzon ...(220792)
              if(debug .or. bdebug.and.hzg.lt.hzoff) then
                      write(6,*) 'setweg 3: ',ie,iweg,hzg
                      write(6,*) znv(nen3v(ii,ie)),hm3v(ii,ie)
              end if
            end do
            if(iweg.eq.0) then
              if(iwegv(ie).gt.0) iwh=iwh+1
              iwegv(ie)=0
            end if
          end do
        end if
c
        iw=iwh
c
        return
   77	continue
	write(6,*) 'drying is not allowed...'
	write(6,*) iweich,ie,iwh,iwegv(ie)
	write(6,*) hzmin,hzoff,hzon,volmin
        do ii=1,3
	  k = nen3v(ii,ie)
          hzg = znv(k) + hm3v(ii,ie)
	  write(6,*) ii,k,znv(k),hzg
	  write(6,*) '     ',zenv(ii,ie),hm3v(ii,ie)
	end do
	stop 'error stop setweg'
        end
c
c****************************************************************
c
        subroutine setuvd
c
c sets velocities in dry areas
c
c ie    element
c dt    time step
c hzmin smallest z allowed
c b,c   form functions
c
c revised ...07.92 by ggu   $$lump  - lumping of matrix
c revised ......92 by ggu   $$eps  - introduction of eps
c revised 12.01.94 by ggu   $$eps0  - use eps only in last control
c revised 12.01.94 by ggu   $$99  - do not jump to 99 in loop
c revised 05.02.94 by ggu   $$azpar - use az to compute velocities
c revised 04.03.94 by ggu   $$azuvdry - one az too much in formula
c revised 27.10.97 by ggu   $$isum - better identification of error 99
c revised 27.10.97 by ggu   $$dpisum - use double prec. for key values
c
        implicit none
c
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it
        integer iwegv(1),nen3v(3,1)
        common /iwegv/iwegv, /nen3v/nen3v
        real zenv(3,1),hm3v(3,1)
        common /zenv/zenv, /hm3v/hm3v
        real uov(1),vov(1),unv(1),vnv(1)
	common /uov/uov, /vov/vov, /unv/unv, /vnv/vnv
	include 'ev.h'
c local
        integer ie,ii,i1,i2,isum,itot
        integer i3,i4,i5,i6,i7,i8
c        real zm,zmed,d1,d2,det
        double precision zm,zmed,d1,d2,det	!$$dpisum
        real adt,axdt,dt,hzmin
	real az,azt,azpar
        real z(3),b(3),c(3),uo,vo,u,v,zz
	real zn(3)
        integer nnn,itmin
c save
        real eps
        save eps
c functions
        real getpar
        integer ieext
        logical iseout
        iseout(ie) = iwegv(ie).gt.0
c data
        data eps /1.e-5/
c
        nnn=-802 	!<0 if not needed
c	itmin=1550100	!0 if not needed
	itmin=0		!0 if not needed
c
	call get_timestep(dt)
        hzmin=getpar('hzmin')
	call getaz(azpar)
	az = azpar
	azt=1.-az
c
        i7=0            !mean value is too low for node
          i8=0          !
            i3=0
              i4=0
                i5=0    !two nodes in element with mean too low
                  i6=0  !node is drying with computed velocity

        do ie=1,nel
        if( iseout(ie) ) then
c
        uo=uov(ie)	!use barotropic velocities
        vo=vov(ie)
	axdt=3.*dt    !$$lump $$azpar
        adt=1./axdt   !$$lump
c
c z average, set b,c
c
        zm=0.
        do ii=1,3
	  zn(ii) = zenv(ii,ie)
          zm=zm+zn(ii)
          b(ii)=ev(3+ii,ie)
          c(ii)=ev(6+ii,ie)
        end do
        zmed=zm/3.
c
c compute final z value for the nodes
c
        itot=0
        isum=0
c        zm=0.
        do ii=1,3
          if(zmed+hm3v(ii,ie).lt.hzmin) then  !$$eps0
            z(ii)=hzmin-hm3v(ii,ie)
            itot=itot+1
            isum=isum+ii
            i7=i7+1
          else
            z(ii)=zmed
          end if
          zm=zm-z(ii)
        end do
c
        if(itot.eq.0) then
          !ok
        else if(itot.eq.1) then
          i1=mod(isum,3)+1
          i2=mod(i1,3)+1
          if(z(i1)+zm*0.5+hm3v(i1,ie).lt.hzmin) then  !$$eps0
            i8=i8+1
            z(i1)=hzmin-hm3v(i1,ie)
            z(i2)=3.*zmed-z(i1)-z(isum)
c            if(z(i2)+hm3v(i2,ie).lt.hzmin-eps) goto 99 !$$99
          else if(z(i2)+zm*0.5+hm3v(i2,ie).lt.hzmin) then !$$eps0
            i3=i3+1
            z(i2)=hzmin-hm3v(i2,ie)
            z(i1)=3.*zmed-z(i2)-z(isum)
c            if(z(i1)+hm3v(i1,ie).lt.hzmin-eps) goto 99  !$$99
          else
            i4=i4+1
            z(i1)=z(i1)+0.5*zm
            z(i2)=z(i2)+0.5*zm
          end if
        else if(itot.eq.2) then
          i5=i5+1
          i1=6-isum
          z(i1)=z(i1)+zm
c          if(z(i1)+hm3v(i1,ie).lt.hzmin-eps) goto 99  !$$99
        else
          goto 99
        end if
c
c control, leave in any case
c
        isum=-1
        zm=0.
        do ii=1,3
          zm=zm+z(ii)
          if(z(ii)+hm3v(ii,ie).lt.hzmin-eps) goto 99
        end do
        isum=-2		!$$isum
        if(abs(zm-3.*zmed).gt.eps) goto 99
c
c now compute velocities
c						!$$azpar
        d1 = (z(1)-zn(1))*adt
     +          - azt*( b(1)*uo + c(1)*vo )
        d2 = (z(2)-zn(2))*adt
     +          - azt*( b(2)*uo + c(2)*vo )
        det=1./(b(1)*c(2)-b(2)*c(1))
c
        u = det * ( c(2)*d1 - c(1)*d2 )
        v = det * (-b(2)*d1 + b(1)*d2 )
c
c with this velocity is there some node drying out ?
c
	isum=-1
        itot=0
        do ii=1,3				!$$azpar
          zz = axdt * ( b(ii)*(azt*uo+az*u) + c(ii)*(azt*vo+az*v) ) 
     +			+ zn(ii)
          if(zz+hm3v(ii,ie).lt.hzmin-eps) itot=itot+1 !$$eps
        end do
c
        if(itot.gt.0) then    !node is drying, set next z to above values
	  isum=-2
          i6=i6+1
c						!$$azpar
          d1 = (z(1)-zn(1))*adt
     +          - azt*( b(1)*uo + c(1)*vo )
          d2 = (z(2)-zn(2))*adt
     +          - azt*( b(2)*uo + c(2)*vo )
          det=1./( az * (b(1)*c(2)-b(2)*c(1)) )	!$$azuvdry
c
c the formula should be   det=1./( az*az * (b(1)*c(2)-b(2)*c(1)) )
c and in the next two lines   u = det * az * (...)   and   v = ...
c but we devide by az both equations
c
          u = det * ( c(2)*d1 - c(1)*d2 )
          v = det * (-b(2)*d1 + b(1)*d2 )
        end if
c
c now set u/v/z
c
        itot=0
        zm=0.          !only for control
        do ii=1,3				!$$azpar
          zz = axdt * ( b(ii)*(azt*uo+az*u) + c(ii)*(azt*vo+az*v) ) 
     +			+ zn(ii)
          if(zz+hm3v(ii,ie).lt.hzmin-eps) itot=itot+1 !$$eps
          zenv(ii,ie)=zz
          zm=zm+zz
        end do
        isum=-3		!$$isum
        if(abs(zm-3.*zmed).gt.eps) goto 99
c
        if(itot.gt.0) goto 97

        unv(ie)=u	!put back to new barotropic velocities
        vnv(ie)=v
c
        end if
        end do

c        write(88,*) i7,i8,i3,i4,i5,i6
c
        return

   99   continue  !use isum to decide which branch has been taken
        write(6,*) '---------------------'
	write(6,*) 'error log setuvd (z computation)'
        write(6,*) 'time: ',it
        write(6,*) 'element number: ',ieext(ie)
        write(6,*) ie,itot,isum
        write(6,*) zmed,zm,hzmin
        write(6,*) (hm3v(ii,ie),ii=1,3)
        write(6,*) (zenv(ii,ie),ii=1,3)
        write(6,*) (z(ii),ii=1,3)
        write(6,*) '---------------------'
        stop 'error stop setuvd : z computation'
   97   continue
        write(6,*) '---------------------'
	write(6,*) 'error log setuvd (u/v computation)'
        write(6,*) 'time: ',it
        write(6,*) 'element number: ',ieext(ie)
        write(6,*) ie,itot,isum
        write(6,*) zmed,zm,hzmin
        write(6,*) uo,vo,u,v
        write(6,*) (hm3v(ii,ie),ii=1,3)
        write(6,*) (zenv(ii,ie),ii=1,3)
        write(6,*) (z(ii),ii=1,3)
        write(6,*) (zn(ii),ii=1,3)
        write(6,*) '---------------------'
        stop 'error stop setuvd : u/v computation'
        end
c
c****************************************************************

        subroutine setzev

c sets array zenv from znv

        implicit none

c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1),iwegv(1)
        common /nen3v/nen3v, /iwegv/iwegv
        real zenv(3,1)
        common /zenv/zenv
	real znv(1)
	common /znv/znv
c local
        integer ie,ii

        do ie=1,nel
          if(iwegv(ie).eq.0) then !element is in system
            do ii=1,3
              zenv(ii,ie)=znv(nen3v(ii,ie))
            end do
          end if
	end do

        end

c****************************************************************

        subroutine setznv

c sets array znv from zenv

        implicit none

c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real eps1,eps2,pi,flag,high,hihi
	common /mkonst/ eps1,eps2,pi,flag,high,hihi

        integer nen3v(3,1),iwegv(1)
        common /nen3v/nen3v, /iwegv/iwegv
        real zenv(3,1)
        common /zenv/zenv
	real znv(1)
	common /znv/znv
	include 'ev.h'
	real v1v(1), v2v(1)
	common /v1v/v1v, /v2v/v2v
c local
        integer ie,ii,k
        integer ntot
	real z,area

c-------------------------------------------------------------
c initialize znv and counters
c-------------------------------------------------------------

        ntot = 0

	do k=1,nkn
	  znv(k) = flag
	  v1v(k) = 0.
	  v2v(k) = 0.
	end do

c-------------------------------------------------------------
c set znv and accumulate
c-------------------------------------------------------------

        do ie=1,nel
          if( iwegv(ie) .eq. 0 ) then		!element is in system
            do ii=1,3
	      k = nen3v(ii,ie)
              z = zenv(ii,ie)
	      if( znv(k) .ne. flag .and. znv(k) .ne. z ) then   !restart
                ntot = ntot + 1
                write(6,*) 'n,ie,ii,k,z,znv(k) ',ntot,ie,ii,k,z,znv(k)
                z=max(z,znv(k))         !using higher value
              end if
	      znv(k) = z
            end do
	  else					!out of system
	    area = 4. * ev(10,ie)
            do ii=1,3
	      k = nen3v(ii,ie)
              z = zenv(ii,ie)
	      v1v(k) = v1v(k) + z * area
	      v2v(k) = v2v(k) + area
            end do
          end if
	end do

c-------------------------------------------------------------
c compute znv for dry areas
c-------------------------------------------------------------

	do k=1,nkn
	  if( znv(k) .eq. flag ) then		!out of system
	    znv(k) = v1v(k) / v2v(k)
	  end if
	end do

c-------------------------------------------------------------
c write debug status
c-------------------------------------------------------------

        if( ntot .gt. 0 ) then
          write(6,*) ntot, ' nodes with different z value'
        end if

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	return
   99	continue
	write(6,*) 'ie,ii,k,z,znv(k) ',ie,ii,k,z,znv(k)
	stop 'error stop setznv: nodal value not unique'
        end
c
c****************************************************************
c
        subroutine zuniq(zv,av)
c
c makes z values in dynamic system unique
c
c works only for lumped mass matrix
c
c zv    aux vector for z value
c av    aux vector for weighting factors (areas)
c
c written ...07.92 by ggu   $$lump  - lumping of matrix
c
        implicit none
c
c arguments
        real zv(1),av(1)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer iwegv(1),nen3v(3,1)
        real zenv(3,1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /iwegv/iwegv, /nen3v/nen3v
        common /zenv/zenv
	include 'ev.h'
c local
        integer ie,i,k
        real aomega
c functions
        logical isein
        isein(ie) = iwegv(ie).eq.0

c initialize aux vectors

        do k=1,nkn
          zv(k)=0.
          av(k)=0.
        end do

c accumulate contributions

        do ie=1,nel
          if(isein(ie)) then
            aomega=ev(10,ie)
            do i=1,3
              k=nen3v(i,ie)
              zv(k)=zv(k)+aomega*zenv(i,ie)
              av(k)=av(k)+aomega
            end do
          end if
        end do

c scaling of z values with weighting functions

        do k=1,nkn
          if(av(k).gt.0.) then
            zv(k)=zv(k)/av(k)
          end if
        end do

c write back to original vector zenv

        do ie=1,nel
          if(isein(ie)) then
            do i=1,3
              zenv(i,ie)=zv(nen3v(i,ie))
            end do
          end if
        end do

        return
        end
c
c*****************************************************************
c