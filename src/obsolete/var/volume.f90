
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-2000,2003-2005,2010,2014-2015  Georg Umgiesser
!    Copyright (C) 2018-2020  Georg Umgiesser
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

! subroutines for computing volumes
!
! contents :
!
!        subroutine invola
!        subroutine rdvola
!        subroutine ckvola
!        subroutine prvola
!        subroutine tsvola
!
!        subroutine wrvola(dtime)			write of vol data
!
!        subroutine volareas(kvol,ivol,n,az,vol)	vol in areas
!        function volarea(kvol,ivol,n,az)		vol in area
!
!	 function voltotal(bz)				vol in total basin
!	 function areatotal()				area in total basin
!
!        subroutine volini			initializes vol routines
!        subroutine volin(kvol,ivol,n)		sets up info structure
!
!	 subroutine volinf(nv,n,kvol,idim,ivolm,ivol)	sets up element info
!	 subroutine wrvolinf(n,ielem,volav,areav,hav)	computes average values
!
!	 subroutine linecl(n,kvol,nnodes,kline,xline,yline)	close line
!
!	 subroutine wrtelem(iunit,n,ielems,ietype)
!	 subroutine wrtline(iunit,n,inodes,il,iltype)
!
! revision log :
!
! 30.04.1998	ggu	newly written routines (subpor deleted)
! 07.05.1998	ggu	check nrdveci on return for error
! 08.05.1998	ggu	restructured with new comodity routines
! 13.09.1999	ggu	type of node computed in own routine voltype
! 20.01.2000	ggu	common block /dimdim/ eliminated
! 31.07.2003	ggu	some comments and little restructuring
! 10.08.2003	ggu	do not call setweg, setnod, setkan
! 02.09.2003	ggu	bug fix $$BUGVOLT in wrvola -> summation on tot vol
! 17.09.2003	ggu	new routine volstats -> statistics to file .vvv
! 26.11.2004	ggu	routine volstats commented (use only if needed)
! 22.02.2005	ggu	iflag substituted by v3v
! 23.03.2010	ggu	changed v6.1.1
! 07.03.2014	ggu	changed VERS_6_1_72
! 18.06.2014	ggu	changed VERS_6_1_77
! 26.11.2014	ggu	changed VERS_7_0_7
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 05.05.2015	ggu	changed VERS_7_1_10
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 30.07.2015	ggu	changed VERS_7_1_83
! 23.09.2015	ggu	changed VERS_7_2_4
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 16.02.2020    ggu     femtime eliminated
! 20.03.2022    ggu     upgraded to da_out
! 28.09.2023    ggu     include substituted by module
!
! notes :
!
! These routines can also be used internally to compute the vol
! over various sections. The following calling sequence must be respected:
!
! call n2int(n,kvol,berror)		converts external to internal nodes
! nvols = klineck(n,kvol)		checks array kvol and computes nvols
! call volin(kvol,ivol,n)		initializes ivol
!
! call volareas(kvol,ivol,n,az,vol)	computes vols and returns in vol()
!
! Initialization can be done anytime.
!
! wrvola
!       volareas
!       voltotal
! volini
!       volin
!               volinf
!       volstats        (commented)
!
!******************************************************************

	module volcomp

	implicit none

	integer, parameter :: nfxdim = 1

        integer, save :: nvols,kvold,kvolm,kvol(nfxdim)
        integer, save :: ivolm,ivol(nfxdim)

	end module volcomp

!******************************************************************

        subroutine mod_vol(mode)
 
        implicit none
 
        integer mode
 
        include 'modules.h'
 
	double precision dtime
 
        if( mode .eq. M_AFTER ) then
	   call get_act_dtime(dtime)
           call wrvola(dtime)
        else if( mode .eq. M_INIT ) then
           call invola
        else if( mode .eq. M_READ ) then
           call rdvola
        else if( mode .eq. M_CHECK ) then
           call ckvola
        else if( mode .eq. M_SETUP ) then
           call volini
        else if( mode .eq. M_PRINT ) then
           call prvola
        else if( mode .eq. M_TEST ) then
           call tsvola
        else if( mode .eq. M_BEFOR ) then
!          nothing
        else
           write(6,*) 'unknown mode : ', mode
           stop 'error stop mod_vol'
        end if
 
        end

!******************************************************************

        subroutine invola

! nvols		total number of areas
! kvold		dimension of kvol
! kvolm		total number of nodes defining areas
! kvol()	node numbers defining areas

	use volcomp

        implicit none

        nvols = -1	!must still be initialized
        kvold = 0	!is set later
        kvolm = 0

	ivolm = 0

        end

!******************************************************************

        subroutine rdvola

	use volcomp

        implicit none

	integer nfxdi
        integer nrdveci

	!write(6,*) 'module not yet converted... do not use'
	!stop 'error stop rdvola: cannot use'

	nfxdi = nfxdim

	kvold = nfxdi
        kvolm = nrdveci(kvol,nfxdi)

        if( kvolm .lt. 0 ) then
          if( kvolm .eq. -1 ) then
            write(6,*) 'dimension error nvodin in section $VOL : '      &
     &                          ,nfxdi
          else
            write(6,*) 'read error in section $VOL'
          end if
          stop 'error stop rdvola'
        end if

        end

!******************************************************************

        subroutine ckvola

	use volcomp

        implicit none

	integer k,ii
        logical berror

	call n2int(kvolm,kvol,berror)

        if( berror ) then
		write(6,*) 'error in section $VOL'
		stop 'error stop: ckvola'
	end if

! the real set up is done in volini
! however, at this stage we do not have all the arrays set up,
! so we have to post-pone it until later

        end

!******************************************************************

	subroutine prvola

	use volcomp

	implicit none

	integer nnode,ifirst,ilast
	integer ntotal,ns
	integer i,ii
	real volav,areav,hav

	integer ipext,ieext
	logical bwrite
	logical nextline
	real voltotal,areatotal

!----------------------------------------------------------

	if( nvols .le. 0 ) return

!----------------------------------------------------------

	bwrite = .true.

!----------------------------------------------------------

	write(6,*)
	write(6,*) 'vol section :'
	write(6,*)
	write(6,*) 'nvols,kvolm ',nvols,kvolm
	write(6,*) 'kvold,ivolm ',kvold,ivolm
	write(6,*)

!----------------------------------------------------------

	ns = 0
	nnode = 0

	if( bwrite ) open(69,file='lines.grd',status='unknown')

	do while( nextline(kvol,kvolm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
	  write(6,*) 'volume (line): ',ns,ntotal
	  write(6,*) (ipext(kvol(i)),i=ifirst,ilast)
	  if( bwrite ) call wrtline(69,ntotal,kvol(ifirst),ns,80+ns)
	end do

	if( bwrite ) close(69)

!----------------------------------------------------------

	ns = 0
	nnode = 0

	if( bwrite ) open(69,file='elems.grd',status='unknown')

	volav = voltotal(.false.)
	areav = areatotal()
	hav = volav / areav
	write(6,*) 'volume (elements): ',0,0
	write(6,*) 'average values: ',volav,areav,hav

	do while( nextline(ivol,ivolm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
	  call wrvolinf(ntotal,ivol(ifirst),volav,areav,hav)
	  write(6,*) 'volume (elements): ',ns,ntotal
	  write(6,*) 'average values: ',volav,areav,hav
	  !write(6,*) (ieext(ivol(i)),i=ifirst,ilast)
	  if( bwrite ) call wrtelem(69,ntotal,ivol(ifirst),80+ns)
	end do

	if( bwrite ) close(69)

!----------------------------------------------------------

	end

!******************************************************************

	subroutine tsvola

	use volcomp

	implicit none

	integer i,ii

	write(6,*) '/kvolc/'
	write(6,*) nvols,kvold,kvolm
	write(6,*) (kvol(i),i=1,kvolm)

	write(6,*) '/ivol/'
	write(6,*) ivolm
	write(6,*) (ivol(i),i=1,ivolm)

	end

!******************************************************************

	subroutine wrvola(dtime)

! write of vol data

	use volcomp

	implicit none

	double precision dtime

        integer iscdim
        parameter(iscdim=500)

	integer idtvol
	integer i
	integer it
	real az,rr

	integer iround,ideffi
	real getpar
	double precision dgetpar
	logical has_output_d,next_output_d,is_over_output_d

	real vol(iscdim)

	real volt(0:iscdim)	!accumulator - could be also double precision
        save volt

	real voltotal
	real voltot

        integer nr
        save nr
	double precision, save :: da_out(4)
        integer icall,nbvol,nvers,idfile
        save icall,nbvol,nvers,idfile
        data icall,nbvol,nvers,idfile /0,0,1,538/

! start of code

        if( icall .eq. -1 ) return

! initialization

        if( icall .eq. 0 ) then

		call init_output_d('itmvol','idtvol',da_out)
		call increase_output_d(da_out)
                if( .not. has_output_d(da_out) ) icall = -1

                if( kvolm .le. 0 ) icall = -1
                if( nvols .le. 0 ) icall = -1
                if( icall .eq. -1 ) return

                if( nvols .gt. iscdim ) then
                  stop 'error stop wrvola: dimension iscdim'
                end if

                nr = 0
                do i=0,nvols
                  volt(i) = 0.
                end do

                nbvol=ideffi('datdir','runnam','.vol','unform','new')
                if(nbvol.le.0) then
        	   stop 'error stop wrvola : Cannot open VOL file'
		end if

		idtvol = nint(da_out(1))

                write(nbvol) idfile,nvers
                write(nbvol) nvols+1,kvolm,idtvol
                write(nbvol) (kvol(i),i=1,kvolm)
                write(nbvol) ivolm
                write(nbvol) (ivol(i),i=1,ivolm)

        end if

	write(6,*) 'module not yet converted... do not use'
	stop 'error stop wrvola: cannot use'

! normal call

        icall = icall + 1
	it = nint(dtime)

        if( .not. is_over_output_d(da_out) ) return

!	accumulate results

        nr = nr + 1

	call volareas(ivolm,ivol,vol)

	do i=1,nvols
	  volt(i) = volt(i) + vol(i)
	end do

	volt(0) = volt(0) + voltotal(.true.)	!total basin $$BUGVOLT

        if( .not. next_output_d(da_out) ) return

!	write results

        rr=1./nr

        do i=1,nvols
          vol(i) = volt(i) * rr
        end do
	voltot = volt(0) * rr

        write(nbvol) it,nvols+1,voltot,(vol(i),i=1,nvols)
!	write(6,*) 'section written: ',nvols,nr,it

!	reset variables

        nr = 0
        do i=0,nvols
          volt(i) = 0.
        end do

	end

!******************************************************************

	subroutine volareas(n,ielem,vol)

! computes vol in all areas and returns them in vol

	implicit none

	integer n
	integer ielem(1)
	real vol(1)

	integer nnode,ifirst,ilast,ntotal
	integer ns
	logical nextline
	real volarea

	nnode = 0
	ns = 0

	do while( nextline(ielem,n,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
	  vol(ns) = volarea(ntotal,ielem(ifirst))
	end do

	end

!******************************************************************

	function volarea(n,ielem)

! computes vol in one area

	implicit none

	real volarea
	integer n
	integer ielem(1)

	integer i,ie
	real volume

	real volele

	volume = 0.

	do i=1,n
	  ie = ielem(i)
	  volume = volume + volele(ie,+1)
	end do

	volarea = volume

	end
	  
!******************************************************************

	function voltotal(bz)

! computes vol in total basin

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real voltotal
	logical bz	!if true use new zeta to compute volume

	integer ie,mode
	real volume

	real volele

	volume = 0.

	if( bz ) then
	  mode = +1
	else
	  mode = 0
	end if

	do ie=1,nel
	  volume = volume + volele(ie,mode)
	end do

	voltotal = volume

	end
	  
!******************************************************************

	function areatotal()

! computes area in total basin

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real areatotal


	integer ie
	real area

	real areaele

	area = 0.

	do ie=1,nel
	  area = area + areaele(ie)
	end do

	areatotal = area

	end
	  
!******************************************************************

	subroutine volstats(n,ielem)

! writes statistics of volumes to file

	implicit none

	integer n
	integer ielem(1)

	integer nnode,ifirst,ilast,ntotal
	integer ns
        integer iunit
        integer i,ie
        real volume,area

	logical nextline
        integer ieext
        integer ifemop
	real volele,areaele

        iunit = 79
        iunit = ifemop('.vvv','form','new')
        if( iunit .le. 0 ) then
          write(6,*) '*** volstats: cannot open file to write'
          return
        end if

	nnode = 0
	ns = 0

	do while( nextline(ielem,n,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1

	  volume = 0.
          area = 0.
          do i=ifirst,ilast
	    ie = ielem(i)
	    volume = volume + volele(ie,0)
	    area = area + areaele(ie)
	  end do

          write(iunit,*) ns,ntotal
          write(iunit,*) volume,area,volume/area
          write(iunit,*) (ieext(ielem(i)),i=ifirst,ilast)
	end do

        close(iunit)

	end

!******************************************************************

	subroutine volini

! initializes vol routines finally

	use volcomp
	use mod_geom

	implicit none

	integer idummy

	integer klineck

! more checks for compatibility

	nvols = klineck(kvolm,kvol)

	if( nvols .lt. 0 ) then
	  write(6,*) 'errors in section $VOL'
	  stop 'error stop : volini'
	end if

! be sure all these arrays are initialized (lenkv,ilinkv must exists)

!	call setweg(-1,idummy)
!	call setnod
!	call setkan(kantv)

! now set info structure for sections

	write(6,*)
	write(6,*) 'setting up section vol... ',kvolm
	write(6,*)

	call volin(kvolm,kvol,kvold,ivolm,ivol)
	!call volstats(ivolm,ivol)       !writes stats to file vvv

	end

!******************************************************************

	subroutine volin(n,kvol,idim,ivolm,ivol)

! sets up element info structure ivol(1) from kvol(1) for all areas

	implicit none

	integer n		!size of kvol			(in)
	integer kvol(1)		!node numbers of lines		(in)
	integer idim		!dimension of ivol		(in)
	integer ivolm		!filling of ivol		(out)
	integer ivol(1)		!element numbers of volumes	(out)

	integer nnode,ifirst,ilast,ntotal
	integer nv
	logical nextline

	nv = 0
	nnode = 0

	do while( nextline(kvol,n,nnode,ifirst,ilast) )
	  nv = nv + 1
	  ntotal = ilast - ifirst + 1
!	  write(6,*) n,nnode,ifirst,ilast,ntotal
	  call volinf(nv,ntotal,kvol(ifirst),idim,ivolm,ivol)
	end do

	end

!******************************************************************

	subroutine volinf(nv,n,kvol,idim,ivolm,ivol)

! sets up element info structure ivol(1) from kvol(1) for one area

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nv		!progressive number of line	(in)
	integer n		!size of kvol			(in)
	integer kvol(1)		!node numbers of line		(in)
	integer idim		!dimension of ivol		(in)
	integer ivolm		!filling of ivol		(out)
	integer ivol(1)		!element numbers of volumes     (out)


	real xline(nkn), yline(nkn)
	integer kline(nkn)

	logical bclosed
	logical bdebug
	integer nnodes,ntot
	integer i,ie
	integer ipext
	real xc,yc

	logical inpoly
	integer ieext

	bdebug = .false.

	if( bdebug ) write(6,*) '0 (volinf) closing line ',nv

! see if we can close it

	call linecl(n,kvol,nnodes,kline,xline,yline)
	bclosed = nnodes .ge. 0

	if( .not. bclosed ) then
	  write(6,*) 'line cannot be closed: ',n
	  write(6,*) (ipext(kvol(i)),i=1,n)
	  write(6,*) 'Please note that first and last node'
	  write(6,*) 'must be the same (no closure needed)'
	  write(6,*) 'or must be on the same boundary line.'
	  stop 'error stop volinf: no closure possible'
	else
	  if( bdebug ) then
	    write(6,*) '0 line has been closed: ',nnodes
	    write(6,'(i1,3i10)') 3,nv,0,nnodes
	    write(6,*) (ipext(kline(i)),i=1,nnodes)
	  end if
	end if

! collect elements in volume

	ntot = 0

	if( bdebug ) write(6,*) ntot,ivolm,nnodes,idim

	do ie=1,nel
	  call baric(ie,xc,yc)
	  if( inpoly(nnodes,xline,yline,xc,yc) ) then
	    ntot = ntot + 1
	    ivolm = ivolm + 1
	    if( ivolm .gt. idim ) goto 99
	    ivol(ivolm) = ie
	  end if
	end do

	if( ntot .eq. 0 ) then
	  stop 'error stop volinf: no elements found'
	end if

	if( bdebug ) write(6,*) '0 total elements: ',ntot

	ivolm = ivolm + 1
	if( ivolm .gt. idim ) goto 99
	ivol(ivolm) = 0

!	write(6,*) ivolm
!	write(6,*) (ivol(ie),ie=1,ivolm)

	return
   99	continue
	write(6,*) 'idim: ',idim
	write(6,*) 'Number of elements in volumes is too big.'
	stop 'error stop volinf: dimension ivol'
	end

!******************************************************************

	subroutine wrvolinf(n,ielem,volav,areav,hav)

! computes average values...

	implicit none

	integer n
	integer ielem(1)
	real volav,areav,hav

	integer i,ie
	real volume,area

	real volele,areaele

	volume = 0.
	area = 0.

	do i=1,n
	  ie = ielem(i)
	  volume = volume + volele(ie,0)
	  area = area + areaele(ie)
	end do

	volav = volume
	areav = area
	hav = volume / area

	end

!******************************************************************

	subroutine linecl(n,kvol,nnodes,kline,xline,yline)

! close line -> nodes will be unique (first & last are different)

	use mod_geom
	use basin

	implicit none

	integer n		!size of kvol
	integer kvol(1)
	integer nnodes
	integer kline(1)
	real xline(1), yline(1)

	integer i,k
	integer kfirst,kstart,knext

! copy original line

	do i=1,n
	  k = kvol(i)
	  kline(i) = k
	  xline(i) = xgv(k)
	  yline(i) = ygv(k)
	end do

	nnodes = n - 1			!make unique

	if( kvol(1) .eq. kvol(n) ) return

! not closed -> try to close it

	nnodes = nnodes + 1
	kfirst = kline(1)
	kstart = kline(n)
	knext = kantv(1,kstart)

        do while( knext .ne. kfirst .and. knext .ne. kstart             &
     &                  .and. knext .gt. 0 )

	  nnodes = nnodes + 1

	  kline(nnodes) = knext
	  xline(nnodes) = xgv(knext)
	  yline(nnodes) = ygv(knext)

	  knext = kantv(1,knext)

	  if( knext .le. 0 ) stop 'error stop linecl: internal error'

	end do
	
	if( knext .le. 0 ) then		!not on boundary -> cannot close
	  nnodes = -1
	else if( knext .eq. kstart ) then	!different boundaries
	  nnodes = -1
	end if

	end

!******************************************************************

	subroutine wrtelem(iunit,n,ielems,ietype)

	use basin

	implicit none

	integer iunit
	integer n
	integer ielems(1)
	integer ietype

	integer k,i,ii,ie,nvert
	integer kn(10)
	real hmed
	real v3v(nkn)

	real depele

	v3v = 0.

	do i=1,n

	  ie = ielems(i)

	  hmed = depele(ie,0)

	  call nindex(ie,nvert,kn)
	  do ii=1,nvert
	    v3v(kn(ii)) = 1.
	  end do

          write(iunit,2000) 2,ipev(ie),ietype,nvert                     &
     &                          ,(ipv(kn(ii)),ii=1,nvert),hmed

	end do

	  
	do k=1,nkn
	  if( v3v(k) .ne. 0. ) then
	    write(iunit,1000) 1,ipv(k),0,xgv(k),ygv(k)
	  end if
	end do

	return
 1000	format(i1,2i10,2e12.4)
 2000	format(i1,6i10,e12.4)
	end

!******************************************************************

	subroutine wrtline(iunit,n,inodes,il,iltype)

	use basin

	implicit none

	integer iunit
	integer n
	integer inodes(1)
	integer il,iltype

	integer k,i
	integer istart

	istart = 1
	if( inodes(1) .eq. inodes(n) ) istart = 2

	do i=istart,n
	  k = inodes(i)
	  write(iunit,1000) 1,ipv(k),0,xgv(k),ygv(k)
	end do

	write(iunit,2000) 3,il,iltype,n
	write(iunit,2001) (ipv(inodes(i)),i=1,n)
	  
	return
 1000	format(i1,2i10,2e12.4)
 2000	format(i1,3i10)
 2001	format((10i7))
	end

!******************************************************************

