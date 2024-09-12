
!--------------------------------------------------------------------------
!
!    Copyright (C) 1988,1990,1994-1998,2001,2005,2009  Georg Umgiesser
!    Copyright (C) 2011-2013,2015-2019  Georg Umgiesser
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

! pre-processing routine
!
! revision log :
!
! 30.08.1988	ggu	(itief in common, no rdpara)
! 22.11.1988	ggu	(new version 3, itief passed as actual param)
! 24.11.1988	ggu	(sp13f., descrr...)
! 30.11.1988	ggu	(back to sp13u.)
! 31.07.1990	ggu	(open all files explicitly)
! 08.10.1994	ggu	(newly designed -> use subroutines)
! 09.10.1994	ggu	(read from .grd files)
! 16.03.1995	ggu	(double precision in clockw)
! 06.03.1996	ggu	renumber also iarv in renel
! 08.09.1997	ggu	introduce raux,neaux for compiler warnings
! 20.03.1998	ggu	automatic optimization of bandwidth introduced
! 08.05.1998	ggu	always process whole file (idepth = 0)
! 18.05.1998	ggu	always process depths elementwise
! 18.05.1998	ggu	dont ask for description anymore
! 17.10.2001	ggu	accept also grd files with some missing data
! 18.10.2005	ggu	some error messages slightly changed
! 06.04.2009	ggu	read param.h
! 24.04.2009	ggu	new call to rdgrd()
! 04.03.2011	ggu	new routine estimate_grade()
! 30.03.2011	ggu	new routine check_sidei(), text in optest()
! 15.07.2011	ggu	calls to ideffi substituted
! 15.11.2011	ggu	new routines for mixed depth (node and elem), hflag
! 09.03.2012	ggu	delete useless error messages, handle nkn/nel better
! 29.03.2012	ggu	use ngr1 to avoid too small dimension for ngr
! 04.10.2013	ggu	in optest better error handling
! 30.07.2015	ggu	vp renamed to shypre
! 18.12.2015	ggu	changed VERS_7_3_17
! 09.09.2016	ggu	changed VERS_7_5_17
! 09.05.2017	ggu	changed VERS_7_5_26
! 09.10.2017	ggu	changed VERS_7_5_33
! 14.11.2017	ggu	changed VERS_7_5_36
! 17.11.2017	ggu	implement output switches (quiet,silent,etc..)
! 24.01.2018	ggu	changed VERS_7_5_41
! 13.04.2018	ggu	accepts partition to write bas file with node partition
! 16.10.2018	ggu	changed VERS_7_5_50
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 21.05.2020    ggu     better handle copyright notice
! 13.07.2020    ggu     honor noopti flag, stack poisoning eliminated
! 06.05.2023    ggu     some enhancements and better error handeling
! 22.05.2023    ggu     locate is defined in module
! 12.09.2024    ggu     flush bas file before closing (intel bug)
!
! notes :
!
! could eliminate de-scrambling of iknot ->
!       no knscr()
!       pass ngr1 to bandop
!       change cmv,rosen
!
!**********************************************************
!**********************************************************
!**********************************************************
! overall call
!**********************************************************
!**********************************************************
!**********************************************************

	subroutine shypre_sub(grdfile,bquiet,bsilent,bauto,binfo &
     &				,bopti,bnepart,eps_area,grdpart)

	use basin
	use evgeom

	implicit none

        character*(*) grdfile
	logical bquiet,bsilent,bauto,binfo,bopti,bnepart
	real eps_area
	character*(*) grdpart

	character*80 name
	character*80 file
        character*80 descrg,descra,basnam
        logical bstop,bwrite,bdebug,bww,bstopall

        integer, save, allocatable :: ipdex(:)		!node sort index
        integer, save, allocatable :: iedex(:)		!element sort index
        integer, save, allocatable :: kphv(:)		!node rank
        integer, save, allocatable :: ierank(:)		!element rank
        real, save, allocatable :: hev(:)		!element depth
        real, save, allocatable :: hkv(:)		!node depth
        integer, save, allocatable :: ng(:)		!node grade

	integer ianz,itief
	integer ie,k
	integer nbout

	integer nk,ne,nl,nn
        integer nknh,nelh,nli,nco
	integer nknddi,nelddi
	integer ngr1
	integer nrec
	integer nlidim,nlndim
	integer nne,nnl
	real hflag

	integer idefbas,ichanm,ifileo

!--------------------------------------------------------
! initialization of parameters
!--------------------------------------------------------

	bstop = .false.
	bstopall = .false.
	bdebug = .false.

        dcorbas=0.
        dirnbas=0.
	sphebas = -1
        descrg=' '
        descrr=' '
        descra=' '
!
	nbout=2

	hflag = -999.

	itief=0		!0=read by element  1=read by node

!--------------------------------------------------------
! set parameters
!--------------------------------------------------------

	bwrite = .not. bquiet
	bww = .not. bsilent

	basnam = grdfile
	call delete_extension(basnam,'.grd')

!--------------------------------------------------------
! read grid
!--------------------------------------------------------

	if( bww ) write(6,*) 'grdfile: ',trim(grdfile)
	call grd_set_write(bww)
	call grd_read(grdfile)

	call grd_get_params(nk,ne,nl,nne,nnl)
	if( bwrite ) write(6,*) 'grid info: ',nk,ne,nl

	if( nk == 0 .or. ne == 0 ) then
	  write(6,*) 'nk,ne: ',nk,ne
	  stop 'error stop shypre: no nodes or elements in basin'
	end if

	call grd_to_basin
	call bas_check_spherical

!--------------------------------------------------------
! set up internal geometrical structure
!--------------------------------------------------------

        !call ev_init(ne)
        call set_ev

!--------------------------------------------------------
! allocate arrays
!--------------------------------------------------------

	allocate(ipdex(nk))
	allocate(iedex(ne))
	allocate(kphv(nk))
	allocate(hev(ne))
	allocate(hkv(nk))
	allocate(ierank(ne))
	allocate(ng(nk))

!--------------------------------------------------------
! handle depth
!--------------------------------------------------------

	call grd_get_depth(nk,ne,hkv,hev)
        call process_shypre_depth(bwrite,bww,nkn,nel,hkv,hev,hflag &
     &                          ,nknh,nelh,itief)

	if( binfo ) stop

!--------------------------------------------------------
! open files
!--------------------------------------------------------

	nbout=idefbas(basnam,'new')
        if(nbout.le.0) stop 'error stop shypre: error opening output file'

!--------------------------------------------------------
! start processing
!--------------------------------------------------------

        call handle_shypre_numbering(bwrite,bdebug,eps_area,nkn,nel,nen3v &
     &                          ,ipv,ipev,ipdex,iedex &
     &                          ,xgv,ygv)

!--------------------------------------------------------
! setup side index, compute mbw and ngr, and do optimization
!--------------------------------------------------------

        !call debug_out('befop ipv',nkn,ipv)
        !call debug_out('befop ipev',nel,ipev)
        !call debug_out_r('befop hev',nel,hev)
        !call debug_out_r('befop hkv',nkn,hkv)

        call handle_shypre_optimization(bwrite,bww,bdebug,bauto,bopti &
     &                          ,nkn,nel,nen3v,ng,ipv &
     &                          ,kphv,iarnv &
     &                          ,xgv,ygv,hkv,mbw,ngr)

        !call debug_out('opti kphv',nkn,kphv)
        !call debug_out_r('afterop hev',nel,hev)
        !call debug_out_r('afterop hkv',nkn,hkv)

!--------------------------------------------------------
! renumber elements
!--------------------------------------------------------

        call renel(bwrite,nel,nen3v,ipev,iarv,hev,ierank)

        !call debug_out('renel ierank',nkn,ierank)
        !call debug_out_r('renel hev',nel,hev)
        !call debug_out_r('renel hkv',nkn,hkv)

!--------------------------------------------------------
! handles depths
!--------------------------------------------------------

        call handle_shypre_depth(bwrite,bdebug &
     &                          ,nkn,nel,nknh,nelh,nen3v &
     &                          ,ipdex,iedex &
     &                          ,ipv,ipev &
     &                          ,itief,hflag,hkv,hev,hm3v)

        !call debug_out('depth ipdex',nkn,ipdex)
        !call debug_out('depth iedex',nkn,iedex)
        !call debug_out('depth ipv',nkn,ipv)
        !call debug_out('depth ipev',nkn,ipev)
        !call debug_out_r('depth hev',nel,hev)
        !call debug_out_r('depth hkv',nkn,hkv)

!--------------------------------------------------------
! process partitions
!--------------------------------------------------------

	call handle_shypre_partition(bwrite,bnepart,grdpart,nkn,nel &
     &				,kphv,ierank)

!--------------------------------------------------------
! write to file
!--------------------------------------------------------

	descrr = descrg

	if( bwrite ) write(6,*) ' ...writing file '
	call basin_write(nbout)
	flush(nbout)
	close(nbout)

	if( bww ) call bas_info

!--------------------------------------------------------
! end of routine
!--------------------------------------------------------

 	end

!********************************************************************
!********************************************************************
!********************************************************************
! high level calls
!********************************************************************
!********************************************************************
!********************************************************************

	subroutine handle_shypre_numbering(bwrite,bdebug,eps_area &
     &				,nkn,nel,nen3v &
     &				,ipv,ipev,ipdex,iedex &
     &				,xgv,ygv)

	use mod_sort

	implicit none

	logical bwrite,bdebug
	real eps_area
	integer nkn,nel
	integer nen3v(3,nel)
	integer ipv(nkn),ipdex(nkn)
	integer ipev(nel),iedex(nel)
	real xgv(nkn),ygv(nkn)

	logical bstop, bstopall
	integer nknddi,nelddi
	integer, allocatable :: iaux(:)

        bstop=.false.
        bstopall=.false.
        nknddi = nkn
        nelddi = nel
	allocate(iaux(nkn))

        call gtest(bdebug,'start numbering',nelddi,nkn,nel,nen3v)

        if( bwrite ) write(6,*) ' ...making external index'

        call sort(nkn,ipv,ipdex)
        call sort(nel,ipev,iedex)

        if( bwrite ) write(6,*) ' ...controlling uniqueness of node numbers'
        call uniqn(bwrite,nkn,ipv,ipdex,bstop)
        if(bstop) goto 99909

        if( bwrite ) write(6,*) ' ...controlling uniqueness of element numbers'
        call uniqn(bwrite,nel,ipev,iedex,bstop)
        if(bstop) goto 99910

        call gtest(bdebug,'end uniqn',nelddi,nkn,nel,nen3v)

        if( bwrite ) write(6,*) ' ...controlling uniqueness of elements'
        call uniqe(bwrite,nel,nen3v,ipev,bstop)
        if( bstop ) then
          bstop = .false.
          bstopall = .true.
          if( bwrite ) write(6,*) 'trying to eliminate non-unique elements'
          call uniqe(bwrite,nel,nen3v,ipev,bstop)
          if(bstop) goto 99915
        end if

        if( bwrite ) write(6,*) ' ...controlling node numbers'
        call needn(nkn,nel,nen3v,ipv,iaux,bstop)
        if(bstop) goto 99918

        call gtest(bdebug,'end sense',nelddi,nkn,nel,nen3v)

        if( bwrite ) write(6,*) ' ...testing sense of nodes in index'
        call clockw(nkn,nel,nen3v,ipev,xgv,ygv,bstop)
        if( bstop ) then
          bstop = .false.
          bstopall = .true.
          if( bwrite ) write(6,*) 'trying to invert sense of elements'
          call clockw(nkn,nel,nen3v,ipev,xgv,ygv,bstop)
          if(bstop) goto 99930
        end if

        if( bwrite ) write(6,*) ' ...testing area of elements'
        call areaz(bwrite,eps_area,nkn,nel,nen3v,ipev,xgv,ygv,bstop)
        if(bstop) goto 99931

        if( bstopall ) goto 99940

	return
99909   write(6,*)' (09) error: no unique definition of nodes'
	stop 'error stop handle_numbering: error'
99910   write(6,*)' (10) error: no unique definition of elements'
	stop 'error stop handle_numbering: error'
99915   write(6,*)' (15) error: no unique definition of elements'
	stop 'error stop handle_numbering: error'
99918   write(6,*)' (18) error: not used nodes present'
	stop 'error stop handle_numbering: error'
99930   write(6,*)' (30) error: elements with clockwise nodes'
	stop 'error stop handle_numbering: error'
99931   write(6,*)' (31) error: zero area elements'
	stop 'error stop handle_numbering: error'
99940   write(6,*)' (40) error: generic error'
	stop 'error stop handle_numbering: error'
	end

!********************************************************************

	subroutine process_shypre_depth(bwrite,bww,nkn,nel,hkv,hev,hflag &
     &				,nknh,nelh,itief)

! checks how depth can be processed

	implicit none

	logical bwrite,bww
	integer nkn,nel
	real hkv(nkn),hev(nel)
	real hflag
	integer nknh,nelh
	integer itief

	integer k

        nknh = 0
        do k=1,nkn
          if( hkv(k) .ne. hflag ) nknh = nknh + 1
        end do

        nelh = 0
        do k=1,nel
          if( hev(k) .ne. hflag ) nelh = nelh + 1
        end do

        if( bww ) then
          write(6,*) 'nkn,nel   : ',nkn,nel
          write(6,*) 'nknh,nelh : ',nknh,nelh
        end if

        if(nkn.le.0 .or. nel.le.0) then
          write(6,*) ' Nothing to process'
          goto 99999
        end if

        itief=0
        if( nel.eq.nelh .and. nkn.eq.nknh ) then
         if( bwrite ) then
          write(6,*) ' Can process depth node or elementwise.'
          write(6,*) ' ...depths are processed elementwise'
         end if
        else if(nel.eq.nelh) then
         if( bwrite ) then
          write(6,*) ' ...depths are processed elementwise'
         end if
        else if(nkn.eq.nknh) then
          itief=1
         if( bwrite ) then
          write(6,*) ' ...depths are processed nodewise'
         end if
        else if(nknh.eq.0.and.nelh.eq.0) then
         if( bwrite ) then
          write(6,*) ' No depth data read. Process anyway'
         end if
        else
          itief=2
         if( bwrite ) then
          write(6,*) '********************************************'
          write(6,*) '********************************************'
          write(6,*) ' Mixed data source for depth. Process anyway'
          write(6,*) '********************************************'
          write(6,*) '********************************************'
         end if
        end if

	return
99999   write(6,*)' (99) error: no grid in grd file'
	stop 'error stop process_depth: no grid'
	end

!********************************************************************

	subroutine handle_shypre_depth(bwrite,bdebug &
     &				,nkn,nel,nknh,nelh,nen3v &
     &				,ipdex,iedex &
     &				,ipv,ipev &
     &				,itief,hflag,hkv,hev,hm3v)

! handles renumbering of depth values

	implicit none

	logical bwrite,bdebug
	integer nkn,nel,nknh,nelh
	integer nen3v(3,nel)
	integer ipdex(nkn),iedex(nel)
	integer ipv(nkn),ipev(nel)
	integer itief
	real hflag
	real hkv(nkn),hev(nel)
	real hm3v(3,nel)

	logical bstop
	integer nmax
	integer, allocatable :: iaux(:)

	bstop = .false.
	nmax = max(nkn,nel)
	allocate(iaux(nmax))

        if( bwrite ) write(6,*) ' ...handling depths'

	hm3v = hflag

        if(itief.eq.0) then
          if( bwrite ) write(6,*) ' ...................elementwise'
          call helem(nel,nelh,iaux,iedex,ipev &
     &                          ,hm3v,hev,bstop)
        else if( itief .eq. 1 ) then
          if( bwrite ) write(6,*) ' ......................nodewise'
          call hnode(nkn,nel,nknh,nen3v,iaux,ipdex,ipv &
     &                      ,hm3v,hkv,bstop)
        else
          if( bwrite ) write(6,*) ' .............first elementwise'
          call helem(nel,nelh,iaux,iedex,ipev &
     &                          ,hm3v,hev,bstop)
          if( bwrite ) write(6,*) ' .................then nodewise'
          call hnode(nkn,nel,nknh,nen3v,iaux,ipdex,ipv &
     &                      ,hm3v,hkv,bstop)
        end if

        if( bstop ) goto 99932

        call check_hm3v(nel,hm3v,hflag)
        call gtest(bdebug,'check hm3v',nel,nkn,nel,nen3v)

	return
99932   write(6,*)' (32) error: handling depth'
	stop 'error stop handle_depth: error'
	end

!********************************************************************
!********************************************************************
!********************************************************************
! low level calls
!********************************************************************
!********************************************************************
!********************************************************************

        subroutine uniqn(bwrite,n,iv,index,bstop)

! controlls uniqueness of numbers in iv

        implicit none

	logical bwrite
        integer n
        logical bstop
        integer iv(n)
        integer index(n)

        integer k1,k2,i

        k1=iv(index(1))
        do i=2,n
          k2=iv(index(i))
          if(k1.eq.k2) then
            if( bwrite ) write(6,*)' warning : ',k1,' not unique'
            bstop=.true.
          end if
          k1=k2
	end do

        return
        end

!*****************************************************************

        subroutine uniqe(bwrite,nel,nen3v,ipev,bstop)

! controlls uniqueness of elements (shell for equale)

	use mod_sort

        implicit none

	logical bwrite
        integer nel
        logical bstop
        integer nen3v(3,nel)
        integer ipev(nel)

        integer ie,ii,ka,km,k1,k,ne,i
        integer ie1,ie2
        integer isum,iold,inew,itest
        integer iaux(nel),index(nel),elist(2,nel)

        isum=0

        do ie=1,nel
          ka=0
          km=1
          k1=1
          do ii=1,3
            k=nen3v(ii,ie)
            ka=ka+k
            km=km*(mod(k,97)+1)
            k1=k1*(mod(k+1,97)+1)
          end do
          iaux(ie)=ka+km+k1
        end do

        call sort(nel,iaux,index)

	iold = 0
	itest = 0
	do ie=1,nel
	  inew = iaux(index(ie))
	  if( inew == iold ) itest = itest + 1
	  iold = inew
	end do

	if( bwrite ) write(6,*) 'testing uniqueness of elements: ',itest

	ne = 0
	elist = 0

        ie1=1
        do while(ie1.le.nel)
          ie2=ie1+1
          do
	    if( ie2 > nel ) exit
	    if( iaux(index(ie1)) /= iaux(index(ie2)) ) exit
            ie2=ie2+1
          end do
          if(ie2-ie1.gt.1) then
            isum=isum+ie2-ie1-1
	    call equale(nel,ie1,ie2-1,nen3v,index,ipev,ne,elist)
          end if
          ie1=ie2
        end do

	if( ne > 0 ) then
	  bstop = .true.
	  if( bwrite ) then
	    write(6,*) 'non-unique elements present: ',ne
	  end if
	  call elime(ne,elist,nel,nen3v,ipev)
	end if

        !write(6,*) 'uniqe (isum) : ',isum

        return
        end

!*****************************************************************

	subroutine elime(ne,elist,nel,nen3v,ipev)

	implicit none

	integer ne
	integer elist(2,ne)
	integer nel
	integer nen3v(3,nel)
	integer ipev(nel)

	integer i,ie

	do i=1,ne
	  ie = elist(2,i)
	  nen3v(:,ie) = nen3v(:,nel)
	  ipev(ie) = ipev(nel)
	  nel = nel - 1
	end do

	end

!*****************************************************************

	subroutine equale(nel,ip1,ip2,nen3v,index,ipev,ne,elist)

! controlls uniqueness of listed elements

        implicit none

        integer nel
        integer ip1,ip2
        integer nen3v(3,nel)
        integer index(nel)
        integer ipev(nel)
	integer ne
        integer elist(2,nel)

        integer ip,i,ie1,ie2
        integer i1,i2,i3,kn1,kn2,kn3

        do ip=ip1,ip2
          ie1=index(ip)
          kn1=nen3v(1,ie1)
          do i=ip+1,ip2
            ie2=index(i)
            do i1=1,3
              if(nen3v(i1,ie2).eq.kn1) then
                i2=mod(i1,3)+1
                i3=mod(i2,3)+1
                kn2=nen3v(2,ie1)
                kn3=nen3v(3,ie1)
                if(nen3v(i2,ie2).eq.kn2 .and. &
     &                  nen3v(i3,ie2).eq.kn3) then
                  write(6,*)' error: element',ipev(ie1) &
     &                  ,'  and',ipev(ie2),'  are identical'
		  ne = ne + 1
		  elist(1,ne) = ie1
		  elist(2,ne) = ie2
                end if
              end if
            end do
          end do
	end do

        return
        end

!*****************************************************************

        subroutine needn(nkn,nel,nen3v,ipv,iaux,bstop)

! controlls if all nodes are needed (use nen3v as one-dim array)

        implicit none

        integer nkn,nel
        logical bstop
        integer nen3v(3,nel)
        integer ipv(nkn)
        integer iaux(nkn)

        integer k,ie,ii

        do k=1,nkn
          iaux(k)=0
        end do

        do ie=1,nel
          do ii=1,3
            iaux(nen3v(ii,ie))=1
          end do
        end do

        do k=1,nkn
          if(iaux(k).eq.0) then
            write (6,*)' error: node ',ipv(k),' not needed'
            bstop=.true.
          end if
        end do

        return
        end

!*****************************************************************

        subroutine areaz(bwrite,eps_area,nkn,nel,nen3v,ipev,xgv,ygv,bstop)

! test for zero area elements

	use mod_sort

        implicit none

	logical bwrite
	real eps_area
        integer nkn,nel
        logical bstop
        integer nen3v(3,nel)
        integer ipev(nel)
        real xgv(nkn),ygv(nkn)

	logical bwarn
        integer ie,ii,ii2,ii3,k1,k2,k3,ke,i
	integer nsmall,nzero
	integer, parameter :: nsmax = 20
        real a,x1,x2,x3,y1,y2,y3,p
        real area(nel),amin,amax
	real eps

	real area_elem

	amin = 1.e+30
	amax = -1.e+30
	eps = eps_area
	if( .not. bwrite ) eps = 0.	!do not write small elements
	bwarn = .false.

	do ie=1,nel
	  ii = 1
          ii2=mod(ii,3)+1
          ii3=mod(ii2,3)+1
          k1=nen3v(ii,ie)
          k2=nen3v(ii2,ie)
          k3=nen3v(ii3,ie)
	  x1=xgv(k1)
	  y1=ygv(k1)
	  x2=xgv(k2)
	  y2=ygv(k2)
	  x3=xgv(k3)
	  y3=ygv(k3)
	  a = area_elem(ie)			!honors spherical
	  a = abs(a)
	  amin = min(amin,a)
	  amax = max(amax,a)
	  area(ie) = a
	end do

	nsmall = 0
	nzero = 0

	do ie=1,nel
	  a = area(ie)
	  p = a / amax
	  ke = ipev(ie)
          if(a.eq.0.) then
	    nzero = nzero + 1
	    if( nzero < nsmax ) then
              write(6,*)' error: area in element ',ipev(ie) &
     &                     ,' is zero'
	    end if
	    if( nzero == nsmax ) then
	      write(6,*) ' ... more zero elements not showing'
	    end if
            bstop=.true.
          else if(p.lt.eps) then
	    nsmall = nsmall + 1
	    if( nsmall < nsmax ) then
              write(6,1000)' warning: element ',ipev(ie) &
     &                     ,' seems too small:',a,p
	    end if
	    if( nsmall == nsmax ) then
	      write(6,*) ' ... more small elements not showing'
	    end if
            bwarn=.true.
          end if
	end do

	if( bwrite .and. bwarn ) then
	  if( nsmall > 0 ) then
	    write(6,*) ' total number of small elements: ',nsmall
	  end if
	  write(6,*) ' area amin,amax: = ',amin,amax,amin/amax
	end if

        return
 1000	format(1x,a,i10,a,2g14.5)
        end

!*****************************************************************

        subroutine clockw(nkn,nel,nen3v,ipev,xgv,ygv,bstop)

! test for anti-clockwise sense of nodes

        implicit none

        integer nkn,nel
        logical bstop
        integer nen3v(3,nel)
        integer ipev(nel)
        real xgv(nkn),ygv(nkn)

        integer ie,ii,iii,k1,k2,ieaux
        double precision a,x1,x2,y1,y2	!new 16.3.95

	do ie=1,nel
          a=0.
          do ii=1,3
            iii=mod(ii,3)+1
            k1=nen3v(ii,ie)
            k2=nen3v(iii,ie)
	    x1=xgv(k1)
	    x2=xgv(k2)
	    y1=ygv(k1)
	    y2=ygv(k2)
            a=a+x1*y2-x2*y1
          end do
          if(a.lt.0.) then
            write(6,*)' error: nodes in element ',ipev(ie) &
     &                     ,' are in clockwise sense... adjusting'
            bstop=.true.
	    ieaux = nen3v(1,ie)
	    nen3v(1,ie) = nen3v(2,ie)
	    nen3v(2,ie) = ieaux
          end if
	end do

        return
        end

!*****************************************************************

        subroutine estimate_grade(nkn,nel,nen3v,ng,ngr)

! estimates ngr

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        integer ng(nkn)

        integer i,ii,ie,k,ngr

	do i=1,nkn
          ng(i)=0
	end do

	do ie=1,nel
          do ii=1,3
	    k = nen3v(ii,ie)
	    ng(k) = ng(k) + 1
	  end do
	end do

	ngr = 0
	do i=1,nkn
          ngr = max(ngr,ng(i))
	end do
	ngr = ngr + 1		!account for boundary nodes

	end

!*****************************************************************

        subroutine check_sidei(nkn,nel,nen3v,ipv,ng,ngr1,bstop)

! set up side index and find grade

        implicit none

        integer nkn,nel,ngr
        integer ngr1
        integer nen3v(3,nel)
        integer ipv(nkn)
        integer ng(nkn)
        !integer iknot(ngr1,nkn)
	logical bstop

        integer iaux(nkn)

        integer ii,ie,k
	integer igr,iel

	bstop = .false.
	iaux = 0

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    iaux(k) = iaux(k) + 1
	  end do
	end do

	do k=1,nkn
	  igr = ng(k)
	  iel = iaux(k)
	  if( igr .gt. iel + 1 .or. igr .lt. iel ) then
	    write(6,*) ' error: irregular connection of node: ' &
     &				,ipv(k),igr,iel
	    bstop = .true.
	  end if
	end do

	if( bstop ) then
	  write(6,*) ' ...please check nodes for connections'
	end if

	end

!*****************************************************************

        subroutine sidei(nkn,nel,nen3v,ng,iknot,ngr1,ngr,bstop)

! set up side index and find grade ngr

        implicit none

        integer nkn,nel,ngr
        integer ngr1
        integer nen3v(3,nel)
        integer ng(nkn)
        integer iknot(ngr1,nkn)
	logical bstop

        integer i,ii,iii,ie,k1,k2

	do i=1,nkn
	  do ii=1,ngr1
	    iknot(ii,i)=0
	  end do
          ng(i)=0
	end do

	do ie=1,nel
          do ii=1,3
            iii=mod(ii,3)+1
            k1=nen3v(ii,ie)
            k2=nen3v(iii,ie)
            do i=1,ng(k1)
              if(iknot(i,k1).eq.k2) goto 1    !already there
            end do
            ng(k1)=ng(k1)+1        !side not yet in index
	    if(ng(k1).gt.ngr1) goto 99
            iknot(ng(k1),k1)=k2
            ng(k2)=ng(k2)+1        !has to be added to index
	    if(ng(k2).gt.ngr1) goto 99
            iknot(ng(k2),k2)=k1
    1       continue
          end do
	end do

	ngr=0
	do i=1,nkn
          if(ng(i).gt.ngr) ngr=ng(i)
	end do

        return
   99	continue
	bstop = .true.
	write(6,*) ' error: dimension of ngr1 is too low: ',ngr1
        end

!**********************************************************

        subroutine bandw(nel,nen3v,mbw)

! determine bandwidth mbw

        implicit none

        integer nel,mbw
        integer nen3v(3,nel)

        integer ie,ii,iii,k,kk
        integer mh,mm

	mh=0

	do ie=1,nel
          do ii=1,3
            k=nen3v(ii,ie)
            do iii=ii+1,3
              kk=nen3v(iii,ie)
              mm=iabs(kk-k)
              if(mm.gt.mh) mh=mm
            end do
          end do
	end do

	mbw=mh

	return
	end

!**********************************************************

	subroutine ininum(nkn,iphv,kphv)

! initializes numbering of nodes

	implicit none

	integer nkn
	integer iphv(nkn), kphv(nkn)

	integer i

	do i=1,nkn
	  iphv(i) = i
	  kphv(i) = i
	end do

	end

!**********************************************************

        subroutine icopy(n,iv,irank)

! copy one array to itself exchanging elements as in irank

        implicit none

        integer n
        integer iv(n)
        integer irank(n)

        integer i
	integer, allocatable :: iauxv(:)

	allocate(iauxv(n))

        do i=1,n
          iauxv(i)=iv(i)
        end do

        do i=1,n
          iv(irank(i))=iauxv(i)
        end do

        return
        end

!**********************************************************

        subroutine rcopy(n,rv,irank)

! copy one array to itself exchanging elements as in irank

        implicit none

        integer n
        real rv(n)
        integer irank(n)

        integer i
	real, allocatable :: rauxv(:)

	allocate(rauxv(n))

        do i=1,n
          rauxv(i)=rv(i)
        end do

        do i=1,n
          rv(irank(i))=rauxv(i)
        end do

        return
        end

!**********************************************************

        subroutine renel(bwrite,nel,nen3v,ipev,iarv,hev,ierank)

! renumbering of elements

! we construct iedex newly and use iaux,iedex as aux arrays
! iedex is also used as a real aux array for rcopy	-> changed
! neaux is probably he3v (real) used as an aux array	-> changed

	use mod_sort

        implicit none

	logical bwrite
        integer nel
        integer nen3v(3,nel)
        integer ipev(nel)
	integer iarv(nel)
        real hev(nel)
	integer ierank(nel)

        integer ie,ii
        integer iedex(nel)
        integer neaux(3,nel)
        integer ival(nel)

        if( bwrite ) write(6,*) ' ...renumbering elements'

        do ie=1,nel
          ival(ie)=min(nen3v(1,ie),nen3v(2,ie),nen3v(3,ie))
	end do

        call sort(nel,ival,iedex)    !iedex is the index table
        call rank(nel,iedex,ierank)   !ierank is the rank table

        call icopy(nel,ipev,ierank)
        call icopy(nel,iarv,ierank)
        call rcopy(nel,hev,ierank)

        do ie=1,nel
          do ii=1,3
            neaux(ii,ie)=nen3v(ii,ie)
          end do
        end do

        do ie=1,nel
          do ii=1,3
            nen3v(ii,ierank(ie))=neaux(ii,ie)
          end do
        end do

        return
        end

!**********************************************************

        subroutine rank(n,index,irank)

! builds rank table from index table

        implicit none

        integer n
        integer index(n),irank(n)

        integer i

        do i=1,n
          irank(index(i))=i
        end do

        return
        end

!**********************************************************

	subroutine check_hm3v(nel,hm3v,hflag)

	implicit none

	integer nel
	real hm3v(3,nel)
	real hflag

	logical bmiss
	integer ie,ii,iflag
	real h

	iflag = 0

	do ie=1,nel
	  bmiss = .false.
	  do ii=1,3
	    h = hm3v(ii,ie)
	    if( h .eq. hflag ) then
	      iflag = iflag + 1
	      bmiss = .true.
	    end if
	  end do
	  if( bmiss ) write(6,*) ie,(hm3v(ii,ie),ii=1,3)
	end do

	if( iflag .gt. 0 ) then
	  write(6,*) '*******************************************'
	  write(6,*) '*******************************************'
	  write(6,*) '*******************************************'
	  write(6,*) 'flag found in depth: ',iflag,iflag/3
	  write(6,*) '*******************************************'
	  write(6,*) '*******************************************'
	  write(6,*) '*******************************************'
	end if

	end

!**********************************************************

        subroutine helem(nel,nhd,iaux,iedex &
     &                    ,ipev,hm3v,hev,bstop)

! depth by elements

	use mod_sort

        implicit none

        integer nel,nhd
        logical bstop
        integer iaux(nel)
        integer iedex(nel)
        integer ipev(nel)
        real hm3v(3,nel),hev(nel)

        integer ie,i,iel,ii
        integer, allocatable :: iphev(:)

        iaux = 0
	allocate(iphev(nel))
	iphev = ipev

        call sort(nel,ipev,iedex)

! it is : ipev(ie) == iphev(i)

        !do i=1,nhd
        do i=1,nel
          iel=iphev(i)
          ie=locate(nel,ipev,iedex,iel)
          if(ie.le.0) then
            write(6,*)' warning : element',iel,' not found'
            bstop=.true.
          else
            if(iaux(ie).ne.0) then
              write(6,*)' for element ',ipev(ie) &
     &                        ,' depth data not unique'
              bstop=.true.
            else
              iaux(ie)=i
            end if
          end if
        end do

        do  ie=1,nel
          if(iaux(ie).eq.0) then
            write(6,*)' for element ',ipev(ie) &
     &                        ,' no depth data found'
            bstop=.true.
	  else
            do ii=1,3
              hm3v(ii,ie)=hev(iaux(ie))
            end do
	  end if
	end do

        return
        end

!**********************************************************

        subroutine hnode(nkn,nel,nhd,nen3v,iaux,ipdex,ipv &
     &                      ,hm3v,hkv,bstop)

! depth by nodes

	use mod_sort

        implicit none

        integer nkn,nel,nhd
        logical bstop
        integer iaux(nkn)
        integer ipdex(nkn)
        integer ipv(nkn)
        integer nen3v(3,nel)
        real hm3v(3,nel),hkv(nkn)

        integer ie,ii,i,k,kn
	real h,hflag
	integer, allocatable :: iphv(:)

	hflag = -999.
        iaux = 0
	allocate(iphv(nkn))
	iphv = ipv

        call sort(nkn,ipv,ipdex)

        do i=1,nhd
          kn=iphv(i)
          k=locate(nkn,ipv,ipdex,kn)
          if(k.le.0) then
            write(6,*)' warning : node',kn,' not found'
            bstop=.true.
          else
            if(iaux(k).ne.0) then
              write(6,*)' for node ',ipv(k) &
     &                        ,' depth data not unique'
              bstop=.true.
            else
              iaux(k)=i
            end if
          end if
        end do

        do  k=1,nkn
          if(iaux(k).eq.0) then
            write(6,*)' for node ',ipv(k) &
     &                        ,' no depth data found'
            bstop=.true.
	  end if
	end do

        do ie=1,nel
          do ii=1,3
            kn=nen3v(ii,ie)
            h = hkv(iaux(kn))
            if( iaux(kn) .gt. 0 .and. h .ne. hflag ) then
              hm3v(ii,ie) = h
            end if
          end do
	end do

        return
        end

!**********************************************************

	subroutine gtest(bdebug,text,nelddi,nkn,nel,nen3v)

	implicit none

	logical bdebug
	character*(*) text
	integer nelddi,nkn,nel
	integer nen3v(3,nelddi)

	integer ie,ii,k
	integer iii
	logical berror
	
	!if( .not. bdebug ) return

	berror = .false.

	if( bdebug ) then
	  write(6,*) 'testing ... ',text
	  write(6,*) nelddi,nkn,nel
	end if

	do ie=1,nel
	  do ii=1,3
	    iii = mod(ii,3) + 1
	    if( nen3v(ii,ie) .eq. nen3v(iii,ie) ) then
	      if( .not. berror ) then
	  	write(6,*) 'internal error found: (there might be more)'
		write(6,*) ie,(nen3v(iii,ie),iii=1,3)
		berror = .true.
	      end if
	    end if
	  end do
	end do

	if( berror ) then
	  write(6,*) 'testing has found errors...'
	  write(6,*) trim(text)
	  stop 'error stop gtest: error'
	end if

	end

!**********************************************************

        subroutine debug_out(text,n,ia)

        implicit none

        character*(*) text
        integer n
        integer ia(n)

        integer iu,nn,i

        iu = 589
        nn = 10

        write(iu,*) trim(text)
        do i=1,n,n/nn
          write(iu,*) ia(i)
        end do

        end

!**********************************************************

        subroutine debug_out_r(text,n,ra)

        implicit none

        character*(*) text
        integer n
        real ra(n)

        integer iu,nn,i

        iu = 589
        nn = 10

        write(iu,*) trim(text)
        do i=1,n,n/nn
          write(iu,*) ra(i)
        end do

        end

!**********************************************************

