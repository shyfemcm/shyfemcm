
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
!
! notes :
!
! could eliminate scrambling of iknot ->
!       no knscr()
!       pass ngr1 to bandop
!       change cmv,rosen
!
!**********************************************************

!==========================================================
	module mod_shypre
!==========================================================

	logical, save :: bopti
	logical, save :: bauto
	logical, save :: bnoopti
	logical, save :: bmanual

	logical, save :: binfo
	logical, save :: bquiet
	logical, save :: bsilent

	real, save :: eps_area = 1.e-4

!==========================================================
	end module mod_shypre
!==========================================================

!**********************************************************

        program shypre

	use basin
	use mod_shypre
	use mod_sort
	use evgeom

	implicit none

	character*80 name
	character*80 file
	character*80 errfil
	character*80 errtex
        character*80 descrg,descra,basnam,grdfile
        logical bstop,bwrite,bdebug,bww,bstopall

        integer, save, allocatable :: ipdex(:), iedex(:)
        integer, save, allocatable :: iphv(:), kphv(:)
        integer, save, allocatable :: iphev(:), iaux(:), ierank(:)
	integer, save, allocatable :: neaux(:,:)
        real, save, allocatable :: raux(:)
        real, save, allocatable :: hev(:)
        real, save, allocatable :: hkv(:)
        integer, save, allocatable :: ng(:)
        integer, save, allocatable :: kvert(:,:)

	integer, save, allocatable :: iknot(:,:)

	integer ianz,idepth,itief
	integer ie,k
	integer nat,net,nb1,nb2,nb3,nb4,nb9

	integer nk,ne,nl,nn
        integer nknh,nelh,nli,nco
	integer nknddi,nelddi
	integer ngr1
	integer nrec
	integer nlidim,nlndim
	integer nne,nnl
	real hflag

	integer idefbas,ichanm,ifileo

	bstop = .false.
	bstopall = .false.
	bdebug = .false.

        dcorbas=0.
        dirnbas=0.
        descrg=' '
        descrr=' '
        descra=' '
!
	nb1=1
	nb2=2
	nb3=3
	nb4=4
	nb9=9
	net=5
	nat=6

	hflag = -999.

	itief=0		!0=read by element  1=read by node

!--------------------------------------------------------
! get name of basin
!--------------------------------------------------------

	call shypre_init(grdfile)

!--------------------------------------------------------
! set parameters
!--------------------------------------------------------

	bwrite = .not. bquiet
	bww = .not. bsilent

	!call pardef
	basnam = grdfile
	call delete_extension(basnam,'.grd')
	!call putfnm('basnam',grdfile)

!--------------------------------------------------------
! always process whole file
!--------------------------------------------------------

	idepth = 0

!--------------------------------------------------------
! read grid
!--------------------------------------------------------

	if( bww ) write(6,*) 'grdfile: ',trim(grdfile)
	call grd_read(grdfile)

	call grd_get_params(nk,ne,nl,nne,nnl)
	if( bwrite ) write(6,*) 'grid info: ',nk,ne,nl

	if( nk == 0 .or. ne == 0 ) then
	  write(6,*) 'nk,ne: ',nk,ne
	  stop 'error stop shypre: no nodes or elements in basin'
	end if

	call grd_to_basin

!--------------------------------------------------------
! set up internal geometrical structure
!--------------------------------------------------------

        !call ev_init(ne)
        call set_ev

!--------------------------------------------------------
! allocate arrays
!--------------------------------------------------------

	allocate(ipdex(nk), iedex(ne))
	allocate(iphv(nk), kphv(nk))
	allocate(iphev(ne), iaux(ne))
	allocate(neaux(3,ne))
	allocate(raux(ne))
	allocate(hev(ne))
	allocate(hkv(nk))
	allocate(ierank(ne))
	allocate(ng(nk))
	allocate(kvert(2,nk))

!--------------------------------------------------------
! handle depth
!--------------------------------------------------------

	call grd_get_depth(nk,ne,hkv,hev)

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
          !write(6,*) 'nli,nco   : ',nli,nco
	end if

        if(nkn.le.0 .or. nel.le.0) then
          write(6,*) ' Nothing to process'
	  goto 99999
        end if

	itief=0
	if( nel.eq.nelh .and. nkn.eq.nknh ) then
	 if( bwrite ) then
	  write(nat,*) ' Can process depth node or elementwise.'
          write(nat,*) ' ...depths are processed elementwise'
	 end if
	else if(nel.eq.nelh) then
	 if( bwrite ) then
          write(nat,*) ' ...depths are processed elementwise'
	 end if
	else if(nkn.eq.nknh) then
	  itief=1
	 if( bwrite ) then
          write(nat,*) ' ...depths are processed nodewise'
	 end if
	else if(nknh.eq.0.and.nelh.eq.0) then
	 if( bwrite ) then
	  write(nat,*) ' No depth data read. Process anyway'
	 end if
	else
	  itief=2
	 if( bwrite ) then
	  write(nat,*) '********************************************'
	  write(nat,*) '********************************************'
	  write(nat,*) ' Mixed data source for depth. Process anyway'
	  write(nat,*) '********************************************'
	  write(nat,*) '********************************************'
	 end if
	end if

	if( binfo ) stop

!--------------------------------------------------------
! open files
!--------------------------------------------------------

	nb2=idefbas(basnam,'new')
        if(nb2.le.0) stop

!--------------------------------------------------------
! start processing
!--------------------------------------------------------

	nknddi = nkn
	nelddi = nel

	call ketest(nel,nen3v)
	if( bdebug ) then
	call gtest_0('end read',nelddi,nkn,nel,nen3v)
	end if

        bstop=.false.

	if( bwrite ) then
        write(nat,*) ' ...making external index'
	end if

        call sort(nkn,ipv,ipdex)
        call sort(nel,ipev,iedex)

	if( bwrite ) then
        write(nat,*) ' ...controlling uniqueness of node numbers'
	end if

        call uniqn(nkn,ipv,ipdex,bstop)
	if(bstop) goto 99909

	if( bwrite ) then
        write(nat,*) ' ...controlling uniqueness of element numbers'
	end if

        call uniqn(nel,ipev,iedex,bstop)

	if( bdebug ) then
	call gtest_0('end uniqn',nelddi,nkn,nel,nen3v)
	end if

	if( bwrite ) then
        write(nat,*) ' ...controlling uniqueness of elements'
	end if

        call uniqe(nel,nen3v,ipev,bstop)
	if( bstop ) then
	  bstop = .false.
	  bstopall = .true.
	  if( bwrite ) then
	    write(6,*) 'trying to eliminate non-unique eliments'
	  end if
          call uniqe(nel,nen3v,ipev,bstop)
          if(bstop) goto 99915
	end if

	!if( bwrite ) then
	!write(nat,*) ' ...changing extern with intern node numbers'
	!end if

        !call chexin(nkn,nel,nen3v,ipv,ipdex,bstop)
	!if(bstop) goto 99920

	if( bwrite ) then
        write(nat,*) ' ...controlling node numbers'
	end if

        call needn(nkn,nel,nen3v,ipv,iaux,bstop)
	if(bstop) goto 99918

	if( bwrite ) then
	write(nat,*) ' ...testing sense of nodes in index'
	end if

	call ketest(nel,nen3v)
	if( bdebug ) then
	call gtest_0('end sense',nelddi,nkn,nel,nen3v)
	end if

        call clockw(nkn,nel,nen3v,ipev,xgv,ygv,bstop)
	if( bstop ) then
	  bstop = .false.
	  bstopall = .true.
          call clockw(nkn,nel,nen3v,ipev,xgv,ygv,bstop)
	  if(bstop) goto 99930
	end if

	if( bwrite ) then
	write(nat,*) ' ...testing area of elements'
	end if

        call areaz(nkn,nel,nen3v,ipev,xgv,ygv,bstop)
	if(bstop) goto 99931

	if( bwrite ) then
	write(nat,*) ' ...setting up side index'
	end if

        call estimate_grade(nkn,nel,nen3v,ng,ngr1)
	allocate(iknot(ngr1,nk))

        call sidei(nkn,nel,nen3v,ng,iknot,ngr1,ngr,bstop)
	if(bstop) goto 99934
        call check_sidei_0(nkn,nel,nen3v,ipv,ng,iknot,ngr1,bstop)
	if(bstop) goto 99935
	call knscr(nkn,ngr,ngr1,iknot)

	if( bww ) write(nat,*) 'Maximum grade of nodes is ',ngr

	call debug_out('befop ipv',nkn,ipv)
	call debug_out('befop ipev',nel,ipev)
        call debug_out_r('befop hev',nel,hev)
        call debug_out_r('befop hkv',nkn,hkv)

!--------------------------------------------------------
! bandwidth optimization
!--------------------------------------------------------

	if( bdebug ) then
	call gtest_0('bandwidth',nelddi,nkn,nel,nen3v)
	end if

        if( bwrite ) write(nat,*) ' ...optimizing band width'

        call bandw(nel,nen3v,mbw)

	if( bdebug ) then
	call gtest_0('bandwidth 1',nelddi,nkn,nel,nen3v)
	end if
	if( bww ) write(nat,*) 'Bandwidth is ',mbw

        call bandop(nkn,ngr,ipv,iphv,kphv,ng,iknot,kvert &
     &			,bopti,bauto,bww)

	!call debug_out('bandop iphv',nkn,iphv)
	!call debug_out('bandop kphv',nkn,kphv)

	if( bdebug ) then
	call gtest_0('bandwidth 2',nelddi,nkn,nel,nen3v)
	end if

        if(bopti) then
          call bandex(nkn,nel,nen3v,kphv,ipv,iarnv,xgv,ygv,hkv) !kphv is n-rank
	  if( bdebug ) then
	    call gtest_0('bandwidth 3',nelddi,nkn,nel,nen3v)
	  end if
          call bandw(nel,nen3v,mbw)
	  if( bdebug ) then
	    call gtest_0('bandwidth 4',nelddi,nkn,nel,nen3v)
	  end if
	  if( bww ) write(nat,*) 'Optimized bandwidth is ',mbw
	end if

	call debug_out('opti kphv',nkn,kphv)
        call debug_out_r('afterop hev',nel,hev)
        call debug_out_r('afterop hkv',nkn,hkv)


!--------------------------------------------------------
! renumber elements
!--------------------------------------------------------

	call ketest(nel,nen3v)
	if( bdebug ) then
	call gtest_0('end bandwidth',nelddi,nkn,nel,nen3v)
	end if

	if( bwrite ) write(nat,*) ' ...renumbering elements'

        call renel(nel,nen3v,ipev,iarv,hev,ierank)	!ierank is element rank

	call debug_out('renel ierank',nkn,ierank)
        call debug_out_r('renel hev',nel,hev)
        call debug_out_r('renel hkv',nkn,hkv)

!--------------------------------------------------------
! save pointers for depth
!--------------------------------------------------------

        if( bwrite ) write(nat,*) ' ...saving pointers'

	do k=1,nkn
	  iphv(k)=ipv(k)
	end do

	do ie=1,nel
	  iphev(ie)=ipev(ie)
	end do

	call ketest(nel,nen3v)
	if( bdebug ) then
	call gtest_0('write',nelddi,nkn,nel,nen3v)
	end if

!--------------------------------------------------------
! process depths
!--------------------------------------------------------

        if( bwrite ) write(nat,*) ' ...processing depths'

	call init_hm3v(nel,hm3v,hflag)

	if(itief.eq.0) then
          if( bwrite ) write(nat,*) ' ...................elementwise'
          call helem(nel,nelh,iphev,iaux,iedex,ipev &
     &				,hm3v,hev,bstop)
	else if( itief .eq. 1 ) then
          if( bwrite ) write(nat,*) ' ......................nodewise'
          call hnode(nkn,nel,nknh,nen3v,iphv,iaux,ipdex,ipv &
     &                      ,hm3v,hkv,bstop)
	else
          if( bwrite ) write(nat,*) ' .............first elementwise'
          call helem(nel,nelh,iphev,iaux,iedex,ipev &
     &				,hm3v,hev,bstop)
          if( bwrite ) write(nat,*) ' .................then nodewise'
          call hnode(nkn,nel,nknh,nen3v,iphv,iaux,ipdex,ipv &
     &                      ,hm3v,hkv,bstop)
	end if

        call debug_out('depth ipdex',nkn,ipdex)
        call debug_out('depth iedex',nkn,iedex)
        call debug_out('depth ipv',nkn,ipv)
        call debug_out('depth ipev',nkn,ipev)
        call debug_out_r('depth hev',nel,hev)
        call debug_out_r('depth hkv',nkn,hkv)


	if( bstop ) goto 99932

	call check_hm3v(nel,hm3v,hflag)
	call ketest(nel,nen3v)

	descrr=descrg

!--------------------------------------------------------
! process partitions
!--------------------------------------------------------

	if( bwrite ) write(6,*) 'handle_partition: ',nkn,nel
	call handle_partition(nkn,nel,kphv,ierank)

!--------------------------------------------------------
! final error check
!--------------------------------------------------------

	if( bstopall ) goto 99940

!--------------------------------------------------------
! write to file
!--------------------------------------------------------

	if( bwrite ) write(nat,*) ' ...writing file '

	call basin_write(nb2)

	close(nb2)

	if( bww ) call bas_info

!--------------------------------------------------------
! end of routine
!--------------------------------------------------------

	stop
99909	write(nat,*)' (09) error: no unique definition of nodes'
	goto 99
99915	write(nat,*)' (15) error: no unique definition of elements'
	goto 99
99918	write(nat,*)' (18) error: not used nodes present'
	goto 99
99920	write(nat,*)' (20) error: extern/intern node numbers'
	goto 99
99930	write(nat,*)' (30) error: elements with clockwise nodes'
	goto 99
99931	write(nat,*)' (31) error: zero area elements'
	goto 99
99932	write(nat,*)' (32) error: processing depth'
	goto 99
99934	write(nat,*)' (34) error: ngr1 is wrong (internal error'
	goto 99
99935	write(nat,*)' (35) error: irregular connections'
	goto 99
99940	write(nat,*)' (40) error: generic error'
	goto 99
99999	write(nat,*)' (99) error: no grid in grd file'
	goto 99

   99	continue

	write(6,*) 'no bas file written... aborting'
	stop ' error stop shypre'
 	end

!*****************************************************************

        subroutine uniqn(n,iv,index,bstop)

! controlls uniqueness of numbers in iv

        implicit none

        integer n
        logical bstop
        integer iv(n)
        integer index(n)

        integer k1,k2,i

        k1=iv(index(1))
        do i=2,n
          k2=iv(index(i))
          if(k1.eq.k2) then
            write(6,*)' warning : ',k1,' not unique'
            bstop=.true.
          end if
          k1=k2
	end do

        return
        end

!*****************************************************************

        subroutine uniqe(nel,nen3v,ipev,bstop)

! controlls uniqueness of elements (shell for equale)

	use mod_shypre
	use mod_sort

        implicit none

        integer nel
        logical bstop
        integer nen3v(3,nel)
        integer ipev(nel)

        integer ie,ii,ka,km,k1,k,ne,i
        integer ie1,ie2
        integer isum,iold,inew,itest
        integer iaux(nel),index(nel),elist(2,nel)

	logical bwrite

	bwrite = .not. bquiet

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
	    do i=1,ne
	      !write(6,*) '  elements: ',ipev(elist(1,i)),ipev(elist(2,i))
	    end do
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

        subroutine chexin(nkn,nel,nen3v,ipv,index,bstop)

! changing extern with intern node numbers in element index

	use mod_sort

        implicit none

        integer nkn,nel
        logical bstop
        integer nen3v(3,nel)
        integer ipv(nkn)
        integer index(nkn)

        integer ie,ii,i,kn
        !integer locate

        do ie=1,nel
          do ii=1,3
            kn=nen3v(ii,ie)
            i=locate(nkn,ipv,index,kn)
            if(i.le.0) then
              write(6,*)' error: node',kn,' not found'
              bstop=.true.
            else
              nen3v(ii,ie)=i
            end if
          end do
        end do

        return
        end

!*****************************************************************

        subroutine areaz(nkn,nel,nen3v,ipev,xgv,ygv,bstop)

! test for zero area elements

	use mod_shypre
	use mod_sort

        implicit none

        integer nkn,nel
        logical bstop
        integer nen3v(3,nel)
        integer ipev(nel)
        real xgv(nkn),ygv(nkn)

	logical bwarn
        integer ie,ii,ii2,ii3,k1,k2,k3,ke,i
	integer nsmall
	integer, parameter :: nsmax = 20
        real a,x1,x2,x3,y1,y2,y3,p
        real area(nel),amin,amax
	real eps

	real get_area, area_elem

	amin = 1.e+30
	amax = -1.e+30
	eps = eps_area
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
          !a = get_area(x1,y1,x2,y2,x3,y3)
	  a = area_elem(ie)			!honors spherical
	  a = abs(a)
	  amin = min(amin,a)
	  amax = max(amax,a)
	  area(ie) = a
	end do

	!call sort(nel,area)
	!do i=1,10
	!  write(6,*) i,area(i),area(i)/amax
	!end do
	!write(6,*)
	!do i=nel-10,nel
	!  write(6,*) i,area(i),area(i)/amax
	!end do
	!call histo_area(nel,area)
	!stop

	nsmall = 0

	do ie=1,nel
	  a = area(ie)
	  p = a / amax
	  ke = ipev(ie)
          if(a.eq.0.) then
            write(6,*)' error: area in element ',ipev(ie) &
     &                     ,' is zero'
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

	if( .not. bquiet .and. bwarn ) then
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

        subroutine check_sidei_0(nkn,nel,nen3v,ipv,ng,iknot,ngr1,bstop)

! set up side index and find grade

        implicit none

        integer nkn,nel,ngr
        integer ngr1
        integer nen3v(3,nel)
        integer ipv(nkn)
        integer ng(nkn)
        integer iknot(ngr1,nkn)
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

	subroutine knscr(nkn,ngr,ngr1,iknot)

! scrambles side index for cmv and rosen

! -> please eliminate need for this

	implicit none

	integer nkn,ngr,ngr1
        integer iknot(ngr1*nkn)

	integer k,j

        do k=1,nkn
          do j=1,ngr
            iknot((k-1)*ngr+j)=iknot((k-1)*ngr1+j)
          end do
        end do

	return
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

        subroutine bandop(nkn,ngr1,ipv,iphv,kphv,ng,iknot,kvert &
     &				,bopti,bauto,bwrite)

! optimize band width

        implicit none

        integer nkn,ngr1
        logical bopti,bauto
        integer ipv(nkn)
        integer iphv(nkn),kphv(nkn)
        integer ng(nkn)
        integer iknot(ngr1,nkn)
        integer kvert(2,nkn)
	logical bwrite,bloop

	integer iantw

	call ininum(nkn,iphv,kphv)

	if( .not. bopti ) return

	if( bauto ) then
	  call ininum(nkn,iphv,kphv)
	  call optest('before optimization: ',nkn,ipv,iphv,kphv)
	  call cmgrade(nkn,ngr1,ipv,iphv,kphv,ng,iknot,1,4,bwrite)
	  call optest('after Cuthill McKee: ',nkn,ipv,iphv,kphv)
	  call revnum(nkn,iphv,kphv,bwrite)
	  call optest('after reversing nodes: ',nkn,ipv,iphv,kphv)
          call rosen(nkn,ngr1,iphv,kphv,ng,iknot,kvert,bwrite)
	  call optest('after Rosen: ',nkn,ipv,iphv,kphv)
	  return
	end if

	bloop = .true.

        do while( bloop )
	  call ininum(nkn,iphv,kphv)

!	  call anneal(nkn,ngr1,kphv,ng,iknot,iphv,kvert)

          if( iantw(' Cuthill McKee algorithm ?') .gt. 0 ) then
            call cmv(nkn,ngr1,ipv,iphv,kphv,ng,iknot,bwrite)
	  end if

          if( iantw(' Reverse numbering of nodes ?') .gt. 0 ) then
	    call revnum(nkn,iphv,kphv,bwrite)
	  end if

          if( iantw(' Rosen algorithm ?') .gt. 0 ) then
            call rosen(nkn,ngr1,iphv,kphv,ng,iknot,kvert,bwrite)
	  end if

	  bloop = iantw(' Repeat optimization of bandwidth ?') .gt. 0
        end do

	bopti = .true.

        end

!**********************************************************

	subroutine revnum(nkn,iphv,kphv,bwrite)

! reverses numbering of nodes

	implicit none

	integer nkn
	integer iphv(nkn), kphv(nkn)
	logical bwrite

	integer i

	if( bwrite ) write(6,*) 'Applying reverse algorithm...'

	do i=1,nkn
          kphv(i) = nkn+1-kphv(i)
          iphv(kphv(i)) = i
	end do

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

	subroutine zernum(nkn,iphv,kphv)

! zeros numbering of nodes

	implicit none

	integer nkn
	integer iphv(nkn), kphv(nkn)

	integer i

	do i=1,nkn
	  iphv(i) = 0
	  kphv(i) = 0
	end do

	end

!**********************************************************

        subroutine bandex(nkn,nel,nen3v,kphv &
     &                      ,ipv,iarnv,xgv,ygv,hkv)

! exchange nodes after optimization

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        integer kphv(nkn)
        integer ipv(nkn)
        integer iarnv(nkn)
        real xgv(nkn),ygv(nkn)
        real hkv(nkn)

        integer ie,ii
	integer neaux(3,nel)

        do ie=1,nel
          do ii=1,3
            neaux(ii,ie)=nen3v(ii,ie)
          end do
        end do

        do ie=1,nel
          do ii=1,3
            nen3v(ii,ie)=kphv(neaux(ii,ie))
          end do
        end do

!       copy arrays with kphv as rank table

        call icopy(nkn,ipv,kphv)
        call icopy(nkn,iarnv,kphv)
        call rcopy(nkn,xgv,kphv)
        call rcopy(nkn,ygv,kphv)
        call rcopy(nkn,hkv,kphv)

        return
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

        subroutine renel(nel,nen3v,ipev,iarv,hev,ierank)

! renumbering of elements

! we construct iedex newly and use iaux,iedex as aux arrays
! iedex is also used as a real aux array for rcopy	-> changed
! neaux is probably he3v (real) used as an aux array	-> changed

	use mod_sort

        implicit none

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

	subroutine init_hm3v(nel,hm3v,hinit)

	implicit none

	integer nel
	real hm3v(3,nel)
	real hinit

	integer ie,ii

	do ie=1,nel
	  do ii=1,3
	    hm3v(ii,ie) = hinit
	  end do
	end do

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

        subroutine helem(nel,nhd,iphev,iaux,iedex &
     &                    ,ipev,hm3v,hev,bstop)

! depth by elements

	use mod_sort

        implicit none

        integer nel,nhd
        logical bstop
        integer iaux(nel)
        integer iedex(nel)
        integer ipev(nel),iphev(nel)
        real hm3v(3,nel),hev(nel)

        integer ie,i,iel,ii
        !integer locate

        do ie=1,nel
          iaux(ie)=0
        end do

	!write(6,*) 'helem: ',nel,nhd

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

        subroutine hnode(nkn,nel,nhd,nen3v,iphv,iaux,ipdex,ipv &
     &                      ,hm3v,hkv,bstop)

! depth by nodes

	use mod_sort

        implicit none

        integer nkn,nel,nhd
        logical bstop
        integer iaux(nkn)
        integer ipdex(nkn)
        integer ipv(nkn),iphv(nkn)
        integer nen3v(3,nel)
        real hm3v(3,nel),hkv(nkn)

        integer ie,ii,i,k,kn
        !integer locate
	real h,hflag

	hflag = -999.

        do k=1,nkn
          iaux(k)=0
        end do

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

	subroutine gtest_0(text,nelddi,nkn,nel,nen3v)

	implicit none

	character*(*) text
	integer nelddi,nkn,nel
	integer nen3v(3,nelddi)

	integer ie,ii,k
	integer iii
	logical berror
	
	berror = .false.

	write(6,*) 'testing ... ',text
	write(6,*) nelddi,nkn,nel

	do ie=1,nel
	  do ii=1,3
	    iii = mod(ii,3) + 1
	    if( nen3v(ii,ie) .eq. nen3v(iii,ie) ) then
	      if( .not. berror ) then
		write(6,*) ie,(nen3v(iii,ie),iii=1,3)
		berror = .true.
	      end if
	    end if
	  end do
	end do

	if( berror ) then
	  write(6,*) 'testing has found errors...'
	end if
	end

!**********************************************************

	subroutine optest(text,nkn,ipv,iphv,kphv)

	implicit none

	character*(*) text
	integer nkn
	integer ipv(nkn)
	integer iphv(nkn)
	integer kphv(nkn)

	integer i

	do i=1,nkn
	  if( iphv(i) .gt. nkn ) goto 99
	  if( kphv(i) .gt. nkn ) goto 99
	  if( iphv(kphv(i)) .ne. i ) goto 99
	  if( kphv(iphv(i)) .ne. i ) goto 99
	end do

	return
   99	continue
	write(6,*) '*** Error in optest: '
	write(6,*) text
	write(6,*) 'error in pointers...'
	write(6,*) 'problem is close to following node...'
	write(6,*) 'node (intern/extern): ',i,ipv(i)
	write(6,*) i,nkn
	write(6,*) iphv(i),kphv(i)
	write(6,*) iphv(kphv(i)),kphv(iphv(i))
	write(6,*) 'maybe the domain is not connected...'
	stop 'error stop optest: pointers'
	end

!**********************************************************

	subroutine ketest(nel,nen3v)

! checks uniquness of nodes in elements

	implicit none

	integer nel
	integer nen3v(3,nel)

	integer ie,ii
	integer kn1,kn2,kn3

	do ie=1,nel
	  kn1 = nen3v(1,ie)
	  kn2 = nen3v(2,ie)
	  kn3 = nen3v(3,ie)
	  if( kn1 .eq. kn2 .or. kn1 .eq. kn3 .or. kn2 .eq. kn3 ) then
	    write(6,*) ie,kn1,kn2,kn3
	    stop 'error stop ketest: not unique nodes in element'
	  end if
	end do

	!write(77,*)
	!do ie=1,10
	!  write(77,*) ie,(nen3v(ii,ie),ii=1,3)
	!end do

	return
	end

!**********************************************************

        function get_area(x1,y1,x2,y2,x3,y3)

! computes area of triangle

        implicit none

        real get_area
        real x1,y1,x2,y2,x3,y3

        get_area = 0.5 * ( (x2-x1) * (y3-y1) - (x3-x1) * (y2-y1) )

        end

!**********************************************************
!**********************************************************
!**********************************************************

	subroutine handle_partition(nn,ne,knrank,ierank)

	use clo
	use grd
	use basin

	implicit none

	integer nn,ne
	integer knrank(nn)
	integer ierank(ne)

	logical bnepart,bgrd
	integer nnpart,nepart
	integer, allocatable :: area_node(:)
	integer, allocatable :: area_elem(:)

	integer i
	character*80 grdfile

        call clo_get_option('nepart',bnepart)
        call clo_get_option('partition',grdfile)
	bgrd = ( grdfile /= ' ' )

	if( .not. bnepart .and. .not. bgrd ) return

	if( bnepart .and. bgrd ) then
	  write(6,*) 'only one of -partition and -nepart can be given'
	  stop 'error stop handle_partition: options'
	end if

	allocate(area_node(nn))
	allocate(area_elem(ne))

	if( bgrd ) then
	  write(6,*) 'reading partitioning file ',trim(grdfile)
	  call grd_read(grdfile)
	  if( nk_grd /= nn ) goto 99
	  if( ne_grd /= ne ) goto 99
	  area_node = ianv
	  area_elem = iaev
	else
	  area_node = iarnv
	  area_elem = iarv
	end if

	call renumber_partition(nn,area_node,nnpart)
	call renumber_partition(ne,area_elem,nepart)

	if( bgrd ) then
          call icopy(nn,area_node,knrank)
          call icopy(ne,area_elem,ierank)
	end if

	call basin_set_partition(nn,ne,nnpart,nepart,area_node,area_elem)

	return
   99	continue
	write(6,*) nk_grd,nn
	write(6,*) ne_grd,ne
	stop 'error stop handle_partition: incompatibility'
	end

!**********************************************************

	subroutine renumber_partition(n,area,npart)

	implicit none

	integer n
	integer area(n)
	integer npart

	integer i,ia,imax
	integer nmin,nmax
	integer, allocatable :: table_in(:),table_out(:)

	nmin = minval(area)
	nmax = maxval(area)

	if( nmin == nmax ) then		!no partition
	  npart = 0
	  area = 0
	  return
	end if

	if( nmin < 0 .or. nmax < 0 ) then
	  write(6,*) nmin,nmax
	  stop 'error stop renumber_partition: nmin,nmax'
	end if

	allocate(table_in(0:nmax),table_out(0:nmax))

	table_in = 0
	table_out = 0

	do i=1,n
	  ia = area(i)
	  table_in(ia) = table_in(ia) + 1
	end do

	imax = -1
	do ia=0,nmax
	  if( table_in(ia) > 0 ) then
	    imax = imax + 1
	    table_out(ia) = imax
	  end if
	end do
	npart = imax

	do i=1,n
	  ia = area(i)
	  area(i) = table_out(ia)
	end do

	end

!**********************************************************
!**********************************************************
!**********************************************************

	subroutine shypre_init(grdfile)

	use clo
	use mod_shypre

	implicit none

	character*(*) grdfile

	call clo_init('shypre','grd-file','3.0')

        call clo_add_info('pre-processes grd file and create bas file')

	call clo_add_sep('general options')
	call clo_add_option('info',.false.,'only give info on grd file')
	call clo_add_option('quiet',.false.,'be quiet in execution')
	call clo_add_option('silent',.false.,'do not write anything')

	call clo_add_sep('optimization options')
	!call clo_add_option('noopti',.false.,'do not optimize bandwidth')
	call clo_add_option('opti',.false.,'do not optimize bandwidth')
        call clo_add_option('manual',.false.,'manual optimization')

	call clo_add_sep('options for partitioning')
	call clo_add_option('partition file',' ' &
     &		,'use file containing partitioning')
	call clo_add_option('nepart',.false. &
     &		,'use node and elem type in file for partitioning')

	call clo_parse_options

        call clo_get_option('info',binfo)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)

        call clo_get_option('opti',bopti)
        !call clo_get_option('noopti',bnoopti)
        call clo_get_option('manual',bmanual)

	!bopti = .not. bnoopti
	bauto = .not. bmanual
	if( bsilent ) bquiet = .true.

	call shyfem_set_short_copyright(bquiet)
	if( .not. bsilent ) then
          call shyfem_copyright('shypre - pre-processing of GRD grid')
	end if

	call clo_check_files(1)
	call clo_get_file(1,grdfile)

	end

!**********************************************************

	subroutine histo_area(n,array)

	implicit none

	integer n
	real array(n)

	integer, parameter :: nbuck = 20
	integer ib,i,ntot
	integer ibuck(0:nbuck)
	real atot,adiff,a1,a2
	real abuck(0:nbuck)
	real x,y,a,y0,z,s,b
	real amin,amax

	amin = minval(array)
	amax = maxval(array)

	b = 0.1
	b = 5.
	b = 10.

	atot = amax - amin
	adiff = atot / nbuck
	a1 = amin + adiff/2.
	a2 = amax - adiff/2.
	y0 = exp(0.)/exp(b*1.0)
	abuck(0) = amin
	do ib=0,nbuck
	  x = ib/float(nbuck)
	  y = exp(b*x)/exp(b*1.0)
	  z = (y-y0)/(1.-y0)		![0-1]
	  s = amin + z*(amax-amin)
	  abuck(ib) = s
	  write(6,*) ib,x,y,z,s
	end do
	ibuck = 0
	
	do i=1,n
	  a = array(i)
	  do ib=1,nbuck
	    if( a > abuck(ib-1) .and. a <= abuck(ib) ) exit
	  end do
	  ib = ib - 1
	  if( ib < 0 ) write(6,*) '*** hist: ',ib,a,a1,a2
	  if( ib > nbuck ) write(6,*) '*** hist: ',ib,a,a1,a2
	  ibuck(ib) = ibuck(ib) + 1
	  !write(6,*) ib,a,abuck(ib)
	end do

	write(6,*) 'histogram areas: ',nbuck
	write(6,*) 'abuck0 = ',abuck(0)
	ntot = 0
	do ib=0,nbuck
	  write(6,*) ib,abuck(ib),ibuck(ib)
	  ntot = ntot + ibuck(ib)
	end do

	write(6,*) ntot,n
	if( n /= ntot ) stop 'error stop histo: ntot/=n'

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

