
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003,2006,2009-2011,2014-2019  Georg Umgiesser
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

! topological set up routines
!
! contents :
!
! subroutine nknel(nnkn,nnel,nnbn,nnli,nnis,nnod)
!			 computes parameters about total numbers
! subroutine newlnk
!			 administrates new topological routines
! function winkk(k)
!			 total angle at node k
! subroutine setnar(nar)
!			 sets number of areas
! subroutine setwnk(wink)
!			 sets total angle at boundary nodes
! subroutine arper
!			 computes statistics about area of submerged zones
!
! revision log :
!
! 01.08.2003	ggu	created from sublnk.f
! 13.08.2003	ggu	new name set_link_info for newlnk
! 17.02.2006	ggu	do not all anymore zvbnds()
! 28.04.2009	ggu	links re-structured
! 23.03.2010	ggu	changed v6.1.1
! 31.05.2011	ggu	changed VERS_6_1_23
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_2
! 19.01.2015	ggu	changed VERS_7_1_3
! 05.05.2015	ggu	changed VERS_7_1_10
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 16.12.2015	ggu	changed VERS_7_3_16
! 18.12.2015	ggu	changed VERS_7_3_17
! 28.04.2016	ggu	changed VERS_7_5_9
! 27.06.2016	ggu	changed VERS_7_5_16
! 05.12.2017	ggu	changed VERS_7_5_39
! 03.04.2018	ggu	changed VERS_7_5_43
! 11.05.2018	ggu	changed VERS_7_5_47
! 16.02.2019	ggu	changed VERS_7_5_60
! 06.11.2019	ggu	eliminated femtime
! 06.06.2023	ggu	do not write newlnk - not working in MPI
!
!*****************************************************************

        subroutine nknel(nnkn,nnel,nnbn,nnli,nnis,nnod)

! computes parameters about total numbers
!
! all numbers are refered to actual elements and nodes inside system
!
! nnkn   total number of nodes
! nnel   total number of elements
! nnbn   total number of boundary nodes
! nnli   total number of links
! nnis   total number of islands
! nnod   total number of not allowed junctions (should be 0)
!
! formulas :
!
! nel = 2*nkn - nbn - 2 + 2*nis - nod = nbn - 2 + 2*nin + 2*nis - nod
! nin = nkn - nbn
! nli = ( 3*nel + nbn + nod ) / 2
!
! the following is always true:
!
! nel < 2*nkn
! nli = nel + (1/2)*nel + nbn < nel + nkn + nbn <= nel + 2*nkn   so
! nli < 2*nkn + nel

	use mod_geom
	use mod_geom_dynamic
	use basin, only : nkn,nel,ngr,mbw
	use shympi

        implicit none

! arguments
        integer nnkn,nnel,nnbn,nnli,nnis,nnod
! local
        integer k,ie,n,i,ne,ntot
	integer elems(maxlnk)
        logical bin
! statement functions
        logical iskbnd,iskins,iseins
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
        iskins(k) = inodv(k).ne.-2
        iseins(ie) = ie.gt.0.and.iwegv(ie).eq.0

        nnod=0
	ntot = nkn_unique
        do k=1,ntot
	  call get_elems_around(k,maxlnk,ne,elems)

          n=0
	  ie = elems(ne)
          bin=iseins(ie)
          do i=1,ne
	    ie = elems(i)
            if( iseins(ie) ) then
              if( .not. bin ) n=n+1
            end if
            bin=iseins(ie)
          end do

          if(n.gt.1) nnod=nnod+n-1
        end do

        nnbn=0
        nnkn=0
	ntot = nkn_unique
        do k=1,ntot
          if( iskbnd(k) ) nnbn=nnbn+1
          if( iskins(k) ) nnkn=nnkn+1
        end do

        nnel=0
	ntot = nel_unique
        do ie=1,ntot
          if( iseins(ie) ) nnel=nnel+1
        end do

        nnel = shympi_sum(nnel)
        nnbn = shympi_sum(nnbn)
        nnod = shympi_sum(nnod)
        nnkn = shympi_sum(nnkn)

        nnis=(nnel+nnbn-2*nnkn+2+nnod)/2
        nnli=(3*nnel+nnbn+nnod)/2

        end

!****************************************************************

        subroutine set_link_info

! administrates new topological routines
!
! ... iwegv has already been set

	use mod_geom
	use basin, only : nkn,nel,ngr,mbw
	use shympi

        implicit none

! local
        character*80 nam,dir,file
        real wink
        integer nnkn,nnel,nnbn,nnli,nnis,nnod,nnar,nnnis,nnnel
! function
        integer iround,ifileo
! save
        integer, save :: iuinfo = 0
	logical, save :: binfo = .false.

	if( .not. binfo ) return			!not working with MPI

        if(iuinfo.eq.0) call getinfo(iuinfo)		!get unit number

        call setnod

        call nknel(nnkn,nnel,nnbn,nnli,nnis,nnod)

! shympi - FIXME - next two calls might not work properly

        call setnar(nnar)	!number of areas
        call setwnk(wink)	!total angle at boundary nodes

        nnnis = iround( (wink-nnbn*180.)/360. ) + nnar
        nnnel=2*nnkn-nnbn+2*(nnnis-nnar)-nnod

        if(shympi_is_master()) then
          write(iuinfo,'(a,10i6)') ' newlnk: ' &
     &          ,nnkn,nnel,nnbn,nnli,nnis,nnod,nnar &
     &          ,nnnis,nnnel-nnel,nel_global-nnel
	end if

        end

!****************************************************************

        function winkk(k)

! total angle at node k
!
! k       node to links
!
	use mod_geom
	use mod_geom_dynamic
	use evgeom
	use basin, only : nkn,nel,ngr,mbw

        implicit none

! arguments
	real winkk
        integer k
! local
        integer i,n,ie
	integer elems(maxlnk)
        double precision w
! functions
        integer ithis

	call get_elems_around(k,maxlnk,n,elems)

        w=0.
        do i=1,n
          ie=elems(i)
          if( iwegv(ie).eq.0 ) then
            w = w + ev(10+ithis(k,ie),ie)
          end if
        end do

        winkk=w

        end

!****************************************************************
!
        subroutine setnar(nar)
!
! sets number of areas
!
! nar   number of areas (return)
!
	use mod_geom
	use mod_geom_dynamic
	use basin
	use shympi

        implicit none
!
! arguments
        integer nar
! local
        integer i,ie,ieo,ien,n1,n2
        logical btest
! functions
        integer inext
! data
        data btest /.false./
!
        do ie=1,nel
          if(iwegv(ie).gt.0) iwegv(ie)=iwegv(ie)+3
        end do
!
        n1=0
        n2=0
        nar=0
        do ie=1,nel
          if(iwegv(ie).eq.0) then
            i=1
            ieo=ie
            do while(iwegv(ieo).lt.3)
              iwegv(ieo)=iwegv(ieo)+1
              i=mod(i,3)+1
              ien=ieltv(i,ieo)
              n1=n1+1
!              write(6,*) ieo,iwegv(ieo),i,ien,iwegv(ien)
!              read(5,'(i10)') n
              if(ien.gt.0) then !possible to enter
               if(iwegv(ien).lt.3) then !possible to enter
                if(iwegv(ien).eq.0.or.iwegv(ieo).eq.3) then !may enter
                  i=inext( nen3v(mod(i,3)+1,ieo) , ien )
                  ieo=ien
                  n2=n2+1
                end if
               end if
              end if
            end do
            nar=nar+1
          end if
        end do
!
        if(btest) then
          write(6,*) nar,n1,n2
          write(6,*)
        end if

        !call shympi_comment('shympi_min(nar)')
        nar = shympi_min(nar)

        do ie=1,nel
          if(iwegv(ie).lt.3) then
            write(6,*) 'we forgot something : ',ie,iwegv(ie)
            stop 'error stop setnar'
          end if
          iwegv(ie)=iwegv(ie)-3
        end do
!
        return
        end
!
!****************************************************************
!
        subroutine setwnk(wink)
!
! sets total angle at boundary nodes
!
! wink  total angle (return)
!
	use mod_geom_dynamic
	use evgeom
	use basin
	use shympi

        implicit none
!
! arguments
        real wink
! local
        integer ie,ii,k
        double precision w
! functions
        logical iskbnd
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
!
! sum angles
!
        w=0.
        do ie=1,nel
          if(iwegv(ie).eq.0) then !element is in system
            do ii=1,3
              k=nen3v(ii,ie)
              if( iskbnd(k) ) w=w+ev(10+ii,ie)
            end do
          end if
	end do

        wink=w

        !call shympi_comment('shympi_sum(wink)')
        wink = shympi_sum(wink)

        return
        end

!****************************************************************
!
        subroutine arper
!
! computes statistics about area of submerged zones
!
! ... iwegv has already been set
!
	use mod_geom_dynamic
	use evgeom
	use basin, only : nkn,nel,ngr,mbw

        implicit none
!
! local
        real arin,arout,artot,area
        integer ie
! functions
        logical isein
        isein(ie) = iwegv(ie).eq.0

        arin=0.
        arout=0.
        artot=0.

        do ie=1,nel
          area=ev(10,ie)
          if(isein(ie)) then
            arin=arin+area
          else
            arout=arout+area
          end if
          artot=artot+area
        end do

        arin=100.*arin/artot
        arout=100.*arout/artot
        artot=12.*artot

!        write(88,'(i8,f12.3,2f12.2,e12.4)') dtime,arin,arout,artot

        end

!****************************************************************

