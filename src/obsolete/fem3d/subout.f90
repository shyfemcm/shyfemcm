
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998,2010,2019  Georg Umgiesser
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

! utility routines to read/write OUT file - file type 81
!
! reads files up to version 6
!
! contents :
!
! function rfout(iunit,nvers,itanf,itend,idt,idtout,href,hzoff,descrp)
! 			reads first record of file OUT
! function wfout(iunit,nvers,itanf,itend,idt,idtout,href,hzoff,descrp)
! 			writes first record of file OUT
! function rdout(iunit,nvers,it,nkn,xv)
!			reads data record of file OUT
! function wrout(iunit,nvers,it,nkn,xv)
!			writes data record of file OUT
!
! revision log :
!
! 08.05.1998	ggu	$$id81 - constant 81 replaced with variable id=81 (bug)
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2019	ggu	changed VERS_7_5_60
!
!************************************************************
!
	function rfout(iunit,nvers,itanf,itend,idt,idtout &
     &				,href,hzoff,descrp)
!
! reads first record of out file
!
! versions (first record) :
!	1	itanf,itend,idt,idtout
!	2-5	nvers
!	6-...	ntype,nvers
!
	implicit none
!
! arguments
	integer rfout
	integer iunit,nvers,itanf,itend,idt,idtout
	real href,hzoff
	character*80 descrp
! local
	integer ntype,ios
!
	rewind(iunit,iostat=ios)
!
	if(ios.ne.0) then
		write(6,*) 'Cannot rewind file for unit :',iunit
		rfout=71
		return
	end if
!
! first record - find out what version
!
	read(iunit,iostat=ios) itanf,itend,idt,idtout
!
	if(ios.eq.0) then	! no read error ==> first version
		nvers=1
	else if(ios.lt.0) then	!eof
		nvers=itanf
        else if(itanf.eq.81) then !new version
          nvers=itend
        else        !no type available --> up to version 5
          nvers=itanf
	end if

!        write(6,*) nvers,itanf,itend,ios
!
! now rewind and read first record
!
	rewind(iunit,iostat=ios)
!
	if(nvers.eq.1) then
		ntype=81
		read(iunit,iostat=ios) itanf,itend,idt,idtout
	else if(nvers.ge.2.and.nvers.le.5) then
		ntype=81
		read(iunit,iostat=ios) nvers
	else if(nvers.ge.6) then
		read(iunit,iostat=ios) ntype,nvers
	end if

!         write(6,*) ntype,nvers,ios
!
	if(ios.gt.0) then
		write(6,*) 'Error encountered while reading'
		write(6,*) 'first record of out file header'
		write(6,*) 'nvers =',nvers
		rfout=22
		return
	else if(ios.lt.0) then
		write(6,*) 'EOF encountered while reading'
		write(6,*) 'first record of out file header'
		write(6,*) 'nvers =',nvers
		rfout=21
		return
	end if
!
! control version number and type of file
!
	if(nvers.le.0.or.nvers.gt.6) then
		write(6,*) 'version not recognized : ',nvers
		rfout=11
		return
	end if
!
	if(ntype.ne.81) then
		write(6,*) 'rfout : Wrong type of file : ',ntype
		write(6,*) 'Expected : 81'
		rfout=15
		return
	end if
!
! second record
!
	if(nvers.eq.1) then
		hzoff=0.05
		href=0.
		descrp=' '
	else if(nvers.ge.2.and.nvers.le.6) then
		read(iunit,iostat=ios)	 itanf,itend,idt,idtout &
     &					,href &
     &					,hzoff &
     &					,descrp
	end if
!
	if(ios.gt.0) then	!error
		write(6,*) 'Error encountered while reading'
		write(6,*) 'second record of out file header'
		write(6,*) 'nvers =',nvers
		rfout=35
		return
	else if(ios.lt.0) then	!eof
		write(6,*) 'EOF encountered while reading'
		write(6,*) 'second record of out file header'
		write(6,*) 'nvers =',nvers
		rfout=36
		return
	end if
!
	rfout=0
!
	return
	end
!
!********************************************************************
!
	function wfout(iunit,nvers,itanf,itend,idt,idtout &
     &				,href,hzoff,descrp)
!
! writes first record of out file
!
! versions (first record) :
!	1	itanf,itend,idt,idtout
!	2-5	nvers
!	6-...	ntype,nvers
!
	implicit none
!
! arguments
	integer wfout
	integer iunit,nvers,itanf,itend,idt,idtout
	integer id
	real href,hzoff
	character*80 descrp
!
	rewind(iunit)
!
	id = 81
	if(nvers.eq.0) nvers=6
!
	if(nvers.eq.1) then
		write(iunit)		 itanf,itend,idt,idtout
	else if(nvers.ge.2.and.nvers.le.5) then
		write(iunit)		 nvers
		write(iunit)		 itanf,itend,idt,idtout &
     &					,href &
     &					,hzoff &
     &					,descrp
	else if(nvers.eq.6) then
		write(iunit)		 id,nvers		!$$id81
		write(iunit)		 itanf,itend,idt,idtout &
     &					,href &
     &					,hzoff &
     &					,descrp
	else
		write(6,*) 'version not recognized : ',nvers
		wfout=11
		return
	end if
!
	wfout=0
!
	return
	end
!
!************************************************************
!
	function rdout(iunit,nvers,it,nkn,nel,xv,zenv,u,v)
!
! reads data record of out file
!
	implicit none
!
! arguments
	integer rdout
	integer iunit,nvers,it,nkn,nel
	real xv(1),zenv(1),u(1),v(1)
! local
	integer ios,nk,ne,i
!
! control version number
!
	if(nvers.le.0.or.nvers.gt.6) then
		write(6,*) 'version not recognized : ',nvers
		rdout=11
		return
	end if
!
! time record
!
	if(nvers.ge.1.and.nvers.le.5) then
		nk=nkn
		ne=nel
		read(iunit,iostat=ios) it
	else if(nvers.eq.6) then
		read(iunit,iostat=ios) it,nk,ne
	end if
!
	if(ios.gt.0) then	!error
		write(6,*) 'Error while reading'
		write(6,*) 'time record of out file'
		rdout=25
		return
	else if(ios.lt.0) then	!eof
		rdout=-1
		return
	end if
!
! check dimensions
!
	if(nkn.eq.0) then
!		ok -> skip
	else if(nk.gt.nkn) then	!array too small
		write(6,*) 'Too much data in node record :',nk
		write(6,*) 'Can read only',nkn,' data'
	else
		nkn=nk
	end if
!
	if(nel.eq.0) then
!		ok -> skip
	else if(ne.gt.nel) then	!array too small
		write(6,*) 'Too much data in element record :',ne
		write(6,*) 'Can read only',nel,' data'
	else
		nel=ne
	end if
!
! data record
!
	if(nvers.ge.1.and.nvers.le.5) then
		read(iunit,iostat=ios)	 (xv(i),i=1,3*nkn)
	else if(nvers.eq.6) then
		read(iunit,iostat=ios)	 (xv(i),i=1,3*nkn) &
     &					,(zenv(i),i=1,3*nel) &
     &					,(u(i),i=1,nel) &
     &					,(v(i),i=1,nel)
	end if
!
	if(ios.gt.0) then	!error
		write(6,*) 'Error while reading'
		write(6,*) 'data record of out file'
		write(6,*) 'it      : ',it
		write(6,*) 'nkn,nel : ',nkn,nel
		rdout=35
		return
	else if(ios.lt.0) then	!eof
		write(6,*) 'EOF encountered while reading'
		write(6,*) 'data record of out file'
		write(6,*) 'it =',it
		write(6,*) 'nkn,nel : ',nkn,nel
		rdout=31
		return
	end if
!
	rdout=0
!
	return
	end
!
!************************************************************
!
	function wrout(iunit,nvers,it,nkn,nel,xv,zenv,u,v)
!
! writes data record of out file
!
	implicit none
!
! arguments
	integer wrout
	integer iunit,nvers,it,nkn,nel
	real xv(1),zenv(1),u(1),v(1)
! local
        integer i
!
! control version number
!
	if(nvers.le.0.or.nvers.gt.6) then
		write(6,*) 'version not recognized : ',nvers
		wrout=11
		return
	end if
!
	if(nvers.ge.1.and.nvers.le.5) then
		write(iunit)	 it
		write(iunit)	 (xv(i),i=1,3*nkn)
	else if(nvers.eq.6) then
		write(iunit)	 it,nkn,nel
		write(iunit)	 (xv(i),i=1,3*nkn) &
     &				,(zenv(i),i=1,3*nel) &
     &				,(u(i),i=1,nel) &
     &				,(v(i),i=1,nel)
	end if
!
	wrout=0
!
	return
	end
!
