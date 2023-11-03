
!--------------------------------------------------------------------------
!
!    Copyright (C) 1995,1997-1999,2001,2003,2010-2017  Georg Umgiesser
!    Copyright (C) 2019  Georg Umgiesser
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

! interactive assignment of files with apn
!
! contents :
!
! subroutine assnam(mode)                       assigns run/bas interactivly
! function iapini(mode,nknddi,nelddi,matdim)    init routine for ap routines
! subroutine pardef                             reads nls and fnm
!
! revision log :
!
! 22.03.1995	ggu	$$EXINEL - call to exinel changed to ieint
! 21.05.1997	ggu	$$NOBAS - in iapini bug if basin is not opened
! 21.01.1998	ggu	$$APNPAR - only look for apn file in local directory
! 01.05.1998	ggu	$$IMEM - close file only for imem > 0
! 01.05.1998	ggu	interactive questions beautified
! 20.06.1998	ggu	general clean up
! 07.08.1998	ggu	handles negative mode -> do not change names
! 12.02.1999	ggu	useless data structures deleted in iapini
! 12.02.1999	ggu	in assnam do not ask for apn file
! 12.02.1999	ggu	honor auto mode
! 12.02.1999	ggu	read runnam from memfil only if not in apn file
! 05.12.2001	ggu	fixed compiler error with -Wall -pedantic
! 20.06.2003	ggu	in iapini statement shiftet for compiler error Origin
! 23.03.2010	ggu	changed v6.1.1
! 15.07.2011	ggu	adjusted, ideffi substituted
! 19.03.2012	ggu	if no basin is given return with "error"
! 27.02.2013	ggu	pass what parameter into nlsa
! 03.05.2013	ggu	changed VERS_6_1_63
! 05.09.2013	ggu	read_apn_file() needs integer now
! 12.09.2013	ggu	changed VERS_6_1_67
! 18.06.2014	ggu	changed VERS_6_1_77
! 12.12.2014	ggu	changed VERS_7_0_9
! 23.12.2014	ggu	changed VERS_7_0_11
! 14.01.2015	ggu	reorganized and cleaned
! 23.01.2015	ggu	changed VERS_7_1_4
! 10.02.2015	ggu	debugged, bcompat gives compatibility with old versions
! 26.02.2015	ggu	changed VERS_7_1_5
! 29.04.2015	ggu	generic changes - now works as anticipated
! 05.06.2015	ggu	changed VERS_7_1_12
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_53
! 20.07.2015	ggu	changed VERS_7_1_81
! 28.04.2016	ggu	changed VERS_7_5_9
! 11.07.2017	ggu	changed VERS_7_5_30
! 14.11.2017	ggu	changed VERS_7_5_36
! 16.02.2019	ggu	changed VERS_7_5_60
!
! notes :
!
!	call iapini
!		->   call pardef
!			->   call nlsina
!			->   call fnminh
!	  		->   call read_apn_file(-1)
!				->   call nlsa
!		->   call assnam
!		->   call read_bas_file(nknddi,nelddi)
!		->   call sp131k(matdim)
!
! .memory is only changed with bask == .true.
! .memory is read with bmem == .true., but not changed
!
!****************************************************************

	subroutine assnam(mode)

! assigns name for run and basin (also interactivly)
!
! mode          1:assign basin 2:assign run 4:assign parameter file
!
! mode < 0	do not ask for names
!
! combinations are possible, all together : 7

	implicit none

	integer mode

	logical bchange,bask
	character*80 runnam,basnam
	character*80 runmem,basmem

!---------------------------------------------------------------------
! read memory file
!---------------------------------------------------------------------

	bask = mode .gt. 0
	bchange = .false.

	call read_memory(runmem,basmem)

!---------------------------------------------------------------------
! if we already have the names through the title section -> use these
!	else use the names just read from memory file
!---------------------------------------------------------------------

	call getfnm('runnam',runnam)
	call getfnm('basnam',basnam)

	runnam = adjustl(runnam)
	basnam = adjustl(basnam)

	if( runnam .eq. ' ' ) runnam = runmem
	if( basnam .eq. ' ' ) basnam = basmem

!---------------------------------------------------------------------
! if necessary ask for new values
!---------------------------------------------------------------------

	if( bask ) then
	  stop 'error stop: bask==.true. not supported any more'
	  !call ask_memory(mode,runnam,basnam,bchange)
	end if

!---------------------------------------------------------------------
! see if we have what we need
!---------------------------------------------------------------------

	if( btest(abs(mode),0) ) then
	  if( basnam == ' ' ) then
	    write(6,*) 'no basin given... exiting'
	    stop
	  end if
	  write(6,*) 'Name of basin      : ',trim(basnam)
	end if

	if( btest(abs(mode),1) ) then
	  if( runnam == ' ' ) then
	    write(6,*) 'no simulation given... exiting'
	    stop
	  end if
	  write(6,*) 'Name of simulation : ',trim(runnam)
	end if

!---------------------------------------------------------------------
! writes new information and new values
!---------------------------------------------------------------------

	call putfnm('runnam',runnam)
	call putfnm('basnam',basnam)

	if(bchange) call write_memory(runnam,basnam)

!---------------------------------------------------------------------
! end of routine
!---------------------------------------------------------------------

	end

!*************************************************************

	subroutine ap_set_names(basin,simul)

! sets basin and simulation

	character*(*) basin,simul

	logical haspar

	if( .not. haspar('runnam') ) call addfnm('runnam',' ')
	if( .not. haspar('basnam') ) call addfnm('basnam',' ')

	if( simul .ne. ' ' ) call putfnm('runnam',simul)
	if( basin .ne. ' ' ) call putfnm('basnam',basin)

	end

!*************************************************************

	subroutine ap_get_names(basin,simul)

! sets basin and simulation

	character*(*) basin,simul

	logical haspar

	basin = ' '
	simul = ' '

	if( haspar('runnam') ) call getfnm('runnam',simul)
	if( haspar('basnam') ) call getfnm('basnam',basin)

	end

!*************************************************************
!*************************************************************
!*************************************************************

	subroutine ap_init_basin

! shell for just reading basin

	call ap_init(.false.,1,0,0)

	end

!*************************************************************

	subroutine ap_init(bask,mode,nknddi,nelddi)

! initializes post processing

	implicit none

	logical bask
	integer mode,nknddi,nelddi

	if( bask ) then
	  call iap_init(mode,nknddi,nelddi,0)
	else
	  call iap_init(-mode,nknddi,nelddi,0)
	end if

	end

!*************************************************************

	function iapini(mode,nknddi,nelddi,matdim)

! init routine for ap routines
!
! calls pardef
! calls assnam
! reads basin
!
! iapini        1:success 0:error
! mode          for assnam 1:assign basin 2:assign run 3:both
! nknddi        dimension for nkn for reading bas file
! nelddi        dimension for nel for reading bas file
! matdim        (probably useless...)
!
! mode negative: do not ask for new basin and simulation

	implicit none

	integer iapini
	integer mode,nknddi,nelddi,matdim

	logical bcompat
	integer nmode,iauto

	real getpar

	bcompat = .false.	!set to .true. for compatibility

	call pardef		!we need this for iauto

	nmode = mode
	iauto = nint(getpar('iauto'))
	if( iauto .ne. 0 ) nmode = -abs(mode)

	if( .not. bcompat ) nmode = -abs(mode)	!never ask

	call iap_init(nmode,nknddi,nelddi,matdim)

	iapini = 1

	end

!*************************************************************

	subroutine iap_init(mode,nknddi,nelddi,matdim)

! init routine for ap routines
!
! calls pardef
! calls assnam
! reads basin
!
! iapini        1:success 0:error
! mode          for assnam 1:assign basin 2:assign run 3:both
! nknddi        dimension for nkn for reading bas file
! nelddi        dimension for nel for reading bas file
! matdim        (probably useless...)
!
! mode negative: do not ask for new basin and simulation

	implicit none

	integer mode,nknddi,nelddi,matdim

!---------------------------------------------------------------------
! assign new parameter file
!---------------------------------------------------------------------

	call pardef

!---------------------------------------------------------------------
! get names of basin and simulation
!---------------------------------------------------------------------

	call assnam(mode)

!---------------------------------------------------------------------
! if no bas file has to be opened we are done
!---------------------------------------------------------------------

	if( .not. btest(abs(mode),0) ) return

!---------------------------------------------------------------------
! open bas file
!---------------------------------------------------------------------

	call read_bas_file(nknddi,nelddi)

!---------------------------------------------------------------------
! set up boundary nodes
!---------------------------------------------------------------------

	!if(matdim.gt.0) call sp131k(matdim)

!---------------------------------------------------------------------
! end of routine
!---------------------------------------------------------------------

	end

!**************************************************************
!**************************************************************
!**************************************************************

	subroutine pardef

! reads default parameters nls and fnm

	implicit none

	character*80 apnnam

	logical bfirst
	save bfirst
	data bfirst /.true./

	if(bfirst) then
	  call nlsina
	  call fnminh

	  call getfnm('apnfil',apnnam)
	  call putfnm('apnnam',apnnam)

	  call read_apn_file(-1)

	  bfirst=.false.
	end if

	end

!**************************************************************

	subroutine read_apn_file(ivar)

	implicit none

	integer ivar

	integer nin
	integer ifileo

	nin=ifileo(-1,'apnstd.str','form','old')
	if( nin .le. 0 ) return

	call nlsa(nin,ivar,.false.)
	close(nin)

	end

!**************************************************************
!**************************************************************
!**************************************************************

	subroutine read_memory(simul,basin)

	implicit none

	character*(*) simul,basin

	integer imem
	character*80 memfil

	integer ifileo

	simul=' '
	basin=' '

	call getfnm('memfil',memfil)

	if(memfil.eq.' ') return

	imem=ifileo(-1,memfil,'form','old')

	if(imem.le.0) return

	read(imem,'(a)') simul
	read(imem,'(a)') basin

	close(imem)

	end

!**************************************************************

	subroutine write_memory(simul,basin)

	implicit none

	character*(*) simul,basin

	integer imem
	character*80 memfil

	integer ifileo

	call getfnm('memfil',memfil)

	if(memfil.eq.' ') return

	imem=ifileo(0,memfil,'form','new')

	if(imem.le.0) return

	write(imem,'(a)') trim(simul)
	write(imem,'(a)') trim(basin)

	close(imem)

	end

!**************************************************************

	subroutine read_bas_file(nknddi,nelddi)

! opens and reads basin file

	use basin

	implicit none

	integer nknddi,nelddi

	integer iunit
	character*80 basnew

	integer idefbas

	character*80 basold
	save basold
	data basold / ' ' /

	call getfnm('basnam',basnew)

	if( basnew .eq. ' ' ) return
	if( basnew .eq. basold ) return

	iunit=idefbas(basnew,'old')
	if(iunit.le.0) then
	  stop 'error stop read_bas_file: error reading basin'
	end if

	call basin_read(iunit)

	close(iunit)

	basold=basnew

	call bas_info

	end

!**************************************************************

