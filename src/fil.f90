
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998,2001,2008-2010,2012,2014-2015  Georg Umgiesser
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

! file opening routines
!
! contents :
!
! function ifileq(iunit,quest,format,status)	asks for filename and opens file
! function ifileo(iunit,file,format,status)	opens file
!
! subroutine useunit(iumax)                     writes usage of units and file
! subroutine filna(iunit,name)			inquires file name of iunit
! function filex(name)				inquires if file name exists
! function ifilun(name)				inquires unit of file name
! subroutine filext(file,ext)			appends file extension
!
! revision log :
!
! 02.05.1998	ggu	new format binary for ifileo
! 07.05.1998	ggu	new subroutine filext
! 21.05.1998	ggu	unit 0 is ok in ifileo()
! 29.06.1998	ggu	filna revisited
! 10.12.2001	ggu	filna still revisited, new routine useunit
! 31.03.2008	ggu	in ifileo remember last unit number
! 09.01.2009	ggu	re-formatted, safeguard units 5 and 6
! 23.03.2010	ggu	changed v6.1.1
! 26.03.2012	ggu	if opening already open file -> better error message
! 25.10.2012	ggu	changed VERS_6_1_59
! 28.04.2014	ggu	if iunit < 0 -> does not complain
! 30.10.2014	ggu	changed VERS_7_0_4
! 30.04.2015	ggu	changed VERS_7_1_9
! 26.04.2018	ggu	changed VERS_7_5_46
! 16.10.2018	ggu	changed VERS_7_5_50
! 16.02.2019	ggu	changed VERS_7_5_60
! 03.03.2020	ggu	better error messaging, if error return -1
! 30.03.2021	ggu	simplified find_unit()
! 20.04.2021	ggu	new routines for flushing and syncing
! 21.04.2021	clr	calling compiler-specific file descriptor function
!
!********************************************************

	function ifileq(iunit,quest,format,status)

! asks for filename and opens file
!
! iunit		unit number to be opened
! quest		text to be written to terminal
! format	f*ormatted or u*nformatted
! status	o*ld, n*ew, s*cratch  or u*nknown
! ifileq	effective unit number the file is opened
!		-1 if no file is opened, 0 if no name is given

	implicit none

	integer ifileq
	integer iunit
	character*(*) quest,format,status

	character*80  name

	integer ifileo

	write(6,*) quest

	read(5,'(a)') name

	if(name.eq.' ') then
		ifileq=0
	else
		ifileq=ifileo(iunit,name,format,status)
	end if

	end

!***************************************************************

	function ifileo(iunit,file,format,status)

! opens file
!
! iunit		unit number to be opened
! file		file name to be opened
! format	f*ormatted, u*nformatted or b*inary
! status	o*ld, n*ew, s*cratch  or u*nknown
! ifileo	effective unit number the file is opened
!
! routine tries to open file at iunit. If this is not possible
! ...higher unit numbers are tested and the file is eventually
! ...opened at one of these units. The effective unit number
! ...is passed back in ifileo. -1 is returned if no file is opened.
! use iunit <= 0 to use standard unit number
! iunit < 0 does not complain if not existing

	implicit none

	integer ifileo
	integer iunit
	character*(*) file,format,status

	integer iu,nfile,ios,iuaux
	logical found,error,opened,ex,od,bquiet
	character*1 cf,cs
	character*15 form,stat,access

	external ichanm0
	integer ichanm0,ifilun

        integer iustd
        save iustd
        data iustd / 20 /

!----------------------------------------------------------------
! initialize basic things
!----------------------------------------------------------------

	ifileo=-1
	iu=iunit
	bquiet = ( iu < 0 )		!if unit negative do not complain
	if( iu .le. 0 ) iu = iustd

!----------------------------------------------------------------
! check key words
!----------------------------------------------------------------

	cf=format(1:1)
	cs=status(1:1)

	access = 'sequential'

	if( cf .eq. 'f' .or. cf .eq. 'F' ) then
		form='formatted'
	else if( cf .eq. 'u' .or. cf .eq. 'U' ) then
		form='unformatted'
	else if( cf .eq. 'b' .or. cf .eq. 'B' ) then
		form='unformatted'
!lahey#		access = 'transparent'
	else
		write(6,*) 'format keyword not recognized :',format
		return
	end if

	if( cs .eq. 'o' .or. cs .eq. 'O' ) then
		stat='old'
	else if( cs .eq. 'n' .or. cs .eq. 'N' ) then
! for VAX      change to	stat='new'
! for DOS/UNIX change to	stat='unknown'
!		stat='new'
		stat='unknown'
	else if( cs .eq. 's' .or. cs .eq. 'S' ) then
		stat='scratch'
	else if( cs .eq. 'u' .or. cs .eq. 'U' ) then
		stat='unknown'
	else
		write(6,*) 'status keyword not recognized :',status
		return
	end if

!----------------------------------------------------------------
! check if file exists (in case of status=old)
!----------------------------------------------------------------

	inquire(file=file,exist=ex)
	if(.not.ex.and.stat.eq.'old') then
	  if( .not. bquiet ) then
	    write(6,*) 'file does not exist: ',trim(file)
	  end if
	  return
	end if

!----------------------------------------------------------------
! find unit where to open file
!----------------------------------------------------------------

	call find_unit(iu)
	found=iu.gt.0

!----------------------------------------------------------------
! open file and check error
!----------------------------------------------------------------

	if(found) then
	  !write(6,*) iu,file,form,stat,access
      open( unit=iu , file=file, form=form , status=stat, access=access, iostat=ios )
	  if(ios.ne.0) then
	    nfile=ichanm0(file)
	    if(nfile.le.0) nfile=1
	    write(6,*) 'error opening file : ',file(1:nfile)
	    write(6,*) 'unit : ',iu,'  iostat : ',ios
	    write(6,*) 'error : ',mod(ios,256)
	    inquire(file=file,opened=opened)
	    if( opened ) then
	      write(6,*) '...the file is already open...'
	      write(6,*) 'If you are using gfortran or pgf90'
	      write(6,*) 'please remember that you can open'
	      write(6,*) 'a file only once. You will have to'
	      write(6,*) 'copy the file to files with different'
	      write(6,*) 'names and open these files instead'
	      iuaux = ifilun(file)
	      write(6,*) iuaux,'  ',trim(file)
	    end if
	    ifileo=-1
	  else
		rewind(iu)
		ifileo=iu
                if( iunit .le. 0 ) iustd=iu ! remember for next time
                !write(10,*) 'ifileo: ',iu,iunit,'  ',file(1:40)
	  end if
	end if

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!*******************************************************

	subroutine useunit(iumax)

! writes usage of units and file names to stdout
!
! iumax		maximum unit number to try

	implicit none

	integer iumax

	integer iu,ium
	character*70 name

        ium = iumax
        if( ium .eq. 0 ) ium = 100

        write(6,*) 'unit usage:'
        do iu=1,ium
          call filna(iu,name)
          if( name .ne. ' ' ) then
                  write(6,'(i5,2x,a70)') iu,name
          end if
        end do

	end

!*******************************************************

	subroutine filna(iunit,name)

! inquires file name of unit number iunit
!
! iunit		unit number
! name		file name (return value)

	implicit none

	integer iunit
	character*(*) name

	logical btest

	name = ' '

	inquire(iunit,opened=btest)
	if( .not. btest ) return

	inquire(iunit,named=btest)
        name = '(no name available)'
	if( .not. btest ) return

	name = ' '
	inquire(unit=iunit,name=name)

	end

!*******************************************************

	function filex(name)

! inquires if file name exists
!
! name		file name
! filex		.true. if file exists

	implicit none

	logical filex
	character*(*) name

	inquire(file=name,exist=filex)

	end

!*******************************************************

	function ifilun(name)
!
! inquires unit of file name if opened
!
! name		file name
! ifilun	unit if file is opened, else 0

	implicit none

	integer ifilun
	character*(*) name

	logical open
	integer iunit

	inquire(file=name,opened=open)
	if(open) then
		inquire(file=name,number=iunit)
	else
		iunit=0
	end if

	ifilun=iunit

	end

!*******************************************************

	subroutine filext(file,ext)

! appends file extension if not there

	implicit none

	character*(*) file,ext

	integer nf,ne,length
	integer nef,nel,nff,nfl
	integer ichanm0

	length = len(file)

	nf = ichanm0(file)
	ne = ichanm0(ext)

	nef = 1
	nel = ne
	if( ext(1:1) .eq. '.' ) then	!in case dot is in extension
	  nef = 2
	  ne = ne - 1
	end if

	if( ne .le. 0 ) return
	if( nf+ne+1 .gt. length ) return

!	now we know that there is enough room for extension

	nff = nf - ne + 1
	nfl = nf

	if( nff .le. 0 .or. file(nff:nfl) .ne. ext(nef:nel) ) then
	  file(nfl+1:) = '.' // ext(nef:nel)
	end if

	end

!*******************************************************

	subroutine find_unit(iunit)

! finds unit to open file - starts to search from iunit
! on return iunit is either the next unit available
! or it is 0 which means there was an error

	implicit none

	integer iunit

	logical opened
	integer, save :: iumax = 1000
	integer iu

	iu = iunit
	if( iu .le. 0 ) iu = 20			!set standard unit

	do
		if( iu .eq. 5 ) iu = 7		!safeguard units 5 and 6
		if( iu > iumax ) then
			write(6,*) 'no unit available to open file'
                        write(6,*) 'unit tried: ',iu
                        call useunit(iu-1)
			iunit = 0
			return
		end if
		inquire(iu,opened=opened)
		if( .not. opened ) exit
		iu=iu+1
	end do

	iunit = iu

	end

!*******************************************************

	subroutine file_flush(iu)

	implicit none

	integer iu

	flush(iu)

	end

!*******************************************************

	subroutine file_sync(iu)

#if defined(__INTEL_COMPILER)
        use ifposix
#endif
        implicit none

        interface
              function fsync (fd) bind(c,name="fsync")
              use iso_c_binding, only: c_int
                integer(c_int), value :: fd
                integer(c_int) :: fsync
              end function fsync
        end interface

        integer iu

        integer ierr,file_descriptor
        character*80 name

        flush(iu)

#if defined(__GFORTRAN__)
	file_descriptor = fnum(iu)
#elif defined(__INTEL_COMPILER)
        call pxffileno (iu,file_descriptor,ierr)
        if( ierr /= 0 ) file_descriptor = -1   
#else
	file_descriptor = -1
#endif

	if( file_descriptor == -1 ) return

        ierr = fsync(file_descriptor)

        if( ierr /= 0 ) then   
          call filna(iu,name)
          write(6,*) iu,'  ',trim(name)
          stop 'error stop file_sync: error syncing file'
        end if

	end

!*******************************************************
!*******************************************************
!*******************************************************
! stub in order to be independent of subsss.f
!*******************************************************
!*******************************************************
!*******************************************************

        function ichanm0(line)

! computes length of line without trailing blanks
!
! line          line of text
! ichanm0       length of line (return value)
!               ... 0 : line is all blank

	implicit none

	integer ichanm0
        character*(*) line

	integer ndim,i

        character*1 blank,tab
        character*1 char
        data blank /' '/

        tab=char(9)

        ndim=len(line)

        do i=ndim,1,-1
          if(line(i:i).ne.blank.and.line(i:i).ne.tab) goto 1
        end do

    1   continue
        ichanm0=i

        end

!*******************************************************

