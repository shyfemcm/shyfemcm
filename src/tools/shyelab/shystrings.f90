
!--------------------------------------------------------------------------
!
!    Copyright (C) 2024  Georg Umgiesser
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

! revision log :
!
! 21.09.2024	ggu	written from scratch
!
!****************************************************************

        program shystrings

! show list of strings

        use shyfem_strings
        use clo

        implicit none

        character(80) :: string

        logical :: berror,bwrite
        logical :: bquiet,bsilent,bverbose
	integer :: is,ioff,ivar,id,idlast,iwrite
	real    :: f(10)
	character*10 :: short
	character*28 :: full,header

	logical is_integer
	integer istof

!---------------------------------------------------------------
! command line options
!---------------------------------------------------------------

        call clo_init('shystrings','[ivar|name]','1.0')
        call clo_add_info('shows strings and variable numbers')

        call clo_add_sep('general options')

        call clo_add_option('verbose',.false.,'write extra information')
        call clo_add_option('quiet',.false.,'be quiet in execution')
        call clo_add_option('silent',.false.,'do not write anything')

        call clo_add_sep('usage:')
        call clo_add_com('  shystrings ivar     shows variable ivar')
        call clo_add_com('  shystrings string   ' // &
     &				'looks for string in strings')
        call clo_parse_options

        call clo_get_option('verbose',bverbose)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)

	if( bsilent ) bquiet = .true.
        call shyfem_set_short_copyright(bquiet)
        if( .not. bsilent ) then
	  call shyfem_copyright('shystrings - show strings')
        end if

	call populate_strings

	header = ' ivar   short        full'

!---------------------------------------------------------------
! read files and set up parameters
!---------------------------------------------------------------

        call clo_check_files(0)
        call clo_get_file(1,string)

	if( bverbose ) write(6,*) 'looking for ',trim(string)

	if( is_integer(string) ) then
	  ioff = 1
	  is = istof(string,f,ioff)
	  ivar = nint(f(1))
	  if( bverbose ) write(6,*) 'looking for variable ',ivar
	  id = strings_get_id_by_ivar(ivar)
	  if( id <= 0 ) then
	    write(6,*) 'no such variable number: ',ivar
	  else
	    call strings_get_ivar_and_names(id,ivar,short,full)
	    write(6,'(a)') header
	    write(6,1000) ivar,'   ',short,'   ',full
	  end if
	else
	  if( bverbose ) then
	    write(6,*) 'looking for string ',trim(string),' in variable names'
	  end if
	  idlast = strings_filling()
	  iwrite = 0
	  do id=1,idlast
	    call strings_get_ivar_and_names(id,ivar,short,full)
	    bwrite = .false.
	    if( index(short,trim(string)) > 0 ) bwrite = .true.
	    if( index(full,trim(string)) > 0 ) bwrite = .true.
	    if( bwrite ) then
	      if( iwrite == 0 ) write(6,'(a)') header
	      write(6,1000) ivar,'   ',short,'   ',full
	      iwrite = iwrite + 1
	    end if
	  end do
	  if( iwrite == 0 ) write(6,*) 'cannot find string ',trim(string)
	end if

 1000	format(i5,6a)

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

        end program

!***************************************************************

