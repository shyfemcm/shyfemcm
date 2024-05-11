
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
! could eliminate de-scrambling of iknot ->
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
	logical, save :: bnepart

	logical, save :: binfo
	logical, save :: bquiet
	logical, save :: bsilent

	real, save :: eps_area = 1.e-4

	character*80, save :: grdpart

!==========================================================
	end module mod_shypre
!==========================================================

!**********************************************************

        program shypre

	use mod_shypre

	implicit none

        character*80 grdfile

!--------------------------------------------------------
! get name of basin and parse command line options
!--------------------------------------------------------

	call shypre_init(grdfile)

	call shypre_sub(grdfile,bquiet,bsilent,bauto,binfo &
     &				,bopti,bnepart,eps_area,grdpart)

	end

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
        call clo_add_option('opti',.false.,'optimize bandwidth' &
     &                             ,'do not optimize bandwidth (default)')
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
        call clo_get_option('manual',bmanual)

        call clo_get_option('nepart',bnepart)
        call clo_get_option('partition',grdpart)

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

