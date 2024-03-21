
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998,2003-2005,2010,2014-2019  Georg Umgiesser
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

! OUS file administration routines
!
! contents :
!
! subroutine wrousa
!
! revision log :
!
! 26.01.1998	ggu	$$ITMOUT - adjust itmout for first write
! 01.09.2003	ggu	new routine wrousa
! 02.09.2003	ggu	bug fix in wrousa: save nbout
! 25.11.2004	ggu	in new file subousa.f
! 18.05.2005	ggu	initial itmout is changed
! 23.03.2010	ggu	changed v6.1.1
! 20.01.2014	ggu	new calls for ous writing implemented
! 28.01.2014	ggu	changed VERS_6_1_71
! 26.11.2014	ggu	changed VERS_7_0_7
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 15.10.2015	ggu	added new calls for shy file format
! 28.04.2016	ggu	changed VERS_7_5_9
! 09.09.2016	ggu	changed VERS_7_5_17
! 30.09.2016	ggu	changed VERS_7_5_18
! 05.12.2017	ggu	changed VERS_7_5_39
! 03.04.2018	ggu	changed VERS_7_5_43
! 11.04.2018	ggu	mpi version is ready and working
! 19.04.2018	ggu	changed VERS_7_5_45
! 11.05.2018	ggu	mpi version working also for zeta levels
! 16.02.2019	ggu	changed VERS_7_5_60
!
!********************************************************

	subroutine wrousa

! writes and administers ous file

	use mod_trace_point
	use mod_depth
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shyfile
	use shympi
	use simul

	implicit none

	logical bdebug
	integer itmout,ierr
	real href,hzoff
	integer date,time
	character*80 title,femver

	integer iround,ideffi
        integer wfout,wrout
	integer nvar,ftype,id
	double precision dtime
	character*80 file

	real getpar
	double precision dgetpar
	integer ifemop
	logical has_output_d,next_output_d

	integer idtout,itout
	integer icall,nbout,nvers
	double precision, save :: da_out(4)
	save idtout,itout
	save icall,nvers,nbout
	data icall,nvers,nbout /0,2,0/

	bdebug = .true.
	bdebug = .false.

	if( icall .eq. -1 ) return

	if( icall .eq. 0 ) then
		call trace_point('init hydro_file in wrousa')
		da_out = 0.
		call init_output_d('itmout','idtout',da_out)

		if( .not. has_output_d(da_out) ) icall = -1
		if( icall .eq. -1 ) return

		if( has_output_d(da_out) ) then
		  call shyfem_init_hydro_file('hydro',.false.,id)
		  da_out(4) = id
		end if
	end if

	icall = icall + 1

	if( next_output_d(da_out) ) then
	  call trace_point('writing hydro_records in wrousa')
	  id = nint(da_out(4))
	  call get_act_dtime(dtime)
	  call shy_write_hydro_records(id,dtime,nlvdi,znv,zenv &
     &					,utlnv,vtlnv)
	end if

	call trace_point('finished wrousa')

	return
   77   continue
	write(6,*) 'Error opening OUS file :'
	stop 'error stop : wrousa'
   78   continue
	write(6,*) 'Error writing first header of OUS file'
	write(6,*) 'unit,err :',nbout,ierr
	stop 'error stop : wrousa'
   75   continue
	write(6,*) 'Error writing second header of OUS file'
	write(6,*) 'unit,err :',nbout,ierr
	stop 'error stop : wrousa'
   79   continue
	write(6,*) 'Error writing data record of OUS file'
	write(6,*) 'unit,err :',nbout,ierr
	stop 'error stop : wrousa'
	end

!********************************************************

