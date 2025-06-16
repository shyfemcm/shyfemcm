
!--------------------------------------------------------------------------
!
!    Copyright (C) 2004,2006,2008-2012,2014-2015,2014-2015  Georg Umgiesser
!    Copyright (C) 2017-2020  Georg Umgiesser
!    Copyright (C) 2012  Andrea Cucco
!    Copyright (C) 2014  Christian Ferrarin
!    Copyright (C) 2018  Marco Bajo
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

! restart routines
!
! contents :
!
! subroutine wrrst(it,iunit)	writes one record of restart data
! subroutine rdrst(itrst,iunit)	reads one record of restart data
!
! subroutine skip_rst(iunit,atime,it,nvers,nrec,nkn,nel,nlv,iflag,ierr)
!				returns info on record in restart file
!
! revision log :
!
! 02.09.2004	ggu	started with routine rdrst()
! 18.10.2006	ggu	included hm3v in restart file -> nvers = 4
! 13.06.2008	ggu	new version 5 -> S/T/rho
! 09.01.2009	ggu	bugfix in inirst - file opened with status=new
! 23.03.2009	ggu	bugfix in admrst - itmrst=itanf if itmrst<itanf
! 07.05.2009	ggu	new parameter ityrst
! 29.05.2009	ggu	use closest record for restart (if ityrst=2)
! 13.11.2009	ggu	keep track of restart: /rstrst/ and has_restart()
! 27.11.2009	ggu	deal with empty file, rdrst() restructured
! 19.01.2010	ggu	initialize also conz, has_restart() is function
! 11.03.2010	ggu	write also vertical velocity
! 23.03.2010	ggu	changed v6.1.1
! 14.07.2011	ggu	changed VERS_6_1_27
! 10.02.2012	ggu	write only last record, restart from last record
! 16.02.2012	aac	write also ecological varibles
! 01.06.2012	ggu	changed VERS_6_1_53
! 28.08.2012	aac	bug fix for restart time = -1 (rdrst_record)
! 07.03.2014	ggu	changed VERS_6_1_72
! 26.11.2014	ggu	changed VERS_7_0_7
! 11.12.2014	ccf	bug fix for atime
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 05.06.2015	ggu	changed VERS_7_1_12
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 20.10.2015	ggu	bug fix to get correct restart time (itanf)
! 30.10.2015	ggu	new names, restructured
! 05.11.2015	ggu	changed VERS_7_3_12
! 09.11.2015	ggu	changed VERS_7_3_13
! 30.11.2015	ggu	allocate cnv/conzv before read
! 16.12.2015	ggu	changed VERS_7_3_16
! 09.05.2017	ggu	changed VERS_7_5_26
! 06.07.2017	ggu	saved hlv
! 04.11.2017	ggu	changed VERS_7_5_34
! 30.05.2018	ggu	some time values now in double
! 31.05.2018	ggu	new version (11), all time values in double
! 28.06.2018	mbj	bug fix for version 11
! 06.07.2018	ggu	changed VERS_7_5_48
! 31.08.2018	ggu	changed VERS_7_5_49
! 25.10.2018	ggu	bug fix with finding desired record
! 16.02.2019	ggu	changed VERS_7_5_60
! 03.05.2019	ggu	return iconz in rst_skip_record()
! 03.05.2019	ggu	new routine to check if rst file (rst_is_rst_file)
! 21.05.2019	ggu	changed VERS_7_5_62
! 09.03.2020	ggu	restart for mercury
! 20.03.2020	ggu	completely restructured
! 17.11.2020	ggu	new routine rst_get_hlv() and array hlvrst
! 18.11.2020	ggu	new version 13 (write ilhv,ilhkv)
! 18.11.2020	ggu	new version 14 (only write vertical for nlv>1)
! 30.03.2021	ggu	new routine init_old_vars()
! 01.04.2021	ggu	save turbulence in restart
! 15.07.2021	ggu	store flags for what has been read
! 11.10.2022	ggu	mpi-version - helper functions defined
! 13.10.2022	ggu	mpi-version finished (no gotm, conz, etc.. yet)
! 14.10.2022	ggu	mpi-version final (gotm, conz)
! 28.04.2023    ggu     update function calls for belem
! 02.05.2023    ggu     fix mpi bug for nlv==1
! 26.01.2025    ggu     fix mpi bug for ibfm /= 0 (only master writes)
! 08.03.2025    ggu     extracted mod_restart in its own file
! 09.03.2025    ggu     call shympi_barrier after finishing restart file
!
! notes :
!
!  |--------|---------|------------|--------------|--------------|
!  | ityrst | no file | file empty ! itrec/=itrst ! itrec==itrst !
!  |--------|---------|------------|--------------|--------------|
!  |   0    |  error  |   error    |     error    |    warm      |
!  |--------|---------|------------|--------------|--------------|
!  |   1    |  cold   |   error    |     error    |    warm      |
!  |--------|---------|------------|--------------|--------------|
!  |   2    |  cold   |   cold     |     warm     |    warm      |
!  |--------|---------|------------|--------------|--------------|
!
! versions :
!
! <3	not supported
! 3	write hydro
! 4	write hm3v
! 5	write saltv,tempv,rhov
! 6	write conv (tracer)
! 7	write wlnv (vertical velocity)
! 8	write ecological variables
! 9	write date,time
! 10	write hlv
! 11	write atime, write idfile (regular header)
! 12	mercury restart
! 13	write ilhv and ilhkv
! 14	write vertical only for nlv > 1
! 15	write gotm arrays
! 16	adapted for mpi
! 17	write bfm restart
!
!	integer, save :: id_hydro_rst = 1	!1		hydro
!	integer, save :: id_depth_rst = 2	!10		depth
!	integer, save :: id_barcl_rst = 3	!100		t/s/rho
!	integer, save :: id_conz_rst  = 4	!1000		tracer
!	integer, save :: id_wvert_rst = 5	!10000		vertical vel.
!	integer, save :: id_eco_rst   = 6	!100000		ecology
!	integer, save :: id_merc_rst  = 7	!1000000	mercury
!	integer, save :: id_gotm_rst  = 8	!10000000	gotm
!	integer, save :: id_bfm_rst   = 9	!100000000	bfm
!
!*********************************************************************

! mod_restart is defined in own file

!*********************************************************************

        subroutine rst_perform_restart

! reads and initializes values from restart

	use mod_restart
	use shympi

        implicit none

	logical blast
	logical bavail,bwant,buse
        integer iunit,ierr,ityrst,id
	integer date,time
	double precision atime0,atime,atrst,ditrst
	double precision dtanf,dtend
        character*80 name
        character*20 aline

        real getpar
	double precision dgetpar
        integer ifileo
	logical rst_is_set

!-----------------------------------------------------------------
! get parameters
!-----------------------------------------------------------------

	call convert_date_d('itrst',ditrst)
        ityrst = nint(dgetpar('ityrst'))
	iflag_want_rst = nint(getpar('flgrst'))
        call getfnm('restrt',name)
        if(name.eq.' ') return

        date = nint(dgetpar('date'))
        time = nint(dgetpar('time'))
	blast = ( ditrst == -1. )

	call dts_to_abs_time(date,time,atime0)
	atrst = atime0 + ditrst
	if( blast ) atrst = -1.

!-----------------------------------------------------------------
! name of restart file given -> open and read
!-----------------------------------------------------------------

        write(6,*) '---------------------------------------------'
        write(6,*) '... performing restart from file:'
        write(6,*) trim(name)
        write(6,*) '---------------------------------------------'

        iunit = ifileo(1,name,'unformatted','old')
        if( iunit .le. 0 ) then
          if( ityrst .le. 0 ) goto 98
          write(6,*) '*** Cannot find restart file ...'
          write(6,*) '*** Continuing with cold start...'
          return
        end if

	atime = atrst
        call rst_read_restart_file(iunit,atime,iflag_avail_rst,ierr)

        if( ierr .gt. 0 ) then
          if( ityrst .le. 1 ) goto 97
          if( ierr .eq. 95 ) then	!file is empty
            write(6,*) '*** No data in restart file ...'
            write(6,*) '*** Continuing with cold start...'
            return
	  else if( ierr == 94 ) then
            write(6,*) '*** hlv not compatible'
            stop 'error stop rst_perform_restart: hlv'
          end if
          write(6,*) '*** Another time record is used for restart'
	  call dts_format_abs_time(atrst,aline)
          write(6,*) '*** Looking for time = ',atrst,aline
	  call dts_format_abs_time(atime,aline)
          write(6,*) '*** Finding time = ',atime,aline
          write(6,*) '*** Continuing with hot start...'
	else
	  call dts_format_abs_time(atime,aline)
	  write(6,*) 'restart time found: ',aline
        end if

	close(iunit)

	if( blast ) then	!reset initial time
	  call dts_format_abs_time(atime,aline)
	  write(6,*) 'setting new initial time: ',aline
	  dtanf = atime-atime0
	  call dputpar('itanf',dtanf)
	  call putfnm('itanf',aline)
	end if

	dtanf = dgetpar('itanf')
	dtend = dgetpar('itend')

        write(6,*) '---------------------------------------------'
        write(6,*) 'A restart has been performed'
	if( blast ) then
          write(6,*) ' requested restart time: last record'
	else
	  call dts_format_abs_time(atrst,aline)
          write(6,*) ' requested restart time = ',aline
	end if
	call dts_format_abs_time(atime,aline)
        write(6,*) ' used restart time =      ',aline
        write(6,*) ' nvers = ',nvers_rst
	write(6,*) '         id bwant bavail buse   description'
	do id=1,nidmax
	  bwant = rst_is_set(id,iflag_want_rst)
	  bavail = rst_is_set(id,iflag_avail_rst)
	  buse = bwant .and. bavail
	  write(6,*) id,'   ',bwant,'   ',bavail,'   ', buse, '  ',descript_rst(id)
	end do
	call dts_format_abs_time(atime0+dtanf,aline)
        write(6,*) ' itanf = ',aline
	call dts_format_abs_time(atime0+dtend,aline)
        write(6,*) ' itend = ',aline
        write(6,*) '---------------------------------------------'

	call init_old_vars	!initializes also old values

	bok_rst = .true.

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

        return
   97   continue
	call dts_format_abs_time(atrst,aline)
        write(6,*) 'no record found for time = ',aline
        stop 'error stop rst_perform_restart: Cannot find time record'
   98   continue
        write(6,*) 'no such file : ',name
        stop 'error stop rst_perform_restart: Cannot read restart file'
        end

!*******************************************************************

        subroutine rst_read_restart_file(iunit,atrst,iflag,ierr)

! reads restart file until it finds atrst

        implicit none

        integer iunit		!unit to read from
        double precision atrst	!absolute time
	integer iflag		!flag of records that are available
        integer ierr            !error code - different from 0 if error

        integer ii,l,ie,k
        integer irec
        logical bloop,blast,bnext
        double precision atime,alast
	character*20 aline

        irec = 0
	ierr = 0
	alast = 0.
	blast = atrst .eq. -1		! take last record

        do
          call rst_read_record(iunit,atime,iflag,ierr)
          if( ierr .gt. 0 ) goto 94
          if( ierr .lt. 0 ) exit
          irec = irec + 1
	  if( .not. blast .and. atime .ge. atrst ) exit
	  alast = atime
        end do

	if( irec .gt. 0 ) then
	  if( blast ) then
	    ierr = 0
	    atrst = alast
	    return
          else if( atime .eq. atrst ) then
	    return
	  end if
	end if

        if( ierr .ne. 0 ) then          !EOF
          if( irec .eq. 0 ) then        !no data found
                goto 95
          else                          !last record in file has smaller time
                goto 97
          end if
        else                            !read past desired time
                goto 97
        end if

        return
   94   continue
        write(6,*) 'error reading restart file: ',ierr
        ierr = 94
	return
   95   continue
        write(6,*) 'reading restart file... '
        write(6,*) 'no records in file (file is empty)'
        ierr = 95
        return
   97   continue
        write(6,*) 'reading restart file... '
	call dts_format_abs_time(atime,aline)
        write(6,*) 'last record read at time = ',atime,aline
	call dts_format_abs_time(atrst,aline)
        write(6,*) 'no record found for time = ',atrst,aline
        ierr = 97
        return
        end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	function rst_has_restart(id)

! gives indication if data from restart is available
!
! id indicates what information is requested

	use mod_restart

	implicit none

	logical rst_has_restart
	integer id

	logical rst_is_set

	if( .not. bok_rst ) then		!no restart file read
	  rst_has_restart = .false.
	else if( id .eq. 0 ) then		!general restart data available
	  rst_has_restart = .true.
	else if( id .ge. 1 .and. id .le. nidmax ) then
	  rst_has_restart = rst_is_set(id,iflag_avail_rst)
	else
	  rst_has_restart = .false.
	end if

	end

!*******************************************************************

	function rst_want_restart(id)

! see if restart for a specific variable is wanted
!
! if id < 0		restart is always wanted
!
! example: iflag = 1011 means that for id 1,2,4 function is true, else false

	use mod_restart

	implicit none

	logical rst_want_restart
	integer id		!number of feature desired

	logical rst_is_set

	rst_want_restart = .true.
	if( iflag_want_rst < 0 ) return

	rst_want_restart = rst_is_set(id,iflag_want_rst)

	end

!*******************************************************************

	function rst_use_restart(id)

! see if restart for a specific variable has been used (available and wanted)
!
! if id < 0		restart is always wanted
!
! example: iflag = 1011 means that for id 1,2,4 function is true, else false

	use mod_restart

	implicit none

	logical rst_use_restart
	integer id		!number of feature desired

	logical rst_has_restart,rst_want_restart

	rst_use_restart = rst_has_restart(id) .and. rst_want_restart(id)

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

        subroutine rst_write_restart

! administers writing of restart file

	use shympi

        implicit none

	logical, parameter :: bdebug = .true.
	integer ierr
        integer iunit
	double precision dtmrst,ddtrst
	double precision atime
	double precision dtanf,dtend

        real getpar
        double precision dgetpar
        integer ifemop
	!integer fsync
	logical has_output_d,next_output_d

	logical, save :: bonce
	double precision, save :: da_out(4) = 0
        integer, save :: icall = 0

        if( icall .le. -1 ) return

!-----------------------------------------------------
! initializing
!-----------------------------------------------------

        if( icall .eq. 0 ) then

          call convert_date_d('itmrst',dtmrst)
          call convert_time_d('idtrst',ddtrst)

	  if( ddtrst .lt. 0. ) then	!only last record saved
	    if( ddtrst .eq. -1. ) then	!only at the end of the simulation
	      dtmrst = -1.
	      call get_first_dtime(dtanf)
	      call get_last_dtime(dtend)
	      ddtrst = -(dtend-dtanf)
	    end if
	    bonce = .true.
	    ddtrst = -ddtrst
	  else
	    bonce = .false.
	  end if

          icall = -1
          call set_output_frequency_d(dtmrst,ddtrst,da_out)
	  !call info_output_d('restart',da_out)
	  !stop
	  !call increase_output_d(da_out)
          if( .not. has_output_d(da_out) ) return
          icall = 1

	  if( .not. bonce ) then
            iunit = ifemop('.rst','unformatted','new')
            if( iunit .le. 0 ) goto 98
	    da_out(4) = iunit
	  end if

        end if

!-----------------------------------------------------
! normal call and writing
!-----------------------------------------------------

        if( .not. next_output_d(da_out) ) return

	call get_absolute_act_time(atime)

	!call check_values	!be sure values of restart are ok

	if( bonce ) then
	  if( bdebug .and. bmpi_master ) then
	    write(6,*) 'writing single restart record'
	  end if
          iunit = ifemop('.rst','unformatted','new')
          if( iunit .le. 0 ) goto 98
          call rst_write_record(atime,iunit)
	  close(iunit)
	else
	  if( bdebug .and. bmpi_master ) then
	    write(6,*) 'writing restart record'
	  end if
	  iunit = nint(da_out(4))
          call rst_write_record(atime,iunit)
	  call file_sync(iunit)
	end if

!-----------------------------------------------------
! end of routine
!-----------------------------------------------------

        return
   98   continue
        stop 'error stop rst_write_restart: Cannot open rst file'
        end

!*******************************************************************
!*******************************************************************
!*******************************************************************

        subroutine rst_write_record(atime,iunit)

! writes one record of restart data

	use mod_geom_dynamic
	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use mod_restart
	use levels, only : nlvdi,nlv,hlv,ilhv,ilhkv
	use basin
	use shympi

        implicit none

        double precision atime
        integer iunit

	logical b3d
        integer it
        integer ii,l,ie,k,i
	integer ibarcl,iconz,ibio,ibfm,ieco,imerc,iturb
        integer nvers
	integer date,time

	logical, parameter :: be = .true.
	logical, parameter :: bn = .false.

	real getpar
        double precision dgetpar

        nvers = nvmax
	b3d = nlv_global > 1
	bmaster = shympi_is_master()

	ibarcl = nint(getpar('ibarcl'))
	iconz = nint(getpar('iconz'))
	ibio = nint(getpar('ibio'))
	ibfm = nint(getpar('ibfm'))
	imerc = nint(getpar('imerc'))
	iturb = nint(getpar('iturb'))
        date = nint(dgetpar('date'))
        time = nint(dgetpar('time'))

	ieco = ibio + ibfm
	ieco = ibio

	if( bmaster ) then
          write(iunit) idfrst,nvers,1
          write(iunit) date,time
          write(iunit) atime
          write(iunit) nkn_global,nel_global,nlv_global
	end if

	if( b3d ) then
          if( bmaster ) write(iunit) hlv_global
	  call restart_write_value(iunit,be,ilhv)
	  call restart_write_value(iunit,bn,ilhkv)
	end if

	call restart_write_value(iunit,be,iwegv)
	call restart_write_value(iunit,bn,znv)
	call restart_write_value(iunit,be,3,zenv)
	call restart_write_value(iunit,be,utlnv)
	call restart_write_value(iunit,be,vtlnv)

	call restart_write_value(iunit,be,3,hm3v)

	call restart_write_value(iunit,ibarcl)
	if( ibarcl .gt. 0 ) then
	  call restart_write_value(iunit,bn,saltv)
	  call restart_write_value(iunit,bn,tempv)
	  call restart_write_value(iunit,bn,rhov)
	end if

	call restart_write_value(iunit,iturb)
	if( iturb .eq. 1 ) then
	  call write_restart_gotm(iunit)
	end if
	
	call restart_write_value(iunit,iconz)
	if( iconz .gt. 0 ) then
	  call write_restart_conz(iunit,iconz)
	end if
	
	call restart_write_value(iunit,nlv_global-1)
	if( nlv_global .gt. 1 ) then
	  call restart_write_value(iunit,bn,wlnv)
	end if

	call restart_write_value(iunit,ieco)
	if( ieco .gt. 0 ) then
	  call write_restart_eco(iunit)
        end if

	call restart_write_value(iunit,imerc)
	if( imerc .gt. 0 ) then
	  call write_restart_mercury(iunit)
        end if

	call restart_write_value(iunit,ibfm)
	if( ibfm .gt. 0 ) then
	  call write_restart_bfm(iunit)
        end if

	if( bmaster ) flush(iunit)
	call shympi_barrier

        end

!*******************************************************************

	subroutine rst_skip_record(iunit,atime,nvers,nrec,nkn,nel,nlv,iflag,ierr)

! returns info on record in restart file and skips data records
!
! in iflag returns availability of specific data

	use mod_restart
	use shympi

	implicit none

	logical brewrite
	integer iuout
	integer iunit,nvers,nrec,nkn,nel,nlv,iflag,ierr
	double precision atime
	integer ibarcl,iconz,iwvert,ieco,imerc,iturb,ibfm
	integer idfile
	integer date,time,it,id

	integer, allocatable :: ival(:)
	real, allocatable :: rval2d(:)
	real, allocatable :: rval2d3(:)
	real, allocatable :: rval3d1(:)
	real, allocatable :: rval3d2(:)
	real, allocatable :: rval3d3(:)

	iuout = 66
	iuout = 0
	brewrite = iuout > 0

	read(iunit,end=2,err=3) idfile,nvers,nrec

        ierr = 0

	nvers_rst = nvers

	date = 0
	time = 0
	iflag = 0

        if( nvers >= 9 ) read(iunit) date,time
        if( nvers >= 11 ) read(iunit) atime

	if( nvers <= 10 ) then
	  it = idfile
	  atime = 0.
	  if( date > 0 ) call dts_to_abs_time(date,time,atime)
	  atime = atime + it
	else if( idfile /= idfrst ) then
	  goto 7
	end if

        read(iunit) nkn,nel,nlv

	allocate(ival(nel))
	allocate(rval2d(nkn))
	allocate(rval2d3(3*nel))
	allocate(rval3d1((nlv+1)*nel))
	allocate(rval3d2((nlv+1)*nel))
	allocate(rval3d3((nlv+1)*nel))

	if( brewrite ) then
	  write(iuout,*) 'rst: new time record'
          write(iuout,*) idfile,nvers,nrec
          write(iuout,*) date,time
          write(iuout,*) atime
          write(iuout,*) nkn,nel,nlv
	end if

	id = id_hydro_rst
	call rst_add_flag(id,iflag)

	call rst_read_vertical(iunit,nvers,nkn,nel,nlv)

	if( brewrite ) then
	  write(iuout,*) 'rst: vertical'
	  write(iuout,*) hlvrst
	  write(iuout,*) ilhrst
	  write(iuout,*) ilhkrst
	end if

	read(iunit) ival(1:nel)
	read(iunit) rval2d(1:nkn)
	read(iunit) rval2d3(1:3*nel)
	read(iunit) rval3d1(1:nlv*nel)
	read(iunit) rval3d2(1:nlv*nel)

	if( brewrite ) then
	  write(iuout,*) 'rst: iwegv'
	  write(iuout,*) ival(1:nel)
	  write(iuout,*) 'rst: znv'
	  write(iuout,*) rval2d(1:nkn)
	  write(iuout,*) 'rst: zenv'
	  write(iuout,*) rval2d3(1:3*nel)
	  write(iuout,*) 'rst: utlnv'
	  write(iuout,*) rval3d1(1:nlv*nel)
	  write(iuout,*) 'rst: vtlnv'
	  write(iuout,*) rval3d2(1:nlv*nel)
	end if

	if( nvers .ge. 4 ) then
	  id = id_depth_rst
	  call rst_add_flag(id,iflag)
	  read(iunit) rval2d3(1:3*nel)
	  if( brewrite ) then
	    write(iuout,*) 'rst: hm3v'
	    write(iuout,*) rval2d3(1:3*nel)
	  end if
	end if

	if( nvers .ge. 5 ) then
	  id = id_barcl_rst
	  read(iunit) ibarcl
	  if( ibarcl .gt. 0 ) then
	    call rst_add_flag(id,iflag)
	    read(iunit) rval3d1(1:nlv*nkn)
	    read(iunit) rval3d2(1:nlv*nkn)
	    read(iunit) rval3d3(1:nlv*nkn)
	  end if
	end if

	if( brewrite ) then
	  write(iuout,*) 'rst: saltv'
	  write(iuout,*) rval3d1(1:nlv*nkn)
	  write(iuout,*) 'rst: tempv'
	  write(iuout,*) rval3d2(1:nlv*nkn)
	  write(iuout,*) 'rst: rhov'
	  write(iuout,*) rval3d3(1:nlv*nkn)
	end if

	if( nvers .ge. 15 ) then
	  id = id_gotm_rst
	  read(iunit) iturb
	  if( iturb .eq. 1 ) then
	    call rst_add_flag(id,iflag)
	    call skip_restart_gotm(iunit)
	  end if
	end if

	if( nvers .ge. 6 ) then
	  id = id_conz_rst
	  read(iunit) iconz
	  if( iconz .gt. 0 ) then
	    call rst_add_flag(id,iflag)
	    call skip_restart_conz(iunit,iconz)
	  end if
	end if

	if( nvers .ge. 7 ) then
	  id = id_wvert_rst
	  read(iunit) iwvert
	  if( iwvert .gt. 0 ) then
	    call rst_add_flag(id,iflag)
	    read(iunit) rval3d1((nlv+1)*nkn)
	    if( brewrite ) then
	      write(iuout,*) 'rst: wlnv'
	      write(iuout,*) rval3d1((nlv+1)*nkn)
	    end if
	  end if
	end if

	if( nvers .ge. 8 ) then
	  id = id_eco_rst
          read(iunit) ieco
          if( ieco .gt. 0 ) then
	    call rst_add_flag(id,iflag)
	    call skip_restart_eco(iunit)
          end if
        end if

	if( nvers .ge. 12 ) then
	  id = id_merc_rst
          read(iunit) imerc
          if( imerc .gt. 0 ) then
	    call rst_add_flag(id,iflag)
	    call skip_restart_mercury(iunit)
          end if
        end if

	if( nvers .ge. 17 ) then
	  id = id_bfm_rst
          read(iunit) ibfm
          if( ibfm .gt. 0 ) then
	    call rst_add_flag(id,iflag)
	    call skip_restart_bfm(iunit)
          end if
        end if

	ierr = 0
	return

    2	continue
	ierr = -1
	return
    3	continue
	write(6,*) 'rst_skip_record: error in reading restart file'
	ierr = 1
	return
    7	continue
	write(6,*) 'idfile,idfrst: ',idfile,idfrst
	write(6,*) 'rst_skip_record: error in idfrst... no restart format'
	ierr = 7
	return
	end

!*******************************************************************

        subroutine rst_read_record(iunit,atime,iflag,ierr)

! reads one record of restart data
!
! iflag is returned, which indicates the available data in the file
! this can be different from the actually read data (if not wanted)

	use mod_geom_dynamic
	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use levels, only : nlvdi,nlv,hlv
	use basin
	use mod_restart
	use shympi

        implicit none

        integer iunit		!unti from which to read
        double precision atime	!absolute time
        integer iflag		!available data, not necessarily read data
        integer ierr            !error code - different from 0 if error

	integer it,idfile,id
        integer ii,l,ie,k,i
        integer nvers,nversaux,nrec
        integer nknaux,nelaux,nlvaux
	integer ibarcl,iconz,iwvert,ieco,imerc,iturb,ibfm
	integer date,time
	real, allocatable :: hlvaux(:)

	logical, parameter :: be = .true.
	logical, parameter :: bn = .false.

	logical rst_want_restart

        read(iunit,end=97) idfile,nvers,nrec
        if( nvers .lt. 3 ) goto 98

        ierr = 0

	nvers_rst = nvers

	date = 0
	time = 0
	iflag = 0

        if( nvers >= 9 ) read(iunit) date,time
        if( nvers >= 11 ) read(iunit) atime

	if( nvers <= 10 ) then
	  it = idfile
	  atime = 0.
	  if( date > 0 ) call dts_to_abs_time(date,time,atime)
	  atime = atime + it
	else if( idfile /= idfrst ) then
	  goto 7
	end if

          read(iunit) nknaux,nelaux,nlvaux
          if( nknaux .ne. nkn_global ) goto 99
          if( nelaux .ne. nel_global ) goto 99
          if( nlvaux .ne. nlv_global ) goto 99

	  call rst_read_vertical(iunit,nvers,nkn,nel,nlv)

	  id = id_hydro_rst
	  call rst_add_flag(id,iflag)
	  if( rst_want_restart(id) ) then
	    call restart_read_value(iunit,be,iwegv)
	    call restart_read_value(iunit,bn,znv)
	    call restart_read_value(iunit,be,3,zenv)
	    call restart_read_value(iunit,be,utlnv)
	    call restart_read_value(iunit,be,vtlnv)
	  else
            read(iunit)
            read(iunit)
            read(iunit)
            read(iunit)
            read(iunit)
	  end if

          if( nvers .ge. 4 ) then
	    id = id_depth_rst
	    call rst_add_flag(id,iflag)
	    if( rst_want_restart(id) ) then
	      call restart_read_value(iunit,be,3,hm3v)
	    else
              read(iunit)
	    end if
          end if

          if( nvers .ge. 5 ) then
	    id = id_barcl_rst
            read(iunit) ibarcl
	    ibarcl_rst = ibarcl
            if( ibarcl .gt. 0 ) then
	      call rst_add_flag(id,iflag)
	      if( rst_want_restart(id) ) then
	        call restart_read_value(iunit,bn,saltv)
	        call restart_read_value(iunit,bn,tempv)
	        call restart_read_value(iunit,bn,rhov)
	      else
                read(iunit)
                read(iunit)
                read(iunit)
	      end if
            end if
          end if

	  if( nvers .ge. 15 ) then
	    id = id_gotm_rst
	    read(iunit) iturb
	    iturb_rst = iturb
	    if( iturb .eq. 1 ) then
	      call rst_add_flag(id,iflag)
	      if( rst_want_restart(id) ) then
	        call read_restart_gotm(iunit)
	      else
	        call skip_restart_gotm(iunit)
	      end if
	    end if
	  end if

          if( nvers .ge. 6 ) then
	    id = id_conz_rst
            read(iunit) iconz
	    iconz_rst = iconz
	    if( iconz > 0 ) then
	      call rst_add_flag(id,iflag)
	      if( rst_want_restart(id) ) then
	        call read_restart_conz(iunit,iconz)
	      else
	        call skip_restart_conz(iunit,iconz)
	      end if
	    end if
	  end if

	  if( nvers .ge. 7 ) then
	    id = id_wvert_rst
	    read(iunit) iwvert
	    iwvert_rst = iwvert
	    if( iwvert .gt. 0 ) then
	      call rst_add_flag(id,iflag)
	      if( rst_want_restart(id) ) then
	        call restart_read_value(iunit,bn,wlnv)
	      else
                read(iunit)
	      end if
	    end if
	  end if

	  if( nvers .ge. 8 ) then
	    id = id_eco_rst
	    read(iunit) ieco
	    ieco_rst = ieco
            if( ieco .gt. 0 ) then
	      call rst_add_flag(id,iflag)
	      if( rst_want_restart(id) ) then
	        call read_restart_eco(iunit)
	      else
	        call skip_restart_eco(iunit)
	      end if
	    end if
          end if

	  if( nvers .ge. 12 ) then
	    id = id_merc_rst
	    read(iunit) imerc
	    imerc_rst = imerc
            if( imerc .gt. 0 ) then
	      call rst_add_flag(id,iflag)
	      if( rst_want_restart(id) ) then
	        call read_restart_mercury(iunit)
	      else
	        call skip_restart_mercury(iunit)
	      end if
	    end if
          end if


	  if( nvers .ge. 17 ) then
	    id = id_bfm_rst
	    read(iunit) ibfm
            if( ibfm .gt. 0 ) then
	      call rst_add_flag(id,iflag)
	      if( rst_want_restart(id) ) then
                call bfm_init_for_restart()
	        call read_restart_bfm(iunit)
	      else
	        call skip_restart_bfm(iunit)
	      end if
	    end if
          end if

        return
    7	continue
	write(6,*) 'rst_read_record: error in idfrst... ' // 'no restart format'
	ierr = 7
	return
   96   continue
        write(6,*) 'error reading header of restart file...'
        stop 'error stop rst_read_record: cannot read header'
   97   continue
        ierr = -1
        return
   98   continue
        write(6,*) 'error reading restart file...'
        write(6,*) 'nvers: ',nvers
        stop 'error stop rst_read_record: cannot read this version'
   99   continue
        write(6,*) 'error reading restart file...'
        write(6,*) 'nkn,nel,nlv:'
        write(6,*) 'shyfem:  ',nkn,nel,nlv
        write(6,*) 'rstfile: ',nknaux,nelaux,nlvaux
        stop 'error stop rst_read_record: incompatible parameters'
        end

!*******************************************************************

	subroutine rst_read_vertical(iunit,nvers,nkn,nel,nlv)

! nkn,... are local values

	use mod_restart
	use shympi

	implicit none

	integer iunit
	integer nvers
	integer nkn,nel,nlv

	logical, parameter :: be = .true.
	logical, parameter :: bn = .false.

	if( .not. allocated(hlvrst) ) then
	  allocate(hlvrst(nlv_global))
	  allocate(ilhrst(nel))
	  allocate(ilhkrst(nkn))
	  hlvrst = 0.
	  ilhrst = 0
	  ilhkrst = 0
	end if

	if( nvers .ge. 14 .and. nlv == 1 ) then
	  hlvrst = 10000.
	  ilhrst = 1
	  ilhkrst = 1
	else
	  if( nvers .ge. 10 ) then
	    read(iunit) hlvrst
	  end if
	  if( nvers .ge. 13 ) then
	    call restart_read_value(iunit,be,ilhrst)
	    call restart_read_value(iunit,bn,ilhkrst)
	  end if
	end if

	return
	end

!*******************************************************************

	subroutine rst_get_vertical(nkn,nel,nlv,hlv,ilhv,ilhkv)

	use mod_restart

	implicit none

	integer nkn,nel,nlv
	real hlv(nlv)
	integer ilhv(nel)
	integer ilhkv(nkn)

	if( nkn <= 0 .or. nel <= 0 .or. nlv <= 0 ) goto 98

	if( .not. allocated(hlvrst) ) then
	  stop 'error stop rst_get_hlv: hlvrst not allocated'
	end if
	if( nlv /= size(hlvrst) ) goto 99
	if( nkn /= size(ilhkrst) ) goto 99
	if( nel /= size(ilhrst) ) goto 99

	hlv = hlvrst
	ilhv = ilhrst
	ilhkv = ilhkrst

	return
   98	continue
	write(6,*) 'nkn,nel,nlv: ',nkn,nel,nlv
	stop 'error stop rst_get_hlv: error in parameters'
   99	continue
	write(6,*) 'nkn: ',nkn,size(ilhkrst)
	write(6,*) 'nel: ',nel,size(ilhrst)
	write(6,*) 'nlv: ',nlv,size(hlvrst)
	stop 'error stop rst_get_hlv: arrays not compatible'
	end 

!*******************************************************************

	subroutine init_old_vars

! this copies vars just read to old so that they are available

	use mod_hydro_vel
	use mod_hydro
	!use mod_hydro_print
	!use mod_hydro_baro

	implicit none

	zeov = zenv
	zov = znv
	utlov = utlnv
	vtlov = vtlnv
	wlov = wlnv

        !call make_new_depth
        !call copy_depth
        !call make_new_depth

	!call ttov
	!call uvint
	!call uvtopr
	!call uvtop0

        !upro  = uprv
        !vpro  = vprv
        !uov   = unv
        !vov   = vnv
        !ulov  = ulnv
        !vlov  = vlnv

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	function rst_is_rst_file(file)

! finds out if file is a restart file

	use mod_restart

	implicit none

	logical rst_is_rst_file
	character*(*) file

	integer idfile,nvers,nrec

	rst_is_rst_file = .false.

	open(1,file=file,status='old',form='unformatted')
	read(1,end=3,err=3) idfile,nvers,nrec
	close(1)

	if( nvers <= 10 ) then
	  if( nrec /= 1 ) return
	else if( idfile /= idfrst ) then
	  return
	end if

	rst_is_rst_file = .true.

	return
    3	continue
	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

        subroutine write_flags(iflag)

        use mod_restart

        implicit none

        integer iflag

        integer id,ivalue
        logical bread
        logical bit10_is_set
        integer bit10_return_value

        write(6,*) 'Meaning of iflag:'
        write(6,*) '         id        flag   set   description'

        do id=1,nidmax
          bread = bit10_is_set(iflag,id)
          ivalue = bit10_return_value(id)
          write(6,*) id,ivalue,'   ',bread,'  ',descript_rst(id)
        end do

        end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	function rst_is_set(id,iflag)
	implicit none
	logical rst_is_set
	integer id,iflag
	logical bit10_is_set
	rst_is_set = .true.
	if( iflag == -1 ) return
	rst_is_set = bit10_is_set(iflag,id)
	end

	subroutine rst_add_flag(id,iflag)
	implicit none
	integer id,iflag
	integer bit10_return_value
	iflag = iflag + bit10_return_value(id)
	end

	function has_restart(id)
	logical has_restart
	integer id
	logical rst_use_restart
	has_restart = rst_use_restart(id)
	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

