
!--------------------------------------------------------------------------
!
!    Copyright (C) 1997-2012,2014-2020  Georg Umgiesser
!    Copyright (C) 2008,2010,2014  Christian Ferrarin
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

! time management routines
!
! contents :
!
!
! revision log :
!
! 23.09.1997	ggu	boundn deleted -> no access to data structure
! 20.03.1998	ggu	minor changes to priout
! 29.04.1998	ggu	new module for semi-implicit time-step
! 07.05.1998	ggu	check for error on return of nrdvecr
! 19.06.1998	ggu	version number is character
! 22.01.1999	ggu	oxygen section added
! 26.01.1999	ggu	new comp3d added
! 11.08.1999	ggu	new compatibility array hlhv initialized
! 19.11.1999	ggu	new routines for section vol
! 20.01.2000	ggu	common block /dimdim/ eliminated
! 04.02.2000	ggu	no priout, dobefor/after, pritime, endtime
! 15.05.2000	ggu	hm3v substituted
! 26.05.2000	ggu	copright statement adjourned
! 21.11.2001	ggu	routines to handle advective index (aix)
! 27.11.2001	ggu	routine to handle info file (getinfo)
! 11.10.2002	ggu	aix routines deleted
! 07.02.2003	ggu	routine added: changeimp, getaz; deleted getaza
! 10.08.2003	ggu	call adjust_chezy instead sp135r
! 14.08.2003	ggu	femver transfered to subver, not called in nlsh2d
! 20.08.2003	ggu	tsmed substituted by ts_shell
! 01.09.2003	ggu	call wrousa
! 03.09.2004	ggu	call admrst, comp3d renamed to init_3d (not used)
! 03.09.2004	ggu	nlv, hlv initialized in nlsh2d (FIXME)
! 28.09.2004	ggu	read lagrangian section
! 01.12.2004	ggu	new routine set_timestep for variable time step
! 17.01.2005	ggu	get_stab_index to newcon.f, error stop in set_timestep
! 14.03.2005	ggu	syncronize idt with end of simulation (set_timestep)
! 07.11.2005	ggu	handle new section sedtr for sediments
! 23.03.2006	ggu	changed time step to real
! 23.05.2007	ggu	recall variable time step pars at every time step
! 02.10.2007	ggu	bug fix in set_timestep for very small rindex
! 10.04.2008	ccf	output in netcdf format
! 28.04.2008	ggu	in set_timestep new call to advect_stability()
! 03.09.2008	ggu	in nlsh2d different error message
! 20.11.2008	ggu	init_3d deleted, nlv initialized to 0
! 18.11.2009	ggu	new format in pritime (write also time step)
! 22.02.2010	ggu	new call to hydro_stability to compute time step
! 22.02.2010	ccf	new routine for tidal pot. (tideforc), locaus deleted
! 26.02.2010	ggu	in set_timestep compute and write ri with old dt
! 22.03.2010	ggu	some comments for better readability
! 29.04.2010	ggu	new routine set_output_frequency() ... not finished
! 04.05.2010	ggu	shell to compute energy
! 22.02.2011	ggu	in pritime() new write to terminal
! 20.05.2011	ggu	changes in set_timestep(), element removal, idtmin
! 31.05.2011	ggu	changes for BFM
! 01.06.2011	ggu	idtmin introduced
! 12.07.2011	ggu	new routine next_output(), revised set_output_frequency
! 14.07.2011	ggu	new routines for original time step
! 13.09.2011	ggu	better error check, rdtitl() more robust
! 23.01.2012	ggu	new section "proj"
! 24.01.2012	ggu	new routine setup_parallel()
! 10.02.2012	ggu	new routines to initialize and access time common block
! 05.03.2014	ggu	code prepared to repeat time step (irepeat) - not ready
! 05.03.2014	ggu	new routines get_last/first_time()
! 10.04.2014	ccf	new section "wrt" for water renewal time
! 29.10.2014	ggu	do_() routines transfered from newpri.f
! 10.11.2014	ggu	time management routines transfered to this file
! 26.11.2014	ggu	changed VERS_7_0_7
! 05.12.2014	ggu	changed VERS_7_0_8
! 12.12.2014	ggu	changed VERS_7_0_9
! 19.12.2014	ggu	accept date also as string
! 23.12.2014	ggu	fractional time step introduced
! 07.01.2015	ggu	fractional time step without rounding (itsplt=3)
! 26.02.2015	ggu	changed VERS_7_1_5
! 30.04.2015	ggu	changed VERS_7_1_9
! 23.09.2015	ggu	time step is now working with dt as double
! 10.10.2015	ggu	use bsync as global to check for syncronization
! 23.09.2016	ggu	cleaned set_timestep()
! 30.09.2016	ggu	changed VERS_7_5_18
! 20.10.2017	ggu	new get_absolute_act_time(),get_absolute_ref_time()
! 04.11.2017	ggu	changed VERS_7_5_34
! 05.12.2017	ggu	changed VERS_7_5_39
! 22.02.2018	ggu	changed VERS_7_5_42
! 23.02.2018	ggu	most parts converted from int to double
! 29.03.2018	ggu	bug fix for syncronization step
! 03.04.2018	ggu	changed VERS_7_5_43
! 03.04.2018	ggu	changed VERS_7_5_44
! 13.04.2018	ggu	hydro_stability includes explicit gravity wave
! 13.04.2018	ggu	set_timestep now is working with mpi
! 16.04.2018	ggu	write warning if time step is over recommended one
! 09.02.2019	ggu	bug fix for syncronization of last time step
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 15.09.2019	ggu	small changes to account for synchorization time step
! 08.02.2020	ggu	utility routines copied to new file
! 16.02.2020	ggu	itunit eliminated
! 16.03.2020	ggu	write also dtime to terminal
! 31.05.2021	ggu	write stability index to inf file
! 15.07.2021	ggu	do not set it (for climatological runs)
! 16.02.2022	ggu	cosmetic changes
! 20.03.2022	ggu	eliminated convert_time -> convert_time_d
! 21.03.2022	ggu	bug in set_timestep: dtmin is real, now use ddtmin as dp
! 28.03.2022	ggu	bug fix: ddtmin was not saved
! 02.04.2023    ggu     only master writes to iuinfo
! 12.12.2023    ggu     introduced dtmax (maximum time to run to)
!
!**********************************************************************
!**********************************************************************
!**********************************************************************

	subroutine print_time

! prints time after time step

	use shympi
	use femtime

	implicit none

        integer nit1,nit2,naver
	integer idtfrac,i
        real perc,dt

	integer year,month,day,hour,min,sec
	integer, save :: isplit,itime
	double precision daux,dtime,ddt,atime

	character*20 dline
	character*9 frac
	character*4 atext
	double precision dgetpar
	logical dts_has_date

	integer, parameter :: MyLongIntType = selected_int_kind (12)
	integer (kind=MyLongIntType) :: inttime

	integer, save :: icall = 0
	integer, save :: iuinfo = 0
	real cpu_time_new
	real, save :: cpu_time_old = 0.
	real, parameter :: cpu_time_max = 1.	!max seconds before flush

!---------------------------------------------------------------
! set parameters and compute percentage of simulation
!---------------------------------------------------------------

	if( icall .eq. 0 ) then
          isplit = nint(dgetpar('itsplt'))
	  call cpu_time(cpu_time_old)
          call getinfo(iuinfo)  !unit number of info file
	end if

        naver = 20
        naver = 0

	ddt = dt_act
	dtime = t_act

        !perc = (100.*(it-itanf))/(itend-itanf)
        perc = (100.*(dtime-dtanf))/(dtend-dtanf)

!---------------------------------------------------------------
! compute total number of iterations
!---------------------------------------------------------------

	nit1 = 0
	idtfrac = 0

	if( bsync ) then	!syncronization - do not count
	  !
	else if( idt .gt. 0 ) then
          nit1 = niter + (dtend-dtime)/ddt
	else if( ddt > 0 ) then	!sub second time step
          daux = (dtend-dtime)/ddt
	  if( daux > 1000000000. ) then
	    write(6,*) '******************************************'
	    write(6,*) '******************************************'
	    write(6,*) '******************************************'
	    write(6,*) t_act,dt_act,dtanf,dtend
	    write(6,*) daux
	    write(6,*) niter,bsync
	    write(6,*) '******************************************'
	    write(6,*) '******************************************'
	    write(6,*) '******************************************'
	    stop 'error stop print_time: internal error'
	  else
            nit1 = niter + nint((dtend-dtime)/ddt)
	    idtfrac = nint(1./ddt)
	  end if
	else
	  write(6,*) 'idt,dt_act: ',idt,ddt
	  write(6,*) 'warning: time step was 0'
	  stop 'error stop print_time: 0 time step'
	end if

	nit2 = nit1
	if( dtime .gt. dtanf ) then
          nit2 = nint(niter*( 1 + (dtend-dtime)/(dtime-dtanf)))
	end if

        nits = nit2
        if( naver .gt. 0 ) nits = ( 1*nit1 + (naver-1)*nit2 ) / naver

	icall = icall + 1

!---------------------------------------------------------------
! write to terminal
!---------------------------------------------------------------

	if( .not. shympi_is_master() ) return

	if( dts_has_date() ) then
	  atime = atime0 + dtime
	  atext = 'date'
	  call dts_format_abs_time(atime,dline)
	else
	  atime = dtime
	  atext = 'time'
	  write(dline,'(f20.4)') dtime
	end if

	if( mod(icall,50) .eq. 0 ) write(6,1004) atext

	  if( isplit == 3 .or. idtorig == 0 ) then
            write(6,1009) dline,dtime,ddt,niter,nits,perc
	  else if( idtfrac == 0 ) then
	    inttime = nint(dtime)
            write(6,1008) dline,dtime,idt,niter,nits,perc
	  else
	    frac = ' '
	    write(frac,'(i9)') idtfrac
	    do i=9,1,-1
	      if( frac(i:i) == ' ' ) exit
	    end do
	    frac(i-1:i) = '1/'
            write(6,1006) dline,frac,niter,nits,perc
	  end if

!---------------------------------------------------------------
! flush if maximum time has passed
!---------------------------------------------------------------

	call cpu_time(cpu_time_new)
	if( cpu_time_new - cpu_time_old > cpu_time_max ) then
	  cpu_time_old = cpu_time_new
	  flush(6)
	  if( iuinfo > 0 ) flush(iuinfo)
	end if

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	return
 1004   format(17x,a4,12x,'dtime',4x,'dt',12x,'iterations',5x,'percent')
 1008   format(1x,a20,1x,f16.2, i6,i10,' /',i10,f10.3,' %')
 1006   format(1x,a20,1x,       a9,i10,' /',i10,f10.3,' %')
 1007   format(1x,a20,1x, f9.2,    i10,' /',i10,f10.3,' %')
 1009   format(1x,a20,1x,f16.2, f9.2,    i10,' /',i10,f8.3,' %')
	end

!********************************************************************

	subroutine print_end_time	!FIXME

! prints stats after last time step

	use shympi
	use femtime

	implicit none

	if( .not. shympi_is_master() ) return

	write(6,*) 'program stop at time = ',aline_act
	write(6,*) 'total iterations = ',niter

	end

!********************************************************************
!********************************************************************
!********************************************************************

	subroutine setup_time

! setup and check time parameters

	use femtime

	implicit none

	integer date,time
	double precision didtm,atime,didt
	character*20 dline

	double precision dgetpar

	call dts_get_date(date,time)
	call dts_to_abs_time(date,time,atime0)

	call convert_date_d('itanf',dtanf)
	call convert_date_d('itend',dtend)
	call convert_time_d('idt',didt)

	if( didt .le. 0 ) then
	  write(6,*) 'Error in compulsory time parameters'
	  write(6,*) 'Time step is not positive'
	  write(6,*) 'idt :',didt
	  stop 'error stop setup_time: idt'
	else if( dtanf+didt .gt. dtend ) then
	  write(6,*) 'Error in compulsory time parameters'
	  write(6,*) 'itend too small, no time step will be performed'
	  write(6,*) 'itanf,itend,idt :',dtanf,dtend,didt
	  atime = atime0 + dtanf
	  call dts_format_abs_time(atime,dline)
	  write(6,*) 'initial time: ',dline
	  atime = atime0 + dtend
	  call dts_format_abs_time(atime,dline)
	  write(6,*) 'final time:   ',dline
	  stop 'error stop setup_time: dtend'
	end if

	niter = 0
	nits = (dtend-dtanf) / didt

	t_act = dtanf
	dt_act = didt
	dt_orig = didt

        atime = atime0 + t_act
        call dts_format_abs_time(atime,aline_act)

	idt = nint(didt)
	idtorig = idt
	itanf = nint(dtanf)
	!itend = nint(dtend)
	it = itanf

	itunit = nint(dgetpar('itunit'))
	if( itunit /= 1. ) then
	  write(6,*) 'itunit = ',itunit
	  write(6,*) 'this parameter is not supported anymore'
	  stop 'error stop setup_time: itunit'
	end if

	end

!********************************************************************

	subroutine setup_date

! setup and check date parameter

	implicit none

	integer date,time
	double precision ddate,dtime
	integer year,month,day,hour,min,sec
	integer ierr
	character*20 text
	double precision dgetpar

	call getfnm('date',text)
	if( text .ne. ' ' ) then	!string value given
	  call dtsunform(year,month,day,hour,min,sec,text,ierr)
	  if( ierr .ne. 0 ) then
	    write(6,*) 'date: ',text
	    stop 'error stop setup_date: cannot parse date'
	  end if
	  date = 10000*year + 100*month + day
	  time = 10000*hour + 100*min + sec
	  ddate = date
	  dtime = time
	  call dputpar('date',ddate)
	  call dputpar('time',dtime)
	  !write(6,*) '===================================='
	  !write(6,*) 'date as string: ',date,time
	  !write(6,*) '===================================='
	end if

	date = nint(dgetpar('date'))
	time = nint(dgetpar('time'))

	call dtsini(date,time)

	if( date >= 0 ) return

	write(6,*) 'The date parameter has not been set.'
	write(6,*) 'The date parameter indicates the absolute date'
	write(6,*) 'to which FEM time is relative to.'
	write(6,*) 'The format for the date parameter is YYYY[MMDD].'
	write(6,*) 'If you want to do without this faeture, please'
	write(6,*) 'explicitly set the date parameter to 0'
	write(6,*) 'in the parameter file.'
	write(6,*) 'Alternatively, you can also specify date as a string.'
	write(6,*) 'The format is ''YYYY-MM-DD[::hh:mm:ss]'''
	write(6,*) 'Examples:'
	write(6,*) '   date = 20120101'
	write(6,*) '   date = 2012                 # same as above'
	write(6,*) '   date = ''2012-01-01''         # same as above'
	write(6,*) '   date = 0                    # no date feature'

	stop 'error stop ckdate: date'
	end

!********************************************************************

	subroutine get_date_time(date,time)

	implicit none

	integer date,time

	double precision dgetpar

	date = nint(dgetpar('date'))
	time = nint(dgetpar('time'))

	end

!********************************************************************
!********************************************************************
!********************************************************************

        subroutine check_timestep(irepeat)

! still to be implemented

	implicit none

	integer irepeat		!on return 1 if time step has to be repeated

	irepeat = 0

	end

!********************************************************************
!********************************************************************
!********************************************************************

        subroutine set_timestep(dtmax)

! controls time step and adjusts it

	use shympi
	use femtime
	use mod_info_output

        implicit none

        double precision :: dtmax	!maximum time for run (normally dtend)

	logical bdebug
        integer idtdone,idtrest,idts
	integer idtfrac
        integer istot,iss
	integer idta(n_threads)
        double precision dt,dtnext,atime,ddts,dtsync,dtime,dt_recom
        double precision, save :: ddtmin
        double precision :: dtbest
	real dtr,dtaux
        real ri,rindex,rindex1,sindex
	real perc,rmax
	real dhpar,chpar,thpar,shpar
	real array(5)
	character*80 format

        real, save :: cmax,tfact,dtmin,zhpar
        integer, save :: idtsync,isplit,idtmin
        integer, save :: iuinfo = 0
        integer, save :: icall = 0

	double precision dgetpar

	bdebug = .false.

        if( icall .eq. 0 ) then
          istot = 0

          isplit = nint(dgetpar('itsplt'))
          cmax   = dgetpar('coumax')
          tfact  = dgetpar('tfact')	!still to be commented

          dhpar = dgetpar('dhpar')
          chpar = dgetpar('chpar')
          thpar = dgetpar('thpar')
          shpar = dgetpar('shpar')
	  zhpar = max(dhpar,chpar,thpar,shpar)	!max scalar diffusion parameter

	  !call convert_time('idtsyn',idtsync)
	  !call convert_time('idtmin',idtmin)
	  call convert_time_d('idtsyn',dtsync)
	  idtsync = nint(dtsync)
	  call convert_time_d('idtmin',ddtmin)
	  idtmin = nint(ddtmin)

	  if( shympi_is_master() ) then
            call getinfo(iuinfo)  !unit number of info file
	  end if

        end if

        icall = icall + 1

	call info_output('new timestep','none',0,array,.false.,'(a)')

!----------------------------------------------------------------------
!        idtsync = 0             !time step for syncronization
!        cmax = 1.0              !maximal Courant number permitted
!        isplit = -1             !mode for variable time step:
!                                ! -1:  time step fixed, 
!                                !      no computation of rindex
!                                !  0:  time step fixed
!                                !  1:  split time step
!                                !  2:  optimize time step (no multiple)
!                                !  3:  optimize time step (fractional)
!	 idtmin = 1		 !minimum time step allowed
!	 tfact = 0		 !factor of maximum decrease of time step
!----------------------------------------------------------------------

	call check_time('in set_timestep start')

        if( isplit .ge. 0 ) then
          dtr = 1.
          call hydro_stability(dtr,rindex)
	  !call scalar_basic_stability(dtr,zhpar,sindex)
	  !write(6,*) 'stability: ',rindex,sindex
        else
          rindex = 0.
        end if

!----------------------------------------------------------------------
! syncronize stability index between domains if running in mpi mode
!----------------------------------------------------------------------

	!write(6,*) 'time domains: ',my_id,rindex,1./rindex,cmax
	rindex = shympi_max(rindex)
	!write(6,*) 'time final: ',my_id,rindex,1./rindex,cmax

!----------------------------------------------------------------------
! split time step
!----------------------------------------------------------------------

	istot = 0
	idtfrac = 0

        if( isplit .le. 0 ) then
          idts = 0
	  dt = dt_act
        else if( isplit .eq. 1 ) then
	  if( idtorig <= 0 ) then
	    stop 'error stop set_timestep: idtorig==0 and isplit==1'
	  end if
          idts = idtorig
	  dt = dt_act
	  call split_equal(rindex,cmax,dt,istot)
        else if( isplit .eq. 2 .or. isplit .eq. 3 ) then
          idts = idtsync
	  idtfrac = 0
	  dt = dt_orig
	  if( rindex > 0 ) dt = cmax / rindex	! maximum allowed time step
	  if( dt >= dt_orig ) then
	    dt = dt_orig
	  else if( dt >= 1. ) then
	    if( isplit .eq. 2 ) dt = int(dt)
	  else
	    idtfrac = ceiling(1./dt)
	    if( isplit .eq. 2 ) dt = 1. / idtfrac
	  end if
        else
          write(6,*) 'isplit = ',isplit
          stop 'error stop set_timestep: value for isplit not allowed'
        end if

        dtaux = 0.
        if( rindex > 0 ) dtaux = cmax / rindex  ! maximum allowed time step
        !if( iuinfo > 0 ) then
	!  write(iuinfo,*) 'stability_hydro: ',rindex,dtaux,dt
	!end if
	array(1:3) = (/rindex,dtaux,real(dt)/)
	call info_output('stability_hydro','none',3,array,.false.)

	if( dt <= 0 ) then
	  write(6,*) 'dt is negative after setting'
	  write(6,*) dtr,isplit,idts,idtfrac,dt_orig
	  write(6,*) dt,cmax,rindex,cmax/rindex
	  stop 'error stop set_timestep: negative time step'
	end if

	call check_time('in set_timestep 1')

!----------------------------------------------------------------------
! dt	 is proposed new time step
! idts   is time step with which to syncronize
! istot  is number of internal time steps (only for isplit = 1)
! rindex is computed stability index (refers to time step == 1)
!----------------------------------------------------------------------

	if( rindex > 0. ) then
	  dt_recom = 1. / rindex
	  if( dt > dt_recom ) then
	    write(6,*) 'warning: time step is bigger than recommended'
	    write(6,*) 'dt recommended: ',dt_recom
	    write(6,*) 'dt used       : ',dt
	    write(6,*) 'this might lead to instability'
	    write(6,*) 'If you know what you are doing'
	    write(6,*) 'then set itsplt = -1'
	  end if
	end if

!----------------------------------------------------------------------
!	syncronize time step
!----------------------------------------------------------------------

	dtbest = dt
	dtime = t_act
	dtsync = idts
        call sync_step(dtanf,dtsync,dtime,dt,bsync)

	if( dt <= 0 ) stop 'error stop set_timestep: dt<=0'

	if( dtime .gt. dtmax ) then	!sync with end of sim
	  dt = dtmax - t_act
	  dtime = dtmax
	  bsync = .true.
	end if

	!write(6,*) 'set_timestep: sync ',rindex,dt

	if( dt <= 0 ) then
	  write(6,*) 'dt is negative after syncronize'
	  write(6,*) dtbest
	  write(6,*) dtime,t_act
	  write(6,*) dtanf,dtend,dtmax
	  write(6,*) dtr,isplit,idtfrac,dt_orig
	  write(6,*) dt,cmax,rindex,cmax/rindex
	  write(6,*) dtnext,bsync
	  stop 'error stop set_timestep: negative time step'
	end if

	call check_time('in set_timestep 2')

!----------------------------------------------------------------------
! ri     is stability index for computed time step
! dtmin  is minimum time step allowed
! dt	 is the computed time step
!----------------------------------------------------------------------

        ri = dt*rindex

	!write(6,*) bsync
	!write(6,*) dt
	!write(6,*) ddtmin
        if( dt .lt. ddtmin .and. .not. bsync ) then    !should never happen
	  dtr = dt
          call error_stability(dtr,rindex)
          write(6,*) 'dt is less than ddtmin'
          !write(6,*) it,itanf,mod(it-itanf,idtorig)
          !write(6,*) idtnew,idtdone,idtrest,idtorig
          write(6,*) idtorig
          write(6,*) idts,idtsync
          write(6,*) dtime,dt,ddtmin
          write(6,*) isplit,istot
          write(6,*) cmax,rindex,ri
	  write(6,*) 'possible computed time step:  dt = ',dt
	  write(6,*) 'minimum time step allowed: ddtmin = ',ddtmin
	  write(6,*) 'please lower dtmin in parameter input file'
          stop 'error stop set_timestep: time step too small'
        end if

!----------------------------------------------------------------------
! set new values
!----------------------------------------------------------------------

        niter=niter+1

	dt_act = dt
	t_act = dtime
	idt = dt
	!it = t_act

	call check_time('in set_timestep end')
	
	atime = atime0 + t_act
	call dts_format_abs_time(atime,aline_act)

	if( bdebug ) then
	  write(107,*) '========================'
	  write(107,*) t_act,dt_orig,dt_act,bsync
	  write(107,*) dtanf,dtend
	  write(107,*) rindex,ri
	  write(107,*) '========================'
	end if

	!perc = (100.*(it-itanf))/(itend-itanf)
	perc = (100.*(t_act-dtanf))/(dtend-dtanf)
	iss = 0
	if( bsync ) iss = 1

!	if( iuinfo > 0 ) then
!          write(iuinfo,1004) 'timestep: ',aline_act &
!     &				,t_act,istot,iss,dt,perc
!	  flush(iuinfo)	
!	end if
	array = (/real(t_act),real(istot),real(iss),real(dt),real(perc)/)
	format = '(a,f18.4,f6.0,f4.0,2f10.2)'
	call info_output(' timestep','none',5,array,.true.,format)

	call shympi_barrier

	if( bsync .or. mod(icall,50) == 0 ) then
	  if( iuinfo > 0 ) flush(iuinfo)
	  flush(6)
	end if

        return
 1004   format(a,a20,f18.4,i5,i3,2f10.2)
        end

!**********************************************************************

        subroutine sync_step(dtanf,dtsync,dtime,dt,bsync)

! syncronizes time step and returns dtime,dt,bsync

        implicit none

        double precision dtanf,dtsync
        double precision dtime,dt
        logical bsync

        double precision dtnext,dtold

        bsync = .false.
        dtold = dtime
        dtime = dtime + dt
        if( dtsync <= 0. ) return

        dtnext = dtanf + dtsync*(1.+aint((dtold-dtanf)/dtsync))

        if( dtime >= dtnext ) then
          dt = dtnext - dtold
          dtime = dtnext
          bsync = .true.
        end if

        end

!**********************************************************************

	subroutine split_equal(rindex,cmax,dt,istot)

! split time step (macro timestep) into equal parts
!
! dt is the new time step, all the rest can be computed from dt

	use femtime

	implicit none

	real rindex,cmax
	double precision dt
	integer istot

	integer idtnew
	real ri

        if( mod(it-itanf,idtorig) .ne. 0 ) return	!inside macro time step

        istot = idtorig*rindex/cmax	!starting point for istot
        if( istot < 1 ) istot = 1

        do 
          idtnew = idtorig/istot
	  if( idtnew == 0 ) exit	!istot too big
          ri = idtnew*rindex
          if( ri .le. cmax .and. istot*idtnew == idtorig ) exit
          istot = istot + 1
        end do

	if( idtnew == 0 ) then
	  write(6,*) rindex,cmax,istot,idtnew
	  stop 'error stop split_equal: cannot find new time step'
	end if

	dt = idtnew

	end

!**********************************************************************
!**********************************************************************
!**********************************************************************

	subroutine check_time(text)

! to be deleted later

	use femtime

	implicit none

	character*(*) text

	integer, save :: nbcheck = 0
	double precision, save :: t_old = 0
	integer ifemop

! -1372895700.0000
! -3520376896.0000

	return

	if( nbcheck == 0 ) then
	  nbcheck=ifemop('.check.txt','form','new')
	  if( nbcheck == 0 ) stop 'error stop check_time'
	end if

	if( t_act < -1372895700.0000 ) then
	  write(6,*) 'error from check_time: ',trim(text)
	  write(6,*) t_act,t_old,dt_act,dtanf,dtend
	  write(nbcheck,*) 'error from check_time: ',trim(text)
	  write(nbcheck,*) t_act,t_old,dt_act,dtanf,dtend
	end if
	
	if( t_act > t_old ) t_old = t_act

	end

!********************************************************************
!********************************************************************
!********************************************************************

        subroutine test_sync

        implicit none

        logical bsync
        integer ierr,nrec
        double precision dtanf,dtsync,dtorig,dt,dtime,r,dtrun,dts
        double precision dtmax,dtnext,dnew,dtn,dyears

        dtanf = 0
        dtanf = -10000000.
        dtmax = dtanf
        dtmax = 1000000000.
        dtsync = 3600
        dtorig = 300

        dt = dtorig
        dtnext = dtanf + dtsync
        ierr = 0
        nrec = 0

        dtime = dtanf
        do
          nrec = nrec + 1
          call random_number(r)
          dt = dtorig/2. + 1.5*(r-0.5)*dtorig
          if( dt > dtorig ) dt = dtorig
          if( dt < 1. ) dt = 1.
          !dt = int(dt)

          call sync_step(dtanf,dtsync,dtime,dt,bsync)

          !write(6,*) dtime,dt
          if( bsync ) then
            write(6,*) 'syn',dtime,dt
            if( dtime /= dtnext ) then
              write(6,*) '**** ',dtnext,dtime
              ierr = ierr + 1
              stop 'error stop: syncronize time wrong'
            end if
            dtnext = dtnext + dtsync
          else
            if( dtime > dtnext ) then
              write(6,*) dtime,dtnext
              stop 'error stop: no syncronization'
            end if
          end if
          if( dtmax > dtanf .and. dtime > dtmax ) exit
        end do

        if( ierr > 0 ) write(6,*) 'errors: ',ierr
        if( nrec > 0 ) write(6,*) 'time steps: ',nrec
        dyears = (dtime-dtanf)/(365*86400)
        write(6,*) 'years simulated: ',dyears

        end

!******************************************************************
!        program test_main
!        call test_sync
!        end
!******************************************************************

