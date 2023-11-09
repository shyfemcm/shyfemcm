
!--------------------------------------------------------------------------
!
!    Copyright (C) 1993,1997-1999,2002-2004,2008,2010  Georg Umgiesser
!    Copyright (C) 2014-2019  Georg Umgiesser
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

! residual currents
!
! contents :
!
! subroutine resid		computes residual currents
! subroutine rmsvel		computes rms currents
! subroutine dischar		computes discharge (old)
! function resi(zov,znv,n)	computes residuum
! subroutine tsmed 		computes average of scalar values
!
! revision log :
!
! 21.12.1993	ggu	written residual currents
! 15.05.1997	ggu	$$IDTRMS - new parameters for RMS velocity
! 21.05.1997	ggu	$$AVERZR - average zr for dry/wet algorithm
! 21.05.1997	ggu	$$HREFBUG - bug in getting href, hzoff
! 16.06.1997	ggu	$$RESCU - not used subroutine rescu deleted
! 28.04.1998	ggu	dischar deleted
! 20.05.1998	ggu	use ifileo to open file
! 26.01.1999	ggu	use 3D NOS routines for rms velocity
! 18.10.1999	ggu	function resi copied from newres
! 08.10.2002	ggu	subroutine tsmed to compute average of T/S
! 10.10.2002	ggu	subroutine tsmed to compute min/max of T/S
! 19.08.2003	ggu	some cleanup in subroutine rmsvel
! 19.08.2003	ggu	new routines cmed_init, cmed_accum, ts_shell
! 04.03.2004	ggu	bug fix in cmed_accum() for array with more variables
! 10.08.2004	ggu	cmed_init, cmed_accum adjusted also for 2D
! 25.11.2004	ggu	resid converted to 3D (not tested)
! 09.10.2008	ggu	new call to confop
! 23.03.2010	ggu	changed v6.1.1
! 20.01.2014	ggu	new calls to ous routines
! 28.01.2014	ggu	changed VERS_6_1_71
! 26.11.2014	ggu	changed VERS_7_0_7
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 05.06.2015	ggu	changed VERS_7_1_12
! 10.07.2015	ggu	changed VERS_7_1_50
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 18.09.2015	ggu	changed VERS_7_2_3
! 23.09.2015	ggu	changed VERS_7_2_4
! 28.04.2016	ggu	changed VERS_7_5_9
! 17.05.2016	ggu	in ts_shell() save isvect/itvect
! 25.05.2016	ggu	changed VERS_7_5_10
! 31.10.2016	ggu	new output format in rms files
! 12.01.2017	ggu	changed VERS_7_5_21
! 03.04.2018	ggu	changed VERS_7_5_43
! 11.05.2018	ggu	changed VERS_7_5_47
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 15.07.2021	ggu	do not use -> it must be converted to dtime
! 21.03.2022	ggu	new version with da_out
! 22.03.2022	ggu	cmed routines removed and put into subconutil.f
!
!********************************************************************
!
	subroutine resid
!
! computes residual currents
!
	use mod_depth
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use simul

	implicit none

! common
! local
	character*80 nam,dir,file
	integer ierr,ii,ie,k,id
        integer nvers,lmax,l
	integer date,time
	integer ftype,nvar
	real href,hzoff,rr,hm
	double precision dtime
	character*80 title,femver
	character*20 aline
! function
	integer iround
	integer wfout,wrout,ifileo
	real getpar
	double precision dgetpar
	integer ifemop
	logical has_output_d,is_over_output_d,next_output_d
! save
	real, save, allocatable :: ur(:,:)
	real, save, allocatable :: vr(:,:)
	real, save, allocatable :: znr(:)
	real, save, allocatable :: zer(:,:)

	integer, save :: icall = 0
	integer, save :: nr = 0
	logical, save :: bfirst = .true.
	logical, save :: bdebug = .false.
	double precision, save :: da_out(4)

	if(icall.eq.-1) return

!---------------------------------------------------------------------
! initialize routine and variables
!---------------------------------------------------------------------

	if(icall.eq.0) then
          call init_output_d('itmres','idtres',da_out)
	  call increase_output_d(da_out)	!itres=itmres+idtres
          if( .not. has_output_d(da_out) ) icall = -1
	  if(icall.eq.-1) return

	  call shyfem_init_hydro_file('resid',.false.,id)
          da_out(4) = id

	  bdebug = iround(getpar('levdbg')) .ge. 1

	  if( bdebug ) then
	    call get_act_timeline(aline)
	    write(6,*) 'rmsvel : res-file opened ',aline
	  end if

	  allocate(ur(nlvdi,nel))
	  allocate(vr(nlvdi,nel))
	  allocate(znr(nkn))
	  allocate(zer(3,nel))
	  nr=0
	  ur = 0.
	  vr = 0.
	  znr = 0.
	  zer = 0.
	end if

!---------------------------------------------------------------------
! already ready for adding?
!---------------------------------------------------------------------

	icall=icall+1

	if( .not. is_over_output_d(da_out) ) return	!before start of accum

	if( bfirst ) then
	   call get_act_timeline(aline)
	   write(6,*) 'resid: starting summing ',aline
	   bfirst = .false.
	end if

!---------------------------------------------------------------------
! sum transports
!---------------------------------------------------------------------

	nr=nr+1
	do ie=1,nel
          lmax = ilhv(ie)
          do l=1,lmax
	    ur(l,ie)=ur(l,ie)+utlnv(l,ie)
	    vr(l,ie)=vr(l,ie)+vtlnv(l,ie)
          end do
	  do ii=1,3
	    zer(ii,ie)=zer(ii,ie)+zenv(ii,ie)
	  end do
	end do
        do k=1,nkn
	  znr(k)=znr(k)+znv(k)
	end do

!---------------------------------------------------------------------
! is it time to write file ?
!---------------------------------------------------------------------

	if( .not. next_output_d(da_out) ) return

!---------------------------------------------------------------------
! write results into file
!---------------------------------------------------------------------

	if( bfirst ) then
	   call get_act_timeline(aline)
	   write(6,*) 'resid: res-file written ',aline,nr
	end if

	  rr=1./nr

	  do ie=1,nel
            lmax = ilhv(ie)
            do l=1,lmax
	      ur(l,ie)=ur(l,ie)*rr
	      vr(l,ie)=vr(l,ie)*rr
            end do
	    do ii=1,3
	      zer(ii,ie)=zer(ii,ie)*rr
	    end do
	  end do
          do k=1,nkn
	    znr(k)=znr(k)*rr
	  end do

          id = nint(da_out(4))
          call get_act_dtime(dtime)
          call shy_write_hydro_records(id,dtime,nlvdi,znr,zer &
     &                                  ,ur,vr)

!---------------------------------------------------------------------
! reset variables
!---------------------------------------------------------------------

	  nr=0
	  ur = 0.
	  vr = 0.
	  znr = 0.
	  zer = 0.

!---------------------------------------------------------------------
! end of routine
!---------------------------------------------------------------------

	return
   97   continue
	write(6,*) 'Error opening residual file : '
	write(6,*) file
	stop 'error stop : resid'
   78   continue
	write(6,*) 'Error writing record 1 of residual file: ',ierr
	stop 'error stop : resid'
   75   continue
	write(6,*) 'Error writing record 2 of residual file: ',ierr
	stop 'error stop : resid'
   79   continue
	write(6,*) 'Error writing data record of residual file: ',ierr
	stop 'error stop : resid'
	end
!
!********************************************************************
!
	subroutine rmsvel
!
! computes rms currents
!
	use mod_hydro_baro
	use mod_hydro
	use levels
	use basin
	use mod_layer_thickness
	use simul

	implicit none

! common
! local
	double precision rr,dtime
	logical bout
	integer ii,ie,k,id,idc
	integer l,lmax
	real hm,u,v
	real v1v(nkn),v2v(nkn)
	real rmse(nlvdi,nel)
	real rmsn(nlvdi,nkn)
	real aux(nlvdi,nkn)
	character*20 aline
! function
	real getpar
	logical has_output_d,is_over_output_d,next_output_d
! save
	double precision, save, allocatable :: rms(:,:)
	logical, save :: bdebug
	double precision, save :: da_out(4)
	integer, save :: icall = 0
	integer, save :: nr = 0

	if(icall.eq.-1) return

!-----------------------------------------------------------
! initialization
!-----------------------------------------------------------

	if(icall.eq.0) then
          call init_output_d('itmrms','idtrms',da_out)
	  call increase_output_d(da_out)	!itres=itmres+idtres

	  if( .not. has_output_d(da_out) ) icall = -1
	  if(icall.eq.-1) return

	  call shyfem_init_scalar_file('rms',1,.false.,id)
	  da_out(4) = id

	  bdebug = nint(getpar('levdbg')) .ge. 1
	  if( bdebug ) then
	    call get_act_timeline(aline)
	    write(6,*) 'rmsvel : rms-file opened ',aline
	  end if

	  allocate(rms(nlvdi,nel))
	  nr=0
	  rms = 0.
	end if

!-----------------------------------------------------------
! see if we have to start summing
!-----------------------------------------------------------

	icall=icall+1

	bout = is_over_output_d(da_out)
	if( .not. bout ) return

!-----------------------------------------------------------
! sum transports
!-----------------------------------------------------------

	nr=nr+1
	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    hm=hdenv(l,ie)
	    u=utlnv(l,ie)/hm
	    v=vtlnv(l,ie)/hm
	    rms(l,ie) = rms(l,ie) + u*u + v*v
	  end do
	end do

!-----------------------------------------------------------
! if time write output to rms file
!-----------------------------------------------------------

	bout = next_output_d(da_out)
	if( bout ) then
	  if( bdebug ) then
	    call get_act_timeline(aline)
	    write(6,*) 'rmsvel : rms-file written ',aline,nr
	  end if

	  rr=1./nr
	  rmse = sqrt(rms*rr)

	  call e2n3d(nlvdi,rmse,rmsn,aux)

	  call get_act_dtime(dtime)
	  idc = 18
          id = nint(da_out(4))
	  call shy_write_scalar_record(id,dtime,idc,nlvdi,rmsn)
	  call shy_sync(id)

	  nr=0
	  rms = 0.
	end if

!-----------------------------------------------------------
! end of routine
!-----------------------------------------------------------

	end

!********************************************************************

        function resi(zov,znv,n)

! computes residuum

        implicit none

	real resi
        integer n
        real zov(n),znv(n)

        integer i
        real res,var,epsr
        data epsr /1.e-5/

        res=0.

        do i=1,n
           var=zov(i)
           if(abs(var).gt.epsr) then
                res=res+abs((var-znv(i))/var)
           else if(abs(znv(i)).gt.epsr) then
                res=res+abs((var-znv(i))/znv(i))
           end if
        end do

        resi=res

        end

!********************************************************************

