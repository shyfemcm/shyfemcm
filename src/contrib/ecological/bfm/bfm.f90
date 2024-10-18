
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016,2018-2020  Georg Umgiesser
!    Copyright (C) 2016,2019  Erik Pascolo
!    Copyright (C) 2016  Leslie Aveytua
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

! routines for BFM module
!
! revision log :
!
! 22.02.2016	ggu&erp	new bfm routines created from newconz
! 22.03.2016	ggu	changed VERS_7_5_6
! 06.06.2016	ggu	initialization from file changed
! 09.06.2016 	laa	modified (Search LAA)
! 10.06.2016	ggu	changed VERS_7_5_13
! 17.06.2016	ggu	changed VERS_7_5_15
! 28.06.2016	ggu	initialize bfmv, new routine bfm_check()
! 09.09.2016	ggu	changed VERS_7_5_17
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 01.05.2019	erp	restructured for new BFM
! 17.10.2019	ggu	restructured and adapted
! 09.04.2020	ggu	new routine bfm_run()

!==================================================================
        module mod_bfm
!==================================================================

      use mod_sinking
      use mod_depth
      use mod_ts
      use mod_diff_visc_fric
      use mod_hydro_print
      use mod_hydro
      use levels
      use basin
      use mod_meteo
      use mod_layer_thickness

        implicit none

#include "BFM_var_list.h"      

        integer,  save :: nst_bfm = 0
        integer,  save :: nkn_bfm = 0
        integer,  save :: nlv_bfm = 0
        integer,  save :: nst_bfm_d = 42 ! number of variables among the 114 variables of d to be saved in bfm_d

        real, allocatable, save :: bfmv(:,:,:)
        real, allocatable, save :: bfm_d(:,:,:)  ! here we save nst_bfm_d variables among the 114 of d
        real, allocatable, save :: ave_bfmv(:,:,:)
        real, allocatable, save :: ave_bfm_d(:,:,:)  ! here we save nst_bfm_d variables among the 114 of d
        real, allocatable, save :: ave_bfm_benthic(:,:,:)  ! here we save nst_bfm_d variables among the 114 of d
        real, allocatable, save :: wsink_bfm(:,:,:)  ! here we save nst_bfm_d variables among the 114 of d
        real, allocatable, save :: rload_bfm(:,:,:)  ! here we save nst_bfm_d variables among the 114 of d
!       Benthic variables
        real, allocatable, save :: bfm_benthic(:,:,:)

	integer, save :: ibfm = 0
	integer, save :: ibfm_state = jptra	!number of state variables
        integer, save :: icall_bfm = 0
        integer, save :: ninfo_bfm = 0
	integer, save :: ibfm_state_benthic = 9!number of benthic state variables

        logical, private, save :: bdebug = .false.
        logical, save :: binfo_bfm = .true.

        real, save :: cref,rkpar,difmol

        real, save :: counter_ave = 0.

        integer, save, allocatable :: idbfm(:)

        double precision, save :: da_out(4)
        double precision, save :: id_out(5)

        real, save, allocatable :: bfminit(:)
        real, save, allocatable :: bfmbound(:)

	character*3, save :: what = 'bfm'
        logical,save :: mod_bfm_initialized=.False. 
        integer, save :: ieco_av=0

!==================================================================
	contains
!==================================================================

        subroutine mod_bfm_init(nst,nkn,nlv)

	implicit none

        integer nst
        integer nkn
        integer nlv
        real getpar

        if( nst == nst_bfm .and. nkn == nkn_bfm .and. nlv == nlv_bfm ) return

        if( nst > 0 .or. nkn > 0 .or. nlv > 0 ) then
          if( nst == 0 .or. nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nst,nkn,nlv: ',nst,nkn,nlv
            stop 'error stop mod_bfm_init: incompatible parameters'
          end if
        end if

        if( nkn_bfm > 0 ) then
          deallocate(bfmv)
          deallocate(bfm_d)
        end if

        nst_bfm = nst
        nkn_bfm = nkn
        nlv_bfm = nlv

        if( nkn == 0 ) return

	write(6,*) 'bfm mod init: ',nst,nkn,nlv
        allocate(bfmv(nlv,nkn,nst))
	bfmv = 0.
        allocate(bfm_d(nlv,nkn,nst_bfm_d))
	bfm_d = 0.
        ieco_av=nint(getpar('ieco_av'))
        if(ieco_av==1)then
          allocate(ave_bfmv(nlv,nkn,nst))
          allocate(ave_bfm_d(nlv,nkn,nst_bfm_d))
        elseif(ieco_av==2) then
          allocate(ave_bfmv(1,nkn,nst))
          allocate(ave_bfm_d(1,nkn,nst_bfm_d))
        endif
        allocate(ave_bfm_benthic(1,nkn,ibfm_state_benthic))
        if(ieco_av>0)then
          ave_bfmv(:,:,:)=0.0
          ave_bfm_d(:,:,:)=0.0
          ave_bfm_benthic(:,:,:)=0.0
        endif
        allocate(wsink_bfm(nlv,nkn,nst))
        wsink_bfm(:,:,:)=0.0
        allocate(rload_bfm(nlv,nkn,nst))
        rload_bfm(:,:,:)=0.0
        
!       Benthic module
        allocate(bfm_benthic(1,nkn,ibfm_state_benthic))
        bfm_benthic(1,:,1)=3.0    ! PO4
        bfm_benthic(1,:,2)=30.0   ! NO3
        bfm_benthic(1,:,3)=200.0  ! NH4
        bfm_benthic(1,:,4)=8.0    ! SiO3
        bfm_benthic(1,:,5)=48000. ! DIC
        bfm_benthic(1,:,6)=.0     ! R6c
        bfm_benthic(1,:,7)=.0     ! R6n
        bfm_benthic(1,:,8)=.0     ! R6p
        bfm_benthic(1,:,9)=.0     ! R6s
        end subroutine mod_bfm_init

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine bfm_run
      
	implicit none
      
	if( ibfm < 0 ) return

	if( bdebug ) write(6,*) 'starting bfm_run'

	if( bdebug ) write(6,*) 'running bfm_compute'
        call bfm_compute

	if( bdebug ) write(6,*) 'running bfm_reactor'
        call bfm_reactor

	if( bdebug ) write(6,*) 'running bfm_write_output_file'
        call bfm_write_output_file

	if( bdebug ) write(6,*) 'finished bfm_run'

	end subroutine bfm_run

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine bfm_init

! initializes bfm computation

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
	
	implicit none

	integer nvar,nbc,nintp,i,id
	integer nmin
	double precision dtime0

	logical has_restart,has_output,has_output_d
	integer nbnds
	real getpar

!-------------------------------------------------------------
! initialization of module
!-------------------------------------------------------------

        if( mod_bfm_initialized==.False. ) then 
	if( ibfm < 0 ) return

        if( ibfm == 0 ) then
          ibfm=nint(getpar('ibfm'))
          if( ibfm <= 0 ) ibfm = -1
          if( ibfm < 0 ) return

          call mod_bfm_init(ibfm_state,nkn,nlvdi)
	  call bfm_init_internal

          write(6,*) 'bfm initialized in bfm_init: ', ibfm,ibfm_state,nkn,nlvdi
          mod_bfm_initialized = .True.
        end if
        end if

!-------------------------------------------------------------
! initialization of parameters
!-------------------------------------------------------------

	rkpar=getpar('chpar')
	difmol=getpar('difmol')

	call get_first_dtime(dtime0)

	nvar = ibfm_state
	allocate(bfminit(nvar))
	allocate(bfmbound(nvar))
	bfminit = 0.				!default initial condition
	bfmbound = 3.				!default boundary condition

	!call bfm_check('before init')

        if( .not. has_restart(4) ) then	!no restart of conzentrations
	  call bfm_init_file(dtime0,nvar,nlvdi,nlv,nkn,bfminit,bfmv)
	end if

	!call bfm_check('after init')

	call bfm_init_output_file

        call getinfo(ninfo_bfm)

        nbc = nbnds()
        allocate(idbfm(nbc))
        idbfm = 0

	nintp = 2
        call bnds_init_new('bfm',dtime0,nintp,nvar,nkn,nlv,bfmbound,idbfm)

	end subroutine bfm_init

!*********************************************************************

	subroutine bfm_compute

	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	logical bfirst
	integer nvar,i
	real wsink
        real one
	double precision dtime
	
!-------------------------------------------------------------
! initialization
!-------------------------------------------------------------

	if( ibfm < 0 ) return

	call is_time_first(bfirst)
	if( bfirst ) stop 'bfm_compute_multi: internal error'

!-------------------------------------------------------------
! normal call
!-------------------------------------------------------------

	nvar = ibfm_state
	call get_act_dtime(dtime)

	!call bfm_check('before advect')

	call bnds_read_new(what,idbfm,dtime)
        one=1.0
	do i=1,nvar

!$OMP TASK FIRSTPRIVATE(i,rkpar,difhv,difv,difmol,idbfm,what,nlvdi) SHARED(bfmv,one,wsink_bfm,rload_bfm)  DEFAULT(NONE)
 
! use scal_adv_fact to pass wsink matrix into advection

          call scal_adv_fact(what,i,one,bfmv(:,:,i), &
                             idbfm,rkpar,            &
                             one,wsink_bfm(:,:,i),   &
                             one,rload_bfm(:,:,i),   &
                             difhv,difv,difmol)

!$OMP END TASK

	end do	

!$OMP TASKWAIT

	!call bfm_check('after advect')

	icall_bfm = icall_bfm + 1

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end subroutine bfm_compute

!*********************************************************************

	subroutine bfm_init_output_file

	implicit none

	integer nvar,id
	logical has_output_d
       

	nvar = ibfm_state

        call init_output_d('itmcon','idtcon',da_out)

        if( has_output_d(da_out) ) then
          call shyfem_init_scalar_file('bfm',nvar,.false.,id)
          id_out(1) = id
          call shyfem_init_scalar_file('bfm_d',nst_bfm_d,.false.,id)
          id_out(2) = id
          call shyfem_init_scalar_file('afm',nvar,.false.,id)
          id_out(3) = id
          call shyfem_init_scalar_file('afm_d',nst_bfm_d,.false.,id)
          id_out(4) = id
          call shyfem_init_scalar_file('bfm_b',nst_bfm_d,.false.,id)
          id_out(5) = id
        end if

	end subroutine

!*********************************************************************

	subroutine bfm_write_output_file
        use levels, only: ilhkv
	use basin, only : nkn,nel,ngr,mbw
        use mod_layer_thickness, only: hdknv

	implicit none

	integer id,nvar,i,k,idbase,idc
        real cmin,cmax,ctot
	real v1v(nkn)
        integer bottom
        real vint,norm
	double precision dtime

	logical next_output,next_output_d

	if( ibfm < 0 ) return

!-------------------------------------------------------------
! write to file
!-------------------------------------------------------------

	idbase = 600
	nvar   = ibfm_state
        if(ieco_av==1)then
          ave_bfmv    = ( ave_bfmv   * counter_ave + bfmv)  /(counter_ave +1.0)
          ave_bfm_d  = ( ave_bfm_d * counter_ave + bfm_d)/(counter_ave +1.0)
          ave_bfm_benthic  = ( ave_bfm_benthic * counter_ave + bfm_benthic)/(counter_ave +1.0)
        elseif(ieco_av==2)then
          DO k=1,nkn

           bottom=ilhkv(k)
           norm=1.0/sum(hdknv(1:bottom,k)) ! Depth in k  

           DO i=1,nvar 
              vint=sum(bfmv(1:bottom,k,i)*hdknv(1:bottom,k))*norm
              ave_bfmv(1,k,i)   = ( ave_bfmv(1,k,i)  * counter_ave + vint)  /(counter_ave +1.0)
           ENDDO

           DO i=1,nst_bfm_d
              vint=sum(bfm_d(1:bottom,k,i)*hdknv(1:bottom,k))*norm
              ave_bfm_d(1,k,i)  = ( ave_bfm_d(1,k,i) * counter_ave + vint)/(counter_ave +1.0)
           ENDDO

           DO i=1,ibfm_state_benthic
              vint=sum(bfm_benthic(1:bottom,k,i)*hdknv(1:bottom,k))*norm
              ave_bfm_benthic(1,k,i)  = ( ave_bfm_benthic(1,k,i) * counter_ave + vint)/(counter_ave +1.0)
           ENDDO

          ENDDO
        endif

        counter_ave = counter_ave  + 1.0 
        write(*,*) 'counter_ave', counter_ave

        if( next_output_d(da_out) ) then
	  call get_act_dtime(dtime)
          id = nint(id_out(1))
	  do i=1,nvar
	    idc = idbase + i
            call shy_write_scalar_record(id,dtime,idc,nlvdi,bfmv(1,1,i))
	  end do
	  idbase = 690 
          id = nint(id_out(2))
	  do i=1,nst_bfm_d
	    idc = idbase + i
           call shy_write_scalar_record(id,dtime,idc,nlvdi,bfm_d(1,1,i))
	  end do
          if(ieco_av>01)then
!           average values
            idbase = 900
            id = nint(id_out(3))
            if(ieco_av==1)then
              do i=1,nvar
                idc = idbase + i
                call shy_write_scalar_record(id,dtime,idc,nlvdi,ave_bfmv(1,1,i))
              end do
            elseif(ieco_av==2)then
              do i=1,nvar
                idc = idbase + i
                call shy_write_scalar_record2d(id,dtime,idc,ave_bfmv(1,1,i))
              end do
            endif

            idbase = 990
            id = nint(id_out(4))
            if(ieco_av==1)then
              do i=1,nst_bfm_d
                idc = idbase + i
                call shy_write_scalar_record(id,dtime,idc,nlvdi,ave_bfm_d(1,1,i))
              end do
            elseif(ieco_av==2)then
              do i=1,nst_bfm_d
                idc = idbase + i
                call shy_write_scalar_record2d(id,dtime,idc,ave_bfm_d(1,1,i))
              end do
            endif

            idbase = 980
            id = nint(id_out(5))
            do i=1,nst_bfm_d
              idc = idbase + i
              call shy_write_scalar_record2d(id,dtime,idc,ave_bfm_benthic(1,1,i))
            end do
           
            ave_bfmv    = 0.0
            ave_bfm_d   = 0.0
            ave_bfm_benthic   = 0.0
          endif
          counter_ave = 0.0 
        end if

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end subroutine

!*********************************************************************

	subroutine bfm_check(string)

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) string

	integer iv,nvar,k,l,lmax
	real vmin,vmax,v
	real bmin,bmax
	character*80 textgen,text,line

	nvar = ibfm_state
	vmin = -1000.
	vmax = +1000.
	textgen = 'NaN check'

	write(6,*) 'bfm_check: ' // string

	do iv=1,nvar
	  write(line,'(i3)') iv
	  text = 'bfmv ' // trim(line) // '  ' // trim(string)
	  call check2Dr(nlvdi,nlv,nkn,bfmv,vmin,vmax,trim(textgen),trim(text))
	  do k=1,nkn
	    lmax = ilhkv(k)
	    do l=1,lmax
	      v = bfmv(l,k,iv)
	      if( v <= 0 ) then
		!write(6,*) 'bfm_check: ',iv,k,l,v
	      end if
	    end do
	  end do

	  call mima3d(nkn,nlvdi,ilhkv,bfmv(:,:,iv),bmin,bmax)
	  write(6,*) iv,bmin,bmax
	end do

	write(66,*) 'bfm_check: ',string
	write(66,*) (bfmv(1,1,iv),iv=1,nvar)

	end subroutine bfm_check

!*********************************************************************

        subroutine bfm_init_file(dtime,nvar,nlvddi,nlv,nkn,val0,val)

! initialization of bfm from file

        implicit none

	double precision dtime
        integer nvar
        integer nlvddi
        integer nlv
        integer nkn
	real val0(nvar)
        real val(nlvddi,nkn,nvar)

        call tracer_file_init('bfm init','bfmini',dtime,nvar,nlvddi,nlv,nkn,val0,val)

	end subroutine bfm_init_file

!*********************************************************************

      subroutine bfm_reactor
      
      implicit none
      
      if( ibfm < 0 ) return

      call bfm_reactor_internal

      call bfm_benthic_internal

      end subroutine bfm_reactor

!==================================================================
        end module mod_bfm
!==================================================================

