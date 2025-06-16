
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2020  Georg Umgiesser
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

! routines for generic concentration
!
! contents :
!
! subroutine mod_conz_init(ncs,nkn,nlv)
!
! revision log :
!
! 10.07.2015	ggu	changed VERS_7_1_50
! 20.07.2015	ggu	changed VERS_7_1_81
! 09.11.2015	ggu	restructured with global values
! 19.02.2016	ggu	set default decay to zero
! 07.06.2016	ggu	changed VERS_7_5_12
! 09.09.2016	ggu	changed VERS_7_5_17
! 13.04.2017	ggu	changed VERS_7_5_25
! 05.12.2017	ggu	changed VERS_7_5_39
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.10.2018	ggu	changed VERS_7_5_50
! 16.02.2019	ggu	changed VERS_7_5_60
! 20.03.2020	ggu	restart routines introduced
! 22.06.2021	ggu	new parameter bage for age computation
! 13.10.2021	ggu	restart for mpi - compatibility with older versions
! 28.04.2023	ggu	update function calls for belem
! 16.04.2025	ggu	contau eliminated, added wsettlv
!
!******************************************************************

!==================================================================
        module mod_conz
!==================================================================

        implicit none

        integer, private, save :: ncs_conz = 0
        integer, private, save :: nkn_conz = 0
        integer, private, save :: nlv_conz = 0

        real, allocatable, save :: conzv(:,:,:)
        real, allocatable, save :: cnv(:,:)

        logical, save :: baccum = .false.
        logical, save :: bage = .false.
	double precision, save :: dtconz_accum
        double precision, allocatable, save :: conz_aver(:,:,:)
        real, allocatable, save :: conz_min(:,:,:)
        real, allocatable, save :: conz_max(:,:,:)

	integer, save :: iconz = 0
        integer, save :: icall_conz = 0
        integer, save :: iuinfo = 0
        integer, save :: level = 0
        integer, save :: iprogr = 0
        integer, save :: idecay = 0

        logical, save :: binfo = .true.

        real, save :: cref,rkpar,difmol

        integer, save, allocatable :: idconz(:)
        double precision, save :: da_out(4)

        real, save, allocatable :: cdefs(:)
        real, save, allocatable :: tauv(:)
        real, save, allocatable :: wsettlv(:)
        real, save, allocatable :: massv(:)

	character*4, save :: what = 'conz'

!==================================================================
	contains
!==================================================================

        subroutine mod_conz_init(ncs,nkn,nlv)

        integer ncs
        integer nkn
        integer nlv

        if( ncs == ncs_conz .and. nkn == nkn_conz .and. nlv == nlv_conz ) return

        if( ncs > 0 .or. nkn > 0 .or. nlv > 0 ) then
          if( ncs == 0 .or. nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'ncs,nkn,nlv: ',ncs,nkn,nlv
            stop 'error stop mod_conz_init: incompatible parameters'
          end if
        end if

        if( nkn_conz > 0 ) then
          deallocate(conzv)
          deallocate(cnv)
          deallocate(conz_aver)
          deallocate(conz_min)
          deallocate(conz_max)
        end if

        ncs_conz = ncs
        nkn_conz = nkn
        nlv_conz = nlv

        if( nkn == 0 ) return

        allocate(conzv(nlv,nkn,ncs))
        allocate(cnv(nlv,nkn))

	cnv = 0.
	conzv = 0.

	if( baccum ) then
          allocate(conz_aver(nlv,nkn,ncs))
          allocate(conz_min(nlv,nkn,ncs))
          allocate(conz_max(nlv,nkn,ncs))
	end if

        end subroutine mod_conz_init

!*****************************************************

        function mod_conz_is_initialized()

        logical mod_conz_is_initialized

        mod_conz_is_initialized = ( nkn_conz > 0 )

        end function mod_conz_is_initialized

!==================================================================
        end module mod_conz
!==================================================================

!*********************************************************************
!*********************************************************************
!*********************************************************************
! restart routines
!*********************************************************************
!*********************************************************************
!*********************************************************************

        subroutine write_restart_conz(iunit,ic)

        use basin
        use levels
        use mod_conz
	use mod_restart

        implicit none

        integer iunit,ic

        integer k,l,i
	logical, parameter :: belem = .false.

        if( ic .eq. 1 ) then
	  call restart_write_value(iunit,belem,cnv)
        else
          !write(iunit) (((conzv(l,k,i),l=1,nlv),k=1,nkn),i=1,iconz)
	  do i=1,ic
	    call restart_write_value(iunit,belem,conzv(:,:,i))
	  end do
        end if

        end

!*********************************************************************

        subroutine skip_restart_conz(iunit,ic)

        use mod_conz
	use mod_restart

        implicit none

        integer iunit,ic

	integer i

	if( nvers_rst < 16 ) then
	  read(iunit)
	else
	  do i=1,ic
            read(iunit)
	  end do
	end if

        end

!*********************************************************************

        subroutine read_restart_conz(iunit,ic)

        use basin
        use levels
        use mod_conz
	use mod_restart
	use shympi

        implicit none

        integer iunit,ic

        integer k,l,i
	logical, parameter :: belem = .false.
	real, allocatable :: conz_aux(:,:,:)

        call mod_conz_init(ic,nkn,nlvdi)

        if( ic .eq. 1 ) then
	  call restart_read_value(iunit,belem,cnv)
        else
          !read(iunit) (((conzv(l,k,i),l=1,nlv),k=1,nkn),i=1,ic)
	  if( nvers_rst < 16 ) then	!hack for compatibility
	    allocate(conz_aux(nlv_global,nkn_global,ic))
	    read(iunit) conz_aux
	    do i=1,ic
	      call shympi_g2l_array(conz_aux(:,:,i),conzv(:,:,i))
	    end do
	  else
	    do i=1,ic
	      call restart_read_value(iunit,belem,conzv(:,:,i))
	    end do
	  end if
        end if

        end

!*********************************************************************

