
!--------------------------------------------------------------------------
!
!    Copyright (C) 2019-2020  Georg Umgiesser
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

! routines for BFM module - dummy routine
!
! revision log :
!
! 02.10.2019	ggu	dummy routine introduced
! 22.04.2020	ggu	error if called
! 23.10.2024	ggu	define ibfm in module

!******************************************************************

!==================================================================
	module mod_bfm
!==================================================================

	implicit none

	integer, save :: ibfm = 0

!==================================================================
	contains
!==================================================================

	subroutine bfm_init
	real getpar
	ibfm = nint(getpar('ibfm'))
	call bfm_init_internal(ibfm)
	end

	subroutine bfm_run
	end

!==================================================================
	end module mod_bfm
!==================================================================

!******************************************************************

	subroutine bfm_init_internal(ibfm)

	implicit none

	integer ibfm

	if( ibfm > 0 ) then
	  write(6,*) 'ibfm = ',ibfm
	  write(6,*) 'bfm routines have not been linked into shyfem'
	  write(6,*) 'please correct Rules.make file'
	  stop 'error stop bfm_init_internal: no bfm routines'
	end if

	end

	subroutine bfm_reactor_internal
	end

        subroutine write_restart_bfm(iunit)
        implicit none
        integer iunit
        end

        subroutine skip_restart_bfm(iunit)
        implicit none
        integer iunit
        end

        subroutine read_restart_bfm(iunit)
        implicit none
        integer iunit
        end

        subroutine bfm_init_for_restart()
        implicit none
        end

!******************************************************************

	subroutine bfm_shyfem_array_interface

! pointers to shyfem variables needed in bfm/vegetation
!
! no area yet

	use mod_ts
	use mod_meteo
	use mod_diff_visc_fric
	use mod_layer_thickness
	use levels

	implicit none

	real, pointer :: temp_bfm(:,:)
	real, pointer :: salt_bfm(:,:)
	real, pointer :: rho_bfm(:,:)		!density anomaly, r = r' + 1025
	real, pointer :: windx_bfm(:)
	real, pointer :: windy_bfm(:)
	real, pointer :: wspeed_bfm(:)		!wind speed
	real, pointer :: solrad_bfm(:)
	real, pointer :: tair_bfm(:)
	real, pointer :: ice_bfm(:)
	real, pointer :: visc_bfm(:,:)		!(0:nlv,nkn)
	real, pointer :: diff_bfm(:,:)		!(0:nlv,nkn)
	real, pointer :: bnstress_bfm(:)	!normalized bottom stress
	real, pointer :: dz_bfm(:,:)		!layer thickness at node
	integer, pointer :: k_bott_bfm(:)	!number of layers at node
	integer, pointer :: k_surf_bfm(:)	!index of surface layer at node

	temp_bfm => tempv
	salt_bfm => saltv
	rho_bfm => rhov
	windx_bfm => wxv
	windy_bfm => wyv
	wspeed_bfm => metws
	solrad_bfm => metrad
	tair_bfm => mettair
	ice_bfm => metice
	visc_bfm => visv
	diff_bfm => difv
	bnstress_bfm => bnstressv
	dz_bfm => hdknv
	k_bott_bfm => ilhkv
	k_surf_bfm => jlhkv

	end subroutine

!******************************************************************

