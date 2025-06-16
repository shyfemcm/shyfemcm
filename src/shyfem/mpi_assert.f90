
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2019  Georg Umgiesser
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
! 07.12.2017	ggu	changed VERS_7_5_40
! 19.04.2018	ggu	changed VERS_7_5_45
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 01.04.2022	ggu	some more checkes
! 01.04.2022	ggu	in shympi_check_hydro_print() check for wprv was wrong

!*****************************************************************
!
! checks if all assumptions on variables are true
!
!*****************************************************************

	subroutine shympi_assert_all

	call shympi_assert_coriolis

	end

!*****************************************************************

	subroutine shympi_assert_coriolis

	use mod_internal
	use shympi

	implicit none

	integer isphe
	integer vals(n_threads)

	call get_coords_ev(isphe)

	call shympi_gather(isphe,vals)

	if( shympi_is_master() ) then
	  if( any(vals/=isphe) ) then
	    write(6,*) 'error in isphe: ',isphe,vals
	    stop 'error stop shympi_assert_coriolis: isphe'
	  end if
	end if

	call shympi_check_2d_elem(fcorv,'fcorv')

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine shympi_check_all(iwhat)

	implicit none

	integer iwhat

	logical, parameter :: bdebug = .false.

	if( bdebug ) then
	  write(6,*) 'shympi_check_all: checking arrays for correctness'
	  write(6,*) 'call at iwhat = ',iwhat
	end if

	call shympi_assert_all
	call shympi_check_all_static
	call shympi_check_all_scalar
	call shympi_check_all_dynamic

	if( bdebug ) then
	  write(6,*) 'arrays ok at iwhat = ',iwhat
	end if

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine shympi_check_all_static

	implicit none

	call shympi_check_depth
	call shympi_check_geom_static
	call shympi_check_levels

	end

!*****************************************************************

	subroutine shympi_check_all_dynamic

	implicit none

	call shympi_check_geom_dynamic

	call shympi_check_hydro
	call shympi_check_hydro_baro
	call shympi_check_hydro_vel
	call shympi_check_hydro_print

	end

!*****************************************************************

	subroutine shympi_check_all_scalar

        use mod_conz
        use mod_ts
        use shympi

	implicit none

        if( mod_conz_is_initialized() ) then
	  call shympi_check_3d_node(cnv,'cnv')
        end if
	call shympi_check_3d_node(saltv,'saltv')
	call shympi_check_3d_node(tempv,'tempv')
	call shympi_check_3d_node(rhov,'rhov')

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine shympi_check_hydro

	use basin
	use mod_hydro
	use mod_layer_thickness
	use shympi

	implicit none

	integer i
	real aux(nel)
	integer ie,iee,ieee,iei,ieext,ieint
	logical bunique

	call shympi_check_2d_node(zov,'zov')
	call shympi_check_2d_node(znv,'znv')

	do i=1,3
	  aux = zeov(i,:)
	  call shympi_check_2d_elem(aux,'zeov')
	  aux = zenv(i,:)
	  call shympi_check_2d_elem(aux,'zenv')
	end do

	!write(6,*) nlv_global,nkn_global,nel_global,nel,my_id
	!write(6,*) 'nel: ',nel,nel_unique,nel_inner,my_id
	!if( my_id == 2 ) then
	!  ie = 132106
	!  iee = ieext(ie)
	!  write(6,*) 'ieext: ',ie,iee
	!end if
	!ieee = 159344
	ieee = 0

	iei = ieint(ieee)
	if( iei /= 0 ) then
	  bunique = shympi_is_unique_elem(iei)
	  write(6,*) ieee,iei,utlnv(1,iei),my_id
	  write(6,*) bunique,id_elem(:,iei),my_id
	end if

	call shympi_check_3d_elem(hdeov,'hdeov')
	call shympi_check_3d_node(hdkov,'hdkov')
	call shympi_check_3d_elem(hdenv,'hdenv')
	call shympi_check_3d_node(hdknv,'hdknv')

	call shympi_check_3d_elem(utlov,'utlov')
	call shympi_check_3d_elem(vtlov,'vtlov')
	call shympi_check_3d_elem(utlnv,'utlnv')
	call shympi_check_3d_elem(vtlnv,'vtlnv')

	end

!*****************************************************************

	subroutine shympi_check_hydro_baro

	use mod_hydro_baro
	use shympi

	implicit none

	call shympi_check_2d_elem(unv,'unv')
	call shympi_check_2d_elem(vnv,'vnv')
	call shympi_check_2d_elem(uov,'uov')
	call shympi_check_2d_elem(vov,'vov')

	end

!*****************************************************************

	subroutine shympi_check_hydro_vel

	use basin
	use levels
	use mod_hydro_vel
	use shympi

	implicit none

	call shympi_check_3d_elem(ulov,'ulov')
	call shympi_check_3d_elem(vlov,'vlov')
	call shympi_check_3d_elem(ulnv,'ulnv')
	call shympi_check_3d_elem(vlnv,'vlnv')
	call shympi_check_3d0_node(wlnv,'wlnv')
	call shympi_check_3d0_node(wlov,'wlov')

	end

!*****************************************************************

	subroutine shympi_check_hydro_print

	use basin
	use levels
	use mod_hydro_print
	use shympi

	implicit none

	integer i
	real aux(nkn)

	call shympi_check_3d_node(uprv,'uprv')
	call shympi_check_3d_node(vprv,'vprv')
	call shympi_check_3d_node(upro,'upro')
	call shympi_check_3d_node(vpro,'vpro')
	call shympi_check_3d_node(wprv,'wprv')
	call shympi_check_2d_node(up0v,'up0v')
	call shympi_check_2d_node(vp0v,'vp0v')

	!do i=1,3
	!  aux = xv(i,:)
	!  call shympi_check_2d_node(aux,'xv')
	!end do

	end

!*****************************************************************

	subroutine shympi_check_depth

	use basin
	use mod_depth
	use shympi

	implicit none

	integer i
	real aux(nel)

	do i=1,3
	  aux = hm3v(i,:)
	  call shympi_check_2d_elem(aux,'hm3v')
	end do

	call shympi_check_2d_elem(hev,'hev')
	call shympi_check_2d_node(hkv,'hkv')
	call shympi_check_2d_node(hkv_min,'hkv_min')
	call shympi_check_2d_node(hkv_max,'hkv_max')

	end

!*****************************************************************

	subroutine shympi_check_geom_dynamic

	use evgeom
	use mod_geom_dynamic
	use mod_diff_visc_fric
	use shympi

	implicit none

	call shympi_check_2d_elem(czv,'czv')
	call shympi_check_2d_elem(iwegv,'iwegv')
	call shympi_check_2d_elem(iwetv,'iwetv')
	call shympi_check_2d_node(inodv,'inodv')	!FIXME - not working

	end

!*****************************************************************

	subroutine shympi_check_geom_static

	use basin
	use evgeom
	use mod_geom
	use shympi

	implicit none

	integer i
	real aux(nel)

	do i=1,evdim
	  aux = ev(i,:)
	  call shympi_check_2d_elem(aux,'ev')
	end do
	call shympi_check_2d_elem(iarv,'iarv')

	end

!*****************************************************************

	subroutine shympi_check_levels

	use levels
	use shympi

	implicit none

	call shympi_check_2d_elem(ilhv,'ilhv')
	call shympi_check_2d_elem(ilmv,'ilmv')
	call shympi_check_2d_node(ilhkv,'ilhkv')
	call shympi_check_2d_node(ilmkv,'ilmkv')

	end

!*****************************************************************

