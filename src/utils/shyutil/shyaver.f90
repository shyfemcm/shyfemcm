
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

! implements average routines for hydro and scalar variables
!
! revision log :
!
! 18.10.2024	ggu	shyaver started

!**************************************************************************

!==================================================================
        module mod_shyaver
!==================================================================

        implicit none

	private

        logical, save :: bdebug = .false.

        type :: scal_entry
	  integer			:: idim
	  double precision		:: time
	  double precision, allocatable :: array1(:)
	  double precision, allocatable :: array2(:,:)
	  double precision, allocatable :: array3(:,:,:)
        end type scal_entry

        type :: hydro_entry
	  double precision		:: time
	  double precision, allocatable :: znv(:)
	  double precision, allocatable :: zenv(:,:)
	  double precision, allocatable :: utlnv(:,:)
	  double precision, allocatable :: vtlnv(:,:)
        end type hydro_entry

        integer, parameter :: idmax = 10
        integer, save :: idlast = 0
        integer, save :: ihlast = 0

	type(scal_entry), save :: pentry(idmax)
	type(hydro_entry), save :: hentry(idmax)

	INTERFACE shyaver_initialize
	MODULE PROCEDURE  shyaver_initialize_1 &
     &			, shyaver_initialize_2 &
     &			, shyaver_initialize_3 &
     &			, shyaver_initialize_hydro
	END INTERFACE

	INTERFACE shyaver_accumulate
	MODULE PROCEDURE  shyaver_accumulate_1 &
     &			, shyaver_accumulate_2 &
     &			, shyaver_accumulate_3 &
     &			, shyaver_accumulate_hydro
	END INTERFACE

	INTERFACE shyaver_average
	MODULE PROCEDURE  shyaver_average_1 &
     &			, shyaver_average_2 &
     &			, shyaver_average_3 &
     &			, shyaver_average_hydro
	END INTERFACE

	!----------------------------------------------------------
	! following are the public calls for shyaver routines
	!----------------------------------------------------------

	public :: shyaver_initialize
	public :: shyaver_accumulate
	public :: shyaver_average

!==================================================================
	contains
!==================================================================

	subroutine allocate_new_scalar(id)

	integer, intent(out) :: id

	idlast = idlast + 1
	if( idlast > idmax ) then
	  write(6,*) idlast,idmax
	  stop 'error stop allocate_new_scalar: idlast>idmax'
	end if

	id = idlast

	end

!******************************************************************

	subroutine allocate_new_hydro(id)

	integer, intent(out) :: id

	ihlast = ihlast + 1
	if( ihlast > idmax ) then
	  write(6,*) ihlast,idmax
	  stop 'error stop allocate_new_hydro: ihlast>idmax'
	end if

	id = ihlast

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shyaver_initialize_1(array,id)

! initialzes 1d array

	use basin

	real, intent(in) ::  array(:)
	integer, intent(out) :: id

	integer n1

	call allocate_new_scalar(id)
	pentry(id)%idim = 1

	n1 = size(array)
	allocate(pentry(id)%array1(n1))

	pentry(id)%time = 0.
	pentry(id)%array1 = 0.

	end

!******************************************************************

	subroutine shyaver_initialize_2(array,id)

! initialzes 2d array

	use basin

	real, intent(in) ::  array(:,:)
	integer, intent(out) :: id

	integer n1,n2

	call allocate_new_scalar(id)
	pentry(id)%idim = 2

	n1 = size(array,1)
	n2 = size(array,2)
	allocate(pentry(id)%array2(n1,n2))

	pentry(id)%time = 0.
	pentry(id)%array2 = 0.

	end

!******************************************************************

	subroutine shyaver_initialize_3(array,id)

! initialzes 2d array

	use basin

	real, intent(in) ::  array(:,:,:)
	integer, intent(out) :: id

	integer n1,n2,n3

	call allocate_new_scalar(id)
	pentry(id)%idim = 3

	n1 = size(array,1)
	n2 = size(array,2)
	n3 = size(array,3)
	allocate(pentry(id)%array3(n1,n2,n3))

	pentry(id)%time = 0.
	pentry(id)%array3 = 0.

	end

!******************************************************************

	subroutine shyaver_initialize_hydro(znv,zenv,utlnv,vtlnv,id)

! initialzes 2d array

	use basin

	real, intent(in) ::  znv(:)
	real, intent(in) ::  zenv(:,:)
	real, intent(in) ::  utlnv(:,:)
	real, intent(in) ::  vtlnv(:,:)
	integer, intent(out) :: id

	integer nk,n1,n2,n3

	call allocate_new_hydro(id)
	pentry(id)%idim = 4

	nk = size(znv,1)
	n1 = size(zenv,1)
	n2 = size(zenv,2)
	n3 = size(utlnv,1)
	allocate(hentry(id)%znv(nk))
	allocate(hentry(id)%zenv(n1,n2))
	allocate(hentry(id)%utlnv(n3,n2))
	allocate(hentry(id)%vtlnv(n3,n2))

	hentry(id)%time = 0.
	hentry(id)%znv = 0.
	hentry(id)%zenv = 0.
	hentry(id)%utlnv = 0.
	hentry(id)%vtlnv = 0.

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shyaver_accumulate_1(dt,array,id)

! accumulate 1d array

	use basin

	real, intent(in) ::  dt
	real, intent(in) ::  array(:)
	integer, intent(in) ::  id

	pentry(id)%time = pentry(id)%time + dt
	pentry(id)%array1 = pentry(id)%array1 + array * dt

	end

!******************************************************************

	subroutine shyaver_accumulate_2(dt,array,id)

! accumulate 2d array

	use basin

	real, intent(in) ::  dt
	real, intent(in) ::  array(:,:)
	integer, intent(in) ::  id

	pentry(id)%time = pentry(id)%time + dt
	pentry(id)%array2 = pentry(id)%array2 + array * dt

	end

!******************************************************************

	subroutine shyaver_accumulate_3(dt,array,id)

! accumulate 2d array

	use basin

	real, intent(in) ::  dt
	real, intent(in) ::  array(:,:,:)
	integer, intent(in) ::  id

	pentry(id)%time = pentry(id)%time + dt
	pentry(id)%array3 = pentry(id)%array3 + array * dt

	end

!******************************************************************

	subroutine shyaver_accumulate_hydro(dt,znv,zenv,utlnv,vtlnv,id)

! accumulate 2d array

	use basin

	real, intent(in) ::  dt
	real, intent(in) ::  znv(:)
	real, intent(in) ::  zenv(:,:)
	real, intent(in) ::  utlnv(:,:)
	real, intent(in) ::  vtlnv(:,:)
	integer, intent(in) ::  id

	hentry(id)%time = hentry(id)%time + dt
	hentry(id)%znv = hentry(id)%znv + znv * dt
	hentry(id)%zenv = hentry(id)%zenv + zenv * dt
	hentry(id)%utlnv = hentry(id)%utlnv + utlnv * dt
	hentry(id)%vtlnv = hentry(id)%vtlnv + vtlnv * dt

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shyaver_average_1(array,id)

! averages 1d array

	use basin

	real, intent(out) :: array(:)
	integer, intent(in) ::  id

	array = pentry(id)%array1 / pentry(id)%time

	pentry(id)%time = 0.
	pentry(id)%array1 = 0.

	end

!******************************************************************

	subroutine shyaver_average_2(array,id)

! averages 2d array

	use basin

	real, intent(out) :: array(:,:)
	integer, intent(in) ::  id

	array = pentry(id)%array2 / pentry(id)%time

	pentry(id)%time = 0.
	pentry(id)%array2 = 0.

	end

!******************************************************************

	subroutine shyaver_average_3(array,id)

! averages 2d array

	use basin

	real, intent(out) ::  array(:,:,:)
	integer, intent(in) ::  id

	array = pentry(id)%array3 / pentry(id)%time

	pentry(id)%time = 0.
	pentry(id)%array3 = 0.

	end

!******************************************************************

	subroutine shyaver_average_hydro(znv,zenv,utlnv,vtlnv,id)

! averages 2d array

	use basin

	real, intent(out) ::  znv(:)
	real, intent(out) ::  zenv(:,:)
	real, intent(out) ::  utlnv(:,:)
	real, intent(out) ::  vtlnv(:,:)
	integer, intent(in) ::  id

	znv = hentry(id)%znv / hentry(id)%time
	zenv = hentry(id)%zenv / hentry(id)%time
	utlnv = hentry(id)%utlnv / hentry(id)%time
	vtlnv = hentry(id)%vtlnv / hentry(id)%time

	hentry(id)%time = 0.
	hentry(id)%znv = 0.
	hentry(id)%zenv = 0.
	hentry(id)%utlnv = 0.
	hentry(id)%vtlnv = 0.

	end

!==================================================================
        end module mod_shyaver
!==================================================================

