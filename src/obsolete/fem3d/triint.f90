
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998,2010,2013,2019  Georg Umgiesser
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

! subroutines for 2 (bi) and 3-dimensional (tri) linear interpolation
!
! contents :
!
! function triint(u,x,y,z)	tri-linear interpolation in cube
! function bilint(u,x,y)	bi-linear interpolation in square
! subroutine exbil_test		exercise bi-linear interpolation
! subroutine extri_test		exercise tri-linear interpolation
! function grand_tri()		random numbers
!
! revision log :
!
! 21.08.1998	ggu	tri and bi linear interpolation
! 23.03.2010	ggu	changed v6.1.1
! 25.10.2013	ggu	changed VERS_6_1_68
! 16.02.2019	ggu	changed VERS_7_5_60
!
!***********************************************************************

	function triint(u,x,y,z)

! tri-linear interpolation in cube
!
! (x,y,z) are in the range [0...1]
! u is array of 8 values distributed in the following way
!
!
!
!                  z                      y
!                  ^                     ^
!                  |                    /
!                  |                   /
!                  |     7            /            8
!                  |       *---------------------*
!                  |      /|        /           /|
!                  |     / |       /           / |
!                  |    /  |      /           /  |
!                  |   /   |     /           /	 |                   z = 1
!                  |  /    |    /           /    |
!                  | /     |   /           /     |
!                  |/      |  /           /      |
!                5 *---------------------* 6     |
!                  |       |/            |       |
!                  |     3 *-------------|-------* 4
!                  |      /              |      /
!                  |     /               |     /
!                  |    /                |    /
!                  |   /                 |   /	                     z = 0
!                  |  /                  |  /
!                  | /                   | /
!                  |/                    |/
!                  *---------------------*-----------------------> x
!                1                         2
!

	implicit none

	real triint
	real u(8)
	real x,y,z	![0...1]

	integer i,j
	real acu,val,xy
	real a(8)

	real atri(8,8)
	save atri
	data atri / &
     &			 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , &
     &			-1 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , &
     &			-1 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , &
     &			 1 ,-1 ,-1 , 1 , 0 , 0 , 0 , 0 , &
     &			-1 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , &
     &			 1 ,-1 , 0 , 0 ,-1 , 1 , 0 , 0 , &
     &			 1 , 0 ,-1 , 0 ,-1 , 0 , 1 , 0 , &
     &			-1 , 1 , 1 ,-1 , 1 ,-1 ,-1 , 1  &
     &		 /

	do i=1,8
	  acu = 0.
	  do j=1,8
	    acu = acu + atri(j,i) * u(j)
	  end do
	  a(i) = acu
	end do

	xy = x * y

	val = a(1) + a(2) * x + a(3) * y + a(4) * xy &
     &		+ z * ( a(5) + a(6) * x + a(7) * y + a(8) * xy )

	triint = val

	end

!*************************************************************

	function bilint(u,x,y)

! bi-linear interpolation in square
!
! (x,y)		normalized coordinates in the range [0...1]
! u		array of 4 values at the vertices of the
!			the square distributed in the following way
! bilint	interpolated value of u onto (x,y) on return
!
!
!
!                  y
!                  ^
!                  |
!                  |
!                  |
!                  |
!                3 |                       4
!                  *---------------------*
!                  |                     |
!                  |                     |
!                  |                     |
!                  |                     |
!                  |                     |
!                  |                     |
!                  |                     |
!                  |                     |
!                  |                     |
!                  *---------------------*-----------------------> x
!                1                         2
!
!

	implicit none

	real bilint
	real u(4)
	real x,y	![0...1]

	integer i,j
	real acu,val
	real a(4)

	real atri(4,4)
	save atri
	data atri / &
     &			 1 , 0 , 0 , 0 , &
     &			-1 , 1 , 0 , 0 , &
     &			-1 , 0 , 1 , 0 , &
     &			 1 ,-1 ,-1 , 1 &
     &		 /

	do i=1,4
	  acu = 0.
	  do j=1,4
	    acu = acu + atri(j,i) * u(j)
	  end do
	  a(i) = acu
	end do

	val = a(1) + a(2) * x + a(3) * y + a(4) * x * y

	bilint = val

	end

!*************************************************************
!*************************************************************
!*************************************************************

	subroutine exbil_test

! exercise bi-linear interpolation

	implicit none

	real u(4)
	real x,y
	real r
	real bilint

	write(6,*) 'Enter values at 4 vertices :'
	write(6,*) '(0,0)  (1,0)  (0,1)  (1,1)'
	read(5,*) u
	write(6,*) u

	x=0.
	y=0.

	do while( x .ge. 0. .and. y .ge. 0. )
	  write(6,*) 'Enter x,y [0..1] (negative to end) :'
	  read(5,*) x,y
	  r = bilint(u,x,y)
	  write(6,*) x,y,r
	end do

	end
	  
!*************************************************************

	subroutine extri_test

! exercise tri-linear interpolation

	implicit none

	integer ndim
	parameter ( ndim = 30 )

	logical berror
	integer n
	integer i,j
	integer ix,iy,iz
	real u(8)
	real uvw(0:2)
	real x,y,z
	real val,compare
	real triint,grand_tri

	real xyz(8,3)
	data xyz / &
     &			 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , &
     &			 0 , 0 , 1 , 1 , 0 , 0 , 1 , 1 , &
     &			 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1   &
     &		 /

! interpolation property

	berror = .false.
	write(6,*) 'Checking interpolation property...'

	do i=1,8
	  u(i) = 0.
	end do

	do i=1,8
	  u(i) = 1.
	  do j=1,8
	    x = xyz(j,1)
	    y = xyz(j,2)
	    z = xyz(j,3)
	    val = triint(u,x,y,z)
	    compare = 0.
	    if( i .eq. j ) compare = 1.
	    if( val .ne. compare ) then
		write(6,*) i,x,y,z,val
		berror = .true.
	    end if
	  end do
	  u(i) = 0.
	end do

	if( berror ) then
	  write(6,*) 'There have been errors in check...'
	  stop 'error stop'
	else
	  write(6,*) '   ...passed'
	end if

! linear interpolation

	do j=0,2

	do i=1,8
	  if( xyz(i,j) .ne. 0. ) u(i) = 1.
	end do

	uvw(j) = 0.5
	do n=1,ndim
	  uvw( mod(j+1,3) ) = grand_tri()
	  uvw( mod(j+2,3) ) = grand_tri()
	  x = uvw(0)
	  y = uvw(1)
	  z = uvw(2)
	  val = triint(u,x,y,z)
	  write(6,*) n,x,y,z,val
	end do

	end do

	end

!*************************************************************

	function grand_tri()

! random numbers

	implicit none

	real grand_tri
	real rand

	integer iseed
	save iseed
	data iseed /37582629/

	grand_tri = rand(iseed)

	end 

!*************************************************************
!
! uncomment next lines to run test routines

!	program hptri_main
!	call extri_test
!	call exbil_test
!	end

!*************************************************************

