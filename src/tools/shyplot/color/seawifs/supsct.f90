
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016,2019  Georg Umgiesser
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

!  color handling routines

! ************************************************************

! revision log :
!
! 17.06.2016	ggu	changed VERS_7_5_15
! 14.02.2019	ggu	changed VERS_7_5_56

! ************************************************************


	subroutine sctini

!  initializes color

	implicit none

	include 'vsect.h'

	logical binit
	save binit
	data binit /.false./

	if( binit ) return

	binit = .true.

        nnode = 0

	end

! ************************************************************

	subroutine sctrd

!  reads color block from str file

	implicit none

	include 'vsect.h'

	character*80 name,text
	logical bdebug
	integer mode
	integer i
	integer ndim,ndim1
	real value

	integer nrdpar

	call sctini

	ndim = isodim
	ndim1 = isodim + 1
	bdebug = .false.
	ncolrd = 0
	nisord = 0

        mode = 1
        do while( mode .ne. 0 )
            mode = nrdpar('color',name,value,text)
            if( mode .eq. 2 ) then
	      if( name .ne. 'vsect' ) then
		write(6,*) 'Variable is no vector : ',name
		stop 'error stop : sctrd'
	      end if
	    end if
            if( name .eq. 'color' ) then
		if( ncolrd .ge. ndim1 ) stop 'error stop sctrd: dimension'
                ncolrd = ncolrd + 1
                ciso(ncolrd) = value
            else if( name .eq. 'isoval' ) then
		if( nisord .ge. ndim ) stop 'error stop sctrd: dimension'
                nisord = nisord + 1
                fiso(nisord) = value
	    end if
        end do

	if( bdebug ) then
	  write(6,*) 'color debug : ',nisord,ncolrd
	  write(6,*) (fiso(i),i=1,nisord)
	  write(6,*) (ciso(i),i=1,ncolrd)
	end if

	end

! ************************************************************
! ************************************************************
! ************************************************************
!  accessor routines
! ************************************************************
! ************************************************************
! ************************************************************

	subroutine colcopy(nval,fisov,col)

!  sets colors and isolevels

	implicit none

	integer nval
	real fisov(1), col(1)
	integer i

	include 'color.h'

	call colini

	if( nval .gt. isopar .or. nval .lt. 0 ) then
		write(6,*) 'nval,isopar: ',nval,isopar
		stop 'error stop colcopy: error in parameters'
	end if

	isoanz = nval

	do i=1,nval
	  fiso(i) = fisov(i)
	  ciso(i) = col(i)
	end do

	ciso(nval+1) = col(nval+1)

	end
 
! ************************************************************

	subroutine colnul(fnul)

!  sets null value

	implicit none

	real fnul

	include 'color.h'

	call colini

	fnull = fnul

	end
 
! ************************************************************

	subroutine colinfo(numiso,fnul)

!  returns info on color table

	implicit none

	integer numiso
	real fnul

	include 'color.h'

	call colini

	numiso = isoanz
	fnul = fnull

	end
 
! ************************************************************

	subroutine colentry(icol,viso,color)

!  returns entry of color table

	implicit none

	integer icol
	real viso,color

	include 'color.h'

	call colini

	viso = 0.
	color = 0.

	if( icol .le. 0 ) then
	  return
	else if( icol .le. isoanz ) then
	  viso = fiso(icol)
	  color = ciso(icol)
	else if( icol .eq. isoanz+1 ) then
	  color = ciso(icol)
	else
	  return
	end if

	end
 
! ************************************************************

	function getcol(value)

!  gets color for value

	implicit none

	real getcol
	real value
	integer i

	include 'color.h'

	do i=1,isoanz
	  getcol = ciso(i)
	  if( fiso(i) .gt. value ) return
	end do

	getcol = ciso(isoanz+1)

	end

! ************************************************************

	function getcolval(rindex)

!  gets value for real index
! 
!  if array has been read -> closest value
!  if regular values (valmin/max) -> use real rindex

	implicit none

	real getcolval
	real rindex

	real ri,val,dval

	include 'color.h'

	ri = rindex
	ri = min(ri,float(isoanz))
	ri = max(ri,1.)

	if( nisord .gt. 0 .or. isoanz .eq. 1 ) then	!values read
	  val = fiso(nint(ri))
	else
	  dval = (fiso(isoanz) - fiso(1)) / (isoanz-1)
	  val = fiso(1) + (ri-1.) * dval
	end if

	getcolval = val

	end
 
! ************************************************************
! ************************************************************
! ************************************************************

	subroutine setcol(nkn,parray,dp)

!  sets new color table -> old , do not use -> subst by colsetup, colauto
! 
!  uses pmin/pmax for extension and ncol colors

	implicit none

	integer ndim
	real eps
	parameter(ndim=100,eps=1.e-4)

	integer nkn
	real dp
	real parray(1)

	logical bdebug
	integer ncol,i
	real dc
	real pmin,pmax
	real colmin,colmax
	real col(ndim),riso(ndim-1)

	real rnexta,rnext
	real getpar

	bdebug = .false.

	call colini

	if( dp .eq. 0. ) return

! 	colmin = getpar('colmin')
! 	colmax = getpar('colmax')
	colmin = 0.3
	colmax = 0.9

	call mima(parray,nkn,pmin,pmax)

	if( bdebug ) then
	  write(6,*) 'new color:'
	  write(6,*) pmin,pmax
	end if

	if( dp .lt. 0. ) then	!dp is aprox. number of isolines
	  dp = - (pmax-pmin) / dp
	  dp = rnext(dp,+3)
	end if

	pmin = rnexta(pmin,-3)
	pmax = rnexta(pmax,+3)

	ncol = (pmax-pmin) / (dp-eps)
	if( ncol .le. 2 ) then
	  dc = 0.
	else 
	  dc = (colmax-colmin) / (ncol-2)
	end if

	if( ncol .gt. ndim ) stop 'error stop setcol: ndim'

	do i = 1,ncol-1
	  col(i) = colmin + (i-1)*dc
	  riso(i) = pmin + i*dp
	end do
	if( ncol .le. 0 ) ncol = 1
	col(ncol) = -1.

	call colcopy(ncol-1,riso,col)

	if( bdebug ) then
	  write(6,*) pmin,pmax
	  write(6,*) ncol,dc,dp
	  write(6,*) (col(i),i=1,ncol)
	  write(6,*) (riso(i),i=1,ncol-1)
	end if
	  
	end
	
! *****************************************************************
! *****************************************************************
! *****************************************************************

	subroutine qchsc( bcolor )

	implicit none

	logical bcolor

	logical bcol
	common /colgry/ bcol
	save /colgry/

	bcol = bcolor

	end

! *****************************************************************

	subroutine qsetc( color )

	implicit none

	real color

	logical bcol
	common /colgry/ bcol
	save /colgry/

	if( bcol ) then
	  call qhue(color)
	else
	  call qgray(color)
	end if

	end

! *****************************************************************
! *****************************************************************
! *****************************************************************

