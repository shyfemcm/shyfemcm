
!--------------------------------------------------------------------------
!
!    Copyright (C) 2001,2008-2010,2012-2014,2016-2019  Georg Umgiesser
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
! 
!  revision log :
! 
!  26.10.2001	ggu	in colrd allow only one array to be given
!  06.06.2008	ggu	new colortable introduced
!  19.11.2008	ggu	new routine coltst()
!  09.01.2009	ggu	some re-formatting
!  09.01.2009	ggu	deleted setcol(), new make_single_isolines()
!  26.01.2009	ggu	some more color tables
!  23.02.2010	ggu	new set_default_color_table(), renamed qchsc()
!  23.03.2010	ggu	make it easier to change between color tables
!  29.09.2010	ggu	changed VERS_6_1_12
!  01.06.2012	ggu	changed VERS_6_1_53
!  12.09.2013	ggu	changed VERS_6_1_67
!  18.07.2014	ggu	changed VERS_7_0_1
!  05.12.2014	ggu	changed VERS_7_0_8
!  17.06.2016	ggu	changed VERS_7_5_15
!  30.09.2016	ggu	changed VERS_7_5_18
!  11.07.2017	ggu	changed VERS_7_5_30
!  07.12.2017	ggu	changed VERS_7_5_40
!  24.01.2018	ggu	changed VERS_7_5_41
!  19.04.2018	ggu	changed VERS_7_5_45
!  25.10.2018	ggu	changed VERS_7_5_51
!  18.12.2018	ggu	changed VERS_7_5_52
!  18.09.2024	ggu	new g/set_default_color(), icolor renamed to icoltab
! 
! ************************************************************

	subroutine colini

!  initializes color

	use color

	implicit none

	logical, save :: binit = .false.

	if( binit ) return

	binit = .true.

	!write(6,*) 'initializing color with colini...'

	isopar = isodim
	isoanz = 0
	nisord = 0
	ncolrd = 0
	icauto = -1
	fnull = -999.0
	ciso(1) = -1.

	end

! ************************************************************

	subroutine colrd

!  reads color block from str file
! 
!  this is not called anymore - the setup is done in colchk

	use color

	implicit none

	character*80 name,text
	logical bdebug
	integer mode
	integer i
	integer ndim
	real value
	double precision dvalue

	integer nrdpar
	real getpar

	call colini

	ndim = isodim
	bdebug = .true.
	bdebug = .false.
	ncolrd = 0
	nisord = 0

        mode = 1
        do while( mode .ne. 0 )
            mode = nrdpar('color',name,dvalue,text)
	    value = dvalue
	    !write(6,*) mode,name,value,text
            if( mode .eq. 2 ) then
	      if( name .ne. 'color' .and. name .ne. 'isoval' ) then
		write(6,*) 'Variable is no vector : ',name
		stop 'error stop : colrd'
	      end if
	    end if
            if( name .eq. 'color' ) then
                ncolrd = ncolrd + 1
		if( ncolrd .gt. ndim+1 ) stop 'error stop colrd: dimension'
                ciso(ncolrd) = value
            else if( name .eq. 'isoval' ) then
                nisord = nisord + 1
		if( nisord .gt. ndim ) stop 'error stop colrd: dimension'
                fiso(nisord) = value
	    end if
        end do

	if( ncolrd .gt. 0 .and. nisord .gt. 0 ) then
	  if( nisord + 1 .ne. ncolrd ) then
		write(6,*) 'ncolrd,nisord: ',ncolrd,nisord
		write(6,*) 'There must be one more color than isovalue'
		stop 'error stop : colrd'
	  end if
	end if

	isoanz = nisord
! 	call putpar('dval',0.)	!if section is given -> no automatic

	if( bdebug ) then
	  write(6,*) 'color debug : ',nisord,ncolrd
	  write(6,*) (fiso(i),i=1,nisord)
	  write(6,*) (ciso(i),i=1,ncolrd)
	end if

	end

! ************************************************************

	subroutine colchk

	use color
	use para

	implicit none

	logical bdebug
	integer i,ndim

	bdebug = .true.
	bdebug = .false.

	call colini

	ndim = isodim

	if( para_has_array('color') ) then
	  call para_get_array_size('color',ncolrd)
	end if
	if( para_has_array('isoval') ) then
	  call para_get_array_size('isoval',nisord)
	end if

	if( ncolrd .gt. ndim+1 ) stop 'error stop colchk: dimension'
	if( nisord .gt. ndim ) stop 'error stop colchk: dimension'

	!call para_info_name('color')
	!call para_info_name('isoval')

	if( ncolrd .gt. 0 .and. nisord .gt. 0 ) then
	  if( nisord + 1 .ne. ncolrd ) then
		write(6,*) 'ncolrd,nisord: ',ncolrd,nisord
		write(6,*) 'There must be one more color than isovalue'
		stop 'error stop : colchk'
	  end if
	end if

	if( ncolrd > 0 ) then
	  call para_get_array_value('color',ndim+1,ncolrd,ciso)
	end if
	if( nisord > 0 ) then
	  call para_get_array_value('isoval',ndim,nisord,fiso)
	end if

	isoanz = nisord

	if( bdebug ) then
	  write(6,*) 'color debug : ',nisord,ncolrd
	  write(6,*) (fiso(i),i=1,nisord)
	  write(6,*) (ciso(i),i=1,ncolrd)
	end if

	end

! ************************************************************

	subroutine coltst

!  debug write of color

	use color

	implicit none

	integer i

	write(6,*) 'coltst: debug write of color global params...'

	write(6,*) 'isodim : ',isodim
	write(6,*) 'isopar : ',isopar
	write(6,*) 'isoanz : ',isoanz
	write(6,*) 'nisord : ',nisord
	write(6,*) 'ncolrd : ',ncolrd
	write(6,*) 'icauto : ',icauto
	write(6,*) 'fnull : ',fnull

	if( nisord .lt. 0 .or. nisord .gt. 1000 ) goto 99
	if( ncolrd .lt. 0 .or. ncolrd .gt. 1000 ) goto 99

	write(6,*) 'fiso : ',(fiso(i),i=1,nisord)
	write(6,*) 'ciso : ',(ciso(i),i=1,ncolrd)

	return
   99	continue
	write(6,*) 'unrealistic values for nisord or ncolrd'
	stop 'error stop coltst: nisord/ncolrd'
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

	use color

	implicit none

	integer nval
	real fisov(nval), col(nval+1)
	integer i

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

	use color

	implicit none

	real fnul

	call colini

	fnull = fnul

	end
 
! ************************************************************

	subroutine colinfo(numiso,fnul)

!  returns info on color table

	use color

	implicit none

	integer numiso
	real fnul

	call colini

	numiso = isoanz
	fnul = fnull

	end
 
! ************************************************************

	subroutine colminmax(vmin,vmax)

!  returns minimum and maximum value in color table

	use color

	implicit none

	real vmin,vmax

	call colini

	vmin = fiso(1)
	vmax = fiso(isoanz)

	end

! ************************************************************

	subroutine colentry(icol,viso,col)

!  returns entry of color table

	use color

	implicit none

	integer icol
	real viso,col

	call colini

	viso = 0.
	col = 0.

	if( icol .le. 0 ) then
	  return
	else if( icol .le. isoanz ) then
	  viso = fiso(icol)
	  col = ciso(icol)
	else if( icol .eq. isoanz+1 ) then
	  col = ciso(icol)
	else
	  return
	end if

	end
 
! ************************************************************

	function getcol(value)

!  gets color for value

	use color

	implicit none

	real getcol
	real value
	integer i

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

	use color

	implicit none

	real getcolval
	real rindex

	real ri,val,dval

	ri = rindex
	ri = min(ri,float(isoanz))
	ri = max(ri,1.)

	!write(6,*) 'ri: ',rindex,ri,nisord,isoanz,iusear
	!write(6,*) 'fiso: ',fiso(1:isoanz)

	if( nisord .gt. 0 .or. isoanz .eq. 1 ) then	!values read
	  val = fiso(nint(ri))
	else if( iusear .ne. 0 ) then			!use array
	  val = fiso(nint(ri))
	else
	  dval = (fiso(isoanz) - fiso(1)) / (isoanz-1)
	  val = fiso(1) + (ri-1.) * dval
	end if

	getcolval = val

	end
 
! ************************************************************

	subroutine make_single_isolines(ntick)

!  computes values for which single isolines are plotted

	use color

	implicit none

	integer ntick

	integer i
	real dtick,rit,value
	real getcolval

	if( ntick .le. 1 ) then
	  nriso = isoanz
	  do i=1,nriso
	    riso(i) = fiso(i)
	  end do
	else
	  nriso = ntick
	  dtick = float(isoanz-1) / (nriso-1)
	  do i=1,nriso
            rit = 1. + (i-1)*dtick
            value = getcolval(rit)
	    riso(i) = value
	  end do
	end if

	!write(6,*) 'make_single_isolines: ',nriso
	!write(6,*) (riso(i),i=1,nriso)

	end

! ************************************************************
! ************************************************************
! ************************************************************

	subroutine set_default_color( rcolor )

!  sets default color of the actual color table

	use color

	implicit none

	real rcolor

	color_def = rcolor
	call qsetc(rcolor)

	end

! ************************************************************

	subroutine get_default_color( rcolor )

!  gets default color of the actual color table

	use color

	implicit none

	real rcolor

	rcolor = color_def

	end

! ************************************************************
! ************************************************************
! ************************************************************

	subroutine set_default_color_table( icoltab )

!  sets default color table

	implicit none

	integer icoltab

	call qtdef(icoltab)

	end

! *****************************************************************

	subroutine get_color_table( icoltab )

!  gets color table used

	use color

	implicit none

	integer icoltab

	icoltab = icoltab_def

	end

! *****************************************************************

	subroutine set_color_table( icoltab )

!  sets color table to be used

	use color

	implicit none

	integer icoltab

	icoltab_def = icoltab
	call qsetc(0.) 		!try if color table exists

	end

! *****************************************************************

	subroutine reset_color_table

!  sets default color table

	use color

	implicit none

	call set_color_table( icoltab_def )

	end

! *****************************************************************

	subroutine write_color_table

!  writes default color table to terminal

	use color

	implicit none

	write(6,*) 'color table: actual = ',icoltab_def
	!write(666,*) 'color table: actual = ',icoltab_def

	end

! *****************************************************************
! 
!  set image colorscale hsb 0.666 0.0 1.0 .min. hsb 0.666 1.0 1.0 .max.
!  Panel "3. White-blue (HSB blending)"
!  
!  set image colorscale rgb 1.0 0.0 0.0 .min. rgb 0.0 0.0 1.0 .max.
!  Panel "4a. Red-blue (RGB blending)"
! 
! *****************************************************************

	subroutine qsetc( col )

!  changes color using the actual color table

	use color

	implicit none

	real col

	integer icoltab

	icoltab = icoltab_def

	if( icoltab .eq. -1 ) then
	  call qcolor(col)
	else if( icoltab .eq. 0 ) then
	  call qgray(col)
	else if( icoltab .eq. 1 ) then
	  call qhue(col)
	else if( icoltab .eq. 2 ) then
	  call white_blue( col )
	else if( icoltab .eq. 3 ) then
	  call white_red( col )
	else if( icoltab .eq. 4 ) then
	  call blue_white_red( col )
	else if( icoltab .eq. 5 ) then
	  call blue_black_red( col )
	else if( icoltab .eq. 6 ) then
	  call hue_sqrt( col )
	else if( icoltab .eq. 7 ) then
	  call hue_pow2( col )
	else if( icoltab .eq. 8 ) then
	  call col_cust( col )
	else
	  write(6,*) 'icoltab = ',icoltab
	  stop 'error stop qsetc: no such color table'
	end if

	end

! *****************************************************************
! *****************************************************************
! *****************************************************************

	subroutine white_blue( color )
	implicit none
	real color
	call qhsb(0.666,color,1.)
	end

	subroutine white_red( color )
	implicit none
	real color
	call qhsb(1.,color,1.)
	end

	subroutine blue_white_red( color )
	implicit none
	real color
	if( color .le. 0.5 ) then
	  call qhsb(0.666,1.-2.*color,1.)
	else
	  call qhsb(1.,2.*(color-0.5),1.)
	end if
	end

	subroutine blue_black_red( color )
	implicit none
	real color
	if( color .le. 0.5 ) then
	  call qhsb(0.666,1.,1.-2.*color)
	else
	  call qhsb(1.,1.,2.*(color-0.5))
	end if
	end

	subroutine hue_sqrt( color )
	implicit none
	real color
	call qhue(sqrt(color))
	end

	subroutine hue_pow2( color )
	implicit none
	real color
	call qhue(color*color)
	end

	subroutine col_cust( col )
	use color
	implicit none
	real col
	real c
	integer ic
	c = col
	c = min(1.,c)
	c = max(0.,c)
	ic = nint((icmax-1)*c) + 1
	call qrgb(coltab(1,ic),coltab(2,ic),coltab(3,ic))
	!write(6,*) 'color... : ',color,c,ic,coltab(:,ic)
	end

! *****************************************************************

