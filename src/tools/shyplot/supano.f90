
!--------------------------------------------------------------------------
!
!    Copyright (C) 2000-2001,2003-2019  Georg Umgiesser
!    Copyright (C) 2016  Christian Ferrarin
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

!  routines for anotating plot
! 
!  contents :
! 
!  revision log :
! 
!  09.02.2000	ggu	use inboxdim to compute box to plot
!  07.04.2000	ggu	new format for legend read
!  12.06.2000	ggu	ndim renamed into legdim
!  04.02.2001	ggu	color bar legend now works with nctick
!  21.05.2003	ggu	various changes in legend
!  21.08.2003	ggu	uses occupy to decide where to put color bar
!  04.09.2003	ggu	new routines to read legend, legend.h
!  05.10.2003	ggu	new routines to plot vectors autonomously
!  08.03.2004	ggu	legdate to write date to plot (called in legplo)
!  03.09.2004	ggu	legwind to write wind data (preliminary)
!  22.09.2004	ggu	small bug fix in legwind, qlwidth(0.) for reset
!  05.10.2004	ggu	changes in legdate and legwind
!  02.12.2004	ggu	legwind and legdate prepared for str input
!  02.03.2005	ggu	new format for exffil
!  11.03.2005	ggu	date and wind are now specified in STR file
!  27.05.2005	ggu	new routines to make absolute coordinates (for text)
!  18.10.2006	ggu	may give relative coords for date/time
!  06.06.2007	ggu	new mode for idate (4)
!  06.12.2008	ggu	new routine blank_window(), in velsh handle wind/waves
!  06.12.2008	ggu	make colbar smaller than box size
!  09.01.2009	ggu	allow any factor in velsh, no traref etc..
!  21.04.2009	ggu	allow just one isoline (bug)
!  09.10.2009	ggu	sclvel for scaling arrow, lots of changes in color bar
!  13.10.2009	ggu	changes to colbar, new inbox routines
!  23.02.2010	ggu	for colorbar switch to generic color table
!  23.03.2010	ggu	changed v6.1.1
!  09.12.2011	ggu	changed VERS_6_1_38
!  30.03.2012	ggu	changed VERS_6_1_51
!  01.06.2012	ggu	new circle item for legend, auto determ for ndec
!  12.09.2012	ggu	changed VERS_6_1_57
!  19.06.2013	ggu	changed VERS_6_1_66
!  05.09.2013	ggu	adjust for scale distortion in ref and wind arrows
!  12.09.2013	ggu	changed VERS_6_1_67
!  05.03.2014	ggu	in annotation use date information
!  30.05.2014	ggu	changed VERS_6_1_76
!  18.07.2014	ggu	changed VERS_7_0_1
!  16.10.2014	ggu	annotes doesnt need it anymore
!  26.11.2014	ggu	changed VERS_7_0_7
!  05.12.2014	ggu	changed VERS_7_0_8
!  23.12.2014	ggu	changed VERS_7_0_11
!  19.01.2015	ggu	changed VERS_7_1_3
!  26.02.2015	ggu	changed VERS_7_1_5
!  01.04.2015	ggu	changed VERS_7_1_7
!  05.05.2015	ggu	changed VERS_7_1_10
!  06.05.2015	ggu	prepared for logarithmic color bar
!  21.05.2015	ggu	changed VERS_7_1_11
!  05.06.2015	ggu	new keyword vart in legend for variable text
!  10.07.2015	ggu	changed VERS_7_1_50
!  17.07.2015	ggu	changed VERS_7_1_80
!  20.07.2015	ggu	changed VERS_7_1_81
!  19.02.2016	ggu	changed VERS_7_5_2
!  17.06.2016	ccf	include kn units and correct date/time for tzshow
!  13.09.2016	ggu	old legend deleted, new keyword ctxt for centering
!  30.09.2016	ggu	changed VERS_7_5_18
!  11.07.2017	ggu	changed VERS_7_5_30
!  14.11.2017	ggu	changed VERS_7_5_36
!  22.02.2018	ggu	changed VERS_7_5_42
!  06.07.2018	ggu	changed VERS_7_5_48
!  18.12.2018	ggu	changed VERS_7_5_52
!  21.05.2019	ggu	changed VERS_7_5_62
!  18.09.2024	ggu	new parameter rfaccol, new log colorbar
!  09.01.2025	ggu	avoid divide by zero in scale_legend(): 10 -> 10.
! 
!  notes :
! 
!  for format of legend please see routine newleg()
! 
! ***************************************************************

	module legend

	implicit none

        integer legdim			!maximum number of legend entries
        parameter(legdim=400)

        integer, save :: nleg,nlegdi,iplotleg

        real, save :: xleg(2,legdim), yleg(2,legdim)

        character*80, save :: legleg(legdim)

        character*4, save :: whatleg(legdim)

        integer, save :: legsiz(legdim)

        real, save :: aleg(legdim), cleg(legdim)

	end module legend

! ***************************************************************

	subroutine annotes(var)

!  writes annotation for simulation

	use simul

	implicit none

	character*(*) var

	character*80 line
        logical debug
	integer it
        integer iaux
	integer date,time
	real x,y,dx,wx
	real xmin,ymin,xmax,ymax
        real aux
	double precision atime

	integer getlev,ialfa
	logical dts_has_date

        debug = .true.
        debug = .false.

! --------------------------------------------------------------
!  set up variables
! --------------------------------------------------------------

	call qgetvp(xmin,ymin,xmax,ymax)

	call qtxts(6)
        call qfont('Courier')

! --------------------------------------------------------------
!  determine total width to use
! --------------------------------------------------------------

	dx = (xmax - xmin)
	wx = 0.9 * dx

        if( debug ) then
	  write(6,*) 'annoting...'
	  write(6,*) var
	  write(6,*) it,getlev()
	  write(6,*) xmin,ymin,xmax,ymax
        end if

! --------------------------------------------------------------
!  write description of simulation
! --------------------------------------------------------------

	y = ymin + 0.6
	x = xmin + 0.1 * dx
	call qtext(x,y,descrp)

! --------------------------------------------------------------
!  write time in seconds
! --------------------------------------------------------------

	y = ymin + 1.1
	x = xmin + 0.1 * dx

	call ptime_get_itime(it)
	write(line,'(i12,a)') it,' sec'
	call qtext(x,y,line)

! --------------------------------------------------------------
!  write time/date
! --------------------------------------------------------------

	x = x + wx/4.
	if( dts_has_date() ) then
	  call ptime_get_atime(atime)
	  call dts_format_abs_time(atime,line)
	  call dts_get_date(date,time)
	  !write(6,*) 'time: ',date,time,atime,line(1:20)
	else
	  call makehm(it,line)    !makes time and date
	end if

	call qtext(x,y,line)

! --------------------------------------------------------------
!  write level
! --------------------------------------------------------------

	x = x + wx/4.
        line = 'level = '
        aux = getlev()
        iaux = ialfa(aux,line(9:),-1,-1)
	call qtext(x,y,line)

! --------------------------------------------------------------
!  write description of what is plotted
! --------------------------------------------------------------

	x = x + wx/4.
	call qtext(x,y,var)

! --------------------------------------------------------------
!  reset plotting area -> is done in annote (basin)
! --------------------------------------------------------------

! 	call qsetvp(xmin,ymin+1.5,xmax,ymax)

! --------------------------------------------------------------
!  end of routine
! --------------------------------------------------------------

	end

! ***************************************************************

	subroutine annote

!  writes annotation

	implicit none

	character*80 line
	real x,y,dx,wx
	real xmin,ymin,xmax,ymax

	call qgetvp(xmin,ymin,xmax,ymax)

	call qtxts(6)
        call qfont('Courier')

	y = ymin + 0.1
	dx = (xmax - xmin)
	wx = 0.9 * dx

	x = 0.1 * dx
	call getfnm('basnam',line)
	call qtext(x,y,line)

	x = x + wx/3.
	call getfnm('runnam',line)
	call qtext(x,y,line)

	x = x + wx/3.
	call getfnm('apnnam',line)
	call qtext(x,y,line)

	call qsetvp(xmin,ymin+1.5,xmax,ymax)

	end

! ****************************************************************

	subroutine makehm(it,line)

!  converts it (seconds) to ddd:hh:mm (alpha string)

	implicit none

!  arguments
	integer it
	character*(*) line
!  local
	integer is,ih,im,id,ii
	character sec*2, min*2, hour*2, day*5, sign*1
!  functions
	integer ialfa

	is = it
	sign = ' '
	if( is .lt. 0 ) then
		is = -is
		sign = '-'
	end if

	id = is/86400
	ih = (is-id*86400)/3600
        im = (is-id*86400-ih*3600)/60
	is = (is-id*86400-ih*3600-im*60)

	
        ii = ialfa(float(is),sec,-1,+1)
	if( ii .eq. 1 ) sec(1:1) = '0'

        ii = ialfa(float(im),min,-1,+1)
	if( ii .eq. 1 ) min(1:1) = '0'

        ii = ialfa(float(ih),hour,-1,+1)
	if( ii .eq. 1 ) hour(1:1) = '0'

        ii = ialfa(float(id),day,-1,-1)

	line = sign//day(1:ii)//' '//hour//':'//min//':'//sec

	return
	end

! *******************************************************************

	subroutine scalbar(x,y,dist,h,ntack,scale,ndec,unit)

!  plots scale

	implicit none

	real x,y,dist,h
	real scale
	integer ndec
	integer ntack
	character*(*) unit

	integer ii,n
	real dh,xa,ya,xe,ye
	real z
	real width,height
	integer ialfa
	character*10 line

	dh = 0.5 * h
	call qgray(0.)

!  outer box

	xe = x + dist*ntack
	ye = y + h

	call qline(x,y,xe,y)
	call qline(xe,y,xe,ye)
	call qline(xe,ye,x,ye)
	call qline(x,ye,x,y)

!  internal division

	ye = y + dh
	call qline(x,ye,xe,ye)

	ye = y + h
	do n = 1,ntack-1
	  xe = x + dist*n
	  call qline(xe,y,xe,ye)
	end do

!  black boxes

	ii = 0
	do n=1,ntack
	  ii = mod(ii+1,2)
	  xa = x + (n-1)*dist
	  xe = x + n*dist
	  ya = y + ii * dh
	  ye = ya + dh
	  call qrfill(xa,ya,xe,ye)
	end do

!  legend

        call qtxts(12)
        call qfont('Times-Roman')
	call qtsize(' ',width,height)

	ye = y - 2. * h
	ye = y - 1.5 * h		!changed 27.11.1997
	ye = y - 1.5 * height		!changed 02.12.1997
	do n=0,ntack
	  xe = x + n*dist
	  z = n * dist * scale
	  ii = ialfa(z,line,ndec,-1)	!ialfa(z,line,ndec,mode)
	  call qtext(xe,ye,line(1:ii))
	end do

	call qtext(999.,999.,unit)
	
	end

! 
! *************************************************************
! 
        subroutine dnord0(xm,ym,radius,angle)
! 
!  plots northern direction
! 
!  xm,ym         coordinates of centre of circle
!  radius        radius of circle
!  angle         angle of inclination of true north
!                ...positive if true north is pointing to the east
! 
	implicit none

	real xm,ym,radius,angle

	real radian
        parameter (radian=3.14159/180.)

	integer i
	real hsymb
	real ri,x,y
	real xtop,ytop,angbot
	real xbot1,ybot1,xbot2,ybot2
	real xh,yh,xr,yr
! 
        hsymb=radius/4.
! 
!  plot circle
! 
        call qmove(xm+radius,ym)
! 
        do i=1,360
        ri=i*radian
        x=xm+radius*cos(ri)
        y=ym+radius*sin(ri)
        call qplot(x,y)
        end do
! 
!  plot arrow
! 
        xtop=xm+radius*sin(angle*radian)
        ytop=ym+radius*cos(angle*radian)
! 
        angbot=amod(angle+180.,360.)
! 
        xbot1=xm+radius*sin((angbot+10.)*radian)
        ybot1=ym+radius*cos((angbot+10.)*radian)
! 
        xbot2=xm+radius*sin((angbot-10.)*radian)
        ybot2=ym+radius*cos((angbot-10.)*radian)
! 
        call qmove(xtop,ytop)
        call qplot(xbot1,ybot1)
        call qplot(xbot2,ybot2)
        call qplot(xtop,ytop)
! 
!  find lower left corner for 'N' and plot it
! 
        xh=-hsymb*2./6.
        yh=radius+0.5*hsymb
! 
        xr= xh*cos(angle*radian)+yh*sin(angle*radian)
        yr=-xh*sin(angle*radian)+yh*cos(angle*radian)
! 
!         call qsymb(xr+xm,yr+ym,hsymb,'N',-angle)
! 
	return
	end
! 
! *************************************************************
! 
        subroutine dnord(xm,ym,radius,angle)
! 
!  plots northern direction
! 
!  xm,ym         coordinates of centre of circle
!  radius        radius of circle
!  angle         angle of inclination of true north
!                ...positive if true north is pointing to the east

	implicit none

	real radian
        parameter (radian=3.14159/180.)

	real xm,ym,radius,angle
	real r,r2,r4
	real xp(4),yp(4)

	r = radius
	r2 = r / 2.
	r4 = r / 4.

	call qgray(0.)

	call qmove(xm,ym+r)
	call qplot(xm+r4,ym+r4)
	call qplot(xm+r,ym)
	call qplot(xm+r4,ym-r4)
	call qplot(xm,ym-r)
	call qplot(xm-r4,ym-r4)
	call qplot(xm-r,ym)
	call qplot(xm-r4,ym+r4)
	call qplot(xm,ym+r)

	call qline(xm-r,ym,xm+r,ym)
	call qline(xm,ym-r,xm,ym+r)

	call qmove(xm+r2,ym+r2)
	call qplot(xm+r4,ym)
	call qplot(xm+r2,ym-r2)
	call qplot(xm,ym-r4)
	call qplot(xm-r2,ym-r2)
	call qplot(xm-r4,ym)
	call qplot(xm-r2,ym+r2)
	call qplot(xm,ym+r4)
	call qplot(xm+r2,ym+r2)

	call qline(xm+r2,ym+r2,xm-r2,ym-r2)
	call qline(xm+r2,ym-r2,xm-r2,ym+r2)

	xp(1) = xm+0.
	yp(1) = ym+0.
	xp(2) = xm-r4
	yp(2) = ym+r4
	xp(3) = xm+0.
	yp(3) = ym+r
	xp(4) = xm+r4
	yp(4) = ym+r4
	call qafill(4,xp,yp)

! 	call qfont('Times-Roman')
! 	call qtxts()

        end

! *****************************************************************

	subroutine scale_legend(x0,y0,x1,y1)

!  creates legend (north wheel and scale bar)

	implicit none

	real x0,y0,x1,y1	!dimension of legend on plot

	character*3 unit
	real dx,dy
	real xr,yr,r		!for north wheel
	real x0s,y0s		!starting point of scale bar
	real h			!height of scale bar
	real dd			!length of scale bar
	real lblank		!factor for blanking
	real scale
	integer is,i0
	integer ndec,ntics

	integer iround,istell
	real rnext
	real getpar

	logical bblank

	bblank = .true.		!blanking for scale legend
	lblank = getpar('lblank')
	bblank = lblank > 0

	!write(6,*) '===================== legend =============='

!  blanking of window

	if( bblank ) call blank_window(lblank,x0,y0,x1,y1)

!  extension of window

	dx = x1 - x0
	dy = y1 - y0

!  dimensions for north wheel

	xr = x0 + 0.5 * dx
	yr = y0 + 0.75 * dy
	r = 0.2 * dy

!  find length of scale bar - dd must not be larger than 80% of full length

	dd = rnext(dx,-3)			!find lower 1,2,3,4,5,8
	if( dd .gt. 0.8*dx ) then
	  dd = rnext(0.8*dx,-3)
	end if

!  dimensions for scale bar (starting point and height)

	x0s = x0 + 0.5 * ( dx - dd )
	y0s = y0 + dy * 0.35
	h = 0.1 * dy

! 	write(6,*) x0,y0,dx,dy
! 	write(6,*) x0s,y0s,h

!  find leading digit (must be 1,2,3,4,5,8)

	is = istell(dd)
	i0 = iround( dd / 10.0**is )

	if( i0 .eq. 2 .or. i0 .eq. 4 .or. i0 .eq. 8 ) then
		ntics = 4
	else if( i0 .eq. 3 ) then
		ntics = 3
	else
		ntics = 5
	end if

	if( dd .lt. 1000 ) then	!less than 1 km
		scale = 1.
		unit = " m "
	else
		scale = 1.e-3
		unit = " km"
	end if

	if( dd .gt. 850. .and. dd .lt. 2500. ) then
		ndec = 1
	else 
		ndec = -1
	end if

	!write(6,*) 'Plotting legend: ',ntics,ndec
	!write(6,*) x0,y0,x1,y1
	!write(6,*) x0s,y0s,dd/ntics,h,ntics,scale,ndec,unit
	!write(6,*) xr,yr,r

        call qcomm('Plotting legend')
        call qtxts(12)
        call qfont('Times-Roman')

        call qcomm('Plotting scale bar')
	call scalbar(x0s,y0s,dd/ntics,h,ntics,scale,ndec,unit)
        call qcomm('Plotting north')
        call dnord(xr,yr,r,0.)

! 	call qmove(x0,y0)
! 	call qplot(x1,y0)
! 	call qplot(x1,y1)
! 	call qplot(x0,y1)
! 	call qplot(x0,y0)

        call qcomm('end of plotting legend')

	end

! ******************************************************************

	subroutine colsh

!  shell for color bar

	use color

	implicit none

	character*80 line
	integer ndec
	integer nctick,ipllog
	logical bhoriz
	real fact,rfact
	real x0,y0,x1,y1

	integer iround
	real getpar
	logical inboxdim

	if( .not. inboxdim('col',x0,y0,x1,y1) ) then

          call occupy_dim(x0,y0,x1,y1)
	  call make_absolute(x0,y0,x1,y1)

	end if

	bhoriz = .true.
	fact = getpar('faccol')
        rfact = getpar('rfaccol')
        if( fact == 1. .and. rfact /= 1. ) fact = 1. / rfact
	ndec = iround(getpar('ndccol'))
	nctick = iround(getpar('nctick'))
	ipllog = iround(getpar('ipllog'))
	call getfnm('legcol',line)

	!write(6,*) 'colsh : ',nctick,ndec,fact,x0,y0,x1,y1
	call colbar(line,bhoriz,nctick,ipllog,ndec,fact,x0,y0,x1,y1)

	end

! ******************************************************************

	subroutine colbar(title,bhoriz,nticks,ipllog,ndec,fact &
     &				,xmin,ymin,xmax,ymax)

!  writes color bar

	use plot_fonts

	implicit none

	character*(*) title
	logical bhoriz		!horizontal?
	integer nticks		!number of ticks to comment (0 -> all colors)
	integer ipllog
	integer ndec
	real fact
	real xmin,ymin,xmax,ymax

!  color
	integer isoanz
	real fnull

	integer ndim
	parameter (ndim=100)
	real aval(ndim)
	real rval(ndim)

	character*80 line
	integer i,imax,nw
	integer ntk,ntks
	integer nmin,nmax
	integer icsave
	integer idec
	logical bblank,brotate
	logical blowtck
	real dx,dy
	real x,y
	real x0,x1,x2,y0,y1,y2
	real dxc
	real value
	real col,xlow,xhigh
	real width, height
	real dtick,dt,rit
	real xperc,yperc
	real ylow,yhigh,dyc,dxy
	real xcm,ycm,dxtick,dytick
	real amin,amax
	real cblank			!factor for blanking
	integer idiv,ntkl
	integer it

	logical bdebug,blog
	integer iround,isoget,ialfa,ideci
	real getcolval
	real getpar

	call qcomm('Plotting colorbar')
        call qfont('Times-Roman')
        !call qtxts(12)
	call qtxts(fs_cbar_num)
		
	bdebug = .true.
	bdebug = .false.
	bblank = .true.		!blanking for color bar
	cblank = getpar('cblank')
	bblank = cblank > 0
	!bhoriz = .false.
	!bhoriz = .true.
	brotate = .true.	!rotate numbers for vertical color bar

	if( bdebug ) then
	  write(6,*) 'colbar: ',bhoriz,brotate
	  write(6,*) 'title: ',trim(title)
	  write(6,*) nticks,ipllog,ndec
	  write(6,*) fact
	  write(6,*) xmin,ymin,xmax,ymax
	  call write_color_table
	end if

!  blanking of window

	if( bblank ) call blank_window(cblank,xmin,ymin,xmax,ymax)

!  determine min/max

	dx = xmax - xmin	! total width of colorbar
	dy = ymax - ymin	! height of color bar (including legend)

	call colinfo(isoanz,fnull)
	imax = isoanz + 1	! total number of colors

	if( bdebug ) then
	  write(6,*) 'colbar:  colors = ',imax,'  isolines = ',isoanz
	end if

!  position of colorbar on plot (colorbar is a little smaller than box)
!  (x0,y0) and (x1,y1) are the extremes of the color bar
!  dx,dy is width of box

	if( bhoriz ) then
	  xperc = 0.10
	  yperc = 0.30
	else
	  xperc = 0.30
	  yperc = 0.10
	end if

        x0 = xmin + xperc * dx
        x1 = xmin + (1.-xperc) * dx
	y0 = ymin + yperc * dy
	y1 = ymin + (1.-yperc) * dy

!  plot single boxes - dxc/dyc are width/height of single color box

	if( bdebug ) call write_color_table
	!call set_auto_color_table
	!if( bdebug ) call write_color_table
	if( bdebug ) then
	  write(6,*) 'isoanz,imax: ',isoanz,imax
	end if

	xlow = x0
	ylow = y0
	dxc = (x1-x0)
	dyc = (y1-y0)
	if( bhoriz ) then
	  dxc = dxc / imax
	else
	  dyc = dyc / imax
	end if

	do i=1,imax
	  call colentry(i,value,col)
	  if( bhoriz ) then
	    xlow = x0 + (i-1) * dxc
	  else
	    ylow = y0 + (i-1) * dyc
	  end if
	  xhigh = xlow + dxc
	  yhigh = ylow + dyc
	  if( col .ge. 0. ) then
	    call qsetc(col)
	    call qrfill(xlow,ylow,xhigh,yhigh)
	  end if
	end do

	if( bdebug ) call write_color_table
	!call reset_auto_color_table
	if( bdebug ) call write_color_table

!  write legend

	call qcm(xcm,ycm)
	dxtick = 0.1*xcm
	dytick = 0.1*ycm

	if( bhoriz ) then
	  call qtxtr(0.)
	  call qtxtcr(0.,+2.5)
	else
	  if( brotate ) then
	    call qtxtr(90.)
	    call qtxtcr(0.,+2.0)
	  else
	    call qtxtcr(-1.,0.)
	  end if
	end if

	dtick = 1.
	ntk = nticks
	blog = ipllog .ne. 0
	if( blog ) ntk = -1

	if( ntk .eq. 0 .or. ntk .gt. isoanz ) then	!real boxes, all colors
	  ntk = isoanz
	  dtick = 1.
	  y2 = y1
	  x2 = x0
	else if( ntk .gt. 0 ) then			!only ticks, skip some
	  if( ntk .eq. 1 ) ntk = 2			!just in case...
	  dtick = float(isoanz-1) / (ntk-1)
	  y2 = y0-dytick
	  x2 = x1+dxtick
	else						!negative -> log
	  idiv = ipllog
	  y2 = y0-dytick
	  x2 = x1+dxtick
	  amin = getcolval(1.)
	  amax = getcolval(float(isoanz))
	  ntk = ndim
	  ntks = nticks
	  call logvals(amin,amax,idiv,ntk,ntks,rval,aval)
	  !write(6,*) 'amin,amax: ',amin,amax
	  !write(6,*) 'idiv,ntk: ',idiv,ntk
	end if

	if( bhoriz ) then
	  dxy = x1-x0-2.*dxc
	else
	  dxy = y1-y0-2.*dyc
	end if

	if( ntk .eq. 1 ) then				!bug fix for 1 isoline
	  dt = 0.
	else
	  dt = dxy / (ntk-1)				!distance between ticks
	end if

	call qgray(0.)

	do i=1,ntk
	  rit = 1. + (i-1)*dtick
	  value = getcolval(rit)
	  if( value > 2.00E+09 ) goto 99
	  if( blog ) value = aval(i)
	  value = fact * value
	  idec = ndec
	  if( idec .eq. 0 ) idec = ideci(value)
	  nw = ialfa(value,line,idec,-1)
	  x = x0 + dxc + (i-1) * dt
	  y = y0 + dyc + (i-1) * dt
	  if( blog ) then
	    x = x0 + dxc + rval(i)*dxy
	    y = y0 + dyc + rval(i)*dxy
	  end if
	  if( bhoriz ) then
	    call qline(x,y0,x,y2)
	    call qtext(x,y0,line)
	  else
	    call qline(x1,y,x2,y)
	    call qtext(x1+2.*dxtick,y,line)
	  end if
	end do

!  plot box around colorbar

        call qgray(0.)
	call pbox(x0,y0,x1,y1)
	!call pbox(xmin,ymin,xmax,ymax)		!just for debug

!  write legend over colorbar

        !call qtxts(15)
	call qtxts(fs_cbar_text)
        
	if( bhoriz ) then
	  call qtxtcr(0.,-1.9)
	  call qtext((x0+x1)/2.,y1,title)
	else
	  call qtxtr(90.)
	  call qtxtcr(0.,-2.0)
	  call qtext(x0,(y0+y1)/2.,title)
	end if

	call qtxtr(0.)
	call qtxtcr(-1.,-1.)

	return
 1000	format(i5,4f14.4,5x,a10)
   99	continue
	write(6,*) 'value = ',value
	stop 'error stop colbar: value too big to handle'
	end

! ************************************************************
! ************************************************************
! ************************************************************

        subroutine legini

!  initializes legend

	use legend

        implicit none

	integer icall
	save icall
	data icall /0/

	if( icall .ne. 0 ) return
	icall = 1

	nleg = 0
	nlegdi = legdim
	iplotleg = 1		!plot legend

	end

! ************************************************************

	subroutine set_plotleg(ipl)

	use legend

	implicit none

	integer ipl

	iplotleg = ipl

	end

! ************************************************************

        subroutine legerr

!  error message for legend

	use legend

        implicit none

	if( nlegdi .ne. legdim ) then
                write(6,*) 'internal error: nlegdi .ne. legdim'
                stop 'error stop legerr: nlegdi'
        else
                write(6,*) 'too many legends'
                write(6,*) 'increase legdim in legend.h'
                stop 'error stop legerr: legdim'
        end if

        end

! ************************************************************

        subroutine legrd

!  reads legend from str file

	use legend

        implicit none

        character*80 line
        integer i

	logical is_letter
	integer nrdlin
	integer ichafs

	call legini
	if( nlegdi .ne. legdim ) call legerr

        do while( nrdlin(line) .gt. 0 )
	  i = ichafs(line)			!first non blank char
	  if( line(i:i) .ne. '#' .and. i .ne. 0 ) then	!trash comments

	    if( is_letter( line(i:i) ) ) then
		call newleg(line(i:))
	    else
		write(6,*) 'first part of legend must be keyword:'
		write(6,*) trim(line(i:))
		write(6,*) 'error processing legend'
		write(6,*) '(old legend format has been discontinued)'
		stop 'error stop legrd: no keyword found'
	    end if
	  end if
	    
        end do

	return
        end

! ************************************************************

	subroutine newleg(line)

!  reads legend - new version
! 
!  text	30500 11800     15      'Chioggia'	#text, 15pt
!  vart	30500 11800     15      'a'		#text, 15pt, change on each call
!  ctxt	0 0					#text centering (horiz, vert)
!  ignu	1 					#ignore underscore (0/1)
!  line	30500 11800 35000 15000			#line
!  vect	30500 11800 35000 15000	0.1		#arrow, tipsize 0.1
!  rect	30500 11800 35000 15000	0.1		#rectangle, fill color 0.1
!  rect	30500 11800 35000 15000	-1		#rectangle (outline, no fill)
!  circ	30500 11800 5000 -1			#circle (outline, no fill)
!  wid	5					#set line width to 5
!  col	0.5					#set color to 0.5
! 
!  values for ctxt:
! 		-1:left/bottom  0:center  +1:right/top
!  possible text for vart:
! 		a	gives a b c d ...
! 		A	gives A B C D ...
! 		1	gives 1 2 3 4 ...

	use legend

	implicit none

	character*(*) line

	character*80 text,keyword
	logical bdebug
	integer ioff
	integer iston,istof,istos

	character*4 what
	real x,y,x2,y2,size,arrow,color
	ioff = 1

	bdebug = .false.

	x = 0
	y = 0
	x2 = 0
	y2 = 0
	size = 0
	arrow = 0
	color = 0
	text = ' '
	what = 'none'

	if( iston(line,keyword,ioff) .le. 0 ) goto 99

	what = keyword(1:4)

	if( keyword .eq. 'text' ) then
	  if( istof(line,x,ioff) .le. 0 ) goto 99
	  if( istof(line,y,ioff) .le. 0 ) goto 99
	  if( istof(line,size,ioff) .le. 0 ) goto 99
	  if( istos(line,text,ioff) .le. 0 ) goto 99
	else if( keyword .eq. 'vart' ) then
	  if( istof(line,x,ioff) .le. 0 ) goto 99
	  if( istof(line,y,ioff) .le. 0 ) goto 99
	  if( istof(line,size,ioff) .le. 0 ) goto 99
	  if( istos(line,text,ioff) .le. 0 ) goto 99
	else if( keyword .eq. 'ctxt' ) then
	  if( istof(line,x,ioff) .le. 0 ) goto 99
	  if( istof(line,y,ioff) .le. 0 ) goto 99
	else if( keyword .eq. 'ignu' ) then
	  if( istof(line,size,ioff) .le. 0 ) goto 99
	else if( keyword .eq. 'rect' ) then
	  if( istof(line,x,ioff) .le. 0 ) goto 99
	  if( istof(line,y,ioff) .le. 0 ) goto 99
	  if( istof(line,x2,ioff) .le. 0 ) goto 99
	  if( istof(line,y2,ioff) .le. 0 ) goto 99
	  if( istof(line,color,ioff) .le. 0 ) goto 99
	else if( keyword .eq. 'circ' ) then
	  if( istof(line,x,ioff) .le. 0 ) goto 99
	  if( istof(line,y,ioff) .le. 0 ) goto 99
	  if( istof(line,arrow,ioff) .le. 0 ) goto 99	!is radius
	  if( istof(line,color,ioff) .le. 0 ) goto 99
	else if( keyword .eq. 'vect' ) then
	  if( istof(line,x,ioff) .le. 0 ) goto 99
	  if( istof(line,y,ioff) .le. 0 ) goto 99
	  if( istof(line,x2,ioff) .le. 0 ) goto 99
	  if( istof(line,y2,ioff) .le. 0 ) goto 99
	  if( istof(line,arrow,ioff) .le. 0 ) goto 99
	else if( keyword .eq. 'line' ) then
	  if( istof(line,x,ioff) .le. 0 ) goto 99
	  if( istof(line,y,ioff) .le. 0 ) goto 99
	  if( istof(line,x2,ioff) .le. 0 ) goto 99
	  if( istof(line,y2,ioff) .le. 0 ) goto 99
	else if( keyword .eq. 'col' ) then
	  if( istof(line,color,ioff) .le. 0 ) goto 99
	else if( keyword .eq. 'wid' ) then
	  if( istof(line,size,ioff) .le. 0 ) goto 99
	else
	  goto 98
	end if

	nleg = nleg + 1
	if( nleg .gt. legdim ) call legerr

	if( bdebug ) then
	      write(6,*) '------------'
	      write(6,*) line
	      write(6,*) x,y,size,x2,y2
	      write(6,*) arrow,color
	      write(6,*) text
	      write(6,*) '------------'
	end if
            
	xleg(1,nleg) = x
	yleg(1,nleg) = y
	xleg(2,nleg) = x2
	yleg(2,nleg) = y2
	aleg(nleg) = arrow
	cleg(nleg) = color
	legsiz(nleg) = nint(size)
	legleg(nleg) = text
	whatleg(nleg) = what

	return
   98	continue
	write(6,*) 'unknown keyword ',keyword
	write(6,*) line
	stop 'error stop newleg: read error'
   99	continue
	write(6,*) 'error reading legend: '
	write(6,*) line
	stop 'error stop newleg: read error'
	end

! ************************************************************

        subroutine legplo

!  plots legend

	use legend

        implicit none

        logical bdebug
	character*4 what,what1
	character*80 text
	integer i,isize
        real color,size,rad
	real dcolor		!default color
        real scale

	integer getisec

	integer icall
	save icall
	data icall /0/

	call legini
	if( nlegdi .ne. legdim ) call legerr

	if( iplotleg == 0 ) return	!do not plot legend

        bdebug = .true.
        bdebug = .false.

	dcolor = 0.

	call qcomm('Plotting user defined legend')
        call qwhite(.true.)
        call qfont('Times-Roman')
        call qgray(dcolor)
	call qtxts(12)
	call qlwidth(0.01)

	if( bdebug ) then
	  write(6,*) 'debug legend: ',nleg,nlegdi
	  write(4,*) 'debug legend: ',nleg,nlegdi
	end if

	do i=1,nleg

	  what = whatleg(i)
	  size = aleg(i)
	  isize = legsiz(i)
          color = cleg(i)
          rad   = aleg(i)
	  text = legleg(i)

	  if( what .eq. 'text' ) then
	    if( isize .gt. 0 ) call qtxts(isize)
	    call make_absolute1(xleg(1,i),yleg(1,i))
	    call qtext(xleg(1,i),yleg(1,i),text)
	  else if( what .eq. 'vart' ) then
	    call adjust_vartext(icall,text)
	    if( isize .gt. 0 ) call qtxts(isize)
	    call make_absolute1(xleg(1,i),yleg(1,i))
	    call qtext(xleg(1,i),yleg(1,i),text)
	  else if( what .eq. 'ctxt' ) then
	    call qtxtcr(xleg(1,i),yleg(1,i))
	  else if( what .eq. 'ignu' ) then
	    call qignu(isize)
	  else if( what .eq. 'vect' ) then
	    if( isize .gt. 0 ) call qlwidth(isize/100.)
            if( color .gt. 0. ) call qgray(color)
	    call make_absolute2(xleg(1,i),yleg(1,i),xleg(2,i),yleg(2,i))
	    call fpfeil(xleg(1,i),yleg(1,i),xleg(2,i),yleg(2,i),size)
            if( color .gt. 0. ) call qgray(dcolor)
	  else if( what .eq. 'line' ) then
	    if( isize .gt. 0 ) call qlwidth(isize/100.)
            if( color .gt. 0. ) call qgray(color)
	    call make_absolute2(xleg(1,i),yleg(1,i),xleg(2,i),yleg(2,i))
	    call qline(xleg(1,i),yleg(1,i),xleg(2,i),yleg(2,i))
            if( color .gt. 0. ) call qgray(dcolor)
	  else if( what .eq. 'rect' ) then
	    call make_absolute2(xleg(1,i),yleg(1,i),xleg(2,i),yleg(2,i))
            if( color .ge. 0. ) then
	      call qgray(color)
	      call qrfill(xleg(1,i),yleg(1,i),xleg(2,i),yleg(2,i))
	      call qgray(dcolor)
	    else
	      call pbox(xleg(1,i),yleg(1,i),xleg(2,i),yleg(2,i))
	    end if
	  else if( what .eq. 'circ' ) then
	    call make_absolute1(xleg(1,i),yleg(1,i))
            if( color .ge. 0. ) then
	      call qgray(color)
	      call qarcf(xleg(1,i),yleg(1,i),rad,0.,360.)
	      call qgray(dcolor)
	    else
	      call qarc(xleg(1,i),yleg(1,i),rad,0.,360.)
	    end if
	  else if( what .eq. 'col' ) then
	    dcolor = color
	    call qgray(color)
	  else if( what .eq. 'wid' ) then
	    call qlwidth(isize/100.)
	  else
	    goto 99
	  end if

	  if( bdebug ) then
	    write(6,*) '------------'
	    write(6,*) xleg(1,i),yleg(1,i),xleg(2,i),yleg(2,i)
	    write(6,*) isize
	    write(6,*) size,color
	    write(6,*) text
	    write(6,*) '------------'
	  end if

	end do

	call qgray(0.)
        call qwhite(.false.)
	call qlwidth(0.)              !reset

	if( getisec() == 0 ) call legdate	!not a section plot
	call legwind

	icall = icall + 1

	return
   99	continue
	write(6,*) 'Unknown keyword ',what
	stop 'error stop legplo: keyword error'
	end

! *****************************************************************

	subroutine adjust_vartext(iadd,text)

	implicit none

	integer iadd
	character*(*) text

	integer al,au,zl,zu,one,nine
	parameter (al=ichar('a'),au=ichar('A'))
	parameter (zl=ichar('z'),zu=ichar('Z'))
	parameter (one=ichar('1'),nine=ichar('9'))

	integer i,ic,add
	character*1 c

	do i=1,len(text)
	  c = text(i:i)
	  ic = ichar(c)
	  add = 0
	  if( ic >= al .and. ic <= zl ) add = iadd
	  if( ic >= au .and. ic <= zu ) add = iadd
	  if( ic >= one .and. ic <= nine ) add = iadd
	  ic = ic + add
	  text(i:i) = char(ic)
	end do

	write(6,*) 'adjusting vartext: ',iadd,trim(text)

	end

! *****************************************************************

        subroutine legwind

!  plots wind vector

        implicit none

        real x,y
        real scale
        integer isize

        integer ndim
        parameter(ndim=20)

        real u,v,s,d
	real xt,yt
	real fact,afact
        real t,sc,r
        real wind(2)
	double precision dtime
        integer it,nintp,nvar,nread
        character*80 file,text

	integer, save :: idwind = 0

        real array(ndim)
        save array
        real xwind,ywind,lwwind,scwind,xtwind,ytwind
        save xwind,ywind,lwwind,scwind,xtwind,ytwind
        integer iwtype,stwind
        save iwtype,stwind
        character*40 wtext,wunit
        save wtext,wunit

	real getpar

        integer icall
        save icall
        data icall / 0 /

	if( icall .eq. -1 ) return

        if( icall .eq. 0 ) then
          icall = 1
          file = 'wind.dat'
	  call getfnm('wfile',file)
	  if( file .eq. ' ' ) then
	    icall = -1
	    return
	  end if
          nintp = 2
          nvar = 2
	  nread = 0
	  call ptime_get_dtime(dtime)
	  call iff_ts_init(dtime,file,nintp,nvar,idwind)
          !call exffil(file,nintp,nvar,nread,ndim,array)
          xwind = getpar('xwind')
          ywind = getpar('ywind')
          iwtype = nint(getpar('iwtype'))
          lwwind = getpar('lwwind')
          scwind = getpar('scwind')
	  call make_absolute1(xwind,ywind)

          xtwind = getpar('xtwind')
          ytwind = getpar('ytwind')
          stwind = nint(getpar('stwind'))
	  call getfnm('wtext',wtext)
	  call getfnm('wunit',wunit)
	  call make_absolute1(xtwind,ytwind)
        end if

        if( iwtype .le. 0 ) return

	x = xwind
	y = ywind

	call qcomm('Plotting wind vector')
        call qwhite(.true.)
        call qfont('Times-Roman')
        call qgray(0.)
	call qtxts(12)
	call qlwidth(lwwind)

	!call ptime_get_itime(it)
        !t = it
	call ptime_get_dtime(dtime)
	call iff_ts_intp(idwind,dtime,wind)
        !call exfintp(array,t,wind)

        if( iwtype .eq. 1 ) then
          u = wind(1)
          v = wind(2)
          s = sqrt(u*u+v*v)
          d = 0.
        else if( iwtype .eq. 3 .or. iwtype .eq. 4 ) then
          s = wind(1)
          d = wind(2)
          if( iwtype .eq. 4 ) s = s * 1852. / 3600.
          call convert_wind(s,d,u,v)
        else
          write(6,*) 'iwtype = ',iwtype
          stop 'error stop legwind: impossible value for iwtype'
        end if

        write(6,*) 'wind legend: ',x,y,s,d,u,v,scwind
	call spherical_fact(fact,afact)
	!u = u / fact
        !call pfeil(x,y,u,v,scwind)
	!r = s*scwind
	!call pcircle(x,y,r)
	xt = x + u*scwind/fact
	yt = y + v*scwind
        call fcmpfeil(x,y,xt,yt,0.3)
	call qlwidth(-1.)      !FIXME -> use negative number to reset

	if( stwind .gt. 0 ) call qtxts(stwind)
        call make_wind_text(wtext,wunit,s,text)
	call qtext(xtwind,ytwind,text)

        end

! *****************************************************************

        subroutine make_wind_text(wtext,wunit,s,text)

        implicit none

        character*(*) wtext,wunit,text
        real s

        integer i
        integer ichanm

        text = wtext
        i = ichanm(text)

        write(text(i+2:),'(f4.1)') s
        i = ichanm(text)

        text(i+2:) = wunit

        end

! *****************************************************************

        subroutine legdate

!  plots date legend

        implicit none

        integer it,iday,ihour
        integer jd,year,month,day
        integer date,time
        character*25 line
        character*3 name

        real xdate,ydate
        save xdate,ydate
        integer sdate,idate
        save sdate,idate

	real, save :: tzshow
	integer itl

        integer icall
        save icall
        data icall /0/

	real getpar
	double precision dgetpar

	if( icall .eq. -1 ) return

        if( icall .eq. 0 ) then

          idate = 1
          date = 20010101
          time = 0

          idate = nint(getpar('idate'))

	  if( idate .le. 0 ) then
	    icall = -1
	    return
	  end if

          date = nint(dgetpar('date'))
          time = nint(dgetpar('time'))
	  if( date .ne. -1 ) then	!if given overwrites date in sim
	    write(6,*) 'initializing date and time: ',date,time
            call dtsini(date,time)
	  end if

          xdate = getpar('xdate')
          ydate = getpar('ydate')
          sdate = nint(getpar('sdate'))
          tzshow = getpar('tzshow')

	  call make_absolute1(xdate,ydate)
	  !write(6,*) '%%%%%%%%%%%%% ',xdate,ydate,idate,date

          icall = 1
        end if

        if( idate .le. 0 ) return

	call qcomm('Plotting date legend')
        call qwhite(.true.)
        call qfont('Times-Roman')
        call qgray(0.)
	call qtxts(12)
	call qlwidth(-1.)

	call ptime_get_itime(it)

	itl = it + nint(tzshow*3600)		!correct for time zone

        if( idate .eq. 1 ) then
          call dtsgf(itl,line)
	  !write(6,*) 'date/time for plot: ',itl,'  ',line
        else if( idate .eq. 2 ) then
          iday = itl / 86400
          ihour = (itl - iday*86400 ) / 3600       !not yet finished
          year = 2002
          jd = iday
          if( jd .le. 0 ) jd = 1
          !write(6,*) 'legdate: ',itl,iday,jd,ihour
          call j2date(jd,year,month,day)
          call month_name(month,name)
          !write(line,'(a,i2,1x,a3,1x,i4)') 'data ',day,name,year
          write(line,'(i2,1x,a3,1x,i4)') day,name,year
        else if( idate .eq. 3 ) then
          call dtsgf(itl,line)
	  line(11:12) = '  '
	  !write(6,*) 'date/time for plot: ',itl,'  ',line
        else if( idate .eq. 4 ) then
          call dtsgf(itl,line)
	  line(11:12) = '  '
	  line(23:25) = 'GMT'
	  !write(6,*) 'date/time for plot: ',itl,'  ',line
        else if( idate .eq. 5 ) then
          call dtsgf(itl,line)
	  line(11:) = '  '
        else
          write(6,*) 'idate = ',idate
          stop 'error stop legdate: impossible value for idate'
        end if

	if( sdate .gt. 0 ) call qtxts(sdate)
	call qtext(xdate,ydate,line)

        end

! *****************************************************************

	subroutine velarr(title,unit,scale,ndec,fact &
     &				,arrlen,sclvel,velmax,x0,y0,x1,y1)

!  creates legend for velocity arrow and plots arrow

	use mod_bash

	implicit none

	character*(*) title,unit
	integer ndec
	real scale		!factor from world to physical coordinates
	real fact		!factor for representation purpose
	real arrlen		!length of arrow in physical units
	real sclvel		!factor for arrow
	real velmax		!maximum velocity in plot
	real x0,y0,x1,y1	!dimension of legend on plot

	character*80 line
	integer nw
	real dx,dy
	real x0s,y0s		!starting point of velocity arrow
	real dd			!length of velocity arrow
	real val,val1,val2,value
	real width,height
	real sfact,afact,xt
	real ablank			!factor for blanking

	integer iround,istell,ialfa
	real rnext
	real getpar

	logical bblank,bdebug

	bblank = .true.		!blanking for velocity arrow
	ablank = getpar('ablank')
	bblank = ablank > 0
	bdebug = .true.
	bdebug = .false.

!  blanking of window

	if( bblank ) call blank_window(ablank,x0,y0,x1,y1)

!  extension of window

	dx = x1 - x0
	dy = y1 - y0

!  find length of arrow
! 
!  dd	length of velocity vector in world coordinates
!  val	length of velocity vector in velocity (transport) coordinates

	if( arrlen .le. 0. ) then	!automatic length
	  dd = 0.5 * dx			!maximum length is half of space
	  val1 = dd / scale		!same in velocity coords.
	  val1 = rnext(val1,-3)		!find lower 1,2,3,4,5,8
	  val2 = rnext(velmax,-3)	!find lower 1,2,3,4,5,8
	  val = sclvel * min(val1,val2)
	else
	  val = arrlen
	end if
	dd = val * scale		!dd is length in world coords.

        call qcomm('Plotting arrow')
        call qtxts(12)
        call qfont('Times-Roman')
	call qtxtr(0.)
	call qtxtcr(-1.,-1.)

	x0s = x0 + 0.5 * ( dx - dd )
	y0s = y0 + 0.5 * dy

	!call pfeil(x0s,y0s,val,0.,scale)

	call spherical_fact(sfact,afact)
	xt = x0s + val*scale/sfact
        !call fcmpfeil(x0s,y0s,xt,0.,0.3)
        call fcmpfeil(x0s,y0s,xt,y0s,0.3)

	call qtsize(title,width,height)
	x0s = x0 + 0.5 * dx - 0.5 * width
	call qtext(x0s,y0s+1.*height,title)

	value = val*fact
	if( bbverb ) write(6,*) 'arrow: ',val,value,fact
	nw = ialfa(value,line,ndec,-1)
	line(nw+1:) = unit
	call qtsize(line,width,height)
	x0s = x0 + 0.5 * dx - 0.5 * width
	call qtext(x0s,y0s-2.*height,line)

	if( bdebug) then
	  write(6,*) 'start velarr ================================='
	  write(6,*) dd,x0s,y0s,xt
	  write(6,*) '  end velarr ================================='
	end if

        call qcomm('end of plotting arrow')

	end

! ******************************************************************

	subroutine velsh(ivel,scale,velmax)

!  shell for velocity arrow
! 
!  ivel = 1      velocities
!  ivel = 2      transports
!  ivel = 3      wind
!  ivel = 4      waves

	implicit none

	integer ivel
	real scale
	real velmax		!maximum velocity in plot

	character*80 title,unit
	integer ndec
	real fact,arrlen,sclvel,eps
	real x0,y0,x1,y1
	logical bdebug

	integer iround
	real getpar
	logical inboxdim

	bdebug = .true.
	bdebug = .false.

	eps = 1.e-3

	if( .not. inboxdim('arr',x0,y0,x1,y1) ) return

	fact = getpar('facvel')
	arrlen = getpar('arrvel')
	sclvel = getpar('sclvel')
	ndec = nint(getpar('ndcvel'))
	call getfnm('legvel',title)

	if( ivel .eq. 2 ) then
	  unit = ' m**2/s'
	  if( fact .ne. 1 ) then
	    call make_unit(fact,unit)
	    !write(6,*) 'Cannot handle factor ',fact,' with ivel ',ivel
	    !stop 'error stop velsh: facvel'
	  end if
	else
	  if( fact .eq. 1. ) then
	    unit = ' m/s'
	  else if( fact .eq. 100. ) then
	    unit = ' cm/s'
	  else if( fact .eq. 1000. ) then
	    unit = ' mm/s'
          else if( abs(fact-1.9438445) < eps ) then   !ccf
            unit = ' Kn'
	  else
	    unit = ' m/s'
	    call make_unit(fact,unit)
	    !write(6,*) 'Cannot handle factor ',fact,' with ivel ',ivel
	    !stop 'error stop velsh: facvel'
	  end if
	end if

	if( bdebug) then
	  write(6,*) 'start velsh ================================='
	  write(6,*) fact,'   ',unit
	  write(6,*) '  end velsh ================================='
	end if

	call velarr(title,unit,scale,ndec,fact,arrlen,sclvel &
     &				,velmax,x0,y0,x1,y1)

	end

! ******************************************************************

	subroutine make_unit(fact,unit)

	implicit none

	real fact
	character*(*) unit

	integer ie
	real rstell
	character*80 line,stell

	integer ialfa

	rstell = nint( log10(fact) )
	rstell = -rstell

	ie = ialfa(rstell,stell,-1,-1)

	line = ' 10**' // stell(1:ie) // unit
	line = ' 10E' // stell(1:ie) // unit
	unit = line

	end

! ******************************************************************

	subroutine setboxdim(name,x0,y0,x1,y1)

!  sets dimension of box

	implicit none

	character*(*) name
	real x0,y0,x1,y1

	character*10 x0var,y0var,x1var,y1var

	x0var = 'x0' // name
	y0var = 'y0' // name
	x1var = 'x1' // name
	y1var = 'y1' // name

	call putpar(x0var,x0)
	call putpar(y0var,y0)
	call putpar(x1var,x1)
	call putpar(y1var,y1)

	end

! ******************************************************************

	subroutine getboxdim(name,x0,y0,x1,y1)

!  gets dimension of box

	implicit none

	character*(*) name
	real x0,y0,x1,y1

	character*10 x0var,y0var,x1var,y1var

	real getpar

	x0var = 'x0' // name
	y0var = 'y0' // name
	x1var = 'x1' // name
	y1var = 'y1' // name

	x0 = getpar(x0var)
	y0 = getpar(y0var)
	x1 = getpar(x1var)
	y1 = getpar(y1var)

	end

! ******************************************************************

	function inboxdim(name,x0,y0,x1,y1)

	implicit none

	logical inboxdim
	character*(*) name
	real x0,y0,x1,y1

	logical inboxdim_intern

	inboxdim = inboxdim_intern(name,.true.,x0,y0,x1,y1)

	end

! ******************************************************************

	function inboxdim_noabs(name,x0,y0,x1,y1)

	implicit none

	logical inboxdim_noabs
	character*(*) name
	real x0,y0,x1,y1

	logical inboxdim_intern

	inboxdim_noabs = inboxdim_intern(name,.false.,x0,y0,x1,y1)

	end
	
! ******************************************************************

	function inboxdim_intern(name,babs,x0,y0,x1,y1)

!  computes size of box inside plot

	implicit none

	logical inboxdim_intern
	character*(*) name
	logical babs
	real x0,y0,x1,y1

	inboxdim_intern = .false.

	call getboxdim(name,x0,y0,x1,y1)

! 	write(6,*) '------ ',name,' ----',x0,y0,x1,y1

	if( x0 .eq. x1 .or. y0 .eq. y1 ) return	!nothing to do

	if( x0 .gt. x1 ) call swapr(x0,x1)
	if( y0 .gt. y1 ) call swapr(y0,y1)

	inboxdim_intern = .true.

! 	we are nearly done -> see if units are relative and convert

	if( babs ) call make_absolute(x0,y0,x1,y1)

	end

! ******************************************************************

	subroutine swapr(a,b)
	c=a
	a=b
	b=c
	end

! ******************************************************************

	function is_box_given(name)

!  checks if box is given

	implicit none

	logical is_box_given
	character*(*) name

	real x0,y0,x1,y1	!we are not interested in these values
	logical inboxdim_noabs

	is_box_given = inboxdim_noabs(name,x0,y0,x1,y1)

	end

! ******************************************************************

	function is_relative(x,y)

!  checks if coordinates are relative or absolute

	implicit none

	logical is_relative
	real x,y

	real relmin,relmax

! 	relmin =  0.
! 	relmax = +1.
	relmin = -1.
	relmax = +2.

	is_relative = .false.

	if( x .lt. relmin .or. x .gt. relmax ) return
	if( y .lt. relmin .or. y .gt. relmax ) return

	is_relative = .true.

	end

! ******************************************************************

	subroutine make_abs(x,y)

!  makes coordinates absolute

	implicit none

	real x,y

	real xmin,ymin,xmax,ymax

	call getbas(xmin,ymin,xmax,ymax)

	x = xmin + x * (xmax-xmin)
	y = ymin + y * (ymax-ymin)

	end

! ******************************************************************

	subroutine make_absolute1(x,y)

!  makes coordinates absolute (1 pair)

	implicit none

	real x,y

	logical is_relative

	if( is_relative(x,y) ) then
	  call make_abs(x,y)
	end if

	end

! ******************************************************************

	subroutine make_absolute2(x0,y0,x1,y1)

!  makes coordinates absolute (2 pairs)

	implicit none

	real x0,y0,x1,y1

	logical is_relative

	if( is_relative(x0,y0) .and. is_relative(x1,y1) ) then
	  call make_abs(x0,y0)
	  call make_abs(x1,y1)
	end if

	end

! ******************************************************************

	subroutine make_absolute(x0,y0,x1,y1)

!  makes absolute coordinates -> old version

	implicit none

	real x0,y0,x1,y1

	real relmin,relmax
	real xmin,ymin,xmax,ymax

! 	relmin =  0.
! 	relmax = +1.
	relmin = -1.
	relmax = +2.

	if( x0 .lt. relmin .or. x0 .gt. relmax ) return
	if( y0 .lt. relmin .or. y0 .gt. relmax ) return
	if( x1 .lt. relmin .or. x1 .gt. relmax ) return
	if( y1 .lt. relmin .or. y1 .gt. relmax ) return

! 	units are relative -> convert

	call getbas(xmin,ymin,xmax,ymax)

! 	write(6,*) '------ minmax ----',xmin,ymin,xmax,ymax

	x0 = xmin + x0 * (xmax-xmin)
	y0 = ymin + y0 * (ymax-ymin)
	x1 = xmin + x1 * (xmax-xmin)
	y1 = ymin + y1 * (ymax-ymin)

	end

! ******************************************************************

	subroutine occupy_dim(x0,y0,x1,y1)

!  uses info from occupy to set up x/y for color bar

	use mod_bash

	implicit none

	real x0,y0,x1,y1

	integer iquad		!quadrante that is free

	call occupy(iquad)

	if( iquad .eq. 1 ) then
	  x0 = 0.1
	  x1 = 0.4
	  y0 = 0.1
	  y1 = 0.2
	else if( iquad .eq. 2 ) then
	  x0 = 0.6
	  x1 = 0.9
	  y0 = 0.1
	  y1 = 0.2
	else if( iquad .eq. 3 ) then
	  x0 = 0.6
	  x1 = 0.9
	  y0 = 0.7
	  y1 = 0.8
	else if( iquad .eq. 4 ) then
	  x0 = 0.1
	  x1 = 0.4
	  y0 = 0.7
	  y1 = 0.8
	end if

        if( bbverb ) write(6,*) 'occupy_dim : ',iquad,x0,y0,x1,y1

	end

! ******************************************************************

	subroutine occupy(iquad)

!  computes free area for legend -> must still be used by another routine

!  quadrante:
! 
! 	+---+---+
! 	| 4 | 3 |
! 	+---+---+
! 	| 1 | 2 |
! 	+---+---+

	use basin
	use mod_bash

	implicit none

	integer iquad		!quadrante that is free

	integer ndim
	parameter (ndim=20)

	integer ia(ndim,ndim)

	integer iqc(4)

	logical bwrite
	character*80 line
	integer i,j,k,n
	integer iq,jq,ip
	integer iqcmin,iqmin
	real x,y,dx,dy
	real x0,y0,x1,y1

	integer ijq(2,2)	!qudrante
	save ijq
	data ijq /1,2,4,3/

	bwrite = bbverb
	bwrite = .false.

	do j=1,ndim
	  do i=1,ndim
	    ia(i,j) = 0
	  end do
	end do

	call getbas(x0,y0,x1,y1)

	dx = x1 - x0
	dy = y1 - y0

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  i = 1 + ndim*(x-x0)/dx
	  j = 1 + ndim*(y-y0)/dy
	  if( i .ge. 1 .and. i .le. ndim ) then
	    if( j .ge. 1 .and. j .le. ndim ) then
	      ia(i,j) = 1
	    end if
	  end if
	end do

	n = min(ndim,80)

	if( bwrite ) then
	 write(6,*) 'occupy... ',n,ndim,x0,y0,x1,y1
	 do j=n,1,-1
	  do i=1,n
	    if( ia(i,j) .ne. 0 ) then
	      line(i:i) = '*'
	    else
	      line(i:i) = '.'
	    end if
	  end do
	  write(6,*) line(1:n)
	 end do
	end if

	do i=1,4
	  iqc(i) = 0
	end do

	do j=1,ndim
	  jq = 1 + 2*(j-1)/ndim
	  do i=1,ndim
	    iq = 1 + 2*(i-1)/ndim
	    ip = ijq(iq,jq)
	    !write(6,*) i,j,iq,jq,ip
	    if( ia(i,j) .ne. 0 ) iqc(ip) = iqc(ip) + 1
	  end do
	end do
	
	iqcmin = ndim*ndim
	iqmin = 0
	do i=1,4
	  if( iqc(i) .lt. iqcmin ) then
	    iqcmin = iqc(i)
	    iqmin = i
	  end if
	end do

	!write(6,*) iqmin,iqc

	iquad = iqmin

	end

! ******************************************************************

        subroutine blank_window(rf,x0,y0,x1,y1)

!  blanks window given by coordinates

        implicit none

	real rf
        real x0,y0,x1,y1

	real dx,dy,xm,ym
        real x0aux,y0aux,x1aux,y1aux

	if( rf <= 0 ) return

	if( rf /= 1. ) then
	  xm = 0.5 * (x1 + x0)
	  ym = 0.5 * (y1 + y0)
	  dy = y1 - y0
	  dx = x1 - x0
	  dy = y1 - y0
	  dx = 0.5 * dx * rf
	  dy = 0.5 * dy * rf
	  x0aux = xm - dx
	  x1aux = xm + dx
	  y0aux = ym - dy
	  y1aux = ym + dy

	  !write(6,*) x0,y0,x1,y1
	  !write(6,*) x0aux,y0aux,x1aux,y1aux
	end if

        call qcomm('Start blanking window')
	call qwhite(.true.)
	call qgray(1.)	!color is white
	if( rf /= 1. ) then
	  call qrfill(x0aux,y0aux,x1aux,y1aux)
	else
	  call qrfill(x0,y0,x1,y1)
	end if
	call qgray(0.)
	call qwhite(.false.)
        call qcomm('End blanking window')

        end

! ******************************************************************


