
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-1999,2009-2012,2014-2015,2014-2015  Georg Umgiesser
!    Copyright (C) 2017-2020  Georg Umgiesser
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

! general utility routines
!
! contents :
!
! function izahl(z,ndec)		computes ciphers of number
! function laezi(l,ioff)		computes length of number
! function iantw(l)			gets answer from terminal
! function inquir(l,f)			gets numbers from terminal
! function getnum(prompt)		gets number from terminal
! function igetxt(prompt,text)		gets text from terminal
! function ialfa(zahl,lh,ndec,mode)	converts number into alphanumeric chars
! subroutine uplow(text,type)		translates text to upper/lower case 
! subroutine tablnc(linold,linnew)	converts tabs into blanks
! subroutine prilin(line,iunit)		writes a line without trailing blanks
! function ichanm(line)			computes length of line
! function ichafs(line)			first occurrence of non blank character 
! function icompr(line)			eliminates blank characters from line
! function ifndch(line,iocc,cha)	finds the iocc occurrence of cha 
! function ifndoc(line,string)		finds number of occurences of string
! function xzentr(x,height,f,ndec)	funtion centers number
! function ideflt(k,text)		writes default and gets value (integer)
! function fdeflt(r,text)		writes default and gets value (real)
! function iround(f)			rounds real value
! function itypch(char)			returns typ of character
! subroutine tabula(tab)		creates tab character
! function rnext(r,mode)		finds closest value to r
! function rnexta(r,mode)		finds closest value to r (absolute)
! function istell(r)			computes exponent of r
! subroutine scplo(xmin,xmax,ymin,ymax)	scales plot
! subroutine swapr(a,b)			swaps values
! function rlen(x1,y1,x2,y2)		length of line
!
! uplow,itypch,tabula depend on the ASCII character set
!
! revision log :
!
! 01.06.1998	ggu	new routines iscan... into own file (subscn.f)
! 17.06.1998	ggu	old routines with new name (iscan0,ialfa0)
! 12.02.1999	ggu	new routine triml
! 26.01.2009	ggu	minor changes to avoid compiler warnings
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2011	ggu	new routine trimline()
! 30.03.2012	ggu	changed VERS_6_1_51
! 30.05.2014	ggu	new routine rnextsub()
! 12.12.2014	ggu	changed VERS_7_0_9
! 06.05.2015	ggu	new routine logvals() for logarithmic values
! 21.05.2015	ggu	changed VERS_7_1_11
! 10.07.2015	ggu	changed VERS_7_1_50
! 14.09.2015	ggu	new routine i2s0()
! 25.05.2017	ggu	changed VERS_7_5_28
! 14.11.2017	ggu	changed VERS_7_5_36
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.02.2020	ggu	rounding routines into new file subround.f
! 11.12.2020	ggu	ichafs() and ichanm() moved to subscn.f
!
!***********************************************************
!
	function izahl(z,ndec)
!
! computes ciphers of number
!
! z		number for which ciphers have to be determined
! ndec		positions after the decimal point
!		(-1 : no decimal point)
! izahl		ciphers of number z with minus sign and
!		...decimal point if necessary (return value)
!
	if(ndec.lt.-1) then
		ndech=iabs(ndec)-1
		zh=z
	else if(ndec.ge.1) then
		ndech=ndec
		zh=z
	else
		ndech=0
		zh=z
	end if
!
	iz=abs(zh)+0.01
!
	if(iz.eq.0) then
		istel=1
	else
		istel=0
	end if
!
	do while (iz.gt.0)
	iz=iz/10
	istel=istel+1
	end do
!
	if(z.lt.0.) istel=istel+1
	if(ndec.ge.0) then
		istel=istel+1+ndec
	else if(ndec.lt.-1) then
		if(istel+1.lt.ndech) then
			istel=istel+1
		else
			istel=ndech
		end if
	end if
!
	izahl=istel
!
	return
	end
!
!******************************************
!
	function iantw(l)
!
! gets answer from terminal
!
! l		text written to terminal
! iantw		1 : true
!		0 : false
!
! the following answers can be given as true or false :
!
!		y j s t 1	==>	true
!		n f 0 'blank'	==>	false
!
	character*(*) l
!
	character l1*1
	data net,nat /5,6/
!
    1	continue
!
	write(nat,*) l
	read(net,1000) l1
 1000	format(a)
!
	call to_lower(l1)
!
	if( l1 .eq. 'y' .or. l1 .eq. 'j' .or. l1 .eq.'s' .or. l1 .eq. 't' .or. l1 .eq. '1' ) then
		iantw=1
	else if(l1 .eq. 'n' .or. l1 .eq. 'f' .or. l1 .eq. '0' .or. l1 .eq. ' ') then
		iantw=0
	else
		write(nat,*) 'Incorrect answer. Try again.'
		goto 1
	end if
!
	return
	end

!********************************************

	function inquire_numbers(l,f,max)

! gets numbers from terminal
!
! do not use -> no way to know if we are out of bounds
!
! l		text written to terminal
! f		array in which the values are stored (return value)
! inquire_numbers	total number of values read in

	implicit none

	integer inquire_numbers
	character*(*) l
	real f(max)
	integer max

	character*80 lh
	integer net,nat
	integer iscanf
	data net,nat /5,6/

	lh=' '

	write(nat,1000) l
 1000	format(1x,a)

	read(net,2000) lh
 2000	format(a)

	inquire_numbers=iscanf(lh,f,max)

	end

!********************************************
!
	function getnum(prompt)
!
! gets number from terminal
!
! prompt	text written to terminal
! getnum	value of number 
!
	character*(*) prompt
!
	character*80 lh
	real f(1)
	data net,nat /5,6/
!
    1	continue
!
	write(nat,1000) prompt
 1000	format(1x,a)
!
	lh=' '
	read(net,2000) lh
 2000	format(a)
!
	if(iscanf(lh,f,1).ne.1) then
		write(nat,1000) 'Erroneous input. Repeat.'
		goto 1
	end if
!
	getnum=f(1)
!
	return
	end
!
!********************************************
!
	subroutine tablnc(linold,linnew)
!
! converts tabs into blanks
!
! linold	old line with tabs
! linnew	new converted line
!
! itab		number of blanks equals one tab (for fortran minimum 6)
!
	character*(*) linold,linnew
!
	character*1 tab
	data itab   / 8 /
!
	call tabula(tab)
	lengo=len(linold)
	lengn=len(linnew)
!
	iold=1
	inew=1
!
	do while(iold.le.lengo.and.inew.le.lengn)
!
	if(linold(iold:iold).ne.tab) then
		linnew(inew:inew)=linold(iold:iold)
		iold=iold+1
		inew=inew+1
	else
		ianf=inew
		iend=itab*(1+((inew-1)/itab))
		if(iend.gt.lengn) iend=lengn
		do i=ianf,iend
		   linnew(inew:inew)=' '
		   inew=inew+1
		end do		
		iold=iold+1
	end if
!
	end do
!
	return
	end
!
!*********************************************
!
	function ichanm_0(line)
!
! computes length of line without trailing blanks
!
! line		line of text
! ichanm	length of line (return value)
!		... 0 : line is all blank
!
	character*(*) line
!
	character*1 blank,tab
	data blank /' '/
!
	call tabula(tab)
!
	ndim=len(line)
!
	do i=ndim,1,-1
	if(line(i:i).ne.blank.and.line(i:i).ne.tab) goto 1
	end do
!
    1	continue
	ichanm_0=i
!
	return
	end
!
!*********************************************
!
	function ichafs_0(line)
!
! computes first occurrence of a non-blank character in line
!
! line		line of text
! ichafs	position of first non-blank character (return value)
!		... 0 : no non-blank character found
!
	character*(*) line
!
	character*1 blank,tab
	data blank /' '/
!
	call tabula(tab)
!
	ndim=len(line)
!
	do i=1,ndim
	if(line(i:i).ne.blank.and.line(i:i).ne.tab) goto 1
	end do
!
	i=0
!
    1	continue
	ichafs_0=i
!
	return
	end
!
!*********************************************
!
	function ifndch(line,iocc,cha)
!
! finds the iocc occurrence of cha in line from offset ioff on
!
! line		line of text
! iocc		occurrence of character
!		... >0 : iocc occurrence is found
!		...  0 : last occurrence is found
!		... -1 : last but one occurrence is found...
! cha		character(s) to be found
! ifndch	position of cha in line (return value)
!		...  0 : cha not found
!		... -1 : error in ioff or ndim
!
	character*(*) line,cha
!
	n=0
!
	ncha=len(cha)-1
	ndim=len(line)
!
	if(iocc.gt.0) then
		iend=ndim-ncha
		do i=1,iend
		if(line(i:i+ncha).eq.cha(1:1+ncha)) then
			n=n+1
			if(n.eq.iocc) then
				ifndch=i
				return
			end if
		end if
		end do
	else
		ianf=ndim-ncha
		do i=ianf,1,-1
		if(line(i:i+ncha).eq.cha(1:1+ncha)) then
			n=n+1
			if(1-n.eq.iocc) then
				ifndch=i
				return
			end if
		end if
		end do
	end if
!
	ifndch=0
!
	return
	end
!
!************************************************************
!
	function fdeflt(r,text)
!
! writes default and gets value (real)
!
! r		default value (real)
! text		text written to terminal
! fdeflt	real value chosen
!
	character*(*) text
!
	character line*80
	real f(1)
!
    1	continue
!
	write(6,*) text,' (default = ',r,' )'
	read(5,'(a)') line
!
	ianz=iscanf(line,f,1)
!
	if(ianz.eq.1) then
		fdeflt=f(1)
	else if(ianz.eq.0) then
		fdeflt=r
	else if(ianz.gt.1) then
		write(6,*) 'Only one number at a time.'
		goto 1
	else
		write(6,*) 'Read error'
		goto 1
	end if
!
	return
	end
!
!***********************************************************
!
	function iround(f)
!
! rounds real value
!
! f		real value
! iround	rounded integer value
!
	if(f.ge.0) then
		iround=f+0.5
	else
		iround=f-0.5
	end if
!
	return
	end
!
!****************************************************************
!
	subroutine tabula(tab)
!
! creates tab character
!
! works only for ASCII set of characters
!
! tab		tab character (return value)
!
	character*1 tab,char
!
	tab=char(9)
!
	return
	end

!*********************************************************

	function istell(r)
!
! computes exponent of r
!
! -2     0.01 <= |r| < 0.1
! -1      0.1 <= |r| < 1
!  0        1 <= |r| < 10
!  1       10 <= |r| < 100
!  2      100 <= |r| < 1000
!
	implicit none

	integer istell
	real r

	integer i
	real eps,ar

	ar = abs(r)
	eps = 0.0001*ar

	i = log10(ar+eps)

	if(ar+eps.lt.1.) i=i-1

	istell = i

	return
	end

!************************************************************

	subroutine i2s0(number,string)

! converts integer to string with blanks substituted by zeros

	implicit none

	integer number
	character*(*) string

	integer i,l
	character*10 format

	l = len(string)
	write(format,'(a,i2,a)') '(i',l,')'
	write(string,format) number

	do i=1,l
	  if( string(i:i) == ' ' ) string(i:i) = '0'
	end do

	end

!************************************************************

