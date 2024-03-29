
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2019  Georg Umgiesser
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
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

! geometrical routines
!
! contents :
!
! function vproj(u1,v1,u2,v2)				projection of vector
! function icut(x1,y1,x2,y2,x3,y3,x4,y4,xout,yout)	intersects two lines
! function intria(x,y,xp,yp)				point in triangle or not
! subroutine mirr(x1,y1,x2,y2,xp,yp)			reflection of point
! subroutine intrfk(x,y,xp,yp,fk)			interpolates point (fk)
! function rintrz(x,y,xp,yp,zp)				interpolates point (zp)
! function icutri(x,y,x1,y1,x2,y2,ds)			intersects line & tri.
!
!***************************************************************
!
	function vproj(u1,v1,u2,v2)
!
! projects a vector (u1,v1) onto another vector (u2,v2)
!
! (u1,v1)	vector to project
! (u2,v2)	vector on whom to project (u1,v1)
! vproj		return value, giving the length of the projected
!		...vector in terms of the second one. The
!		...projected vector can be obtained multiplying
!		...the second vector by vproj :
!		...   (up,vp) = vproj * (u2,v2)
!
! formulas :	u is to be projected onto v
!		... |uproj| = |u| * cos(alfa)    and
!		...  uproj  = |uproj| * v / |v|    with
!		... cos(alfa) = (u * v) / (|u| * |v|)
!
	real vproj

	uv2 = u2*u2 + v2*v2
!
	if(uv2.lt.1.e-6) then
		write(6,*) 'Cannot project on nill vector'
		write(6,*) 'Length of second vector : ',uv2
		vproj = 0.	!$$ALPHA
		return
	end if
!
	vproj = (u1*u2 + v1*v2) / uv2
!
	return
	end
!
!*************************************************************
!
	function icut(x1,y1,x2,y2,x3,y3,x4,y4,xout,yout)
!
! intersects two lines and gives back intersection point
!
! x1,y1		starting point of first line
! x2,y2		ending point of first line
! x3,y3		starting point of second line
! x4,y4		ending point od second line
! xout,yout	point of intersection (if any)
! icut		return code :
!		-1 : lines parallel (but not equal)
!		-2 : lines equal
!		0  : lines intersect
!		1  : intersection inside 1 but outside 2
!		2  : intersection inside 2 but outside 1
!		3  : intersection outside of both lines
!		...(not-negative return code returns intersection
!		...point in xout,yout)
!
! formulas :	linear system
!				a*r1 + b*r2 = e
!				c*r1 + d*r2 = f
!
	parameter (eps1=1.e-6,eps2=1.e-6)
!	parameter (rmin=0.,rmax=1.)
	parameter (rmin=0.-eps1,rmax=1.+eps1)
!
	a=x2-x1
	b=x3-x4
	c=y2-y1
	d=y3-y4
	e=x3-x1
	f=y3-y1
!
	det = a*d-c*b
!
	if(abs(det).lt.eps2) then	!lines parallel
		if(abs(a*f-c*e).lt.eps2) then	!lines equivalent
			icut=-2
		else
			icut=-1
		end if
		return
	else				!compute intersection
		r1=(e*d-f*b)/det
		r2=(a*f-c*e)/det
!
		if(r1.ge.rmin.and.r1.le.rmax) then
			if(r2.ge.rmin.and.r2.le.rmax) then
				icut=0
			else
				icut=1
			end if
			xout=x1+r1*(x2-x1)
			yout=y1+r1*(y2-y1)
		else
			if(r2.ge.rmin.and.r2.le.rmax) then
				icut=2
				xout=x3+r2*(x4-x3)
				yout=y3+r2*(y4-y3)
			else
				icut=3
				if(abs(r1).lt.abs(r2)) then
					xout=x1+r1*(x2-x1)
					yout=y1+r1*(y2-y1)
				else
					xout=x3+r2*(x4-x3)
					yout=y3+r2*(y4-y3)
				end if
			end if
		end if
	end if
!
	return
	end
!
!***********************************************************
!
	function intria(x,y,xp,yp)
!
! point in triangle or not
!
! x,y		array of coordinates of vertices of triangle
! xp,yp		coordinates of point
! intria	1: point is in triangle  0: point outside (return value)
!
	dimension x(3),y(3)
	data eps /1.e-7/
!
	xs=0.
	ys=0.
	do i=1,3
	   xs=xs+x(i)
	   ys=ys+y(i)
	end do
	xs=xs/3.
	ys=ys/3.
!
	intria=0
!
	do k1=1,3
	   k2=mod(k1,3)+1
	   x12=x(k1)-x(k2)
	   y12=y(k1)-y(k2)
	   det=(xs-xp)*y12-(ys-yp)*x12
	   if(abs(det).ge.eps) then
		detlam=(x(k1)-xp)*y12-(y(k1)-yp)*x12
		rlamb=detlam/det
		if(rlamb.gt.0..and.rlamb.lt.1.) return	!outside
	   end if
	end do
!
	intria=1	!inside
!
	return
	end
!
!****************************************************************
!
	function mirr(x1,y1,x2,y2,xp,yp)
!
! reflection of one point on a line
!
! x1,y1,x2,y2	points defining line
! xp,yp		on entry : point to be reflected
! 		on return : reflected point
! mirr		0 : ok     -1 : error
!
	u=x2-x1
	v=y2-y1
!
	uv = u*u + v*v
!
	if(uv.lt.1.e-6) then
		write(6,*) 'Cannot reflect on nill vector'
		write(6,*) 'Length of second vector : ',uv
		mirr=-1
		return
	end if
!
	r=(u*(xp-x1)+v*(yp-y1))/uv
!
	xp=2.*(r*u+x1)-xp
	yp=2.*(r*v+y1)-yp
!
	mirr=0
!
	return
	end
!
!*********************************************************************
!
	subroutine intrfk(x,y,xp,yp,fk)
!
! interpolates point in triangle and returns fk
!
! x,y		vertices of triangle
! xp,yp		coordinates of point
! fk		form functions (return value)
!
! if c(x(i),y(i)),i=1,3 are the values at vertices,
! ... then c(xp,yp) can be obtained by :
!
!		c(xp,yp) = c(1)*fk(1) + c(2)*fk(2) + c(3)*fk(3)
!
	dimension x(3),y(3),fk(3)
!
	f=0.
	do i=1,3
	   ii=mod(i,3)+1
	   iii=mod(ii,3)+1
	   a=x(ii)*y(iii)-x(iii)*y(ii)
	   b=y(ii)-y(iii)
	   c=x(iii)-x(ii)
	   fk(i) = a + xp*b + yp*c
	   f=f+a
	end do
!
	do i=1,3
	   fk(i) = fk(i)/f
	end do
!
	return
	end
!
!*********************************************************************
!
	function rintrz(x,y,xp,yp,zp)
!
! interpolates point in triangle and returns value
!
! x,y		vertices of triangle
! xp,yp		coordinates of point
! zp		values at vertices
! rintrz	interpolated value at (xp,yp)
!
	dimension x(3),y(3),zp(3)
!
	f=0.
	z=0.
	do i=1,3
	   ii=mod(i,3)+1
	   iii=mod(ii,3)+1
	   a=x(ii)*y(iii)-x(iii)*y(ii)
	   b=y(ii)-y(iii)
	   c=x(iii)-x(ii)
	   fk = a + xp*b + yp*c
	   z = z + zp(i)*fk
	   f=f+a
	end do
!
	rintrz = z/f
!
	return
	end
!
!********************************************************************
!
	function icutri(x,y,x1,y1,x2,y2,ds)
!
! intersects line with sides of triangle
!
! x,y		vertices of triangle
! x1,y1		coordinates of starting point of line
! x2,y2		coordinates of ending point of line
! ds		distance of intersection with vertex
!		...ds(i) = [0...1] intersection with line i
!		...ds(i) = -1 no intersection with line i
!		...(line i is opposit to vertex i)
!		...if ds(i) = [0...1] and xp,yp is intersection point
!		...then ds(i) = sqrt ( (xp-x(i+1))**2 + (yp-y(i+1))**2 )
! icutri	number of sides intersect with line
!
	dimension x(3),y(3),ds(3)
!
	icutri=0
!
	do i=1,3
	   ii=mod(i,3)+1
	   iii=mod(ii,3)+1
	   x3=x(ii)
	   y3=y(ii)
	   x4=x(iii)
	   y4=y(iii)
	   isw=icut(x1,y1,x2,y2,x3,y3,x4,y4,xp,yp)
	   if(isw.eq.0) then
		rpunt = (xp-x3)*(xp-x3) + (yp-y3)*(yp-y3)
		rside = (x4-x3)*(x4-x3) + (y4-y3)*(y4-y3)
		ds(i) = sqrt( rpunt/rside )
		icutri = icutri +1
	   else
		ds(i)=-1.
	   end if
	end do
!
	return
	end
