c
c $Id: lagrange_util.f,v 1.1 2009-02-13 17:22:44 georg Exp $
c
c general utilities for lagrangian model
c
c revision log :
c
c 05.02.2009    ggu     copied from other files
c
c*******************************************************************

        function dist_node(k1,k2)

c returns distance between two nodes

        implicit none

        real dist_node
        integer k1,k2

        real xgv(1), ygv(1)
        common /xgv/xgv, /ygv/ygv

        real x1,y1,x2,y2,dx,dy

        x1 = xgv(k1)
        y1 = ygv(k1)
        x2 = xgv(k2)
        y2 = ygv(k2)

        dx = x1 - x2
        dy = y1 - y2

        dist_node = sqrt( dx*dx + dy*dy )

        end

c*******************************************************************

	subroutine xy_minmax(n,x,y,xmin,xmax,ymin,ymax)

c returns min/max coordinates of polygon

	implicit none

	integer n
	real x(1),y(1)
	real xmin,xmax
	real ymin,ymax

	integer i

	xmin = x(1)
	xmax = x(1)
	ymin = y(1)
	ymax = y(1)

	do i=2,n
	  xmin = min(xmin,x(i))
	  xmax = max(xmax,x(i))
	  ymin = min(ymin,y(i))
	  ymax = max(ymax,y(i))
	end do

	end

c*******************************************************************

	subroutine basin_center(xm,ym)

c returns center of gravity of total basin

	implicit none

	real xm,ym

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real xgv(1)
	common /xgv/xgv
	real ygv(1)
	common /ygv/ygv

	call xy_center(nkn,xgv,ygv,xm,ym)

	end
	
c*******************************************************************

	subroutine xy_center(n,x,y,xm,ym)

c returns center of gravity of polygon

	implicit none

	integer n
	real x(1),y(1)
	real xm,ym

	integer i

	xm = 0.
	ym = 0.

	do i=1,n
	  xm = xm + x(i)
	  ym = ym + y(i)
	end do

	xm = xm / n
	ym = ym / n

	end

c*******************************************************************

	subroutine xy_and_minmax(ie,x,y,xmin,xmax,ymin,ymax)

c returns x,y and min/max coordinates of vertices of element ie

	implicit none

	integer ie
	real x(3),y(3)
	real xmin,xmax
	real ymin,ymax

	integer nen3v(3,1)
	common /nen3v/nen3v
	real xgv(1)
	common /xgv/xgv
	real ygv(1)
	common /ygv/ygv

	integer ii,k

	do ii=1,3
          k=nen3v(ii,ie)
          x(ii)=xgv(k)
          y(ii)=ygv(k)
        end do

        xmin=min(x(1),x(2),x(3))
        xmax=max(x(1),x(2),x(3))
        ymin=min(y(1),y(2),y(3))
        ymax=max(y(1),y(2),y(3))

	end

c*******************************************************************

	subroutine compute_total_area(area)

	implicit none

	real area

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	include 'ev.h'

	integer ie
	double precision a,tot_area

	tot_area = 0.

	do ie=1,nel
	  a = ev(10,ie)
	  tot_area = tot_area + a
	end do

	area = 12. * tot_area

	end

c*******************************************************************
