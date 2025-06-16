
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-1999,2004-2005,2009-2020  Georg Umgiesser
!    Copyright (C) 2011  Marco Bajo
!    Copyright (C) 2011  Christian Ferrarin
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

! routines for interpolation onto regular grid
!
! contents :
!
! subroutine setgeo(x0,y0,dx,dy,flag)
!		sets grid values for interpolation
! subroutine getgeo(x0,y0,dx,dy,flag)
!		gets grid values for interpolation
! subroutine getgeoflag(flag)
!		gets flag value for interpolation
!
! subroutine av2am(av,am,ip,jp)
!		interpolation of av onto a regular net
! subroutine av2amk(bwater,av,am,ip,jp)
!		interpolation of av onto a regular net
! function intri(x,y,xp,yp)
!		point in triangle or not
! function intrid(x,y,xp,yp)
!		point in triangle or not (double precision version)
! subroutine am2av(am,av,ip,jp)
!		interpolation of am onto finite element mesh
! function am2val(am,ip,jp,xx,yy)
!		interpolation of am onto finite element mesh
! subroutine ave2am(av,am,ip,jp)
!		interpolation of av (elementwise) onto a regular net
!
! subroutine mkmask(bwater,zv,href,hzoff)
!		makes mask in element
!
! subroutine mimareg(am,ip,jp,amin,amax)
!		computes min/max of regular matrix (without flag values)
! subroutine a2char(am,ac,ip,jp)
!		creates 1 char representation of matrix
! subroutine prchar(ac,ip,jp)
!		prints 1 char representation of matrix
!
! subroutine femintp(ie,z,xp,yp,zp)
!               interpolation in element (with ev)
! subroutine elemintp(x,y,z,xp,yp,zp)
!               interpolation in element (no ev)
!
! subroutine find_elem_from_old(ieold,xp,yp,ielem)
!		finds element for point (xp,yp) starting from ieold
! subroutine find_element(xp,yp,ielem)
!		finds element for point (xp,yp)
! function in_element(ie,xp,yp)
!		checks if point (xp,yp) is in element ie
! subroutine get_xy_elem(ie,x,y)
!		returns x,y of vertices of element ie
!
! revision log :
!
! 18.11.1998	ggu	routine commented
! 18.11.1998	ggu	routine setgeo introduced
! 19.11.1998	ggu	routines a2char, prchar added
! 19.10.1999	ggu	routine mkmask added from subutl
! 25.11.2004	ggu	new routines femintp and elemintp for interpolation
! 14.03.2005	ggu	new routines for interpolation in element
! 11.03.2009	ggu	new helper routine getgeoflag()
! 12.06.2009	ggu	passing to double precision, intrid, bug bug_f_64bit
! 23.03.2010	ggu	changed v6.1.1
! 26.01.2011	ggu&mbj	handling extrapolation in am2av()
! 27.01.2011	ggu&ccf	bug fix in find_elem_from_old() BUG_27.01.2011
! 31.03.2011	ggu	new routine elemmask()
! 14.04.2011	ggu	changed VERS_6_1_22
! 07.06.2011	ggu	changed VERS_6_1_25
! 24.11.2011	ggu	new routine find_close_elem()
! 09.12.2011	ggu	changed VERS_6_1_38
! 24.01.2012	ggu	changed VERS_6_1_41
! 30.03.2012	ggu	changed VERS_6_1_51
! 20.06.2012	ggu	new routine get_scal_elem()
! 07.10.2012	ggu	new routine av2fm()
! 10.10.2012	ggu	new routine fm2am2d() and fm2am3d()
! 25.10.2012	ggu	changed VERS_6_1_59
! 26.10.2012	ggu	bug fix: do not access not existing storage
! 05.11.2012	ggu	changed VERS_6_1_60
! 25.01.2013	ggu	changed VERS_6_1_62
! 30.05.2014	ggu	in av2amk() do not interpolate for flag values
! 18.06.2014	ggu	changed VERS_6_1_77
! 27.06.2014	ggu	changed VERS_6_1_78
! 07.07.2014	ggu	new routine intp_reg()
! 26.11.2014	ggu	changed VERS_7_0_7
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 05.05.2015	ggu	changed VERS_7_1_10
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 25.09.2015	ggu	new routines intp_reg_nodes(), intp_reg_elems()
! 16.12.2015	ggu	changed VERS_7_3_16
! 18.12.2015	ggu	changed VERS_7_3_17
! 19.02.2016	ggu	changed VERS_7_5_2
! 28.04.2016	ggu	changed VERS_7_5_9
! 05.05.2016	ggu	file restructured (module)
! 14.05.2016	ggu	allow for extension of grid -> bregextend
! 25.05.2016	ggu	changed VERS_7_5_10
! 10.06.2016	ggu	changed VERS_7_5_13
! 14.06.2016	ggu	changed VERS_7_5_14
! 17.06.2016	ggu	changed VERS_7_5_15
! 23.06.2016	ggu	allow for eps in computing box
! 23.09.2016	ggu	allow for eps in computing box and reg intp
! 30.09.2016	ggu	changed VERS_7_5_18
! 11.10.2016	ggu	changed VERS_7_5_20
! 20.01.2017	ggu	changed VERS_7_5_22
! 23.04.2017	ggu	new routine intp_reg_single_nodes()
! 09.05.2017	ggu	changed VERS_7_5_26
! 23.05.2017	ggu	file split into subreg, submask and subfind
! 11.07.2017	ggu	changed VERS_7_5_30
! 02.09.2017	ggu	changed VERS_7_5_31
! 17.11.2017	ggu	changed VERS_7_5_37
! 04.12.2017	ggu	check t,u values and correct if out of bounds
! 18.05.2018	ggu	more checks on fr routines, introduced ierr, bextend
! 26.05.2018	ggu	even more checks and debug output
! 06.07.2018	ggu	changed VERS_7_5_48
! 13.07.2018	ggu	changed VERS_7_4_1
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 14.05.2019	ggu	in fm_extra_setup() use double precision
! 21.05.2019	ggu	changed VERS_7_5_62
! 01.04.2020	ggu	new routine make_reg_box()
! 11.04.2022	ggu	new routine condense_valid_coordinates for single nodes
! 25.11.2022	ggu	new routine recollocate_nodes() and ipg array
! 20.02.2025	ggu	new routine create_reg()
!
! notes :
!
! pxareg,pyareg         coordinates of lower left point of matrix
! pxdreg,pydreg         grid size of matrix
! pzlreg                value of z for land points
!
!******************************************************

!==================================================================
        module regular
!==================================================================

	implicit none

	logical, save :: bregextend = .false.

	real, save :: pxareg = 0.	!x coordinate of lower,left point (x0)
	real, save :: pyareg = 0.	!y coordinate of lower,left point (y0)
	real, save :: pxdreg = 0.	!grid spacing in x direction (dx)
	real, save :: pydreg = 0.	!grid spacing in y direction (dy)
	real, save :: pzlreg = -999.	!flag for land points

!==================================================================
        end module regular
!==================================================================

	subroutine setgeo(x0,y0,dx,dy,flag)

! sets grid values for interpolation
!
! pxareg,pyareg         coordinates of lower left point of matrix
! pxdreg,pydreg         grid size of matrix
! pzlreg                value of z for land points

	use regular

	implicit none

	real x0,y0	!coordinates of lower,left point
	real dx,dy	!grid spacing in x/y direction
	real flag	!flag for land points

	pxareg = x0
	pyareg = y0
	pxdreg = dx
	pydreg = dy
	pzlreg = flag

	end

!******************************************************

	subroutine getgeo(x0,y0,dx,dy,flag)

! gets grid values for interpolation
!
! pxareg,pyareg         coordinates of lower left point of matrix
! pxdreg,pydreg         grid size of matrix
! pzlreg                value of z for land points

	use regular

	implicit none

	real x0,y0	!coordinates of lower,left point
	real dx,dy	!grid spacing in x/y direction
	real flag	!flag for land points

	x0   = pxareg
	y0   = pyareg
	dx   = pxdreg
	dy   = pydreg
	flag = pzlreg

	end

!******************************************************

	subroutine getgeoflag(flag)

! gets flag value for interpolation
!
! flag                value for land points

	use regular

	implicit none

	real flag	!flag for land points

	flag = pzlreg

	end

!******************************************************

	subroutine setregextend(bextend)

! sets flag to decide if extend interpolated grid

	use regular

	implicit none

	logical bextend

	bregextend = bextend

	end

!******************************************************

	subroutine getregextend(bextend)

! gets flag to decide if extend interpolated grid

	use regular

	implicit none

	logical bextend

	bextend = bregextend

	end

!******************************************************
!******************************************************
!******************************************************

	subroutine create_reg(dx,dy,xmin,ymin,xmax,ymax,regpar)

	implicit none

	real, intent(in) :: dx,dy
	real, intent(in) :: xmin,ymin,xmax,ymax
	real, intent(out) :: regpar(7)

	integer nx,ny
	real x0,y0,x1,y1
	real flag

        if( dx <= 0. .or. dy <= 0. ) goto 96
        if( xmin >= xmax ) goto 96
        if( ymin >= ymax ) goto 96

        x0 = int(xmin/dx)*dx
        y0 = int(ymin/dy)*dy

        nx = 2 + (xmax - xmin) / dx
        x1 = x0 + (nx-1)*dx
        if( x1 < xmax ) nx = nx + 1
        x1 = x0 + (nx-1)*dx

        ny = 2 + (ymax - ymin) / dy
        y1 = y0 + (ny-1)*dy
        if( y1 < ymax ) ny = ny + 1
        y1 = y0 + (ny-1)*dy

        if( x0 > xmin .or. x1 < xmax ) goto 95
        if( y0 > ymin .or. y1 < ymax ) goto 95

	call getgeoflag(flag)
	call setreg(regpar,nx,ny,x0,y0,dx,dy,flag)

	return
   95   continue
        write(6,*) x0,xmin,xmax,x1
        write(6,*) y0,ymin,ymax,y1
        stop 'error stop create_reg: internal error (1)'
   96   continue
        write(6,*) 'dx,dy: ',dx,dy
        write(6,*) 'xmin,xmax: ',xmin,xmax
        write(6,*) 'ymin,ymax: ',ymin,ymax
        stop 'error stop create_reg: error in input parameters'
	end

!******************************************************

	subroutine setreg(regpar,nx,ny,x0,y0,dx,dy,flag)

	implicit none

	real regpar(7)
	integer nx,ny
	real x0,y0,dx,dy
	real flag

	regpar(1) = nx
	regpar(2) = ny
	regpar(3) = x0
	regpar(4) = y0
	regpar(5) = dx
	regpar(6) = dy
	regpar(7) = flag

	end

!******************************************************

	subroutine getreg(regpar,nx,ny,x0,y0,dx,dy,flag)

	implicit none

	real regpar(7)
	integer nx,ny
	real x0,y0,dx,dy
	real flag

	nx = nint(regpar(1))
	ny = nint(regpar(2))
	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)
	flag = regpar(7)

	end

!******************************************************

	subroutine printreg(regpar)

	implicit none

	real regpar(7)
	integer nx,ny
	real x0,y0,dx,dy,x1,y1
	real flag

        call getreg(regpar,nx,ny,x0,y0,dx,dy,flag)

        x1 = x0 + (nx-1)*dx
        y1 = y0 + (ny-1)*dy
        write(6,'(4x,a,2i16)') 'nx,ny: ',nx,ny
        write(6,'(4x,a,2f16.6)') 'x0,y0: ',x0,y0
        write(6,'(4x,a,2f16.6)') 'x1,y1: ',x1,y1
        write(6,'(4x,a,2f16.6)') 'dx,dy: ',dx,dy
        write(6,'(4x,a,2f16.6)') 'flag : ',flag

	end

!******************************************************

	subroutine extend_reg(next,nx,ny,zreg)

	implicit none

	integer next
	integer nx,ny
	real zreg(nx,ny)

	integer n,nflag
	real flag

	n = 0
	do
	  n = n + 1
	  if( next > 0 .and. n > next ) exit		!done with extension
	  call extend_reg_by_one(nx,ny,zreg,nflag)
	  if( nflag == 0 ) exit				!no more flagged values
	end do

	end

!******************************************************

	subroutine extend_reg_by_one(nx,ny,zreg,nflag)

	implicit none

	integer nx,ny
	real zreg(nx,ny)
	integer nflag

	integer n,ix,iy,ixx,iyy
	double precision accum
	real, allocatable :: zaux(:,:)
	real flag

	call getgeoflag(flag)

	allocate(zaux(nx,ny))
	zaux = zreg
	nflag = 0

	do iy=1,ny
	  do ix=1,nx
	    if( zreg(ix,iy) /= flag ) cycle
	    nflag = nflag + 1
	    n = 0
	    accum = 0.
	    ixx = ix-1
	    if( ixx >= 1 .and. zreg(ixx,iy) /= flag ) then
	      n = n + 1
	      accum = accum + zreg(ixx,iy)
	    end if
	    ixx = ix+1
	    if( ixx <= nx .and. zreg(ixx,iy) /= flag ) then
	      n = n + 1
	      accum = accum + zreg(ixx,iy)
	    end if
	    iyy = iy-1
	    if( iyy >= 1 .and. zreg(ix,iyy) /= flag ) then
	      n = n + 1
	      accum = accum + zreg(ix,iyy)
	    end if
	    iyy = iy+1
	    if( iyy <= ny .and. zreg(ix,iyy) /= flag ) then
	      n = n + 1
	      accum = accum + zreg(ix,iyy)
	    end if
	    if( n > 0 ) zaux(ix,iy) = accum / n
	  end do
	end do

	zreg = zaux

	end

!******************************************************
!******************************************************
!******************************************************

	subroutine find_position_to_coord(x,y,ix,iy)

! finds closest position (ix,iy) to coordinate (x,y) in reg grid

	implicit none

	real x,y
	integer ix,iy

	real x0,y0,dx,dy,flag

	call getgeo(x0,y0,dx,dy,flag)

	ix = nint( (x-x0)/dx + 1. )
	iy = nint( (y-y0)/dy + 1. )

	end

!******************************************************

	subroutine find_coord_to_position(ix,iy,x,y)

! finds coordinate (x,y) to given position (ix,iy) in reg grid

	implicit none

	integer ix,iy
	real x,y

	real x0,y0,dx,dy,flag

	call getgeo(x0,y0,dx,dy,flag)

	x = x0 + (ix-1)*dx
	y = y0 + (iy-1)*dy

	end

!******************************************************
!******************************************************
!******************************************************

	subroutine av2am(av,am,ip,jp)

! interpolation of av onto a regular net (nodal values)
!
! av                    array to be interpolated
! am                    matrices of interpolated values (u,v,z,h)
! ip,jp                 dimension of matrices
!
! pxareg,pyareg         coordinates of lower left point of matrix
! pxdreg,pydreg         grid size of matrix
! pzlreg                value of z for land points

	use basin

	implicit none

	integer ip,jp
	real av(nkn)
	real am(ip,jp)

	logical bwater(nel)

	bwater = .true.

	call av2amk(bwater,av,am,ip,jp)

	end

!******************************************************

	subroutine av2amk(bwater,av,am,ip,jp)

! interpolation of av onto a regular net (nodal values) with mask
!
! bwater		mask for water points
! av                    array to be interpolated
! am                    matrices of interpolated values (u,v,z,h)
! ip,jp                 dimension of matrices
!
! pxareg,pyareg         coordinates of lower left point of matrix
! pxdreg,pydreg         grid size of matrix
! pzlreg                value of z for land points

	use basin

	implicit none

! arguments
	integer ip,jp
	real av(nkn)
	real am(ip,jp)
	logical bwater(nel)
! parameter
	double precision eps
	parameter ( eps = 1.d-14 )
! local
	real pxareg,pyareg,pxdreg,pydreg,pzlreg
	integer i,j,ii,iii,ie,k,kn,iin
	integer imin,imax,jmin,jmax
	integer iflag
	double precision x(3),y(3),z(3),a(3),b(3),c(3)
	double precision zh,fh,f,xp,yp
	double precision xmin,xmax,ymin,ymax
! function
	integer intrid

	call getgeo(pxareg,pyareg,pxdreg,pydreg,pzlreg)

	am=pzlreg

	do ie=1,nel
	  if( bwater(ie) ) then	!wet
	    iflag = 0
	    do i=1,3
		kn=nen3v(i,ie)
		x(i)=xgv(kn)
		y(i)=ygv(kn)
		z(i)=av(kn)
	        if( z(i) > pzlreg ) iflag = iflag + 1
	    end do
	    if( iflag .ne. 3 ) cycle

	    !f=0.
	    do i=1,3
		ii=mod(i,3)+1
		iii=mod(ii,3)+1
		a(i)=x(ii)*y(iii)-x(iii)*y(ii)
		b(i)=y(ii)-y(iii)
		c(i)=x(iii)-x(ii)
		!f=f+a(i)
	    end do
	    f = c(3)*b(2) - c(2)*b(3)		!bug_f_64bit
	    if( f .le. eps ) goto 99

	    xmax=max(x(1),x(2),x(3))
	    xmin=min(x(1),x(2),x(3))
	    ymin=min(y(1),y(2),y(3))
	    ymax=max(y(1),y(2),y(3))

	    imin=(xmin-pxareg)/pxdreg+1.99
	    imax=(xmax-pxareg)/pxdreg+1.01
	    jmin=(ymin-pyareg)/pydreg+1.99
	    jmax=(ymax-pyareg)/pydreg+1.01

	    if(imin.lt.1) imin=1
	    if(imax.gt.ip)imax=ip
	    if(jmin.lt.1) jmin=1
	    if(jmax.gt.jp)jmax=jp

	    do i=imin,imax
		do j=jmin,jmax
		    xp=(i-1)*pxdreg+pxareg
		    yp=(j-1)*pydreg+pyareg

		    iin=intrid(x,y,xp,yp)

		    if(iin.ne.0) then
			zh=0.
			do k=1,3
			   fh=(a(k)+xp*b(k)+yp*c(k))/f
			   zh=zh+z(k)*fh
			end do
			am(i,j)=zh
		    end if
		end do
	    end do
	  end if
	end do

	return
   99	continue
	write(6,*) ie,f
	write(6,*) x
	write(6,*) y
	write(6,*) a
	write(6,*) b
	write(6,*) c
	stop 'error stop av2amk: area of element'
	end

!************************************************
!************************************************
!************************************************

	subroutine av2fm(fm,ip,jp)

! computation of interpolation matrix (nodal values to regular grid) with mask

	use basin

	implicit none

	real fm(4,ip,jp)	!values for interpolation (fm(4,i,j) = ie)
	integer ip,jp		!dimension of matrices

	logical bwater(nel)	!wet mask for each element

	bwater = .true.

	call av2fmk(bwater,fm,ip,jp)

	end

!************************************************

	subroutine av2fmk(bwater,fm,ip,jp)

! computation of interpolation matrix (nodal values to regular grid) with mask
!
! the interpolation can be carried out as
!
!	do j=1,jp
!	  do i=1,ip
!	    ie = nint(fm(4,i,j))
!	    if( ie .gt. 0 ) then
!	      a = 0.
!	      do ii=1,3
!	        k = nen3v(ii,ie)
!		a = a + val(k) * fm(ii,i,j)
!	      end do
!	    else
!	      a = flag
!	    end if
!	    am(i,j) = a
!	  end do
!	end do
	        
	use basin

	implicit none

! arguments
	logical bwater(nel)	!wet mask for each element
	real fm(4,ip,jp)	!values for interpolation (fm(4,i,j) = ie)
	integer ip,jp		!dimension of matrices
! parameter
	double precision eps
	parameter ( eps = 1.d-14 )
! local
	real pxareg,pyareg,pxdreg,pydreg,pzlreg
	integer i,j,ii,iii,ie,k,kn,iin
	integer imin,imax,jmin,jmax
	logical bok
	double precision x(3),y(3),z(3),a(3),b(3),c(3)
	double precision zh,fh,f,xp,yp
	double precision xmin,xmax,ymin,ymax
! function
	integer intrid

	call getgeo(pxareg,pyareg,pxdreg,pydreg,pzlreg)

	fm = 0.

	do ie=1,nel
	  bok = bwater(ie)
	  if( bok ) then			!wet
	    do i=1,3
		kn=nen3v(i,ie)
		x(i)=xgv(kn)
		y(i)=ygv(kn)
	    end do

	    !f=0.
	    do i=1,3
		ii=mod(i,3)+1
		iii=mod(ii,3)+1
		a(i)=x(ii)*y(iii)-x(iii)*y(ii)
		b(i)=y(ii)-y(iii)
		c(i)=x(iii)-x(ii)
		!f=f+a(i)
	    end do
	    f = c(3)*b(2) - c(2)*b(3)		!bug_f_64bit
	    if( f .le. eps ) goto 99

	    xmax=max(x(1),x(2),x(3))
	    xmin=min(x(1),x(2),x(3))
	    ymin=min(y(1),y(2),y(3))
	    ymax=max(y(1),y(2),y(3))

	    imin=(xmin-pxareg)/pxdreg+1.99
	    imax=(xmax-pxareg)/pxdreg+1.01
	    jmin=(ymin-pyareg)/pydreg+1.99
	    jmax=(ymax-pyareg)/pydreg+1.01

	    if(imin.lt.1) imin=1
	    if(imax.gt.ip)imax=ip
	    if(jmin.lt.1) jmin=1
	    if(jmax.gt.jp)jmax=jp

	    do i=imin,imax
		do j=jmin,jmax
		    xp=(i-1)*pxdreg+pxareg
		    yp=(j-1)*pydreg+pyareg

		    iin=intrid(x,y,xp,yp)

		    if(iin.ne.0) then
			do ii=1,3
			   fh=(a(ii)+xp*b(ii)+yp*c(ii))/f
			   fm(ii,i,j) = fh
			end do
			fm(4,i,j) = ie
		    end if
		end do
	    end do
	  end if
	end do

	return
   99	continue
	write(6,*) ie,f
	write(6,*) x
	write(6,*) y
	write(6,*) a
	write(6,*) b
	write(6,*) c
	stop 'error stop av2fm: area of element'
	end

!************************************************

	subroutine av2fm_single(fm,np,xx,yy)

	use basin

	implicit none

! arguments
	real fm(4,np,1)	!values for interpolation (fm(4,i,j) = ie)
	integer np		!dimension of matrices
	real xx(np),yy(np)
! parameter
	double precision eps
	parameter ( eps = 1.d-14 )
! local
	integer i,j,ii,iii,ie,k,kn,iin
	integer imin,imax,jmin,jmax
	logical bok
	double precision x(3),y(3),z(3),a(3),b(3),c(3)
	double precision zh,fh,f,xp,yp
	double precision xmin,xmax,ymin,ymax
! function
	integer intrid

	fm = 0.

	do ie=1,nel
	    do i=1,3
		kn=nen3v(i,ie)
		x(i)=xgv(kn)
		y(i)=ygv(kn)
	    end do

	    do i=1,3
		ii=mod(i,3)+1
		iii=mod(ii,3)+1
		a(i)=x(ii)*y(iii)-x(iii)*y(ii)
		b(i)=y(ii)-y(iii)
		c(i)=x(iii)-x(ii)
	    end do
	    f = c(3)*b(2) - c(2)*b(3)		!bug_f_64bit
	    if( f .le. eps ) goto 99

	    xmax=max(x(1),x(2),x(3))
	    xmin=min(x(1),x(2),x(3))
	    ymin=min(y(1),y(2),y(3))
	    ymax=max(y(1),y(2),y(3))

	    j = 1
	    do i=1,np
		    xp=xx(i)
		    yp=yy(i)

		    iin=intrid(x,y,xp,yp)

		    if(iin.ne.0) then
			do ii=1,3
			   fh=(a(ii)+xp*b(ii)+yp*c(ii))/f
			   fm(ii,i,j) = fh
			end do
			fm(4,i,j) = ie
		    end if
	     end do
	end do

	return
   99	continue
	write(6,*) ie,f
	write(6,*) x
	write(6,*) y
	write(6,*) a
	write(6,*) b
	write(6,*) c
	stop 'error stop av2fm_single: area of element'
	end

!************************************************

        subroutine fm2am2d(femval,nx,ny,fm,am)

! interpolation 2d of fem values to regular grid using fm matrix

	use basin

        implicit none

        real femval(nkn)		!values of fem array
        integer nx,ny			!dimension of regular matrix
        real fm(4,nx,ny)		!interpolation matrix
        real am(nx,ny)			!interpolated values (return)

	integer nlvdi,nlv
	integer, allocatable :: ilhv(:)

	allocate(ilhv(nel))

	nlvdi = 1
	nlv = 1
	ilhv = 1

        call fm2am3d(nlvdi,ilhv,femval,nlv,nx,ny,fm,am)

	end

!************************************************

        subroutine fm2am3d(nlvdi,ilhv,femval,nlv,nx,ny,fm,am)

! interpolation 3d of fem values to regular grid using fm matrix

	use basin

        implicit none

	integer nlvdi			!vertical dimension of fem array
        integer ilhv(nel)		!vertical discretization (element!!)
        real femval(nlvdi,nkn)		!values of fem array
        integer nlv,nx,ny		!dimension of regular matrix
        real fm(4,nx,ny)		!interpolation matrix
        real am(nlv,nx,ny)		!interpolated values (return)

	logical bflag
        integer i,j,l,lmax,ie,ii,k
        real a
        real flag

	call getgeoflag(flag)

        do j=1,ny
          do i=1,nx
            ie = nint(fm(4,i,j))
            lmax = 0
            if( ie .gt. 0 ) lmax = ilhv(ie)
	    lmax = min(lmax,nlv)
            do l=1,lmax
              a = 0.
	      bflag = .false.
              do ii=1,3
                k = nen3v(ii,ie)
                a = a + femval(l,k) * fm(ii,i,j)
		if( femval(l,k) == flag ) bflag = .true.
              end do
	      if( bflag ) a = flag
              am(l,i,j) = a
            end do
            do l=lmax+1,nlv
              am(l,i,j) = flag
            end do
          end do
        end do

        end

!************************************************
!************************************************
!************************************************

	subroutine fm_extra_setup(nx,ny,fmextra)

! sets up fmextra structure to allow interpolation from fem nodes to reg grid

	use basin

	implicit none

        integer nx,ny			!dimension of regular matrix
	real fmextra(6,nkn)

	logical bout,berror
	integer ix,iy,iix,iiy
	integer j,k
	real xr0,yr0,drx,dry,flag
	double precision x0,y0,dx,dy
	double precision xmax,ymax
	double precision x,y,x1,y1
	double precision eps,tueps
	double precision t,u
	double precision d,w,d2

	integer, save :: jx(4) = (/0,1,1,0/)
	integer, save :: jy(4) = (/0,0,1,1/)
	double precision fmweight(nx,ny)

	eps = 0.01
	tueps = 0.0001
	fmextra = 0.
	fmweight = 0.

	call getgeo(xr0,yr0,drx,dry,flag)
	x0 = xr0
	y0 = yr0
	dx = drx
	dy = dry

	xmax = x0 + (nx-1)*dx
	ymax = y0 + (ny-1)*dy

!	---------------------------------------------------------
!	set up contribution from each fem node to regular grid
!	---------------------------------------------------------

	berror = .false.
	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  ix = (x-x0)/dx+1.
	  iy = (y-y0)/dy+1.
	  if( ix < 1 .or. ix >= nx ) cycle
	  if( iy < 1 .or. iy >= ny ) cycle
	  x1 = x0+(ix-1)*dx
	  y1 = y0+(iy-1)*dy
	  t = (x-x1)/dx
	  u = (y-y1)/dy
	  bout = .false.
	  if( t-1. > tueps .or. t < -tueps ) bout = .true.
	  if( u-1. > tueps .or. u < -tueps ) bout = .true.
	  !if( t.gt.1. .or. t.lt.0. ) bout = .true.
	  !if( u.gt.1. .or. u.lt.0. ) bout = .true.
	  if( bout ) then
	    write(6,*) 'out of domain: '
	    write(6,*) '... ',k,nx,ny,dx,dy
	    write(6,*) '... ',x0,y0,xmax,ymax
	    write(6,*) '... ',x1,y1,x,y
	    write(6,*) '... ',ix,iy,t,u
	    berror = .true.
	  end if
	  fmextra(1,k) = ix
	  fmextra(2,k) = iy
	  do j=1,4
	    iix = ix+jx(j)
	    iiy = iy+jy(j)
	    !if( fm(4,iix,iiy) == 0 ) then
	      x1 = x0+(iix-1)*dx
	      y1 = y0+(iiy-1)*dy
	      d2 = ((x1-x)/dx)**2 + ((y1-y)/dy)**2	!normalized distance
	      d = sqrt( d2 )
	      if( d2 > 2. ) then
		write(6,*) 'distance too large: ',ix,iy,iix,iiy,d,d2
		berror = .true.
	      end if
	      !w = 2. - d			!weight - could be gaussian
	      w = exp(-d2/2.)			!sigma is 1
	      fmextra(2+j,k) = w
	      fmweight(iix,iiy) = fmweight(iix,iiy) + w
	    !end if
	  end do
	end do

	if( berror ) then
	  stop 'error stop fm_extra_setup: internal error (3)'
	end if

!	---------------------------------------------------------
!	scale weight to 1
!	---------------------------------------------------------

	do k=1,nkn
	  ix = nint(fmextra(1,k))
	  iy = nint(fmextra(2,k))
	  if( ix == 0 .or. iy == 0 ) cycle
	  do j=1,4
	    iix = ix+jx(j)
	    iiy = iy+jy(j)
	    !if( fm(4,iix,iiy) == 0 ) then
	      w = fmweight(iix,iiy)
	      if( w > 0. ) fmextra(2+j,k) = fmextra(2+j,k) / w
	    !end if
	  end do
	end do

!	---------------------------------------------------------
!	check if weight sums up to 1
!	---------------------------------------------------------

	fmweight = 0.

	do k=1,nkn
	  ix = nint(fmextra(1,k))
	  iy = nint(fmextra(2,k))
	  if( ix == 0 .or. iy == 0 ) cycle
	  do j=1,4
	    iix = ix+jx(j)
	    iiy = iy+jy(j)
	    !if( fm(4,iix,iiy) == 0 ) then
	      w = fmextra(2+j,k)
	      fmweight(iix,iiy) = fmweight(iix,iiy) + w
	    !end if
	  end do
	end do

	berror = .false.
	do iy=1,ny
	  do ix=1,nx
	    w = fmweight(ix,iy)
	    if( w > 0 ) then
	      if( abs(w-1.) > eps ) then
		berror = .true.
		write(6,*) 'error... ',ix,iy,w
	      end if
	    end if
	  end do
	end do

	if( berror ) then
	  stop 'error stop fm_extra_setup: internal error (2)'
	end if

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end

!************************************************

	subroutine fm_extra_3d(nlvdi,nlv,il,nx,ny,fmextra,femdata,regdata)

! interpolates from fem to reg grid using fmextra structure
!
! interpolation is done only in points that have flag set

	use basin

	implicit none

	integer nlvdi,nlv
	integer il(nkn)
        integer nx,ny			!dimension of regular matrix
	real fmextra(6,nkn)
	real femdata(nlvdi,nkn)
	real regdata(nlvdi,nx,ny)

	integer k,l,lmax
	real flag
	real fem2d(nkn)
	real reg2d(nx,ny)

	call getgeoflag(flag)

!	---------------------------------------------------------
!	make sure fem data below bottom is flag
!	---------------------------------------------------------

	do k=1,nkn
	  lmax = il(k)
	  femdata(lmax+1:nlvdi,k) = flag
	end do

!	---------------------------------------------------------
!	interpolate layer by layer
!	---------------------------------------------------------

	do l=1,nlv
	  fem2d(:) = femdata(l,:)
	  reg2d(:,:) = regdata(l,:,:)
	!write(6,*) l,nx,ny,nx*ny
	!write(6,*) (fem2d(k),k=1,nkn,nkn/20)
	!write(6,*) 'before'
	!write(6,*) reg2d
	  call fm_extra_2d(nx,ny,fmextra,fem2d,reg2d)
	  regdata(l,:,:) = reg2d(:,:)
	!write(6,*) 'after'
	!write(6,*) reg2d
	end do

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end

!************************************************

	subroutine fm_extra_2d(nx,ny,fmextra,femdata,regdata)

! interpolates from fem to reg grid using fmextra structure
!
! interpolation is done only in points that have flag set

	use basin

	implicit none

        integer nx,ny			!dimension of regular matrix
	real fmextra(6,nkn)
	real femdata(nkn)
	real regdata(nx,ny)

	integer ix,iy,iix,iiy
	integer j,k
	real x0,y0,dx,dy
	real eps,flag
	real regval,femval
	double precision d,w

	integer, save :: jx(4) = (/0,1,1,0/)
	integer, save :: jy(4) = (/0,0,1,1/)
	double precision fmweight(nx,ny)
	double precision fmdata(nx,ny)

	eps = 0.01
	fmweight = 0.
	fmdata = 0.

	call getgeoflag(flag)

!	---------------------------------------------------------
!	accumulate on regular grid (only where flag is set)
!	---------------------------------------------------------

	do k=1,nkn
	  ix = nint(fmextra(1,k))
	  iy = nint(fmextra(2,k))
	  if( ix == 0 .or. iy == 0 ) cycle
	  femval = femdata(k)
	  if( femval == flag ) cycle
	  do j=1,4
	    iix = ix+jx(j)
	    iiy = iy+jy(j)
	    regval = regdata(iix,iiy)
	    if( regval /= flag ) cycle
	    w = fmextra(2+j,k)
	!write(6,*) ix,iy,iix,iiy,femval,regval,w
	!write(6,*) ix,iy,iix,iiy,femval,regval
	    fmweight(iix,iiy) = fmweight(iix,iiy) + w
	    fmdata(iix,iiy) = fmdata(iix,iiy) + w * femval
	  end do
	end do
	!write(6,*) 'fmweight'
	!write(6,*) fmweight
	!write(6,*) regdata

!	---------------------------------------------------------
!	correct for weight and set where flag
!	---------------------------------------------------------

	where ( fmweight > 0. ) 
	  fmdata = fmdata / fmweight
	else where
	  fmdata = flag
	end where
	!write(6,*) fmdata
	where ( regdata == flag ) regdata = fmdata

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end

!************************************************
!************************************************
!************************************************

	function intri(x,y,xp,yp)

! point in triangle or not
!
! x,y		array of coordinates of vertices of triangle
! xp,yp		coordinates of point
! intri		1: point is in triangle  0: point outside (return value)

	implicit none

! arguments
	integer intri
	real x(3),y(3),xp,yp
! local
	integer k1,k2
	double precision x21,y21,xn,yn
	double precision scal,eps
! save
	save eps
	data eps /1.e-13/

	intri=0

	do k1=1,3
	   k2=mod(k1,3)+1
	   x21=x(k2)-x(k1)
	   y21=y(k2)-y(k1)
	   yn = x21
	   xn = -y21
	   scal=(xp-x(k1))*xn+(yp-y(k1))*yn
	   if(scal.lt.0.) return
	end do

	intri=1	!inside

	end

!************************************************

	function intrid(x,y,xp,yp)

! point in triangle or not (double precision version)
!
! x,y		array of coordinates of vertices of triangle
! xp,yp		coordinates of point
! intri		1: point is in triangle  0: point outside (return value)

	implicit none

! arguments
	integer intrid
	double precision x(3),y(3),xp,yp
! local
	integer k1,k2
	double precision x21,y21,xn,yn
	double precision scal,eps
! save
	save eps
	data eps /1.e-13/

	intrid=0

	do k1=1,3
	   k2=mod(k1,3)+1
	   x21=x(k2)-x(k1)
	   y21=y(k2)-y(k1)
	   yn = x21
	   xn = -y21
	   scal=(xp-x(k1))*xn+(yp-y(k1))*yn
	   if(scal.lt.0.) return
	end do

	intrid=1	!inside

	end

!****************************************************************
!
	function intri0(x,y,xp,yp)
!
! point in triangle or not
!
! x,y		array of coordinates of vertices of triangle
! xp,yp		coordinates of point
! intri		1: point is in triangle  0: point outside (return value)
!
	implicit none
!
! arguments
	integer intri0
	real x(3),y(3),xp,yp
! local
	integer i,k1,k2
	double precision xs,ys,x12,y12
	double precision det,detlam,rlamb
	double precision eps
! save
	save eps
	data eps /1.e-13/
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
	intri0=0
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
	intri0=1	!inside
!
	return
	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine intp_reg_nodes( nx, ny, x0, y0, dx, dy, flag, regval, femval, ierr )

! interpolates regular grid to FEM grid - values are on nodes

	use basin

	implicit none

	integer nx,ny
	real x0,y0,dx,dy
	real flag
	real regval(nx,ny)
	real femval(nkn)	!interpolated values on fem grid (return)
	integer ierr		!error code (return)

	call intp_reg( nx, ny, x0, y0, dx, dy, flag, regval, nkn, xgv, ygv, femval, ierr )

	end

!****************************************************************

	subroutine intp_reg_elems( nx, ny, x0, y0, dx, dy, flag, regval, femval, ierr )

! interpolates regular grid to FEM grid - values are on elements

	use basin

	implicit none

	integer nx,ny
	real x0,y0,dx,dy
	real flag
	real regval(nx,ny)
	real femval(nel)	!interpolated values on fem grid (return)
	integer ierr		!error code (return)

	real xp(nel),yp(nel)

	call bas_get_elem_coordinates(xp,yp)

	call intp_reg( nx, ny, x0, y0, dx, dy, flag, regval, nel, xp, yp, femval, ierr )

	end

!****************************************************************

	subroutine intp_reg_single_nodes( nx, ny, x0, y0, dx, dy, flag, regval, np, nodes, femval, ierr )

! interpolates regular grid to single nodes

	implicit none

	integer nx,ny
	real x0,y0,dx,dy
	real flag
	real regval(nx,ny)
	integer np		!total number of nodes
	integer nodes(np)	!node numbers
	real femval(np)		!interpolated values on nodes (return)
	integer ierr		!error code (return)

	integer nc
	integer nodesc(np),ipg(np)
	real xp(np),yp(np)

	call condense_valid_coordinates(np,nodes,nc,nodesc,ipg)
	call bas_get_special_coordinates(nc,nodesc,xp,yp)

	call intp_reg(nx,ny,x0,y0,dx,dy,flag,regval,nc,xp,yp,femval,ierr)

	call recollocate_nodes(np,nc,ipg,femval)

	end

!****************************************************************

	subroutine condense_valid_coordinates(np,nodes,nc,nodesc,ipg)

	implicit none

	integer np,nc
	integer nodes(np)
	integer nodesc(np)
	integer ipg(np)		!pointer to good nodes

	integer i,k

	nc = 0
	do i=1,np
	  k = nodes(i)
	  if (k == 0 ) cycle
	  nc = nc + 1
	  nodesc(nc) = k
	  ipg(nc) = i
	end do

	end

!****************************************************************

	subroutine recollocate_nodes(np,nc,ipg,femval)

	implicit none

	integer np		!size of original nodes
	integer nc		!size of condensed nodes
	integer ipg(nc)		!pointer to nodes
	real femval(np)		!data (2D)

	integer i
	real femaux(nc)
	real, parameter :: flag = -777.

	if( np == nc ) return	!nothing to do

	femaux(1:nc) = femval(1:nc)
	femval = flag

	do i=1,nc
	  femval(ipg(i)) = femaux(i)
	end do
	
	end

!****************************************************************

	subroutine intp_reg( nx, ny, x0, y0, dx, dy, flag &
     &			, regval, np, xp, yp, femval, ierr )

! interpolation of regular array onto fem grid - general routine
!
! ierr:
!		= 0	no errors
!		< 0	interpolation out of domain (extrapolation)
!		> 0	flag found in interpolation data

	implicit none

	integer nx,ny
	real x0,y0,dx,dy
	real flag
	real regval(nx,ny)
	integer np		!number of fem points
	real xp(np)
	real yp(np)
	real femval(np)		!interpolated values on fem grid (return)
	integer ierr		!error code (return)

	logical bextra,bout,bflag
	logical bintpout,bintpflag,bextend
	logical bdebug
	integer k
	integer imin,jmin
	integer iflag,iout
	real xx,yy,z1,z2,z3,z4,x1,y1,t,u
	real xn,yn
	real zz(4)
 
	real, parameter :: eps = 3.e-4
	!real, parameter :: eps = 0.
	logical outbox
	outbox(t) = ( t-1. > eps .or. t < -eps )

	!bintpout = .false.	!interpolate even if outside
	!bintpout = .true.	!interpolate even if outside
	!bintpflag = .false.	!interpolate even if flag
	!bintpflag = .true.	!interpolate even if flag
	bdebug = .false.	

	call getregextend(bextend)
	bintpout = bextend
	bintpflag = bextend

	iflag = 0	!used flag for interpolation (no data)
	iout = 0	!used outside point for interpolation

	imin = 0
	jmin = 0

	if( dx < 0. .or. dy < 0. ) goto 97

	xn = x0 + (nx-1)*dx
	yn = y0 + (ny-1)*dy

	do k=1,np
	    xx = xp(k)
	    yy = yp(k)
 
	    femval(k) = flag
 
	    if( xx .le. x0 ) then
	      imin = 1
	    else if( xx .ge. xn ) then
	      imin = nx-1
	    else
	      imin=1+(xx-x0)/dx
	    end if
	    if( yy .le. y0 ) then
	      jmin = 1
	    else if( yy .ge. yn ) then
	      jmin = ny-1
	    else
	      jmin=1+(yy-y0)/dy
	    end if

	    if( imin.lt.1 .or. jmin.lt.1 ) goto 99
	    if( imin+1.gt.nx .or. jmin+1.gt.ny ) goto 99

	    x1 = x0+(imin-1)*dx
	    y1 = y0+(jmin-1)*dy
	    t = (xx-x1)/dx
	    u = (yy-y1)/dy

	    bout = .false.
	    if( outbox(t) ) bout = .true.
	    if( outbox(u) ) bout = .true.
	    if( bout ) then
	      if( bintpout ) then
	        if( t .le. 2. ) t = min(1.,t)
	        if( u .le. 2. ) u = min(1.,u)
	        if( t .ge. -1. ) t = max(0.,t)
	        if( u .ge. -1. ) u = max(0.,u)
	      end if
	      bout = .false.
	      if( outbox(t) ) bout = .true.
	      if( outbox(u) ) bout = .true.
	      if( bout ) then
		iout = iout + 1
		if( bdebug ) then
		  write(6,*) 'reg intp: ',t,u
	          call reg_debug_1('out',k,xx,yy)
		end if
		cycle
	      end if
	    end if

	    z1 = regval(imin,jmin)
	    z2 = regval(imin+1,jmin)
	    z3 = regval(imin+1,jmin+1)
	    z4 = regval(imin,jmin+1)
	    zz = (/z1,z2,z3,z4/)

	    if( any(zz == flag) .and. bintpflag ) then
	      call recover_flag(zz,z1,z2,z3,z4,flag)
	    end if

	    bflag = .false.
	    if( z1.eq.flag .or. z2.eq.flag ) bflag = .true.
	    if( z3.eq.flag .or. z4.eq.flag ) bflag = .true.
	    if( bflag ) then
	      iflag = iflag + 1
	      if( bdebug ) call reg_debug_1('flag',k,xx,yy)
	      cycle
	    end if

	    femval(k)=(1-t)*(1-u)*z1+t*(1-u)*z2+t*u*z3+(1-t)*u*z4
	end do
 
	ierr = 0
	if( iout .gt. 0 ) ierr = - iout - iflag
	if( iflag .gt. 0 ) ierr = iflag

	!write(6,*) 'intp_reg: ierr = ',ierr,iout,iflag

	!if( ierr /= 0 ) then
	!  write(6,*) 'reg: ',ierr,np
	!  write(6,*) femval
	!  stop
	!end if

	return
   97	continue
	write(6,*) 'dx,dy: ',dx,dy
	stop 'error stop intp_reg: dx or dy are negative'
   99	continue
	write(6,*) imin,jmin,nx,ny
	stop 'error stop intp_reg: internal error (1)'
	end

!****************************************************************

	subroutine reg_debug_1(what,k,x,y)

	implicit none

	character*(*) what
	integer k
	real x,y

	write(166,*) trim(what),k,x,y

	end

!****************************************************************

	subroutine recover_flag(zz,z1,z2,z3,z4,flag)

	implicit none

	real zz(4)
	real z1,z2,z3,z4
	real flag

	integer i,ic
	real zt

	!write(6,*) 'recovering flag: ',zz

	ic = count(zz == flag)
	if( ic == 4 ) return

	zt = 0.
	do i=1,4
	  if( zz(i) /= flag ) zt = zt + zz(i)
	end do
	zt = zt / (4-ic)
	do i=1,4
	  if( zz(i) == flag ) zz(i) = zt
	end do

	z1 = zz(1)
	z2 = zz(2)
	z3 = zz(3)
	z4 = zz(4)

	!write(6,*) 'recovered flag: ',zz

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine intp_reg_setup_fr(nx,ny,x0,y0,dx,dy,np,xp,yp,fr,ierr)

! interpolation of regular array onto fem grid - general routine
!
! produces array fr that can be used to interpolate

	use basin

	implicit none

	integer nx,ny
	real x0,y0,dx,dy
	integer np		!total number of sparse points
	real xp(np),yp(np)	!coordinates of sparse grid
	real fr(4,np)		!interpolation array on sparse grid (return)
	integer ierr		!number of points with no interpolation

	logical bextend,bout,bdebug
	integer k
	integer imin,jmin
	real xx,yy,x1,y1,t,u
	real xn,yn
	real, parameter :: eps = 1.e-4
	real, parameter :: fr_flag = -999.
 
	bdebug = ( ierr /= 0 )
	ierr = 0
	imin = 0
	jmin = 0
	fr = 0.

	call getregextend(bextend)

	xn = x0 + (nx-1)*dx
	yn = y0 + (ny-1)*dy

	do k=1,np
	    xx = xp(k)
	    yy = yp(k)
 
	    if( xx .le. x0 ) then
	      imin = 1
	    else if( xx .ge. xn ) then
	      imin = nx-1
	    else
	      imin=1+(xx-x0)/dx
	    end if
	    if( yy .le. y0 ) then
	      jmin = 1
	    else if( yy .ge. yn ) then
	      jmin = ny-1
	    else
	      jmin=1+(yy-y0)/dy
	    end if

	    if( imin.lt.1 .or. jmin.lt.1 ) goto 99
	    if( imin+1.gt.nx .or. jmin+1.gt.ny ) goto 99

	    x1 = x0+(imin-1)*dx
	    y1 = y0+(jmin-1)*dy
	    t = (xx-x1)/dx
	    u = (yy-y1)/dy

	    bout = .false.
	    if( t < 0. ) then
	      if( t >= -eps ) t = 0.
	      if( t >= -1. .and. bextend ) t = 0.
	      if( t < 0. ) bout = .true.
	    end if
	    if( u < 0. ) then
	      if( u >= -eps ) u = 0.
	      if( u >= -1. .and. bextend ) u = 0.
	      if( u < 0. ) bout = .true.
	    end if
	    if( t > 1. ) then
	      if( t <= 1.+eps ) t = 1.
	      if( t <= 2. .and. bextend ) t = 1.
	      if( t > 1. ) bout = .true.
	    end if
	    if( u > 1. ) then
	      if( u <= 1.+eps ) u = 1.
	      if( u <= 2. .and. bextend ) u = 1.
	      if( u > 1. ) bout = .true.
	    end if
	    if( bout .and. bdebug ) then
	      call reg_debug_1('out',k,xx,yy)
	      t = fr_flag
	      u = fr_flag
	      ierr = ierr + 1
	    end if

	    fr(1,k) = imin
	    fr(2,k) = jmin
	    fr(3,k) = t
	    fr(4,k) = u
	end do
 
	return
   99	continue
	write(6,*) imin,jmin,nx,ny
	stop 'error stop intp_reg_setup_fr: internal error (1)'
	end

!****************************************************************

	subroutine intp_reg_intp_fr(nx,ny,flag,regval,np,fr,femval,ierr)

! interpolation of regular array onto fem grid - general routine
!
! ierr:
!		= 0	no errors
!		< 0	interpolation out of domain (extrapolation)
!		> 0	flag found in interpolation data

	implicit none

	integer nx,ny
	real flag
	real regval(nx,ny)
	integer np		!number of fem points
	real fr(4,np)
	real femval(np)		!interpolated values on fem grid (return)
	integer ierr		!error code (return)

	logical bdebug,bflag,bextend,bout
	integer k,ks
	integer imin,jmin
	integer iflag_tot,iout_tot
	real, parameter :: fr_flag = -999.
	real z1,z2,z3,z4,t,u
	real zz(4)
 
	bdebug = .true.
	bdebug = ( ierr /= 0 )

        call getregextend(bextend)

	iflag_tot = 0	!used flag for interpolation
	iout_tot = 0	!used outside point for interpolation

	do k=1,np
 
	    femval(k) = flag
 
	    imin = nint(fr(1,k))
	    jmin = nint(fr(2,k))
	    t = fr(3,k)
	    u = fr(4,k)

	    bout = ( t == fr_flag )	!pre-computed
	    if( bout ) then
	      if( bdebug ) then
	        write(6,*) 'debug intp_reg_intp_fr: ',k,u,t
	        call reg_debug_1('out',k,0.,0.)
	      end if
	      iout_tot = iout_tot + 1
	      cycle
	    end if

	    z1 = regval(imin,jmin)
	    z2 = regval(imin+1,jmin)
	    z3 = regval(imin+1,jmin+1)
	    z4 = regval(imin,jmin+1)
            zz = (/z1,z2,z3,z4/)

            if( any(zz == flag) .and. bextend ) then
              call recover_flag(zz,z1,z2,z3,z4,flag)
            end if

	    bflag = any(zz == flag)
            if( bflag .and. bdebug ) then
	      call reg_debug_1('flag',k,real(ierr),0.)
	      iflag_tot = iflag_tot + 1
              cycle
            end if

	    femval(k)=(1-t)*(1-u)*z1+t*(1-u)*z2+t*u*z3+(1-t)*u*z4
	end do
 
	ierr = 0
	if( iout_tot .gt. 0 ) ierr = - iout_tot
	if( iflag_tot .gt. 0 ) ierr = iflag_tot

	return
   99	continue
	write(6,*) imin,jmin,nx,ny
	stop 'error stop intp_reg_intp_fr: internal error (1)'
	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine am2av(am,av,ip,jp)

! compatibility for old calls - from regular to fem

	use basin

	implicit none

	integer ip,jp
	real av(nkn)
	real am(ip,jp)

	integer ierr
	real pxareg,pyareg,pxdreg,pydreg,pzlreg

	call getgeo(pxareg,pyareg,pxdreg,pydreg,pzlreg)

	call intp_reg(ip,jp,pxareg,pyareg,pxdreg,pydreg,pzlreg,am,nkn,xgv,ygv,av,ierr)

	if( ierr /= 0 ) then
	  write(6,*) 'interpolation am2av: ierr = ',ierr
	  !stop 'error stop am2av: error in interpolation'
	end if

	end

!******************************************************
!******************************************************
!******************************************************

        subroutine elemintp(x,y,z,xp,yp,zp)

! interpolation in element (no ev)
!
! interpolates in element given by x,y nodal values z to point xp,yp
! result is in zp
!
! needs no other vectors but needs x,y of nodes

        real x(3),y(3)  !coordinates of nodes
        real z(3)       !values on nodes
        real xp,yp      !coordinates of point
        real zp         !interpolated value (return)

        double precision eps
        parameter ( eps = 1.d-14 )

        integer i,ii,iii
        double precision zh,f,fh
        double precision a(3),b(3),c(3)

        f = 0.
        do i=1,3
           ii=mod(i,3)+1
           iii=mod(ii,3)+1
           a(i)=x(ii)*y(iii)-x(iii)*y(ii)
           b(i)=y(ii)-y(iii)
           c(i)=x(iii)-x(ii)
           f=f+a(i)
        end do
        f = c(3)*b(2) - c(2)*b(3)               !bug_f_64bit
        if( f .le. eps ) goto 99

        zh=0.
        do i=1,3
           fh=(a(i)+xp*b(i)+yp*c(i))/f
           zh=zh+z(i)*fh
        end do

        zp = zh

        return
   99   continue
        write(6,*) 0,f
        write(6,*) x
        write(6,*) y
        write(6,*) a
        write(6,*) b
        write(6,*) c
        stop 'error stop elemintp: area of element'
        end

!******************************************************
!******************************************************
!******************************************************

	subroutine make_reg_box(dreg,regpar)

	use basin

	implicit none

	real dreg
        real regpar(7)

        real xmin,ymin,xmax,ymax
        integer nx,ny
        logical, save :: bdebug = .false.
        real, save :: flag = -999.
	real rround

        call bas_get_minmax(xmin,ymin,xmax,ymax)

        if( bdebug ) write(6,*) xmin,xmax
        xmin = rround(xmin,dreg,-1)
        xmax = rround(xmax,dreg,+1)
        nx = 1 + (xmax-xmin)/dreg
        if( xmin+(nx-1)*dreg < xmax ) nx = nx + 1
        if( bdebug ) write(6,*) xmin,xmax,nx,xmin+(nx-1)*dreg

        if( bdebug ) write(6,*) ymin,ymax
        ymin = rround(ymin,dreg,-1)
        ymax = rround(ymax,dreg,+1)
        ny = 1 + (ymax-ymin)/dreg
        if( ymin+(ny-1)*dreg < ymax ) ny = ny + 1
        if( bdebug ) write(6,*) ymin,ymax,ny,ymin+(ny-1)*dreg

	call setreg(regpar,nx,ny,xmin,ymin,dreg,dreg,flag)

	end

!******************************************************

