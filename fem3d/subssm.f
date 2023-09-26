
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999,2010,2015,2019  Georg Umgiesser
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

! utility routines for operations on matrices
!
! contents :
!
! subroutine matinv(a,ip,hv,n,nd)	inversion of matrix (no band matrix)
!
! subroutine matpri(iunit,a,n,nd)	print of matrix (no band matrix)
! subroutine matstru(iunit,a,n,nd)	prints structure of non band matrix
!
! function ldiag(a,b,ip,n,m,nd,md)	diagonalization of linear system a*x=b
! function ltria(a,b,ip,n,m,nd,md)	triangolization of linear system a*x=b
! subroutine lsolvd(a,b,ip,n,m,nd,md)	solves diagonal linear system
! function lsqua(a,b,ip,hv,ah,bh,iph,n,m,nd,md,nmdh)
!				minimum distance to origin of point on plane
!
! revision log :
!
! 04.06.1999	ggu	new routine matstru
! 22.06.1999	ggu	headers in new routine matstru
! 23.03.2010	ggu	changed v6.1.1
! 10.02.2015	ggu	some double precision routines introduced
! 26.02.2015	ggu	changed VERS_7_1_5
! 16.02.2019	ggu	changed VERS_7_5_60
!
!************************************************

	subroutine dmatinv(a,ip,hv,n,nd)

! inversion of square matrix (no band matrix)
!
! a		square matrix
! ip		aux vector for pivot search
! hv		aux vector for resolution of matrix
! n		actual dimension of matrix
! nd		formal dimension of matrix

	implicit none

	double precision eps
	parameter (eps=1.e-30)

	integer n,nd
	double precision a(nd,nd),hv(nd)
	integer ip(nd)

	integer i,j,k,ir
	integer ihi,ipk
	double precision amax,hr

	do j=1,n
	  ip(j)=j
	end do

	do j=1,n

! look for pivot

	amax=abs(a(j,j))
	ir=j
	do i=j+1,n
	  if(abs(a(i,j)).gt.amax) then
	    amax=abs(a(i,j))
	    ir=i
	  end if
	end do
	if(amax.lt.eps) goto 99

! change row

	if(ir.gt.j) then
	  do k=1,n
	    hr=a(j,k)
	    a(j,k)=a(ir,k)
	    a(ir,k)=hr
	  end do
	  ihi=ip(j)
	  ip(j)=ip(ir)
	  ip(ir)=ihi
	end if

! transformation

	hr=1./a(j,j)
	do i=1,n
	  a(i,j)=hr*a(i,j)
	end do
	a(j,j)=hr
	do k=1,n
	  if(k.ne.j) then
	    do i=1,n
	      if(i.ne.j) a(i,k)=a(i,k)-a(i,j)*a(j,k)
	    end do
	    a(j,k)=-hr*a(j,k)
	  end if
	end do

	end do

! change column

	do i=1,n
	  do k=1,n
	    ipk=ip(k)
            hv(ipk)=a(i,k)
	  end do
	  do k=1,n
	    a(i,k)=hv(k)
	  end do
	end do

	return

   99	continue
	write(6,*) 'matrix singular. cannot invert matrix'
	write(6,*) 'i,j,amax(i,j) :',ir,j,amax
	stop 'error stop : matinv'
	end

!**********************************************

	subroutine matinv(a,ip,hv,n,nd)
!
! inversion of square matrix (no band matrix)
!
! a		square matrix
! ip		aux vector for pivot search
! hv		aux vector for resolution of matrix
! n		actual dimension of matrix
! nd		formal dimension of matrix
!
	parameter (eps=1.e-30)
	real a(nd,nd),hv(nd)
	integer ip(nd)
!
	do j=1,n
	  ip(j)=j
	end do

	do j=1,n

! look for pivot

	amax=abs(a(j,j))
	ir=j
	do i=j+1,n
	  if(abs(a(i,j)).gt.amax) then
	    amax=abs(a(i,j))
	    ir=i
	  end if
	end do
	if(amax.lt.eps) goto 99

! change row

	if(ir.gt.j) then
	  do k=1,n
	    hr=a(j,k)
	    a(j,k)=a(ir,k)
	    a(ir,k)=hr
	  end do
	  ihi=ip(j)
	  ip(j)=ip(ir)
	  ip(ir)=ihi
	end if

! transformation

	hr=1./a(j,j)
	do i=1,n
	  a(i,j)=hr*a(i,j)
	end do
	a(j,j)=hr
	do k=1,n
	  if(k.ne.j) then
	    do i=1,n
	      if(i.ne.j) a(i,k)=a(i,k)-a(i,j)*a(j,k)
	    end do
	    a(j,k)=-hr*a(j,k)
	  end if
	end do

	end do

! change column

	do i=1,n
	  do k=1,n
	    ipk=ip(k)
            hv(ipk)=a(i,k)
	  end do
	  do k=1,n
	    a(i,k)=hv(k)
	  end do
	end do

	return

   99	continue
	write(6,*) 'matrix singular. cannot invert matrix'
	write(6,*) 'i,j,amax(i,j) :',ir,j,amax
	stop 'error stop : matinv'

	end

!**********************************************

	subroutine dmatnorm(anorm,a,n,nd)

! computes norm of square matrix (column sum norm)

	implicit none

	double precision anorm		!computed norm (return)
	double precision a(nd,nd)	!matrix
	integer n			!size of matrix
	integer nd			!formal dimension of matrix
	
	integer i,j
	double precision acu

	anorm = 0.

	do j=1,n
	  acu = 0.
	  do i=1,n
	    acu = acu + abs(a(i,j))
	  end do
	  anorm = max(anorm,acu)
	end do

	end

!**********************************************

	subroutine dmatmult(a,b,c,n,nd)

! multiplies two square matrices: c = a * b

	implicit none

	double precision a(nd,nd)	!matrix
	double precision b(nd,nd)	!matrix
	double precision c(nd,nd)	!matrix
	integer n			!size of matrix
	integer nd			!formal dimension of matrix
	
	integer i,j,k
	double precision acu

	do i=1,n
	  do j=1,n
	    acu = 0.
	    do k=1,n
	      acu = acu + a(i,k) * b(k,j)
	    end do
	    c(i,j) = acu
	  end do
	end do

	end

!**********************************************

	subroutine matpri(iunit,a,n,nd)

! prints non band matrix in rectangular form
!
! iunit		output unit number
! a		matrix
! n		actual dimension of matrix
! nd		formal dimension of matrix

	implicit none

	integer iunit
	integer n,nd
	real a(nd,nd)

	integer i,j

	do i=1,n
	  write(iunit,'(1x,10e13.5)')(a(i,j),j=1,n)
	end do

	end

!**********************************************

	subroutine matstru(iunit,a,n,nd)

! prints structure of non band matrix in rectangular form
!
! iunit		output unit number
! a		matrix
! n		actual dimension of matrix
! nd		formal dimension of matrix

	implicit none

	integer iunit
	integer n,nd
	real a(nd,nd)

	integer ndim
	parameter (ndim=120)

	character*120 line,lineh,linen
	integer i,j,nc

	character*10 numbers
	save numbers
	data numbers /'0123456789'/

	nc = min(n,ndim)

	do i=1,nc
	  lineh(i:i) = '-'
	  j = mod(i,10) + 1
	  linen(i:i) = numbers(j:j)
	end do

	write(iunit,*)
	write(iunit,'(a,a,a)') '     ',linen(1:nc),' '
	write(iunit,'(a,a,a)') '    +',lineh(1:nc),'+'

	do j=1,n
	  do i=1,nc
	    if( a(i,j) .ne. 0. ) then
	      line(i:i) = '*'
	    else
	      line(i:i) = '.'
	    end if
	  end do
	  write(iunit,'(i3,a,a,a,i3)') j,' |',line(1:nc),'| ',j
	end do

	write(iunit,'(a,a,a)') '    +',lineh(1:nc),'+'
	write(iunit,'(a,a,a)') '     ',linen(1:nc),' '
	write(iunit,*)

	end

!************************************************
!
	function ldiag(a,b,ip,n,m,nd,md)
!
! diagonalization of linear system a*x=b (not necessary square system)
!
! a		matrix of linear system a(n,m)
! b		constant vector
! ip		aux vector for column exchange, contains the
!		...effective variable number for each column on return
! n,m		actual dimension of matrix (n equations, m variables)
! nd,md		formal dimension of matrix
! ldiag		rank of matrix (positive if solvable, negative if not)
!
	integer ldiag
	parameter (eps=1.e-5)
	real a(nd,md),b(nd)
	integer ip(md)
!
!	lingl=0
	ldiag=0		!$$ALPHA
	if(n.le.0.or.m.le.0) return
!
	do j=1,m
	  ip(j)=j
	end do
!
	do j=1,n
! look for pivot in rectangle (j...n,j...m)
	amax=abs(a(j,j))
	jr=j
	jc=j
	do i=j,n
	  do k=j,m
	    if(abs(a(i,k)).gt.amax) then
	      amax=abs(a(i,k))
	      jr=i
	      jc=k
	    end if
	  end do
	end do
	if(amax.lt.eps) goto 99
! change row
	if(jr.le.j) goto 40
	do k=1,m
	  hr=a(j,k)
	  a(j,k)=a(jr,k)
	  a(jr,k)=hr
	end do
	hr=b(j)
	b(j)=b(jr)
	b(jr)=hr
   40	continue
! change column
	if(jc.le.j) goto 41
	do k=1,n
	  hr=a(k,j)
	  a(k,j)=a(k,jc)
	  a(k,jc)=hr
	end do
	ipr=ip(jc)
	ip(jc)=ip(j)
	ip(j)=ipr
   41   continue
! transformation of row j
	hr=1./a(j,j)
	do i=j+1,m
	  a(j,i)=hr*a(j,i)
	end do
	b(j)=hr*b(j)
	a(j,j)=1.
! transformation of rows >j
	do k=j+1,n
	  hr=a(k,j)
	  do i=j+1,m
	    a(k,i)=a(k,i)-hr*a(j,i)
	  end do
	  b(k)=b(k)-hr*b(j)
	  a(k,j)=0.
	end do

	end do
!
	j=n+1
   99	continue
! diagonalization
	do l=j-1,2,-1
	  do k=1,l-1
	    hr=a(k,l)
	    do i=j,m
	      a(k,i)=a(k,i)-hr*a(l,i)
	    end do
	    b(k)=b(k)-hr*b(l)
	    a(k,l)=0.
	  end do
	end do
! look for solvability
	do i=j,n
	  if(abs(b(i)).gt.eps) goto 85
	end do
!
! solvable
	ldiag=j-1
	return
!
   85	continue
! not solvable
	ldiag=-(j-1)
	return
!
	end
!
!************************************************
!
	function ltria(a,b,ip,n,m,nd,md)
!
! triangolization of linear system a*x=b (not necessary square system)
!
! a		matrix of linear system a(n,m)
! b		constant vector b(n)
! ip		aux vector for column exchange, contains the
!		...effective variable number for each column on return
! n,m		actual dimension of matrix (n equations, m variables)
! nd,md		formal dimension of matrix
! ltria		rank of matrix (positive if solvable, negative if not)
!
	integer ltria
	parameter (eps=1.e-5)
	real a(nd,md),b(nd)
	integer ip(md)
!
!	lingl=0
	ltria=0		!$$ALPHA
	if(n.le.0.or.m.le.0) return
!
	do j=1,m
	  ip(j)=j
	end do
!
	do j=1,n
! look for pivot in rectangle (j...n,j...m)
	amax=abs(a(j,j))
	jr=j
	jc=j
	do i=j,n
	  do k=j,m
	    if(abs(a(i,k)).gt.amax) then
	      amax=abs(a(i,k))
	      jr=i
	      jc=k
	    end if
	  end do
	end do
	if(amax.lt.eps) goto 99
! change row
	if(jr.le.j) goto 40
	do k=1,m
	  hr=a(j,k)
	  a(j,k)=a(jr,k)
	  a(jr,k)=hr
	end do
	hr=b(j)
	b(j)=b(jr)
	b(jr)=hr
   40	continue
! change column
	if(jc.le.j) goto 41
	do k=1,n
	  hr=a(k,j)
	  a(k,j)=a(k,jc)
	  a(k,jc)=hr
	end do
	ipr=ip(jc)
	ip(jc)=ip(j)
	ip(j)=ipr
   41   continue
! transformation of row j
	hr=1./a(j,j)
	do i=j+1,m
	  a(j,i)=hr*a(j,i)
	end do
	b(j)=hr*b(j)
	a(j,j)=1.
! transformation of rows >j
	do k=j+1,n
	  hr=a(k,j)
	  do i=j+1,m
	    a(k,i)=a(k,i)-hr*a(j,i)
	  end do
	  b(k)=b(k)-hr*b(j)
	  a(k,j)=0.
	end do

	end do
!
	j=n+1
   99	continue
! look for solvability
	do i=j,n
	  if(abs(b(i)).gt.eps) goto 85
	end do
!
! solvable
	ltria=j-1
	return
!
   85	continue
! not solvable
	ltria=-(j-1)
	return
!
	end
!
!************************************************
!
	subroutine lsolvd(a,b,ip,n,m,nd,md)
!
! solves linear system that is already diagonalizied 
!
! a		matrix of linear system a(n,m)
! b		constant vector for indices [1...n]
!		...variable values for indices [n+1...m]
! ip		aux vector for column exchange, contains the
!		...effective variable number for each column
!		...is used to sort the variables, on return
!		...the variables are in the original order
! n,m		actual dimension of matrix (n equations, m variables)
!		...n is the number of the effective equations
! nd,md		formal dimension of matrix
!
	parameter (eps=1.e-5)
	real a(nd,md),b(nd)
	integer ip(md)

! solve for rest of variables

	do k=1,n
	  hr=b(k)
	  do j=n+1,m
	    hr=hr-a(k,j)*b(j)
	  end do
	  b(k)=hr
	end do

! sorting

	do j=1,m
	  do i=j,m
	    if(ip(i).eq.j) goto 5
	  end do
	  goto 99
    5	  continue
	  hr=b(j)
	  b(j)=b(i)
	  b(i)=hr
	  ipr=ip(j)
	  ip(j)=ip(i)
	  ip(i)=ipr
	end do

	return
   99	continue
	write(6,*) 'error in sorting'
	stop 'error stop : lsolvd'
	end
!
!************************************************
!
	function lsqua(a,b,ip,hv,ah,bh,iph,n,m,nd,md,nmdh)
!
! minimum distance to origin of point on plane
!
! the plane is decribed by n linear equations, there are
! ...in total m variables, the matrix a is diagonalized
! ...the rank is lll (lll =< n) thus only lll equations
! ...are considered
! ...the first lll variables are dependent on the last m-lll,
! ...the least square method is used
! ... the function to minimize is :
!
!	 f = x(1)**2/hv(i) + x(2)**2/hv(2) + ... + x(m)**2/hv(m)
!
! be sure that the dimension    nmdh >= m-lll
!
! a		matrix of linear system a(n,m)
! b		constant vector b(n)
! ip		aux vector for column exchange, contains the
!		...effective variable number for each column ip(m)
! hv		on call vector containing weigths for distance function hv(m)
!		...if one of the hv is 0., weigth 1. is used for all hv
!		...on return contains the solution of the system
! ah		auxiliar matrix to solve for minimum ah(m-lll,m-lll)
! bh		auxiliar vector to solve for minimum bh(m-lll)
! iph		auxiliar vector used for pivot search iph(m-lll)
! n,m		actual dimension of matrix a (n equations, m variables)
! nd,md		formal dimension of matrix a 
! nmdh		formal dimension of matrix ah
! lsqua		rank of solution
!
	parameter (eps=1.e-5)
	real a(nd,md),b(md)
	integer ip(md)
	real ah(nmdh,nmdh),bh(nmdh),hv(nmdh)
	integer iph(nmdh)
!
! set hv
!
	do i=1,m
	if(hv(i).eq.0.) goto 45
	end do
!
   45	continue
!
	if( i .gt. m ) then
	  do i=1,m
	    hv(i)=1.
	  end do
	end if
!
! diagonalize matrix a
!
	lll=ldiag(a,b,ip,n,m,nd,md)
!
	if(lll.lt.0) goto 99
	if(m-lll.gt.nmdh) goto 98
	if(lll.eq.0) goto 97
!
! get distance system (lll equations (n), m variables)
!
	do k=lll+1,m
	  kmat=k-lll
!	  diagonal elements
	  hr=1./hv(k)
	  do i=1,lll
	    hr=hr+a(i,k)*a(i,k)/hv(i)
	  end do
	  ah(kmat,kmat)=hr
!	  out of diagonal elements
	  do j=k+1,m
	    jmat=j-lll
	    hr=0.
	    do i=1,lll
	      hr=hr+a(i,k)*a(i,j)/hv(i)
	    end do
	    ah(kmat,jmat)=hr
	    ah(jmat,kmat)=hr
	  end do
	  hr=0.
	  do i=1,lll
	    hr=hr+b(i)*a(i,k)/hv(i)
	  end do
	  bh(kmat)=hr
	end do
!
! invert matrix
!
	call matinv(ah,iph,hv,m-lll,nmdh)
!
! solve for variables
!
	do k=1,m-lll
	  hr=0.
	  do j=1,m-lll
	    hr=hr+ah(k,j)*bh(j)
	  end do
	  hv(k+lll)=hr
	end do
!
	do i=1,lll
	  hv(i)=b(i)
	end do
!
	call lsolvd(a,hv,ip,lll,m,nd,md)
!
	lsqua=lll
!
	return
   99	continue
	write(6,*) 'system not solvable'
	stop 'error stop : lsqua'
   98	continue
	write(6,*) 'dimension error m-lll,nmdh :',m-lll,nmdh
	stop 'error stop : lsqua'
   97	continue
	do i=1,m
	  b(i)=0.
	end do
	lsqua=0
	return
	end

!************************************************

	subroutine mat_test

	implicit none

	integer i,ndim
	real r

	do i=1,100
	  call random_number(r)
	  ndim = 10.*r + 1
	  call mat_test_n(ndim)
	end do

	end

!************************************************

	subroutine mat_test_n(ndim)

	implicit none

	integer ndim

	integer n,nd
	integer i,j
	real r
	double precision an,aninv,anind,cond
	double precision a(ndim,ndim)
	double precision ainv(ndim,ndim)
	double precision aind(ndim,ndim)
	double precision hv(ndim)
	integer ip(ndim)

	n = ndim
	nd = ndim

	do i=1,n
	  do j=1,n
	    call random_number(r)
	    a(i,j) = 10.*(r-0.5)
	  end do
	end do

	ainv = a
	call dmatinv(ainv,ip,hv,n,nd)

	call dmatnorm(an,a,n,nd)
	call dmatnorm(aninv,ainv,n,nd)
	cond = an * aninv

	call dmatmult(a,ainv,aind,n,nd)
	do i=1,n
	  aind(i,i) = aind(i,i) - 1.
	end do
	call dmatnorm(anind,aind,n,nd)

	write(6,'(i5,3f12.4,e12.4)') ndim,an,aninv,cond,anind

	end

!************************************************

!	program mat_test_main
!	call mat_test
!	end

!************************************************

