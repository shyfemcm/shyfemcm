c
c $Id: adjgrd.f,v 1.7 2004/12/02 09:29:50 georg Exp $
c
c description :
c
c grading routines
c
c contents :
c
c subroutine maxgrd(nkn,nel,nen3v,ngrade,nmax)
c			computes maximum grade (except boundary)
c subroutine statgrd(nkn,nmax,nsg,ngrade)
c			computes statistics on grade
c subroutine mkbound(nkn,nel,ngrdim,nen3v,ngrade,ngri)
c			marks boundary nodes and determines exact grade
c subroutine sortgr(n,iv)
c			compacts array iv (throws out double numbers)
c subroutine compat(n,iv)
c			compacts array iv (throws out double numbers)
c
c function nextgr(k,kgr,ngrdim,ngrade,ngri)
c	               finds next node in grade index
c subroutine exchgr(k,kold,knew,ngrdim,ngrade,ngri)
c			exchanges node in grade index
c subroutine delgr(k,kold,ngrdim,ngrade,ngri)
c	      	       deletes grade in index (kold from index of node k)
c subroutine insgr(k,kafter,knew,ngrdim,ngrade,ngri)
c	               inserts grade in index 
c subroutine prgr(k,ngrdim,ngrade,ngri)
c	               prints grade index
c
c revision log :
c
c 19.05.2003    ggu     some changes in stats
c
c***********************************************************

	subroutine maxgrd(nkn,nel,nen3v,ngrade,nmax)

c computes maximum grade (except boundary -> is one too low)

	implicit none

	integer nkn,nel,nmax
	integer nen3v(3,1)
	integer ngrade(1)

	integer k,ie,ii

	do k=1,nkn
	  ngrade(k) = 0
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ngrade(k) = ngrade(k) + 1
	  end do
	end do

	nmax = 0
	do k=1,nkn
	  if( ngrade(k) .gt. nmax ) nmax = ngrade(k)
	end do

	end

c***********************************************************

	subroutine statgrd(text,nkn,nmax,ngrade,nbound)

c computes statistics on grade

	implicit none

	character*(*) text
	integer nkn,nmax
	integer ngrade(1)
	integer nbound(1)

	integer nsa(0:30)
	integer nsi(0:30)
	integer nsb(0:30)

	integer i,k,n,nb
	integer naux,ndiff
	integer nin,nbn,nto,nis

	do i=0,nmax
	  nsa(i) = 0
	  nsi(i) = 0
	  nsb(i) = 0
	end do

	nin = 0			! total number of internal nodes
	nbn = 0			! total number of boundary nodes
	nto = 0			! total grade

	do k=1,nkn
	  n = ngrade(k)
	  nb = nbound(k)	! > 0 if on boundary
	  nto = nto + n
	  if( n .lt. 0 ) stop 'error stop statgrd: internal error'
	  nsa(n) = nsa(n) + 1
	  if( nb .ne. 0 ) then
	    nsb(n) = nsb(n) + 1
	  else
	    nsi(n) = nsi(n) + 1
	  end if

	  if( nb .ne. 0 ) then
		nbn = nbn + 1
	  else
		nin = nin + 1
	  end if
	end do

c	ntheo = 6*nin + 4*nbn + 6*nis - 6	!theoretical grade

	naux = 6*nin + 4*nbn
	ndiff = nto-naux	! must be divisible by 6, -6 for no island
	nis = ndiff/6 + 1	! total number of islands

	if( mod(ndiff,6) .ne. 0 ) write(6,*) 'error in element layout'

	write(6,'(a,a)') 'statistics on grades: ',text
	write(6,'(5x,a)') '  internal  boundary   islands    total grade'
	write(6,'(5x,3i10,i15)') nin,nbn,nis,nto
	write(6,'(3x,a)') 'grade       all  internal  boundary'
	do i=0,nmax
	  write(6,'(3x,i5,3i10)') i,nsa(i),nsi(i),nsb(i)
	end do

	end

c***********************************************************

	subroutine mkbound(nkn,nel,ngrdim,nen3v,ngrade,nbound,ngri)

c marks boundary nodes and determines exact grade

	implicit none

	integer nkn,nel,ngrdim
	integer nen3v(3,1)
	integer ngrade(1)
	integer nbound(1)
	integer ngri(2*ngrdim,1)

	integer k,ie,ii,n,nb
	integer i1,i2

	write(6,*) 'grading nodes...'

	do k=1,nkn
	  ngrade(k) = 0
	  nbound(k) = 0
	end do

	do ie=1,nel
	  do ii=1,3
	    i1 = mod(ii,3) + 1
	    i2 = mod(i1,3) + 1
	    k  = nen3v(ii,ie)
	    n = ngrade(k)
c	    nodes are inserted in counter-clockwise sense
	    ngri(n+1,k) = nen3v(i1,ie)
	    ngri(n+2,k) = nen3v(i2,ie)
	    ngrade(k) = ngrade(k) + 2
	  end do
	end do

	nb = 0
	do k=1,nkn
	  n = ngrade(k)
	  call sortgr(n,ngri(1,k))
	  if( ngrade(k)/2 .ne. n ) then		!boundary node
	    nb = nb + 1
	    nbound(k) = 1
	  end if
	  ngrade(k) = n
	end do

	write(6,*) 'grading nodes done ... boundary nodes: ',nb

	end

c***********************************************************

	subroutine sortgr(n,iv)

c compacts array iv (throws out double numbers)

	implicit none

	integer n
	integer iv(n)

        integer ndim
        parameter (ndim=60)

	logical bnew
	integer ib,ie,new
	integer i
	integer iaux(ndim)      !aux array

c initialize

	ib = ndim/2
	ie = ib + 1

	iaux(ib) = iv(1)
	iaux(ie) = iv(2)
	iv(1) = 0
	iv(2) = 0

c look forward

	bnew = .true.
	do while( bnew )
	  new = 0
	  do i=1,n,2		!only odd nodes
	    if( iv(i) .eq. iaux(ie) ) new = i
	  end do
	  if( new .eq. 0 ) then
	    bnew = .false.
	  else
	    ie = ie + 1
            if( ie .gt. ndim ) stop 'error stop sortgr: ndim'
	    iaux(ie) = iv(new+1)
	    iv(new) = 0
	    iv(new+1) = 0
	  end if
	end do

c look backward

	bnew = .true.
	do while( bnew )
	  new = 0
	  do i=2,n,2		!only even nodes
	    if( iv(i) .eq. iaux(ib) ) new = i
	  end do
	  if( new .eq. 0 ) then
	    bnew = .false.
	  else
	    ib = ib - 1
            if( ib .lt. 1 ) stop 'error stop sortgr: ndim'
	    iaux(ib) = iv(new-1)
	    iv(new) = 0
	    iv(new-1) = 0
	  end if
	end do

c adjust for internal node

	if( iaux(ie) .eq. iaux(ib) ) ie = ie - 1

c copy to original array

	n = 0
	do i=ib,ie
	  n = n + 1
	  iv(n) = iaux(i)
	end do

	end

c***********************************************************

	subroutine compat(n,iv)

c compacts array iv (throws out double numbers)
c
c not used

	implicit none

	integer n
	integer iv(1)

	integer ivold,i,idiff,nmax

c set to 0 all double nodes

	ivold = iv(1)
	do i=2,n
	  if( ivold .eq. iv(i) ) then
	    iv(i) = 0
	  else
	    ivold = iv(i)
	  end if
	end do

c copy remaining nodes to front of array (leave no holes)

	idiff = 0
	nmax = 0
	do i=1,n
	  if( iv(i) .eq. 0 ) then
	    idiff = idiff + 1
	  else
	    iv(i-idiff) = iv(i)
	    nmax = nmax + 1
	  end if
	end do

c delete useless entries

	do i=nmax+1,n
	  iv(i) = 0
	end do

	n = nmax

	end

c***********************************************************
c***********************************************************
c***********************************************************

	function nextgr(k,kgr,ngrdim,ngrade,ngri)

c finds next node in grade index
c
c in grade index of node k find node after kgr
c does not take into account boundary nodes

	implicit none

	integer nextgr
	integer k,kgr,ngrdim
	integer ngrade(1)
	integer ngri(2*ngrdim,1)

	integer i,n,kn,kext

	nextgr = 0

	n = ngrade(k)

	do i=1,n
	  kn = ngri(i,k)
	  if( kn .eq. kgr) goto 1
	end do

	write(6,*) k,kgr,n,ngrdim
        call nint2ext(k,kext)
        write(6,*) 'int/ext k: ',k,kext
        call nint2ext(kgr,kext)
        write(6,*) 'int/ext kgr: ',kgr,kext
        write(6,*) (ngri(i,k),i=1,n)
	stop 'error stop nextgr: no such node'

    1	continue
	i = i + 1
	if( i .gt. n ) i = 1
	nextgr = ngri(i,k)

	end

c**************************************************************

	subroutine exchgr(k,kold,knew,ngrdim,ngrade,ngri)

c exchanges node in grade index

	implicit none

	integer k,kold,knew,ngrdim
	integer ngrade(1)
	integer ngri(2*ngrdim,1)

	integer n,i

	n = ngrade(k)

	do i=1,n
	  if( ngri(i,k) .eq. kold ) then
		ngri(i,k) = knew
		return
	  end if
	end do

	write(6,*) k,kold,knew
	stop 'error stop exchgr: no such node'

	end

c**************************************************************

	subroutine delgr(k,kold,ngrdim,ngrade,ngri)

c deletes grade in index (kold from index of node k)

	implicit none

	integer k,kold,ngrdim
	integer ngrade(1)
	integer ngri(2*ngrdim,1)

	integer n,i,idiff

	n = ngrade(k)
	ngrade(k) = ngrade(k) - 1

	idiff = 0
	do i=1,n
	  if( idiff .eq. 0 .and. ngri(i,k) .eq. kold ) then
	    idiff = 1
	  else
	    ngri(i-idiff,k) = ngri(i,k)
	  end if
	end do

	if( idiff .eq. 0 ) then
	  write(6,*) k,kold
	  stop 'error stop delgr: no such node'
	end if

	end

c**************************************************************

	subroutine insgrb(k,kbefor,knew,ngrdim,ngrade,ngri)

c inserts grade in index (knew in index of node k befor node kbefor)

	implicit none

	integer k,kbefor,knew,ngrdim
	integer ngrade(1)
	integer ngri(2*ngrdim,1)

	integer n,i

	n = ngrade(k)
	ngrade(k) = ngrade(k) + 1

	do i=n,1,-1			!backwards
	  if( ngri(i,k) .eq. kbefor ) then
	    ngri(i+1,k) = ngri(i,k)
	    ngri(i,k) = knew
	    return
	  end if
	  ngri(i+1,k) = ngri(i,k)
	end do

	write(6,*) k,kbefor,knew
	stop 'error stop insgrb: no such node'

	end 

c**************************************************************

	subroutine insgr(k,kafter,knew,ngrdim,ngrade,ngri)

c inserts grade in index (knew in index of node k after node kafter)

	implicit none

	integer k,kafter,knew,ngrdim
	integer ngrade(1)
	integer ngri(2*ngrdim,1)

	integer n,i

	n = ngrade(k)
	ngrade(k) = ngrade(k) + 1

	do i=n,1,-1			!backwards
	  if( ngri(i,k) .eq. kafter ) then
	    ngri(i+1,k) = knew
	    return
	  end if
	  ngri(i+1,k) = ngri(i,k)
	end do

	write(6,*) 'insgr: ',k,kafter,knew
	write(6,*) 'No such node : ',kafter
	stop 'error stop insgr'

	end 

c**************************************************************

	subroutine prgr(k,ngrdim,ngrade,ngri)

c prints grade index

	implicit none

	integer k,ngrdim
	integer ngrade(1)
	integer ngri(2*ngrdim,1)

	integer n,i

	n = ngrade(k)

	write(6,*) 'grade index : ',k,n
	write(6,*) (ngri(i,k),i=1,n)

	end 

c**************************************************************
