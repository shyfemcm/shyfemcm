c
c $Id: nosextr_records.f,v 1.1 2008-07-16 15:41:39 georg Exp $
c
c 18.11.1998    ggu     check dimensions with dimnos
c
c**********************************************************

	program nosextr_records

c extracts whole records from nos file
c
c records have to be specified on stdin

        implicit none

	include 'param.h'

        integer ndim
	parameter ( ndim = 10000 )

	integer iu(ndim)
	character*80 name,file

        integer nrdim
	parameter ( nrdim = 2000 )

	integer irec(nrdim)

	character*80 title
	real cv(nkndim)
	real cv3(nlvdim,nkndim)

        double precision accum(nlvdim,nkndim)
	real amin(nlvdim,nkndim)
	real amax(nlvdim,nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

        integer nread,ivarold,nextr
        integer l,k,nin,nb
        integer nkn,nel,nlv,nvar
        integer it,ivar
        integer ierr
        integer nvers
        real r,rnull
	real conz,high

        integer iapini,ideffi,ifileo

c-------------------------------------------------------------------
c initialize params
c-------------------------------------------------------------------

	nread=0
	nextr=0
	rnull=0.
        ivarold = 0
	high = 1.e+30

c-------------------------------------------------------------------
c get simulation
c-------------------------------------------------------------------

	if(iapini(2,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

	nin=ideffi('datdir','runnam','.nos','unform','old')
	if(nin.le.0) goto 100

c-------------------------------------------------------------------
c open NOS file and read header
c-------------------------------------------------------------------

        nvers=3
	call rfnos(nin,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title

        call dimnos(nin,nkndim,neldim,nlvdim)

	call rsnos(nin,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

	write(6,*) 'Available levels: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

c-------------------------------------------------------------------
c get records to extract from STDIN
c-------------------------------------------------------------------

	call get_records_from_stdin(nrdim,irec)

c-------------------------------------------------------------------
c open NOS output file
c-------------------------------------------------------------------

	call mkname(' ','extract','.nos',file)
	write(6,*) 'writing file ',file(1:50)
	nb = ifileo(55,file,'unform','new')
	if( nb .le. 0 ) goto 98
	call wfnos(nb,3,nkn,nel,nlv,1,title,ierr)
	if( ierr .ne. 0 ) goto 99
	call wsnos(nb,ilhkv,hlv,hev,ierr)
	if( ierr .ne. 0 ) goto 99

c-------------------------------------------------------------------
c loop on input records
c-------------------------------------------------------------------

  300   continue

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100
        if( ivarold .eq. 0 ) ivarold = ivar
        if( ivar .ne. ivarold ) goto 91

	nread=nread+1
	write(6,*) 'time : ',nread,it,ivar

	if( nread .le. nrdim .and. irec(nread) .ne. 0 ) then
	  call wrnos(nb,it,ivar,nlvdim,ilhkv,cv3,ierr)
	  if( ierr .ne. 0 ) goto 99
	  nextr = nextr + 1
	end if

	goto 300

  100	continue

c-------------------------------------------------------------------
c end of loop
c-------------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*) nextr,' records written to file extract.nos'
	write(6,*)

        if( nextr .le. 0 ) stop 'no file written'

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	stop
   91	continue
	write(6,*) 'file may have only one type of variable'
	write(6,*) 'error ivar : ',ivar,ivarold
	stop 'error stop nosextr_records: ivar'
   98	continue
	write(6,*) 'error opening file'
	stop 'error stop nosextr_records'
   99	continue
	write(6,*) 'error writing file'
	stop 'error stop nosextr_records'
	end

c***************************************************************

        subroutine mimar(xx,n,xmin,xmax,rnull)

c computes min/max of vector
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull		invalid value

        implicit none

        integer n,i,nmin
        real xx(n)
        real xmin,xmax,x,rnull

	do i=1,n
	  if(xx(i).ne.rnull) goto 1
	end do
    1	continue

	if(i.le.n) then
	  xmax=xx(i)
	  xmin=xx(i)
	else
	  xmax=rnull
	  xmin=rnull
	end if

	nmin=i+1

        do i=nmin,n
          x=xx(i)
	  if(x.ne.rnull) then
            if(x.gt.xmax) xmax=x
            if(x.lt.xmin) xmin=x
	  end if
        end do

        end

c***************************************************************

	subroutine get_records_from_stdin(ndim,irec)

c gets records to extract from stdin

	implicit none

	integer ndim
	integer irec(ndim)

	integer i,ir

	do i=1,ndim
	  irec(i) = 0
	end do

	write(6,*) 'Please enter the record numbers to be extracted.'
	write(6,*) 'Enter every record on a single line.'
	write(6,*) 'Finish with 0 on the last line.'
	write(6,*) 'example:'
	write(6,*) '   5'
	write(6,*) '  10'
	write(6,*) '  15'
	write(6,*) '  0'
	write(6,*) ' '

	do while(.true.)
	  write(6,*) 'Enter record to extract (0 to end): '
	  ir = 0
	  read(5,'(i10)') ir

	  if( ir .le. 0 ) then
	    return
	  else if( ir .gt. ndim ) then
	    write(6,*) 'Cannot extract records higher than ',ndim
	    write(6,*) 'Please change ndim and recompile.'
	  else
	    irec(ir) = 1
	  end if
	end do

	end

c***************************************************************
