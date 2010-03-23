c
c $Id: subnosa.f,v 1.2 2006/09/28 08:56:59 georg Exp $
c
c auxiliar nos routines
c
c contents :
c
c wrnos2d(name,title,value)             write 2d nos file
c
c revision log :
c
c 05.01.2005    ggu     new routine wrnos2d() to write 2d nos file
c 28.09.2006    ggu     new routines wrnos2d_it, wrnos2d_index, extract_level
c
c***************************************************************************

        subroutine wrnos2d(name,title,value)

c write 2d nos file

        implicit none

        character*(*) name,title
        real value(1)

	integer it

	it = 0
	call wrnos2d_it(it,name,title,value)

	end

c***************************************************************************

        subroutine wrnos2d_it(it,name,title,value)

c write 2d nos file

        implicit none

	integer it
        character*(*) name,title
        real value(1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real hev(1)
        common /hev/hev

        character*80 pfile
        character*80 ptitle
        integer nvers,nvar,ivar,nlv
        integer nb,ierr
        integer ie,ii
        integer ilhkv(1)
        real hlv(1)
        real h

        integer ifileo

c-----------------------------------------------------------------
c initialize variables
c-----------------------------------------------------------------

        nvers = 3
        nvar = 1
        ivar = 99
        nlv = 1

        ptitle = title

c-----------------------------------------------------------------
c writing file
c-----------------------------------------------------------------

        pfile = name
        call filext(pfile,'.nos')
        write(6,*) 'writing file ',pfile(1:50)
        nb = ifileo(55,pfile,'unform','unknown')
        if( nb .le. 0 ) goto 98

        call wfnos(nb,nvers,nkn,nel,nlv,nvar,ptitle,ierr)
        if(ierr.ne.0) goto 99
        call wsnos(nb,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 99
        call wrnos(nb,it,ivar,nlv,ilhkv,value,ierr)
        if( ierr .ne. 0 ) goto 99

	call delnos(nb)
        close(nb)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

        return
   98   continue
        write(6,*) pfile
        write(6,*) nb
        stop 'error stop wrnos2d: opening file'
   99   continue
        write(6,*) pfile
        write(6,*) ierr
        stop 'error stop wrnos2d: writing file'
        end

c***************************************************************************

        subroutine wrnos2d_index(it,index,name,title,value)

	implicit none

	integer it
	integer index
        character*(*) name,title
        real value(1)

	integer i
	character*5 number
	character*80 file

	write(number,'(i5)') index
	do i=1,5
	  if( number(i:i) .eq. ' ' ) number(i:i) = '0'
	end do

	file = name//number
        call wrnos2d_it(it,file,title,value)

	end

c***************************************************************************

	subroutine extract_level(nlvdim,nkn,level,v3,v2)

	implicit none

	integer nlvdim
	integer nkn
	integer level
	real v3(nlvdim,1)
	real v2(1)

	integer k

	if( level .gt. nlvdim ) stop 'error stop extract_level: level'

	do k=1,nkn
	  v2(k) = v3(level,k)
	end do

	end

c***************************************************************************
