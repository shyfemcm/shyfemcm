c
c $Id: subrst.f,v 1.11 2010-03-11 15:36:39 georg Exp $
c
c restart routines
c
c contents :
c
c subroutine inirst		reads and initializes values from restart
c subroutine admrst		administers writing of restart file

c subroutine wrrst(it,iunit)	writes one record of restart data
c subroutine rdrst(itrst,iunit)	reads one record of restart data
c
c subroutine skip_rst(iunit,it,nvers,nrec,nkn,nel,nlv,ierr)
c				returns info on record in restart file
c
c revision log :
c
c 02.09.2004    ggu     started with routine rdrst()
c 18.10.2006    ggu     included hm3v in restart file -> nvers = 4
c 13.06.2008    ggu     new version 5 -> S/T/rho
c 09.01.2009    ggu     bugfix in inirst - file opened with status=new
c 23.03.2009    ggu     bugfix in admrst - itmrst=itanf if itmrst<itanf
c 07.05.2009    ggu     new parameter ityrst
c 29.05.2009    ggu     use closest record for restart (if ityrst=2)
c 13.11.2009    ggu     keep track of restart: /rstrst/ and has_restart()
c 27.11.2009    ggu     deal with empty file, rdrst() restructured
c 19.01.2010    ggu     initialize also conz, has_restart() is function
c 11.03.2010    ggu     write also vertical velocity
c 10.02.2012    ggu     write only last record, restart from last record
c 16.02.2011    aac     write also ecological varibles
c
c*****************************************************************

        subroutine inirst

c reads and initializes values from restart

        implicit none

	integer iokrst,nvers,ibarcl,iconz,iwvert,ieco
	common /rstrst/ iokrst,nvers,ibarcl,iconz,iwvert,ieco
	save /rstrst/

        integer itrst,iunit,ierr,ityrst,it
	double precision dit
        character*80 name

        real getpar
	double precision dgetpar
        integer ifileo

c-----------------------------------------------------------------
c get parameters
c-----------------------------------------------------------------

	iokrst = 0	!if different from 0, restart has been performed

        itrst = nint(getpar('itrst'))
        ityrst = nint(getpar('ityrst'))
        call getfnm('restrt',name)
        if(name.eq.' ') return

c-----------------------------------------------------------------
c name of restart file given -> open and read
c-----------------------------------------------------------------

        write(6,*) '---------------------------------------------'
        write(6,*) '... performing restart from file'
        write(6,*) name
        write(6,*) '---------------------------------------------'

        iunit = ifileo(1,name,'unformatted','old')
        if( iunit .le. 0 ) then
          if( ityrst .le. 0 ) goto 98
          write(6,*) '*** Cannot find restart file ...'
          write(6,*) '*** Continuing with cold start...'
          return
        end if

	it = itrst
        call rdrst(it,iunit,ierr)
        if( ierr .gt. 0 ) then
          if( ityrst .le. 1 ) goto 97
          if( ierr .eq. 95 ) then
            write(6,*) '*** No data in restart file ...'
            write(6,*) '*** Continuing with cold start...'
            return
          end if
          write(6,*) '*** Another time record is used for restart'
          write(6,*) '*** Looking for time = ',itrst
          write(6,*) '*** Finding time = ',it
          write(6,*) '*** Continuing with hot start...'
        end if

	close(iunit)

	if( itrst .eq. -1 ) then	!reset initial time
	  write(6,*) 'setting new initial time: ',it
	  dit = it
	  call dputpar('itanf',dit)
	  write(6,*) 'new itanf: ',nint(dgetpar('itanf'))
	end if

        write(6,*) '---------------------------------------------'
        write(6,*) 'A restart has been performed'
        write(6,*) ' requested restart time = ',itrst
        write(6,*) ' used restart time = ',it
        write(6,*) ' nvers,ibarcl,iconz,ieco ',nvers,ibarcl,iconz,ieco
        write(6,*) ' iwvert ',iwvert
        write(6,*) '---------------------------------------------'

	iokrst = 1

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

        return
   97   continue
        write(6,*) 'no record found for time = ',itrst
        stop 'error stop inirst: Cannot find time record'
   98   continue
        write(6,*) 'no such file : ',name
        stop 'error stop inirst: Cannot read restart file'
        end

c*******************************************************************

	function has_restart(icode)

c gives indication if restart has been performed
c
c icode indicates what information is requested
c
c icode = 0	general restart
c icode = 1	basic restart (hydro values)
c icode = 2	depth values
c icode = 3	t/s/rho values
c icode = 4	conz values
c icode = 5	vertical velocity
c icode = 6	ecological variables

	implicit none

	logical has_restart
	integer icode

	integer iokrst,nvers,ibarcl,iconz,iwvert,ieco
	common /rstrst/ iokrst,nvers,ibarcl,iconz,iwvert,ieco
	save /rstrst/

	if( iokrst .le. 0 ) then
	  has_restart = .false.
	else if( icode .eq. 0 .or. icode .eq. 1 ) then
	  has_restart = .true.
	else if( icode .eq. 2 ) then
	  has_restart = nvers .ge. 4
	else if( icode .eq. 3 ) then
	  has_restart = ibarcl .gt. 0
	else if( icode .eq. 4 ) then
	  has_restart = iconz .gt. 0
	else if( icode .eq. 5 ) then
	  has_restart = iwvert .gt. 0
	else if( icode .eq. 6 ) then
	  has_restart = ieco .gt. 0
	else
	  has_restart = .false.
	end if

	end

c*******************************************************************

        subroutine admrst

c administers writing of restart file

        implicit none

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        character*80 file
        real getpar
        integer ifemop

	logical bonce
	save bonce
        integer idtrst,itmrst,itnext,iunit
        save idtrst,itmrst,itnext,iunit
        integer icall
        save icall
        data icall / 0 /

        if( icall .le. -1 ) return

c-----------------------------------------------------
c initializing
c-----------------------------------------------------

        if( icall .eq. 0 ) then

          idtrst = nint(getpar('idtrst'))
          itmrst = nint(getpar('itmrst'))

          icall = -1
          if( idtrst .eq. 0 ) return
          if( itmrst .gt. itend ) return

	  if( idtrst .lt. 0 ) then	!only last record saved
	    bonce = .true.
	    idtrst = -idtrst
	  else
	    bonce = .false.
	  end if

          icall = 1
          if( idtrst .le. idt ) idtrst = idt
	  if( itmrst .lt. itanf ) itmrst = itanf
          itnext = itmrst
	  if( itmrst .eq. itanf ) itnext = itnext + idtrst

	  if( .not. bonce ) then
            iunit = ifemop('rst','unformatted','new')
            if( iunit .le. 0 ) goto 98
	  end if

        end if

c-----------------------------------------------------
c normal call and writing
c-----------------------------------------------------

        if( it .lt. itnext ) return

	if( bonce ) then
          iunit = ifemop('rst','unformatted','new')
          if( iunit .le. 0 ) goto 98
          call wrrst(it,iunit)
	  close(iunit)
	else
          call wrrst(it,iunit)
	end if

        itnext = itnext + idtrst

c-----------------------------------------------------
c end of routine
c-----------------------------------------------------

        return
   98   continue
        stop 'error stop admrst: Cannot open rst file'
        end

c*******************************************************************

        subroutine wrrst(it,iunit)

c writes one record of restart data

        implicit none

        integer it
        integer iunit

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        real zenv(3,1)
        common /zenv/zenv
        real znv(1)
        common /znv/znv
        integer iwegv(1)
        common /iwegv/iwegv
        real utlnv(nlvdim,1), vtlnv(nlvdim,1)
        common /utlnv/utlnv, /vtlnv/vtlnv
	real hm3v(3,1)
	common /hm3v/hm3v
        real saltv(nlvdim,1),tempv(nlvdim,1),rhov(nlvdim,1)
        common /saltv/saltv, /tempv/tempv, /rhov/rhov
        real wlnv(0:nlvdim,1)
	common /wlnv/wlnv

        real cnv(nlvdim,nkndim)
        common /cnv/cnv
        real conzv(nlvdim,nkndim,ncsdim)
        common /conzv/conzv

        integer ii,l,ie,k,i
	integer ibarcl,iconz,ieco
        integer nvers

	real getpar

        nvers = 8

	ibarcl = nint(getpar('ibarcl'))
	iconz = nint(getpar('iconz'))
	ieco = nint(getpar('ibfm'))

        write(iunit) it,nvers,1
        write(iunit) nkn,nel,nlv

        write(iunit) (iwegv(ie),ie=1,nel)
        write(iunit) (znv(k),k=1,nkn)
        write(iunit) ((zenv(ii,ie),ii=1,3),ie=1,nel)
        write(iunit) ((utlnv(l,ie),l=1,nlv),ie=1,nel)
        write(iunit) ((vtlnv(l,ie),l=1,nlv),ie=1,nel)

        write(iunit) ((hm3v(ii,ie),ii=1,3),ie=1,nel)

        write(iunit) ibarcl
	if( ibarcl .gt. 0 ) then
          write(iunit) ((saltv(l,k),l=1,nlv),k=1,nkn)
          write(iunit) ((tempv(l,k),l=1,nlv),k=1,nkn)
          write(iunit) ((rhov(l,k),l=1,nlv),k=1,nkn)
	end if

        write(iunit) iconz
	if( iconz .eq. 1 ) then
          write(iunit) ((cnv(l,k),l=1,nlv),k=1,nkn)
	else if( iconz .gt. 1 ) then
          write(iunit) (((conzv(l,k,i),l=1,nlv),k=1,nkn),i=1,iconz)
	end if
	
        write(iunit) nlv-1
	if( nlv .gt. 1 ) then
          write(iunit) ((wlnv(l,k),l=0,nlv),k=1,nkn)
	end if

	write(iunit) ieco
	if( ieco .gt. 0 ) then
	  call write_restart_eco(iunit)
        end if

        end

c*******************************************************************

        subroutine rdrst(itrst,iunit,ierr)

c reads restart file until it finds itrst

        implicit none

        integer itrst
        integer iunit
        integer ierr            !error code - different from 0 if error

        integer ii,l,ie,k
        integer itaux
        integer irec
        logical bloop,blast,bnext

        irec = 0
	ierr = 0
        bloop = .true.
	blast = itrst .eq. -1		! take last record

        do while( bloop )
          call rdrst_record(itaux,iunit,ierr)
          if( ierr .eq. 0 ) irec = irec + 1
	  bnext = itaux .lt. itrst .or. blast	!look for more records
          bloop = ierr .eq. 0 .and. bnext
        end do

	if( irec .gt. 0 ) then
	  if( blast ) then
	    ierr = 0
	    itrst = itaux
	    return
          else if( itaux .eq. itrst ) then
	    return
	  end if
	end if

        if( ierr .ne. 0 ) then          !EOF
          if( irec .eq. 0 ) then        !no data found
                goto 95
          else                          !last record in file has smaller time
                goto 97
          end if
        else                            !read past desired time
                goto 97
        end if

        return
   95   continue
        write(6,*) 'reading restart file... '
        write(6,*) 'no records in file'
        ierr = 95
        return
   97   continue
        write(6,*) 'reading restart file... '
        write(6,*) 'last record read at time = ',itaux
        write(6,*) 'no record found for time = ',itrst
        ierr = 97
        return
        end

c*******************************************************************

	subroutine skip_rst(iunit,it,nvers,nrec,nkn,nel,nlv,ierr)

c returns info on record in restart file and skips data records

	implicit none

	integer iunit,it,nvers,nrec,nkn,nel,nlv,ierr
	integer ibarcl,iconz,iwvert,ieco

	read(iunit,end=2,err=3) it,nvers,nrec
        read(iunit) nkn,nel,nlv

	read(iunit)
	read(iunit)
	read(iunit)
	read(iunit)
	read(iunit)
	if( nvers .ge. 4 ) read(iunit)

	if( nvers .ge. 5 ) then
	  read(iunit) ibarcl
	  if( ibarcl .gt. 0 ) then
	    read(iunit)
	    read(iunit)
	    read(iunit)
	  end if
	end if

	if( nvers .ge. 6 ) then
	  read(iunit) iconz
	  if( iconz .gt. 0 ) then
	    read(iunit)
	  end if
	end if

	if( nvers .ge. 7 ) then
	  read(iunit) iwvert
	  if( iwvert .gt. 0 ) then
	    read(iunit)
	  end if
	end if

	if( nvers .ge. 8 ) then
          read(iunit) ieco
          if( ieco .gt. 0 ) then
	    call skip_restart_eco(iunit)
          end if
        end if

	ierr = 0
	return

    2	continue
	ierr = -1
	return
    3	continue
	write(6,*) 'skip_rst: error in reading restart file'
	ierr = 1
	return
	end

c*******************************************************************

        subroutine rdrst_record(it,iunit,ierr)

c reads one record of restart data

        implicit none

        integer it
        integer iunit
        integer ierr            !error code - different from 0 if error

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        real zenv(3,1)
        common /zenv/zenv
        real znv(1)
        common /znv/znv
        integer iwegv(1)
        common /iwegv/iwegv
        real utlnv(nlvdim,1), vtlnv(nlvdim,1)
        common /utlnv/utlnv, /vtlnv/vtlnv
        real hm3v(3,1)
        common /hm3v/hm3v
        real saltv(nlvdim,1),tempv(nlvdim,1),rhov(nlvdim,1)
        common /saltv/saltv, /tempv/tempv, /rhov/rhov
        real wlnv(0:nlvdim,1)
	common /wlnv/wlnv

        real cnv(nlvdim,nkndim)
        common /cnv/cnv
        real conzv(nlvdim,nkndim,ncsdim)
        common /conzv/conzv

	integer iokrst,nvers,ibarcl,iconz,iwvert,ieco
	common /rstrst/ iokrst,nvers,ibarcl,iconz,iwvert,ieco
	save /rstrst/

        integer ii,l,ie,k,i
        integer nversaux,nrec
        integer nknaux,nelaux,nlvaux

        ierr = 0
	ibarcl = 0
	iconz = 0
	ieco = 0

          read(iunit,end=97) it,nvers,nrec
          if( nvers .lt. 3 ) goto 98

          read(iunit) nknaux,nelaux,nlvaux
          if( nknaux .ne. nkn ) goto 99
          if( nelaux .ne. nel ) goto 99
          if( nlvaux .ne. nlv ) goto 99

          read(iunit) (iwegv(ie),ie=1,nel)
          read(iunit) (znv(k),k=1,nkn)
          read(iunit) ((zenv(ii,ie),ii=1,3),ie=1,nel)
          read(iunit) ((utlnv(l,ie),l=1,nlv),ie=1,nel)
          read(iunit) ((vtlnv(l,ie),l=1,nlv),ie=1,nel)

          if( nvers .ge. 4 ) then
            read(iunit) ((hm3v(ii,ie),ii=1,3),ie=1,nel)
          end if

          if( nvers .ge. 5 ) then
            read(iunit) ibarcl
            if( ibarcl .gt. 0 ) then
              read(iunit) ((saltv(l,k),l=1,nlv),k=1,nkn)
              read(iunit) ((tempv(l,k),l=1,nlv),k=1,nkn)
              read(iunit) ((rhov(l,k),l=1,nlv),k=1,nkn)
            end if
          end if

          if( nvers .ge. 6 ) then
            read(iunit) iconz
	    if( iconz .eq. 1 ) then
              read(iunit) ((cnv(l,k),l=1,nlv),k=1,nkn)
	    else if( iconz .gt. 1 ) then
              read(iunit) (((conzv(l,k,i),l=1,nlv),k=1,nkn),i=1,iconz)
	    end if
	  end if

	  if( nvers .ge. 7 ) then
	    read(iunit) iwvert
	    if( iwvert .gt. 0 ) then
              read(iunit) ((wlnv(l,k),l=0,nlv),k=1,nkn)
	    end if
	  end if

	  if( nvers .ge. 8 ) then
	    read(iunit) ieco
            if( ieco .gt. 1 ) then
	      call read_restart_eco(iunit)
	    end if
          end if

        return
   97   continue
        ierr = -1
        return
   98   continue
        write(6,*) 'reading restart file...'
        write(6,*) 'nvers: ',nvers
        stop 'error stop rdrst: cannot read this version'
   99   continue
        write(6,*) 'reading restart file...'
        write(6,*) nkn,nel,nlv,nlvdim
        write(6,*) nknaux,nelaux,nlvaux
        stop 'error stop rdrst: dimensions'
        end

c*******************************************************************

