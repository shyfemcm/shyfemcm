
!--------------------------------------------------------------------------
!
!    Copyright (C) 2009-2014,2019  Georg Umgiesser
!    Copyright (C) 2014  Petras Zemlys
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

! regular file time handling
!
! revision log :
!
! 10.03.2009	ggu	finished coding
! 27.03.2009	ggu	bug fix ITACT (itact was not adjourned)
! 31.03.2009	ggu	changed header of regular files (x0,.., description)
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2011	ggu	completely restructured - store geo information
! 18.02.2011	ggu	bug fix in rgf_check_header() -> get instead of set
! 22.11.2011	ggu	changed VERS_6_1_37
! 23.02.2012	ggu	new routines for time series (ts) and regular check
! 16.03.2012	ggu	changed VERS_6_1_48
! 02.05.2013	ggu	bug fix: call to rgf_intp() was wrong
! 15.05.2014	ggu&pzy	bug fix: compare reals in rgf_check_header()
! 16.02.2019	ggu	changed VERS_7_5_60
!
! notes :
!
! ifidim	dimension of array (bigger than needed)
! ifimax	elements used for integer information
! ifi_geobase	base from where geo vars start
!
! structure of ifile:
!
!	basic variables

!	  ifile(1)   iunit
!	  ifile(2)   itold
!	  ifile(3)   itnew
!	  ifile(4)   itact
!	  ifile(5)   nvar
!	  ifile(6)   nx
!	  ifile(7)   ny
!	  ifile(8)   n (=nvar*nx*ny)
!
!	from here new variables...
!
!	  ifile(9)	mode (file type)
!			 0 if mode still has to be determined
!			 1 unformatted
!			 2 direct
!			 3 regular
!			 4 time series
!			-1 file is closed or is not existing
!	  ifile(10)	nkn
!
!	at the end geo variables
!
!	  ifile(19)  x0			ifi_geobase = 18
!	  ifile(20)  y0
!	  ifile(21)  dx
!	  ifile(22)  dy
!	  ifile(23)  flag
! 
!	mode: 0 unknown  1 unformatted  2 regular  3 time series  4 direct
!
! in dfile are three records:
!	1	actual timestep
!	2	old timestep
!	3	new timestep
!
! every record may contain more than one variable
!
!	dile(ix,iy,ivar)
!	first data point is lower left corner, then rowwise
!	after this new variable
!
! format of regular meteo forcing file (one time record):
!
! (data is rowwise, lower left corner is first point)
!
!	read(iunit,*,end=1) it,nvar,nx,ny,x0,y0,dx,dy,flag
!	do ivar=1,nvar
!	  read(iunit,'(a)') description
!	  read(iunit,*) ((data(ix,iy,ivar),ix=1,nx),iy=1,ny)
!	end do
!
!*********************************************************************

	subroutine rgf_admin(it,file,nvar0,ndim,ifile,dfile)

! administers interpolation of regular fields

! can be called any time
! at first call ifile(1) must be -1 (to indicate need for initialization)
! ifile and dfile are returned and contain all needed information

	implicit none

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

	integer it		!time for interpolation
	character*(*) file	!name of file (only needed for first call)
	integer nvar0		!how many variables are expected
	integer ndim		!dimension of dfile
	integer ifile(ifidim)	!info on regular field data
	real dfile(ndim,3)	!data of regular field
	
	integer iunit
	integer itact,itold,itnew
	integer i
	integer nvar
	integer nx,ny,n

	integer ifileo

	itact = it

!-------------------------------------------------------------
! initialize regular data files
!-------------------------------------------------------------

	iunit = ifile(1)
	if( iunit .eq. 0 ) return	!file has encountered EOF before

	if( iunit .lt. 0 ) then

	  do i=1,ifidim
	    ifile(i) = 0
	  end do

	  if( file .eq. ' ' ) return	!no file given

	  write(6,*) 'initialization of regular file: ',file

	  iunit = 0
	  iunit = ifileo(0,file,'form','old')
	  if( iunit .le. 0 ) goto 99

	  iunit = -iunit	!force initialization

	  call rgf_intp(iunit,ndim,itact,itold,itnew,ifile &
     &				,dfile(1,1),dfile(1,2),dfile(1,3))

	  nvar  = ifile(5)
	  nx    = ifile(6)
	  ny    = ifile(7)
	  n     = ifile(8)
	  if( nvar .ne. nvar0 ) goto 97

	  write(6,*) 'regular file opened :',iunit,nvar,nx,ny

	  ifile(1) = iunit
	  ifile(2) = itold
	  ifile(3) = itnew
	  ifile(4) = itact
	end if

!-------------------------------------------------------------
! interpolate in time
!-------------------------------------------------------------

	itold = ifile(2)
	itnew = ifile(3)
	itact = ifile(4)

	if( it .ne. itact ) then
	  itact = it			!bugfix ITACT
	  call rgf_intp(iunit,ndim,itact,itold,itnew,ifile &
     &				,dfile(1,1),dfile(1,2),dfile(1,3))
	end if

	ifile(1) = iunit
	ifile(2) = itold
	ifile(3) = itnew
	ifile(4) = itact

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	return
   97	continue
	write(6,*) nvar,nvar0
	stop 'error stop rgf_admin: number of variables'
   99	continue
	write(6,*) iunit,file
	stop 'error stop rgf_admin: cannot open file'
	end

!*********************************************************************

	subroutine rgf_intp(iunit,ndim,it,itold,itnew &
     &                          ,ifile,dact,dold,dnew)

! handles time interpolation of regular fields
!
! file must alread be open on unit iunit
! the first time must be called with negative iunit
! sets iunit = 0 if EOF has been found

	implicit none

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

	integer iunit
	integer ndim
	integer it,itold,itnew
	integer ifile(ifidim)
	real dact(1)
	real dold(1)
	real dnew(1)

	logical bdata

!--------------------------------------------------------------
! read first time
!--------------------------------------------------------------

	if( iunit .eq. 0 ) return	!EOF found before

	if( iunit .lt. 0 ) then
	  iunit = -iunit
	  call rgf_read(iunit,ndim,itnew,ifile,dnew,bdata)
	  if( .not. bdata ) goto 99

	  itold = itnew
	  call rgf_copy(ifile,dnew,dold)
	end if

!--------------------------------------------------------------
! read new record
!--------------------------------------------------------------

	do while( it .gt. itnew )
	  itold = itnew
	  call rgf_copy(ifile,dnew,dold)
	  call rgf_read(iunit,ndim,itnew,ifile,dnew,bdata)
	  if( .not. bdata ) then
	    itnew = it
	    iunit = 0
	  end if
	end do

!--------------------------------------------------------------
! interpolate
!--------------------------------------------------------------

	call rgf_time_interpolate(it,itold,itnew &
     &				,ifile,dact,dold,dnew)

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	return
   99	continue
	stop 'error stop rgf_intp: no data'
	end

!*********************************************************************

	subroutine rgf_time_interpolate(it,itold,itnew &
     &				,ifile,dact,dold,dnew)

! time interpolation of regular field data

	implicit none

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

	integer it,itold,itnew
	integer ifile(ifidim)
	real dact(1)
	real dold(1)
	real dnew(1)

	integer i,n
	real rit,do,dn

        if( itnew .gt. itold ) then
          rit=float(it-itold)/float(itnew-itold)
        else
          rit = 0.
        end if

	n = ifile(8)

	do i=1,n
	  do = dold(i)
	  dn = dnew(i)
	  dact(i) = rit*(dn-do) + do
	end do

	end

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine rgf_copy(ifile,dsource,dtarget)

! copy regular field record

	implicit none

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

	integer ifile(ifidim)
	real dsource(1)
	real dtarget(1)

	integer i,n

	n = ifile(8)

	do i=1,n
	  dtarget(i) = dsource(i)
	end do

	end

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine rgf_read(iunit,ndim,it,ifile,data,bdata)

! read complete record of regular field

	implicit none

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

	integer iunit
	integer ndim
	integer it
	integer ifile(ifidim)
	real data(ndim)
	logical bdata

	integer n,i

!------------------------------------------------------------
! read header
!------------------------------------------------------------

	call rgf_read_header(iunit,it,ifile,bdata)

!------------------------------------------------------------
! check for errors
!------------------------------------------------------------

	if( .not. bdata ) return		!no more data

	n = ifile(8)
	if( n .gt. ndim ) then
	  write(6,*) 'nvar,nx,ny,n:  ',(ifile(i),i=5,8)
	  write(6,*) 'ndim (provided): ',ndim
	  write(6,*) 'ndim (needed): ',n
	  stop 'error stop rgf_read: ndim'
	end if

!------------------------------------------------------------
! read data
!------------------------------------------------------------

	call rgf_read_data(iunit,ifile,data)

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	end

!*********************************************************************

	subroutine rgf_read_header(iunit,it,ifile,bdata)

! reads header of regular field

	implicit none

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

	integer iunit
	integer it
	integer ifile(ifidim)
	logical bdata

	integer nvar
	integer n,nx,ny
	integer i
	real x0,y0,dx,dy,flag

	bdata = .true.
        n = ifile(8)	!total number of data expected - if 0 -> initialize

	read(iunit,*,end=1) it,nvar,nx,ny,x0,y0,dx,dy,flag

	if( n .le. 0 ) then	!initialize
	  call rgf_set_header(ifile,nvar,nx,ny,x0,y0,dx,dy,flag)
	else
	  call rgf_check_header(ifile,it,nvar,nx,ny,x0,y0,dx,dy,flag)
	end if

        n = ifile(8)	!total number of data read
	if( n .le. 0 ) goto 98

	return
    1	continue
	bdata = .false.
	return
   98	continue
	write(6,*) 'nvar,nx,ny,n: ',(ifile(i),i=5,8)
	stop 'error stop rgf_intp: parameters are zero'
	end

!*********************************************************************

	subroutine rgf_read_data(iunit,ifile,data)

! reads data of regular field

	implicit none

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

	integer iunit
	integer ifile(ifidim)
	real data(1)

	integer nvar,nx,ny,nxy
	integer ivar,j,ip
	character*80 description

        nvar  = ifile(5)
        nx    = ifile(6)
        ny    = ifile(7)
	nxy = nx * ny

	do ivar=1,nvar
	  ip = (ivar-1) * nxy
	  read(iunit,'(a)') description
	  read(iunit,*) (data(ip+j),j=1,nxy)
	end do

	end

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine rgf_set_header(ifile,nvar,nx,ny,x0,y0,dx,dy,flag)

! sets header information

	implicit none

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

	integer ifile(ifidim)
	integer nvar
	integer nx,ny
	real x0,y0,dx,dy,flag

	ifile(5) = nvar
	ifile(6) = nx
	ifile(7) = ny
	ifile(8) = nx*ny*nvar

	call rgf_set_geo(ifile,x0,y0,dx,dy,flag)

	end

!*********************************************************************

	subroutine rgf_check_header(ifile,it,nvar,nx,ny,x0,y0,dx,dy,flag)

! checks header information

	implicit none

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

	integer ifile(ifidim)
	integer it
	integer nvar
	integer nx,ny
	real x0,y0,dx,dy,flag

	integer i
	real x0_a,y0_a,dx_a,dy_a,flag_a
	real eps

	eps = 1.e-5

	call rgf_get_geo(ifile,x0_a,y0_a,dx_a,dy_a,flag_a)

	if( ifile(5) .ne. nvar ) goto 99
	if( ifile(6) .ne. nx ) goto 99
	if( ifile(7) .ne. ny ) goto 99
	if( ifile(8) .ne. nx*ny*nvar ) goto 99

        if( abs(x0-x0_a) .gt. eps ) goto 99
        if( abs(y0-y0_a) .gt. eps ) goto 99
        if( abs(dx-dx_a) .gt. eps ) goto 99
        if( abs(dy-dy_a) .gt. eps ) goto 99
        if( abs(flag-flag_a) .gt. eps ) goto 99

	return
   99	continue
	write(6,*) nvar,nx,ny,nx*ny*nvar
	write(6,*) (ifile(i),i=5,8)
	write(6,*) x0,y0,dx,dy,flag
	write(6,*) x0_a,y0_a,dx_a,dy_a,flag_a
	write(6,*) 'some of the above parameters are differing'
	write(6,*) 'file opened at unit: ',ifile(1)
	write(6,*) 'time where inconsistency is happening: ',it
	stop 'error stop rgf_check_header: parameter mismatch'
	end

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine rgf_set_geo(ifile,x0,y0,dx,dy,flag)

! set geo information

	implicit none

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

	integer ifile(ifidim)
	real x0,y0,dx,dy,flag

	integer i
	integer ifi_geobase
	integer ia(5)
	real ra(5)
	equivalence(ra(1),ia(1))

	ifi_geobase = 18

	ra(1) = x0
	ra(2) = y0
	ra(3) = dx
	ra(4) = dy
	ra(5) = flag

	do i=1,5
	  ifile(ifi_geobase+i) = ia(i)
	end do

	end

!*********************************************************************

	subroutine rgf_get_geo(ifile,x0,y0,dx,dy,flag)

! get geo information

	implicit none

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

	integer ifile(ifidim)
	real x0,y0,dx,dy,flag

	integer i
	integer ifi_geobase
	integer ia(5)
	real ra(5)
	equivalence(ra(1),ia(1))

	ifi_geobase = 18

	do i=1,5
	  ia(i) = ifile(ifi_geobase+i)
	end do

	x0 = ra(1)
	y0 = ra(2)
	dx = ra(3)
	dy = ra(4)
	flag = ra(5)

	end

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine rgf_info(nout,ifile)

! prints info of regular field to unit nout

	implicit none

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

	integer nout
	integer ifile(ifidim)

	integer unit,i
	real x0,y0,dx,dy,flag

	unit = nout
	if( nout .le. 0 ) unit = 6

	call rgf_get_geo(ifile,x0,y0,dx,dy,flag)

	write(unit,*) 'rgf_info: ',(ifile(i),i=1,ifimax)
	write(unit,*) 'x0,y0,dx,dy,flag: ',x0,y0,dx,dy,flag

	end

!*********************************************************************

	subroutine rgf_print(nout,ifile,dfile)

! prints data of regular field to unit nout

	implicit none

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

	integer nout
	integer ifile(ifidim)
	real dfile(1)

	integer unit,i
	integer nvar,nx,ny,n
	real x0,y0,dx,dy,flag

	unit = nout
	if( nout .le. 0 ) unit = 6

        n     = ifile(8)
	call rgf_get_geo(ifile,x0,y0,dx,dy,flag)

	write(unit,*) 'rgf_print: ',(ifile(i),i=1,ifimax)
	write(unit,*) 'x0,y0,dx,dy,flag: ',x0,y0,dx,dy,flag
	write(unit,*) (dfile(i),i=1,n)

	end

!*********************************************************************

	subroutine rgf_get_val(ifile,dfile,ivar,i,j,value)

! gets value of regular field

	implicit none

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

	integer ifile(ifidim)
	real dfile(1)
	integer ivar		!variable to get
	integer i,j		!coordinates of variable in grid
	real value		!value of point (i,j) (return)

	integer nvar,nx,ny,n,ip

        nvar  = ifile(5)
        nx    = ifile(6)
        ny    = ifile(7)
        n = nx * ny

	if( i .le. 0 .or. i .gt. nx ) goto 99
	if( j .le. 0 .or. j .gt. ny ) goto 99
	if( ivar .le. 0 .or. ivar .gt. nvar ) goto 99

	ip = (ivar-1)*n + (j-1)*nx + i
	value = dfile(ip)

	return
   99	continue
	write(6,*) 'nvar,nx,ny: ',nvar,nx,ny
	write(6,*) 'ivar,i,j:   ',ivar,i,j
	stop 'error stop rgf_get_val: input parameters'
	end

!*********************************************************************
!*********************************************************************
!*********************************************************************

	function is_meteo_ts(file,nvar_expected)

	implicit none

	logical is_meteo_ts
        character*(*) file
        integer nvar_expected

	integer it,nvar

	call ts_test_file(file,it,nvar_expected,nvar)

        if( nvar_expected .eq. 0 ) then
          is_meteo_ts = nvar .gt. 0
        else
          is_meteo_ts = nvar .eq. nvar_expected
        end if

	end

!*********************************************************************

	subroutine ts_test_file(file,it,nvar_expected,nvar)

! tests file if it is regular meteo forcing file

	implicit none

        character*(*) file
        integer it,nvar_expected,nvar
        integer ierr

        integer iunit
	integer it1,i
	real dummy

	integer ifileo

	it = 0
	nvar = 0

	iunit = ifileo(0,file,'form','old')
	if( iunit .le. 0 ) return

	read(iunit,*,err=2) it,(dummy,i=1,nvar_expected)
	read(iunit,*,err=2) it1,(dummy,i=1,nvar_expected)

	if( it1 .le. it ) goto 2

	nvar = 1	!we should really be able to check number of variables
	if( nvar_expected .gt. 0 ) nvar = nvar_expected

    2	continue
	close(iunit)

	end

!*********************************************************************

	subroutine ts_read(it,ifile,data,bdata)

! reads one line of time series

	implicit none

        integer it
        integer ifile(9)
        real data(1)
        logical bdata

        integer iunit,nvar,i

	it = 0
	bdata = .false.

	iunit = ifile(1)
	nvar = ifile(5)

	read(iunit,*,end=2) it,(data,i=1,nvar)

	bdata = .true.
    2	continue

	end

!*********************************************************************
!*********************************************************************
!*********************************************************************

	function is_meteo_regular(file,id_expected)

	implicit none

	logical is_meteo_regular
        character*(*) file
        integer id_expected

	integer it,id,nvar,nx,ny

	call rgf_test_file(file,it,id,nvar,nx,ny)

        if( id_expected .eq. 0 ) then
          is_meteo_regular = id .gt. 0
        else
          is_meteo_regular = id .eq. id_expected
        end if

	end

!*********************************************************************

	subroutine rgf_test_file(file,it,id,nvar,nx,ny)

! tests file if it is regular meteo forcing file

	implicit none

        character*(*) file
        integer it,id,nvar,nx,ny
        integer ierr

        integer iunit
	integer it1,nvar1,nx1,ny1
	integer ix,iy,ivar
	real x0,y0,dx,dy,flag
	real dummy
	character*80 description

	integer ifileo

	id = 0
	it = 0
	nvar = 0
	nx = 0
	ny = 0

	iunit = ifileo(0,file,'form','old')
	if( iunit .le. 0 ) return

	read(iunit,*,err=2) it,nvar,nx,ny,x0,y0,dx,dy,flag

	if( nvar .lt. 1 .or. nvar .gt. 4 ) goto 2
	if( nx .lt. 1 .or. nx .gt. 10000 ) goto 2
	if( ny .lt. 1 .or. ny .gt. 10000 ) goto 2

	do ivar=1,nvar
	  read(iunit,'(a)',err=2) description
	  read(iunit,*,err=2) ((dummy,ix=1,nx),iy=1,ny)
	end do

	read(iunit,*,err=2) it1,nvar1,nx1,ny1,x0,y0,dx,dy,flag

	if( nvar .ne. nvar1 ) goto 2
	if( nx .ne. nx1 ) goto 2
	if( ny .ne. ny1 ) goto 2
	if( it1 .le. it ) goto 2

	id = 2000	!this is the regular version

    2	continue
	close(iunit)

	end

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine test_regular

        implicit none

        integer ndim
        parameter (ndim=10000)

	integer ifidim,ifimax
	parameter ( ifidim = 30 , ifimax = 10 )

        character*80 file
        integer iunit,irec,ierr,it,id
        integer n,nvar,i,nx,ny
        real data(ndim)
	integer ifile(ifidim)
	logical bdata
	logical is_meteo_regular

        irec = 0
        iunit = 1
        file = 'meteodat.win'
        file = 'wind.win'
        file = 'heat_skadar_20all.tmp'

	do i=1,ifidim
	  ifile(i) = 0
	end do

	if( is_meteo_regular(file,0) ) then
	  write(6,*) 'file is regular meteo forcing'
	else
	  write(6,*) 'file is not regular meteo forcing'
	end if

	call rgf_test_file(file,it,id,nvar,nx,ny)
        write(6,*) 'file info: ',it,id,nvar,nx,ny
        if( id .eq. 0 ) stop

        open(iunit,file=file,form='formatted',status='old')

        do while(.true.)
          n = 0
	  call rgf_read(iunit,ndim,it,ifile,data,bdata)
          if( .not. bdata ) goto 1
          irec = irec + 1
          write(6,*) irec,it
        end do

    1   continue
	close(iunit)

        end

!***************************************************************
!       program test_regular_main
!       call test_regular
!       end
!***************************************************************

