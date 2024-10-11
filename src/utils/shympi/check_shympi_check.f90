
!*************************************************************

	program check_shympi_check

! simplified version of check_shympi_debug

	use clo

	implicit none

	logical btwo
	integer ierr,iu,iu2,ios
	integer l,ll,lmax
	integer ll2,lmax2
	integer nrec,icall,isact,nsize
	integer nrec2,icall2,isact2,nsize2
	logical belem,belem2
	integer itime
	integer nc
	double precision dtime,dtime2
	character*80 what,what2
	character*80 name_one,name_two
	integer, allocatable :: ies(:,:)
	real, allocatable :: vals(:)
	integer, allocatable :: ies2(:,:)
	real, allocatable :: vals2(:)

	iu = 1
	iu2 = 2
	btwo = .false.
	call check_shympi_check_init

	nc = clo_number_of_files()
	if( nc == 0 ) stop 'no file given'
	if( nc > 2 ) stop 'error stop: not yet ready for 2 files'
	call clo_get_file(1,name_one)
	open(iu,file=name_one,status='old',form='unformatted',iostat=ios)
	if( ios /= 0 ) stop 'error stop check500: cannot open file'
	write(6,*) 'file opened: ',trim(name_one)
	if( nc == 2 ) then
	  btwo = .true.
	  call clo_get_file(2,name_two)
	  open(iu2,file=name_two,status='old',form='unformatted',iostat=ios)
	  if( ios /= 0 ) stop 'error stop check500: cannot open file 2'
	  write(6,*) 'file opened: ',trim(name_two)
	end if

	do
	  call read_time_header(iu,dtime,isact,nsize,lmax,belem,what,ierr)
	  if( ierr /= 0 ) exit
	  if( btwo ) then
	    call read_time_header(iu2,dtime2,isact2,nsize2 &
     &					,lmax2,belem2,what2,ierr)
	    if( ierr /= 0 ) exit
	    if( dtime /= dtime2 ) goto 99
	    if( isact /= isact2 ) goto 99
	    if( nsize /= nsize2 ) goto 99
	    if( lmax /= lmax2 ) goto 99
	    if( lmax /= lmax2 ) goto 99
	    if( belem /= belem2 ) goto 99
	    if( what /= what2 ) goto 99
	  end if
	  write(6,*) dtime,isact,nsize,lmax,trim(what)
	  do l=1,lmax
	    call read_level_header(iu,ll,nsize,nrec,icall)
	    if( l /= ll ) stop 'error stop check500: l/=ll'
	    if( btwo ) then
	      call read_level_header(iu2,ll2,nsize2,nrec2,icall2)
	      if( ll /= ll2 ) goto 98
	      if( nsize /= nsize2 ) goto 98
	      if( nrec /= nrec2 ) goto 98
	      if( icall /= icall2 ) goto 98
	    end if
	    if( .not. btwo ) write(6,*) l,ll,nsize,nrec,icall
	    allocate(ies(2,nsize),vals(nsize))
	    call read_data(iu,nsize,ies,vals)
	    if( btwo ) then
	      allocate(ies2(2,nsize),vals2(nsize))
	      call read_data(iu2,nsize,ies2,vals2)
	    end if
	    if( btwo ) call compare_files(l,nsize,ies,vals,ies2,vals2)
	    deallocate(ies,vals)
	    if( btwo ) deallocate(ies2,vals2)
	  end do
	end do

	if( ierr > 0 ) stop 'error stop check500: read error'

	stop
   98	continue
	write(6,*) 'level records are not comparable:'
	write(6,*) ll,nsize,nrec,icall
	write(6,*) ll2,nsize2,nrec2,icall2
	stop 'error stop 98'
   99	continue
	write(6,*) 'time records are not comparable:'
	write(6,*) dtime,isact,nsize,lmax,belem,trim(what)
	write(6,*) dtime2,isact2,nsize2,lmax2,belem2,trim(what2)
	stop 'error stop 99'
	end

!*************************************************************

	subroutine compare_files(l,nsize,ies,vals,ies2,vals2)

	implicit none

	integer l
	integer nsize
	integer ies(2,nsize)
	real vals(nsize)
	integer ies2(2,nsize)
	real vals2(nsize)

	integer i,ierror
	integer, save :: ierror_max = 10

	ierror = 0

	do i=1,nsize
	  if( any(ies(:,i) /= ies2(:,i) ) ) then
	    ierror = ierror + 1
	    if( ierror > ierror_max ) cycle
	    write(6,*) l,ies(:,i),ies2(:,i)
	  end if
	end do
	if( ierror > 0 ) then
	  write(6,*) ierror,' errors found, max shown ',ierror_max
	  stop 'error stop: error in ies'
	end if

	do i=1,nsize
	  if( vals(i) /= vals2(i) ) then
	    ierror = ierror + 1
	    if( ierror > ierror_max ) cycle
	    write(6,*) ierror,i,l,ies(:,i)
	    write(6,*) vals(i),vals2(i)
	  end if
	end do
	if( ierror > 0 ) then
	  write(6,*) ierror,' errors found, max shown ',ierror_max
	  stop 'error stop: error in vals'
	end if

	end

!*************************************************************

	subroutine read_time_header(iu,dtime,isact,nsize,lmax,belem,what,ierr)

	implicit none

	integer ierr,iu
	integer lmax
	integer isact
	integer nsize
	logical belem
	double precision dtime
	character*80 what

	integer ios

	what = 'unknown'
	read(iu,iostat=ios) dtime,isact,nsize,lmax,belem,what
	ierr = ios

	end

!*************************************************************

	subroutine read_level_header(iu,ll,nsize,nrec,icall)

	implicit none

	integer iu
	integer ll,nsize,nrec,icall

	read(iu) ll,nsize,nrec,icall

	end

!*************************************************************

        subroutine check_shympi_check_init

        use clo

        implicit none

        logical baux
        character*80 version

        version = '2.0'

        call clo_init('check_shympi_check','dbg-file(s)',trim(version))

        call clo_add_info('checks files fort.501')

        call clo_parse_options

	end

!*************************************************************

	subroutine read_data(iu,nsize,ies,vals)

        implicit none

	integer iu
	integer nsize
	integer ies(2,nsize)
	real vals(nsize)

	integer i

	read(iu) ies(2,:)
	read(iu) vals(:)

	do i=1,nsize
	  ies(1,i) = i
	end do

	end

!*************************************************************

