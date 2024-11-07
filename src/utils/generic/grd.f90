
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998,2001,2005,2009-2010,2012,2012  Georg Umgiesser
!    Copyright (C) 2014-2019  Georg Umgiesser
!    Copyright (C) 2018  Christian Ferrarin
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

! rdgrd routines - read GRD files
!
! contents :
!
!       subroutine rdgrd(...)
!				reads grd file
!       subroutine rdcom(...)
!				reads comment from .grd file
!       subroutine rdnode(...)
!				reads node from .grd file
!       subroutine rdelem(...)
!				reads element from .grd file
!       subroutine rdline(...)
!				reads line from .grd file
!	subroutine read_node_list(...)
!				reads node list
!	subroutine rdunknown(iwhat,berr)
!				handles unknown type
!
!       function ifstch(line)
!				finds first char of line that is not blank
!       subroutine fempar(line)
!				read parameters for fem model
!
!	subroutine ex2in(nkn,ne,nl,ipnv,ipaux,nen3v,inodlv,berr)
!				changing extern with intern node numbers
!	subroutine chex2in(nkn,n,nodes,ipnv,index,berr)
!				changing extern with intern node numbers (list)
!
!	internal routines:
!
!	subroutine grd_internal_init(file)
!	function grd_next_line()
!	subroutine grd_nvals(nvals)
!	subroutine grd_vals(nvals,vals)
!	subroutine grd_ival(ival,val)
!	subroutine grd_line_info(iline_grd,line_grd)
!	subroutine grd_write_line_info
!
!	writing routines
!
!	subroutine grd_write_grid(
!				writes grd file
!	subroutine grd_write_node_list(nout,n,nodes,ipnv,depth)
!				writes out node list
!
! revision log :
!
! 20.03.1998	ggu	declare f(6) in subroutines to avoid compiler warnings
! 20.05.1998	ggu	open file with ifileo()
! 22.05.1998	ggu	bug fix (ifileo not declared)
! 27.03.2001	ggu	assign depth -999. to depth not set
! 09.10.2001	ggu	read from stdin
! 09.10.2001	ggu	read node type (ianv)
! 18.10.2005	ggu	error messages slightly changed
! 22.04.2009	ggu	changes from spline integrated (read lines)
! 24.04.2009	ggu	newly restructured
! 23.03.2010	ggu	changed v6.1.1
! 14.02.2012	ggu	changed VERS_6_1_44
! 09.03.2012	ggu	handle dimension error more gracefully
! 16.03.2012	ggu	changed VERS_6_1_48
! 23.12.2014	ggu	changed VERS_7_0_11
! 08.01.2015	ggu	common blocks in include file
! 05.05.2015	ggu	changed VERS_7_1_10
! 10.07.2015	ggu	changed VERS_7_1_50
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 30.07.2015	ggu	changed VERS_7_1_83
! 02.10.2015	ggu	new routine is_grd_file()
! 02.10.2015	ggu	now stopping on first error
! 10.10.2015	ggu	changed VERS_7_3_2
! 18.12.2015	ggu	changed VERS_7_3_17
! 25.05.2016	ggu	changed VERS_7_5_10
! 06.06.2016	ggu	bstop substituted with berr, new accessor routines
! 14.06.2016	ggu	changed VERS_7_5_14
! 09.09.2016	ggu	changed VERS_7_5_17
! 10.02.2017	ggu	bug fix: do not allocate at least 1 array element
! 14.11.2017	ggu	changed VERS_7_5_36
! 24.01.2018	ggu	changed VERS_7_5_41
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.04.2018	ggu	changed VERS_7_5_43
! 06.07.2018	ggu	new handling of line reading routines: read_all_lines()
! 25.10.2018	ccf	grid output in gr3 and msh formats
! 16.02.2019	ggu	changed VERS_7_5_60
! 28.05.2020	ggu	bgrdwrite and grd_set_write() implemented
! 12.02.2022	ggu	new variable bcdepth
! 09.03.2022	ggu	new routine grd_write_node()
! 08.06.2022	ggu	new calling sequence in grd_write_item()
! 16.06.2022	ggu	new routine write_grd_file() for simplified writing
! 13.12.2022	ggu	new routine write_grd_file_with_depth()
! 18.10.2023	ggu	deleted routines dealing with gr3 and msh files
! 10.05.2024	ggu	handles isphe given in grd file (FEM-SPHER)
! 07.11.2024	ggu	in write_grd_file_with_detail() also write zoom window
!
!**********************************************************

!==============================================================
	module grd
!==============================================================

	implicit none

        real, save :: xscale_grd = 1.
        real, save :: yscale_grd = 1.
        real, save :: zscale_grd = 1.
	character*80, save :: title_grd = ' '
	real, save :: dcor_grd = 0.
	real, save :: dirn_grd = 0.
	integer, save :: sphe_grd = -1

        integer, save :: nin_grd,iline_grd,ianz_grd
        real, save :: f_grd(81)
        character*132, save :: line_grd

	integer, save :: nk_grd = 0
	integer, save :: ne_grd = 0
	integer, save :: nl_grd = 0
	integer, save :: nne_grd = 0
	integer, save :: nnl_grd = 0

	integer, save :: nk_grd_alloc = 0
	integer, save :: ne_grd_alloc = 0
	integer, save :: nl_grd_alloc = 0
	integer, save :: nne_grd_alloc = 0
	integer, save :: nnl_grd_alloc = 0

	integer, save, allocatable :: ippnv(:),ippev(:),ipplv(:)
	integer, save, allocatable :: ianv(:),iaev(:),ialv(:)
        real, save, allocatable :: hhnv(:),hhev(:),hhlv(:)
        real, save, allocatable :: xv(:),yv(:)

        integer, save, allocatable :: ipntev(:),inodev(:)
        integer, save, allocatable :: ipntlv(:),inodlv(:)

	logical, save :: bcdepth = .true.	!needs complete set of depth

	logical, save :: bgrderror = .true.	!writes error if found
	logical, save :: bgrdwrite = .true.	!writes information messages

!==============================================================
	contains
!==============================================================

	subroutine grd_init(nkk,nee,nll,nnee,nnll)

	integer nkk,nee,nll,nnee,nnll

	integer nk,ne,nl,nne,nnl

	logical :: bdebug = .false.
	logical :: balloc

	if( nkk == 0 .and. nee == 0 .and. nll == 0 .and. &
     &				nnee == 0 .and. nnll == 0 ) then
	  balloc = .false.
	else
	  balloc = .true.
	end if

	nk = nkk
	ne = nee
	nl = nll
	nne = nnee
	nnl = nnll

	if( bdebug ) write(6,*) 'nk: ',nk,nk_grd,nk_grd_alloc

	if( nk .ne. nk_grd ) then
	  if( allocated(ippnv) ) then
	    deallocate(ippnv,ianv,hhnv,xv,yv)
	  end if
	  nk_grd_alloc = nk
	  if( balloc ) then
	    allocate(ippnv(nk),ianv(nk),hhnv(nk),xv(nk),yv(nk))
	  end if
	end if
	nk_grd = nk

	if( bdebug ) write(6,*) 'ne: ',ne,ne_grd,ne_grd_alloc

	if( ne .ne. ne_grd ) then
	  if( allocated(ippev) ) then
	    deallocate(ippev,iaev,hhev)
	  end if
	  ne_grd_alloc = ne
	  if( balloc ) then
	    allocate(ippev(ne),iaev(ne),hhev(ne))
	  end if
	end if
	ne_grd = ne
	if( allocated(ipntev) ) deallocate(ipntev)
	allocate(ipntev(0:ne))
	ipntev = 0

	if( bdebug ) write(6,*) 'nl: ',nl,nl_grd,nl_grd_alloc

	if( nl .ne. nl_grd_alloc ) then
	  if( allocated(ipplv) ) then
	    deallocate(ipplv,ialv,hhlv)
	  end if
	  nl_grd_alloc = nl
	  if( balloc ) then
	    allocate(ipplv(nl),ialv(nl),hhlv(nl))
	  end if
	end if
	nl_grd = nl
	if( allocated(ipntlv) ) deallocate(ipntlv)
	allocate(ipntlv(0:nl))
	ipntlv = 0

	if( bdebug ) write(6,*) 'nne: ',nne,nne_grd,nne_grd_alloc

	if( nne .ne. nne_grd_alloc ) then
	  if( allocated(inodev) ) then
	    deallocate(inodev)
	  end if
	  nne_grd_alloc = nne
	  if( balloc ) then
	    allocate(inodev(nne))
	  end if
	end if
	nne_grd = nne

	if( bdebug ) write(6,*) 'nnl: ',nnl,nnl_grd,nnl_grd_alloc

	if( nnl .ne. nnl_grd_alloc ) then
	  if( allocated(inodlv) ) then
	    deallocate(inodlv)
	  end if
	  nnl_grd_alloc = nnl
	  if( balloc ) then
	    allocate(inodlv(nnl))
	  end if
	end if
	nnl_grd = nnl

	end subroutine grd_init

!==============================================================
	end module grd
!==============================================================

!*****************************************************************

	subroutine grd_set_write(bset)

	use grd

	implicit none

	logical bset

	bgrdwrite = bset

	end subroutine grd_set_write

!*****************************************************************

	subroutine grd_set_error(berr)

	use grd

	implicit none

	logical berr

	bgrderror = berr

	end subroutine grd_set_error

!*****************************************************************

	function grd_write_error()

	use grd

	implicit none

	logical grd_write_error

	grd_write_error = bgrderror

	end function grd_write_error

!*****************************************************************

	subroutine grd_init_fake

	use grd

	implicit none

	call grd_init(1,1,1,1,1)

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine grd_get_params(nk,ne,nl,nne,nnl)

	use grd

	integer nk,ne,nl,nne,nnl

	nk = nk_grd
	ne = ne_grd
	nl = nl_grd
	nne = nne_grd
	nnl = nnl_grd

	end subroutine grd_get_params

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine grd_read(file)

	use grd

	implicit none

	character(*) file
	integer ner,nco
	logical berr,bwrite
	integer nk,ne,nl,nne,nnl
	integer nkndi,neldi,nlidi,nendi,nlndi

        ner = 6
        berr = .false.
	bwrite = .false.

        call grd_info(file,nk,ne,nl,nne,nnl,berr)
	if( berr ) goto 99

	if( bwrite ) write(6,*) 'grd_info: ',nk,ne,nl,nne,nnl

	nkndi = nk
	neldi = ne
	nlidi = nl
	nendi = nne
	nlndi = nnl

	call grd_init(nkndi,neldi,nlidi,nendi,nlndi)

	call rdgrd( &
     &			 file &
     &			,berr &
     &			,nco,nk,ne,nl,nne,nnl &
     &			,nkndi,neldi,nlidi,nendi,nlndi &
     &			,ippnv,ippev,ipplv &
     &			,ianv,iaev,ialv &
     &			,hhnv,hhev,hhlv &
     &			,xv,yv &
     &                  ,ipntev,inodev &
     &                  ,ipntlv,inodlv &
     &			)
	if( berr ) goto 99

	call ex2in(nk,nne,nnl,ippnv,inodev,inodlv,berr)
	if( berr ) goto 99

	if( bwrite ) then
	  write(6,*) 'total number of lines read: ',iline_grd
	end if

	return
   99	continue
	stop 'error stop grd_read: error reading grd file'
	end

!*****************************************************************

        subroutine grd_info(gfile,nk,ne,nl,nne,nnl,berr)

! reads grd file to obtain basic parameters

        implicit none

        character*(*) gfile
        integer nk              !total number of nodes
        integer ne              !total number of elems
        integer nl              !total number of lines
        integer nne             !total number of nodes in elems
        integer nnl             !total number of nodes in lines
        logical berr		!true if error

        integer ner,nco
        integer nkndi0,neldi0,nlidi0,nendi0,nlndi0

	integer ippnv(1),ippev(1),ipplv(1)
	integer ianv(1),iaev(1),ialv(1)
        real hhnv(1),hhev(1),hhlv(1)
        real xv(1),yv(1)

        integer ipntev0(0:0),inodev(1)
        integer ipntlv0(0:0),inodlv(1)

!-----------------------------------------------------------------
! initialize parameters
!-----------------------------------------------------------------

        ner = 6
        berr = .false.

        nkndi0 = 0
        neldi0 = 0
        nlidi0 = 0
        nendi0 = 0
        nlndi0 = 0

	call grd_init_fake

!-----------------------------------------------------------------
! read grd file
!-----------------------------------------------------------------

        call rdgrd( &
     &                   gfile &
     &                  ,berr &
     &                  ,nco,nk,ne,nl,nne,nnl &
     &                  ,nkndi0,neldi0,nlidi0,nendi0,nlndi0 &
     &                  ,ippnv,ippev,ipplv &
     &                  ,ianv,iaev,ialv &
     &                  ,hhnv,hhev,hhlv &
     &                  ,xv,yv &
     &                  ,ipntev0,inodev &
     &                  ,ipntlv0,inodlv &
     &                  )

        end

!**********************************************************

	function is_grd_file(gfile)

	implicit none

	logical is_grd_file
        character*(*) gfile

        logical berr
        integer nk,ne,nl,nne,nnl
        integer ner,nco
        integer nkndi0,neldi0,nlidi0,nendi0,nlndi0

	integer ippnv(1),ippev(1),ipplv(1)
	integer ianv(1),iaev(1),ialv(1)
        real hhnv(1),hhev(1),hhlv(1)
        real xv(1),yv(1)

        integer ipntev0(0:0),inodev(1)
        integer ipntlv0(0:0),inodlv(1)

!-----------------------------------------------------------------
! initialize parameters
!-----------------------------------------------------------------

        ner = 6
        berr = .false.

        nkndi0 = 0
        neldi0 = 0
        nlidi0 = 0
        nendi0 = 0
        nlndi0 = 0

	call grd_set_error(.false.)	!do not write errors to terminal

	call grd_init_fake

!-----------------------------------------------------------------
! read grd file
!-----------------------------------------------------------------

        call rdgrd( &
     &                   gfile &
     &                  ,berr &
     &                  ,nco,nk,ne,nl,nne,nnl &
     &                  ,nkndi0,neldi0,nlidi0,nendi0,nlndi0 &
     &                  ,ippnv,ippev,ipplv &
     &                  ,ianv,iaev,ialv &
     &                  ,hhnv,hhev,hhlv &
     &                  ,xv,yv &
     &                  ,ipntev0,inodev &
     &                  ,ipntlv0,inodlv &
     &                  )

	call grd_set_error(.true.)	!original behavior

	is_grd_file = .not. berr

        end

!**********************************************************

	subroutine grd_close

	use grd

	implicit none

	call grd_init(0,0,0,0,0)

	end

!**********************************************************

	subroutine grd_write(file)

	use grd

	implicit none

	character*(*) file

	integer nk,ne,nl,nne,nnl
	integer nco
	integer i,iu
	logical bdebug

	bdebug = .false.
	nco = 0

	call grd_get_params(nk,ne,nl,nne,nnl)

	if( bdebug ) then
	iu = 91
	write(iu,*) '========================================'
	write(iu,*) 'grd_write: ',nk,ne,nl,nne,nnl
	write(iu,*) trim(file)
	write(iu,*)
	write(iu,'(5i10)') (ippnv(i),i=1,nk)
	write(iu,*)
	write(iu,'(5i10)') (ippev(i),i=1,ne)
	write(iu,*)
	write(iu,'(5i10)') (ipntev(i),i=0,ne)
	write(iu,*)
	write(iu,'(5i10)') (inodev(i),i=1,nne)
	write(iu,*)
	write(iu,*) 'grd_write: ',nk,ne,nl,nne,nnl
	write(iu,*) '========================================'
	end if

	call grd_write_grid( &
     &			 file &
     &			,nco,nk,ne,nl,nne,nnl &
     &			,ippnv,ippev,ipplv &
     &			,ianv,iaev,ialv &
     &			,hhnv,hhev,hhlv &
     &			,xv,yv &
     &                  ,ipntev,inodev &
     &                  ,ipntlv,inodlv &
     &			)

	end

!**********************************************************
!**********************************************************
!**********************************************************

	subroutine rdgrd( &
     &			 file &
     &			,berr &
     &			,nco,nkn,nel,nli,nne,nnl &
     &			,nknddi,nelddi,nliddi,nenddi,nlnddi &
     &			,ipnv,ipev,iplv &
     &			,ianv,iaev,ialv &
     &			,hnv,hev,hlv &
     &			,xgv,ygv &
     &                  ,ipntev,inodev &
     &                  ,ipntlv,inodlv &
     &			)

! reads grd file
!
! after read nen3v and inodlv contain still external node numbers
! use ex2in to convert external to internal node numbers
!
! works only with triangles as elements

	!use grd

	implicit none

	character*(*) file	!file name

	logical berr		!true on return if error

	integer nco		!total number of comments read
	integer nkn		!total number of nodes read
	integer nel		!total number of elements read
	integer nli		!total number of lines read
	integer nne		!total number of nodes in elems
	integer nnl		!total number of nodes in lines

	integer nknddi		!dimension for number of nodes
	integer nelddi		!dimension for number of elements
	integer nliddi		!dimension for number of lines
	integer nenddi		!dimension for node numbers of elems
	integer nlnddi		!dimension for node numbers of lines

	integer ipnv(nknddi)	!external node number
	integer ipev(nelddi)	!external element number
	integer iplv(nliddi)	!external line number

	integer ianv(nknddi) 	!node type
	integer iaev(nelddi)	!element type
	integer ialv(nliddi)	!line type

	real hnv(nknddi)	!depth of node
	real hev(nelddi)	!depth of element
	real hlv(nliddi)	!depth of line

	real xgv(nknddi)	!x coordinate of node
	real ygv(nknddi)	!y coordinate of node

	integer ipntev(0:nelddi)!pointer into inodev
	integer inodev(nenddi)	!node numbers of elems
	integer ipntlv(0:nliddi)!pointer into inodlv
	integer inodlv(nlnddi)	!node numbers of lines

	logical berrwrite
	integer iwhat,ner
	real value

	logical grd_next_line

!--------------------------------------------------------------------
! initialize variables
!--------------------------------------------------------------------

	
	berrwrite = berr	!write read errors to terminal
	berr = .false.
	ner = 6

	nco=0
        nkn=0
        nel=0
        nli=0
	nne=0
	nnl=0

	ipntev(0) = 0
	ipntlv(0) = 0

!--------------------------------------------------------------------
! open file or STDIN
!--------------------------------------------------------------------

	call grd_internal_init(file)

!--------------------------------------------------------------------
! loop on lines and read
!--------------------------------------------------------------------

        do while( grd_next_line() )

	  call grd_ival(1,value)
	  iwhat = nint(value)

	  berr = berrwrite
          if( iwhat.eq.0 ) then  !comment or error
		call rdcom(nco,berr)
          else if(iwhat.eq.1) then              !node
        	call rdnode(nkn,nknddi,berr &
     &				,ipnv,ianv,xgv,ygv,hnv)
          else if(iwhat.eq.2) then              !element
        	call rdelem(nel,nne,nelddi,nenddi,berr &
     &				,ipev,iaev,ipntev,inodev,hev)
          else if(iwhat.eq.3) then              !line
        	call rdline(nli,nnl,nliddi,nlnddi,berr &
     &				,iplv,ialv,ipntlv,inodlv,hlv)
          else
	  	call rdunknown(iwhat,berr)
          end if

	  if( berr ) exit
        end do

	call grd_internal_close

	if( nkn == 0 ) berr = .true.	!we must really read something

!--------------------------------------------------------------------
! end of routine
!--------------------------------------------------------------------

	end

!**********************************************************

	subroutine rdcom(nco,berr)

! reads comment

	implicit none

	integer nco
	logical berr

	integer ner,iline
	integer i,n
	character*132 line

	integer ifstch
	logical grd_write_error

	call grd_line_info(iline,line)

	i=ifstch(line)
	n=len(line)

	ner = 6

	if(i.le.0) then
	  !blank line
	else if(line(i:i).eq.'0') then
	  nco=nco+1
	  if(i.lt.n) then
	    call fempar(line(i+1:))
	  end if
	else
	  if( grd_write_error() ) then
	    write(ner,*) 'Read error on line ',iline
	    write(ner,*) line
	  end if
          berr=.true.
	end if

	end

!**********************************************************

        subroutine rdnode(nkn,nknddi,berr &
     &				,ipnv,ianv,xgv,ygv,hnv)

! reads nodes from .grd file

	use grd, only : xscale_grd,yscale_grd,zscale_grd

        implicit none

	integer nkn,nknddi
	logical berr
        integer ipnv(nknddi)
        integer ianv(nknddi)         !node type
	real xgv(nknddi),ygv(nknddi)
	real hnv(nknddi)

	logical bread
	integer ner
	integer ianz
        real f(6)
	real depth

	ner = 6

	bread = nknddi .gt. 0		!read nodes?

	call grd_nvals(ianz)
        if(ianz.lt.5) goto 88
        if(ianz.gt.6) goto 87
	call grd_vals(ianz,f)

        nkn=nkn+1
	if( .not. bread ) return
	if(bread .and. nkn.gt.nknddi) then
	  bread = .false.
	  berr = .true.
	  if( nkn .eq. nknddi+1 ) then		!just one time
	    write(ner,*) 'dimension of nknddi too low: ',nknddi
	  end if
	end if
	if( .not. bread ) return

        ipnv(nkn)=nint(f(2))
        ianv(nkn)=nint(f(3))
        xgv(nkn)=f(4)*xscale_grd
        ygv(nkn)=f(5)*yscale_grd

       	depth = -999.
	if( ianz .ge. 6 ) depth = f(6)*zscale_grd
	hnv(nkn) = depth

        return
   87	continue
	write(ner,*) 'Too much data on line'
	call grd_write_line_info
	berr=.true.
	return
   88	continue
	write(ner,*) 'Not enough data on line'
	call grd_write_line_info
	berr=.true.
	return
   99	continue
	write(ner,*) 'dimension of nknddi too low: ',nknddi
	stop 'error stop rdnode: nknddi'
        end

!**********************************************************

        subroutine rdelem(nel,nne,nelddi,nenddi,berr &
     &				,ipev,iaev,ipntev,inodev,hev)

! reads elements from .grd file

        implicit none

	integer nel,nne,nelddi,nenddi
	logical berr
        integer ipev(nelddi),iaev(nelddi)
	integer ipntev(0:nelddi)
	integer inodev(max(1,nenddi))
	real hev(nelddi)

	logical bread
	integer ner,ii
        real f(4)
	integer ilist(10)
	integer inum,itype,ianz,ipnt
	integer ivert,nvert,istart
	real depth

	ner = 6

	bread = nelddi .gt. 0		!read elements?

	call grd_nvals(ianz)
        if(ianz.lt.4) goto 88
	call grd_vals(4,f)

        nel=nel+1
	if(bread .and. nel.gt.nelddi) then
	  bread = .false.
	  berr = .true.
	  if( nel .eq. nelddi+1 ) then		!just one time
	    write(ner,*) 'dimension of nelddi too low: ',nelddi
	  end if
	end if

        inum=nint(f(2))
        itype=nint(f(3))
        nvert=nint(f(4))

        if(nvert.gt.3) goto 87

        istart=4
        ivert=nvert
	if( .not. bread ) ivert = -ivert

        istart=4
        ivert=nvert
	if( bread ) then
	  ipnt = ipntev(nel-1)
	  if( ipnt + nvert .gt. nenddi ) goto 98
	else
	  ipnt = 0
	  ivert = -ivert
	end if

	call read_node_list(ivert,istart,inodev(ipnt+1),depth)

	if(ivert.lt.nvert) goto 86

	if( bread ) then
          ipev(nel) = inum
          iaev(nel) = itype
	  hev(nel)  = depth
	  ipntev(nel) = ipntev(nel-1) + nvert
	end if
	nne = nne + nvert

	return
   86	continue
	write(ner,*) 'Could not read all vertices for element'
	write(ner,*) '   internal = ',nel,'    external = ', inum
	call grd_write_line_info
	berr=.true.
	return
   87	continue
	write(ner,*) 'Not a triangle on line'
	call grd_write_line_info
	berr=.true.
	return
   88	continue
	write(ner,*) 'Not enough data on line'
	call grd_write_line_info
	berr=.true.
	return
   98   continue
        write(6,*) 'dimension of nenddi too low: ',nenddi
        stop 'error stop rdelem: nenddi'
   99	continue
	write(ner,*) 'dimension of nelddi too low: ',nelddi
	stop 'error stop rdelem: nelddi'
	end

!**********************************************************

        subroutine rdline(nli,nnl,nliddi,nlnddi,berr &
     &				,iplv,ialv,ipntlv,inodlv,hlv)

! reads lines from .grd file

        implicit none

	integer nli,nnl,nliddi,nlnddi
	logical berr
        integer iplv(nliddi),ialv(nliddi)
	integer ipntlv(0:nliddi)		!pointer into inodlv
	integer inodlv(max(1,nlnddi))
	real hlv(nliddi)

	logical bread
	integer ner
	integer inum,itype,ianz,ipnt
	integer ivert,nvert,istart
        real f(4)
	real depth

	ner = 6

	bread = nliddi .gt. 0		!read lines?

	call grd_nvals(ianz)
        if(ianz.lt.4) goto 88
	call grd_vals(4,f)

        nli=nli+1
	if(bread .and. nli.gt.nliddi) goto 99
        inum=nint(f(2))
        itype=nint(f(3))
        nvert=nint(f(4))

        istart=4
        ivert=nvert
	if( bread ) then
	  ipnt = ipntlv(nli-1)
	  if( ipnt + nvert .gt. nlnddi ) goto 98
	else
	  ipnt = 0
	  ivert = -ivert
	end if

	call read_node_list(ivert,istart,inodlv(ipnt+1),depth)

	if(ivert.lt.nvert) goto 86

	if( bread ) then
          iplv(nli) = inum
          ialv(nli) = itype
	  hlv(nli)  = depth
	  ipntlv(nli) = ipntlv(nli-1) + nvert
	end if
	nnl = nnl + nvert

	return
   86	continue
	write(ner,*) 'Could not read all nodes for line ',nli
	call grd_write_line_info
	berr=.true.
	return
   88	continue
	write(ner,*) 'Not enough data on line'
	call grd_write_line_info
	berr=.true.
	return
   98	continue
	write(6,*) 'dimension of nlnddi too low: ',nlnddi
	stop 'error stop rdline: nlnddi'
   99	continue
	write(6,*) 'dimension of nliddi too low: ',nliddi
	stop 'error stop rdline: nliddi'
	end

!**********************************************************

	subroutine read_node_list(nvert,istart,nodes,depth)

! reads node list

	use grd, only : zscale_grd

	implicit none

	integer nvert		!number of vertices to read
	integer istart
	integer nodes(abs(nvert))
	real depth

	logical bread,bline
	integer i,ivert,ianz
	real value

	logical grd_next_line

	bread = nvert .gt. 0
	nvert = abs(nvert)

        i=istart
        ivert=0
	bline = .true.

	call grd_nvals(ianz)

        do while( bline .and. ivert .lt. nvert )    !read vertices

          if( i .ge. ianz ) then			!read new line
	    i = 0
	    bline = grd_next_line()
	    call grd_nvals(ianz)
	  else						!get new node number
            i=i+1
            ivert=ivert+1
	    call grd_ival(i,value)
            if( bread ) nodes(ivert) = nint(value)
	  end if

	end do

       	depth = -999.
	if( i .lt. ianz ) then
	  call grd_ival(i+1,value)
	  depth = value*zscale_grd
	end if

	nvert = ivert		!pass back number of vertices read

	end

!******************************************************************************

        subroutine rdunknown(iwhat,berr)

! handles unknown type

	logical berr

	integer iwhat
	integer iline
	character*132 line

	integer ner

	ner = 6

	if( berr ) then
	  call grd_line_info(iline,line)
          write(ner,*) 'Type ',iwhat,' at line ',iline,' not recognized'
          write(ner,*) line
	end if

        berr=.true.

	end

!******************************************************************************

	function ifstch(line)

! finds first char of line that is not blank or tab

        implicit none

	integer ifstch
	character*(*) line

	integer length,i

	length=len(line)

	do i=1,length
	  if(line(i:i).ne.' ' .and. line(i:i).ne.'	' ) then
		ifstch=i
		return
	  end if
	end do

	ifstch=-1

	return
	end

!******************************************************************************

	subroutine fempar(gline)

! read parameters for fem model 
!
! the following special symbols are recognized on a comment line :
!
! (FEM-TITLE)	title for basin			default: first comment line
! (FEM-SCALE)	scale in x/y/z for basin	default: 1.,1.,1.
! (FEM-LATID)	latitude of basin (degrees)	default: 0.
! (FEM-NORTH)	true north of basin (degrees) 	default: 0.
! (FEM-SPHER)	coordinates are lat/lon 	default: -1 (must be determined)
!
! all entries are optional
!
! example (! of fortran comment must be deleted, line starts with 0) :
!
! 0 (FEM-TITLE) test basin
! 0 (FEM-SCALE) 0.5 0.5 2.
! 0 (FEM-LATID) 45.0
! 0 (FEM-NORTH) 90.0
! 0 (FEM-SPHER) 0
!
	use grd
	!use basin

        implicit none

	character*(*) gline

	integer i,j,n
	integer ifstch,iscanf
	logical, save :: btitle = .false.
	logical berror

!----------------------------------------------------
! initialize parameters
!----------------------------------------------------

	berror = .false.

	i=ifstch(gline)
	n=len(gline)

!----------------------------------------------------
! check comments for FEM information
!----------------------------------------------------

	if( i.gt.0 .and. i+10.lt.n ) then
	  if( gline(i:i+10) .eq. '(FEM-TITLE)' ) then
		title_grd=gline(i+11:)
		btitle=.true.
	  else if( gline(i:i+10) .eq. '(FEM-SCALE)' ) then
		j=iscanf(gline(i+11:),f_grd,4)
		if(j.eq.3) then
		  xscale_grd=f_grd(1)
		  yscale_grd=f_grd(2)
		  zscale_grd=f_grd(3)
		else
		  write(6,*) 'error reading (FEM-SCALE) :',j
		  write(6,*) gline
		  write(6,*) gline(i+11:)
		  berror = .true.
		end if
	  else if( gline(i:i+10) .eq. '(FEM-LATID)' ) then
		j=iscanf(gline(i+11:),f_grd,2)
		if(j.eq.1) then
		  dcor_grd=f_grd(1)
		else
		  write(6,*) 'error reading (FEM-LATID) :'
		  write(6,*) gline
		  berror = .true.
		end if
	  else if( gline(i:i+10) .eq. '(FEM-NORTH)' ) then
		j=iscanf(gline(i+11:),f_grd,2)
		if(j.eq.1) then
		  dirn_grd=f_grd(1)
		else
		  write(6,*) 'error reading (FEM-NORTH) :'
		  write(6,*) gline
		  berror = .true.
		end if
	  else if( gline(i:i+10) .eq. '(FEM-SPHER)' ) then
		j=iscanf(gline(i+11:),f_grd,2)
		if(j.eq.1) then
		  sphe_grd=f_grd(1)
		else
		  write(6,*) 'error reading (FEM-SPHER) :'
		  write(6,*) gline
		  berror = .true.
		end if
	  end if
	end if

!----------------------------------------------------
! check parameters read
!----------------------------------------------------

	if( dcor_grd < -90. .or. dcor_grd > 90. ) then
	  write(6,*) 'value for FEM-LATID out of range: ',dcor_grd
	  berror = .true.
	else if( dirn_grd < -360. .or. dirn_grd > 360. ) then
	  write(6,*) 'value for FEM-NORTH out of range: ',dirn_grd
	  berror = .true.
	else if( sphe_grd < -1 .or. sphe_grd > 1. ) then
	  write(6,*) 'value for FEM-SPHER out of range: ',sphe_grd
	  berror = .true.
	end if
	
!----------------------------------------------------
! error handling
!----------------------------------------------------

	if( berror ) then
	  stop 'error stop fempar: error reading FEM comments'
	end if

!----------------------------------------------------
! use first comment as title
!----------------------------------------------------

	if( i.gt.0 .and. .not.btitle ) then
		title_grd=gline(i:)
		btitle=.true.
	end if

!----------------------------------------------------
! end of routine
!----------------------------------------------------

	end

!******************************************************************************
!******************************************************************************
!******************************************************************************

	subroutine ex2in(nkn,ne,nl,ippnv,inodev,inodlv,berr)

! changing extern with intern node numbers in elements and lines
!
! if no elements or lines are given, set ne or nl to 0

	implicit none

        integer nkn,ne,nl
        logical berr
        integer ippnv(nkn)
        integer inodev(ne)
        integer inodlv(nl)

        integer ipaux(nkn)	!local

	call isort(nkn,ippnv,ipaux)
	call chex2in(nkn,ne,inodev,ippnv,ipaux,berr)
	call chex2in(nkn,nl,inodlv,ippnv,ipaux,berr)

	end

!*****************************************************************
 
        subroutine chex2in(nkn,n,nodes,ipnv,index,berr)
 
! changing extern with intern node numbers node list
 
        implicit none
 
        integer nkn,n
        logical berr
        integer nodes(n)
        integer ipnv(nkn)
        integer index(nkn)
 
        integer k,i,kn
	integer ner
        integer locate
 
	ner = 6

        do i=1,n
            kn=nodes(i)
            k=locate(nkn,ipnv,index,kn)
            if(k.le.0) then
              write(ner,*)' warning : node',kn,' not found'
              berr=.true.
            else
              nodes(i)=k
            end if
        end do
 
        end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine grd_internal_init(file)

! initializes reading from grid file

	use grd

	implicit none

	character*(*) file

	integer ner
	integer ifileo

	nin_grd = 0
	iline_grd = 0
	ianz_grd = 0
	ner = 6

	xscale_grd=1.
	yscale_grd=1.
	zscale_grd=1.

        if( file .eq. ' ' ) then
          write(ner,*) 'no file given...'
          write(ner,*) 'cannot read GRD file from stdin...'
	  stop 'error stop grd_internal_init: cannot read file'
        else
	  nin_grd = ifileo(nin_grd,file,'formatted','old')
        end if

	if(nin_grd.le.0) then
	  write(ner,*) 'error opening file'
	  write(ner,*) file
	  stop 'error stop grd_internal_init: cannot open file'
	end if

	end

!*****************************************************************

	subroutine grd_internal_close

	use grd

	implicit none

	if( nin_grd .ne. 5 ) close(nin_grd)
	nin_grd = 0

	end

!*****************************************************************

	function grd_next_line()

! reads next line from file

	use grd

	implicit none

	logical grd_next_line

	integer ner,ios
	integer iscanf

	ner = 6
	grd_next_line = .false.

        read(nin_grd,'(a)',iostat=ios) line_grd

        if( ios .lt. 0 ) then
          !write(ner,*) iline,' lines read'
	  close(nin_grd)
	  return
        else if( ios .gt. 0 ) then
          write(ner,*) 'read error close to line ',iline_grd
          write(ner,*) line_grd
          stop 'error stop grd_next_line: reading line'
        end if

	iline_grd = iline_grd + 1
        ianz_grd=iscanf(line_grd,f_grd,81)
	if( ianz_grd .gt. 80 ) goto 99

	grd_next_line = .true.

	return
   99	continue
	write(ner,*) 'dimension of array f too small ',80,ianz_grd
	stop 'error stop grd_next_line: ianz_grd'
	end

!*****************************************************************

	subroutine grd_nvals(nvals)

! returns number of values on line

	use grd

	implicit none

	integer nvals

	nvals = ianz_grd

	end

!*****************************************************************

	subroutine grd_vals(nvals,vals)

! returns nvals in vals

	use grd

	implicit none

	integer nvals
	real vals(nvals)

	integer i,n,nmin,nmax

	n = min(ianz_grd,nvals)

	do i=1,n
	  vals(i) = f_grd(i)
	end do

!	if nvals is greater than ianz_grd return zero

	nmin = max(1,n+1)
	nmax = min(80,nvals)

	do i=nmin,nmax
	  vals(i) = 0.
	end do

	end

!*****************************************************************

	subroutine grd_ival(ival,val)

! returns value at position ival

	use grd

	implicit none

	integer ival
	real val

	val = 0.
	if( ival .ge. 1 .and. ival .le. ianz_grd ) val = f_grd(ival)

	end

!*****************************************************************

	subroutine grd_line_info(iline_gr,line_gr)

! returns info on line

	use grd

	implicit none

	integer iline_gr
	character*(*) line_gr

	iline_gr = iline_grd
	line_gr = line_grd

	end

!*****************************************************************

	subroutine grd_write_line_info

! write info on line

	use grd

	implicit none

	integer ner

	ner = 6

	write(ner,*) 'line number: ',iline_grd
	write(ner,*) line_grd

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine grd_write_grid( &
     &			 file &
     &			,nco,nk,ne,nl,nne,nnl &
     &			,ippnv,ippev,ipplv &
     &			,ianv,iaev,ialv &
     &			,hhnv,hhev,hhlv &
     &			,xv,yv &
     &                  ,ipntev,inodev &
     &                  ,ipntlv,inodlv &
     &			)

! writes grd file

	implicit none

	character*(*) file	!file name

	integer nco		!total number of comments read
	integer nk		!total number of nodes read
	integer ne		!total number of elements read
	integer nl		!total number of lines read
	integer nne		!total number of nodes in elems 
	integer nnl		!total number of nodes in lines

	integer ippnv(nk)	!external node number
	integer ippev(ne)	!external element number
	integer ipplv(nl)	!external line number

	integer ianv(nk) 	!node type
	integer iaev(ne)	!element type
	integer ialv(nl)	!line type

	real hhnv(nk)		!depth of node
	real hhev(ne)		!depth of element
	real hhlv(nl)		!depth of line

	real xv(nk)		!x coordinate of node
	real yv(nk)		!y coordinate of node

	integer ipntev(0:ne)	!pointer into inodev
	integer inodev(nne)	!node numbers of elems
	integer ipntlv(0:nl)	!pointer into inodlv
	integer inodlv(nnl)	!node numbers of lines

	integer nout
	integer i
	integer n,ib,nmax
	integer k,ie,il,kext,itype,in
	integer ia,ieext,ilext
	real x,y,depth
	logical bsort,bextern
	integer, allocatable :: ipdex(:)
	integer, allocatable :: nextern(:)
	integer, allocatable :: nodes(:)

	integer ifileo

	bsort = .true.		!sort nodes on external numbering
	bsort = .false.

	bextern = .false.
	bextern = .true.	!use external nodes

	if( .not. bextern ) bsort = .false.

	nout = ifileo(1,file,'formatted','unknown')
	if( nout .le. 0 ) goto 99

	nmax = max(nk,ne,nl)
	allocate(ipdex(nmax))
	do i=1,nmax
	  ipdex(i) = i
	end do

	allocate(nextern(nk),nodes(nk))
	do k=1,nk
	  if( bextern ) then
	    nextern(k) = ippnv(k)		!use external numbering
	  else
	    nextern(k) = k			!use internal numbering
	  end if
	end do

	write(nout,*)

	if( bsort ) call isort(nk,ippnv,ipdex)

	do i=1,nk
	  k = ipdex(i)
	  kext = nextern(k)
	  itype = ianv(k)
	  x = xv(k)
	  y = yv(k)
	  depth = hhnv(k)
	  call grd_write_node(nout,kext,itype,x,y,depth)
	end do

	write(nout,*)

	if( bsort ) call isort(ne,ippev,ipdex)

	do i=1,ne
	  ie = ipdex(i)
	  ieext = ie
	  if( bextern ) ieext = ippev(ie)
	  ia = iaev(ie)
	  depth = hhev(ie)
	  n = ipntev(ie) - ipntev(ie-1)
	  ib = ipntev(ie-1)
	  do in=1,n
	    nodes(in) = nextern(inodev(ib+in))
	  end do
	  call grd_write_item(nout,2,ieext,ia,n,nodes,depth)
	end do

	write(nout,*)

	if( bsort ) call isort(nl,ipplv,ipdex)

	do i=1,nl
	  il = ipdex(i)
	  ilext = il
	  if( bextern ) ilext = ipplv(il)
	  ia = ialv(il)
	  depth = hhlv(il)
	  n = ipntlv(il) - ipntlv(il-1)
	  ib = ipntlv(il-1)
	  do in=1,n
	    nodes(in) = nextern(inodlv(ib+in))
	  end do
	  call grd_write_item(nout,3,ilext,ia,n,nodes,depth)
	end do

	write(nout,*)

	close(nout)

	deallocate(ipdex)
	deallocate(nextern,nodes)

	return
   99	continue
	write(6,*) 'error opening output file'
	write(6,*) file
	stop 'error stop grd_write_grid: cannot open file'
	end

!*****************************************************************

	subroutine grd_write_node(nout,number,itype,x,y,depth)

	implicit none

	integer nout
	integer number
	integer itype
	real x,y
	real depth

	real, parameter :: flag = -999.

	if( depth .eq. flag ) then
	  write(nout,1000) 1,number,itype,x,y
	else
	  write(nout,1000) 1,number,itype,x,y,depth
	end if

	return
 1000	format(i1,2i10,3e16.8)
	end

!*****************************************************************

	subroutine grd_write_item(nout,iwhat,number,itype,n, &
     &				nodes,depth)

	implicit none

	integer nout
	integer iwhat
	integer number
	integer itype
	integer n
	integer nodes(n)
	real depth

	integer i
	real, parameter :: flag = -999.

	if( n > 3 ) then
	  write(nout,3000) iwhat,number,itype,n
	  call grd_write_node_list(nout,n,nodes,depth)
	else if( depth == flag ) then
	  write(nout,2000) iwhat,number,itype,n, &
     &			(nodes(i),i=1,n)
	else
	  if( n == 2 ) then
	    write(nout,2002) iwhat,number,itype,n, &
     &			(nodes(i),i=1,n),depth
	  else
	    write(nout,2000) iwhat,number,itype,n, &
     &			(nodes(i),i=1,n),depth
	  end if
	end if

	return
 2002	format(i1,5i10,e16.8)
 2000	format(i1,6i10,e16.8)
 3000	format(i1,3i10)
	end

!*****************************************************************

	subroutine grd_write_node_list(nout,n,nodes,depth)

! writes out node list

	implicit none

	integer nout
	integer n
	integer nodes(1)
	real depth

	integer i,iend,ii,nend
	integer lnodes(6)
	character*30 format

	do i=1,n,6
	  iend = min(n,i+5)
	  nend = iend - i + 1
	  !write(6,*) i,n,iend,nend
	  do ii=1,nend
	    lnodes(ii) = nodes(i+ii-1)
	  end do
	  if( iend .eq. n .and. depth .ne. -999. ) then
	      write(format,'(a4,i1,a10)') '(1x,',nend,'i10,e16.8)'
	      write(nout,format) (lnodes(ii),ii=1,nend),depth
	  else
	      write(nout,3001) (lnodes(ii),ii=1,nend)
	  end if
	end do

	return
 3001	format(1x,6i10)
	end

!*****************************************************************

	subroutine grd_get_depth(nk,ne,hkv,hev)

	use grd

	implicit none

	integer nk,ne
	real hkv(nk)
	real hev(ne)

	if( ne > 0 ) then
	  if( ne .ne. ne_grd ) then
	    write(6,*) 'ne,ne_grd: ',ne,ne_grd
	    stop 'error stop grd_get_depth: dimension mismatch'
	  end if
	  hev = hhev
	end if

	if( nk > 0 ) then
	  if( nk .ne. nk_grd ) then
	    write(6,*) 'nk,nk_grd: ',nk,nk_grd
	    stop 'error stop grd_get_depth: dimension mismatch'
	  end if
	  hkv = hhnv
	end if

	end

!*****************************************************************

	subroutine grd_get_nodes(np,xp,yp,hp)

	use grd

	implicit none

	integer np
	real xp(np)
	real yp(np)
	real hp(np)

	if( np < nk_grd ) then
	  write(6,*) 'np,nk_grd: ',np,nk_grd
	  stop 'error stop grd_get_nodes: np'
	end if

	np = nk_grd
	xp(1:nk_grd) = xv(1:nk_grd)
	yp(1:nk_grd) = yv(1:nk_grd)
	hp(1:nk_grd) = hhnv(1:nk_grd)

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine grd_get_total_lines(nl)

	use grd

	implicit none

	integer nl

	nl = nl_grd

	end

!*****************************************************************

	subroutine grd_get_line_params(il,inum,itype,nvert,depth)

	use grd

	implicit none

	integer il
	integer inum
	integer itype
	integer nvert
	real depth

        inum = ipplv(il)
        itype = ialv(il)
	depth = hhlv(il)
	nvert = ipntlv(il) - ipntlv(il-1)

	end

!*****************************************************************

	subroutine grd_get_line_array(il,nvert,nodes,x,y)

	use grd

	implicit none

	integer il
	integer nvert
	integer nodes(nvert)
	real x(nvert)
	real y(nvert)

	integer ibase,i,node

	ibase = ipntlv(il-1)
	nvert = ipntlv(il) - ipntlv(il-1)
	nodes(1:nvert) = inodlv(ibase+1:ibase+nvert)

	do i=1,nvert
	  node = nodes(i)
	  x(i) = xv(node)
	  y(i) = yv(node)
	end do

	end

!*****************************************************************

	subroutine grd_get_total_nodes(nk)

	use grd

	implicit none

	integer nk

	nk = nk_grd

	end

!*****************************************************************

	subroutine grd_get_node_params(ik,inum,itype,x,y,depth)

	use grd

	implicit none

	integer ik
	integer inum
	integer itype
	real x,y
	real depth

        inum = ippnv(ik)
        itype = ianv(ik)
	x = xv(ik)
	y = yv(ik)
	depth = hhnv(ik)

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine grd_get_node_arrays(n,ipv,iav,hv,x,y)

	use grd

	implicit none

	integer n
	integer ipv(n),iav(n)
	real hv(n),x(n),y(n)

	if( n /= nk_grd ) stop 'error stop grd_get_node_arrays: n'

	ipv = ippnv
	iav = ianv
	hv = hhnv
	x = xv
	y = yv

	end 

!*****************************************************************

	subroutine grd_get_elem_arrays(n,ipv,iav,hv)

	use grd

	implicit none

	integer n
	integer ipv(n),iav(n)
	real hv(n)

	if( n /= ne_grd ) stop 'error stop grd_get_node_arrays: n'

	ipv = ippev
	iav = iaev
	hv = hhev

	end 

!*****************************************************************
!*****************************************************************
!*****************************************************************
! line routines
!*****************************************************************
!*****************************************************************
!*****************************************************************
!
! all following routines read lines from GRD, BND or XY format
! they return n,x,y,ifl
! if n == 0 on entry only counting of x/y points
! if n == 0 on return an error has occured
! ifl == 1 for start of line and 0 for following points
! XY format can read only one line
!
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine read_all_lines(file,n,x,y,ifl)

! reads file with lines and general format (grd,bnd,xy)

	implicit none

	character*(*) file
	integer n		!0 for counting points, >0 for reading them
	real x(n),y(n)
	integer ifl(n)

	logical is_grd_file,is_bnd_file,is_xy_file

	if( is_grd_file(file) ) then
	  call read_grd_lines(file,n,x,y,ifl)
	else if( is_bnd_file(file) ) then
	  call read_bnd_lines(file,n,x,y,ifl)
	else if( is_xy_file(file) ) then
	  call read_xy_lines(file,n,x,y,ifl)
	else
	  write(6,*) 'cannot determine file format: ',trim(file)
	  stop 'error stop read_all_lines: file format'
	end if

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine read_grd_lines(grdfile,n,x,y,ifl)

! reads grd file for lines to be treated (plot or particle release)

	implicit none

	character*(*) grdfile
	integer n		!0 for counting points, >0 for reading them
	real x(n),y(n)
	integer ifl(n)

	integer nl,il,inum,itype,nvert,ndim,ip
	real depth

	logical is_grd_file

	ndim = n
	n = 0

	if( .not. is_grd_file(grdfile) ) then
	  write(6,*) 'file is not a grd file: ',trim(grdfile)
	  write(6,*) 'cannot read old format of bnd files'
	  write(6,*) 'use bnd2grd to convert from bnd to grd file'
	  !stop 'error stop read_bnd_lines: no bnd file read'
	  return
	end if

	call grd_read(grdfile)

	call grd_get_total_lines(nl)

	if( ndim == 0 ) then
	  do il=1,nl
	    call grd_get_line_params(il,inum,itype,nvert,depth)
	    n = n + nvert
	    !write(6,*) 'counting line: ',nl,il,nvert
	  end do
	  !write(6,*) 'total number of points in line: ',n
	else
	  ip = 1
	  do il=1,nl
	    call grd_get_line_array(il,nvert,ifl(ip),x(ip),y(ip))
	    !write(6,*) 'reading line: ',nl,il,nvert,ip
	    ifl(ip:ip+nvert-1) = 0
	    ifl(ip) = 1
	    ip = ip + nvert
	    if( ip-1 > ndim ) goto 99
	  end do
	  n = ip - 1
	  !write(6,*) 'number of points in line read: ',n
	end if

	return
   99	continue
	write(6,*) 'ndim = ',ndim
	stop 'error stop read_grd_lines: ndim'
	end

!*****************************************************************

	subroutine read_bnd_lines(bndfile,n,x,y,ifl)

! reads old bnd file format

	implicit none

	character*(*) bndfile
	integer n		!0 for counting points, >0 for reading them
	real x(n),y(n)
	integer ifl(n)

	integer ndim,ios,iflag
	real xx,yy

	ndim = n
	n = 0

	open(1,iostat=ios,file=bndfile,status='old',form='formatted')
	if( ios /= 0 ) return

	do
	  read(1,*,iostat=ios) xx,yy,iflag
	  if( ios /= 0 ) exit
	  n = n + 1
	  if( ndim > 0 ) then
	    if( n > ndim ) goto 99
	    x(n) = xx
	    y(n) = yy
	    ifl(n) = iflag
	  end if
	end do

	close(1)

	if( ios > 0 ) n = 0

	return
   99	continue
	write(6,*) 'ndim = ',ndim
	stop 'error stop read_bnd_lines: ndim'
	end

!*****************************************************************

	subroutine read_xy_lines(xyfile,n,x,y,ifl)

! reads xy file format - only one line allowed

	implicit none

	character*(*) xyfile
	integer n		!0 for counting points, >0 for reading them
	real x(n),y(n)
	integer ifl(n)

	integer ndim,ios,iflag
	real xx,yy

	ndim = n
	n = 0
	iflag = 1	!only first point flagged with 1, rest 0

	open(1,iostat=ios,file=xyfile,status='old',form='formatted')
	if( ios /= 0 ) return

	do
	  read(1,*,iostat=ios) xx,yy
	  if( ios /= 0 ) exit
	  n = n + 1
	  if( ndim > 0 ) then
	    if( n > ndim ) goto 99
	    x(n) = xx
	    y(n) = yy
	    ifl(n) = iflag
	  end if
	  iflag = 0
	end do

	close(1)

	if( ios > 0 ) n = 0

	return
   99	continue
	write(6,*) 'ndim = ',ndim
	stop 'error stop read_xy_lines: ndim'
	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	function is_bnd_file(file)

	implicit none

	logical is_bnd_file
	character*(*) file

	character*80 line
	integer ios,ianz
	real f(4)
	integer iscanf

	is_bnd_file = .false.

	open(1,iostat=ios,file=file,status='old',form='formatted')
	if( ios /= 0 ) return
	read(1,*,iostat=ios) line
	close(1)

	if( ios /= 0 ) return
	ianz = iscanf(line,f,4)
	if( ianz /= 3 ) return
	if( nint(f(3)) /= 0 .and. nint(f(3)) /= 1 ) return

	is_bnd_file = .true.

	end

!*****************************************************************

	function is_xy_file(file)

	implicit none

	logical is_xy_file
	character*(*) file

	character*80 line
	integer ios,ianz
	real f(4)
	integer iscanf

	is_xy_file = .false.

	open(1,iostat=ios,file=file,status='old',form='formatted')
	if( ios /= 0 ) return
	read(1,*,iostat=ios) line
	close(1)

	if( ios /= 0 ) return
	ianz = iscanf(line,f,4)
	if( ianz /= 2 ) return

	is_xy_file = .true.

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************
! simplified routines for grd writing
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine write_grd_file(file,text,nk,ne,xg,yg,index &
     &			,inext,ieext,intype,ietype)

! writes grd file with no depth information

	implicit none

	character*(*) file		!file name
	character*(*) text		!text string
	integer nk,ne			!total nodes and elements
	real xg(nk),yg(nk)		!coordinates
	integer index(3,ne)		!element index (internal numbers)
	integer inext(nk),ieext(ne)	!external nodes and element numbers
	integer intype(nk),ietype(ne)	!node and element type

	integer nout
	integer k,ie,ii,itype,kext,eext,n
	integer nen3v(3)
	real depth,x,y
	real, parameter :: flag = -999.

	nout = 1
	open(nout,file=file,status='unknown',form='formatted')

	depth = flag

	if( text /= ' ' ) then
	  write(nout,*)
	  write(nout,'(a,a)') '0 ',trim(text)
	  write(nout,*)
	end if

	do k=1,nk
	  kext = inext(k)
	  itype = intype(k)
	  x = xg(k)
	  y = yg(k)
	  call grd_write_node(nout,kext,itype,x,y,depth)
	end do

	n = 3
	do ie=1,ne
	  eext = ieext(ie)
	  itype = ietype(ie)
	  do ii=1,n
	    k = index(ii,ie)
	    nen3v(ii) = inext(k)
	  end do
          call grd_write_item(nout,2,eext,itype,n,nen3v,depth)
	end do

	close(nout)

	end

!*****************************************************************

	subroutine write_grd_file_with_depth(file &
     &			,text,nk,ne,xg,yg,index &
     &			,inext,ieext,intype,ietype &
     &			,rndepth,redepth)

! writes grd file with depth information

	implicit none

	character*(*) file		!file name
	character*(*) text		!text string
	integer nk,ne			!total nodes and elements
	real xg(nk),yg(nk)		!coordinates
	integer index(3,ne)		!element index (internal numbers)
	integer inext(nk),ieext(ne)	!external nodes and element numbers
	integer intype(nk),ietype(ne)	!node and element type
	real rndepth(nk),redepth(ne)	!node and element depths

	integer nout
	integer k,ie,ii,itype,kext,eext,n
	integer nen3v(3)
	real depth,x,y
	real, parameter :: flag = -999.

	nout = 1
	open(nout,file=file,status='unknown',form='formatted')

	depth = flag

	if( text /= ' ' ) then
	  write(nout,*)
	  write(nout,'(a,a)') '0 ',trim(text)
	  write(nout,*)
	end if

	do k=1,nk
	  kext = inext(k)
	  itype = intype(k)
	  depth = rndepth(k)
	  x = xg(k)
	  y = yg(k)
	  call grd_write_node(nout,kext,itype,x,y,depth)
	end do

	n = 3
	do ie=1,ne
	  eext = ieext(ie)
	  itype = ietype(ie)
	  depth = redepth(ie)
	  do ii=1,n
	    k = index(ii,ie)
	    nen3v(ii) = inext(k)
	  end do
          call grd_write_item(nout,2,eext,itype,n,nen3v,depth)
	end do

	close(nout)

	end

!*****************************************************************

	subroutine write_grd_file_with_detail(file &
     &			,text,nk,ne,xg,yg,index &
     &			,inext,ieext,intype,ietype &
     &			,kcext,window)

! writes grd file with depth information

	implicit none

	character*(*) file		!file name
	character*(*) text		!text string
	integer nk,ne			!total nodes and elements
	real xg(nk),yg(nk)		!coordinates
	integer index(3,ne)		!element index (internal numbers)
	integer inext(nk),ieext(ne)	!external nodes and element numbers
	integer intype(nk),ietype(ne)	!node and element type
	integer kcext			!central node (external number)
	real window			!window around node (fact)

	integer nout
	integer k,ie,ii,itype,kext,eext,n,i
	integer nen3v(3)
	real depth,x,y
	real, parameter :: flag = -999.
	character*80 zoom

	integer, parameter :: ndim = 100
	integer kn(3*ndim)
	real xn(3*ndim),yn(3*ndim)	!be sure the array is big enough
	integer kc,ksum
	real xmin,ymin,xmax,ymax,dx,dy
	integer, allocatable :: kin(:),ein(:)

!---------------------------------------------------------
! determine min/max of detail
!---------------------------------------------------------

	do k=1,nk
	  if( inext(k) == kcext ) exit
	end do
	kc = k
	if( kc > nk ) then
	  write(6,*) 'kcext = ',kcext
	  stop 'error stopwrite_grd_file_with_detail: no such node'
	end if

	!write(6,*) 'node (ext/int) : ',kcext,kc

	n = 0
	do ie=1,ne
	  if( all( index(:,ie) /= kc ) ) cycle
	  if( n+3 > ndim ) stop 'error stop write_grd_file_with_detail: ndim'
	  do ii=1,3
	    n = n + 1
	    k = index(ii,ie)
	    kn(n) = k
	    xn(n) = xg(k)
	    yn(n) = yg(k)
	  end do
	end do

	xmin = minval(xn(1:n))
	ymin = minval(yn(1:n))
	xmax = maxval(xn(1:n))
	ymax = maxval(yn(1:n))
	dx = xmax - xmin
	dy = ymax - ymin
	xmin = xmin - window*dx
	ymin = ymin - window*dy
	xmax = xmax + window*dx
	ymax = ymax + window*dy

	allocate(kin(nk))
	allocate(ein(ne))
	kin = 0
	ein = 0

	do k=1,nk
	  x = xg(k)
	  y = yg(k)
	  if( x < xmin .or. x > xmax ) cycle
	  if( y < ymin .or. y > ymax ) cycle
	  kin(k) = 1
	end do

	do ie=1,ne
	  ksum = 0
	  do ii=1,3
	    k = index(ii,ie)
	    ksum = ksum + kin(k)
	  end do
	  if( ksum == 3 ) ein(ie) = 1
	end do

!---------------------------------------------------------
! now write grd
!---------------------------------------------------------
	
	depth = flag

	nout = 1
	open(nout,file=file,status='unknown',form='formatted')

	depth = flag

	if( text /= ' ' ) then
	  write(nout,*)
	  write(nout,'(a,a)') '0 ',trim(text)
	  write(nout,*)
	end if

	do k=1,nk
	  if( kin(k) == 0 ) cycle
	  kext = inext(k)
	  itype = intype(k)
	  x = xg(k)
	  y = yg(k)
	  call grd_write_node(nout,kext,itype,x,y,depth)
	end do

	n = 3
	do ie=1,ne
	  if( ein(ie) == 0 ) cycle
	  eext = ieext(ie)
	  itype = ietype(ie)
	  do ii=1,n
	    k = index(ii,ie)
	    nen3v(ii) = inext(k)
	  end do
          call grd_write_item(nout,2,eext,itype,n,nen3v,depth)
	end do

	close(nout)

	zoom = 'zoom_' // trim(file)
	open(nout,file=zoom,status='unknown',form='formatted')
	call make_rect(nout,1.,xmin,ymin,xmax,ymax)
	call make_rect(nout,0.5,xmin,ymin,xmax,ymax)
	call make_rect(nout,0.25,xmin,ymin,xmax,ymax)
	close(nout)

	end

!*****************************************************************

	subroutine make_rect(nout,fact,xmin,ymin,xmax,ymax)

	implicit none

	integer nout
	real fact
	real xmin,ymin,xmax,ymax

	integer, save :: n = 0
	integer, save :: l = 0
	integer, save :: itype = 0
	real, parameter :: flag = -999.
	integer nodes(5)

	real dx,dy,x0,y0
	real xmin0,ymin0,xmax0,ymax0
	real depth

	depth = flag

	dx = xmax-xmin
	dy = ymax-ymin
	x0 = 0.5*(xmin+xmax)
	y0 = 0.5*(ymin+ymax)

	xmin0 = x0 - 0.5*fact*dx
	ymin0 = y0 - 0.5*fact*dy
	xmax0 = x0 + 0.5*fact*dx
	ymax0 = y0 + 0.5*fact*dy

	call grd_write_node(nout,n+1,itype,xmin0,ymin0,depth)
	call grd_write_node(nout,n+2,itype,xmin0,ymax0,depth)
	call grd_write_node(nout,n+3,itype,xmax0,ymax0,depth)
	call grd_write_node(nout,n+4,itype,xmax0,ymin0,depth)
	
	l = l + 1
	nodes = (/n+1,n+2,n+3,n+4,n+1/)
        call grd_write_item(nout,3,l,itype,5,nodes,depth)

	n = n + 4

	end

!*****************************************************************

