
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

! gr3 and msh routines - read grid files in gr3 and msh formats
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
! 18.10.2023	ggu	extracted gr3 and msh routines to this file
!
!**********************************************************

	subroutine gr3_write(file)

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
	write(iu,*) 'gr3_write: ',nk,ne,nl,nne,nnl
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

	call gr3_write_grid( &
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

	subroutine msh_write(file)

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
	write(iu,*) 'msh_write: ',nk,ne,nl,nne,nnl
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

	call msh_write_grid( &
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

	subroutine gr3_write_grid( &
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

	integer nout,nout1
	integer i
	integer n,ib,nmax
	integer k,ie,il
	integer ia,ieext,ilext
	real depth
	logical bsort,bextern
	integer, allocatable :: ipdex(:)
	integer, allocatable :: nextern(:)

	integer ifileo

	bsort = .false.
	bsort = .true.		!sort nodes on external numbering

	bextern = .true.	!use external nodes
	bextern = .false.

	!if( .not. bextern ) bsort = .false.

	nout = ifileo(1,file,'formatted','unknown')
	if( nout .le. 0 ) goto 99
	nout1 = ifileo(1,'bas_bnd.gr3','formatted','unknown')
	if( nout1 .le. 0 ) goto 99

	nmax = max(nk,ne,nl)
	allocate(ipdex(nmax))
	do i=1,nmax
	  ipdex(i) = i
	end do

	allocate(nextern(nk))
	do k=1,nk
	  if( bextern ) then
	    nextern(k) = ippnv(k)		!use external numbering
	  else
	    nextern(k) = k			!use internal numbering
	  end if
	end do

	write(nout,*)file
	write(nout,*)ne,nk
	write(nout1,*)file
	write(nout1,*)ne,nk

	call isort(nk,ippnv,ipdex)

	do k=1,nk
	  depth = hhnv(k)
	  write(nout,1000) nextern(k),xv(k),yv(k),depth
	  depth = 0.
	  write(nout1,1000) nextern(k),xv(k),yv(k),depth
	end do

	!write(nout,*)

	if( bsort ) call isort(ne,ippev,ipdex)

	do ie=1,ne
	  ieext = ie
	  if( bextern ) ieext = ippev(ie)
	  ia = iaev(ie)
	  depth = hhev(ie)
	  n = ipntev(ie) - ipntev(ie-1)
	  ib = ipntev(ie-1)
	  write(nout,2000) ieext,n,(nextern(inodev(ib+k)),k=1,n)
	  write(nout1,2000) ieext,n,(nextern(inodev(ib+k)),k=1,n)
	end do

	close(nout)

	deallocate(ipdex)
	deallocate(nextern)

	return
 1000	format(i10,3f16.8)
 2000	format(5i10)
   99	continue
	write(6,*) 'error opening output file'
	write(6,*) file
	stop 'error stop gr3_write_grid: cannot open file'
	end

!**********************************************************

	subroutine msh_write_grid( &
     &			 file &
     &			,nco,nk,ne,nl,nne,nnl &
     &			,ippnv,ippev,ipplv &
     &			,ianv,iaev,ialv &
     &			,hhnv,hhev,hhlv &
     &			,xv,yv &
     &                  ,ipntev,inodev &
     &                  ,ipntlv,inodlv &
     &			)

! writes grd file in msh format (gmsh 2.2)

        use mod_geom_dynamic

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

	integer nb		!number of boundary nodes
	integer nit		!number of items (boundary nodes + elem)
	integer nout
	integer i,idit
	integer n,ib,nmax
	integer k,ie,il
	integer ia,flag,ilext
	real depth
	logical bsort,bextern,bbound
	integer, allocatable :: ipdex(:)
	integer, allocatable :: nextern(:)

	integer ifileo

	bbound = .false.
	bbound = .true.		!include boundary nodes
	bextern = .true.	!use external nodes
	bextern = .false.

	!if( .not. bextern ) bsort = .false.

	nout = ifileo(1,file,'formatted','unknown')
	if( nout .le. 0 ) goto 99

	nmax = max(nk,ne,nl)
	allocate(ipdex(nmax))
	do i=1,nmax
	  ipdex(i) = i
	end do

	allocate(nextern(nk))
	do k=1,nk
	  if( bextern ) then
	    nextern(k) = ippnv(k)		!use external numbering
	  else
	    nextern(k) = k			!use internal numbering
	  end if
	end do

	! write header
	write(nout,'((a14))')'$MeshFormat   '
	write(nout,'((a14))')'2.2 0 8       '
	write(nout,'((a14))')'$EndMeshFormat'
	write(nout,'((a14))')'$Nodes        '
	write(nout,*)nk

	! write nodes
	nb = 0
	do k=1,nk
	  depth = hhnv(k)
	  write(nout,1000) nextern(k),xv(k),yv(k),depth
	  if (inodv(k) < 0 ) nb = nb + 1
	end do
	write(nout,'((a14))')'$EndNodes     '

	! write items (boundary nodes + elements)
	nit = ne
	if ( bbound ) nit = nit + nb

	write(nout,'((a14))')'$Elements     '
	write(nout,*)nit

	! write boundary nodes
        idit = 0
        flag = 15		!for boundary node
	if ( bbound ) then
  	  do k=1,nk
	    if (inodv(k) < 0 ) then
              idit = idit + 1
              ia = ianv(k)
	      write(nout,2000) idit,flag,1,ia,k
            end if
	  end do
        end if

	! write elements
        flag = 2		!for 3-node triangle
	do ie=1,ne
          idit = idit + 1
	  ia = iaev(ie)
	  depth = hhev(ie)
	  n = ipntev(ie) - ipntev(ie-1)
	  ib = ipntev(ie-1)
	  write(nout,2000) idit,flag,1,ia,(nextern(inodev(ib+k)),k=1,n)
	end do

	write(nout,'((a12))')'$EndElements  '
	close(nout)

	deallocate(ipdex)
	deallocate(nextern)

	return
 1000	format(i10,3f16.8)
 2000	format(7i10)
   99	continue
	write(6,*) 'error opening output file'
	write(6,*) file
	stop 'error stop msh_write_grid: cannot open file'
	end

!**********************************************************
!**********************************************************
!**********************************************************

