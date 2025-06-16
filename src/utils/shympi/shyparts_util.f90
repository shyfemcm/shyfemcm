
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2019  Georg Umgiesser
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

! revision log :
!
! 19.05.2020	ccf	started from scratch
! 24.05.2020	ggu	debug option added
! 28.05.2020	ggu	some more informational messages
! 15.07.2020	ggu	bug fix for counting elements
! 22.04.2021	ggu	resolve bound check error (not yet finished)
! 12.04.2022	ggu	preapred for online partitioning
! 10.04.2024	ggu	new routine write_single_nodes()
! 18.11.2024	ggu	writes partition.np.txt
! 21.11.2024	ggu	use make_name_with_number() to create file names
! 23.11.2024	ggu	new routine info_partition_quality()
! 29.11.2024	ggu	compute number of non contiguous areas
! 30.11.2024	ggu	use color as aux array for flood_fill
! 03.01.2025	ggu	new routine check_tripple_points()
!
!****************************************************************

	subroutine info_partition(nparts,area_node,pquality)

! write partition information to terminal

	use basin

	implicit none

	integer, intent(in) :: nparts
	integer, intent(in) :: area_node(nkn)
	real, intent(out) :: pquality

	integer ic,k,ia
	integer netot,neint
	integer iareas
	integer min,max
	character*80 name

        integer, allocatable  :: nc(:)    !total number of nodes with color ic
        integer, allocatable  :: ne(:)    !total number of elements of ic
        integer, allocatable  :: ni(:)    !internal number of elements of ic
        integer, allocatable  :: na(:)    !number of areas
        integer, allocatable  :: color(:) !aux array for area_node

	write(6,*) 'writing information on partion to terminal...'

        allocate(nc(0:nparts))
        allocate(ne(0:nparts))
        allocate(ni(0:nparts))
        allocate(na(0:nparts))
	allocate(color(nkn))

        nc = 0
        ne = 0
        ni = 0
        na = 0
	netot = 0
	neint = 0
	min = minval(area_node)
	max = maxval(area_node)
	write(6,*) 'min/max: ',min,max
        if( min < 0 .or. max > nparts ) then
          write(6,*) 'nparts: ',nparts
          write(6,*) 'min/max: ',min,max
          stop 'error stop info_partition: internal error (1)'
        end if
	do k=1,nkn
          ic = area_node(k)
          nc(ic) = nc(ic) + 1
	end do
	color = area_node	!to not change area_node
	do ic=min,max
	  call count_elements(nkn,nel,nen3v,ic,color,netot,neint)
	  iareas = -1
	  call check_part_color(ic,nkn,color,iareas)
	  if( iareas == 0 ) iareas = 1
	  na(ic) = iareas
	  ne(ic) = netot
	  ni(ic) = neint
        end do
        write(6,*) 
        write(6,*) 'total number of nodes: ',nkn
        write(6,*) 'total number of elems: ',nel
        write(6,*) 
        write(6,*) 'Information on domains: ',nparts
        write(6,*) 
        write(6,*) '   domain      area     nodes   percent' &
     &				// '  elements     ghost' &
     &				// '   percent   areas'
        do ic=min,max
	  ia = ic
	  if( min == 0 ) ia = ic + 1
          write(6,'(3i10,f10.2,2i10,f10.2,i8)')  &
     &		 ia-1,ia,nc(ic),(100.*nc(ic))/nkn &
     &		,ne(ic),ne(ic)-ni(ic),(100.*(ne(ic)-ni(ic)))/ne(ic),na(ic)
        end do
        write(6,*) 

	call make_name_with_number('partition',nparts,'txt',name)
	open(1,file=name,form='formatted',status='unknown')
	write(1,'(a)') '        np       nkn       nel'
	write(1,'(3i10)') nparts,nkn,nel
        write(1,*) '   domain      area     nodes   percent' &
     &				// '  elements     ghost' &
     &				// '   percent'
        do ic=min,max
	  ia = ic
	  if( min == 0 ) ia = ic + 1
          write(1,'(3i10,f10.2,2i10,f10.2)')  &
     &		 ia-1,ia,nc(ic),(100.*nc(ic))/nkn &
     &		,ne(ic),ne(ic)-ni(ic),(100.*(ne(ic)-ni(ic)))/ne(ic)
        end do
	close(1)

	call check_tripple_points(area_node)

	call info_partition_quality(nparts,ne,ni,pquality)

	write(6,*) 'partition quality: ',pquality

	end

!****************************************************************

	subroutine check_tripple_points(area_node)

	use basin

	implicit none

	integer area_node(nkn)

	integer ie,ii,k
	integer ic,ids
	integer id_elem(3)
	integer icount(0:3)
	integer, allocatable :: ides(:)

	allocate(ides(nel))
	ides = 0

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ic = area_node(k)
	    id_elem(ii) = ic
	  end do
	  ids = 0
	  if( id_elem(1) /= id_elem(2) ) ids = ids + 1
	  if( id_elem(2) /= id_elem(3) ) ids = ids + 1
	  if( id_elem(3) /= id_elem(1) ) ids = ids + 1
	  if( ids == 0 ) ids = 1
	  ides(ie) = ids
	end do

	icount = 0
	do ie=1,nel
	  ids = ides(ie)
	  icount(ids) = icount(ids) + 1
	end do

	write(6,*) 'total ghost elements:'
	write(6,*) '(elements with 3 colors are tripple points)'
	write(6,'(a)') '     colors   elements'

	do ids=1,3
	  ic = icount(ids)
	  write(6,'(2i11)') ids,ic
	end do
	ic = icount(3)
	if( ic == 0 ) then
	  write(6,*) 'no tripple points found'
	else
	  write(6,*) 'total number of tripple points: ',ic
	end if
	write(6,*)

	end

!****************************************************************

	subroutine info_partition_quality(np,ne,ni,pquality)

	use basin

	implicit none

	integer np
	integer ne(0:np)
	integer ni(0:np)
	real pquality

	integer ic
	real gmax,g,gtheo
	real base,logfact

	gmax = 0
	do ic=0,np
	  if( ne(ic) == 0 ) cycle
	  g = (ne(ic)-ni(ic))/real(ne(ic))
	  gmax = max(gmax,g)
	end do

        base = 2
        logfact = log(real(np))/log(real(base))
        gtheo = logfact*sqrt(real(nkn/np))/real(nel/np)

	pquality = gmax / gtheo

	end

!****************************************************************

	subroutine write_partition_to_grd(grdname,bdebug &
     &			,nparts,npart,epart)

! write grd files

	use basin
	use grd

	implicit none

	character*(*) grdname
	logical bdebug
	integer nparts
	integer npart(nkn)
	integer epart(nel)

	logical bldebug
	integer i,nmax,ntot,ic
	character*80 name

	bldebug = .true.
	bldebug = .false.

	write(6,*) 'writing grd-file...'
	call basin_to_grd
	call grd_set_write(.false.)

        ianv = npart
	call make_name_with_number(grdname,nparts,'node.grd',name)
        call grd_write(name)
	write(6,*) 'Grid with partition on nodes in file: ',trim(name)

	if( bldebug ) then
	  nmax = maxval(epart)
	  ntot = 0
	  write(6,*) 'max colors in epart:'
	  do i=-1,nmax
	    ic = count( epart == i )
	    write(6,*) i,ic
	    ntot = ntot + ic
	  end do
	  write(6,*) 'total elements: ',ntot,nel
	end if

        iaev = epart
	call make_name_with_number(grdname,nparts,'elem.grd',name)
        call grd_write(name)
	write(6,*) 'Grid with partition on elements in file: ',trim(name) 

	if( bdebug ) then
	  call grd_write_debug(grdname,nparts,npart)
	end if

        end

!****************************************************************

	subroutine grd_write_debug(grdname,nparts,npart)

	use basin
        use grd

	implicit none

	character*(*) grdname
	integer nparts
	integer npart(nkn)

	logical bhasnode
	integer i,k,ie,ii
	character*80 name,pre,post,numb
	integer epart(nel)

	call make_name_with_number(grdname,nparts,'',pre)

	write(6,*) 'writing debug grd files...'

	do i=1,nparts
	  !write(numb,'(i3)') i
          !numb = adjustl(numb)
          !name = trim(pre)//trim(numb)//trim(post)
	  call make_name_with_number(pre,i,'grd',name)
	  write(6,*) 'writing debug file: ',trim(name)
          ianv = npart
	  where( ianv /= i ) ianv = -1
	  iaev = 0
	  do ie=1,nel
	    bhasnode = .false.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( ianv(k) == i ) bhasnode = .true.
	    end do
	    if( bhasnode ) iaev(ie) = i
	  end do
          call grd_write(name)
	end do

	end

!*******************************************************************

        subroutine write_single_nodes(name,n,ks,itype)

	use basin

	implicit none

	character*(*) name
	integer n
	integer ks(n)
	integer itype

	integer, save :: iu = 100
	integer i,k

	write(6,*) 'writing single nodes grd files...'

        open(iu,file=name,status='unknown',form='formatted')

	do i=1,n
	  k = ks(i)
	  write(iu,'(i1,2i10,2f16.8)') 1,i,itype,xgv(k),ygv(k)
	end do

	close(iu)

	end

!*******************************************************************

