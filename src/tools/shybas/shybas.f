
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999,2005,2007,2009-2012,2015-2020  Georg Umgiesser
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

c revision log :
c
c 06.04.1999	ggu	completely restructured
c 04.06.1999	ggu	new statistics are computed
c 28.11.2005	ggu	new call to makehkv
c 31.05.2007	ggu	added area and volume frequency curve
c 24.08.2007	ggu	added new routine write_grd_from_bas
c 06.04.2009	ggu	read param.h
c 12.06.2009	ggu	areatr in double precision - new algorithm
c 01.03.2010	ggu	new routine basqual() to compute grid quality
c 22.03.2010	ggu	write external element number in basqual()
c 17.05.2011	ggu	changes in freqdep()
c 12.07.2011	ggu	better treatment of freqdep()
c 16.11.2011	ggu	basin.h introduced
c 10.02.2012	ggu	use angles in quality of basin (basqual)
c 30.09.2015	ggu	shybas started
c 01.10.2015	ggu	shybas nearly finished
c 02.10.2015	ggu	only basproj is missing
c 10.10.2015	ggu	changed VERS_7_3_2
c 12.10.2015	ggu	changed VERS_7_3_3
c 16.12.2015	ggu	changed VERS_7_3_16
c 17.03.2016	ggu	new routine write_depth_from_bas()
c 15.04.2016	ggu	changed VERS_7_5_8
c 25.05.2016	ggu	changed VERS_7_5_10
c 11.10.2016	ggu	changed VERS_7_5_20
c 21.03.2017	ggu	new routine to compute area/vol on area code
c 31.03.2017	ggu	changed VERS_7_5_24
c 09.10.2017	ggu	changed VERS_7_5_33
c 17.11.2017	ggu	changed VERS_7_5_37
c 07.12.2017	ggu	changed VERS_7_5_40
c 24.01.2018	ggu	changed VERS_7_5_41
c 22.02.2018	ggu	changed VERS_7_5_42
c 04.04.2018	ggu	new code for real area (m) and explicit stability
c 04.04.2018	ggu	new code for nodal partition check and write
c 13.04.2018	ggu	new routine to elab partition and write to file
c 26.04.2018	ggu	changed VERS_7_5_46
c 16.10.2018	ggu	changed VERS_7_5_50
c 25.10.2018	ccf	grid output in gr3 and msh formats
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c 13.02.2020	ggu	new routine write_regular_depth() with -reg
c 01.04.2020    ggu     new option -custom (bcustom) 
c 28.05.2020    ggu     implement bquiet and bsilent
c 12.12.2020    ggu     compute transport CFL
c 22.04.2021    ggu     initialize levels for bounds check
c 10.11.2021    ggu     avoid warning for stack size
c 16.02.2022    ggu     new call to basboxgrd() to re-create grd file from index
c 12.10.2022    ggu     new routine code_count called with -detail
c 12.01.2023    ggu     correct statistics of area also for lat/lon
c 29.01.2023    ggu     more on correct area computation (eliminated areatr)
c 10.05.2024    ggu     new routine write_basin_txt() (bbastxt)
c 03.10.2024    ggu     new call to test_fast_find()
c 21.11.2024    ggu     renamed call to some routines
c
c todo :
c
c reading grd file ngr is 1 too high
c
c****************************************************************

        program shybas

c writes information and manipulates basin

	use mod_depth
	use mod_geom
	use mod_geom_dynamic
	use evgeom
	use basin
	use levels
	use clo
	use basutil
	use shympi

	implicit none

	integer nc
	logical bwrite,bbastxt
	character*80 file

	bbastxt = .false.	! make this a command line option

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

	call basutil_init('BAS')

	bwrite = .not. bsilent

	call clo_get_file(1,file)
        if( file == ' ' ) call clo_usage
	call read_command_line_file(file,bwrite)

	call shympi_init(.false.)

c-----------------------------------------------------------------
c initialize modules
c-----------------------------------------------------------------

	call levels_init_2d(nkn,nel)

	call ev_set_verbose(bwrite)
	call set_ev

	call set_geom

	call mod_geom_dynamic_init(nkn,nel)
	call setnod_bas

	call mod_depth_init(nkn,nel)

        call makehev(hev)
        call makehkv(hkv)

c-----------------------------------------------------------------
c info on basin read
c-----------------------------------------------------------------

	if( .not. bquiet ) then
	  call bas_info
	  call basstat(bnomin,bdetail)
	  if( barea ) call basstat_area
	  call bas_stabil
	  if( bbastxt ) call write_basin_txt
	end if

        call node_test				!basic check
	if( bcheck ) call bascheck		!extra check

c-----------------------------------------------------------------
c transformations and extra info on basin
c-----------------------------------------------------------------

	if( bfreq ) call freqdep
	if( bquality ) call basqual		!grid quality
	if( bresol ) call bas_resolution	!grid resolution

	if( bcompare ) call bascompare		!compares 2 basins
	if( hsigma > 0 ) call bashsigma(hsigma)	!creates hybrid layers
	if( bfile /= ' ' ) call basbathy	!interpolate bathymetry
	if( bsmooth ) call bas_smooth		!limit and smooth
	if( binvert ) call invert_depth		!inverts depth values
	if( bbox ) call basbox			!creates box index
	if( bboxgrd ) call basboxgrd		!creates grd from index
	if( slayers /= ' ' ) call bas_layers(slayers)	!compute area/vol

	if( bcustom ) call bas_custom
	if( bfastfind ) call test_fast_find

c-----------------------------------------------------------------
c loop for interactive information on nodes and elems
c-----------------------------------------------------------------

	if( binter ) call basin_interactive

c-----------------------------------------------------------------
c write output files
c-----------------------------------------------------------------

	if( bgrd ) call write_grd_from_bas
        if( bxyz ) call write_xy('bas.xyz',nkn,ipv,xgv,ygv,hkv)
        if( bdepth ) call write_depth_from_bas
        if( breg ) call write_regular_depth(dreg)
	if( bunique ) call write_grd_with_unique_depth !for sigma levels
	if( bdelem ) call write_grd_with_elem_depth !for zeta levels
	if( bnpart ) call write_nodal_partition(bwrite)	!nodal partition
	if( lfile /= ' ' ) call bas_partition(lfile)	!creates partition file
	if( bgr3 ) call write_gr3_from_bas
	if( bmsh ) call write_msh_from_bas

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine read_command_line_file(file,bwrite)

	use basin
	use basutil

	implicit none

	character*(*) file
	logical bwrite
	logical is_grd_file,filex

	call grd_set_write(bwrite)

	if( .not. filex(file) ) then
	  write(6,*) 'file not existing: ',trim(file)
	  stop 'error stop read_command_line_file: no such file'
	else if( basin_is_basin(file) ) then
	  if( bwrite ) write(6,*) 'reading BAS file: ',trim(file)
	  call basin_read(file)
	  breadbas = .true.
	else if( is_grd_file(file) ) then
	  if( bwrite ) write(6,*) 'reading GRD file: ',trim(file)
	  call grd_read(file)
	  call grd_to_basin
	  call bas_check_spherical
	  call estimate_ngr(ngr)
	  breadbas = .false.
	else
	  write(6,*) 'Cannot read this file: ',trim(file)
	  stop 'error stop read_command_line_file: unknown format'
	end if

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine basin_interactive

	logical bnode,belem

	bnode = .true.
	belem = .true.

	do while( bnode .or. belem )

	   call nodeinfo(bnode)
	   if( .not. bnode .and. .not. belem ) exit

	   call eleminfo(belem)
	   if( .not. bnode .and. .not. belem ) exit

	end do

	end

c*******************************************************************

	subroutine nodeinfo(bnode)

c info on node number

	use mod_depth
	use basin
	use evgeom

	implicit none

	logical bnode

	integer ie,ii,in,ios
	integer kext,kint
	integer ipext,ipint
	logical bloop

	real area_node

	bnode = .false.
	bloop = .true.

        write(6,*) 

c look for node and give info

	do while( bloop )

        write(6,*) 'Enter node number : '
        read(5,'(i10)',iostat=ios) kext
	if( ios /= 0 ) kext = 0
	if( kext .gt. 0 ) bnode = .true.
	kint = ipint(kext)

	if( kext .le. 0 ) then
	   bloop = .false.
        else if(kint.le.0) then
           write(6,*) ' no node number : ',kext
           if(kext.le.nkn) then
              write(6,*) '(intern : ',kext
     +                          ,' extern : ',ipext(kext),')'
           end if
           write(6,*)
	else
           write(6,*) ' extern : ',kext,' intern : ',kint
           write(6,*) '(intern : ',kext,' extern : ',ipext(kext),')'
           write(6,*) ' (x,y)  : ',xgv(kint),ygv(kint)
           write(6,*) ' depth  : ',hkv(kint)
           write(6,*) ' area   : ',area_node(kint)
           write(6,2200)

           do ie=1,nel
              in=0
              do ii=1,3
                 if(nen3v(ii,ie).eq.kint) in=ii
              end do
              if(in.gt.0) then
                 write(6,2000)   ipev(ie)
     +                          ,(ipext(nen3v(ii,ie)),ii=1,3)
     +                          ,(hm3v(ii,ie),ii=1,3)
     +                          ,iarv(ie),hm3v(in,ie)
              end if
	   end do
           write(6,*)
	end if

	end do

	return
 2200   format(/1x,'element',8x,'nodes',14x,'depth in element'
     +          ,3x,'  area',4x,'depth of node')
 2000   format(1x,i6,2x,3i6,2x,3f8.2,2x,i5,5x,f8.2)
	end

c*****************************************************************

	subroutine eleminfo(belem)

c info on element number

	use basin

	implicit none

	logical belem

	integer ie,ii,k,ios
	integer eext,eint
	integer ipext,ieext,ieint
	logical bloop
	real area_elem

	belem = .false.
	bloop = .true.

	do while( bloop )

        write(6,*) 'Enter element number : '
        read(5,'(i10)',iostat=ios) eext
	if( ios /= 0 ) eext = 0
        if(eext.gt.0) belem = .true.
	eint = ieint(eext)

	if( eext .le. 0 ) then
	  bloop = .false.
        else if(eint.le.0) then
           write(6,*) ' no element number : ',eext
           if(eext.le.nel) then
              write(6,*) '(intern : ',eext
     +                     ,' extern : ',ieext(eext),')'
           end if
           write(6,*)
	else
           write(6,*) ' extern : ',eext,' intern : ',eint
           write(6,*) '(intern : ',eext,' extern : ',ieext(eext),')'
           write(6,*)

	   ie = eint
           do ii=1,3
              k=nen3v(ii,ie)
              write(6,*) ' (x,y) : ',xgv(k),ygv(k)
	   end do

           write(6,3200)
           write(6,3000) ieext(ie)
     +                  ,(ipext(nen3v(ii,ie)),ii=1,3)
     +                  ,(hm3v(ii,ie),ii=1,3)
     +                  ,area_elem(ie)
     +                  ,iarv(ie)
           write(6,*)
	end if

	end do

        return
 3200   format(/1x,'element',8x,'nodes',14x,'depth in element'
     +          ,3x,'   area',3x,'  area code')
 3000           format(1x,i6,2x,3i6,2x,3f8.2,2x,e10.2,2x,i5)
        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine basstat_area

c writes statistics on basin for each area code

	use basin

	implicit none

	integer ie,ia,na
	integer imin,imax
	real area,vol

c-----------------------------------------------------------------
c area code
c-----------------------------------------------------------------

	imin = iarv(1)
	imax = imin

	do ie=1,nel
	  imin = min(imin,iarv(ie))
	  imax = max(imax,iarv(ie))
	end do

c-----------------------------------------------------------------
c area
c-----------------------------------------------------------------

	write(6,'(2a)') '      ia    ntot         area       '
     +			,'         volume              depth'

	do ia=imin,imax
	  call areavol(ia,na,area,vol)
	  if( na == 0 ) cycle
	  write(6,'(2i8,3g20.8)') ia,na,area,vol,vol/area
	end do

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*****************************************************************

	subroutine areavol(ia,na,area,vol)

c computes area and volume of area code ia

	use basin

	implicit none

	integer ia	!area code to use
	integer na	!number of elements found
	real area	!area of code ia
	real vol	!volume of code ia

	logical bflag
	integer ie,ii
	double precision a,atot,vtot,h
	real hk
	real, parameter :: hflag = -999.
	real area_elem

	na = 0
	atot = 0.
        vtot = 0.

	do ie=1,nel
	  if( ia /= iarv(ie) ) cycle
	  na = na + 1
	  a = area_elem(ie)
	  atot = atot + a
          h = 0.
	  bflag = .false.
          do ii=1,3
	    hk = hm3v(ii,ie)
            h = h + hk
	    bflag = hk == hflag
          end do
	  h = h / 3.
	  if( bflag ) h = 0.
          vtot = vtot + a * h
	end do

	area = atot
	vol = vtot

	end

c*****************************************************************

	subroutine basstat(bnomin,bdetail)

c writes statistics on basin

	use basin
	use evgeom

	implicit none

	logical bnomin		!do not compute minimum distance
	logical bdetail		!write details on area and node code

	integer ie,ii,k,ii1
	integer imin,imax
	real area,amin,amax,atot
        real vtot
	real aptot,vptot
	real x(3),y(3)
	real xmin,xmax,ymin,ymax
	real dxmax,dymax
	real h
        real xx,yy,dx,dy
        real dist,distmin
	logical bflag
	real dtot,dptot
	real hk
	real, parameter :: hflag = -999.
        integer i,k1,k2,km1,km2

	real area_elem

c-----------------------------------------------------------------
c area code
c-----------------------------------------------------------------

	imin = iarv(1)
	imax = imin

	do ie=1,nel
	  imin = min(imin,iarv(ie))
	  imax = max(imax,iarv(ie))
	end do

	write(6,*) 'Area code min/max:      ',imin,imax

	if( bdetail ) call code_count(imax,nel,iarv)

	imin = iarnv(1)
	imax = imin

	do k=1,nkn
	  imin = min(imin,iarnv(k))
	  imax = max(imax,iarnv(k))
	end do

	write(6,*) 'Node code min/max:      ',imin,imax

	if( bdetail ) call code_count(imax,nkn,iarnv)

c-----------------------------------------------------------------
c node numbers
c-----------------------------------------------------------------

	imin = ipv(1)
	imax = imin

	do k=1,nkn
	  imin = min(imin,ipv(k))
	  imax = max(imax,ipv(k))
	end do

	write(6,*) 'Node number min/max:    ',imin,imax

c-----------------------------------------------------------------
c element numbers
c-----------------------------------------------------------------

	imin = ipev(1)
	imax = imin

	do ie=1,nel
	  imin = min(imin,ipev(ie))
	  imax = max(imax,ipev(ie))
	end do

	write(6,*) 'Element number min/max: ',imin,imax

c-----------------------------------------------------------------
c area
c-----------------------------------------------------------------

	amin = area_elem(1)
	amax = amin
	atot = 0.
        vtot = 0.
	aptot = 0.
        vptot = 0.

	do ie=1,nel
	  area = area_elem(ie)
	  atot = atot + area
	  amin = min(amin,area)
	  amax = max(amax,area)
          h = 0.
	  bflag = .false.
          do ii=1,3
	    hk = hm3v(ii,ie)
            h = h + hk
	    bflag = hk == hflag
          end do
	  h = h / 3.
	  if( bflag ) h = 0.
          vtot = vtot + area * h
	  if( h .gt. 0. ) then		!only positive depths
	    aptot = aptot + area
            vptot = vptot + area * h
	  end if
	end do

	write(6,*) 'Area min/max:           ',amin,amax
	write(6,*) 'Total area (all/h>0):   ',atot,aptot
	write(6,*) 'Total volume (all/h>0): ',vtot,vptot

c-----------------------------------------------------------------
c coordinates
c-----------------------------------------------------------------

	xmin = xgv(1)
	ymin = ygv(1)
	xmax = xgv(1)
	ymax = ygv(1)

	do k=1,nkn
	  xmin = min(xmin,xgv(k))
	  ymin = min(ymin,ygv(k))
	  xmax = max(xmax,xgv(k))
	  ymax = max(ymax,ygv(k))
	end do

	write(6,*) 'X-Coordinates min/max:  ',xmin,xmax
	write(6,*) 'Y-Coordinates min/max:  ',ymin,ymax

c-----------------------------------------------------------------
c size of elements
c-----------------------------------------------------------------

	dxmax = 0.
	dymax = 0.

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    x(ii) = xgv(k)
	    y(ii) = ygv(k)
	  end do
	  xmin = min(x(1),x(2),x(3))
	  xmax = max(x(1),x(2),x(3))
	  ymin = min(y(1),y(2),y(3))
	  ymax = max(y(1),y(2),y(3))
	  call compute_distance(xmin,ymin,xmax,ymax,dx,dy)
	  dxmax = max(dxmax,dx)
	  dymax = max(dymax,dy)
	end do

	write(6,*) 'Element dxmax/dymax:    ',dxmax,dymax

c-----------------------------------------------------------------
c depth
c-----------------------------------------------------------------

	amin = 999999.
	amax = -amin

	do ie=1,nel
	  h = 0
	  do ii=1,3
	    h = h + hm3v(ii,ie)
	  end do
	  h = h / 3.
	  amin = min(amin,h)
	  amax = max(amax,h)
	end do

	dtot = vtot/atot
	dptot = 0.
	if( aptot > 0. ) dptot = vptot/aptot

	write(6,*) 'Depth min/max:          ',amin,amax
	write(6,*) 'Depth average (all/h>0):',dtot,dptot

c-----------------------------------------------------------------
c minimum distance of nodes
c-----------------------------------------------------------------

        !distmin = (xmax-xmin)**2 + (ymax-ymin)**2
	call compute_distance(xmin,ymin,xmax,ymax,dx,dy)
        distmin = 2.*(dx**2+dy**2)
        km1 = 0
        km2 = 0

	if( bnomin ) then
	  distmin = 0.
	else
	  do ie=1,nel
	    do ii=1,3
	      ii1 = mod(ii,3) + 1
	      k1 = nen3v(ii,ie)
	      k2 = nen3v(ii1,ie)
	      call compute_distance(xgv(k1),ygv(k1),xgv(k2),ygv(k2),dx,dy)
              dist = dx**2 + dy**2
              if( dist .lt. distmin ) then
                km1 = k1
                km2 = k2
                distmin = dist
              end if
	    end do
	  end do
          distmin = sqrt(distmin)
	  write(6,*) 'min node distance:      ',distmin,ipv(km1),ipv(km2)
	end if

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	!write(6,*)

	end

c*******************************************************************

	subroutine bas_stabil

	use basin

	implicit none

	integer ie,ii,ii1,k1,k2,iemin,ietmin
	real dist2,distmin,dx,dy,hmax
	real dt2,dtt2,dt,dtt
	real, parameter :: high = 1.e+30
	real, parameter :: grav = 9.81

	dt2 = high
	dtt2 = high
	iemin = 0
	ietmin = 0

	do ie=1,nel
	  distmin = high
	  do ii=1,3
	    ii1 = mod(ii,3) + 1
	    k1 = nen3v(ii,ie)
	    k2 = nen3v(ii1,ie)
	    call compute_distance(xgv(k1),ygv(k1),xgv(k2),ygv(k2),dx,dy)
            dist2 = dx**2 + dy**2
	    distmin = min(distmin,dist2)
	  end do
	  !write(6,*) ie,distmin
	  hmax = maxval(hm3v(:,ie))
	  hmax = max(hmax,0.10)			!handle salt marshes
	  if( distmin/hmax < dt2 ) then		!gravity wave CFL condition
	    dt2 = distmin/hmax
	    iemin = ie
	  end if
	  if( distmin < dtt2 ) then		!transport CFL condition
	    dtt2 = distmin
	    ietmin = ie
	  end if
	end do

! dt2 = dist*dist/h = dist2/h
! dt < dist/sqrt(g*h) = sqrt(dist2/(g*h)) = sqrt(dt2/g)
! transport CFL: dt < dmin/c, with c=1m/s -> dt < sqrt(dtt2)

	dt = sqrt( dt2/grav )
	dtt = sqrt( dtt2 )

	write(6,*) 'explicit time-step                  ' 
     +				// 'time-step [s] '
     +				// '  element (ext)'
	write(6,*) 'max explicit time-step (gravity):   ',dt,ipev(iemin)
	write(6,*) 'max explicit time-step (transport): ',dtt,ipev(ietmin)

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine freqdep

c writes frequency distribution of depth

	use mod_depth
	use basin

	implicit none

	integer ndim
	parameter (ndim=10000)

	integer ie,i
	integer imax,ih
	real area,vol,z
	real hmin,hmax,dh,h,fr
	real fa,fv,fap

	double precision, allocatable :: freqa(:)
	double precision, allocatable :: freqv(:)

	real area_elem

c-----------------------------------------------------------------
c area code
c-----------------------------------------------------------------

	allocate(freqa(0:ndim))
	allocate(freqv(0:ndim))

	dh = 0.01

	hmin = 1.e+30
	hmax = -1.e+30
	do ie=1,nel
	  h = hev(ie)
	  hmin = min(hmin,h)
	  hmax = max(hmax,h)
	end do

	imax = (hmax-hmin)/dh
	if( imax .gt. ndim ) then
	  dh = (hmax-hmin)/ndim
	  imax = ndim
	end if

	write(6,*) 'computing depth frequency curve ',imax
	write(6,*) '  hmin,hmax,dh: ',hmin,hmax,dh

	do i=0,imax
	  freqa(i) = 0.
	  freqv(i) = 0.
	end do

	do ie=1,nel
	  h = hev(ie)
	  area = area_elem(ie)
	  vol = area * h
	  ih = (hmax-h)/dh
	  if( ih .lt. 0 .or. ih .gt. imax ) then
	    if( ih .lt. 0 ) then
		ih = 0
	    else if( ih .gt. imax ) then
		ih = imax
	    else		!never
	        write(6,*) ih,imax,hmin,hmax,h
	        stop 'error stop freqdep: index out of range'
	    end if
	  end if
	  do i=ih,imax
	    freqa(i) = freqa(i) + area
	    freqv(i) = freqv(i) + area * (i-ih)*dh
	  end do
	end do

	area = freqa(imax)	!total area
	vol  = freqv(imax)	!total volume

	write(6,*) 'total area/vol for freq curve: ',area,vol

	open(1,file='areavol.dat',status='unknown',form='formatted')

	do i=0,imax
	  z = -hmax + i*dh
	  fap = freqa(i) * 100. / area
	  fa = freqa(i)
	  fv = freqv(i)
	  write(1,*) z,fap,fa,fv
	end do

	close(1)

	write(6,*) 'frequency curve written to file areavol.dat'

	end

c*******************************************************************

        function areatr(ie)

c determination of area of element (not working for lat/lon)
c
c ie            number of element (internal)
c areatr        element area (return value)

	use basin

	real areatr
	integer ie

	integer ii,i1,i2,k1,k2
	double precision f,x(3),y(3)

        do ii=1,3
          k=nen3v(ii,ie)
	  x(ii) = xgv(k)
	  y(ii) = ygv(k)
        end do

	f = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))

        areatr = f / 2.D0

        end

c*******************************************************************

	subroutine make_unique_depth

c makes unique depth 

	use evgeom
	use basin

	implicit none

	integer k,ie,ii
	double precision h,area
	double precision hkv(nkn)
	double precision weight(nkn)

	hkv = 0.
	weight = 0.

	do ie=1,nel
	  area = ev(10,ie)
	  do ii=1,3
	    h = hm3v(ii,ie)
	    k = nen3v(ii,ie)
	    hkv(k) = hkv(k) + h*area
	    weight(k) = weight(k) + area
	  end do
	end do

	do k=1,nkn
	  if( weight(k) > 0. ) hkv(k) = hkv(k) / weight(k)
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    hm3v(ii,ie) =  hkv(k)
	  end do
	end do

	end

c*******************************************************************

	subroutine make_constant_depth

c makes constant depth 

	use evgeom
	use basin

	implicit none

	integer ie,ii
	double precision h,hm

	do ie=1,nel
	  hm = 0.
	  do ii=1,3
	    h = hm3v(ii,ie)
	    hm = hm + h
	  end do
	  hm = hm / 3.
	  hm3v(:,ie) = hm
	end do

	end

c*******************************************************************

	subroutine round_depth(hd)

c rounds depth values to nearest given value

	use basin

	implicit none

	real hd

	integer ie,ii,i
	real h

	do ie=1,nel
	  do ii=1,3
	    h = hm3v(ii,ie)
	    i = nint(h/hd)
	    if( i <= 0 ) i = 1
	    h = i * hd
	    hm3v(ii,ie) = h
	  end do
	end do

	end

c*******************************************************************

	subroutine invert_depth

c inverts depth values

	use basin
	use basutil

	implicit none

	hm3v = -hm3v

	bgrd = .true.		!ensure writing of file

	end

c*******************************************************************

	subroutine write_grd_with_unique_depth

c writes grd file extracting info from bas file

	implicit none

        write(6,*) 'making unique depth...'
	call make_unique_depth

        call basin_to_grd
        call grd_write('basunique.grd')

        write(6,*) 'The basin has been written to basunique.grd'

	end

c*******************************************************************

	subroutine write_grd_with_elem_depth

c writes grd file extracting info from bas file

	implicit none

        write(6,*) 'making constant elem depth...'
	call make_constant_depth
	!call round_depth(2.)

        call basin_to_grd

        call grd_write('baselem.grd')
        write(6,*) 'The basin has been written to baselem.grd'

	end

c*******************************************************************

	subroutine write_grd_from_bas

c writes grd file extracting info from bas file

	implicit none

        call basin_to_grd

        call grd_write('bas.grd')
        write(6,*) 'The basin has been written to bas.grd'

	end

c*******************************************************************

	subroutine write_gr3_from_bas

c writes grid in gr3 format from bas file
c the gr3 format is used by WWMIII wave model

	use basutil

	implicit none

	if (.not. breadbas) then
	  write(*,*)
	  write(*,*)'You need a bas file for creating a gr3 grid'
	  stop 'error: write_gr3_from_bas'
	end if

        write(6,*) 'making unique depth...'
	call make_unique_depth

        call basin_to_grd
        call gr3_write('bas.gr3')

        write(6,*) 'The basin has been written to bas.gr3'
        write(6,*) 'The boundary file been written to bas_bnd.gr3'

	end

c*******************************************************************

	subroutine write_msh_from_bas

c writes grid in msh (gmsh v.2) format from bas file
c the msh format is used by WWM3 wave model

	use basutil

	implicit none

	if (.not. breadbas) then
	  write(*,*)
	  write(*,*)'You need a bas file for creating a msh grid'
	  stop 'error: write_msh_from_bas'
	end if

        write(6,*) 'making unique depth...'
	call make_unique_depth

        call basin_to_grd
        call msh_write('bas.msh')

        write(6,*) 'The basin has been written to bas.msh'

	end

c*******************************************************************

	subroutine basqual

c writes statistics on grid quality

	use evgeom
	use basin

	implicit none

	real areav(nkn)

	integer ie,ii,k
	integer ia,ic,ilow,ihigh
	integer imin,imax
	real area,amin,amax,atot
        real vtot
	real aptot,vptot
	real x(3),y(3)
	real xmin,xmax,ymin,ymax
	real dxmax,dymax
	real h
        real xx,yy
        real dist,distmin
	real fmax,fmin,a,f
        integer i,k1,k2
	integer iemin,iemax

	integer iangle(0:18)

	integer ieext,ipext

c-----------------------------------------------------------------
c initialization
c-----------------------------------------------------------------

	do k=1,nkn
	  areav(k) = 0.
	end do

c-----------------------------------------------------------------
c compute fraction
c-----------------------------------------------------------------

	do ie=1,nel
	  area = 4.*ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    areav(k) = areav(k) + area
	  end do
	end do

	fmax = 0
	fmin = 10.
	do ie=1,nel
	  area = 12.*ev(10,ie)
	  amax = 0.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    a = areav(k)
	    amax = max(amax,a)
	  end do
	  f = amax/area
	  if( f .gt. 10. ) then
	    write(6,*) 'bad quality of element: ',ie,ieext(ie),f
	  end if
	  if( f .gt. fmax ) then
	    fmax = f
	    iemax = ie
	  end if
	  if( f .lt. fmin ) then
	    fmin = f
	    iemin = ie
	  end if
	end do

	write(6,*) 'Grid quality: (internal/external element number)'
	write(6,*) '   minimum: ',iemin,ieext(iemin),fmin
	write(6,*) '   maximum: ',iemax,ieext(iemax),fmax
	write(6,*)

c-----------------------------------------------------------------
c compute angles
c-----------------------------------------------------------------

	do i=0,18
	  iangle(i) = 0
	end do

	do ie=1,nel
	  !area = 4.*ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    a = ev(10+ii,ie)
	    ia = a / 10.
	    if( ia .ge. 16 ) then
		write(6,*) 'bad angle found: ',k,ipext(k),a
	    end if
	    iangle(ia) = iangle(ia) + 1
	  end do
	end do

	do i=0,18
	  ic = iangle(i)
	  ilow = i * 10
	  ihigh = (i+1) * 10
	  write(6,*) 'angles between ',ilow,ihigh,ic
	end do

	write(6,*)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

        subroutine write_xy(nfile,nkn,ipv,xgv,ygv,hkv)

c writes xy as data to file (format is gis format)

        implicit none

        character*(*) nfile
        integer nkn
        integer ipv(nkn)
        real xgv(nkn)
        real ygv(nkn)
        real hkv(nkn)

        integer k

        open(1,file=nfile,status='unknown',form='formatted')
        write(1,*) 0,nkn,5,'0'
        do k=1,nkn
          write(1,*) k,xgv(k),ygv(k),1,hkv(k)
        end do
        close(1)
        write(6,*) 'The coordinates have been written to ',nfile

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

        subroutine bascheck

c writes statistics on basin

        use evgeom
        use basin

        implicit none

        integer ie,ii,k,i
        integer imin,imax
        real area,amin,amax
        real x(3),y(3)
        real xmin,xmax,ymin,ymax
        real dxmax,dymax
        real h,w,eps
        integer iang,ic

	integer, parameter :: ndim = 20
        integer icount(ndim)
        integer count(nkn)

        integer ipext

        eps = 1.e-5

c-----------------------------------------------------------------
c area code
c-----------------------------------------------------------------

	write(6,*) 'checking basin...'
	write(6,*) '(node numbers are external)'

        do k=1,nkn
          count(k) = 0
        end do

        do ie=1,nel
          do ii=1,3
            k = nen3v(ii,ie)
            count(k) = count(k) + 1
          end do
        end do

        do k=1,nkn
          if( count(k) .le. 1 ) then
            write(6,*) 'low count for node ',ipext(k),' : ',count(k)
          end if
        end do

c-----------------------------------------------------------------
c angle
c-----------------------------------------------------------------

	icount = 0

        do ie=1,nel
          do ii=1,3
            k = nen3v(ii,ie)
            w = ev(10+ii,ie)
            iang = (w-90.)/10. + 1. - eps
	    if( iang > ndim ) stop 'error stop bascheck: iang'
            if( iang .gt. 0 ) then
              icount(iang) = icount(iang) + 1
c             write(6,*) k,ii,w,iang
            end if
          end do
        end do

        do i=1,ndim
          ic = icount(i)
          if( ic .gt. 0 ) then
            write(6,*) 'angle > ',80+10*i,' : ',ic
          end if
        end do

        do ie=1,nel
          do ii=1,3
            k = nen3v(ii,ie)
            w = ev(10+ii,ie)
            if( w .gt. 120. ) then
              write(6,*) 'big angle ',ipext(k),' : ',w
            end if
          end do
        end do

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

        end

c*******************************************************************

	subroutine bascompare

c compares two basins and writes delta depths to file

	use basin
	use evgeom
	use clo

	implicit none

	integer nel_aux
	real hm3v_aux(3,nel)
	character*80 file

	nel_aux = nel
	hm3v_aux = hm3v

	call clo_get_file(2,file)
	if( file == ' ' ) then
	  write(6,*) 'for -compare we need two bas files'
	  stop 'error stop bascomp: missing second file'
	end if
	call basin_read(file)

	if( nel /= nel_aux ) then
	  write(6,*) 'dimension of basins incompatible: ',nel,nel_aux
	  stop 'error stop bascomp: nel'
	end if

	!call ev_init(nel)
	call set_ev

	hm3v = hm3v - hm3v_aux

        call basin_to_grd
        call grd_write('bascomp.grd')
        write(6,*) 'The basin has been written to bascomp.grd'

	end

c*******************************************************************

        subroutine node_test

        use basin

        implicit none

        logical bstop
        integer ie,ii,iii,k,k1

        bstop = .false.

        !write(6,*) 'node_test ... ',nel,nkn

        do ie=1,nel
          do ii=1,3
            k = nen3v(ii,ie)
            if( k .le. 0 ) then
                write(6,*) ie,ii,k
                bstop = .true.
            end if
            iii = mod(ii,3) + 1
            k1 = nen3v(iii,ie)
            if( k .eq. k1 ) then
                write(6,*) ie,(nen3v(iii,ie),iii=1,3)
                bstop = .true.
            end if
          end do
        end do

        !write(6,*) 'end of node_test ... '

        if( bstop ) stop 'error stop node_test: errors'

        end

c*******************************************************************

        subroutine write_depth_from_bas

! writes depth values of elements to file

        use basin

        implicit none

        integer ie,ii,k
        double precision xm,ym,hm

	open(1,file='depth.grd',status='unknown',form='formatted')
	open(2,file='depth.xyz',status='unknown',form='formatted')

        do ie=1,nel

          xm = 0.
          ym = 0.
          hm = 0.
          do ii=1,3
            k = nen3v(ii,ie)
            xm = xm + xgv(k)
            ym = ym + ygv(k)
            hm = hm + hm3v(ii,ie)
          end do
          xm = xm / 3.
          ym = ym / 3.
          hm = hm / 3.

          write(2,*) xm,ym,hm
          write(1,1000) 1,ie,0,xm,ym,hm
 1000     format(i1,2i10,3g18.8)

        end do

	close(1)
	close(2)

        write(6,*) 'The depth values have been written to depth.grd/xyz'

        end

c*******************************************************************

	subroutine write_nodal_partition(bwrite)

	use basin

	implicit none

	logical bwrite

	integer ierr1,ierr2
	integer ie,k,ii,ic,nc,n
	integer, allocatable :: nic(:)
	integer icolor(nkn)

        if( bwrite ) write(6,*) 'making nodal partition...'

	nc = maxval(iarnv)
	allocate(nic(0:nc))
	nic = 0
	do k=1,nkn
	  ic = iarnv(k)
	  nic(ic) = nic(ic) + 1
	end do
	do ic=0,nc
	  n = nic(ic)
	  if( n > 0 ) then
	    if( bwrite ) then
	      write(6,'(2i10,f10.2,a)') ic,n,(100.*n)/nkn,' %'
	    end if
	  end if
	end do
	deallocate(nic)

	do ie=1,nel
	  k = nen3v(1,ie)
	  ic = iarnv(k)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( iarnv(k) /= ic ) ic = 0
	  end do
	  iarv(ie) = ic
	end do

	call grd_set_write(bwrite)
        call basin_to_grd

        call grd_write('npart.grd')
	if( bwrite ) write(6,*) 'The basin has been written to npart.grd'

	return

        call link_set_stop(.false.)     !do not stop after error
        call link_set_write(.false.)    !do not write error

	call check_bas_connectivity(ierr1)
	call check_bas_connections(ierr2)

	if( ierr1 /= 0 .or. ierr2 /= 0 ) then
	  write(6,*) 'there were errors in link structure:'
	  write(6,*) 'connectivity: ',ierr1
	  write(6,*) 'connections:  ',ierr2
	end if

        end

c*******************************************************************

        subroutine setnod_bas

c sets (dynamic) array inodv
c
c inodv 
c        0: internal node  
c       >0: open boundary node
c       -1: boundary node  
c       -2: out of system
c
c if open boundary node, inodv(k) is number of boundary (ggu 15.11.2001)

        use mod_geom_dynamic
        use evgeom
        use basin

        implicit none

        double precision, parameter :: winmax = 359.99
        integer ie,ii,k,n
        integer ibc,ibtyp
        integer nbc,ndry
        double precision winkv(nkn)

c initialize array to hold angles

        ndry = 0
        winkv = 0.

c sum angles

        do ie=1,nel
          if(iwegv(ie).eq.0) then !element is in system
            do ii=1,3
              k=nen3v(ii,ie)
              winkv(k)=winkv(k)+ev(10+ii,ie)
            end do
          else
            ndry = ndry + 1
          end if
        end do

c set up inodv

        do k=1,nkn
          if(winkv(k).gt.winmax) then     !internal node
            inodv(k)=0
          else if(winkv(k).eq.0.) then  !out of system
            inodv(k)=-2
          else                          !boundary node
            inodv(k)=-1
          end if
        end do

        end

c****************************************************************

	subroutine sort_basin

! try sort on basin points - only for testing

	use basin

	implicit none

	integer k
	integer ip1,ip2,ia,ju1,ju2,ju,iaa
	integer nerror,nequal
	real x,y,ra
	logical bwrite
	real rk(nkn)
	integer index(nkn)

	integer locater

	write(6,*) 'testing basin...',nkn

	bwrite = .true.
	bwrite = .false.
	nerror = 0
	nequal = 0

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  rk(k) = x + y
	end do

	call isortr(nkn,rk,index)
	!call sort_invert(nkn,index)

	do k=2,nkn
	  ip1 = index(k-1)
	  ip2 = index(k)
	  if( rk(ip1) > rk(ip2) ) then
	    write(6,*) '*** error ',k,rk(ip1),rk(ip2)
	    write(6,*) xgv(ip1),xgv(ip2),ygv(ip1),ygv(ip2)
	    nerror = nerror + 1
	  else if( rk(ip1) == rk(ip2) ) then
	    if( bwrite ) then
	      write(6,*) 'equal ',k,rk(ip1),rk(ip2)
	      write(6,*) xgv(ip1),xgv(ip2),ygv(ip1),ygv(ip2)
	    end if
	    nequal = nequal + 1
	  end if
	end do

	write(6,*) 'sorted basin tested...',nkn,nerror,nequal
	if( nerror>0 ) stop 'error stop sort_basin: sort coords'

	write(6,*) 'locate nodes...'

	nerror = 0
	nequal = 0

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  ra = x + y
	  call locater_all(nkn,rk,index,ra,ju1,ju2)
	  if( ju1 == 0 ) then
	    ia = 0
	    write(6,*) '*** no entry found: ',k,ra
	  else if( ju1 /= ju2 ) then
	    if( bwrite ) write(6,*) k,ju1,ju2,ju2-ju1+1
	    nequal = nequal + 1
	    ia = 0
	    do ju=ju1,ju2
	      iaa = index(ju)
	      if( xgv(iaa) == x .and. ygv(iaa) == y ) ia = iaa
	    end do
	  else
	    ia = index(ju1)
	  end if
	  if( ia /= k ) then
	    write(6,*) '*** cannot find node: ',k,ia,ju1,ju2
	    nerror = nerror + 1
	  end if
	end do

	write(6,*) 'all nodes found...',nkn,nerror,nequal

	end

c****************************************************************

        subroutine write_regular_depth(dreg)

! interpolates bathymetry on regular grid

        use basin
        use mod_depth

	implicit none

	real dreg

	integer ie,k
	integer nx,ny
	integer iunit,iformat,datetime(2),np
	double precision dtime
	integer ilhkv(1)
	real xmin,ymin,xmax,ymax,dx,dy
	real hlv(1),hd(1)
	real regpar(7)
	character*80 file,string
	real, allocatable :: hreg(:)
	real, save :: flag = -999.

	if( dreg <= 0. ) return

!-------------------------------------------------------------
! find min/max for regular grid
!-------------------------------------------------------------

	call make_reg_box(dreg,regpar)
	call getreg(regpar,nx,ny,xmin,ymin,dx,dy,flag)

	write(6,*) 'xmin,xmax: ',xmin,xmax
	write(6,*) 'ymin,ymax: ',ymin,ymax
	write(6,*) 'dx,dy: ',dx,dy
	write(6,*) 'nx,ny: ',nx,ny

!-------------------------------------------------------------
! interpolate on regular grid
!-------------------------------------------------------------

	allocate( hreg(nx*ny) )

	if( dx /= dy ) stop 'error stop write_regular_depth: dx/=dy'
	if( dx /= dreg ) stop 'error stop write_regular_depth: dx/=dreg'
	call setgeo(xmin,ymin,dx,dy,flag)

	call av2am(hkv,hreg,nx,ny)

!-------------------------------------------------------------
! write fem file
!-------------------------------------------------------------

	iunit = 1
	file = 'regbathy.fem'
	string = 'bathymetry'
	np = nx*ny

	open(iunit,file=file,status='unknown',form='formatted')
	call write_regular_2d_1var_record(iunit,string,regpar,np,hreg)
	close(iunit)

	write(6,*) 'The regular bathymetry has been written to file '
     +			,trim(file)

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

c****************************************************************

	subroutine code_count(imax,n,iarray)

	implicit none

	integer imax,n
	integer iarray(n)

	integer i,ia,ic
	integer, allocatable :: count(:)

	allocate(count(0:imax))
	count = 0

	write(6,*) '  details:'
	write(6,*) '       code       count'

	do i=1,n
	  ia = iarray(i)
	  count(ia) = count(ia) + 1
	end do

	do ia=0,imax
	  ic = count(ia)
	  write(6,*) ia,ic
	end do

	end

c****************************************************************

	subroutine write_basin_txt

! writes basin in txt form - sorted by external numbers

	use basin
	use mod_depth
	use mod_sort

	implicit none

	integer iu,nmax
	integer ie,ii,k,i,ind,kext,iext
	integer, allocatable :: index(:)
	double precision darea_n(nkn)
	real area_n(nkn)
	double precision area

	real area_elem
	integer ipint,ieint

	iu = 567
	nmax = max(nkn,nel)
	allocate(index(nmax))

!-------------------------------------------------------
! compute area of nodes
!-------------------------------------------------------

	darea_n = 0.

	do ie=1,nel
	  area = area_elem(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    darea_n(k) = darea_n(k) + area
	  end do
	end do
	    
	area_n = darea_n / 3.

!-------------------------------------------------------
! write general info
!-------------------------------------------------------

	write(iu,*) nkn,nel,ngr,mbw

!-------------------------------------------------------
! write info on nodes
!-------------------------------------------------------

	call sort(nkn,ipv,index)

	write(iu,*) 'nodes: ',nkn

	do i=1,nkn
	  ind = index(i)
	  kext = ipv(ind)
	  k = ipint(kext)
	  write(iu,2300) i,kext,xgv(k),ygv(k),hkv(k),area_n(k)
	end do

!-------------------------------------------------------
! write info on elems
!-------------------------------------------------------

	call sort(nel,ipev,index)

	write(iu,*) 'elems: ',nel

	do i=1,nel
	  ind = index(i)
	  iext = ipev(ind)
	  ie = ieint(iext)
	  write(iu,2400) i,iext,hm3v(:,ie),area_elem(ie)
	end do

!-------------------------------------------------------
! end of routine
!-------------------------------------------------------

	write(6,*) 'info on basin written to unit ',iu

	return
 2300   format(i6,i8,2f18.6,f10.2,e14.6)
 2400   format(i6,i8,3f10.2,e14.6)
 2200   format(/1x,'element',8x,'nodes',14x,'depth in element'
     +          ,3x,'  area',4x,'depth of node')
 2000   format(1x,i6,2x,3i6,2x,3f8.2,2x,i5,5x,f8.2)
	end

c****************************************************************

