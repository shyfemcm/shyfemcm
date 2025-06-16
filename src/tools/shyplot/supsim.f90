
!--------------------------------------------------------------------------
!
!    Copyright (C) 2002-2005,2007-2019  Georg Umgiesser
!    Copyright (C) 2011  Debora Bellafiore
!    Copyright (C) 2012  Christian Ferrarin
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

!  utility routines for plotsim
! 
!  contents :
! 
!  subroutine plofem(type,ivarin)
!  subroutine plonos(type)
!  subroutine mkvarline(ivar,line)       makes variable description in nos file
!  subroutine plozet
!  subroutine plobar			plots barene
!  subroutine plosim(bvel)
!  subroutine ploval(nkn,parray,title)
!  subroutine plozbar
!  subroutine plovel(ivel)
!  subroutine plo2vel(bvel)
!  subroutine plo3vel(bvel)
!  subroutine por2vel(n,up,vp,uv,vv,ht)
!  subroutine plobas
! 
!  subroutine elkdep(nkn,hkv)
!  subroutine eltype0(nel,iarv)
!  subroutine eldepth(nel,hm3v,auxv,color,title)
!  subroutine eltype(nel,iarv,color,title)
! 
!  subroutine plotel( ie , color )
!  subroutine basplc( nel , color , line , descrr )
! 
!  subroutine bnd2val(a,val)			sets values on boundary to val
!  subroutine moduv(u,v,uv,n,uvmax)		computes modulus and maximum
!  subroutine normuv(u,v,vv,n)			normalize horizontal velocities
!  subroutine apply_dry_mask(bmask,value,n,flag)	aplies mask to nodes
! 
!  revision log :
! 
!  30.01.2002	ggu	better debug info, allow for vertical plotting
!  30.01.2002	ggu	new routines mkvarline, aspecial, traxy, rdspecial
!  05.10.2003	ggu	changes to allow for automatical plotting of vector
!  31.10.2003	ggu	subroutine plowind()
!  17.03.2004	ggu	use okvar to decide which var to plot (0->all)
!  05.10.2004	ggu	use znv instead xv
!  05.10.2004	ggu	inorm instead ivert, overlay over z/h, comp_scale()
!  16.12.2004	ggu	plot also regular net
!  04.03.2005	ggu	gray plot over bathy
!  14.03.2007	ggu	wave plotting with plowave
!  17.09.2008	ggu	plot last layer
!  06.12.2008	ggu	new routine get_minmax_flag(), bvel -> ivel for velsh
!  09.01.2009	ggu	deleted traref, tramin
!  27.01.2009	ggu	bug in plonos (input variable ivar has been changed)
!  07.05.2009	ggu	bug fix for grid and color plotting (plobas)
!  14.09.2009	ggu	section plot for scalars in plonos
!  09.10.2009	ggu	new routine plopres() for atmos. pressure
!  13.10.2009	ggu	section plot for velocities in plosim
!  23.03.2010	ggu	changed v6.1.1
!  26.03.2010	ggu	section plot for velocities in plosim adapted
!  17.12.2010	ggu	substituted hv with hkv
!  31.03.2011	ggu	no plotting in dry nodes implemented - read fvl file
!  14.04.2011	ggu	changed VERS_6_1_22
!  17.05.2011	ggu	in plobas may plot node and element numbers
!  31.05.2011	ggu	changed VERS_6_1_23
!  12.07.2011	ggu	eliminated all references to out routines
!  31.08.2011	ggu	new eos plotting (pleos,ploeval)
!  01.09.2011	ggu	changed VERS_6_1_32
!  07.10.2011	ggu&dbf	error calling extelev with nkn, and not nel
!  18.10.2011	ggu	changed VERS_6_1_33
!  09.12.2011	ggu	changed VERS_6_1_38
!  24.01.2012	ggu	changed VERS_6_1_41
!  10.02.2012	ggu	belem in plobas to plot bathymetry on elements
!  26.03.2012	ccf&ggu	call mkht only for bvel .or. btrans (plo2vel)
!  01.06.2012	ggu	changed VERS_6_1_53
!  03.05.2013	ggu	changed VERS_6_1_63
!  13.06.2013	ggu	new routine plofem()
!  05.09.2013	ggu	endtime() and nplot introduced
!  12.09.2013	ggu	changed VERS_6_1_67
!  28.01.2014	ggu	changed VERS_6_1_71
!  30.05.2014	ggu	flag no data points and do not plot
!  18.07.2014	ggu	changed VERS_7_0_1
!  20.10.2014	ggu	new time management
!  30.10.2014	ggu	changed VERS_7_0_4
!  26.11.2014	ggu	changed VERS_7_0_7
!  05.12.2014	ggu	changed VERS_7_0_8
!  23.12.2014	ggu	changed VERS_7_0_11
!  15.01.2015	ggu	changed VERS_7_1_1
!  19.01.2015	ggu	changed VERS_7_1_2
!  19.01.2015	ggu	changed VERS_7_1_3
!  05.06.2015	ggu	some plotting routines adjourned (flag)
!  10.07.2015	ggu	changed VERS_7_1_50
!  17.07.2015	ggu	changed VERS_7_1_52
!  17.07.2015	ggu	changed VERS_7_1_80
!  20.07.2015	ggu	changed VERS_7_1_81
!  14.09.2015	ggu	prepared for plotting velocities given in fem file
!  06.11.2015	ggu	set valref to 1 if vel field == 0
!  19.02.2016	ggu	changed VERS_7_5_2
!  28.04.2016	ggu	changed VERS_7_5_9
!  25.05.2016	ggu	changed VERS_7_5_10
!  10.06.2016	ggu	changed VERS_7_5_13
!  14.06.2016	ggu	changed VERS_7_5_14
!  27.06.2016	ggu	changed VERS_7_5_16
!  12.01.2017	ggu	changed VERS_7_5_21
!  20.01.2017	ggu	changed VERS_7_5_22
!  21.03.2017	ggu	new parameter valmax introduced
!  31.03.2017	ggu	changed VERS_7_5_24
!  13.04.2017	ggu	changed VERS_7_5_25
!  09.05.2017	ggu	changed VERS_7_5_26
!  16.05.2017	ggu	changed VERS_7_5_27
!  11.07.2017	ggu	changed VERS_7_5_30
!  14.11.2017	ggu	changed VERS_7_5_36
!  17.11.2017	ggu	changed VERS_7_5_37
!  17.11.2017	ggu	changed VERS_7_5_38
!  24.01.2018	ggu	changed VERS_7_5_41
!  13.02.2018	ggu	new routine for plotting boxes (plobox)
!  22.02.2018	ggu	changed VERS_7_5_42
!  19.04.2018	ggu	changed VERS_7_5_45
!  16.10.2018	ggu	changed VERS_7_5_50
!  25.10.2018	ggu	changed VERS_7_5_51
!  18.12.2018	ggu	changed VERS_7_5_52
!  13.03.2019	ggu	in plobas use parameters from STR file
!  21.05.2019	ggu	changed VERS_7_5_62
!  13.07.2021	ggu	better box plotting
!  21.10.2021	ggu	allow for vertical velocity overlay
!  01.02.2023	ggu	use real area in plobox()
!  08.03.2025	ggu	plot nodal code
! 
!  notes :
! 
!  customize belem in plobas() to plot bathymetry on elements
!  customize bnumber in plobas() to write node and element numbers
! 
! **********************************************************
! **********************************************************
! **********************************************************

	subroutine plofem(type,ivarin)

!  3D concentrations

	use mod_hydro_plot
	use levels, only : nlvdi,nlv,ilhkv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) type
	integer ivarin

        character*80 line
	integer nrec,ivel,nplot
	integer ivaria,ivar
	integer level,isect
	logical femnext,ptime_ok,ptime_end
	integer getlev,getisec
	double precision atime

	nrec = 0
	nplot = 0
	level = getlev()
        isect = getisec()

	call femopen(type)
	ivar = ivarin
	write(6,*) 'ivar: ',ivar
	call checkvar(ivar)

        call mkvarline(ivar,line)
	ivaria = ivar

	do while( femnext(atime,ivaria,nlvdi,nkn,p3) )
	  nrec = nrec + 1
	  write(6,*) 'new record: ',nrec,ivaria,ivar,atime
	  if( ivar .ne. ivaria ) goto 99
	  !call ptime_info
	  if( ptime_end() ) exit	!FIXME
	  if( ptime_ok() ) then
	    nplot = nplot + 1
            write(6,*) '==================================='
            write(6,*) 'creating new plot with plofem... ',nplot
            write(6,*) '==================================='
            if( isect .eq. 0 ) then
	      write(6,*) '..........horizontal plotting nodes '
	      if( ivar .eq. 21 .or. ivar .eq. 2 ) then	!wind or vel
		ivel = 3		!wind
		!if( ivar .eq. 2 ) ivel = 3
		if( ivar .eq. 2 ) ivel = 5	!plot 2d vel field
		call set_uvfem(nlvdi,nkn,level,ilhkv,p3)
	        call reset_dry_mask
	        call plovel(ivel)
	      else
	        call extnlev(level,nlvdi,nkn,p3,parray)
	        call prepare_dry_mask
	        call ploval(nkn,parray,line)
	      end if
            else
	      write(6,*) '..........vertical plotting nodes '
              call plot_sect(.false.,p3)
            end if
            write(6,*) '==================================='
            write(6,*) 'end of plot'
            write(6,*) '==================================='
	  end if
	end do

	call femclose

	write(6,*) 'Total number of plots: ',nplot

	return
   99	continue
	write(6,*) 'ivar,ivaria: ',ivar,ivaria
	stop 'error stop plofem: no such variable in file'
	end

! **********************************************************

	subroutine set_uvfem(nlvddi,nkn,level,ilhkv,p)

	use mod_hydro_plot

	implicit none

	integer nlvddi
	integer nkn
	integer level
	integer ilhkv(nkn)
	real p(nlvddi,nkn,2)

	integer k,lev

	lev = max(1,level)

	do k=1,nkn
	  if( level .le. ilhkv(k) ) then
	    uvnode(k) = p(lev,k,1)
	    vvnode(k) = p(lev,k,2)
	  else
	    uvnode(k) = 0.
	    vvnode(k) = 0.
	  end if
	end do

	end

! **********************************************************

	subroutine set_uv(nlvddi,nkn,p)

	use mod_hydro_plot

	implicit none

	integer nlvddi
	integer nkn
	real p(nlvddi,nkn,2)

	integer k

	do k=1,nkn
	  uvnode(k) = p(1,k,1)
	  vvnode(k) = p(1,k,2)
	end do

	end

! **********************************************************

	subroutine plonos(type,ivar_in)

!  3D concentrations

	use mod_hydro_plot
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) type	!default extension for file
	integer ivar_in		!desired variable id


        character*80 line
	integer nrec,it,nplot
	integer ivaria
	integer level,ivar,isect
	logical nosnext,okvar,ptime_ok,ptime_end
	integer getlev,getisec

	nrec = 0
	nplot = 0
	level = getlev()
        isect = getisec()
	ivar = ivar_in

	call nosopen(type)
	call fvlopen('.fvl')
	call checkvar(ivar)

        call mkvarline(ivar,line)

	do while( nosnext(it,ivaria,nlvdi,p3) )
	  nrec = nrec + 1
	  call fvlnext(it,nlvdi,fvlv)
	  write(6,*) 'new record: ',nrec,it,ivaria,ivar
	  !call ptime_info()
	  if( ptime_end() ) exit
	  if( ptime_ok() .and. okvar(ivaria) ) then
	    nplot = nplot + 1
            write(6,*) '==================================='
            write(6,*) 'creating new plot with plonos... ',nplot
            write(6,*) '==================================='
            if( isect .eq. 0 ) then
	      write(6,*) '..........horizontal plotting nodes '
	      call extnlev(level,nlvdi,nkn,p3,parray)
	      call prepare_dry_mask
	      call ploval(nkn,parray,line)
            else
	      write(6,*) '..........vertical plotting nodes '
              call plot_sect(.false.,p3)
            end if
            write(6,*) '==================================='
            write(6,*) 'end of plot'
            write(6,*) '==================================='
	  end if
	end do

	call nosclose

	write(6,*) 'Total number of plots: ',nplot

	end

! **********************************************************

	subroutine ploeos(type,ivar_in)

!  3D concentrations (element values)

	use mod_hydro_plot
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) type	!default extension for file
	integer ivar_in		!desired variable id


        character*80 line
	integer nrec,it,nplot
	integer ivaria
	integer level,ivar,isect
	logical eosnext,okvar,ptime_ok,ptime_end
	integer getlev,getisec

	nrec = 0
	nplot = 0
	level = getlev()
        isect = getisec()
	ivar = ivar_in

	call eosopen(type)
	call checkvar(ivar)

        call mkvarline(ivar,line)

	do while( eosnext(it,ivaria,nlvdi,p3) )
	  nrec = nrec + 1
	  write(6,*) nrec,it,ivaria,ivar
	  if( ptime_end() ) exit
	  if( ptime_ok() .and. okvar(ivaria) ) then
	    nplot = nplot + 1
            write(6,*) '==================================='
            write(6,*) 'creating new plot with ploeos... ',nplot
            write(6,*) '==================================='
            if( isect .eq. 0 ) then
	      write(6,*) '..........horizontal plotting elements '
	      !call extelev(level,nlvdi,nkn,p3,parray)
	      call extelev(level,nlvdi,nel,p3,parray)
	      call prepare_dry_mask
	      call ploeval(nel,parray,line)
            else
	      write(6,*) '..........no vertical plotting for elements'
	      nplot = nplot - 1
            end if
            write(6,*) '==================================='
            write(6,*) 'end of plot'
            write(6,*) '==================================='
	  end if
	end do

	call eosclose

	write(6,*) 'Total number of plots: ',nplot

	end

! **********************************************************

        subroutine mkvarline(ivar,line)

!  makes variable description in nos file

        implicit none

        integer ivar            !variable number
        character*(*) line      !description on return

        integer iaux,isub
        integer ialfa

	call ivar2string(ivar,line,isub)

	if( line .eq. ' ' ) then
          line = 'variable = '
          iaux = ialfa(float(ivar),line(12:),-1,-1)
	end if

        end

! **********************************************************

	subroutine plozet

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nrec,it,k,nplot
	logical ousnext,ptime_ok,ptime_end

	nrec = 0
	nplot = 0

	call ousopen

	do while( ousnext(it) )
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  if( ptime_end() ) exit
	  if( ptime_ok() ) then
	    nplot = nplot + 1
	    write(6,*) '..........plotting ',nplot
	    call prepare_dry_mask
	    call ploval(nkn,znv,'water levels')
	  end if
	end do

	call ousclose

	write(6,*) 'Total number of plots: ',nplot

	end

! **********************************************************

	subroutine plobar

!  plots barene

	use mod_hydro_plot
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nrec,it,k,nplot
	logical ousnext,ptime_ok,ptime_end

	nrec = 0
	nplot = 0

	call ousopen

	do while( ousnext(it) )
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  if( ptime_end() ) exit
	  if( ptime_ok() ) then
	    nplot = nplot + 1
	    write(6,*) '..........plotting ',nplot
	    call plozbar
	  end if
	end do

	call ousclose

	write(6,*) 'Total number of plots: ',nplot

	end

! **********************************************************

	subroutine plowave

	implicit none

	integer nrec,it,nplot
        integer ivel
	logical wavenext,ptime_ok,ptime_end

	nrec = 0
	nplot = 0
        ivel = 4	!wave

	call waveopen

	do while( wavenext(it) )
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  if( ptime_end() ) exit
	  if( ptime_ok() ) then
	    nplot = nplot + 1
	    write(6,*) '..........plotting ',nplot
	    call reset_dry_mask
	    call plovel(ivel)
	  end if
	end do

	call waveclose

	write(6,*) 'Total number of plots: ',nplot

	end

! **********************************************************

	subroutine plosim(bvel)

	use mod_hydro_plot

	implicit none

	logical bvel

	integer nrec,it,nplot
        integer ivel,isect


	logical velnext,ptime_ok,ptime_end
	integer getisec

        if( bvel ) then
          ivel = 1
        else
          ivel = 2
        end if
        isect = getisec()

	nrec = 0
	nplot = 0

	call velopen

	do while( velnext(it) )
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  if( ptime_end() ) exit
	  if( ptime_ok() ) then
	    nplot = nplot + 1
            write(6,*) '==================================='
            write(6,*) 'creating new plot with plosim... ',nplot
            write(6,*) '==================================='
            if( isect .eq. 0 ) then
	      write(6,*) '..........horizontal plotting '
	      call prepare_dry_mask
	      call plovel(ivel)
            else
	      write(6,*) '..........vertical plotting '
	      call prepare_vel(p3)
              call plot_sect(.true.,p3)
            end if
            write(6,*) '==================================='
            write(6,*) 'end of plot'
            write(6,*) '==================================='
	  end if
	end do

	call velclose

	write(6,*) 'Total number of plots: ',nplot

	end

! **********************************************************

	subroutine ploreg(nreg,preg,regpar,title,bintp)

!  plots scalar values from regular grid

	use basin
	use mod_hydro_plot

	implicit none

	integer nreg		!total number of regular values
	real preg(nreg)		!regular scalar values
	real regpar(7)		!description of regular data
        character*(*) title	!title for plot
	logical bintp		!interpolate on fem grid

	integer nflag
	integer nx,ny
	integer ierr
	real x0,y0,dx,dy
	real pmin,pmax,flag
	integer np
	real, allocatable :: pa(:)

	if( bintp ) then
	  np = nkn
	  allocate(pa(np))
	  call getreg(regpar,nx,ny,x0,y0,dx,dy,flag)
	  call setregextend(.false.)
	  call setregextend(.true.)
	  call set_flag(flag)
          call intp_reg_nodes(nx,ny,x0,y0,dx,dy,flag,preg &
     &                          ,pa,ierr)
	  if( ierr /= 0 .and. berrintp ) then
	    write(6,*) 'intp_reg_nodes ierr = ',ierr
	    !stop 'error stop ploreg: interpolation'
	  end if
	else
	  flag = regpar(7)
	  call set_flag(flag)
	  np = nreg
	  allocate(pa(np))
	  pa = preg
	end if

	call setreg_grid(regpar)

	call qstart
        call annotes(title)
	call bash(0)

	call get_minmax_flag(pa,np,pmin,pmax,flag)
	call count_flag(pa,np,flag,nflag)
	call colauto(pmin,pmax)
        if( bminmax ) write(6,*) 'min/max: ',np,nflag,pmin,pmax

	if( bintp ) then
          call qcomm('Plotting interpolated regular grid')
          call isoline(pa,np,0.,2)
	else
          call qcomm('Plotting regular grid')
          call isoreg(pa,np,regpar,0.,2)
	end if

        call colsh
	call bash(4)	! overlays grid
	call bash(2)

	call qend

	end

! **********************************************************

	subroutine plo_scal_val(n,pa,title)

	use basin

	implicit none

	integer n
	real pa(n)
        character*(*) title

	if( n == nkn ) then
	  call ploval(nkn,pa,title)
	else if( n == nel ) then
	  call ploeval(nel,pa,title)
	else
	  write(6,*) 'no nodal or elemental values: ',n,nkn,nel
	  stop 'error stop plo_scal_val: cannot plot'
	end if

	end

! **********************************************************

	subroutine ploval(nkn,pa,title)

!  plots node values

	use mod_hydro_plot

	implicit none

	integer nkn
	real pa(nkn)
        character*(*) title

	integer nflag
	real pmin,pmax,flag
	real getpar

        call get_flag(flag)

	call qstart
        call annotes(title)
	call bash(0)

	call get_minmax_flag(pa,nkn,pmin,pmax,flag)
	call apply_dry_mask(bkwater,pa,nkn,flag)	!flags dry nodes
	call count_flag(pa,nkn,flag,nflag)
	call colauto(pmin,pmax)
        if( bminmax ) write(6,*) 'min/max: ',nkn,nflag,pmin,pmax

        call qcomm('Plotting isolines')
        call isoline(pa,nkn,0.,2)
	call plot_dry_areas
        call colsh

	call bash(4)	! overlays grid

	call bash(2)
	call qend

	end

! **********************************************************

	subroutine ploeval(nel,pa,title)

!  plots element values

	use mod_hydro_plot

	implicit none

	integer nel
	real pa(nel)
        character*(*) title

	real pmin,pmax,flag
	real getpar

        call get_flag(flag)

	call qstart
        call annotes(title)
	call bash(0)

	call get_minmax_flag(pa,nel,pmin,pmax,flag)
        if( bminmax ) write(6,*) 'min/max: ',nel,pmin,pmax
	call apply_dry_mask(bwater,pa,nel,flag)		!flags dry nodes
	call colauto(pmin,pmax)

        call qcomm('Plotting element values')
        call isoline(pa,nel,0.,3)			!plot on elements
	call plot_dry_areas
        call colsh

	call bash(4)	! overlays grid

	call bash(2)
	call qend

	end

! **********************************************************

	subroutine ploarea(n,area,mode,title)

!  plots element type values

	use mod_hydro_plot

	implicit none

	integer n
	integer area(n)
	integer mode		!2: node  3: elem
        character*(*) title

	logical bmincol
	integer iamax,iamin
	real pmin,pmax,flag
	integer color(n)
	real pa(n)
	real getpar

	bmincol = .true.

	call qstart
        call annotes(title)
	call bash(0)

	color = area
	bmincol = ( nint(getpar('icolmin')) /= 0 )
	if( bmincol ) then
	  call transform_color_to_min(color)
	end if

	iamin = minval(color)
	iamax = maxval(color)
        if( bminmax ) write(6,*) 'min/max: ',n,iamin,iamax
	call colset_reg(iamin,iamax)

	pa = color

        call qcomm('Plotting area of '//trim(title))
        call isoline(pa,n,0.,mode)			!plot on elements
        !call colsh

	call bash(4)	! overlays grid

	call bash(2)
	call qend

	end

! **********************************************************

	subroutine plobox(title)

!  interprets element type as boxes and plots

	use basin
	use mod_hydro_plot

	implicit none

        character*(*) title

	integer iamax,iamin
	integer ia,ie,ii,k,k2,i2
	integer i,n,k1
	real x,y,a
	real pmin,pmax,flag
	real pa(nel)
	double precision :: area
	character*80 string
	integer, allocatable :: na(:)
	integer, allocatable :: connect(:,:,:)
	double precision, allocatable :: aa(:)
	double precision, allocatable :: xa(:)
	double precision, allocatable :: ya(:)

	real getpar,area_elem

	call set_plotleg(0)

	call qstart
        call annotes(title)
	call bash(0)

	iamin = minval(iarv)
	iamax = maxval(iarv)
        if( bminmax ) write(6,*) 'min/max: ',nel,iamin,iamax
	call colset_reg(iamin,iamax)

	allocate(connect(0:ngr,nkn,2))
	allocate(na(0:iamax))
	allocate(aa(0:iamax))
	allocate(xa(0:iamax))
	allocate(ya(0:iamax))
	connect = 0
	na = 0
	aa = 0.
	xa = 0.
	ya = 0.
	area = 0.
	do ie=1,nel
	  ia = iarv(ie)
	  area = area_elem(ie)
	  aa(ia) = aa(ia) + 3*area
	  do ii=1,3
	    k = nen3v(ii,ie)
	    x = xgv(k)
	    y = ygv(k)
	    na(ia) = na(ia) + 1
	    !xa(ia) = xa(ia) + x
	    !ya(ia) = ya(ia) + y
	    xa(ia) = xa(ia) + area*x
	    ya(ia) = ya(ia) + area*y
	    i2 = mod(ii,3) + 1
	    k2 = nen3v(i2,ie)
	    call insert_connect(ngr,nkn,k,k2,ia,connect)
	  end do
	end do
	do ia=0,iamax
	  if( na(ia) > 0 ) then
	    !xa(ia) = xa(ia) / na(ia)
	    !ya(ia) = ya(ia) / na(ia)
	    xa(ia) = xa(ia) / aa(ia)
	    ya(ia) = ya(ia) / aa(ia)
	  end if
	end do

	pa = iarv

        call qcomm('Plotting box areas of '//trim(title))

	do k=1,nkn
	  n = connect(0,k,1)
	  do i=1,n
	    k1 = connect(i,k,1)
	    ia = connect(i,k,2)
	    if( ia >= 0 ) then
	      call qline(xgv(k),ygv(k),xgv(k1),ygv(k1))
	    end if
	  end do
	end do

	call qtxts(12)
	call qtxtcc(0,0)

	write(6,*) trim(title)
	write(6,*) '      index     x-center         y-center' // &
     &			'             area     index'
	do ia=0,iamax
	  if( na(ia) > 0 ) then
	    x = xa(ia)
	    y = ya(ia)
	    a = aa(ia)
	    write(string,'(i10)') ia
	    string = adjustl(string)
	    write(6,*) ia,x,y,a,trim(string)
	    call qtext(x,y,trim(string))
	  end if
	end do
	
	call bash(4)	! overlays grid

	call bash(2)
	call qend

	call set_plotleg(1)

	end

! **********************************************************

	subroutine insert_connect(ngr,nkn,k1,k2,ia,connect)

	implicit none

	integer ngr,nkn,k,ia
	integer connect(0:ngr,nkn,2)

	integer i,n,k1,k2

	if( k1 == 0 .or. k2 == 0 ) then
	  write(6,*)k1,k2,ia
	  stop 'error stop insert_connect: k == 0'
	end if

	n = connect(0,k1,1)
	do i=1,n
	  if( connect(i,k1,1) == k2 ) then
	   if( connect(i,k1,2) == ia ) then
	    connect(i,k1,2) = -1
	   end if
	   exit
	  end if
	end do
	if( i > n ) then
	  n = n + 1
	  if( n > ngr ) goto 99
	  connect(0,k1,1) = n
	  connect(n,k1,1) = k2
	  connect(n,k1,2) = ia
	end if

	n = connect(0,k2,1)
	do i=1,n
	  if( connect(i,k2,1) == k1 ) then
	   if( connect(i,k2,2) == ia ) then
	    connect(i,k2,2) = -1
	   end if
	   exit
	  end if
	end do
	if( i > n ) then
	  n = n + 1
	  if( n > ngr ) goto 99
	  connect(0,k2,1) = n
	  connect(n,k2,1) = k1
	  connect(n,k2,2) = ia
	end if

	return
   99	continue
	write(6,*) k1,k2,ia,n,ngr
	write(6,*) connect(:,k1,1)
	write(6,*) connect(:,k2,1)
	stop 'error stop insert_connect: ngr'
	end

! **********************************************************

	subroutine plozbar

	use mod_hydro
	use basin

	implicit none

	logical bdry
	integer ie,ii,k

	call qstart
	call bash(0)

        call qcomm('Plotting barene')
	do ie=1,nel
	  bdry = .false.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( zenv(ii,ie) .ne. znv(k) ) bdry = .true.
	  end do
	  if( bdry ) then
	    call plotel(ie,0.1)
	  else
	    call plotel(ie,0.8)
	  end if
	end do

	call bash(2)
	call qend

	end

! **********************************************************

	subroutine plovel(ivel)

	implicit none

	integer ivel

	if( ivel == 3 ) then
	  write(6,*) 'plotting wind...'
	  call plo2vel(ivel,'2D ')
	else if( ivel == 4 ) then
	  write(6,*) 'plotting wave...'
	  call plo2vel(ivel,'2D ')
	else if( ivel == 5 ) then
	  write(6,*) 'plotting velocities...'
	  call plo2vel(ivel,'3D ')
	else
	  write(6,*) 'plotting 3d...'
	  call plo3vel(ivel)
	end if

	end

! **********************************************************

	subroutine plo2vel(ivel,title)

	implicit none

	integer ivel			!what to plot
        character*(*) title

	call plovect(ivel,title,.false.)

	end

! **********************************************************

	subroutine plovect(ivel,title,bregdata)

!  plots entity with vectors
! 
!  ivel = 1	velocities
!  ivel = 2	transports
!  ivel = 3	wind
!  ivel = 4	waves
!  ivel = 5	velocities already prepared
! 
!  for ivel == 1,2:	utrans,vtrans must be set
!  for other values:	uvnode,vvnode must be set
! 
!  bonelem,bistrans are set in mod_hydro_plot

	use mod_hydro_plot
	use mod_depth
	use mod_hydro
	use basin
	use plotutil

	implicit none

	integer ivel			!what to plot
        character*(*) title		!title for annotation
	logical bregdata		!is data on regular grid?

	real uvmod(nkn)

	logical boverl,bnode,bbound,bnorm
        logical bspecial,bhasbasin
	logical bdebug
	logical bregplot
	logical bvel,btrans,bwind,bwave,bvelok
        character*80 anoline
        character*80 spcvel
	integer ie,k
	integer i,j,nnn
	integer level
	integer ioverl,inorm
	real u,v,xm,ym
	real x0,y0,dx,dy,flag
	real ut,vt
	real href,scale
	real velref,valref
	real velmin,valmin
	real velmax,valmax
	real typls,typlsf
	real color
	real xmin,ymin,xmax,ymax
	real pmin,pmax
	real uvmax,uvmed
	real typsca
	real fscale
	real regpar(7)
	integer nx,ny

	character*13, save :: varname(5) =	(/ &
     &						 'velocity     ' &
     &						,'transport    ' &
     &						,'wind velocity' &
     &						,'waves        ' &
     &						,'velocity     ' &
     &						/)

	real getpar,getcol
	integer getlev

! ------------------------------------------------------------------
!  set up parameters
! ------------------------------------------------------------------

! 	-----------------------------------------------------------
! 	locally defined parameters (please change here)
! 	-----------------------------------------------------------

	bdebug = .true.
	bdebug = .false.
	bnode = .true.		!plot arrows on nodes
	bnode = .false.		!plot arrows on nodes
	bbound = .true. 	!plot arrows on boundary nodes
	bbound = .false. 	!plot arrows on boundary nodes
	fscale = 3.		!extra empirical factor for computing arrows

! 	-----------------------------------------------------------
! 	global parameters from STR file
! 	-----------------------------------------------------------

        bvel   = ivel .eq. 1		!want vel but have trans
        btrans = ivel .eq. 2		!want trans
        bwind  = ivel .eq. 3
        bwave  = ivel .eq. 4
        bvelok = ivel .eq. 5		!want vel and have vel

	if( ivel < 1 .or. ivel > 5 ) then
          stop 'error stop plo2vel: internal error (3)'
        end if

	velref = getpar('velref')
	velmin = getpar('velmin')
	velmax = getpar('velmax')

	typlsf = getpar('typlsf')               !additional factor for length
	href = getpar('href')
        call getfnm('spcvel',spcvel)

	ioverl = nint(getpar('ioverl')) 	!mode of overlay plotting
	inorm  = nint(getpar('inorm'))	 	!normalization of arrows

	bspecial = spcvel .ne. ' '
	if( bspecial ) ioverl = 0

	boverl = ioverl .ne. 0			!overlay color
	bnorm  =  inorm .eq. 1			!normalize vectors

        anoline = title//varname(ivel)

! 	-----------------------------------------------------------
! 	check some settings
! 	-----------------------------------------------------------

	if( bvel .and. .not. bistrans ) goto 99
	if( btrans .and. .not. bistrans ) goto 99
	if( bistrans .and. .not. bonelem ) goto 99
	
! 	-----------------------------------------------------------
! 	start plotting
! 	-----------------------------------------------------------

	call qstart
	call annotes(anoline)
	call bash(0)

! ------------------------------------------------------------------
!  see if regular grid -> set bregplot and nx,ny if needed
! ------------------------------------------------------------------

	call getgeoflag(flag)
	if( bregdata ) then		!data is on regular grid (global value)
	  bregplot = .true.
	else				!see if we have to plot regular
	  call prepare_regular(nx,ny,bregplot)
	end if
	!if( bregplot ) write(6,*) 'plotting on regular grid'
	!write(6,*) 'bregdata,bregplot: ',bregdata,bregplot

! ------------------------------------------------------------------
!  prepare for velocity or transport
! ------------------------------------------------------------------

	valref = velref
	valmin = velmin
	valmax = velmax

	if( bvel ) then
	  call mkht(hetv,href)
	  call por2vel(nel,utrans,vtrans,uvelem,vvelem,hetv)
        else if( btrans ) then
	  uvelem = utrans
	  vvelem = vtrans
        else if( bwind .or. bwave .or. bvelok ) then
	  !ok
	else
	  stop 'error stop plo2vel: internal error (77)'
	end if

! ------------------------------------------------------------------
!  compute values on nodes -> 0 for dry nodes
! ------------------------------------------------------------------

	if( bonelem ) then
	  call intp2node(uvelem,uvnode,bwater)
	  call intp2node(vvelem,vvnode,bwater)
	else if( .not. bregdata ) then
	  call intp2elem(uvnode,uvelem,bwater)
	  call intp2elem(vvnode,vvelem,bwater)
	else
	  uvelem = 0.
	  vvelem = 0.
        end if

	nnn = 0
        if( nnn > 0 ) then
          write(112,*) 'writing uvelem ',nel
          write(112,*) (uvelem(i),i=1,nel,nnn)
          write(112,*) 'writing vvelem ',nel
          write(112,*) (vvelem(i),i=1,nel,nnn)
          write(112,*) 'writing utrans ',nel
          write(112,*) (utrans(i),i=1,nel,nnn)
          write(112,*) 'writing vtrans ',nel
          write(112,*) (vtrans(i),i=1,nel,nnn)
          write(112,*) 'writing uvnode ',nkn
          write(112,*) (uvnode(i),i=1,nkn,nnn)
          write(112,*) 'writing vvnode ',nkn
          write(112,*) (vvnode(i),i=1,nkn,nnn)
        end if

! ------------------------------------------------------------------
!  regular interpolation
! ------------------------------------------------------------------

	if( bregdata ) then	!is already regular grid - we just copy
	  call mod_hydro_get_regpar(regpar)	!regular data description
	  nx = nint(regpar(1))
	  ny = nint(regpar(2))
	  x0 = regpar(3)
	  y0 = regpar(4)
	  dx = regpar(5)
	  dy = regpar(6)
	  flag = regpar(7)
	  call setgeo(x0,y0,dx,dy,flag)
	  call mod_hydro_plot_regular_init(nx,ny)
	  ureg = reshape(uvnode(:),(/nx,ny/))
	  vreg = reshape(vvnode(:),(/nx,ny/))
	  bhasbasin = basin_has_read_basin()
	  if( bhasbasin ) then		!FIXME - might still be wrong
	    call am2av(ureg,uvnode,nx,ny)
	    call am2av(vreg,vvnode,nx,ny)
	  else
	    uvnode = 0.
	    vvnode = 0.
	  end if
	else if( bregplot ) then
	  call av2amk(bwater,uvnode,ureg,nx,ny)
	  call av2amk(bwater,vvnode,vreg,nx,ny)
	end if

! ------------------------------------------------------------------
!  compute modulus and maximum of velocity
! ------------------------------------------------------------------

	call moduv(uvnode,vvnode,uvmod,nkn,uvmax,uvmed) !compute mod/max/aver

	if( .not. bbound ) then		!set boundary vectors to 0
	  call bnd2val(uvnode,0.)	
	  call bnd2val(vvnode,0.)	
	end if

! ------------------------------------------------------------------
!  underlying color 
! ------------------------------------------------------------------

	if( boverl ) then
	  if( ioverl .eq. 1 ) then		!horizontal velocities
	    uvover = uvmod
	  else if( ioverl .eq. 2 ) then		!vertical velocities
	    uvover = wsnv
          else if( ioverl .eq. 3 ) then		!water levels
	    uvover = znv
          else if( ioverl .eq. 4 ) then		!bathymetry
	    uvover = hkv
	  else
            write(6,*) 'ioverl = ',ioverl
            stop 'error stop plo2vel: value not allowed for ioverl'
	  end if
	  call get_minmax_flag(uvover,nkn,pmin,pmax,flag)
	  call apply_dry_mask(bkwater,uvover,nkn,flag)
	  call colauto(pmin,pmax)
	  if( bverb ) write(6,*) 'overlay color... ',pmin,pmax
          call qcomm('Plotting overlay')
          call isoline(uvover,nkn,0.,2)
	end if

	call plot_dry_areas

! ------------------------------------------------------------------
!  get scale (if not given, typls is computed in bastlscale)
! ------------------------------------------------------------------

	typls  = getpar('typls')		!is available only now...
	if( valref .le. 0. ) then		!use maximum
	    if( bnorm ) then			!normalize vectors
	      valref = uvmax
	    else
	      valref = uvmax / fscale		!use empirical smaller value
	    end if
	end if
	typsca = typls * typlsf
	if( valref <= 0. ) valref = 1		!if vel field == 0

	if( bminmax ) then
	  if( .not. boverl ) then	!FIXME
	    call get_minmax_flag(uvmod,nkn,pmin,pmax,flag)
	  end if
	  write(6,*) 'max/med: ',uvmax,uvmed
	  write(6,*) 'min/max: ',pmin,pmax
	end if

	if( bverb ) then
	  write(6,*) 'scale (type): ',typls,typlsf,ioverl
	  write(6,1100) ' scale (value): ',valref,valmin,uvmax,uvmed
 1100	 format(a,4f14.3)
	end if

! ------------------------------------------------------------------
!  plotting of arrow
! ------------------------------------------------------------------

	call qgray(0.)

	if( bspecial ) then
	  scale = typsca / valref
          call aspecial(spcvel,xgv,ygv,uvnode,vvnode,scale)
        else if( bregplot ) then
	  call getgeo(x0,y0,dx,dy,flag)
	  do j=1,ny
	    do i=1,nx
	      u = ureg(i,j)
	      v = vreg(i,j)
	      if( u > flag .and. v > flag ) then
		xm = x0 + (i-1) * dx
		ym = y0 + (j-1) * dy
		call comp_scale(inorm,typsca,valmin,valmax,valref,u,v,scale)
	        call pfeil(xm,ym,u,v,scale)
	      end if
	    end do
	  end do
	else if( bnode ) then
	  do k=1,nkn
	    u = uvnode(k)
	    v = vvnode(k)
	    if( u > flag .and. v > flag ) then
	      call comp_scale(inorm,typsca,valmin,valmax,valref,u,v,scale)
	      call pfeil(xgv(k),ygv(k),u,v,scale)
	    end if
	  end do
	else
	  do ie=1,nel
	    call baric(ie,xm,ym)
	    u = uvelem(ie)
	    v = vvelem(ie)
	    if( u > flag .and. v > flag ) then
	      call comp_scale(inorm,typsca,valmin,valmax,valref,u,v,scale)
	      call pfeil(xm,ym,u,v,scale)
	    end if
	  end do
	end if

! ------------------------------------------------------------------
!  plotting of color bar or scale arrow
! ------------------------------------------------------------------

	if( boverl ) then			!plot color scale
	  call colsh
	end if
	!if( .not. bnorm ) then	                !vel arrows are still in scale
	if( inorm == 0 ) then	                !vel arrows are still in scale
	  scale = typsca / valref		!arrows are in scale
	  call velsh(ivel,scale,valref)
	end if

! ------------------------------------------------------------------
!  end of plot
! ------------------------------------------------------------------

	call bash(4)
	call bash(2)
	call qend

! ------------------------------------------------------------------
!  end of routine
! ------------------------------------------------------------------

	return
   99	continue
	write(6,*) 'problems... ',bvel,btrans,bistrans,bonelem
	stop 'error stop plo2vel: internal error (1)'
	end

! **********************************************************

	subroutine plo3vel(ivel)

	use mod_hydro_plot
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ivel

	integer ie,ii,k,l
	integer level,lmax,lact
        real utot,vtot
	logical bdebug

	integer getlev

	bdebug = .true.
	bdebug = .false.

!  make vertical velocities

	call make_vertical_velocity

!  handle level

	level = getlev()
	write(6,*) 'plo3vel: level = ',level,' ivel = ',ivel
	write(6,*) 'nlvdi: ',nlvdi,'  nlv: ',nlv

	if( ivel .eq. 3 .or. ivel .eq. 4 ) then	!wave or wind
	  call plo2vel(ivel,'3D ')
	  return
	end if

	if( level .lt. -1 ) stop 'internal error (1)'

	if( level .eq. 0 ) then		!barotropic
          do ie=1,nel
            utot = 0.
            vtot = 0.
            lmax = ilhv(ie)
            do l=1,lmax
              utot = utot + utlnv(l,ie)
              vtot = vtot + vtlnv(l,ie)
            end do
            utrans(ie) = utot
            vtrans(ie) = vtot
          end do
	  do k=1,nkn
	    wsnv(k) = wlnv(level,k)
	  end do
	else
	  do ie=1,nel
	    lmax = ilhv(ie)
	    lact = level
	    if( lact .eq. -1 ) lact = lmax

	    if( level .gt. lmax ) then	!no such layer
	      utrans(ie) = 0.
	      vtrans(ie) = 0.
	    else
	      utrans(ie) = utlnv(lact,ie)
	      vtrans(ie) = vtlnv(lact,ie)
	    end if
	  end do
	  do k=1,nkn
	    lmax = ilhkv(k)
	    lact = level
	    if( lact .eq. -1 ) lact = lmax
	    wsnv(k) = wlnv(lact,k)
	  end do
	end if

!  debug

	if( bdebug ) then
	  ie = 100
	  write(6,*) 'debugging 3d routine...'
	  write(6,*) ie,level,ilhv(ie),utrans(ie),vtrans(ie)
	  do l=1,ilhv(ie)
	    write(6,*) l,utlnv(l,ie),vtlnv(l,ie)
	  end do
	end if

!  call 2d routine

	call plo2vel(ivel,'3D ')
	
	end

! **********************************************************

	subroutine por2vel(n,up,vp,uv,vv,ht)

	implicit none

	integer n
	real up(n),vp(n),uv(n),vv(n),ht(n)

	where( ht > 0. )
	  uv = up / ht
	  vv = vp / ht
	else where
	  uv = 0.
	  vv = 0.
	end where

	end

! **********************************************************

	subroutine plobas

	use mod_depth
	use basin
	use mod_hydro_plot

	implicit none

	logical bnumber,belem,bbox
        real pmin,pmax

	real getpar

	!bnumber = .true.	! plot node and element numbers
	!bnumber = .false.
	!belem = .true.		! plot bathymetry on elements
	!belem = .false.
	!bbox = .true.		! plot boxes
	!bbox = .false.		!

	belem = ( nint(getpar('isoinp')) == 2 )
	bnumber = ( nint(getpar('inumber')) /= 0 )
	bbox = ( nint(getpar('ibox')) /= 0 )

	call reset_dry_mask
	call adjust_no_plot_area
	call make_dry_node_mask(bwater,bkwater)

!  only boundary line (plot 1)

	call qstart
	call bash(0)
	call bash(2)
	call qend

!  grid with black (plot 2)

	call qstart
	call bash(0)
	call bash(1)
	call bash(2)
	call qend

!  grid with gray (plot 3)

	call qstart
	call bash(0)
	call bash(3)
	call bash(2)
	call qend

!  bathymetry (gray or color) (plot 4)

	if( belem ) then
	  call ploeval(nkn,hev,'basin')
	else
	  call ploval(nkn,hkv,'basin')
	end if

!  bathymetry with overlay grid (black) (plot 5)

	call qstart
	call bash(0)

	call mima(hkv,nkn,pmin,pmax)
	call colauto(pmin,pmax)
        call isoline(hkv,nkn,0.,2)

        call qcomm('before basin overlay');
        call qlwidth(0.001)
	call bash(1)    !black grid
        call colsh
        call qcomm('after basin overlay');
        call qlwidth(0.01)
	call bash(2)
	call qend

!  bathymetry with overlay grid (gray) (plot 6)

	call qstart
	call bash(0)

	call mima(hkv,nkn,pmin,pmax)
	call colauto(pmin,pmax)
        call isoline(hkv,nkn,0.,2)

        call qcomm('before basin overlay');
        call qlwidth(0.001)
	call bash(3)    !gray grid
        call colsh
        call qcomm('after basin overlay');
        call qlwidth(0.01)
	call bash(2)
	call qend

!  element code (plot 7)

	if( maxval(iarv) > 0 ) then
	  write(6,*) 'plotting element code... ',maxval(iarv)
	  call ploarea(nel,iarv,3,'element code')
	  if( bbox ) call plobox('box information')
	end if

!  nodal code (plot 8)

	if( maxval(iarnv) > 0 ) then
	  write(6,*) 'plotting nodal code... ',maxval(iarnv)
	  call ploarea(nkn,iarnv,2,'nodal code')
	end if

	!call plot_partition

!  boundary line with regular overview grid (plot 9)

	call qstart
	call bash(0)
	call bash(2)
        call reggrid(10,0.,0.5)
	call qend

!  here only debug (node and element numbers)

	if( bnumber ) then

	call qstart
	call bash(0)
	call bash(3)
	call basin_number(1)	!external node
	call bash(2)
	call qend

	call qstart
	call bash(0)
	call bash(3)
	call basin_number(2)	!external element
	call bash(2)
	call qend

	call qstart
	call bash(0)
	call bash(3)
	call basin_number(-1)	!internal node
	call bash(2)
	call qend

	call qstart
	call bash(0)
	call bash(3)
	call basin_number(-2)	!internal element
	call bash(2)
	call qend

	end if

!  end of routine

	end

! **************************************************************

	subroutine plot_partition

	use basin

	implicit none

	integer nnp,nep
	integer ie,ii,k,col
	logical bunique
	integer area_node(nkn),area_elem(nel)
	integer values(nel)

	call basin_get_partition(nkn,nel,nnp,nep,area_node,area_elem)

	if( nep > 0 ) then
	  write(6,*) 'plotting element partition...'
	  call ploarea(nel,area_elem,3,'element partition')
	end if

	if( nnp > 0 ) then
	  values = 0
	  do ie=1,nel
	    bunique = .true.
	    k = nen3v(1,ie)
	    col = area_node(k)
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( col /= area_node(k) ) bunique = .false.
	    end do
	    if( bunique ) values(ie) = 1
	  end do
	  write(6,*) 'plotting node partition...'
	  call ploarea(nel,values,3,'node partition')
	end if

	end

! **************************************************************

	subroutine elkdep(nkn,hkv)

	implicit none

	integer nkn
	real hkv(nkn)

	call qstart

	call bash(0)

        call qcomm('Plotting isolines')
	call isoline(hkv,nkn,0.,2)
	call colsh

        call bash(2)

	call qend

	end

! **************************************************************

	subroutine eltype0(nel,iarv)

	implicit none

	integer nel
	integer iarv(nel)

	integer ie,ia

	call qstart

	call bash(0)

        call qcomm('Plotting element types')

	do ie=1,nel
	  ia = iarv(ie)
! 	  if( ia .eq. 1 ) then					!channels
! 	  if( ia .ge. 3 .and. ia .le. 5 .or. ia .eq. 9 ) then	!inlets
! 	  if( ia .ne. 0 ) then					!special type
! 	  if( ia .ge. 3 .and. ia .le. 5 .or. ia .eq. 9 ) then	!inlets
          if( ia .ge. 12 ) then                                 !valli pesca
            call plotel(ie,0.4)
          end if
	end do

        call bash(2)

	call qend

	end

! **************************************************************

	subroutine eldepth(nel,hm3v,auxv,color,title)

	implicit none

	integer nel
	real hm3v(3,nel)
	real auxv(nel)
	real color(nel)
	character*(*) title

	integer ii,ie
	real hm

	do ie=1,nel
	  hm = 0.
	  do ii=1,3
	    hm = hm + hm3v(ii,ie)
	  end do
	  auxv(ie) = hm / 3.
	end do

	do ie=1,nel
	  hm = auxv(ie)
	  color(ie) = 1.
	  if( hm .le. 0. ) color(ie) = 0.7
	end do
	call basplc(nel,color,'depth = 0.',title)

	do ie=1,nel
	  hm = auxv(ie)
	  color(ie) = 1.
	  if( hm .le. 0.2 ) color(ie) = 0.7
	end do
	call basplc(nel,color,'depth = 0.2',title)

	do ie=1,nel
	  hm = auxv(ie)
	  color(ie) = 1.
	  if( hm .le. 0.4) color(ie) = 0.7
	end do
	call basplc(nel,color,'depth = 0.4',title)

	do ie=1,nel
	  hm = auxv(ie)
	  color(ie) = 1.
	  if( hm .le. 0.6 ) color(ie) = 0.7
	end do
	call basplc(nel,color,'depth = 0.6',title)

	end

! **************************************************************

	subroutine eltype(nel,iarv,color,title)

	implicit none

	integer nel,iarv(nel)
	real color(nel)
	character*(*) title

	integer ia,ie

	do ie=1,nel
	  ia = iarv(ie)
	  if( ia .eq. 0 ) ia = 1	!treat barene and bassifondi the same
	  color(ie) = 0.1 * ( 10 - ia )
	end do
	call basplc(nel,color,'type of elements (1=0)',title)

	do ie=1,nel
	  ia = iarv(ie)
	  if( ia .eq. 0 ) ia = 1	!treat barene and bassifondi the same
	  if( ia .ne. 1 ) ia = 5
	  color(ie) = 0.1 * ( 10 - ia )
	end do
	call basplc(nel,color,'channels and lagoon',title)

	do ie=1,nel
	  ia = iarv(ie)
	  if( ia .gt. 9 ) then
		ia = 1
	  else
		ia = 5
	  end if
	  color(ie) = 0.1 * ( 10 - ia )
	end do
	call basplc(nel,color,'strange types and islands',title)

	do ie=1,nel
	  ia = iarv(ie)
	  if( ia .eq. 1 ) then
		ia = 5
	  else
		ia = 1
	  end if
	  color(ie) = 0.1 * ( 10 - ia )
	end do
	call basplc(nel,color,'barene',title)

	end

! **************************************************************

	subroutine plotel( ie , color )

	use basin

	implicit none

	integer ie
	real color

	integer ii,k
	real xp(3),yp(3)

	if( color .lt. 0. ) return

	do ii=1,3
	  k = nen3v(ii,ie)
	  xp(ii) = xgv(k)
	  yp(ii) = ygv(k)
	end do

	call qsetc(color)
	call qafill(3,xp,yp)

	end

! **************************************************************

	subroutine basplc( nel , color , line , descrr )

	implicit none

	integer nel
	real color(nel)
	character*(*) line,descrr

	integer ie

	call qstart

        call plot_basin(0)
	call frame(0)

        call qcomm('Plotting elements')
	do ie=1,nel
	  call plotel(ie,color(ie))
	end do

        call qcomm('Plotting basin')
        call plot_basin(2)

	call qend

	end

! *****************************************************************

	subroutine bnd2val(a,val)	

!  sets values on boundary to val

	use mod_geom
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real a(nkn)
	real val

	integer k

	logical isbound
	isbound(k) = kantv(1,k) .gt. 0

	do k=1,nkn
	  if( isbound(k) ) then
	    a(k) = val
	  end if
	end do

	end

! *****************************************************************

	subroutine comp_scale(mode,amax,vmin,vmax,vref,u,v,scale)

!  computes scale to apply for plot
! 
!  mode:
!        0       uses real length
!        1       normalizes vector (all vectors have same length)
!        2       uses linear scale
!        3       uses logarithmic scale

	implicit none

	integer mode		!type of scaling
	real amax		!velocity independent scale
	real vmin		!minimum velocity
	real vmax		!maximum velocity
	real vref		!reference velocity -> gives amax length
	real u,v		!velocity values
	real scale		!computed scale (return)

	real vmod,aux

	vmod = sqrt(u*u+v*v)

	if( vmod .le. vmin ) then
	  scale = 0.
	  return
	else if( vmax > 0. .and. vmod >= vmax ) then
	  scale = 0.
	  return
	end if

	if( mode .eq. 0 ) then
	  scale = amax / vref
	else if( mode .eq. 1 ) then
	  scale = amax / vmod
	else if( mode .eq. 2 ) then
	  aux = (vmod-vmin) / (vref-vmin)                       !aux in [0-1]
	  scale = ( amax / vmod ) * aux
	else if( mode .eq. 3 ) then
	  aux = log10( 1 + 9 * (vmod-vmin)/(vref-vmin) )        !aux in [0-1]
	  scale = ( amax / vmod ) * aux
	else
	  write(6,*) 'mode = ',mode
	  stop 'error stop comp_scale: no such mode'
	end if

	end

! *****************************************************************

	subroutine moduv(u,v,uv,n,uvmax,uvmed)

!  computes modulus for velocity vector and maximum and average

	implicit none

	integer n
	real u(n), v(n), uv(n)
	real uvmax
	real uvmed

	integer i,nn
	real um,vm,rmod,rmax,rmed
	real flag

        call get_flag(flag)
	rmax = 0.
	rmed = 0.
	nn = 0

	do i=1,n
	  um = u(i)
	  vm = v(i)
	  if( um > flag .and. vm > flag ) then
	    nn = nn + 1
	    rmod = sqrt( um*um + vm*vm )
	    rmax = max(rmax,rmod)
	    rmed = rmed + rmod
	    uv(i) = rmod
	  else
	    uv(i) = flag
	  end if
	end do

	uvmax = rmax
        uvmed = 0.
	if( nn > 0 ) uvmed = rmed / nn

	end

! *****************************************************************

	subroutine normuv(u,v,vv,n)

!  normalize horizontal velocities

	implicit none

	integer n
	real u(n), v(n), vv(n)

	where( vv > 0. )
	  u = u / vv
	  v = v / vv
	end where

	end

! *****************************************************************

	subroutine intp2elem(vnv,vev,bwater)

!  compute elemental values vev()

	use basin

	implicit none

        real vnv(nkn)
	real vev(nel)
	logical bwater(nel)

	integer ie,ii,k,iflag
	real sum,flag

        call get_flag(flag)

	do ie=1,nel
	  if( bwater(ie) ) then
            sum = 0.
	    iflag = 0
	    do ii=1,3
	      k = nen3v(ii,ie)
	      sum = sum + vnv(k)
	      if( vnv(k) > flag ) iflag = iflag + 1
	    end do
            vev(ie) = sum / 3.
	    if( iflag .ne. 3 ) vev(ie) = flag
          else
            vev(ie) = 0.
	  end if
	end do

	end

! *****************************************************************

	subroutine intp2node(vev,vnv,bwater)

!  compute nodal values vnv()

	use basin

	implicit none

	real vev(nel)
        real vnv(nkn)
	logical bwater(nel)

	integer ie,ii,k
	real r,flag
	real v2v(nkn)

        call get_flag(flag)

	do k=1,nkn
	  vnv(k) = 0.
	  v2v(k) = 0.
	end do

	do ie=1,nel
	 if( bwater(ie) ) then
	  do ii=1,3
	    if( vev(ie) > flag ) then
	      k = nen3v(ii,ie)
	      vnv(k) = vnv(k) + vev(ie)
	      v2v(k) = v2v(k) + 1.
	    end if
	  end do
	 end if
	end do

	do k=1,nkn
	  r = v2v(k)
	  if( r .gt. 0. ) then
	    vnv(k) = vnv(k) / v2v(k)
	  else
	    vnv(k) = 0.
	    vnv(k) = flag
	  end if
	end do

	end

! *****************************************************************

	subroutine count_flag(value,n,flag,nflag)

!  counts total number of flagged nodes

	implicit none

	integer n
	real flag
	real value(n)
	integer nflag

	integer i,iflag

	iflag = 0
	do i=1,n
	  if( value(i) .eq. flag ) iflag = iflag + 1
	end do

	!write(6,*) 'flagged values: (',flag,')',iflag,' of ',n

	nflag = iflag

	end

! *****************************************************************

	subroutine apply_dry_mask(bmask,value,n,flag)

!  aplies mask to nodal value

	implicit none

	integer n
	real flag
	real value(n)
	logical bmask(n)

	integer i

	do i=1,n
	  if( .not. bmask(i) ) value(i) = flag
	end do

	end

! *****************************************************************

        subroutine aspecial(file,xgv,ygv,uv,vv,scale)

!  plots arrow on special points

	use basin, only : nkn

        implicit none

        character*80 file
        real scale
        real xgv(nkn)
        real ygv(nkn)
        real uv(nkn)
        real vv(nkn)

        integer ndim
        parameter (ndim=100)

        integer ns
        integer ks(ndim)
        real us(ndim),vs(ndim)
        real xs(ndim),ys(ndim)
        save ns,ks,us,vs,xs,ys

        logical berror
        integer i,k
        real xn,yn,dist

        integer icall
        save icall
        data icall / 0 /

        if( icall .lt. 0 ) return
        if( icall .eq. 0 ) then
          write(6,*) 'reading special velocity file...'
          write(6,*) file(1:79)
          ns = ndim
          call rdspecial(file,ns,ks,us,vs)
          call n2int(ns,ks,berror)
          if( berror ) stop 'error stop aspecial: no such node'
          do i=1,ns
            xs(i) = xgv(ks(i))
            ys(i) = ygv(ks(i))
          end do
          icall = 1
        end if

        !call qhue(1.)
        call qgray(0.5)
        dist = 50.
	do i=1,ns
          k = ks(i)
          call traxy(xgv(k),ygv(k),us(i),vs(i),dist,xn,yn)
	  call pfeil(xn,yn,us(i),vs(i),scale)
	end do

        call qgray(0.)
	do i=1,ns
          k = ks(i)
	  call pfeil(xgv(k),ygv(k),uv(k),vv(k),scale)
	end do

        call qgray(0.)

        end

! *****************************************************************

        subroutine traxy(xo,yo,u,v,dist,xn,yn)

!  translates x/y to new coordinate to the left of direction (u,v) by dist

	implicit none

	real xo,yo,u,v,dist,xn,yn

	real sp,un,vn

          sp = sqrt(u**2+v**2)
          un = - v/sp
          vn = + u/sp
          xn = xo + dist * un
          yn = yo + dist * vn

          end

! *****************************************************************

        subroutine rdspecial(file,n,ks,us,vs)

        implicit none

        character*80 file
        integer n
        integer ks(n)
        real us(n)
        real vs(n)

        integer ndim,i,k
        real pi,rad
        real sp,dir
        real f(10)
        character*80 line

        integer iscanf

        ndim = n
        n = 0
        pi = 4. * atan(1.)
        rad = pi / 180.

        open(9,file=file,form='formatted',status='old')

    1   continue
          !read(9,*,end=2) i,k,sp,dir
          read(9,'(a)',end=2) line
          if( line(1:1) .eq. '#' ) goto 1       !comment
          i = iscanf(line,f,10)
          if( i .eq. 0 ) goto 1                 !empty
          if( i .lt. 4 ) goto 99
          n = n + 1
          if( n .gt. ndim ) stop 'error stop rdspecial: ndim'
          k = nint(f(2))
          if( i .eq. 4 ) then
            sp = f(3)
            dir = f(4)
          else if( i .ge. 6 ) then
            sp = f(5)
            dir = f(6)
          end if
          ks(n) = k
          write(6,*) i,k,sp,dir
          sp = sp / 100.
          dir = 90. - dir
          if( dir .lt. 0 ) dir = dir + 360.
          dir = dir * rad
          us(n) = sp * cos(dir)
          vs(n) = sp * sin(dir)
          goto 1
    2   continue
        close(9)

        return
   99   continue
        write(6,*) 'i = ',i
        write(6,*) line
        stop 'error stop rdspecial: line'
        end

! *****************************************************************

	subroutine get_minmax_flag(pa,nkn,pmin,pmax,flag)

        implicit none

        integer nkn
        real pa(nkn)
        real pmin,pmax
	real flag

        integer k
        real high,val

        high = 1.e+30

        pmin = high
        pmax = -high
        
        do k=1,nkn
          val = pa(k)
          if( val .ne. flag ) then
            pmin = min(pmin,val)
            pmax = max(pmax,val)
          end if
        end do

        end

! *****************************************************************

	subroutine plot_dry_areas

!  plots dry areas with gray shading

	use mod_hydro_plot
	use basin

	implicit none

	integer ie
	real dgray,hgray,h
	real x(3),y(3)

	double precision dgetpar

	dgray = dgetpar('dgray')
	hgray = dgetpar('hgray')
	if( dgray < 0 ) return

	call qgray(dgray)

	do ie=1,nel
	  h = sum(hm3v(:,ie))/3.
	  if( bwater(ie) .and. h > hgray ) cycle
	  if( .not. bplot(ie) ) cycle
	  call set_xy(ie,x,y)
	  call qafill(3,x,y)
	end do

	call qgray(0.)

	end

! *****************************************************************
! *****************************************************************
! *****************************************************************

	subroutine init_regular

	implicit none

	integer nx,ny
	real dxygrd,x,y
	real xmin,ymin,xmax,ymax
	real x0,y0,x1,y1

	logical inboxdim
	real getpar

        call setgeo(0.,0.,0.,0.,-999.)

        dxygrd = getpar('dxygrd')
	if( dxygrd <= 0 ) return

        call basmima(0)
        if( inboxdim(' ',x0,y0,x1,y1) ) then  !size of plotting is given
          call setbas(x0,y0,x1,y1)
        end if

	call getbas(xmin,ymin,xmax,ymax)

        x = xmin + dxygrd/2.
        y = ymin + dxygrd/2.
	nx = nint((xmax-xmin)/dxygrd)
	ny = nint((ymax-ymin)/dxygrd)

        call setgeo(x,y,dxygrd,dxygrd,-999.)

	end

! *****************************************************************

	subroutine prepare_regular(nx,ny,bregplot)

! checks if we have to plot results on regular grid

	use mod_hydro_plot

	implicit none

	integer nx,ny
	logical bregplot

	real x0,y0,dx,dy,flag
	real xmin,ymin,xmax,ymax

	nx = 0
	ny = 0
	call getgeo(x0,y0,dx,dy,flag)
	bregplot = dx .gt. 0. .and. dy .gt. 0.
	if( .not. bregplot ) return
	
	call getbas(xmin,ymin,xmax,ymax)
	nx = nint((xmax-xmin)/dx)
	ny = nint((ymax-ymin)/dy)

	call mod_hydro_plot_regular_init(nx,ny)	!nx,ny may be changed here

	!write(6,*) 'ggu: x0,y0,dx,dy ',x0,y0,dx,dy,nx,ny
	!write(6,*) 'ggu: ',xmin,ymin,xmax,ymax

	end

! *****************************************************************

	subroutine info_regular(bregplot,nx,ny,dx,dy)

	implicit none

	logical bregplot
	integer nx,ny
	real x0,y0,dx,dy,flag

	call prepare_regular(nx,ny,bregplot)
	if( .not. bregplot ) return

	call getgeo(x0,y0,dx,dy,flag)

	end

! *****************************************************************

	subroutine transform_color_to_min(ids)

	use basin
	use mod_color
	use mod_connect

	implicit none

	integer ids(nel)

	integer nmax,ncol,i,ierr
	integer, allocatable :: matrix(:,:)
	integer color(nel)
	integer ecv(3,nel)

	call connect_init(nkn,nel,nen3v,ierr)
	if( ierr /= 0 ) then
	  stop 'error stop transform_color_to_min: ierr/=0'
	end if
	call connect_get_ecv(nel,ecv)
	call make_elem_matrix(nel,ecv,ids,nmax,matrix)
	call color_graph(nel,ids,nmax,matrix,ncol,color)
	call release_color_matrix(matrix)
	call connect_release

	write(6,*) 'needed colors for coloring: ',ncol

	ids = color

	end

! *****************************************************************

