
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016-2020  Georg Umgiesser
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

! routines for shy post processing output
!
! revision log :
!
! 25.05.2016	ggu	changed VERS_7_5_10
! 30.05.2016	ggu	changed VERS_7_5_11
! 07.06.2016	ggu	changed VERS_7_5_12
! 14.06.2016	ggu	changed VERS_7_5_14
! 17.06.2016	ggu	changed VERS_7_5_15
! 27.06.2016	ggu	changed VERS_7_5_16
! 09.09.2016	ggu	changed VERS_7_5_17
! 23.09.2016	ggu	expand routine written
! 30.09.2016	ggu	changed VERS_7_5_18
! 05.10.2016	ggu	routines dealing with regular field copied to own file
! 11.10.2016	ggu	changed VERS_7_5_20
! 31.03.2017	ggu	changed VERS_7_5_24
! 09.05.2017	ggu	changed VERS_7_5_26
! 04.11.2017	ggu	changed VERS_7_5_34
! 14.11.2017	ggu	changed VERS_7_5_36
! 24.01.2018	ggu	changed VERS_7_5_41
! 24.05.2018	ccf	add outformat option off
! 06.07.2018	ggu	changed VERS_7_5_48
! 25.10.2018	ggu	changed VERS_7_5_51
! 10.02.2019	ggu	write final message for -averbas
! 16.02.2019	ggu	changed VERS_7_5_60
! 14.05.2019	ggu	for nc_output use interpolated values, allow ivar=1
! 16.05.2019	ggu	write regular grid information
! 28.01.2020	ggu	changes for vorticity
! 28.04.2023    ggu     update function calls for belem
! 03.10.2024    ggu     added var_dim
!
!***************************************************************
!
!===============================================================
        module shyelab_out
!===============================================================

	implicit none

	integer, save :: nwrite = 0
	integer, save :: nwtime = 0

	integer, save :: iformat = 1
	integer, save :: ntype = 1

	logical, save :: breg = .false.
	integer, save :: nxreg = 0
	integer, save :: nyreg = 0
	integer, save :: lmaxreg = 0
	real, save :: regpar(7) = 0.

	real, save, allocatable :: s2dvalue(:)
	real, save, allocatable :: zvalue(:)
	real, save, allocatable :: svalue(:,:)
	real, save, allocatable :: uvalue(:,:)
	real, save, allocatable :: vvalue(:,:)
	integer, save, allocatable :: ilcoord(:)
	real, save, allocatable :: xcoord(:)
	real, save, allocatable :: ycoord(:)
	real, save, allocatable :: hcoord(:)
	real, save, allocatable :: fmreg(:,:)
	real, save, allocatable :: fmextra(:,:)
	real, save, allocatable :: xlon(:)
	real, save, allocatable :: ylat(:)

	! arrays for nc routines

	integer, save :: iwrite_nc = 0
	!integer, save, allocatable :: var_ids(:)
	real, save, allocatable :: var3d(:)
	real, save, allocatable :: value2d(:,:)
	real, save, allocatable :: value3d(:,:,:)
	real, save, allocatable :: vnc3d(:,:,:)

!===============================================================
        end module shyelab_out
!===============================================================

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine shyelab_init_output(id,idout,ftype,nvar,ivars)

! initializes in case of btrans, bout, bsplit, bsumvar
!
! in case of bout also distinguishes between output formats

	use basin
	use levels
	use elabutil
	use shyelab_out
	use shyfile
        use mod_offline

	implicit none

	integer id,idout
	integer ftype
	integer nvar
	integer ivars(nvar)

	integer ierr,iunit,np,ncid
	integer nxy
	integer ftype_out
	logical bshy
	character*60 file,form
	character*80 string,title
	integer ifileo

	idout = 0

	if( .not. boutput ) return

	ftype_out = ftype
	if( bvorticity ) ftype_out = 2	!scalar

	bshy = ( outformat == 'shy' .or. outformat == 'native' )

	string = regstring
	call fem_regular_parse(string,regpar,nxreg,nyreg)
	if( nxreg > 0 .and. nyreg > 0 ) then
	  nxy = nxreg * nyreg
	  allocate(svalue(nlvdi,nxy),s2dvalue(nxy))
	  allocate(zvalue(nxy),uvalue(nlvdi,nxy),vvalue(nlvdi,nxy))
	  allocate(ilcoord(nxy),xcoord(nxy),ycoord(nxy),hcoord(nxy))
	  allocate(fmreg(4,nxy),fmextra(6,nkn))
	  allocate(xlon(nxreg),ylat(nyreg))
	  call fem_regular_setup(nxreg,nyreg,regpar,ilhv &
     &				,fmreg,fmextra &
     &				,ilcoord,xcoord,ycoord,hcoord &
     &				,xlon,ylat)
	  if( .not. bquiet ) call fem_reg_print(regpar)
	  breg = .true.
	  ntype = ntype + 10
	  !write(6,*) 'regular setup: ',regpar
	else
	  allocate(svalue(nlvdi,nkn),s2dvalue(nkn)) !only for nodal values
	  allocate(zvalue(nkn),uvalue(nlvdi,nkn),vvalue(nlvdi,nkn))
	  allocate(ilcoord(nkn),xcoord(nkn),ycoord(nkn),hcoord(nkn))
	  xcoord = xgv
	  ycoord = ygv
	  ilcoord = ilhkv
	  call makehkv_minmax(hcoord,+1)
	end if

	if( breg .and. ( bsplit .or. bshy ) ) then
	  write(6,*) 'regular output only for format gis, fem, nc'
	  stop 'error stop shyelab_init_output: regular output'
	end if

	!if( breg .and. outformat == 'gis' ) then
	!  write(6,*) 'regular output for gis not yet ready'
	!  stop 'error stop shyelab_init_output: gis not ready'
	!end if

	!if( breg .and. outformat == 'nc' ) then
	!  write(6,*) 'regular output for nc not yet ready'
	!  stop 'error stop shyelab_init_output: nc not ready'
	!end if

	if( bsplit ) then
	  ! nothing to initialize ... is done later
	else
	  if( outformat == 'shy' .or. outformat == 'native') then
	    file = 'out.shy'
            idout = shy_init(file)
	    if( idout <= 0 ) goto 74
            call shy_clone(id,idout)
            call shy_set_ftype(idout,ftype_out)
            call shy_set_nvar(idout,nvar)
            if( b2d ) call shy_convert_2d(idout)
	    if( bvorticity ) call shy_set_npr(idout,1)
            call shy_write_header(idout,ierr)
            if( ierr /= 0 ) goto 75
	  else if( outformat == 'gis' ) then
	    call gis_write_connect
	  else if( outformat == 'fem' ) then
	    file = 'out.fem'
	    call fem_file_write_open(file,iformat,iunit)
	    if( iunit <= 0 ) goto 74
	    idout = iunit
	  else if( outformat == 'nc' ) then
	    call shy_get_title(id,title)
	    call nc_set_quiet(bquiet)
	    call nc_output_set_vars(breg,nxreg,nyreg,hcoord,fmreg,xlon,ylat)
	    call nc_output_init(ncid,title,nvar,ivars,b2d,sncglobal)
	    idout = ncid
	  else if( outformat == 'off' ) then
	    file = 'out.off'
            iunit = ifileo(60,file,'unformatted','new')
            if( iunit <= 0 ) goto 74
	    idout = iunit
	  else
	    write(6,*) 'outformat = ',trim(outformat)
	    stop 'error stop: outformat not recognized'
	  end if
	end if

	return
   74	continue
        write(6,*) 'error opening file ',trim(file)
        stop 'error stop shyelab_init_output: opening file'
   75	continue
        write(6,*) 'error writing header, ierr = ',ierr
        write(6,*) 'file = ',trim(file)
        stop 'error stop shyelab_init_output: writing header'
	end

!***************************************************************

	subroutine shyelab_header_output(idout,ftype,dtime,nvar)

! writes header of record

	use basin
	use levels
	use elabutil
	use elabtime
	use shyelab_out
	use shyfile

	implicit none

	integer idout,ftype
	double precision dtime
	integer nvar,ncid
	logical bscalar,bhydro

	integer ierr,iunit,np,nvers,lmax

	if( .not. boutput ) return

        bhydro = ftype == 1
        bscalar = ftype == 2

	if( bsplit ) then
	  ! nothing to be done
	else
	  if( outformat == 'shy' .or. outformat == 'native') then
	    ! nothing to be done
	  else if( outformat == 'gis' ) then
	    ! nothing to be done
	  else if( outformat == 'fem' ) then
	    if( bhydro ) return
	    iunit = idout
	    nvers = 0
	    np = nkn
	    if( breg ) np = nxreg*nyreg
	    lmax = nlv
	    if( b2d ) lmax = 1
            call fem_file_write_header(iformat,iunit,dtime &
     &                          ,nvers,np,lmax &
     &                          ,nvar,ntype &
     &                          ,nlvdi,hlv,datetime_elab,regpar)
	  else if( outformat == 'nc' ) then
	    ncid = idout
	    call nc_output_time(ncid,dtime)
	  else if( outformat == 'off' ) then
	    if( bhydro ) return
            write(6,*) 'off format valid only on hydro files'
	    stop 'error stop: off outformat'
	  else
	    write(6,*) 'outformat = ',trim(outformat)
	    stop 'error stop: outformat not recognized'
	  end if
	end if

	end

!***************************************************************

	subroutine shyelab_record_output(id,idout,dtime,ivar,iv &
     &					,belem,n,m &
     &					,lmax,nlvddi,cv3)

! writes data record

	use basin
	use levels
	use elabutil
	use shyelab_out
	use shyfile

	implicit none

	integer id,idout
	double precision dtime
	logical belem
	integer ivar,iv,n,m,ncid,var_id,var_dim
	integer lmax,nlvddi
	real cv3(nlvddi,n*m)

	integer ierr,iunit,nvers
	integer np
	integer id_out
	integer ftype
	logical bscalar,bhydro,bshy
	character*60 string

	if( .not. boutput ) return
	if( bsumvar ) return

	call shy_get_ftype(id,ftype)
        bhydro = ftype == 1
        bscalar = ftype == 2
	bshy = ( outformat == 'shy' .or. outformat == 'native' )

	if( bhydro .and. .not. bshy ) return

	if( bverb ) write(6,*) 'writing to output: ',ivar

	if( breg ) then
	  if( n /= nkn ) goto 98
	  np = nxreg * nyreg
	  call fem_regular_interpolate_shell(regexpand,nlvdi,cv3,svalue)
	else if( .not. belem ) then
	  np = nkn
	  svalue(:,1:nkn) = cv3(:,1:nkn)
	else if( belem ) then
	  if( .not. bhydro ) then
	    !convert from nel to nkn
	    stop 'error stop shyelab_record_output: ' // &
     &				'belem==.true. not ready'
	  end if
	else if( .not. bsplit .and. .not. bshy ) then
	  goto 98
	end if

	if( .not. bsplit .and. .not. bshy ) then
	  !if( ivar == 1 .or. ivar == 3 ) goto 99
	  if( ivar == 3 ) goto 99
	  if( m /= 1 ) goto 98
	end if

!	the value to be written is in svalue - already interpolated

	if( bsplit .and. .not. bhydro ) then
	  call shy_split_id(ivar,id,id_out)
	  call shy_write_output_record(id_out,dtime,ivar &
     &						,belem,n,m &
     &						,lmax,nlvddi,cv3)
	else
	  if( bshy ) then
	    call shy_write_output_record(idout,dtime,ivar &
     &						,belem,n,m &
     &						,lmax,nlvddi,cv3)
	  else if( outformat == 'gis' ) then
            call gis_write_record(dtime,ivar,np,nlvddi,ilcoord &
     &					,svalue,xcoord,ycoord)
	  else if( outformat == 'fem' ) then
	    iunit = idout
	    nvers = 0
	    call get_string_description(ivar,string)
            call fem_file_write_data(iformat,iunit &
     &                          ,nvers,np,lmax &
     &                          ,string &
     &                          ,ilcoord,hcoord &
     &                          ,nlvddi,svalue)
	  else if( outformat == 'nc' ) then
	    ncid = idout
	    call nc_output_get_var_id(iv,var_id)
	    call nc_output_get_var_dim(iv,var_dim)
	    call nc_output_record(ncid,var_id,var_dim,np,svalue)
	  else if( outformat == 'off' ) then
	    ! nothing to be done
	  else
	    write(6,*) 'outformat = ',trim(outformat)
	    stop 'error stop: outformat not recognized'
	  end if
	end if

	nwrite = nwrite + 1
	if( iv == 1 ) nwtime = nwtime + 1 

	return
   98   continue
        write(6,*) 'output format = ',trim(outformat)
        write(6,*) 'n,m = ',n,m
        write(6,*) 'nkn,nel = ',nkn,nel
        stop 'error stop shyelab_record_output: cannot handle'
   99   continue
        write(6,*) 'output format = ',trim(outformat)
        write(6,*) 'ivar = ',ivar
        stop 'error stop shyelab_record_output: cannot handle'
	end

!***************************************************************

	subroutine shyelab_post_output(id,idout,dtime,nvar,n,m,nndim &
     &					,lmax,nlvddi,cv3all)

! writes complete time record after loop

	use basin
	use levels
	use elabutil
	use shyelab_out
	use shyfile

	implicit none

	integer id,idout
	double precision dtime
	integer nvar,n,m,nndim
	integer lmax,nlvddi
	real cv3all(nlvddi,nndim,0:nvar)

	integer ierr,iunit,nvers
	integer it,nb,ivar,np,ncid,irx
	integer id_out
	integer ftype
	logical bscalar,bhydro,bshy,belem
	character*60 string

        real, allocatable :: znv(:)
        real, allocatable :: uprv(:,:)
        real, allocatable :: vprv(:,:)
        real, allocatable :: sv(:,:)
        real, allocatable :: dv(:,:)
        real, allocatable :: cvaux(:,:)

	if( .not. boutput ) return

	call shy_get_ftype(id,ftype)
        bhydro = ftype == 1
        bscalar = ftype == 2
	bshy = ( outformat == 'shy' .or. outformat == 'native' )

	if( bscalar .and. .not. bsumvar ) return

	if( bhydro ) then
	  allocate(znv(nkn),uprv(nlvddi,nkn),vprv(nlvddi,nkn))
	  allocate(sv(nlvddi,nkn),dv(nlvddi,nkn),cvaux(nlvddi,nkn))
          call prepare_hydro(.true.,nndim,cv3all,znv,uprv,vprv)
          call convert_to_speed(uprv,vprv,sv,dv)
	  if( breg ) then
	    np = nxreg * nyreg
	    irx = regexpand
	    cvaux(1,:) = znv(:)		!we must have nlvddi vertical dimension
	    call fem_regular_interpolate_shell(irx,1,cvaux,svalue)
	    zvalue(:) = svalue(1,:)
	    call fem_regular_interpolate_shell(irx,nlvddi,uprv,uvalue)
	    call fem_regular_interpolate_shell(irx,nlvddi,vprv,vvalue)
	  else
	    np = nkn
	    zvalue = znv
	    uvalue = uprv
	    vvalue = vprv
	  end if
	end if	

	if( bsplit ) then
	  if( bhydro ) then
            call shy_split_hydro(id,dtime,znv,uprv,vprv,sv,dv)
	  end if
	else
	  if( bshy ) then
	    if( bhydro ) then
	      ! nothing to be done
	    else if( bsumvar ) then
              cv3all(:,:,0) = 0.
              cv3all(:,:,0) = sum(cv3all,dim=3)
              ivar = 10
	      belem = .false.
	      call shy_write_output_record(idout,dtime,ivar &
     &					,belem,n,m &
     &					,lmax,nlvddi,cv3all(:,:,0))
	    else
	      stop 'error stop shyelab_post_output: internal error (1)'
	    end if
	  else if( outformat == 'gis' ) then
	    call gis_write_hydro(dtime,np,nlvddi,ilcoord &
     &				,zvalue,uvalue,vvalue,xcoord,ycoord)
	  else if( outformat == 'fem' ) then	!also covers breg
	    call fem_write_hydro(idout,dtime,np,nlvddi &
     &					,zvalue,uvalue,vvalue)
	  else if( outformat == 'nc' ) then
	    ncid = idout
	    call nc_output_hydro(ncid,znv,uprv,vprv)
	  else if( outformat == 'off' ) then
            call off_output_hydro(idout,dtime,nndim,cv3all)
	  else
	    write(6,*) 'outformat = ',trim(outformat)
	    stop 'error stop: outformat not recognized'
	  end if
          if( .not. (bshy.and.bhydro) ) then
	    nwrite = nwrite + 1
	    nwtime = nwtime + 1
	  end if
	end if

	if( bhydro ) deallocate(znv,uprv,vprv,sv,dv,cvaux)

	end

!***************************************************************

	subroutine shyelab_final_output(id,idout,nvar,ivars)

! writes info on end of routine

	use basin
	use levels
	use elabutil
	use shyelab_out
	use shyfile
	use shyfem_strings

	implicit none

	integer id,idout
	integer nvar
	integer ivars(nvar)

	integer ierr,iunit,nvers
	integer it,nb,ivar,ncid,iv
	integer id_out
	integer ftype
	logical bscalar,bhydro,bshy
	character*80 file
	character*5 format
	character*10 short
	character*40 full

	if( .not. boutput .and. .not. baverbas ) return
	if( bsilent ) return

	write(6,*) 'output written to following files: '

	call shy_get_ftype(id,ftype)
        bhydro = ftype == 1
        bscalar = ftype == 2
	bshy = ( outformat == 'shy' .or. outformat == 'native' )

	if( bnodes ) then
	  write(6,*) '  what.dim.node'
          write(6,*) 'what is one of the following:'
	  do iv=1,nvar
	    ivar = ivars(iv)
            call strings_get_short_name(ivar,short)
            call strings_get_full_name(ivar,full)
            write(6,*) '  ',short,full
          end do
          write(6,*) 'dim is 2d or 3d'
          write(6,*) '  2d for depth averaged variables'
          write(6,*) '  3d for output at each layer'
          write(format,'(i5)') nnodes
          format = adjustl(format)
	  write(6,'(a)') ' node is consecutive node numbering: ' &
     &				//'1-'//format
	end if

	if( bsplit ) then

	  if( bhydro ) then
	    write(6,*) 'zeta.shy'
	    write(6,*) 'velx.shy'
	    write(6,*) 'vely.shy'
	    write(6,*) 'speed.shy'
	    write(6,*) 'dir.shy'
	  else if( bscalar ) then
	    ivar = 0
	    do
	      ivar = ivar - 1
	      call shy_split_id(ivar,id,idout)
	      if( idout < 0 ) exit		!no more variables
	      if( idout > 0 ) then
		call shy_get_filename(idout,file)
	        write(6,*) trim(file)
	      end if
	    end do
	  end if

	else if( boutput ) then

	  if( bshy ) then
	    write(6,*) 'out.shy'
	  else if( outformat == 'gis' ) then
	    write(6,*) 'connectivity.gis'
	    write(6,*) 'extract_*.gis'
	  else if( outformat == 'fem' ) then
	    write(6,*) 'out.fem'
	  else if( outformat == 'nc' ) then
	    ncid = idout
	    call nc_output_final(ncid)
	    write(6,*) 'out.nc'
	  else if( outformat == 'off' ) then
	    write(6,*) 'out.off'
	  else
	    write(6,*) 'outformat = ',trim(outformat)
	    stop 'error stop: outformat not recognized'
	  end if

	else if( baverbas ) then
	  write(6,*) '  what.dim.area'
          write(6,*) 'what is one of the following:'
	  do iv=1,nvar
	    ivar = ivars(iv)
            call strings_get_short_name(ivar,short)
            call strings_get_full_name(ivar,full)
            write(6,*) '  ',short,'       ',full
          end do
          write(6,*) '  ','volume_and_area','  volume and area'
          write(6,*) 'dim is 0d'
          write(6,*) 'area is 0'
	  if( barea ) then
	    write(6,*) 'the average is relative to the area given'
	    write(6,*) 'by the line(s) contained in ',trim(areafile)
	  else
	    write(6,*) 'the average is relative to the whole basin'
	  end if
	end if

	end

!***************************************************************
!***************************************************************
!***************************************************************


!***************************************************************
!***************************************************************
!***************************************************************

	subroutine fem_write_hydro(idout,dtime,np,nlvddi,zv,uv,vv)

! writes hydro records

	use levels
	use elabtime
	use shyelab_out

	implicit none

	integer idout
	double precision dtime
	integer np
	integer nndim
	integer nlvddi
	real zv(np)
	real uv(nlvddi,np)
	real vv(nlvddi,np)

	integer iunit,nvers,lmax,nvar,ivar
	character*60 string,stringx,stringy

	iunit = idout
	nvers = 0
	lmax = nlv
	nvar = 3

        call fem_file_write_header(iformat,iunit,dtime &
     &                          ,nvers,np,lmax &
     &                          ,nvar,ntype &
     &                          ,nlvdi,hlv,datetime_elab,regpar)

	ivar = 1
	lmax = 1
	call get_string_description(ivar,string)

        call fem_file_write_data(iformat,iunit &
     &                          ,nvers,np,lmax &
     &                          ,string &
     &                          ,ilcoord,hcoord &
     &                          ,lmax,zv)

	ivar = 2
	lmax = nlv
	call get_string_description(ivar,string)
	stringx = trim(string) // ' x'
	stringy = trim(string) // ' y'

        call fem_file_write_data(iformat,iunit &
     &                          ,nvers,np,lmax &
     &                          ,stringx &
     &                          ,ilcoord,hcoord &
     &                          ,nlvddi,uv)

        call fem_file_write_data(iformat,iunit &
     &                          ,nvers,np,lmax &
     &                          ,stringy &
     &                          ,ilcoord,hcoord &
     &                          ,nlvddi,vv)

	end

!***************************************************************

	subroutine get_string_description(ivar,string)

	implicit none

	integer ivar
	character*(*) string

	integer isub

	call ivar2string(ivar,string,isub)

	if( string == ' ' ) then
	  stop 'error stop get_string_description: unknown ivar'
	end if

	end

!***************************************************************
!***************************************************************
!***************************************************************

        subroutine fem_regular_interpolate_shell(regexpand,lmax,cv3,am)

        use basin
        use levels
        use shyelab_out

        implicit none

        integer regexpand
        integer lmax                   !vertical dimension of regular array
        real cv3(nlvdi,nkn)            !values of fem array
        real am(nlvdi,nxreg*nyreg)     !interpolated values (return)

	integer k,l,lm
	real flag

	call getgeoflag(flag)

	do k=1,nkn
	  lm = min(lmax,ilhkv(k))
	  !do l=1,lm
	  !  if( cv3(l,k) == flag ) then
	  !    write(6,*) 'warning: flag in scals ',l,k
	  !  end if
	  !end do
	  cv3(lm+1:nlvdi,k) = flag
	end do

	call fem_regular_interpolate(nxreg,nyreg,regexpand,lmax &
     &                  ,fmreg,fmextra,ilcoord,cv3,am)

	end

!***************************************************************

	subroutine shyelab_increase_nwrite

        use shyelab_out

	implicit none

	nwrite = nwrite + 1

	end

!***************************************************************

	subroutine shyelab_get_nwrite(nw,nt)

        use shyelab_out

	implicit none

	integer nw,nt

	nw = nwrite
	nt = nwtime

	end

!***************************************************************

