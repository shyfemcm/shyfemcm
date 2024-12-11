
!==========================================================================
        module mod_dbg_flx
!==========================================================================

        implicit none

        logical, save :: bsilent = .false.      !be silent
        logical, save :: bquiet = .false.       !be quiet
        logical, save :: bverbose = .false.     !be verbose
        logical, save :: bcheck = .true.        !check for differences
        logical, save :: bstop = .true.         !stop on error
        logical, save :: bnostop = .false.      !do not stop on differences
        logical, save :: bnodiff = .true.       !do not show differences
        logical, save :: bsummary = .true.      !only show summary
        logical, save :: bbalance = .true.      !balance time step
        integer, save :: maxdiff = 0.           !max difference allowed

!==========================================================================
        end module mod_dbg_flx
!==========================================================================

        program dbg_flx

        use clo
        use mod_dbg_flx

        implicit none

        integer nc,ierr

        call dbg_flx_init

        nc = clo_number_of_files()

        if( .not. bquiet ) write(6,*) 'running dbg_flx...'

        if( nc == 0 ) then
          call clo_usage
        else if( nc == 1 ) then
          call read_dbg_flx_file
        else if( nc == 2 ) then
          !call compare_files(ierr)
        else
          write(6,*) 'nc = ',nc
          stop 'error stop dbg_flx: wrong number of files'
        end if

        if( ierr > 0 ) then
          if( ierr == 99 ) ierr = 100   !terrible hack - FIXME
          call exit(ierr)
        else
          call exit(99)
        end if

        end

!**************************************************************************

        subroutine read_dbg_flx_file

! reads one file and outputs info

        use clo
        use mod_dbg_flx

        implicit none

        integer nc
        integer ntime,nrec
        integer nvers,nsect,kfluxm,idtflx
        integer nlmax,nvar,ivar,iv,ierr
	integer is
	integer nlvddi,lmax,l
	integer nh,nv,nt
	integer iu1
        integer ios
        double precision dtime,atime0,atime
        character*60 name_one,text
        character*80 title,femver

	integer, save, allocatable :: kflux(:)
	integer, save, allocatable :: nlayers(:)
	real, save, allocatable :: fluxes(:,:,:)
	character*80, save, allocatable :: strings(:)

        call clo_get_file(1,name_one)

	iu1 = 1
        open(iu1,file=name_one,status='old',form='unformatted',iostat=ios)

        if( ios /= 0 ) then
          write(6,*) 'no such file: ',trim(name_one)
          stop 'error opening file'
        end if

        if( .not. bquiet ) write(6,*) 'file: ',trim(name_one)

	call flx_is_flx_file(iu1,nvers)
	if( nvers == 0 ) then
          write(6,*) 'file is not a flx file: ',trim(name_one)
          stop 'error opening file'
	end if

        ntime = 0
        call flx_read_header(iu1,nvers,nsect,kfluxm,idtflx &
     &                                  ,nlmax,nvar,ierr)
	if( ierr /= 0 ) goto 99

	allocate(kflux(kfluxm))
	allocate(nlayers(nsect))
	allocate(strings(nsect))
        call flx_read_header2(iu1,nvers,nsect,kfluxm &
     &                          ,kflux,nlayers &
     &                          ,atime0,title,femver,strings,ierr)
	if( ierr /= 0 ) goto 98

	write(6,*) nsect,kfluxm,idtflx,nlmax,nvar
	if( bverbose ) then
	  do is=1,nsect
	    write(6,*) is,nlayers(is),trim(strings(is))
	  end do
	end if

	iv = 0
	nlvddi = nlmax
	allocate(fluxes(0:nlvddi,3,nsect))

        do while(.true.)

          call flx_read_record(iu1,nvers,atime &
     &                  ,nlvddi,nsect,ivar &
     &                  ,nlayers,fluxes,ierr)

	  if( ierr < 0 ) exit
	  if( ierr > 0 ) goto 97

	  iv = iv + 1
          if( mod(iv,nvar) == 1 ) ntime = ntime + 1
          if( .not. bquiet ) write(6,*) 'time = ',atime,ntime,ivar
	  if( bverbose ) then
	    do is=1,nsect
	      lmax = nlayers(is)
	      do l=1,lmax
	        write(6,*) l,fluxes(l,:,is)
	      end do
	    end do
	  end if

        end do

        if( .not. bsilent ) write(6,*) 'time records read: ',ntime

	return
   97	continue
	write(6,*) 'error reading record: ',ierr
	stop 'error stop read_dbg_flx_file: error reading record'
   98	continue
	write(6,*) 'error reading header2: ',ierr
	stop 'error stop read_dbg_flx_file: error reading header2'
   99	continue
	write(6,*) 'error reading header: ',ierr
	stop 'error stop read_dbg_flx_file: error reading header'
        end


!**************************************************************************

        subroutine dbg_flx_init

        use clo
        use mod_dbg_flx

        implicit none

        logical baux
        character*80 version

        version = '1.0'

        call clo_init('dbg_flx_init','flx-file(s)',trim(version))

        call clo_add_info('reads flx files')

        call clo_add_sep('general options:')
        call clo_add_option('silent',.false.,'be silent')
        call clo_add_option('quiet',.false.,'be quiet')
        call clo_add_option('nodiff',.false.,'do not show differences')
        call clo_add_option('verbose',.false.,'be verbose')
        call clo_add_option('nostop',.false.,'do not stop at error')
        call clo_add_option('summary',.false.,'do only summary')
        call clo_add_option('balance',.false.,'balance time records')
        call clo_add_option('maxdiff',0.,'maximum tolerated difference')

        call clo_parse_options

        call clo_get_option('silent',bsilent)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('nodiff',bnodiff)
        call clo_get_option('verbose',bverbose)
        call clo_get_option('nostop',baux)
        call clo_get_option('summary',bsummary)
        call clo_get_option('balance',bbalance)
        call clo_get_option('maxdiff',maxdiff)

        if( baux ) bstop = .false.
        if( bsilent ) bquiet = .true.
        if( bquiet ) bverbose = .false.

        end

!**************************************************************************


