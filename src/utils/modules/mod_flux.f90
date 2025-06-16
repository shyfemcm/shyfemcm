
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-2000,2002-2003,2006-2007,2006-2007  Georg Umgiesser
!    Copyright (C) 2009-2020  Georg Umgiesser
!    Copyright (C) 2017  Christian Ferrarin
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

! definition of flux module
!
! revision log :
!
! 30.04.1998	ggu	newly written routines (subpor deleted)
! 07.05.1998	ggu	check nrdveci on return for error
! 08.05.1998	ggu	restructured with new comodity routines
! 13.09.1999	ggu	type of node computed in own routine flxtype
! 19.11.1999	ggu	iskadj into sublin
! 20.01.2000	ggu	old routines substituted, new routine extrsect
! 20.01.2000	ggu	common block /dimdim/ eliminated
! 20.01.2000	ggu	common block /dimdim/ eliminated
! 01.09.2002	ggu	ggu99 -> bug in flx routines (how to reproduce?)
! 26.05.2003	ggu	in flxnov substituted a,b with b,c
! 26.05.2003	ggu	new routine make_fluxes (for lagrangian)
! 10.08.2003	ggu	do not call setweg, setnod, setkan
! 23.03.2006	ggu	changed time step to real
! 28.09.2007	ggu	use testbndo to determine boundary node in flxtype
! 28.04.2009	ggu	links re-structured
! 23.03.2010	ggu	changed v6.1.1
! 23.02.2011	ggu	new routine call write_node_fluxes() for special output
! 01.03.2011	ggu	changed VERS_6_1_20
! 01.06.2011	ggu	documentation to flxscs() changed
! 14.07.2011	ggu	changed VERS_6_1_27
! 21.09.2011	ggu	some lower-level subroutines copied to subflx.f
! 07.10.2011	ggu	adjusted for 3d flux routines
! 18.10.2011	ggu	changed VERS_6_1_33
! 19.10.2011	ggu	added T/S variables, created fluxes_*() routines
! 19.10.2011	ggu	added conz variables, created fluxes_template()
! 09.12.2011	ggu	changed VERS_6_1_38
! 01.06.2012	ggu	changed VERS_6_1_53
! 10.05.2013	ggu	introduced subflxa.h, common routines to subflxu.f
! 25.10.2013	ggu	changed VERS_6_1_68
! 07.03.2014	ggu	changed VERS_6_1_72
! 26.11.2014	ggu	changed VERS_7_0_7
! 19.12.2014	ggu	changed VERS_7_0_10
! 19.01.2015	ggu	changed VERS_7_1_3
! 20.05.2015	ggu	modules introduced
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 09.11.2015	ggu	changed VERS_7_3_13
! 12.04.2016	ggu	fluxes_template adjourned
! 15.04.2016	ggu	fluxes_template debugged and finished
! 09.09.2016	ggu	changed VERS_7_5_17
! 31.03.2017	ggu	changed VERS_7_5_24
! 22.09.2017	ccf	added total sediment concentration
! 26.10.2017	ggu	reads itable, chflx and section description
! 04.11.2017	ggu	changed VERS_7_5_34
! 14.11.2017	ggu	changed VERS_7_5_36
! 17.11.2017	ggu	sediment output adapted to new framework
! 17.11.2017	ggu	changed VERS_7_5_38
! 27.03.2018	ggu	new code for generic flux handling (fluxes_generic())
! 03.04.2018	ggu	changed VERS_7_5_43
! 06.07.2018	ggu	changed VERS_7_5_48
! 16.02.2019	ggu	changed VERS_7_5_60
! 06.03.2020	ggu	new flux0d, get_barotropic_flux()
! 27.05.2022	ggu	changed to be used with mpi, fluxes now double
! 30.05.2022	ggu	more changes for mpi
! 18.05.2023	ggu	in flx_write() call flx_collect_3d()
! 22.05.2023	ggu	need fluxes_r for write
! 23.10.2024	ggu	definition taken out from flux.f90
! 10.12.2024	ggu	allocation taken out from flux.f90, new name mod_flux
!
!******************************************************************

!==================================================================
        module mod_flux
!==================================================================

        implicit none

	integer, save :: nl_flux = 0
	integer, save :: ns_flux = 0

        logical, save :: bflxinit = .false.
        logical, save :: bflxalloc = .false.
        integer, save :: nsect = -1
        integer, save :: kfluxm = 0
        integer, save, allocatable :: kflux(:)
        integer, save, allocatable :: kflux_ext(:)
        integer, save, allocatable :: iflux(:,:)
        integer, save, allocatable :: itable(:,:)
        character*80, save, allocatable :: chflx(:)
        double precision, save :: da_out(4) = 0.

        integer, save :: nlmax
        integer, save, allocatable :: nlayers(:)
        integer, save, allocatable :: nlayers_global(:)

        double precision, save, allocatable :: fluxes(:,:,:)
        double precision, save, allocatable :: flux0d(:)

        double precision, save, allocatable :: masst(:,:,:)
        double precision, save, allocatable :: saltt(:,:,:)
        double precision, save, allocatable :: tempt(:,:,:)
        double precision, save, allocatable :: conzt(:,:,:)
        double precision, save, allocatable :: ssctt(:,:,:)

        real, save, allocatable :: fluxes_r(:,:,:)

!==================================================================
        contains
!==================================================================

        subroutine flx_alloc_arrays(nl,ns)

        implicit none

        integer nl      !layers
        integer ns      !sections

        if( nl == nl_flux .and. ns == ns_flux ) return

        if( nl > 0 .or. ns > 0 ) then
          if( nl == 0 .or. ns == 0 ) then
            write(6,*) 'nl,ns: ',nl,ns
            stop 'error stop flx_alloc_arrays: incompatible parameters'
          end if
        end if

        !write(6,*) 'flx_alloc_arrays: ',nl,ns

        if( ns_flux > 0 ) then
          deallocate(nlayers)
          deallocate(nlayers_global)
          deallocate(fluxes)
          deallocate(fluxes_r)
          deallocate(flux0d)
          deallocate(masst)
          deallocate(saltt)
          deallocate(tempt)
          deallocate(conzt)
          deallocate(ssctt)
        end if

        nl_flux = nl
        ns_flux = ns

        if( ns == 0 ) return

        allocate(nlayers(ns))
        allocate(nlayers_global(ns))
        allocate(fluxes(0:nl,3,ns))
        allocate(fluxes_r(0:nl,3,ns))
        allocate(flux0d(ns))

        allocate(masst(0:nl,3,ns))
        allocate(saltt(0:nl,3,ns))
        allocate(tempt(0:nl,3,ns))
        allocate(conzt(0:nl,3,ns))
        allocate(ssctt(0:nl,3,ns))

        flux0d = 0.

        bflxalloc = .true.

        end subroutine flx_alloc_arrays

!==================================================================
        end module mod_flux
!==================================================================

