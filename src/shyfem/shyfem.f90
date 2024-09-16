
!--------------------------------------------------------------------------
!
!    Copyright (C) 1995,1997-2020  Georg Umgiesser
!    Copyright (C) 2004,2008  Andrea Cucco
!    Copyright (C) 2006,2008,2011-2012,2014-2015,2014-2015  Christian Ferrarin
!    Copyright (C) 2019  Christian Ferrarin
!    Copyright (C) 2007,2013  Debora Bellafiore
!    Copyright (C) 2012  Aaron Roland
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

! finite element model shyfem (version 3D)
!
! original version from march 1991
!
! revision log :
!
! 30.08.1995	ggu	$$AUST - austausch coefficient introduced
! 11.10.1995	ggu	$$BCLBND - boundary condition for baroclinic runs
! 04.08.1997	ggu	$$ZEONV - new arrays for water level elementwise
! 19.03.1998	ggu	$$IPCCV close data items commented or deleted
! 03.04.1998	ggu	$$DESCRP BUG overwriting descrp with boundary vals.
! 30.04.1998	ggu	finally eliminated /semimp/, /trock/, /ffloat/
! 05.05.1998	ggu	izeit eliminated
! 28.05.1998	ggu	call to sp131g changed
! 14.08.1998	ggu	call to biocos (biological model) introduced
! 20.08.1998	ggu	iextpo finally eliminated, momentum arrays added
! 20.08.1998	ggu	spb11 -> sp111
! 21.08.1998	ggu	xv eliminated
! 03.09.1998	ggu	biological reactor integratated
! 06.11.1998	ggu	hv renamed into hkv
! 08.04.1999	ggu	equilibrium tide (tidal potential) introduced
! 19.04.1999	ggu	subroutine custom called
! 22.10.1999	ggu	igrv, ngrv eliminated (subst by ilinkv, lenkv)
! 20.01.2000	ggu	common block /dimdim/ eliminated
! 07.03.2000	ggu	arrays for vertical turbulent diffusion coefficient
! 20.06.2000	ggu	alv -> visv    slv -> difv
! 19.11.2001	ggu	h1hv eliminated, ulv, vlv eliminated
! 03.12.2001	ggu	new arrays [cst]difhv
! 21.08.2002	ggu	new arrays /kvolc/, /ivol/ copied from hp
! 31.07.2003	ggu	array winv eliminated
! 10.08.2003	ggu	big restructuring
! 13.08.2003	ggu	some more restructuring
! 14.08.2003	ggu	even more restructuring and cleaning up
! 04.12.2003	ggu	integration of wave and sediment module
! 06.03.2004	aac	lagrangian trajectories computation module
! 03.09.2004	ggu	restart now in ht
! 15.10.2004	ggu	gotm substituted by general turbulence closure
! 02.12.2004	ggu	variable time step implemented
! 17.01.2005	ggu	new routine diff_h_set, new difhv, ausv deleted
! 24.02.2005	ggu	new routine smagorinski and common smagv
! 03.03.2005	ggu	smagorinski uses difhv, inserted in subdif.f
! 03.03.2005	ggu	new 3d boundary arrays for C/T/S
! 12.08.2005	ggu	minimum level index, TVD arrays
! 07.11.2005	ggu	new array sed2dn introduced (file name for sediments)
! 16.02.2006	ggu	new routine atoxi3d (TOXI)
! 23.03.2006	ggu	changed time step to real
! 18.10.2006	ccf	radiation stress and waves included (with pipe)
! 10.11.2006	ggu	initialize depth values after restart
! 16.11.2006	ggu	turbulence values included
! 02.04.2007	ggu	changes in algorithm (ASYM)
! 31.05.2007	dbf	new arrays bpresxv, bclevvar (debora)
! 26.09.2007	ggu	deleted arrays rcv,rtv,rsv
! 27.09.2007	ggu	deleted call to tstvol,tstvol1
! 20.03.2008	aac	new call for ERSEM ecological model (BFM MODULE)
! 07.04.2008	aac	new array bfm*bc introduced (file name for ersem)
! 10.04.2008	ggu&ccf	upro, waveov, stokes, z0bk
! 16.04.2008	ggu	evaporation mass flux (evapv)
! 22.04.2008	ggu	gradx/yv non global due to parallelization
! 23.04.2008	ggu	deleted r3{c|t|s}v
! 28.04.2008	ggu	new routine init_stability(), new conzm3sh()
! 29.04.2008	ggu&aac	new BMF ERSEM model
! 12.11.2008	ggu	new sigma level initialization
! 10.12.2008	ggu	new array rfricv for bottom friction
! 18.12.2008	ggu	new routine debug_output(), might be removed later
! 04.03.2009	ggu	matrix amat into solver routines
! 11.03.2009	ggu	new arrays for meteo data
! 06.04.2009	ggu	renamed nlidim to nlkdim
! 30.11.2009	ggu	write output file for successfull completion
! 19.02.2010	ggu	init_stability() changed to reset_stability()
! 22.02.2010	ggu	new arrays wxv, wyv
! 26.02.2010	ggu	new arrays sauxe1/2
! 29.04.2010	ggu	write volumes (wrfvla)
! 26.01.2011	ggu	new arrays for observations and nudging
! 16.02.2011	ggu	new iarnv, call to aquabc
! 17.02.2011	ccf	new radiation stress in 3D
! 23.03.2011	ggu	new call to adjust_spherical()
! 31.03.2011	ggu	write finite volumes at initial time step
! 20.05.2011	ggu	iwetv introduced, wet and dry from main
! 25.10.2011	ggu	hlhv eliminated
! 18.11.2011	ggu	new routine handle_projection
! 24.01.2012	ggu	new call to setup_parallel()
! 23.02.2012	ggu&ccf	meteo arrays adjusted (3*nkn)
! 09.03.2012	ggu	call to residence time added
! 21.06.2012	ggu&aar	fluid mud variables integrated
! 05.08.2012	ggu	bug because lam2dn and dmfd2n not defined
! 10.05.2013	dbf	initialization for non hydrostatic routines
! 13.06.2013	ggu	set/copydepth simplified, offline version
! 05.09.2013	ggu	changed order of adjust depth and barene structures
! 29.10.2013	ggu	nudging implemented
! 25.03.2014	ggu	new offline
! 25.06.2014	ggu	new arrays hkv_min, hkv_max
! 05.12.2014	ccf	new interface for waves
! 30.07.2015	ggu	routine renamed from ht to shyfem
! 18.09.2015	ggu	new routine scalar, call to hydro()
! 23.09.2015	ggu	changed VERS_7_2_4
! 29.09.2015	ccf	inverted set_spherical() and handle_projection()
! 30.09.2015	ggu	changed VERS_7_2_6
! 10.10.2015	ggu	fluid mud routines handled differently
! 12.10.2015	ggu	changed VERS_7_3_3
! 13.10.2015	ggu	changed VERS_7_3_5
! 22.10.2015	ggu	changed VERS_7_3_7
! 23.10.2015	ggu	changed VERS_7_3_9
! 05.11.2015	ggu	changed VERS_7_3_12
! 09.11.2015	ggu	changed VERS_7_3_13
! 20.11.2015	ggu	changed VERS_7_3_15
! 16.12.2015	ggu	changed VERS_7_3_16
! 19.02.2016	ggu	changed VERS_7_5_2
! 22.02.2016	ggu	changed VERS_7_5_4
! 07.06.2016	ggu	changed VERS_7_5_12
! 14.06.2016	ggu	changed VERS_7_5_14
! 13.04.2017	ggu	changed VERS_7_5_25
! 09.05.2017	ggu	changed VERS_7_5_26
! 16.05.2017	ggu	changed VERS_7_5_27
! 05.10.2017	ggu	command line options introduced, subs rearranged
! 14.11.2017	ggu	changed VERS_7_5_36
! 17.11.2017	ggu	changed VERS_7_5_37
! 05.12.2017	ggu	changed VERS_7_5_39
! 07.12.2017	ggu	changed VERS_7_5_40
! 03.04.2018	ggu	changed VERS_7_5_43
! 19.04.2018	ggu	changed VERS_7_5_45
! 26.04.2018	ggu	changed VERS_7_5_46
! 11.05.2018	ggu	changed VERS_7_5_47
! 06.07.2018	ggu	changed VERS_7_5_48
! 25.10.2018	ggu	changed VERS_7_5_51
! 18.12.2018	ggu	changed VERS_7_5_52
! 12.02.2019	ccf	bottom shear stress in substress.f
! 16.02.2019	ggu	changed VERS_7_5_60
! 12.03.2019	ccf	include new computation of tide potential/analysis
! 04.07.2019	ggu	new ww3 routines introduced
! 15.09.2019	ggu	subroutine to test only forcing
! 02.10.2019	ggu	delete include files
! 17.10.2019	ggu	no call to bfm_write, is done inside subroutine
! 06.11.2019	ggu	femelab eliminated
! 03.04.2020	ggu	write real start and end time of simulation
! 09.04.2020    ggu     run bfm through bfm_run()
! 21.05.2020    ggu     better handle copyright notice
! 04.06.2020    ggu     debug_output() substituted with shympi_debug_output()
! 30.03.2021    ggu     more on debug, call sp111(2) outside time loop
! 01.04.2021    ggu     turbulence cleaned
! 02.04.2022    ggu     new option -mpi_debug -> calls shympi_check_all()
! 02.04.2022    ggu     new routine shympi_write_debug_special()
! 03.04.2022    ggu     timing problems in handle_debug_output() solved
! 11.04.2022    ggu     no -mpi switch necessary anymore
! 12.04.2022    ggu     message to show if mpi support is available
! 18.05.2022    ggu     cpu_time routines introduced
! 10.03.2023    ggu     do not use bmpirun anymore
! 28.04.2023    ggu     update function calls for belem
! 22.05.2023    ggu     new names for closing: close_init, close_handle
! 05.06.2023    lrp     introduce z-star
! 11.12.2023    ggu     f90 style format
! 11.12.2023    ggu     create subroutines for init/run/finalize
! 12.12.2023    ggu     introduce dtmax to make limited run
! 10.05.2024    ggu     set spherical just after basin read
! 06.09.2024    lrp     nuopc-compliant
!
!*****************************************************************

	program shyfem
	call shyfem_main
	end program

!*****************************************************************
