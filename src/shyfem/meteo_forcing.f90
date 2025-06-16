
!--------------------------------------------------------------------------
!
!    Copyright (C) 2009-2012,2014-2020  Georg Umgiesser
!    Copyright (C) 2009  Debora Bellafiore
!    Copyright (C) 2012,2017  Christian Ferrarin
!    Copyright (C) 2016  Ivan Federico
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

! handle meteo files with new fem format
!
! revision log :
!
! 10.03.2009	ggu	finished coding
! 24.03.2009	ggu	use metrain instead rqdsv
! 07.05.2009	ggu	new call to init_coords()
! 18.06.2009	ggu&dbf	bug fix -> wind speed not interpolated on metws
! 23.02.2010	ggu	call to wstress changed (new wxv,wyv)
! 26.01.2011	ggu	write wind field to debug output (iumetw)
! 05.02.2011	ggu	changed order of records in dfile (actual is first)
! 16.02.2011	ggu	pass idata to files, use geo info from files
! 18.11.2011	ggu	deleted projection code from subroutines
! 10.02.2012	ggu	limit cc and rh to acceptable range
! 16.02.2012	ggu	new routine meteo_get_solar_radiation()
! 22.02.2012	ggu	new routines for regular and ts reading of meteo
! 23.02.2012	ggu&ccf	bug fix meteo_copy_to_old and meteo_interpolate_in_time
! 20.05.2014	ggu	new routines for new file format
! 30.05.2014	ggu	changed VERS_6_1_76
! 18.06.2014	ggu	changed VERS_6_1_77
! 27.06.2014	ggu	changed VERS_6_1_78
! 07.07.2014	ggu	changed VERS_6_1_79
! 18.07.2014	ggu	changed VERS_7_0_1
! 30.10.2014	ggu	changed VERS_7_0_4
! 05.11.2014	ggu	changed VERS_7_0_5
! 07.11.2014	ggu	changed VERS_7_0_6
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 09.01.2015	ggu	changed VERS_7_0_12
! 19.01.2015	ggu	changed VERS_7_1_3
! 26.02.2015	ggu	changed VERS_7_1_5
! 30.04.2015	ggu	ice integrated
! 04.05.2015	ggu	bug in ice eliminated
! 12.05.2015	ggu	introduced ia_icefree for icefree elements
! 21.05.2015	ggu	changed VERS_7_1_11
! 05.06.2015	ggu	changed VERS_7_1_12
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 29.09.2015	ggu	changed VERS_7_2_5
! 08.01.2016	ggu	bug fix in meteo_convert_wind_data() - no wind bug
! 19.02.2016	ggu	changed VERS_7_5_2
! 10.03.2016	ggu	check for pressure to be in reasonable bounds
! 01.04.2016	ggu	changed VERS_7_5_7
! 15.04.2016	ggu	changed VERS_7_5_8
! 10.06.2016	ggu	changed VERS_7_5_13
! 17.06.2016	ggu	changed VERS_7_5_15
! 27.06.2016	ggu	changed VERS_7_5_16
! 23.07.2016	ivn	new heat formulation for iheat==8
! 09.09.2016	ggu	new variable ihtype to choose between rh, wbt, dpt
! 16.09.2016	ggu	allow for Pa and mbar in pressure
! 30.09.2016	ggu	changed VERS_7_5_18
! 12.01.2017	ccf	bug fix in determining pressure units
! 13.06.2017	ggu	changed VERS_7_5_29
! 02.09.2017	ggu	changed VERS_7_5_31
! 05.12.2017	ggu	changed VERS_7_5_39
! 24.01.2018	ggu	changed VERS_7_5_41
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.04.2018	ggu	changed VERS_7_5_43
! 06.07.2018	ggu	changed VERS_7_5_48
! 31.08.2018	ggu	changed VERS_7_5_49
! 03.10.2018	ggu	better output of meteo variables (output_meteo_data)
! 16.10.2018	ggu	changed VERS_7_5_50
! 27.12.2018	ggu	changed VERS_7_5_54
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 17.10.2019	ggu	check number of meteo variables written
! 09.12.2019	ggu	more documentation
! 27.01.2020	ggu	code to use full ice cover (ballcover)
! 17.04.2020	ggu	wind conversion routines out of this file
! 11.11.2020	ggu	new routine meteo_has_ice_file()
! 03.06.2021	mbj	added Hersbach wind stress formulation
! 26.01.2022	ggu	bug in short name of icecover fixed
! 01.02.2022	ggu	automatically convert cloudcover from % to fraction
! 21.03.2022	ggu	new calls for write in debug mode
! 03.06.2022	ggu	in meteo_convert_heat_data() save read vapor to array
! 08.07.2022	ggu	avoid divide by zero when computing dice
! 06.12.2022	ggu	rfact for rain introduced
! 02.04.2023    ggu     only master writes to iuinfo
! 06.09.2024    lrp     nuopc-compliant
! 13.09.2024    lrp     iatm and coupling with atmospheric model
! 22.09.2024    ggu     read meteo output times and from str file
! 13.10.2024    ggu     deal with INTEL_BUG
! 13.11.2024    ggu     marked old code with INTEL_BUG_OLD
! 03.12.2024    lrp     new parameter irain for the coupled model
! 25.01.2025    ggu     tentative cubic interpolation for wind (bbspline)
!
! notes :
!
! info on file format can be found in subfemfile.f
!
! contents :
!
! meteo_set_wind_data
! meteo_convert_wind_data
! meteo_set_rain_data
! meteo_convert_rain_data
! meteo_set_ice_data
! meteo_convert_ice_data
! meteo_set_heat_data
! meteo_convert_heat_data
!
!*********************************************************************
!
! DOCS  START   S_wind
!
! In this section the wind data can be given directly without
! the creation of an external file. Note, however, that
! a wind file specified in the |name| section takes precedence
! over this section. E.g., if both a section |wind| and a
! wind file in |name| is given, the wind data from the file is used.
!
! The format of the wind data in this section is the same as the
! format in the ASCII wind file, i.e., three columns, with
! the first specifying the time in seconds and the other two columns
! giving the wind data. The interpretation of the wind data
! depends on the value of |iwtype|. For more information please
! see the description of |iwtype| in section |para|.
!
! DOCS  END

! DOCS  START   P_wind
!
! DOCS  WIND            Wind parameters
!
! The next two parameters deal with the wind stress to be
! prescribed at the surface of the basin.
!
! The wind data can either be specified in an external file (ASCII
! or binary) or directly in the parameter file in section |wind|.
! The ASCII file or the wind section contain three columns, the first
! giving the time in seconds, and the others the components of
! the wind speed. Please see below how the last two columns are
! interpreted depending on the value of |iwtype|. For the format
! of the binary file please see the relative section.
! If both a wind file and section |wind| are given, data from the
! file is used.
!
! The wind stress is normally computed with the following formula
! \beq
! \tau^x = \rho_a c_D \vert u \vert u^x \quad
! \tau^y = \rho_a c_D \vert u \vert u^y
! \eeq
! where $\rho_a,\rho_0$ is the density of air and water respectively,
! $u$ the modulus of wind speed and $u^x,u^y$ the components of
! wind speed in $x,y$ direction. In this formulation $c_D$ is a
! dimensionless drag coefficient that varies between 1.5 \ten{-3} and
! 3.2 \ten{-3}. The wind speed is normally the wind speed measured
! at a height of 10 m.
!
! |iwtype|      The type of wind data given (default 1):
!               \begin{description}
!               \item[0] No wind data is processed
!               \item[1] The components of the wind is given in [m/s]
!               \item[2] The stress ($\tau^x,\tau^y$) is directly specified
!               \item[3] The wind is given in speed [m/s] and direction
!                        [degrees]. A direction of 0\degrees{} specifies
!                        a wind from the north, 90\degrees{} a wind
!                        from the east etc.
!               \item[4] As in 3 but the speed is given in knots
!               \end{description}
!
! |itdrag|	Formula to compute the drag coefficient. 
!		\begin{description}
!		\item[0] constant value given in |dragco|. 
!		\item[1] Smith and Banke (1975) formula
!		\item[2] Large and Pond (1981) formula
!		\item[3] Spatial/temporal varying in function of wave. Need
!		the coupling with WWMIII.
!		\item[4] Spatial/temporal varying in function of heat flux. 
!		Only for |iheat| = 6. 
!		\item[5] Hersbach (2011) formula. This is a fit of kinematic
!		viscosity for light wind and a Charnock coefficient for strong 
!		wind. Neutral wind should be used (small differences).
!		See the paper and/or the ECMWF report.
!		\end{description}
!		(Default 0)
!
! |dragco|	Drag coefficient used in the above formula. (Default 2.5E-3). 
!		If |itdrag| = 5 this is the Charnock parameter and you should 
!		use values from 0.01 (swell) to 0.04 (steep young waves). 
!		(Default 0.025).
!		Please note that in case of |iwtype| = 2 this parameter
!		is of no interest, since the stress is specified directly.
!
! |wsmax|	Maximum wind speed allowed in [m/s]. This is in order to avoid
!		errors if the wind data is given in a different format
!		from the one specified by |iwtype|. (Default 50)
!
! |wslim|	Limit maximum wind speed to this value [m/s]. This provides
!		an easy way to exclude strong wind gusts that might
!		blow up the simulation. Use with caution. 
!		(Default -1, no limitation)
!
! |rfact|	Precipitation (rain) has to be given in mm/day. If the
!		input data is in a different unit, |rfact| specifies 
!		the conversion factor. E.g., if the data is in mm/hour,
!		|rfact = 24| converts it to mm/day.
!		(Default 1)
!
! DOCS  END

!================================================================
        module meteo_forcing_module
!================================================================

	use intp_fem_file

	implicit none

	real, parameter :: pstd = 101325.	!standard pressure
	real, parameter :: tkelvin = 273.15	!0C in Kelvin
	real, parameter :: dstd = 2.5e-3	!standard drag coefficient
	real, parameter :: cstd = 0.025		!standard Charnock coefficient
	real, parameter :: nmile = 1852.	!nautical mile in m

	real, parameter :: amice = 1.	!use momentum reduction due to ice
					!1: use  0: do not reduce momentum

	integer, save :: idwind,idheat,idrain,idice

	integer, save :: ndbgfreq = 0			!debug output
	double precision, save, private :: da_out(4) = 0
	double precision, save, private :: da_met(4) = 0

	integer, save :: iatm = 0
	integer, save :: iheat = 0
	integer, save :: iwtype,itdrag
        integer, save :: ievap,irain
	integer, save :: irtype
	integer, save :: ihtype
	integer, save :: ictype
	integer, save :: ia_icefree		!area type which is ice free
	real, save :: wsmax,wslim,dragco,roluft,rowass
	real, save :: pfact = 1.
	real, save :: wfact = 1.
	real, save :: sfact = 1.
	real, save :: rfact = 1. !factor for rain if other units than mm/day

	logical, save :: has_pressure = .false.

	logical, save, private :: bdebug = .true.

	integer, save, private :: icall = 0

	character*80, save :: wxss = 'wind stress - x [N/m**2]'
	character*80, save :: wyss = 'wind stress - y [N/m**2]'
	character*80, save :: wxms = 'wind velocity - x [m/s]'
	character*80, save :: wyms = 'wind velocity - y [m/s]'
	character*80, save :: wsms = 'wind speed [m/s]'
	character*80, save :: wskn = 'wind speed [knots]'
	character*80, save :: wdir = 'wind direction [deg]'

	character*80, save :: papa = 'pressure (atmospheric) [Pa]'
	character*80, save :: pamb = 'pressure (atmospheric) [mbar]'

	character*80, save :: rain = 'rain [mm/day]'
	character*80, save :: ice  = 'ice cover [0-1]'

        character*80, save :: srad = 'solar radiation [W/m**2]'
        character*80, save :: tair = 'air temperature [C]'
        character*80, save :: rhum = 'humidity (relative) [%]'
        character*80, save :: shum = 'humidity (specific)'
        character*80, save :: ccov = 'cloud cover [0-1]'
        character*80, save :: wbtm = 'wet bulb temperature [C]'
        character*80, save :: dewp = 'dew point temperature [C]'    

!================================================================
        contains
!================================================================

!*********************************************************************

	subroutine meteo_forcing_fem

! administers meteo files (still to be cleaned)

	use mod_meteo
	use mod_depth
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

	character*60 windfile,heatfile,rainfile,icefile
	character*4 what

	integer nvar,lmax
	integer nintp
	integer modehum
	integer nvarm,nid,nlev
	integer i
	!integer it0
	double precision dtime0,dtime
	real flag
	real val0,val1,val2
	real metaux(nkn)

	real vconst(4)
	integer nodes(1)
	logical batm

	logical bbspline	!interpolate wind with 4th order bspline
	logical bbwrite		!write spline output to file fort.77
	real getpar

!------------------------------------------------------------------
! batm is true if first ocean timestep of coupled atm-oce timestep
!------------------------------------------------------------------

	bbspline = .false.
	bbspline = .true.
	bbwrite = .true.
	bbwrite = .false.

	batm = iatm == 1 .and. icall_nuopc == 1

	if( icall .lt. 0 ) return

!------------------------------------------------------------------
! initialization
!------------------------------------------------------------------

	call shympi_barrier

	if( icall .eq. 0 ) then

	  write(6,*) 'initialization of meteo forcing fem'

!	  ---------------------------------------------------------
!	  initialization of data files
!	  ---------------------------------------------------------

	  call getfnm('wind',windfile)
	  call getfnm('qflux',heatfile)
	  call getfnm('rain',rainfile)
	  call getfnm('ice',icefile)

	  call get_first_dtime(dtime0)

!	  ---------------------------------------------------------
!	  initializing wind file
!	  ---------------------------------------------------------

	  write(6,'(a)') 'opening wind file...'
          nvar = 0      !not sure if 2 or 3
          call iff_get_file_nvar(windfile,nvar)
          if( nvar <= 0 ) nvar = 3      !if no file fake 3
          nintp = 2
          if( bbspline .and. nint(getpar('ibwind')) > 0 ) nintp = 4
          what = 'wind'
          vconst = (/ 0., 0., pstd, 0. /)
	  call iff_init(dtime0,windfile,nvar,nkn,0,nintp,nodes,vconst,idwind)
	  call iff_set_description(idwind,0,'meteo wind')

	  call meteo_set_wind_data(idwind,nvar)

!	  ---------------------------------------------------------
!	  initializing ice file
!	  ---------------------------------------------------------

	  write(6,'(a)') 'opening ice file...'
	  nvar = 1
	  nintp = 2
	  what = 'ice'
	  vconst = (/ 0., 0., 0., 0. /)
	  dtime = -1.
	  call iff_init(dtime,icefile,nvar,nkn,0,nintp,nodes,vconst,idice)
	  call iff_set_description(idice,0,'meteo ice')
	  call iff_flag_ok(idice)	!we do not need all ice data
	  if( iff_has_file(idice) ) then
	    call iff_populate_records(idice,dtime0)  !only now because of flag
	  end if

	  call meteo_set_ice_data(idice,nvar)

!	  ---------------------------------------------------------
!	  initializing rain file
!	  ---------------------------------------------------------

	  write(6,'(a)') 'opening rain file...'
	  nvar = 1
	  nintp = 2
	  what = 'rain'
	  vconst = (/ 0., 0., 0., 0. /)
	  call iff_init(dtime0,rainfile,nvar,nkn,0,nintp,nodes,vconst,idrain)
	  call iff_set_description(idrain,0,'meteo rain')

	  call meteo_set_rain_data(idrain,nvar)

!	  ---------------------------------------------------------
!	  initializing heat file
!	  ---------------------------------------------------------

	  write(6,'(a)') 'opening heat flux file...'
	  nvar = 4
	  nintp = 2
	  what = 'heat'
	  vconst = (/ 0., 0., 50., 0. /)
	  call iff_init(dtime0,heatfile,nvar,nkn,0,nintp,nodes,vconst,idheat)
	  call iff_set_description(idheat,0,'meteo heat')

	  call meteo_set_heat_data(idheat,nvar)

	end if

!------------------------------------------------------------------
! end of initialization
!------------------------------------------------------------------

	icall = icall + 1

!------------------------------------------------------------------
! time interpolation
!------------------------------------------------------------------

	call get_act_dtime(dtime)
	lmax = 1

!	---------------------------------------------------------
!	read and interpolate ice file
!	---------------------------------------------------------

	if( .not. iff_is_constant(idice) .or. icall == 1 ) then
	  call iff_read_and_interpolate(idice,dtime)
	  metice = 0.	!assume no ice cover in areas not covered by data
	  call iff_time_interpolate(idice,dtime,1,nkn,lmax,metice)
	end if

!	---------------------------------------------------------
!	read and interpolate wind file
!	---------------------------------------------------------

	if( .not. iff_is_constant(idwind) .or. icall == 1 ) then
	  call iff_read_and_interpolate(idwind,dtime)
	  call iff_time_interpolate(idwind,dtime,1,nkn,lmax,wxv)
	  call iff_time_interpolate(idwind,dtime,2,nkn,lmax,wyv)
	  if( iff_get_nvar(idwind) == 3 ) then
	    call iff_time_interpolate(idwind,dtime,3,nkn,lmax,ppv)
	  else
	    ppv = pstd
	  end if
	end if

!	---------------------------------------------------------
!	read and interpolate rain file
!	---------------------------------------------------------

        if( .not. iff_is_constant(idrain) .or. icall == 1 ) then
          call iff_read_and_interpolate(idrain,dtime)
          call iff_time_interpolate(idrain,dtime,1,nkn,lmax,metrain)
        end if

!	---------------------------------------------------------
!	read and interpolate heat file
!	---------------------------------------------------------

        if( .not. iff_is_constant(idheat) .or. icall == 1 ) then
          call iff_read_and_interpolate(idheat,dtime)
          call iff_time_interpolate(idheat,dtime,1,nkn,lmax,metrad)
          call iff_time_interpolate(idheat,dtime,2,nkn,lmax,mettair)
          call iff_time_interpolate(idheat,dtime,3,nkn,lmax,metaux)
          call iff_time_interpolate(idheat,dtime,4,nkn,lmax,metcc)
        end if

!------------------------------------------------------------------
! extra treatment of data
!------------------------------------------------------------------

!	---------------------------------------------------------
!	treat ice data
!	---------------------------------------------------------

	if( .not. iff_is_constant(idice) .or. icall == 1 ) then
	  call meteo_convert_ice_data(idice,nkn,metice)
	end if

!	---------------------------------------------------------
!	treat wind data
!	---------------------------------------------------------

	if( .not. iff_is_constant(idwind) .or. icall == 1 .or. batm ) then
	  call meteo_convert_wind_data(idwind,nkn,wxv,wyv     &
      &			,windcd,tauxnv,tauynv,metws,ppv,metice)
	end if

!	---------------------------------------------------------
!	treat heat data
!	---------------------------------------------------------

        if( .not. iff_is_constant(idheat) .or. icall == 1 .or. batm ) then
          call meteo_convert_heat_data(idheat,nkn                 &
      &                       ,metaux,mettair,metcc,ppv,methum)
        end if

!	---------------------------------------------------------
!	treat rain data
!	---------------------------------------------------------

	if( .not. iff_is_constant(idrain) .or. icall == 1 .or. batm ) then
	  call meteo_convert_rain_data(idrain,nkn,metrain)
	end if

	if( bbspline ) then
	  if( bbwrite .and. icall > 1 ) then
	    write(77,*) dtime,metws(1),wxv(1),wyv(1)
	  end if
	end if

!------------------------------------------------------------------
! output
!------------------------------------------------------------------

	call output_debug_data
	call output_meteo_data

	if( bextra_exchange ) then
	  call shympi_exchange_2d_node(tauxnv)
	  call shympi_exchange_2d_node(tauynv)
	end if

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

	return
   99	continue
	stop 'error stop meteo_regular: need wind data for heat module'
	end subroutine meteo_forcing_fem

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine output_debug_data

	use mod_meteo

	integer, parameter :: nvarm = 4	!total number of meteo vars written
	integer, save :: icall = 0
	integer, save :: id = 0
	double precision dtime

	if( ndbgfreq .le. 0 ) return

	if( icall == 0 ) then
	  call shyfem_init_scalar_file('meteo',nvarm,.true.,id)
	end if

	icall = icall + 1

        if( mod(icall,ndbgfreq) .eq. 0 ) then
	  call get_act_dtime(dtime)
	  call shy_write_scalar_record2d(id,dtime,20,ppv)
	  call shy_write_scalar_record2d(id,dtime,28,metws)
	  call shy_write_scalar_record2d(id,dtime,23,mettair)
	  call shy_write_scalar_record2d(id,dtime,85,metice)
        end if

	end subroutine output_debug_data

!*********************************************************************

	subroutine output_meteo_data

	use basin
	use mod_meteo
	use shympi

	integer			:: id
	integer			:: nvar_act
	integer, save		:: imetout
	double precision 	:: dtime
	integer, save		:: nvar = 0
	logical, save 		:: b2d = .true.
	logical, save		:: bwind,bheat,brain,bice
        real, parameter		:: zconv = 86400. / 1000. !convert m/s to mm/day
	real, allocatable	:: maux(:)

	logical has_output_d,next_output_d
	integer bit10_extract_value
	real getpar

	if( da_met(4) < 0 ) return

	if( da_met(4) == 0 ) then
	  call init_output_d('itmmet','idtmet',da_met)
	  if( .not. has_output_d(da_met) ) da_met(4) = -1	!no output
	  if( da_met(4) < 0 ) return

	  imetout = getpar('imetout')	!what type of meteo output
	  bwind = ( bit10_extract_value(imetout,1) > 0 )
	  bheat = ( bit10_extract_value(imetout,2) > 0 )
	  brain = ( bit10_extract_value(imetout,3) > 0 )
	  bice  = ( bit10_extract_value(imetout,4) > 0 )
	  bwind = bwind .and. iff_has_file(idwind)
	  bheat = bheat .and. iff_has_file(idheat)
	  brain = brain .and. iff_has_file(idrain)
	  bice  = bice  .and. iff_has_file(idice)

	  nvar = 0
	  if( bwind ) nvar = nvar + 4
	  if( bheat ) nvar = nvar + 4
	  if( brain ) nvar = nvar + 1
	  if( bice  ) nvar = nvar + 1
	  if( nvar == 0 ) da_met(4) = -1
	  if( da_met(4) < 0 ) return

          call shyfem_init_scalar_file('meteo',nvar,b2d,id)
          da_met(4) = id

	  write(6,*) 'meteo output: ',imetout,nvar,bwind,bheat,brain,bice,my_id

	end if

        if( .not. next_output_d(da_met) ) return

	call get_act_dtime(dtime)
        id = nint(da_met(4))

	call shy_reset_nvar_act(id)

	if( bwind ) then
	  call shy_write_scalar_record2d(id,dtime,21,wxv)
	  call shy_write_scalar_record2d(id,dtime,21,wyv)
	  call shy_write_scalar_record2d(id,dtime,28,metws)
	  call shy_write_scalar_record2d(id,dtime,20,ppv)
	end if
	if( bheat ) then
          call shy_write_scalar_record2d(id,dtime,22,metrad)
          call shy_write_scalar_record2d(id,dtime,23,mettair)
          call shy_write_scalar_record2d(id,dtime,24,methum)
          call shy_write_scalar_record2d(id,dtime,25,metcc)
	end if
	if( brain ) then
	  allocate(maux(nkn))
	  maux = metrain * zconv
          call shy_write_scalar_record2d(id,dtime,26,maux)
	end if
	if( bice ) then
          call shy_write_scalar_record2d(id,dtime,85,metice)
	end if

	call shy_get_nvar_act(id,nvar_act)

	if( nvar /= nvar_act ) then
	  write(6,*) 'nvar,nvar_act: ',nvar,nvar_act
	  write(6,*) 'number of variables written differs from nvar'
	  stop 'error stop output_meteo_data: nvar /= nvar_act'
	end if

	call shy_sync(id)

	end subroutine output_meteo_data

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine meteo_set_wind_data(id,nvar)

! iwtype = 0 no wind data processed
! iwtype = 1 wind data in format (wx,wy)
! iwtype = 2 wind data in format (tx,ty) (stress)
! iwtype = 3 wind data in format (speed,dir) (speed in m/s)
!          dir = 0  -> wind from north   dir = 90 -> wind from east
! iwtype = 4 wind data in format (speed,dir) (speed in nodes)

	integer id
	integer nvar

	integer il
	character*80 string,string1,string2
	character*10 dir,unit

	logical string_is_this_short
	real, external :: getpar

!	---------------------------------------------------------
!	check nvar and get parameters
!	---------------------------------------------------------

	if( iff_has_file(id) ) then
	 if( nvar /= 2 .and. nvar /= 3 ) then
	  write(6,*) 'no support for nvar = ',nvar
	  stop 'error stop meteo_set_wind_data: wind'
	 end if
	end if

        iatm = nint(getpar('iatm'))
        iwtype = nint(getpar('iwtype'))
        itdrag = nint(getpar('itdrag'))
        wsmax = getpar('wsmax')
        wslim = getpar('wslim')
        dragco = getpar('dragco')
        roluft = getpar('roluft')
        rowass = getpar('rowass')

	if( dragco .lt. 0 ) dragco = dstd
	if( dragco == dstd .and. itdrag .eq. 5 ) dragco = cstd

!	---------------------------------------------------------
!	handle wind
!	---------------------------------------------------------

	call iff_get_var_description(id,1,string1)
	call iff_get_var_description(id,2,string2)

	if( .not. iff_has_file(id) ) then	!no wind file

	  if ( iatm == 1 ) then			!no file but coupling atm-oce
            if( iwtype .ne. 0 .and. iwtype .ne. 2 ) then
              write(6,*) 'atmosphere-ocean coupling but no wind type'
              write(6,*) 'set iwtype = 0 or iwtype = 2'
	      stop 'error stop meteo_set_wind_data: iwtype'
            end if
	  else
            iwtype = 0                          !no file and no coupling atm-oce
	  end if

	else if( string1 == ' ' ) then		!TS file or constant

	  if( iff_has_file(id) ) then
	    if( iwtype .le. 0 ) then
	      write(6,*) 'wind file given but no wind type available'
	      write(6,*) 'either set wind description or set iwtype'
	      stop 'error stop meteo_set_wind_data: no wind type'
	    else if( iwtype > 4 ) then
	      write(6,*) 'not supported value for iwtype = ',iwtype
	      stop 'error stop meteo_set_wind_data: iwtype'
	    end if
	  else				!no file opened
	    iwtype = 0
	  end if

	  if( iwtype == 1 ) then
	    call iff_set_var_description(id,1,wxms)
	    call iff_set_var_description(id,2,wyms)
	  else if( iwtype == 2 ) then
	    call iff_set_var_description(id,1,wxss)
	    call iff_set_var_description(id,2,wyss)
	  else if( iwtype == 3 ) then
	    call iff_set_var_description(id,1,wsms)
	    call iff_set_var_description(id,2,wdir)
	  else if( iwtype == 4 ) then
	    call iff_set_var_description(id,1,wskn)
	    call iff_set_var_description(id,2,wdir)
	  end if

	else	!description given -> check and set iwtype

          if( string_is_this_short('wind',string1) .and.             &
      &		string_is_this_short('wind',string2) ) then
	    iwtype = 1
          else if( string_is_this_short('wstress',string1) .and.     &
      &		string_is_this_short('wstress',string2) ) then
	    iwtype = 2
          else if( string_is_this_short('windspeed',string1) .and.   &
      &		string_is_this_short('winddir',string2) ) then
	    iwtype = 3
	  else
	    write(6,*) 'description string for wind not recognized: '
	    write(6,*) trim(string1)
	    write(6,*) trim(string2)
	    stop 'error stop meteo_set_wind_data: wind description'
	  end if

	  if( iwtype == 3 ) then
	    call string_direction_and_unit(string1,dir,unit)
	    if( unit == 'knots' ) iwtype = 4
	  end if

	end if

	wfact = 1. / rowass
	if( iwtype /= 2 ) wfact = roluft / rowass
	sfact = 1.
	if( iwtype == 4 ) sfact = nmile / 3600.

!	---------------------------------------------------------
!	handle pressure
!	---------------------------------------------------------

	if( nvar == 3 ) then
	  has_pressure = .true.
	  call iff_get_var_description(id,3,string)
          if( string_is_this_short('airp',string) ) then
	    call string_direction_and_unit(string1,dir,unit)
	    if( unit == 'Pa' ) then
	      pfact = 1.
	    else if( unit == 'mbar' ) then
	      pfact = 100.
	    else if( unit == ' ' ) then		!must determine later
	      pfact = 0.			!check later
	    end if
	  else if( string == ' ' ) then
	    pfact = 1.
	  else
	    write(6,*) 'description for pressure not recognized: '
	    write(6,*) trim(string)
	    stop 'error stop meteo_set_wind_data: press description'
	  end if
	  if( bdebug ) then
	    write(6,*) 'pressure initialized: ',pfact,trim(string)
	    write(6,*) 'wind type: ',iwtype
	  end if
	end if

!	---------------------------------------------------------
!	remember values and write to monitor
!	---------------------------------------------------------

        call putpar('iwtype',real(iwtype))

	if( iwtype == 0 ) then
	  write(6,*) 'no wind file opened'
	else
	  write(6,*) 'wind file opened: ',iwtype
	  call iff_get_var_description(id,1,string1)
	  call iff_get_var_description(id,2,string2)
	  write(6,*) 'content: '
	  write(6,*) ' 1    ',trim(string1)
	  write(6,*) ' 2    ',trim(string2)
	  if( nvar == 3 ) then
	    call iff_get_var_description(id,3,string)
	    write(6,*) ' 3    ',trim(string)
	  end if
	end if

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end subroutine meteo_set_wind_data

!*********************************************************************

	subroutine meteo_convert_wind_data(id,n,wx,wy,cdv,tx,ty,ws,pp,cice)

	use basin
	use shympi

	integer id
	integer n
	real wx(n),wy(n)
	real cdv(n)
	real tx(n),ty(n)
	real ws(n)
	real pp(n)
	real cice(n)

	logical bnowind,bstress,bspeed
	integer k
	real cd,wxymax,txy,wspeed,wdir,fact,fice,aice,ach
	real pmin,pmax
	real wparam
	double precision dwxy		!INTEL_BUG
	character*80 string
	double precision dtime

	bnowind = iwtype == 0
	bstress = iwtype == 2
	bspeed = iwtype > 2
	cd = dragco
	ach = dragco       !Charnock parameter (ireib=5)
	wxymax = 0.
	aice = amice       !ice cover for momentum: 1: use  0: do not use
	
	call get_act_dtime(dtime)
	
!	---------------------------------------------------------
!	convert wind
!	---------------------------------------------------------

        if( bnowind ) then              !no wind
	  ws = 0.
	  wx = 0.
	  wy = 0.
	  tx = 0.
	  ty = 0.
	  pp = pstd
        else if( bstress ) then         !data is stress -> normalize it
          if( cd .le. 0 ) cd = dstd
          do k=1,n
	    fice = 1. - aice*cice(k)
            tx(k) = fice * wfact * wx(k)
            ty(k) = fice * wfact * wy(k)
            txy = sqrt( tx(k)**2 + ty(k)**2 )
            wspeed = sqrt(txy/(cd*wfact*roluft))
            wxymax = max(wxymax,wspeed)
            wx(k) = tx(k) / (cd*wspeed)
            wy(k) = ty(k) / (cd*wspeed)
	    ws(k) = wspeed
          end do
        else 
	  if( bspeed ) then             !data is speed/dir
            do k=1,n
	      wspeed = sfact * wx(k)
	      wdir = wy(k)
              call convert_sd_uv(wspeed,wdir,wx(k),wy(k),.true.)
	      ws(k) = wspeed
	    end do
	  else				!data is wind velocity [m/s]
            do k=1,n
	      dwxy = dble(wx(k))**2 + dble(wy(k))**2
              wspeed = sqrt( dwxy )			!INTEL_BUG
	      !wspeed = sqrt( wx(k)**2 + wy(k)**2 )	!INTEL_BUG_OLD
	      ws(k) = wspeed
	    end do
	  end if

          do k=1,n
	    fice = 1. - aice*cice(k)
	    !if( k .eq. 1000 ) write(6,*) 'ice: ',k,fice
            wspeed = ws(k)
            wxymax = max(wxymax,wspeed)
	    cd = cdv(k)
            if( itdrag .gt. 0 .and. ( itdrag == 1           & 
      &          .or. itdrag == 2 .or. itdrag == 5 )) then
		call get_drag(itdrag,wspeed,cd,ach)
	    end if
	    wparam = fice * wfact * cd * wspeed		!INTEL_BUG_OLD
            tx(k) = wparam * wx(k)
            ty(k) = wparam * wy(k)
          end do
        end if

!	---------------------------------------------------------
!	limit wind
!	---------------------------------------------------------

	if( wslim > 0 .and. wxymax > wslim ) then !artificially limit wind speed
          do k=1,n
	    fice = 1. - aice*cice(k)
            wspeed = ws(k)
	    if( wspeed <= wslim ) cycle
	    ws(k) = wslim
	    fact = wslim/wspeed
	    wx(k) = fact*wx(k)
	    wy(k) = fact*wy(k)
            wspeed = ws(k)
	    cd = cdv(k)
            if( itdrag .gt. 0 .and. ( itdrag == 1           &
      &          .or. itdrag == 2 .or. itdrag == 5 )) then
		call get_drag(itdrag,wspeed,cd,ach)
	    end if
            tx(k) = fice * wfact * cd * wspeed * wx(k)
            ty(k) = fice * wfact * cd * wspeed * wy(k)
          end do
	  wxymax = wslim
	end if

!	---------------------------------------------------------
!	check wind speed
!	---------------------------------------------------------

        if( wxymax .gt. wsmax ) then
          write(6,*) 'maximum wind speed: ',wxymax
          write(6,*) 'maximum allowed wind speed: ',wsmax
          write(6,*) 'Are you sure the wind is in the correct format?'
          write(6,*) 'If no, please set iwtype, else increase wsmax.'
          stop 'error stop meteo_convert_wind_data: wind speed too high'
        end if

!	---------------------------------------------------------
!	convert pressure
!	---------------------------------------------------------

	if( .not. has_pressure ) return

	call iff_get_var_description(id,3,string)

	if( string == ' ' .or. pfact == 0. ) then !only if not yet determined
	  pmin = minval(pp)
	  pmax = maxval(pp)
	  if( pmin /= 0 .and. pmax /= 0. ) then
	    if( pmin > 85000 .and. pmax < 110000 ) then
	      pfact = 1.
	      string = papa
	    else if( pmin > 850 .and. pmax < 1100 ) then
	      pfact = 100.
	      string = pamb
	    else
	      pfact = 1.
	      string = 'unknown'
	    end if
	    !call iff_set_var_description(id,3,string)
	  end if
	end if

	if( pfact /= 1. ) pp = pfact * pp

	pmin = minval(pp)
	pmax = maxval(pp)

	if( pmin /= 0 .or. pmax /= 0. ) then
	  if( pmin < 85000 .or. pmax > 110000 ) then
	    write(6,*) 'pmin,pmax: ',pmin,pmax
	    write(6,*) 'pressure values out of range'
	    write(6,*) 'pressure should be given in Pascal'
	    stop 'error stop meteo_convert_wind_data: pressure'
	  end if
	end if

!	---------------------------------------------------------
!	end of routine wind speed
!	---------------------------------------------------------

	end subroutine meteo_convert_wind_data

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine meteo_set_rain_data(id,nvar)

	integer id
	integer nvar

	character*60 string

        logical string_is_this_short
	real getpar

!	---------------------------------------------------------
!	check nvar and get parameters
!	---------------------------------------------------------

	if( nvar /= 1 ) then
	  write(6,*) 'no support for nvar = ',nvar
	  stop 'error stop meteo_set_rain_data: rain'
	end if

        rfact = nint(getpar('rfact'))
        ievap = nint(getpar('ievap'))
        irain = nint(getpar('irain'))

!	---------------------------------------------------------
!	handle rain
!	---------------------------------------------------------

        if ( iatm == 1 ) then                 !no file but coupling atm-oce
          if( ievap .eq. 1 ) then
            write(6,*) 'with atmosphere-ocean coupling'
            write(6,*) 'take precip/evap from atmosphere model'
            write(6,*) 'but ievap=1 (compute evap internally)'
            write(6,*) 'pleae set ievap equal to 0'
            write(6,*) 'to avoid considering evaporation twice'
            stop 'error stop meteo_set_rain_data: ievap'
          end if
        endif

	call iff_get_var_description(id,1,string)

	if( string == ' ' ) then	!TS file or constant
	  irtype = 0
	  if( iff_has_file(id) ) irtype = 1

	  if( irtype == 1 ) then
	    call iff_set_var_description(id,1,rain)
	  end if
	else
          if(string_is_this_short('rain',string)) then
	    irtype = 1
	  else
	    write(6,*) 'description string for rain not recognized: '
	    write(6,*) string
	    stop 'error stop meteo_set_rain_data: rain description'
	  end if
	end if

!	---------------------------------------------------------
!	remember values and write to monitor
!	---------------------------------------------------------

	if( irtype == 0 ) then
	  write(6,*) 'no rain file opened'
	else
	  write(6,*) 'rain file opened: ',irtype
	  call iff_get_var_description(id,1,string)
	  write(6,*) 'content: '
	  write(6,*) ' 1    ',string
	end if

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end subroutine meteo_set_rain_data

!*********************************************************************

	subroutine meteo_convert_rain_data(id,n,r)

! convert rain from mm/day to m/s

	integer id
	integer n
	real r(n)

	integer i
	real fact
        real, parameter :: zconv = 1.e-3 / 86400.	!convert mm/day to m/s

	fact = zconv * rfact				!rfact is extra factor

	do i=1,n
	  r(i) = r(i) * fact
	end do

	end subroutine meteo_convert_rain_data

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine meteo_set_ice_data(id,nvar)

	integer id
	integer nvar

	character*60 string

        logical string_is_this_short
	real getpar

!	---------------------------------------------------------
!	check nvar and get parameters
!	---------------------------------------------------------

	if( nvar /= 1 ) then
	  write(6,*) 'no support for nvar = ',nvar
	  stop 'error stop meteo_set_ice_data: ice'
	end if

!	---------------------------------------------------------
!	handle ice
!	---------------------------------------------------------

	call iff_get_var_description(id,1,string)
	string = adjustl(string)

	if( string == ' ' ) then	!TS file or constant
	  ictype = 0
	  if( iff_has_file(id) ) ictype = 1

	  if( ictype == 1 ) then
	    call iff_set_var_description(id,1,ice)
	  end if
	else
	  if(string_is_this_short('icecover',string)) then
	    ictype = 1
	  else
	    write(6,*) 'description string for ice not recognized: '
	    write(6,*) string
	    write(6,*) 'expecting: ',trim(ice)
	    stop 'error stop meteo_set_ice_data: ice description'
	  end if
	end if

!	---------------------------------------------------------
!	handle ice free areas
!	---------------------------------------------------------

	ia_icefree = -1		!this does not change ice cover

	if( ictype /= 0 ) then	!ice file has been opened
	  ia_icefree = nint(getpar('iaicef'))
	  if( ia_icefree == -99 ) then
	    write(6,*) 'ice file has been opened'
	    write(6,*) 'but parameter iaicef has not been set'
	    write(6,*) 'please set iaicef to the area code that'
	    write(6,*) 'indicates ice free conditions'
	    write(6,*) 'if no such areas exist please set iaicef=-1'
	    stop 'error stop meteo_set_ice_data: iaicef'
	  end if
	end if

!	---------------------------------------------------------
!	remember values and write to monitor
!	---------------------------------------------------------

	if( ictype == 0 ) then
	  write(6,*) 'no ice file opened'
	else
	  write(6,*) 'ice file opened: ',ictype
	  call iff_get_var_description(id,1,string)
	  write(6,*) 'content: '
	  write(6,*) ' 1    ',string
	end if

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end subroutine meteo_set_ice_data

!*********************************************************************

	subroutine meteo_convert_ice_data(id,n,r)

! convert ice data (delete ice in ice free areas, compute statistics)

	use evgeom
	use basin
	use shympi
	use mod_info_output

	integer id
	integer n
	real r(n)	!ice concentration

	integer k,ie,ii,ia,nflag,nice
	real rarea,rnodes
	real array(3)
	double precision dacu,dice,darea,area
	character*20 aline

	integer, save :: iuinfo = 0
	real, parameter :: flag = -999.
	logical, parameter :: ballcover = .false.  !cover whole area with ice

        if( iuinfo == 0 ) then
          iuinfo = -1
          if(shympi_is_master()) call getinfo(iuinfo)
        end if

	if( n /= nkn_local ) then
	  write(6,*) n,nkn_local
	  stop 'error stop meteo_convert_ice_data: n/=nkn_local'
	end if

	nflag = 0
	nice = 0

	do k=1,nkn_inner
	  if( r(k) == flag ) then
	    nflag = nflag + 1
	    r(k) = 0
	  else if( r(k) /= 0. ) then
	    nice = nice + 1
	  end if
	end do

	if( ballcover .and. nice > 0 ) r = 1.

	dice = 0.
	darea = 0.
	do ie=1,nel
	  ia = iarv(ie)
	  area = 4. * ev(10,ie)
	  if( ia == ia_icefree ) then
	    do ii=1,3
	      k = nen3v(ii,ie)
	      r(k) = 0.
	    end do
	  else
	    do ii=1,3
	      k = nen3v(ii,ie)
	      dice = dice + area*r(k)
	      darea = darea + area
	    end do
	  end if
	end do
	if( darea > 0. ) dice = dice / darea

	dacu = 0.
	darea = 0.
	do ie=1,nel
	  ia = iarv(ie)
	  area = 4. * ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    dacu = dacu + area*r(k)
	    darea = darea + area
	  end do
	end do
	dacu = dacu / darea

	rarea = dice
	rnodes = dacu
	!call get_act_timeline(aline)
	!if( iuinfo > 0 ) then
	!  write(iuinfo,*) 'ice: ',aline,rarea,rnodes,nflag
	!end if

	array = (/rarea,rnodes,real(nflag)/)
	call info_output('ice','sum',3,array,.true.)

	end subroutine meteo_convert_ice_data

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine meteo_set_heat_data(id,nvar)

	use shyfem_strings

	integer id
	integer nvar

	character*80 string,strings(4)
	integer i,ierr
	character*80 vapor,vshort

	logical string_is_this_short
        real getpar  

!	---------------------------------------------------------
!	check nvar and get parameters
!	---------------------------------------------------------

	if( nvar == 5 ) then
	  write(6,*) 'old format for heatflux not supported any more'
	  write(6,*) 'time series should have only 4 columns: '
	  write(6,*) 'sol_rad t_air humidity cloud_cover'
	  write(6,*) 'no support for nvar = ',nvar
	  stop 'error stop meteo_forcing_fem: heat'
	else if( nvar /= 4 ) then
	  write(6,*) 'no support for nvar = ',nvar
	  stop 'error stop meteo_set_heat_data: heat'
	end if

!	---------------------------------------------------------
!	handle heat
!	---------------------------------------------------------


	if ( iatm == 1 ) then			!no file but coupling atm-oce
          if( iheat .ne. 0 .and. iheat .ne. 7 ) then
            write(6,*) 'atmosphere-ocean coupling but the heat type is not available'
            write(6,*) 'set iheat = 0 or iheat = 7'
	    stop 'error stop meteo_set_heat_data: heat'
          end if
	end if

	do i=1,nvar
	  call iff_get_var_description(id,i,strings(i))
	end do
	call adjust_humidity_string(strings(3))		!FIXME

        ihtype = nint(getpar('ihtype'))
	if( ihtype == 1 ) then
	  vapor = rhum
	else if( ihtype == 2 ) then
	  vapor = wbtm
	else if( ihtype == 3 ) then
	  vapor = dewp
	else if( ihtype == 4 ) then
	  vapor = shum
	else
	  write(6,*) 'ihtype = ',ihtype
	  stop 'error stop meteo_set_heat_data: erroneous ihtype'
	end if

	call strings_get_short_name(vapor,vshort)

	if( strings(1) == ' ' ) then	!TS file or constant
	  ierr = 0
	  if( .not. iff_has_file(id) ) ihtype = 0	!no heat

          if( ihtype > 0 ) then
            call iff_set_var_description(id,1,srad)
            call iff_set_var_description(id,2,tair)
            call iff_set_var_description(id,3,vapor)
            call iff_set_var_description(id,4,ccov)
          end if
        else
	  ierr = 0
          if(.not.string_is_this_short('srad',strings(1))) then
            ierr = 1
          else if(.not.string_is_this_short('airt',strings(2))) then
            ierr = 2
	  else if(.not.string_is_this_short(vshort,strings(3))) then
            ierr = 3
          else if(.not.string_is_this_short('cc',strings(4))) then
            ierr = 4
          end if
        end if

	if( ierr /= 0 ) then
	  write(6,*) 'description string for heat not recognized: '
	  write(6,*) 'possible number of string with error: ',ierr
	  if( ierr == 3 ) then
	    write(6,*) 'ihtype = ',ihtype,'  ',trim(vapor)
	    write(6,*) 'description in file: ',trim(strings(3))
	  end if
	  do i=1,nvar
	    call iff_get_var_description(id,i,string)
	    write(6,*) i,'    ',string
	  end do
	  stop 'error stop meteo_set_heat_data: heat description'
	end if

!	---------------------------------------------------------
!	remember values and write to monitor
!	---------------------------------------------------------

	if( ihtype == 0 ) then
	  write(6,*) 'no heat file opened'
	else
	  write(6,*) 'heat file opened: ',ihtype
	  write(6,*) 'content: '
	  do i=1,nvar
	    call iff_get_var_description(id,i,string)
	    if( i == 3 ) call adjust_humidity_string(string)		!FIXME
	    write(6,*) i,'    ',trim(string)
	  end do
	end if

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end subroutine meteo_set_heat_data

!*********************************************************************
	subroutine adjust_humidity_string(string)
	implicit none
	character*(*) string
	if( string == 'humidity [%]' ) string = rhum
	end subroutine adjust_humidity_string
!*********************************************************************

        subroutine meteo_convert_heat_data(id,n,metaux,mettair,metcc,ppv,methum)

	integer id
	integer n
	real metaux(n)		!this is the vapor information read
	real mettair(n)
	real metcc(n)
	real ppv(n)
	real methum(n)		!return

	logical bnoheat

	bnoheat = ihtype == 0
	
        if( bnoheat ) then              !no heat
	  !nothing to be done
        else
	  call meteo_convert_temperature(n,mettair)	!kelvin -> C
	  call meteo_convert_cloudcover(n,metcc)	!% -> [0-1]
	  call meteo_convert_vapor(ihtype,n,metaux,mettair,ppv,methum) !rhum

	  !if( ihtype == 1 ) methum = metaux	!done in meteo_convert_vapor()
	  !if( ihtype == 2 ) metwbt = metaux
	  !if( ihtype == 3 ) metdew = metaux
	  !if( ihtype == 4 ) metshum = metaux
	end if

	end subroutine meteo_convert_heat_data

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine meteo_convert_cloudcover(n,cc)

! converts percent to fraction

	use mod_meteo
	use shympi

	implicit none

	integer n
	real cc(n)	!cloud cover

	real maxcc
	integer, save :: icall = 0
	integer, save :: ncall = 1
	integer, save :: ncall_max = 100

	if( any( cc > 1.5 ) ) then
	  maxcc = maxval(cc)
	  cc = cc / 100.
	  icall = icall + 1
	  if( mod(icall,ncall) == 0 ) then
	    ncall = ncall * 2
	    if( bmpi_master .and. ncall <= ncall_max ) then
	      write(6,*) 'cloudcover in %... converting to fraction',maxcc
	    end if
	  end if
	end if

	end subroutine meteo_convert_cloudcover

!*********************************************************************

	subroutine meteo_convert_temperature(n,tav)

! computes temperture to Celsius if in Kelvin

	use mod_meteo

	implicit none

	integer n
	real tav(n)	!dry air temperature

	where( tav > 200. ) tav = tav - tkelvin

	end subroutine meteo_convert_temperature

!*********************************************************************

	subroutine meteo_convert_vapor(mode,n,aux,tav,pav,rhv)

! computes relative humidity from other vapor values

	use mod_meteo

	implicit none

	integer mode
	integer n
	real aux(n)	!value read
	real tav(n)	!dry air temperature
	real pav(n)	!atmospheric pressure
	real rhv(n)	!relative humidity (return)

	integer i
	real ta,pp,rh,wb,dp,sh,val

	if( mode < 1 .or. mode > 4 ) then
	      write(6,*) 'mode = ',mode
	      stop 'error stop meteo_convert_hum: mode'
	end if

	do i=1,n
	  val = aux(i)
	  ta = tav(i)
	  pp = pav(i)/100.			!pressure in mbar
	  if( mode .eq. 1 ) then		!val is relative humidity
	      rh = val
	  else if( mode .eq. 2 ) then		!val is wet bulb
	    call wb2rh(ta,pp,val,rh)
	  else if( mode .eq. 3 ) then		!val is dew point
	    call dp2rh(ta,pp,val,rh)
	  else if( mode .eq. 4 ) then		!val is specific humidity
	    call sh2rh(ta,pp,val,rh)
	  end if
	  rhv(i) = rh
	end do

	end subroutine meteo_convert_vapor

!*********************************************************************
!*********************************************************************
!*********************************************************************

        subroutine get_drag(itdrag,wxy,dragco,ach)

! computes drag coefficient

        implicit none

        integer itdrag          !type of formula
        real wxy                !wind speed
        real dragco             !computed drag coefficient
	real ach		!Charnock parameter

        if( itdrag .le. 0 ) then
          !nothing
        else if( itdrag .eq. 1 ) then   !Smith and Banke (1975)
          dragco = 0.001 * (0.63 + 0.066*wxy)
        else if( itdrag .eq. 2 ) then   !Large and Pond (1981)
          if ( wxy .gt. 11. ) then
            dragco = 0.001 * (0.49 + 0.066*wxy)
          else
            dragco = 0.001 * 1.2
          end if
        else if( itdrag .eq. 5 ) then   !Hersbach
	  call hersbach_cd(ach,wxy,dragco)
        else
          write(6,*) 'erroneous value for itdrag = ',itdrag
          stop 'error stop get_drag: itdrag'
        end if

        end subroutine get_drag

!*********************************************************************

        subroutine hersbach_cd(ach,wxy,dragco)

! Computes drag coefficient with a formulation proposed in:
! Hersbach, H., 2011: Sea Surface Roughness and Drag Coefficient as 
! Functions of Neutral Wind Speed. J. Phys. Oceanogr., 41, 247–251, 
! https://doi.org/10.1175/2010JPO4567.1

	implicit none

	real, intent(in) :: ach,wxy
	real, intent(out) :: dragco

	double precision, parameter :: z = 10 !Height of the wind[m]
	double precision, parameter :: k = 0.4 !Von Karman constant
	double precision, parameter :: g = 9.8182 !gravity[m s**-2]
	double precision, parameter :: am = 0.11 !kinematic parameter
	double precision, parameter :: nu = 1.5*10.**(-5) !molecular viscosity[m**2 s**-1]
	double precision, parameter :: p = -12 !coefficient best representing the regime transition
	double precision :: bfit,bnu,bch
	double precision :: R,A

	if ( ach < 0.0001 .or. ach > 0.5 ) then
	  write(6,*) 'Strange value of the Charnock parameter = ',ach
	  write(6,*) 'Try values from 0.01 to 0.04.'
          stop 'error stop hersbach_cd: dragco'
	end if

	R = (z / (am * nu)) * (k * wxy)
	A = (ach / (g * z)) * (k * wxy)**2
	bnu = -1.47 + 0.93 * log(R)
	bch = 2.65 - 1.44 * log(A) - 0.015 * (log(A))**2

	bfit = (bnu**p + bch**p)**(1/p)

	dragco = (k / bfit)**2
	
	end subroutine hersbach_cd

!================================================================
        end module meteo_forcing_module
!================================================================

!*********************************************************************

	subroutine meteo_get_heat_values(k,qs,ta,rh,twb,uw,cc,p)

! returns meteo parameters for one node
!
! pressure is returned in [mb]

	use mod_meteo
	use meteo_forcing_module, only: iheat

	implicit none

        integer k                       !node number
        real qs                         !solar radiation [W/m**2]
        real ta                         !air temperature [Celsius]
        real rh                         !relative humidity [%, 0-100]
        real twb                        !wet bulb temperature
        real uw                         !wind speed [m/s]
        real cc                         !cloud cover [0-1]
        real p                          !atmospheric pressure [mbar, hPa]

	qs = metrad(k)
	ta = mettair(k)
	rh = methum(k)
	uw = metws(k)
	cc = metcc(k)

	if (iheat .ne. 7) then				  !skip, if we read fluxes
	  cc = max(0.,cc)
	  cc = min(1.,cc)
	  rh = max(0.,rh)
	  rh = min(100.,rh)

	  p = ppv(k)
	  p = 0.01 * p					  !Pascal to mb

	  call rh2wb(ta,p,rh,twb)
	endif

	end subroutine meteo_get_heat_values

!*********************************************************************

        subroutine meteo_get_heat_extra(k,dp,uuw,vvw)
 
! returns iextra meteo parameters for one node (iheat == 8)

        use mod_meteo

        implicit none

        integer k                       !node number
        real dp                         !dew point temperature [Celsius]  
        real uuw                        !u-component wind speed [m/s]
        real vvw                        !v-component wind speed [m/s]

        !dp  = metdew(k)  
	call rh2dp(mettair(k),ppv(k),methum(k),dp)
        uuw = wxv(k)
        vvw = wyv(k)

        end subroutine meteo_get_heat_extra

!*********************************************************************

	subroutine get_pe_values(k,r,e,eeff)

! returns precipitation and evaporation values
!
! eeff is evaporation used in model, if ievap==0 => eeff==0.

	use mod_meteo

	implicit none

	integer k
	real r			!rain [m/s]
	real e			!evaporation [m/s]
	real eeff		!effective evaporation [m/s] - used in model

	integer ievap
	save ievap
	data ievap /-1/

	real getpar

        if( ievap .eq. -1 ) ievap = nint(getpar('ievap'))

	r = metrain(k)
	e = evapv(k)
	eeff = e*ievap

	end subroutine get_pe_values

!*********************************************************************

	subroutine set_evap(k,e)

! sets evaporation

	use mod_meteo

	implicit none

	integer k
	real e			!evaporation [m/s]

	evapv(k) = e

	end subroutine set_evap

!*********************************************************************

	subroutine meteo_get_solar_radiation(k,qs)

	use mod_meteo

	implicit none

        integer k                       !node number
        real qs                         !solar radiation [W/m**2]

	qs = metrad(k)

	end subroutine meteo_get_solar_radiation

!*********************************************************************

        subroutine get_wind(k,wx,wy)

! helper function -> return wind for node k

	use mod_meteo

        implicit none

        integer k
        real wx,wy

        wx = wxv(k)
        wy = wyv(k)

        end subroutine get_wind

!*********************************************************************

        subroutine meteo_set_matrix(qs,ta,rh,uw,cc)

! interpolates files spatially - to be deleted

	use mod_meteo
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        real qs,ta,rh,wb,uw,cc

        integer k

        do k=1,nkn
          metrad(k) = qs
          mettair(k) = ta
          methum(k) = rh
          metws(k) = uw
          metcc(k) = cc
        end do

        end

!*********************************************************************

        subroutine meteo_get_ice_usage(uice)

! get parameter for ice usage

        use meteo_forcing_module

        implicit none

	real uice

	uice = amice

	end

!*********************************************************************

	subroutine meteo_has_ice_file(bice)

        use meteo_forcing_module

        implicit none

	logical bice

	bice = iff_has_file(idice)

	end

!*********************************************************************

	subroutine set_iheat(iheat_local)

        use meteo_forcing_module

        implicit none

	integer iheat_local

	iheat = iheat_local

	end

!*********************************************************************
