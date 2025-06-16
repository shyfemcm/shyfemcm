
!--------------------------------------------------------------------------
!
!    Copyright (C) 2011-2012,2014-2015,2018-2019  Georg Umgiesser
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

! $Id$
!
! meteo files
!
! revision log :
!
! 16.02.2011	ggu	created by copying mainly from subn11.f
! 14.07.2011	ggu	changed VERS_6_1_27
! 17.02.2012	ggu	changed VERS_6_1_45
! 24.02.2012	ggu	changed VERS_6_1_46
! 30.03.2012	ggu	changed VERS_6_1_51
! 30.05.2014	ggu	new imreg == 3
! 10.07.2014	ggu	only new file format allowed
! 18.07.2014	ggu	changed VERS_7_0_1
! 30.10.2014	ggu	changed VERS_7_0_4
! 05.11.2014	ggu	changed VERS_7_0_5
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 09.01.2015	ggu	changed VERS_7_0_12
! 19.01.2015	ggu	changed VERS_7_1_3
! 26.02.2015	ggu	changed VERS_7_1_5
! 30.04.2015	ggu	ice cover implemented
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 23.09.2015	ggu	changed VERS_7_2_4
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 16.02.2019	ggu	get_ice() -> get_ice_cover(), new set_ice_cover()
! 02.05.2022	ggu	in convert_distributed() use ie_mpi
! 13.10.2024	ggu	in convert_distributed() exchange rqdsv
! 03.12.2024    lrp     irain for the coupled model
!
!*********************************************************************

	subroutine convert_distributed

! converts distributed source from [m/s] to [m**3/s]

	use mod_bound_dynamic
	use evgeom
	use basin
	use shympi

	implicit none

	integer k,ie,ii,ie_mpi
	real area3
	real v1v(nkn)

	do k=1,nkn
	  v1v(k) = 0.
	end do

	do ie_mpi=1,nel
	  ie = ip_sort_elem(ie_mpi)
	  area3 = 4. * ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    v1v(k) = v1v(k) + area3
	  end do
	end do

	do k=1,nkn
	  rqdsv(k) = rqdsv(k) * v1v(k)
	end do

	call shympi_exchange_2d_node(rqdsv)

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine evap_init

! initializes evaporation mass flux

	use mod_meteo
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k

	do k=1,nkn
	  evapv(k) = 0.	
	end do

	end

!*******************************************************************

	subroutine rain_evap_set

! adds evaporation mass flux to distributed source

	use mod_meteo
	use mod_bound_dynamic
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,ievap,irain
        double precision racum,eacum
	real getpar

	ievap = nint(getpar('ievap'))
        irain = nint(getpar('irain'))

	do k=1,nkn
	  rqdsv(k) = rqdsv(k) + irain*metrain(k) - ievap*evapv(k)
	end do

        !racum = sum(metrain) / nkn
        !eacum = sum(evapv) / nkn
        !write(345,*) racum,eacum

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine windcd_init

! initializes evaporation mass flux

	use mod_meteo
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k
	real dragco

	real getpar

	dragco = getpar('dragco')

	do k=1,nkn
	  windcd(k) = dragco
	end do

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine meteo_init

! initializes meteo variables

	use meteo_forcing_module

	implicit none

	write(6,*) 'initializing meteo'

	call evap_init
	call windcd_init

	call meteo_forcing_fem	!new file format

	end

!*******************************************************************

	subroutine meteo_force

! update meteo variables and admin rain/evaporation

	use meteo_forcing_module

	implicit none

	call meteo_forcing_fem	!new file format

	!call compute_heat_flux

	call rain_evap_set		!add evaporation
	call convert_distributed	!convert from [m/s] to [m**3/s]

	end

!*******************************************************************

	subroutine compute_heat_flux

! computes heat flux through bulk formulas

	use mod_ts
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	double precision dq,dtime
	real dt

        call get_timestep(dt)
	call get_act_dtime(dtime)
        call qflux3d(dtime,dt,nkn,nlvdi,tempv,dq)	!compute heat flux

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine get_meteo_forcing(k,wx,wy,tauxn,tauyn,p)

! returns wind (wx/y), normalized stress (taux/yn) and pressure (p)

	use mod_meteo

	implicit none

	integer k		!node number
	real wx,wy		!wind velocity [m/s]
	real tauxn,tauyn	!normalized stress [m**2/s**2]
	real p			!pressure [Pa]

	wx = wxv(k)
	wy = wyv(k)
	tauxn = tauxnv(k)
	tauyn = tauynv(k)
	p = ppv(k)

	end

!*******************************************************************

	subroutine get_light(k,rad_light)

! returns light intensity [W/m**2]

	implicit none

        integer k               !node number
        real rad_light          !watt/m**2

	call meteo_get_solar_radiation(k,rad_light)

        end

!*******************************************************************

	subroutine set_ice_cover(k,ice_cover)

! sets ice cover [fraction 0-1]

	use mod_meteo

	implicit none

        integer k               !node number
        real ice_cover          ![0-1]

	metice(k) = ice_cover

	end

!*******************************************************************

	subroutine get_ice_cover(k,ice_cover)

! returns ice cover [fraction 0-1]

	use mod_meteo

	implicit none

        integer k               !node number
        real ice_cover          ![0-1]

	ice_cover = metice(k)

        end

!*******************************************************************

	subroutine get_ice_cover_all(ice_cover)

! returns ice cover array [fraction 0-1]

	use mod_meteo
	use basin, only : nkn,nel,ngr,mbw

	implicit none

        real ice_cover(nkn)          ![0-1]

	integer k

	do k=1,nkn
	  ice_cover(k) = metice(k)
	end do

        end

!*******************************************************************

