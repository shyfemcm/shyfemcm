
!--------------------------------------------------------------------------
!
!    Copyright (C) 2009-2010,2014,2018-2019  Georg Umgiesser
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

! heat flux module (Gill)
!
! contents :
!
! subroutine heatgill(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)
!	computes heat fluxes from formulas in Gill
! subroutine tempgill(dt,dh,qsol,t,p,w,ur,cc,ts,tsnew,rtot,evap)
!	computes heat fluxes from formulas in Gill and adjusts water temp
!
! revision log :
!
! 27.08.2009	ggu	call to heatgill changed (pass ur, and not e,q)
! 23.03.2010	ggu	changed v6.1.1
! 05.11.2014	ggu	changed VERS_7_0_5
! 24.01.2018	ggu	changed VERS_7_5_41
! 16.02.2019	ggu	changed VERS_7_5_60
!
!*******************************************************************

	subroutine heatgill(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)

! computes heat fluxes from formulas in Gill
!
! heat fluxes are positive upward (from sea to atmosphere)

	use heat_const

	implicit none

	real t		!air temperature [C]			- in
	real p		!pressure [mb]				- in
	real w		!wind speed [m/s]			- in
	real ur		!relative humidity [%] ([0-100])	- in
	real cc		!cloud cover [0-1]			- in
	real ts		!sea temperature [C]			- in
	real qsens	!sensible heat flux [W/m**2]		- out
	real qlat	!latent heat flux [W/m**2]		- out
	real qlong	!long wave radiation [W/m**2]		- out
	real evap	!evaporation [kg/m**2/s]		- out

! other variables:
!
!	real q		!specific humidity [0-1]
!	real e		!vapor pressure [mb]
!	real r		!mixing ratio [0-1]

	real rdry,lv0,pascal
	parameter(rdry=287.04,lv0=2.5008e6,pascal=100.)
	real sigma,epsbbb,sigma0
	parameter(sigma=5.67e-8,epsbbb=0.985,sigma0=sigma*epsbbb)

	real rho,cp,lv,tv
	real e,r,q
	real es,rs,qs
	real dt,ch,ce
	real theta,theta2,theta4
	real cloud

        if( p > 10000 ) stop 'error stop heatgill: p not in mbar'

!	------------------------------------------------
!	initialization
!	------------------------------------------------

	call vapor(t,p,ur,e,r,q)	!compute e,r,q

	tv = ( t + kelv ) * ( 1. + 0.6078 * q )
	rho = pascal * p / ( rdry * tv )
	lv = lv0 - 2.3e3 * ts
	cp = 1004.6 * ( 1. + 0.8375 * q )

	call satur(ts,p,es,rs,qs)	!compute saturation values for sea

!	------------------------------------------------
!	sensible heat flux
!	------------------------------------------------

	dt = ts - t
	if( dt .lt. 0 ) then	!stable
	  ch = 0.83e-3
	else			!unstable
	  ch = 1.10e-3
	end if

	qsens = rho * cp * ch * w * dt

!	------------------------------------------------
!	latent heat flux
!	------------------------------------------------

	ce = 1.5e-3

	evap  = rho * ce * w * ( qs - q )
	qlat = lv * evap

!	------------------------------------------------
!	long wave radiation
!	------------------------------------------------

	theta = ts + 273.15
	theta2 = theta * theta
	theta4 = theta2 * theta2
	cloud = 1. - 0.6 * cc * cc

	qlong = sigma0 * theta4 * ( 0.39 - 0.05 * sqrt(e) ) * cloud

!	------------------------------------------------
!	end of routine
!	------------------------------------------------

        if( qlong .le. 0. ) then
          write(6,*) qsens,qlat,qlong
          write(6,*) ts,theta,theta4
          write(6,*) e,sqrt(e),0.05*sqrt(e),cc,cloud
          write(6,*) t,w,q,qs
          write(6,*) tv,lv,cp,rho
          stop 'error stop heatgill: qlong'
        end if

	end

!*******************************************************************

	subroutine tempgill(dt,dh,qsol,t,p,w,ur,cc,ts,tsnew,rtot,evap)

! computes heat fluxes from formulas in Gill and adjusts water temperature
!
! heat fluxes are positive upward (from sea to atmosphere)

	use heat_const

	implicit none

	real dt		!time step [s]				- in
	real dh		!depth of layer [m]			- in
	real qsol	!solar radiation [W/m**2]		- in
	real t		!temperature [C]			- in
	real p		!pressure [mb]				- in
	real w		!wind speed [m/s]			- in
	real ur		!relative humidity [%] ([0-100])	- in
	real cc		!cloud cover [0-1]			- in
	real ts		!sea temperature [C]			- in
	real tsnew	!new sea temperature [C]		- out
	real rtot	!total heat input [W/m**2]		- out
	real evap	!evaporation [kg/m**2/s]		- out

        real ct
        real qsens,qlat,qlong

!	------------------------------------------------
!	constants
!	------------------------------------------------

	ct   = cpw * rhow * dh	!heat capacity / area

!	------------------------------------------------
!	compute total radiation - positive if into water
!	------------------------------------------------

	call heatgill(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)
        write(67,*) qsol,-qsens,-qlat,-qlong

	rtot = qsol - ( qsens + qlat + qlong )

!	------------------------------------------------
!	compute new temperature: dQ = dT * rho * cpw * dh / dt
!	------------------------------------------------

	tsnew = ts + rtot * dt / ct

!	------------------------------------------------
!	end of routine
!	------------------------------------------------

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

