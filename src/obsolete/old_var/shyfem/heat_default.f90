
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998,2001-2002,2006,2010-2011,2014  Georg Umgiesser
!    Copyright (C) 2018-2019  Georg Umgiesser
!    Copyright (C) 1998  Lucia Zampato
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

! heat flux module
!
! contents :
!
! subroutine subtem(dt,dh,qs,uw,cc,ta,tb,ts,tsnew,rtot)
!			compute water temperature in element of depth dh
! subroutine subtem0(dt,dh,qs,uw,cc,ta,tb,ts,tsnew,evap)
!			same as subtem, but use heatlucia
! subroutine heatlucia(t,p,w,tb,cc,ts,qsens,qlat,qlong,evap)
!			computes heat fluxes from Lucia modules
! subroutine longwave(ts,ta,tb,cc,rb)
!			long-wave radiation
! subroutine evcon(ts,ta,tb,uw,p,re,rc)
!			latent heat and convective heat
!
! revision log :
!
! 01.06.1998	ggu&lcz	written from scratch (nearly)
! 24.06.1998	ggu&lcz	subroutines from lucia integrated
! 30.04.2001	ggu	new routine qtotal_tb
! 09.12.2002	ggu	cleaned and re-arranged
! 23.03.2006	ggu	changed time step to real
! 23.03.2010	ggu	changed v6.1.1
! 08.10.2010	ggu	changed VERS_6_1_13
! 16.02.2011	ggu	pstd introduced
! 05.11.2014	ggu	changed VERS_7_0_5
! 24.01.2018	ggu	changed VERS_7_5_41
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
!
! notes :
!
! qs	net solar radiation (reflection already subtracted)
! ta	air temperature
! tb	wet bulb temperature
! uw	wind speed
! cc	cloud cover (0 clear sky, 1 totally covered)
!
!*****************************************************************************

      subroutine subtem(dt,dh,qs,uw,cc,ta,tb,ts,tsnew,evap)

!  Calcola la temperatura dell'acqua in un bacino di profondita' dh,
!  risolvendo l'equazione del bilancio termico all'interfaccia acqua-aria.
!
!  dt = timestep del modello FEM
!  dh = profondita' dello strato
!  qs = radiazione solare (short wave radiation) entrante nella superficie 
!       dell'acqua (radiazione incidente al netto della riflessione)
!  uw = velocita' del vento
!  cc = copertura nuvolosa (tra 0 e 1)
!  ts,ta,tb = temperatura dell'acqua, dell'aria, temp. dell'aria a bulbo
!             bagnato
!  tsnew = nuova temperatura dell'acqua
!
!  pa = pressione atmosferica in mbar in subroutine evcon 
!	(!!!!! verificare la sensibilita' a pa)
!
!  cpw = calore specifico dell'acqua di mare (J/kg C) da Gill, per ts=16 C
!  rhow = densita' dell'acqua di mare (kg/m3) da Gill, per ts=16 C
!
!  rtot = total radiation
!  qs = solar radiation (short wave)
!  rb = long wave radiation
!  re = evaporation term
!  rc = convection term 
!
!  Formulazione del bilancio termico da:
!  Mariotti M. e C. Dejak, 1982: 'Valutazione dei flussi energetici
!  all'interfaccia acqua-aria nella Laguna di Venezia', Ambiente Risorse 10.
!
!  Vengono usate le unita' SI.
!
!  L. Zampato - Dicembre 1997

      implicit none

      include 'subqfxm.h'

      real dt
      real dh
      real p
      real qs, uw, cc
      real ta, tb, ts, tsnew
      real evap

      real pstd
      parameter ( pstd = 1013.25 )

      real rb, re, rc, rtot
      real ct

! constants

      p = pstd

!  long-wave term
!   input: ts,ta,tb,cc
!   output: rb = long wave radiation

      call longwave(ts,ta,tb,cc,rb)

!  evaporation and convection terms
!   input: ts,ta,tb,uw
!   output: re,rc = evaporation and convection terms

      call evcon(ts,ta,tb,uw,p,re,rc)

!  total radiation: rtot (W/m2) (positive if into water)

      rtot = qs+rb+re+rc

!  evaporation [kg/(m**2 s)]

      evap = re / 2.5e+6                !divide by latent heat of evaporation
      evap = evap / rhow                 !in [m/s]
      evap = evap * 1000. * 86400.      !in [mm/day]

!  heat capacity/area

      ct = cpw*rhow*dh

!  new temperature      formula:  dQ = dT * rho * cpw * dh / dt

      tsnew = ts + rtot*dt/ct

      end

!*****************************************************************************

      subroutine subtem0(dt,dh,qs,uw,cc,ta,tb,ts,tsnew,evap)

! same as subtem, but use heatlucia

      implicit none

      include 'subqfxm.h'

      real dt
      real dh
      real qs, uw, cc
      real ta, tb, ts, tsnew
      real evap

      real pstd
      parameter ( pstd = 1013.25 )

      real rtot
      real qsens,qlat,qlong
      real ct, p

! constants

      ct = cpw*rhow*dh	!heat capacity	[ J / (m**2 K) ]
      p = pstd

      call heatlucia(ta,p,uw,tb,cc,ts,qsens,qlat,qlong,evap)

      rtot = qs - ( qlong + qlat + qsens )

!  new temperature      formula:  dQ = dT * rho * cpw * dh / dt

      tsnew = ts + rtot*dt/ct

      evap = -evap			!in [kg/(m**2 s)]
      evap = evap / rhow                 !in [m/s]
      evap = evap * 1000. * 86400.      !in [mm/day]

      end

!*****************************************************************************

        subroutine heatlucia(t,p,w,tb,cc,ts,qsens,qlat,qlong,evap)

! computes heat fluxes from Lucia modules
!
! heat fluxes are positive upward (from sea to atmosphere)

        implicit none

        real t          !air temperature [C]                    - in
        real p          !pressure [mb]                          - in
        real w          !wind speed [m/s]                       - in
        real tb         !wet bulb temperature [C]               - in
        real cc         !cloud cover [0-1]                      - in
        real ts         !sea temperature [C]                    - in
        real qsens      !sensible heat flux [W/m**2]            - out
        real qlat       !latent heat flux [W/m**2]              - out
        real qlong      !long wave radiation [W/m**2]           - out
        real evap       !evaporation [kg/(m**2 s)]              - out

	real ta,uw
	real rb,re,rc

	ta = t		!t air
	uw = w		!wind speed

        if( p > 10000 ) stop 'error stop heatlucia: p not in mbar'

	call longwave(ts,ta,tb,cc,rb)
	call evcon(ts,ta,tb,uw,p,re,rc)

	qlong = -rb
	qsens = -rc
	qlat  = -re
	evap  = -re / 2.5e+6

	end

!*****************************************************************************

      subroutine longwave(ts,ta,tb,cc,rb)

!  termine di radiazione long-wave nel bilancio termico atmosfera-mare
!
!  es, ea = emissivita' del mare e dell'aria
!  v = tensione di vapore alla temperatura t
!  ra, rb = radiazione emessa dall'atmosfera e dal mare
!  rb = termine long-wave in W/m2 da inserire nel bilancio termico
!  rbkj = termine long-wave in kJ/(m2 ora)
!
      implicit none
!
      real ts, ta, tb
      real cc
      real sigma, es, ea
      real alpha, beta
      real va, vd, psi, rr
      real ra, rs, rbkj, rb
!
      sigma=5.67/(10**(8))
      es=0.97
      ea=0.92
!     
      if (ta.ge.0.) then
         alpha=17.27
         beta=237.3
      else
         alpha=21.88
         beta=265.5
      endif
      va=6.11*exp(alpha*ta/(ta+beta))
      vd=6.11*exp(17.27*tb/(tb+237.3))
      psi=0.707 + vd/158.
      rr=1 + (0.25-0.005*(va-vd))*cc**2
!
      ra=ea*sigma*rr*psi*(ta+273)**4
      rs=es*sigma*(ts+273)**4
      rb=ra-rs
      rbkj=rb*3.6
!
      end

!*****************************************************************************

      subroutine evcon(ts,ta,tb,uw,p,re,rc)

!  termini di evaporazione-convezione nel bilancio termico atmosfera-mare
!
!  gam = costante psicrometrica in mbar/K
!  fu = funzione della velocita' del vento in W/(m2 mbar)
!  vpa = pressione parziale di vapore nell'aria  in mbar
!  re, rc = termini di evaporazione e convezione in W/m2
!  rekj, rckj = termine di evaporazione e convezione in kj/(m2 ora)

      implicit none

      real ts, ta, tb, uw, p, re, rc

      real pa, fu
      real gam
      real vd, vs, vpa
      real rekj, rckj

      pa = p
      gam = 0.66

      vd=6.11*exp(17.27*tb/(tb+237.3))
      vs=6.11*exp(17.27*ts/(ts+237.3))
      vpa= vd - (pa-vd)*(ta-tb)/(1540-1.3*tb)

      fu = 4.4 + 1.82*uw

      re = - fu*(vs-vpa)
      rekj = re*3.6

      rc = - fu*gam*(ts-ta)
      rckj = rc*3.6

      end

!*****************************************************************************

