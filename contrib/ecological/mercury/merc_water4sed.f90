!--------------------------------------------------------------------------
!
!    Copyright (C) 2018  Ginevra Rosat
!    Copyright (C) 2018  Donata Melaku Canu
!    Copyright (C) 2019  Celia Laurent
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main
!    directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------
!
! sediment dynamics in water for mercury subroutines (sed4merc_water.f)
!
! revision log :
!
! 30.05.2018	dmc&gir	integration of the 0d module into SHYFEM
! 10.07.2019	clr	3d sinking, debug
!
! notes :
!
! dmc= donata canu; gir= ginevra rosati; clr= célia laurent
! laa= Leslie Aveytua Alcazar
!
!*****************************************************************

        subroutine sed4merc_water(bbottom,dtday,tday,wat_vol, &
     &                          wdepth,k,temp,sal,taub,area, &
     &                          C,Dssink,Dpsink,Vds,Vdp, &
     &                          ds_gm2s, dp_gm2s)


       implicit none

       integer nstate
       parameter (nstate=2)


      real C(nstate)         !sed variables: C(1):Silt, C(2):POM  [mg/L]
      real CD(nstate)        !derivatives CD(1) Dsilt, CD(2) DPOM [g/d]
      integer m              !time indicator for write outputs, debug
      integer k              !node
      logical bbottom        !is on bottom?
      real cold(nstate)      !old state variables
      real dtday,tday        !time step [day], time [day]
      real temp,sal          !temperature [C] and salinity [-]
      real tkel              !temperature [K]
      real tkref             !reference temperature, 293 [K]
      real loads(nstate)     !atmospheric loadings

      real Sw, POMw           !State variables, silt and POM [mg/L]
      real dsilt,sdens, spd   !diameter [m], particle density [g/cm3], particle density [kg/m3] for silt
      real dPOM, podens,ppd   !diameter [m], particle density [g/cm3], particle density [kg/m3] for POM
      real vis, swd           !viscosity [Pa s-1] and density [kg/m3] for seawater
      real g                  !acceleration gravity [m/sec2]

      real wdepth, area, wat_vol   !depth [m], area [m2] and volume [m3] of water elements
      real Vss,Vsp                 !Stoke's settling vel for silt and POM               [m/s]
      real taub                    !bottom stress from subssed.f [Pa]
      real tCDs                    !critical shear stress for deposition                 [Pa]
      real Swm,POMwm               !masses of silt and POM in water [g]
      real ter1                    !intermediate term for calculation [m/s]
      real Pd                      !Pd, probability of deposition [-]
      real Vds, Vdp                !Pd x Stoke's velocity for silt and POM              [m/s]
      real ds_gm2s, dp_gm2s        !deposition Flux for silt and POM [g/m2 s]
      real dep_gm2s                !total deposition flux (silt+POM) [g/m2 s]
      real Dssink,Dpsink,Dsink     !Sink of silt (Dssink), POM (Dpsink)and sum (Dsink) [g/s]

      integer ipext,ipint,kext     !nodes external and internal numbers
      integer fortfilenum, iter
      logical constant_parameters

!	variables
!       _______________________________________________________
!       assigne old value to variables

      Sw=C(1)
      POMw=C(2)        ![mg/l]

! --------------------------------------------------------------
!
      tkel=temp+273

!       ________________________________________________________

      constant_parameters=.False.

        if(constant_parameters)then  !claurent-OGS: enable the use of constant parameters for debug
          vis = 0.0012   ! din seawater viscosity (Pa s-1) or [kg/m-sec] 
          swd = 1030.0     ! Seawater density [kg m-3]! COMCelia: controllare che sed4merc_gas_exchange da stessa dimensione
          area=2000. 
          wdepth = 1.8 
          temp = 15.
          wat_vol = wdepth*area       ! m3

        else
          call sed4merc_gas_exchange(sal,temp,area,vis,swd)
        endif
 
!      write(*,*) 'RhoW', swd, 'Vis', vis, 'main'


! _______________________________________________________________
!     silt and POM particle properties
! _______________________________________________________________
!   ! move to INIT diameters FIXME

        dsilt   = 2./10.**5. ! silt diameter [m]
!       5*10^-5 = coarse silt, 6*10^-6 very fine silt, 1*10^-6 fine clay
        dPOM    = 5./10.**5. ! POM diameter [m] 
!       diatom cell 2*10-5 - 2*10^-4 um, picoplankton < 2*10^-6
        sdens   = 2.65       ! silt particle density    [g/cm3] 
        podens  = 1.25       ! POM particle density     [g/cm3]
        spd = sdens*10.**3.  ! silt particle density    [kg/m3] 
        ppd = podens*10.**3. ! POM particle density     [kg/m3] 
        g  = 9.81            ! acceleration gravity   [m sec-2]    
! _______________________________________________________________
!       Input critical shear for deposition and erosion
! _______________________________________________________________

        tCDs = 0.65 !0.7 (ORIG) !.08  !!da 0.06 a 1         !0.06*g*(spd-swd)*dsilt

! _______________________________________________________________
! Compute Stoke's settling velocities for silt and POM
! _______________________________________________________________

       ter1 = g/(18.*vis)                 ![m s-2]/[kg m-2 s-1]= [m s-1]
       Vss= ter1*(spd-swd)*(dsilt**2.)    ![m s-1]
       Vsp= ter1*(ppd-swd)*(dPOM**2.)     
! ______________________________________________________________
! _____ Deposition Occurrence 
! ______________________________________________________________   
       
       if (taub>1) then
         taub=1
       end if
 
       if (taub < tCDs) then            ! DEPOSITION
          Pd = (1. - taub/tCDs)          ! INVERTITI I  SEGNI
       else
          Pd = 0.
       end if

****** Compute Deposition Flux and rate *****************block
!     
       Vds = Pd*Vss               ![m s-1]
       Vdp = Pd*Vsp
       ds_gm2s = Vds*Sw           !Flux: [m s-1]*[g m-3]-> [g m-2 s-1]
       dp_gm2s = Vdp*POMw         !Flux
!       
       Dssink = ds_gm2s *area     !Sink of silt [g s-1] = [g m-2 s-1] * [m2] =Dsflux di Ginevra 
       Dpsink = dp_gm2s *area      !Sink of POM [g s-1] = Dpflux di Ginevra
       Dsink  = Dssink + Dpsink   ! = Dflux di Ginevra
!        
       dep_gm2s = Dsink/area
      
       kext=ipext(k)
 
!       if (kext .EQ. 2284) then
!       write(487,*) Vdp, Vds, 'sed4MERCw'
!       write(486,*) Dssink, Dpsink, 'sed4MERCw'
!       write (500,*) Sw
!       write (501,*) POMw
!       end if

       Swm   = Sw*wat_vol        ! [g m3]*m3 --> g
       POMwm = POMw *wat_vol  ! masses of soilds in water

        if (POMw .LE. 0.0) then  !if
        write(*,*),'instability - negative POMw in sed4merc_wat kext=', &
     &  ipext(k)
        stop
        else if (Sw .LE. 0.0) then 
        write(*,*),'instability - negative siltw in sed4merc_wat kext=', &
     &  ipext(k)
        stop
        end if
      
C       _________________________________________________________

        C(1)=Sw
        C(2)=POMw
!       __________________________________________________________

!         Dssink=Dssink    ! [g/sec] =Dsflux di Ginevra COMCelia: CANCELLED
!         Dsink=Dpsink     ! [g/sec] =Dflux di Ginevra COMCelia: CANCELLED

        CD(1) = -Dssink *86400   !g/day
        CD(2) = -Dpsink *86400   !g/day

!       call merc_euler (2,dtday,wat_vol,wat_vol,c,cold,cd)    ! c(i)=( c(i)*vol+dt*cd(i) )/vol
        call merc_euler (2,dtday,wat_vol,c,cold,cd)   


        if (POMw .LE. 0.0) then  !if
        write(*,*) 'POMw<=0',POMw,wdepth,ipext(k),'s4m_wat aft'
!        write(449,*) 'POMw<=0',POMw,wdepth,ipext(k),'s4m_wat aft'
!        stop
        else if (Sw .LE. 0.0) then
        write(*,*) 'siltw<=0',Sw,wdepth,ipext(k),'s4m_wat aft'
!        write(444,*) 'siltw<=0',Sw,wdepth,ipext(k),'s4m_wat aft'
!        stop
        end if

      iter=nint(tday*86400.)

      if (MOD (iter,1800) .EQ. 0) then
           kext=ipext(k)
           fortfilenum=-1
           if(kext==3985)then       !ve8 closer
               fortfilenum=250 
           elseif(kext==3986) then  !ve8  
               fortfilenum=251
           elseif(kext==3988) then  !ve8
               fortfilenum=252
           elseif(kext==4007) then  !ve8
               fortfilenum=253
           elseif(kext==3763) then  !ve7 
               fortfilenum=254
           elseif(kext==3765) then  !ve7
               fortfilenum=255
           elseif(kext==3764)then   !ve7
               fortfilenum=256
           elseif(kext==3762) then  !ve7 closer
               fortfilenum=257
           else if(kext==2150)then  !ve6 closer
               fortfilenum=258
           elseif(kext==2009) then  !ve6
               fortfilenum=259
           elseif(kext==2359) then  !ve6
               fortfilenum=260
           elseif(kext==2358) then  !ve6
               fortfilenum=261
           elseif(kext==2341) then  !ve5 closer
               fortfilenum=262
           elseif(kext==2408) then  !ve5
               fortfilenum=263
           elseif(kext==2191)then   !ve5
               fortfilenum=264
           elseif(kext==2192) then  !ve5
               fortfilenum=265
           elseif(kext==2654)then   !ve4 closer
               fortfilenum=266
           elseif(kext==2505) then  !ve4
               fortfilenum=267
           elseif(kext==2655) then  !ve4
               fortfilenum=268
           elseif(kext==2653) then  !ve4
               fortfilenum=269
           else if(kext==1372)then  !ve3 closer
               fortfilenum=270
           elseif(kext==1375) then  !ve3
               fortfilenum=271
           elseif(kext==1331) then  !ve3
               fortfilenum=272
           elseif(kext==1378) then  !ve3
               fortfilenum=273
           elseif(kext==4347) then  !ve3
               fortfilenum=247
             elseif(kext==3216) then  !ve2 closer
               fortfilenum=275
           elseif(kext==3057)then   !ve2
               fortfilenum=276
           elseif(kext==2953) then  !ve2
               fortfilenum=277 
           elseif(kext==3217) then  !ve2
               fortfilenum=278
           elseif(kext==2405) then  !ve1
               fortfilenum=279
           elseif(kext==2407)then   !ve1
               fortfilenum=280
           elseif(kext==2284) then  !ve1 closer
               fortfilenum=281
           elseif(kext==2404) then  !ve1
               fortfilenum=282
           endif
           if(fortfilenum.ge.0)then
!               if(fortfilenum==250)
!     +             write(*,*) 'stamp to file 250... at iter=',iter,
!     +             ', tday=', tday
               write(fortfilenum,"(2(i10,','),4(f15.7,','))")  &
     &         iter,kext,wdepth, Sw, POMw, taub
           endif
         endif


!      write(*,*) 'iteration', iter
!      write(*,*) 'ipint-3985' , ipint(3985), ipext(330)
!      write(*,*)
!      write(*,*) 'ipint-3763', ipint(3763), ipext(929)
!      write(*,*)
!      write(*,*) 'ipint-2009', ipint(2009), ipext(2135)

        end
!---------------------------------------------------------
!       set atmospheric loading on the surface layer
!--------------------------------------------------------
!********************************************************************

      subroutine load0ds(dt,cds,loads,vol)

! integrate loadings

      implicit none

        integer nstate          !total number of state parameters
        parameter( nstate =     2 )

      real cds(nstate)      !source term [g]
      real loads(nstate)      !loading for c [g/ day)on the element]
      real vol            !volume of box [m**3]
      real dt

      integer i

!        write(6,*) cds(1),cds(2),cds(3),dt,'cds before load'
      do i=1,nstate
        cds(i) = cds(i) +  loads(i)*dt
      end do

!        write(6,*) loads(1),loads(2),loads(3),dt,'loads'
!        write(6,*) cds(1),cds(2),cds(3),dt,'cds after load'

      end



!********************************************************************

	subroutine  sed4merc_gas_exchange(salin,temp,area,visc,rhow)
!
! 4.05.2017	dmc	

!	computing mercury exchange between air and the water column
!	Sorensen et al., 2010 An Improved Global Mercury Model
!	for Air Sea Exchange of Mercury: High Concentrations over the North
!	Atlantic


	implicit none
        logical bvis1		!select viscosity calculation bvis1 true
	save bvis1		!bvis true use model output else use Soerensen

!	general constant s

!	 parameters 
!        real AW		!Constant based on the Weibull distribution of 
!			!wind speeds over oceans
!        real mw		!molecular weight of water [g mol-1]\
!        real molHg	!molal volume of mercury at its normal boiling
                        !temperature [cm3 mol-1]
!        real phiw	!solvent association factor introduced to define 
			!the effective molecular weight of the solvent
			! with respect to the diffusion process 
	
!	auxiliar variables	
         real a,b, p,c,d
         real e,g,h,m,n,o
         real v1,v2,v3,v4,v5,v6


!	from the hydrodynamic model
     	real temp	!water temperature °C
     	real salin	!water salinity
!       	real uwind10	!wind speed normalised at 10 m above sea surface
     	real area	!element surface

!  	from mercury module
!      	real Hg0	!Hg0 concentration in water [g/ m3] or mg/L
!       	real Hg0atm	!Hg0 concentration in air   [g/ m3] or mg/L

!	variables and parameters calculated in this routine
	
!     	real mex	!Air_sea Exchange of mercury Hg0 at each time step
!			![kg*s-1]
!     	real mex2	!alternative formulation for test
!     	real flux2	!alternative formulation for test
!     	real flux	!Hg0 air-sea exchange flux [ng m-2 h-1]
!     	real Hlw	!dimensionless Henry's law constant
!     	real kw		!water side mass transfer coefficient for steady winds
!      	real kwb	!water side mass transfer coefficient Borgest et al.2004
!			!kwb Borges et al., 2004 formulation for microtidal syst
!      	real SchHg	!Schmidt number for mercury
!      	real ScCO2	!Schmidt number for CO2
!     	real kvis	!kinematic viscosity [cm2 s-1]
!     	real diff	!diffusivity (Wilke-Cang method) [cm s-1]
     	real visc	 ! water viscosity [cP]
        real rhow 	! density of the (sea)water  (KG/M**3)
      	real tempk	!temperature [K]
    	
!        mw=18.0		!molecular weight of water [g mol-1]
!        molHg= 12.74	!molal volume of mercury at its normal boiling
!     			!temperature [cm3 mol-1]
!        phiw=2.26	!solvent association factor 
!       	AW=0.25		!Weibull Constant based on wind distribution

!       from the hydrodynamic model

	tempk=temp+273.15	!temperature Kelvin

C ======================================================================
C ======================================================================
C
C Compute the density and the dynamic viscosity of water from the temperature
C and the salinity

C compute the dynamic/molecular viscosity
!      VISC0=1.802863d-3 - 6.1086d-5*TEMP + 1.31419d-06*TEMP**2 -
!       &1.35576d-08*TEMP**3 + 2.15123d-06*SALIN + 3.59406d-11*SALIN**2

        a=0.0001529
        b=0.000016826
        p=1.013253
        c=8.3885*(10**(-8))
        d=p**(2)
        e=0.0024727
        g=4.8429*(10**(-5))
        h=4.7172*(10**(-6))
        m=7.5986*(10**(-8))
        n=6.0574*(10**(-6))
        o= 2.676*(10**(-9))
        v1= temp*(0.06144-temp*(0.001451-temp*b))
        v2=a*p
        v3=c*d
        v4=e*salin
        v5=(n*p-o*d)*temp
        v6=((temp*g)-temp*(h-temp*m))*salin
 
        visc=(1.791- v1-v2+v3+v4+ v5+v6)/1000     ![kg m-1* s-1]

! mpute the water density according to Brydon et al. 1999, J. Geoph. Res.
C 104/C1, 1537-1540, equation 2 with Coefficient of Table 4, without pressure
C component. Ranges TEMP -2 - 40øC, S 0-42, surface water.
C      RHOW=9.20601d-2 + 5.10768d-2*TEMP + 8.05999d-1*SALIN
C     &     -7.40849d-3*TEMP**2 - 3.01036d-3*SALIN*TEMP +
C     %     3.32267d-5*TEMP**3 + 3.21931d-5*SALIN*TEMP**2
C      RHOW=RHOW+1000d0

C compute the water density according to EOS80, Fofonoff 198599,
C J. Geoph. Res. 90/C2, 3332-3342, without pressure component.
!	[kg * m-2]

      RHOW=999.842594d0 +6.793952d-2*TEMP -9.095290d-3*TEMP**2
     &   +1.00168d-4*TEMP**3 -1.120083d-6*TEMP**4 +6.536332d-9*TEMP**5
     & +(8.24493d-1 -4.0899d-3*TEMP +7.6438d-5*TEMP**2
     &   -8.2467d-7*TEMP**3 +5.3875d-9*TEMP**4) * SALIN
     & +(-5.72466d-3 +1.0227d-4*TEMP -1.6546d-6*TEMP**2) * SALIN**1.5d0
     & +4.8314d-4*SALIN**2

!	bvis1=.true.
!	if(bvis1)then
!
!	kvis=(VISC/RHOW)*10000	![cm2 s-1]
!	else
!
!	kvis=0.017*exp(-0.025*tempk)
!	
!	end if
!       write(*,*) 'RhoW', RHOW, 'Vis', visc, 'subroutine'

C ======================================================================
C ======================================================================

!        diff=((7.4*0.00000001)*((phiw*mw)**0.5)*tempk)
!        diff=diff/(visc*1000*(molHg**0.6))      !diffusivity

!	Hlw=exp((-2403.3/tempk)+6.92)
!	ScCO2=0.11*temp**2-6.16*temp+644.7	!Soerensen da Poissant et al 2000
!c       ScCO2=2073.1+125.62*temp+3.6276*temp*temp-0.043219*temp**3 !Wanninkhof salt
! c      ScCO2=1911.1+118.11*temp+3.4527*temp*temp-0.04132*temp**3 !Wanninkhof fresh
!        SchHg=kvis/diff
!        kw=Aw/100*24*uwind10**2*(SchHg/ScCO2)**(-0.5)  ![m day-1]
!        flux=(kw)*(Hg0-Hg0atm/Hlw)       !g m-2 day-1] FIXMME
!c        !mex=area*flux/1000    ![kg s-1]

!	alternative equation of 

!        kwb=(0.1+2.26*uwind10*(SchHg/ScCO2)**(-0.5))/100*24  ! [m day-1] Borges et al., 2004
!       flux2=(kwb)*(Hg0-Hg0atm/Hlw)        ![g m-2 day-1]
!       !mex2=area*flux2/1000            ![kg day-1]


!c        write(97,*)visc,RHOW,diff,temp
!c        write(98,*)ScCO2,SchHg, Hlw
!c	write(82,*)kw,kwb
!c	write(77,*)uwind10,Hg0,Hg0atm,flux,flux2 


        !write(6,*) RHOW,visc,temp,salin, 'rhow,visc in gas_exchange'
	end			!end of routine






