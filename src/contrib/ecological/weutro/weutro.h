
!--------------------------------------------------------------------------
!
!    Copyright (C) 2004  Georg Umgiesser
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

! weutro.h: header for weutro.f
!
! revision log :
!
! 24.08.2004    ggu     some new variables (idbox, wdebug etc)
!
!**************************************************************
!
!--------------------------------------------------------------
!--------------------------------------------------------------
!	EUTRO.CMN	
!--------------------------------------------------------------
!--------------------------------------------------------------

	integer segmax
	parameter (segmax=1)

	integer ISEG,graztype	!if graztype=1 simulate zoo if=0 use wasp form.

	integer idbox
	logical wdebug
	logical topseg
	logical botseg
	logical sedseg

	real daytime		!actual time [day]

        integer SYSBY(9)!bypass options for system: 0=simulated, 1=bypassed

	common/general/ISEG,graztype,idbox,SYSBY
	save/general/

	common/general_log/wdebug,topseg,botseg,sedseg
	save/general_log/

	common /eetime/daytime	!changed due to compiler warnings (ggu,13/3/01)
	save /eetime/

!--------------------------------------------------------------
!	state variables
!--------------------------------------------------------------

	real NH3,NO3,OPO4,PHYT,CBOD,DO,ON,OP	!8 state variables
	common /statev/ NH3,NO3,OPO4,PHYT,CBOD,DO,ON,OP
	save /statev/
	
	real CD(9,segmax)		!
	common /delta/ CD
	save /delta/

!--------------------------------------------------------------
!	geometrical and  physical variables set in main1.f
!--------------------------------------------------------------

	real VOL	!volume of box
	real H		!depth of box
	real SA		!surface of box
	real VEL	!velocity
	real STP	!temperature
	real SAL	!salinity
	real STP20	!STP-20

	common /geom/ VOL,H,SA,VEL,STP,SAL,STP20
	save /geom/

! depthg()	-> h
! salin()	-> sal
! bvol()	-> vol

!--------------------------------------------------------------
!	subroutine phyto
!--------------------------------------------------------------

!------ growth and respiration
!--------------variables-----

!	the following variables are calculated in phyto.f and used in the 
!	subroutines written in ()

	real GP1	!saturated growth rate of phytop (ammonia)
	real GPP 	!total growth of phytop (inorganop)
	real RESP	!phytop. respiration term (dissoxyg)
	real KMPHYT	!half saturation constant for phytop. (eutroint)
	real XEMPRC 	!auxiliar variable compute phytop. (organ -cn, -p)
	real GITMX1	!auxiliar var to compute growth value (smith) 
	real GITMAX	!auxiliar var to compute growth value (smith)

!-------------parameters---

        real K1RT	!respiration temperature coefficient unitless 
	real K1RC	!respiration rate constant endogenous term day-1
	real K1T 	!growth rate temperature coeff WHEN SYSBY(4)=1 
	real K1C	!maximum growth rate day-1 WHEN SYSBY(4)=1

	common /phytog/ GP1,GPP,RESP,KMPHYT,XEMPRC,K1C,GITMX1,K1T,GITMAX &
     &                ,K1RT,K1RC
	save /phytog/

!------ decomposition
!--------------variables----

        real DP1	!auxiliar variable to compute decomp value (eutrodmp)
	real DPP        !phytoplancton decomposition term (ammonia,cbodsv,
			!inorganp,organicn,organop)
	real K1D	!K1D=non predatory phytopl. death rate  day-1

!--------------parameters--

	real KPZDC	!decomp rate constant in sediment at 20 C, per day
	real KPZDT	!decomposition temperaure coefficient in sediment day-1

	common /phytod/ KPZDC,KPZDT,DP1,DPP,K1D
	save /phytod/


!--------nutrients limitation-
!------------variables-------

	real PNH3G1	!uptake of nitrogen (adimensional) ammonia preference
			!(ammonia,dissoxyg,nitrate)
	real XEMP1	!Michaelis Menten Limitations for  N (eutrodmp) 
	real XEMP2	!Michaelis Menten Limitations for  P (eutrodmp)
	real RNUTR	!auxiliar var to compute nutrient limitation (eutrodmp)

!	real CN		!nutrients concentration values in Michaelis formula N/L
!			!CN = NH3 + NO3 (local in pphyto)
!	real DOPO4	!nutrients concentr. in Michaelis formula P/L (local)

!-----------parameters----

	real KMNG1     	!nitrogen half saturn const for phyto growth microg N/L
	real KMPG1 	!P  half saturn const for phyt growth microg P/L
	real NUTLIM	!nutrient limitation option: 0=min, 1=multiple
!	real CCHL1  	!carbon to clorophyll ratio used only when LGTHSW=1
			!microg/microg
	real FOPO4	!spatially variable dissolved fraction of inorganic P	
!
	common /phytonp/ PNH3G1,KMNG1,KMPG1,XEMP1,XEMP2,FOPO4, &
     &                   RNUTR,NUTLIM
	save /phytonp/
!
!--------grazing------------
!	set in param.f if zooplankton is a state variable or a forcing:
!       if graztype=1 simulate
	
!       grazing from WASP
!	zooplankton: Time and segment variable herbivorous zooplankton populs
!	 are described as the product of the time variable population ZOO and 
!	segment specific ratios ZOOSG
!
!-----------parameters for forcing zooplankton (graztype not equal 1)--- ---

	real K1G	!K1G=grazing rate on Phytop. per unit Zoop.L/mgC-day
	real ZOO	!ZOO=time funct or const (default when set in param.f)
!			!zooplankton population mg C/L
	real ZOOSG(segmax)!ZOOSG=parameter  !fraction of herbivorous zooplankto

	common /phytograz/ K1G,ZOO,ZOOSG
	save /phytograz/
!
!-----variables and--parameters for variable zooplankton (graztype = 1)----
!-----variables------
	real GRZ	!grazing value
	real zoosk	!sink term of zooplankton, source term for 
			!organicn and organop
!-----parameters----
	real KGRZ	!grazing rate constant
	real KPHYZ	!half saturation constant for phytop in grazing
	real EFF	!grazing efficiency for zooplankton growth
	real KDZ	!death rate for zooplankton
	common /phytograzv/GRZ,KGRZ,KPHYZ,EFF,KDZ,zoosk 
	save /phytograzv/
!--------------------------------------------------------------
!	subroutine organop, inorganp (P=phosphorous)
!--------------------------------------------------------------

!-----------variables---

	real SK58	!sink term of organic P (inorganp,organop)
	real FPO4(segmax)!sediment flux P [mg/(m**2 day)]

!	real SR5P	! Source of organic P from phytop. resp (organop) 
!	real SR8P	!Source of organic P from phytop death (inorganp)
!	real SR8OP	!sink of inorganic P (inorganp) 
!	real SK8P	!sink term for inorganic P: algal uptake (inorganp)

!-----------parameters--
	real PCRB 	!P/C ratio mgP/mgC
	real FOP	!fraction of dead and respired phyto 
			!recycled to the organic P pool unitless
	real K58T	!temperature coeff:  mineralization  of dissolved  org P
	real K58C 	!rate for mineralization of dissolved  organic P
	real KOPDC	!rate for mineralization  of organic P in sediment
	real KOPDT	!temperature coeff for mineralization  of org P in sed
	real KPO4	!half saturation constant for P microgP/L

	common /com_phos/ PCRB,FOP,SK58, &
     &				K58C,K58T,KOPDC,KOPDT,KPO4,FPO4
	save /com_phos/

!--------------------------------------------------------------
!	subroutine organicn, ammonia,nitrate
!--------------------------------------------------------------

!-----------variables----

!	real SR13P	!Source of inorganic nitrogen from phyto(ammonia)
!	real SR13ON	!Source term of ammonia (ammonia)
!	real SK13P1	!sink term for inorganic N: algal uptake (ammonia)
!	real SR1413	!source for nitrate (nitrate)
!	real SK14P1	!sink term for nitretae: algal uptake (nitrate)
!	real SR10P	!Source of organic nitrogen from phytop (organicn)

	real SK1013	!Sink term of organic Nitrogen (ammonia, organicn)
	real SK1314	!nitrification term (ammonia,dissoxyg,nitrate)
	real SK14D	!denitrification term (cbodsv, nitrate)
	real FNH4(segmax)!sediment flux N [mg/(m**2 day)]

!-----------parameters--

	real NCRB	!N/C ratio
	real FON	!fraction of organic Nitrogen
	real K1320C	!rate for nitrification
	real K1320T	!temperature coefficient for nitrification
	real K140C	!rate for denitrification
	real K140T	!temperaure coefficient for denitrification
	real KNIT	!half saturation constant for nitrification
	real KNO3	!half saturation constant for and denitrificatiion
	real K1013C	!rate for mineralization dissolved organic nitrgen
	real K1013T	!temperature coefficient for min. dissolved organic N 
	real KONDC	!rate for mineralization of N in sediment
	real KONDT	!temperature coeff for mineralization of N in sediment

	common /comn/ NCRB,FON,SK1013, &
     &                 SK1314,SK14D,K1320C,K1320T, &
     &                K140C,K140T,KNIT,KNO3,K1013C,K1013T,KONDC, &
     &			KONDT &
     &			,FNH4
	save /comn/ 

!--------------------------------------------------------------
!	subroutine CBOD 	
!--------------------------------------------------------------

!------------variables----

!	real SR18P	!Source: phyto decomp.(if sedseg) or death (if water c.)
!	real SK18D	!Sink term: denitrification (cbodsv) 
!
	real DEATH	!Source term: phytop death (cbodsv,organicn,phyto)
	real SK180	!Sink term: oxidation (cbodsv,dissoxyg)

!----------parameters---

	real OCRB	!O/C ratio
	real KDC	!oxidation rate in water column


	real KDT	!temperature coefficient for oxidation in water column
	real KDSC	!oxidation rate in sediment at 20 C
	real KDST	!temperature coefficient for oxidation in sedim at 20 C
	real KBOD	!half saturation constant for oxidation

	common /comcbod/ OCRB,DEATH,SK180,KDC,KDT,KDSC,KDST,KBOD
	save /comcbod/ 

!--------------------------------------------------------------
!       subroutine dissoxyg
!--------------------------------------------------------------

!-------dissoxyg-physical-----
!
!------------variables---

!	real TK		!Temperature Kelvin, calculated internally (dissoxyg)
!	real RLNCS	!ln dissolved oxygen saturation value (dissoxyg)
!	real CL		!Cl value calc.from salinity value (dissoxyg,eutrodmp)

	integer ITYPE(segmax)	!segment types: (dissoxyg,eutrodmp,main1) 
			!1=surface water segment
			!2=subsurface water segment
			!3=upper bed segment
			!4=lower bed segment
	real KA		!Reareation rate (auxiliar variable) (dissoxyg,eutrodmp)
	real CS	!dissolved oxygen saturation valuie (dissoxyg,eutrodmp)

!-----------parameters---

	real WIND	!wind speed	
	real XICECVR	!time variable ice cover funct% default =1 no ice cover
	real AIRTMP	!air temperature
	real WTYPE	!water body size
			! 1 = small
			! 2 = medium
			! 3 = large

	real K2		!reareation rate 20 C (when not calc from physical vv) 

	common /rearphysic/ ITYPE,K2,WIND, &
     &          AIRTMP,WTYPE,KA, &
     &			CS, XICECVR
    	save /rearphysic/
	
!--------dissoxyg-biological--
!
!-----------variables---

!	real SR190	!Source: reareation term
!	real SR19PA	!Source: phytop. growth using NH3 and CO2
!	real SR19PB	!Source: phyto growth from N03 and CO2
!	real SK1913	!Sink term: nitrificatioin 
!	real SK1918	!Sink term: oxidation of CBOD
!	real SK19S	!Sink term: Sediment Ox demand

	real SR19P	!Source: phytoplankton growth (dissoxyg, eutrodmp)
	real SK19P	!Sink: algal respiration (dissoxyg, eutrodmp)

!-----------parameters--

	real SOD1D(segmax)!sediment oxygen demand for segment
	real SODTA(segmax)!segment specific temperature
			!correction coefficient (theta) for sediment ox. demand 
	
	common /rearphyto/ SR19P,SK19P,SOD1D,SODTA
	save /rearphyto/

!--------------------------------------------------------------
!	light limitation
!--------------------------------------------------------------

!-------------variables---

!	real GIT1	!growth rate (with light limitation)	
	real RLIGHT	!light limitation computed in smith and ditoro

	integer LGHTSW	!LIGHT SWITCH: 1=Di Toro 2=Smith  3=Steele
	real IS2	!optimum light intensity [??] (only for Steele)
	real IINST	!instantaneous light intensity

	common /lighti/ LGHTSW
	common /lightr/ RLIGHT,IS2,IINST
	save /lighti/,/lightr/

!--------------------------------------------------------------
!	subroutine ditoro
!--------------------------------------------------------------

! KESG(ISEG)	-> SKE
! KEFN(ISEG)	-> IKE
! KE(5)		-> SKE = SKE * KE(IKE)

!----------variables---

	INTEGER IKE
!	real SKE	!non algal light att.(SKE=KESG(ISEG)  (ditoro, smith)
	real KE(5)	!time var extin coeff funct unitless  (ditoro, smith)
	real IAVBOT	!auxiliar variable (used in  eutro?) (ditoro, smith)
	real RLGHTS(0:segmax,2)!RLGHTS(0:SG,2)controllare
			!in eutroint.f  RLGHTS (J, 1) = 0.0  RLGHTS (J, 2) = 1.0
			!(ditoro, smith)
	integer ITO	!ITO=IBOTSG(ISEG) segment immediately below ISEG 
			!(ditoro,main1,smith)
	real IAV	!average solar rad during the whole day (main1,smith) 
	real ITOTmp	!incident light intensity during the daylight
			!capire meglio !!!???!!!

!--------parameters---

	integer KEFN(segmax) !IKE=KEFN=pointer designating time variable
!	 ext. coeff. (KE) to be used for the segment in wasp KEFN=1,4.
	real KESG(segmax)!non algal light att.(SKE=KESG(ISEG)

	real FDAY	!fraction of day that is daylight unitless
	real IS1	!(Ly/day) saturation light intensity for phytop.used 
			! when LGHTSW=1 deifined originally in eutroint.f
	real CCHL	!C/Chla ratio

	common /comditoro/ KESG,KE,FDAY,IS1,IKE, &
     &                    IAVBOT,CCHL,IAV,RLGHTS,ITO,KEFN
	save /comditoro/

!--------------------------------------------------------------
!	subroutine smith
!--------------------------------------------------------------

!----------variables---

!	real I0		!time variable incident light intensity just below
	real IAVBOTX(segmax) !(smith)
			! the surface, assumed to follow a half sin function
			! over daylight hours, ly/day
!	real IMAX	!auxiliar variable for I0
!	real IAVSG	!auxiliar variable for RLIGHT
!	real SUM	!auxiliar variable for RLIGHT
	real IS1X(segmax) !=IS1

!--------parameters---

	real PI		!3.14
	real PHIMX	!mg carbon fixed per mole of light quanta absorbed
	real CCHLX(segmax)  !carbon to chlorophyll ratio.
		 !Used only when LGTHSW=1 (mg carbon/mb chl a).20-50 Default=30
	real XKC	!m2/mg chl a, phytop. self light attenuat(smith form)
	real ITOT	!incident solar radiation, ly/day
	real DTDAY	!??????verificare, serve a modulare il seno: day
	real NEWDAY	!unknown????if G.E. 1 calcolo smith (a new day)

	common /comsmith/ PHIMX,CCHLX,XKC,PI,ITOT,DTDAY, &
     &			NEWDAY,IS1X,IAVBOTX
	save /comsmith/

!--------------------------------------------------------------
!	subroutine kawind
!--------------------------------------------------------------

!	common /comwind/ DIFF,TW,VW,VA,TA,PA,PW,WS,RK,WH,SRCD,EF,F1, &
!     &		F2,FP1,FP2,FP3,FP4,N,SRCD2,ERR,US,Z0,RK1,GAMU,RK2, &
!     &		RK3,depth,ut,uc,ze,gam
!	save /comwind/

!--------------------------------------------------------------
!	subroutine kaydra
!--------------------------------------------------------------

!	real CFOREA 	!set in kahydra  (kahydra)      
!	real AVDEPE	!depth, H, (kahydra) 
!	real AVVELE	!abs(vel) water velocity (kahydra)
!	real REAK	!coefficient for reareation set in kahydra.f 
!	real EXPREV	!coefficient for reareation set in kahydra.f
!	real EXPRED	!coefficient for reareation set in kahydra.f
!	real TRANDP	!auxiliar variable in kahydra
!	real DIF	!auxiliar variable in kahydra

!---------variables----
	real K20	!auxiliar variable in the computation of krear(kahydra)
			!(kahydra, dissoxyg)

	common /comhydra/K20
	save /comhydra/

