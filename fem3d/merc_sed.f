
!--------------------------------------------------------------------------
!
!    Copyright (C) 2019  Georg Umgiesser
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

! revision log :
!
! 20.12.2019	ggu	changed VERS_7_5_67

         subroutine mercury_sed_react(dtday, &
     &                          k,temp,area,C,Cw, &
     &             Shgsil, Shgpom, Smhgsil, Smhgpom, &
     &             fdiss1w,fdiss2w,fdoc1w,fdoc2w, &
     &             silt,pom,Vr,bvels,bvelp)

!	
************************************************************************
* 
* !! reaction mass yied coefficients: 1.07 (met); 0.93 (demet) !!
* !! oppure calcolare in mol e poi convertire (?) !!

! -----------------------------------------------------------------------
! ---- Call from bentic sediment model: ---------------------------------
! 	 Vr, Bvels, Bvelp
!      silt, POM, por    
!-----------------------------------------------------------------------
!
! INITIAL CONDITIONS: concentrations of hg and mehg [ mg(hg)/kg(sed) ]
! that are converted to hgit and mehgt [ ug(hg)/m3(s+w) ] 
!
************************************************************************
      implicit none
!	include 'mercury.h'
      integer m,nvmerc !water or sediment box
      parameter (nvmerc=2)

      real, intent(IN)    :: dtday
      real, intent(INOUT) :: C(nvmerc)       !mercury variables C(1):hgit, C(2):mehgt     [ug/m3]
      real, intent(IN)    :: Cw(3)           !mercury in water var : Hg0, Hg2, MeHg       [ug/m3]
      real, intent(IN)    :: temp            !temperature (^C) 
      real, intent(IN)    :: area            ![m2] 
      integer,intent(IN)  :: k               !nodes numeber (internal)       [-]
      real, intent(IN)    :: fdiss1w, fdoc1w !Hg fraction diss and in DOC [-]    
      real, intent(IN)    :: fdiss2w, fdoc2w !MeHg fraction diss and in DOC [-]
      real, intent(IN)    :: Shgsil,Shgpom,Smhgpom,Smhgsil    ![ug/d]
      real, intent(IN)    :: silt, pom       ![mg/l]
      real, intent(IN)    :: Vr,bvels,bvelp ![m/s]

      real                :: CD(nvmerc)     !derivatives CD(1), CD(2)
      real                :: CDw(3)          !mercury in water derivative: Hg0, Hg2, MeHg [ug/s]
 
      real depth, vol
      real dt    !time step, day
      real cold(nvmerc)
      real coldw(3)

!------------------------------------------------------------------
!----- variables - initial COND. - ug(hg)/g(sed)  
!	------------------------------------------------------------------
      real Hg2sed, MeHgsed   ![ug(hg)/g(sed)]
      real hgit, hgp, hgd,hgp1, hgp2     ! [ug m-3] or [ng(hg)/l(w+s)] 
      real mehgt, mehgd_ngl, mehgp1, mehgp2, mehgpi
!
      ! hgp = Hg2sed*Csilt/10**3 !   
      ! kd/por =  Hg2sed/ hgd 
      ! hgit = hgp + hgd	        [ug m-3]or [ng(hg) l(w+s)-1] 
! ----------------------------------------------------------------------
!     ! Hg2x  = hgit   * fx1
!     ! MeHgx = mehgt * fx1     [ug m-3]or [ng(hg) l(w+s)-1] 
! ----------------------------------------------------------------------
      real Hg2Cl, Hg2DOC, Hgpsilt, Hg2POM, Hg2P      ! [ug m-3]or [ng(hg) l(w+s)-1] 
      real MeHgCl, MeHgDOC,MeHgPOM, MeHgsilt, MeHgP  ! [ug m-3]or [ng(hg) l(w+s)-1] 
      real Hg2D, MeHgD ![ng l(w+s)-1]  
      real Hg2dpw, MeHgdpw, Hg2pw, MeHgpw ![ng l(w)-1] hg species in pore-water

      real hgitw, mehgtw,HgDw,MeHgDw   !hg species in water [ng(hg)l-1] 

! 	Hg2dpw, MeHgdpw sono usate x calcolo fllussi diff; Hg2pw, MeHgpw calcolate dopo integrazione
!	real Hg0sed    ! ? ---- reductive demethylation produces hg0 ---- 
!
!----- partition coefficients and fractions *****************************
!	------------------------------------------------------------------
      real K1silt, K1doc,K1POM,K1tp, num1  !Kd of Hg2 to silt,doc !K1sand, 
      real K2silt, K2doc, K2POM, K2tp,num2 !Kd of MeHg to silt, sand, doc, poc
      real partden1, partden2  !auxiliar variable in computation
      real d1,si1,sa1,P1, d2,si2,sa2,P2
      real faq1, fsilt1, fsand1, fPOM1, fdoc1, ftot1 ! fraction of hg sp. in sed
      real faq2, fsilt2, fsand2, fPOM2, fdoc2, ftot2 ! fraction of mehg sp. in sed
! --------------------------------------------------------------------
!----- Hg and MeHg resuspension and burial -----------------------------
!   --------------------------------------------------------------------
      real  DOC                 ! [g  m-3]
      real  Cdoc, Csilt,  Cpom, Ctot       ! [g cm-3] or [kg l-1]
      real  Rhgsil, Rhgpom, Rmhgsil, Rmhgpom !resuspension vel [m s-1] and fluxes [..]         
      real  bursHg, burpHg, bursMHg, burpMHg !burial vel [m s-1] and fluxes
! ----------------------------------------------------------------------
!----- Pore-water diffusion --------------------------------------------
!   --------------------------------------------------------------------	
      real pw_m3,pw_L, por, tor !pore-w volume,liters, porosity and tortuosity 
      real hgd_ugl, hgd_ngl     !conversions for diss hg in pore-water
!	
      real Dchg_t, Dchg_25 !diffusion coefficient [m2 s-1]
      real ft1, JHgD,   num3, num4, Jngm2d !difffusion flux
      real ft2, JMHgD,  num5, JMngm2d !difffusion flux 	
! --------------------------------------------------------------------
!----- Hg methylation and MeHg demethylation ---------------------------
!   --------------------------------------------------------------------    
      real ksmet, ksdem, Qbac
      real Eadem, tkel, tkref, R, Rcal, deltat
      real cksmet,sksme, cksdem, sksdem
      real add1,add2
!
!       deposition terms
!
      integer fortfilenum
      integer ipext,ipint,kext     !nodes external and internal

!	------------------------------------------------------
!------  BOX features --------------------------------------------------
!	------------------------------------------------------
      depth = 0.02             ! m
      vol= depth*area          ! m3
      
      DOC = 15.                  ! [mg/l] di DOC in sediment
      por = 0.7                  ! FIXME read from sed4merc_sed GR 11-02-2020     
!	-----------------------------------------------------------------------------
!	---------REACTIONS--rate constants-------------------------------------------

      ksmet = 0.09      ! d-1
      ksdem = 0.12      ! d-1
      Qbac = 1.5

!	------------------------------------------------------	
! -- partition coefficients [L/kg] for mercury into silt,sand,DOC,POC
        K1silt= 10.**(5)        !Hg2 in silt   
        K1doc=  1.0*10.**5      !Hg2 in doc
        K1POM=  10.**(5)        ! Hg2 to POM particles 
  
        K2silt= 1.0*10.**4.     !MeHg in silt  
        K2doc=  2.0*10.**4     !MeHg  in doc
        K2POM=  10.**5         !MeHg in POC particles 

! ---------------Assign old variables-----------------------

        hgitw=Cw(2)
        mehgtw=Cw(3)
        hgit =C(1)
        mehgt=C(2)

! -----------------------------------------------------------
! ----- Lamborg et al., 2016 --------------------------------
! ------- log(kdPOM) = 6.72+-2%; kdCacO3 = log(6.71)+-4% ----
! -----------------------------------------------------------
! --------log(kd-lithogenic part) =5.84+-9% -----------------
! -----------------------------------------------------------
! --- Soerensen Baltic sea: one Kd proportional to LOI ------
! ---- LOI = 4.2854*p_OC+0.859-------------------------------
! -----kd hg = 2.97 + 0.15*LOI ---(= 10**3 -10**4)----------- 
! -----kd mehg = 1.98 + 0.18*LOI ---(= 10**2.2 -10**3.1)-----

***** call from solids water module********************************
!	
!	Sw, POMw, Vds, Vdp   ! solids water and sink velocities
!--------------------------------------------------------------
****** call from Hg water module***********************************
! 
!     hgitw, MeHgw ! hgit in water [ng l-1] or [μg m-3]    
!     fsilt1w,fpom1w, fdiss1w  !fraction Hg-silt -pom -diss water
!      fsilt2w, fpom2w, fdiss2w  ! fraction MeHg-silt -pom -diss water
!
****** call from solids module***********************************
****** CONCENTRATIONS of silt, POM, DOC in sed*********part 04**********
******* converted from [g m-3] to [g cm-3] or [kg l-1]***********
 
      Cdoc  = por*DOC/10.**6. ! g(DOC)/m3(w) * 0.8 m3(w)/m3(s+w) -> mg(DOC)/L(s+w)/10**6 -> kgl-1
      Csilt = silt/10.**6.  ! from [g m-3] to [g cm-3] or [kg l-1]      
      Cpom  = POM /10.**6.  ! [kg l-1]   
      Ctot  = Csilt+Cpom      
 
      num1 = (K1silt *Csilt) + (K1POM*Cpom) !+(K1DOC*Cdoc) ! L/kg* kg/L
      num2 = (K2silt *Csilt) + (K2POM*Cpom) !+(K2DOC*Cdoc)	
      K1tp = num1/Ctot!  KD to all solids
      K2tp = num2/Ctot

!      write(*,*) 'K1tp',K1tp,'K2tp',K2tp 

      pw_m3 = por*vol    ! [m-3(w)] 
      pw_L  = pw_m3*1000.               ! g m-3 * m3 -> g of pore water
 
!------------------------------------------------------------------------
! ---------- Hg and MeHg Fractions ------------- ------------------------
!------------------------------------------------------part 06-----------
! ---------- kD  [g(hg)/kg(s)] / [g(hg)/l(w)] --> [l(w)/kg(s)] ----------
! ---------- por [l(w)/l(s+w)] ------------------------------------------ 
! ---------- kD/por = [l(w)/kg(s)]*[l(s+w)/l(w)]*Csolid [kg(s) l-1(s+w)]- 
!------------------------------------------------------------------------

      d1  = (K1doc/por) * Cdoc  !   
      si1 = (K1silt/por)* Csilt
      P1  = (K1POM/por) * Cpom   ! kd'*m
           
      d2  = (K2doc/por)  * Cdoc    
      si2 = (K2silt/por) * Csilt   
      P2  = (K2POM/por)  * Cpom    
      
      partden1 = 1. + si1 + P1 + d1  ! partden1=1.0+(d1+si1+sa1+P1)   
      partden2 = 1. + si2 + P2 + d2  ! partden2 =1+(si1+sa2+ P2 + d2)	
 
      faq1= 1. /partden1     !fraction of freely dissolved  Hg2	
      faq2 = 1./partden2     !fraction of freely dissolved  MeHg
 
      fsilt1 = si1/partden1  !(0.000001*k1silt*silt)/partden1 !silt sorbed Hg2
      fsilt2 = si2/partden2  !silt sorbed MeHg

      fPOM1 = P1/partden1   !(0.000001*K1POM*POC)/partden1	!POC matter sorbed Hg2
      fPOM2 = P2/partden2   !POC matter sorbed MeHg

      fdoc1 = d1/partden1
      fdoc2 = d2/partden2 !doc sorbed MeHg

      ftot1 = faq1 + fsilt1 + fPOM1 +fdoc1    !fsand1	
      ftot2 = faq2 + fsilt2 + fPOM2 +fdoc2     !+fsand2+

!      write(*,*) ':::::::::::::: PARTITIONING :::::::::::::::::::::::::'
!      write(*,*) '    ' 
!      write(*,*) 'K1doc :',K1doc , 'K2doc :', K2doc	!---------- 
!      write(*,*) 'K1silt', K1silt, 'K2silt:', K2silt	!------------
!      write(*,*) 'K1POM',  K1POM,  'K2POM:', K2POM	!------------
!      write(*,*) '    '           
!      write(*,*) 'por', por, 'Doc', DOC	!-----------------------
!      write(*,*) 'partden1', partden1, 'partden2', partden2	   
!      write(*,*) 'fsilt1', fsilt1, 'fsilt2', fsilt2	!------------- 
!      write(*,*) 'CDOC', Cdoc, 'Cpom', Cpom	!------------------- 
!      write(*,*) 'fPOM1', fPOM1, 'fPOM2', fPOM2	!------------------- 
!      write(*,*) 'fdoc1', fdoc1, 'fdoc2', fdoc2	!-------------------       
!      write(*,*) 'ftot HgSed', ftot1, 'ftot MeHgSed', ftot2    !----------------                
!      write(*,*) ':::::::::::::::::::::::::::::::::::::::::::::::::::::'
!      write(*,*) ' '       
!------------------------------------------------------------------------
!------------------------------------------------------------------------
! ---------- Benthic dynamics   -----------------------------------------   
!------------------------------------------------------part 07-----------
! ----------Hg MeHg RESUSPENSION  --------------------------------------- 
!   

      kext=ipext(k)
 
      Rhgsil = Vr * hgit * fsilt1 *area  !m s-1 * ug m-3 * m2-> ug s-1  
      Rhgpom = Vr * hgit * fPOM1 *area   !m s-1 * ug m-3 * m2-> ug s-1  
!        
      Rmhgsil = Vr * mehgt * fsilt2 *area  ! m s-1 * ug m-3 -> ug s-1  
      Rmhgpom = Vr * mehgt * fPOM2 *area   ! m s-1 * ug m-3 -> ug s-1   

!      
!      if (hgitw .GT. 50) then
!      write(777,*) hgitw,Vr,hgit, fsilt1,ipext(k)
!      write(778,*) Rhgsil, Rhgpom, Rmhgsil, Rmhgpom,ipext(k)  
!      end if


!      write(771,*) Vr, Rhgsil/area, Rmhgsil/area,ipext(k) 
!      write(772,*) Vr, Rhgpom/area, Rmhgpom/area,ipext(k)


!      write(*,*) '    ' 
!      write(*,*) ':::::::::::::: RESUSPENSION :::::::::::::::::::::::::'                       
!      write(*,*)  'res [ug s-1] hg-silt:', Rhgsil,'hg-pom:', Rhgpom   
!      write(*,*)  'res [ug s-1] mhg-silt:', Rmhgsil,'mhg-pom:', Rmhgpom  
!      write(*,*) '    '  
!      write(*,*)  'Vr:',Vr,'hgit:',hgit,'fsilt1:', fsilt1       
!      write(*,*) ':::::::::::::::::::::::::::::::::::::::::::::::::::::'      
!      write(*,*) ' '  
!------------------------------------------------------------------------    
! --------- Hg MeHg BURIAL  ---------------------------------------------
!------------------------------------------------------------------------
      bursHg  = (hgit*fsilt1 *Bvels*area) ![ug s-1]=[ug m-3]*[m s-1]*[m2]
      burpHg  = (hgit*fpom1  *Bvelp*area)  
       
      bursMHg = (mehgt*fsilt2*Bvels*area) ![ug s-1]=[ug m-3]*[m s-1]*[m2]
      burpMHg = (mehgt*fpom2*Bvelp*area)          
     
!      if (kext .EQ. 2284) then
!      write(665,*) Bvels, Bvelp, 'merc_SED'
!      write(664,*) silt,POM, 'merc_SED'
!      end if     

!      write(*,*) '    ' 
!      write(*,*) '::::::::::::::: BURIAL ::::::::::::::::::::::::'                       
!      write(*,*)  'bur [ug s-1] hg-silt:', bursHg, 'hg-pom',burpHg
!      write(*,*)  'bur [ug s-1] mhg-silt:', bursMHg, 'mhg-pom',burpMHg
!      write(*,*) '    '  
!      write(*,*)  'Bvels:', Bvels,'fsilt1:', fsilt1,'fpom1:', fpom1              
!      write(*,*) '::::::::::::::::::::::::::::::::::::::::::::::'      
!      write(*,*) '    ' 
!-----------------------------------------------------------------------
! --------- Hg MeHg DIFFUSION ------------------------------------------
!-------------------------------------------------------part 08---------
!---------- leggere Oxygen conc from biogeochem. model ? ---------------
! --Soerensen et al., 2016 --- if Oxy<0 :  Dchg_25_ANOX = 10.0/10**10 --
! ----------------------------------------------------------------------
! ---Sunderland et al., 2010 -------------------------------------------
! ---------- Dchg_25 = 2.0/10.**10.- Hg aggregated with macromolecular --
! -------------------------------------- colloidal organic matter -------
!-------- Dchg_25 = 9.5/10.**10.- HgCl4----------------------------------
! ----------------------------------------------------------------------- 
! ---------- Dchg_25 = 1.2/10.**9. -- MeHg sulfides (CH3HgSH0)----------- 
! -----------------------------------------------------------------------
      HgDw  = hgitw*(fdiss1w+fdoc1w)  ! Hg diss acque [ng l-1(w+s)] - da model water
      MeHgDw  = mehgtw*(fdiss2w+fdoc2w)  ! Hg diss acque [ng l-1(w+s)] - da model water
      Dchg_25 = 2.0/(10.**10.)                ! [m2 s-1]
      Dchg_t = Dchg_25/(1.+0.048*(25.-temp))  ! temperature dependent difffusion coefficient
      tor = 1. - log(por**2.)
     
      Hg2dpw  =  hgit  * (fdoc1+faq1)/por    !Soerensen et al., 2016  
      MeHgdpw =  mehgt * (fdoc1+faq1)/por

      num3 = (por*Dchg_t)/tor        
      num4 = (Hg2dpw - HgDw)/(depth/2.)  ![ug m-3] * [m-1] -> [ug m-4]  Hg2dpw = Hg2D/por?     
      num5 = (MeHgdpw - MeHgDw)/(depth/2.)  ![ug m-3] * [m-1] -> [ug m-4]  Hg2dpw = Hg2D/por? 

      ft1 =  - num3*num4    ![ug m-4] * [m2 s-1] -> [ug m-2 s-1]
      ft2 =  - num3*num5     
    
      JHgD     = ft1*area  ![ug s-1] < 0: da sedimenti ad acque
      Jngm2d     = ft1*1000.*86400. ![ng m-2] < 0: da sedimenti ad acque

!    ! calculated:  7 - 570       [ng m-2 d-1] Emili et al. 2012
!    ! measured  : 2.000 - 70.000 [ng m-2 d-1] Emili et al. 2012
      
      JMHgD     = ft2*area  ![ug s-1] < 0: da sedimenti ad acque
      JMngm2d   = ft2*1000.*86400. ![ug m-2] < 0: da sedimenti ad acque
!     
!      JHgD_kgy     = JHgD *365/10**9! kg y-1
!      write(*,*) '    ' 
!      write(*,*) '::::::::::::::: DIFFUSION FLUX ::::::::::::::::::::::'         
!      write(*,*) 'Jngm2d:', Jngm2d,num3,num4,Hg2dpw,MeHgdpw
!      write(*,*) 'hgit,mehgt', hgit,mehgt,fdoc1,faq1
!      write(*,*) 'hgdw,depth', hgdw,depth,Hg2dpw
!      write(*,*) '7 - 570       [ng m-2 d-1] Emili et al. 2012'
!      if (JHgD .LT. 0.0) then 
!      write(*,*) 'diff flux from sediment to water'
!      else if (JHgD .GT. 0.0) then 
!      write(*,*) 'diff flux from water to sediment'
!      else 
!      write(*,*) 'diff flux = 0'
!
!      write(*,*) 'JHgD', JHgD   
!      write(*,*) '::::::::::::::: DIFFUSION FLUX ::::::::::::::::::::::'               
!      write(*,*) '    ' 
!      write(*,*) 'hgit' , hgit, 'mehgt', mehgt

!-----------------------------------------------------------------------
! --------- Hg methylation and MeHg demethylation ---------------------- 
! -----------------------------------------------------------------------
      cksmet=(ksmet*Qbac**((temp-20.)/10.))/86400.       ! [d-1] to [s-1]
      sksme=cksmet*(hgit*(faq1+fdoc1))*vol    ! [s-1] * [ug m-3] = ug s-1

! -----------------------------------------------------------------------
! --------- Soerensen et al. 2016 FOR WATER COLUMN ----------------------
! -- OXIC COND.----- kmet = rmr/100 ! rmr = remineralization rate--------
! -- ANOXIC COND.--- kmet = (PO4 - PO4subox)*0.0005----------------------
! -----------------------------------------------------------------------
! ---met - demet in sedimentfrom observations and ratio MeHg/Hg----------
! -----------------------------------------------------------------------

      Eadem=10.      !WASP manual around 10 kcal/mol
      tkel=temp+273.
      tkref=293.     !reference temperature, K
      R=8.314472        ![J K-1 mol-1] universal gas constant
      Rcal=R/4.184
      deltat = tkel-tkref

      cksdem=(ksdem*exp(Eadem*1000*(deltat/(Rcal*tkel*tkref))))/86400
      sksdem=(cksdem*mehgt*(faq2+fdoc2))*vol ! [s-1] * [ug m-3] = ug s-1

! -------------------------------------------------------------------------------- 
! --------- part 09--------------------------------------------------------------- 
! ------------------------integration--------------------------------------------- 
!       Check Ginevra, 
!        C(1) = hgit              ! ug m-3
!        C(2) = mehgt

      Rhgsil=0.0
      Rhgpom=0.0
      Rmhgsil=0.0
      Rmhgpom=0.0
!      JHgD=0.0
!      JMHgD=0.0    
 
      if (bursHg .GE. 0.) then       
      CD(1) = +Shgsil +Shgpom +(-Rhgsil -Rhgpom +JHgD -bursHg -burpHg 
     & +sksdem-sksme)*86400
      else 
      CD(1) = +Shgsil +Shgpom +(-Rhgsil +JHgD -Rhgpom +sksdem-sksme
     & )*86400
      end if
      if (bursMHg .GE. 0.) then       
      CD(2) = +Smhgsil+Smhgpom+(-Rmhgsil-Rmhgpom+JMHgD-bursMHg-burpMHg
     & -sksdem+sksme)*86400
      else 
      CD(2) = +Smhgsil+Smhgpom+(-Rmhgsil-Rmhgpom+JMHgD -sksdem+sksme
     & )*86400
      end if    

       !write (889,*) (C(m), m=1,nvmerc),k    !'HgSED vars old'

!      if (kext .EQ. 1372) then
!          write (666,*) (CD(m), m=1,nvmerc),k
!      write (666,*) Shgsil/area
!      write (667,*) Shgpom/area
!      write (668,*) Rhgsil/area
!      write (669,*) Rhgpom/area
!      write (670,*) JHgD/area
!      write (671,*) bursHg/area
!      write (672,*) burpHg/area
!      write (673,*) vol
!      write (674,*) area
!      write (675,*) (C(m), m=1,nvmerc),k
!      elseif (kext .EQ. 3216) then 
!          write (666,*) (CD(m), m=1,nvmerc),k
!          write (667,*) (C(m), m=1,nvmerc),k
!      elseif (kext .EQ. 2407) then
!          write (666,*) (CD(m), m=1,nvmerc),k
!          write (667,*) (C(m), m=1,nvmerc),k
!      elseif (kext .EQ. 2654) then
!          write (666,*) (CD(m), m=1,nvmerc),k
!          write (667,*) (C(m), m=1,nvmerc),k
!      elseif (kext .EQ. 2341) then
!          write (666,*) (CD(m), m=1,nvmerc),k
!          write (667,*) (C(m), m=1,nvmerc),k
!      elseif (kext .EQ. 2150) then
!          write (666,*) (CD(m), m=1,nvmerc),k
!          write (667,*) (C(m), m=1,nvmerc),k
!      elseif (kext .EQ. 3762) then
!          write (666,*) (CD(m), m=1,nvmerc),k
!          write (667,*) (C(m), m=1,nvmerc),k
!      elseif (kext .EQ. 3985) then
!          write (666,*) (CD(m), m=1,nvmerc),k
!          write (667,*) (C(m), m=1,nvmerc),k
!      end if 
      
        CDw=0
        CDw(2)=(Rhgsil+Rhgpom-JHgD)*86400   !Resuspension and pw diffusion to the water                                                           
        CDw(3)=(Rmhgsil+Rmhgpom-JMHgD)*86400    !Resuspension and diffusion from ug s-1 to ug d-1                
 

       if(C(1).LT.0) then
       write(*,*) 'Hg2Sed<=0',C(1),'at node',ipext(k),'merc_sed bef'
       stop
       end if

       if(C(2).LT.0) then
       write(*,*) 'MHgSed c<=0',C(2),'at node',ipext(k),'merc_sed bef'
       stop
       end if

       if(Cw(2) .LT.0) then
       write(*,*) 'Hg2W<=0',Cw(2),C(1),'node',ipext(k),'merc_sed bef'
!       write(667,*) 'Hg2W<=0',Cw(2),C(1),'node',ipext(k),'merc_sed bef'
!       stop
       end if

       if(Cw(3).LT.0) then
       write(*,*) 'MHgW<=0',Cw(3),C(2),'node',ipext(k),'merc_sed bef'
!       write(778,*) 'MHgW<=0',Cw(3),C(2),'node',ipext(k),'merc_sed bef'
!       stop
       end if


       call merc_euler_sed(nvmerc,dtday,vol,vol,c,cold,cd) !claurent-OGS:add second volume
!       call merc_euler_sed (3,dtday,vol,vol,cw,coldw,cdw) !claurent-OGS:add second volume

!       call merc_euler (nvmerc,dtday,vol,c,cold,cd) 
       call merc_euler(3,dtday,vol,cw,coldw,cdw) 

       if(C(1) .LT.0) then
       write(*,*) 'Hg2Sed<=0',C(1),'at node',ipext(k),'merc_sed aft'
       stop
       end if

       if(C(2).LT.0) then
       write(*,*) 'MHgSed c<=0',C(2),'at node',ipext(k),'merc_sed aft'
       stop
       end if

       if(Cw(2) .LT.0) then
       write(*,*) 'Hg2W<=0',Cw(2),C(1),'at node',ipext(k),'merc_sed aft'
!       stop
       end if

       if(Cw(3).LT.0) then
       write(*,*) 'MHgW<=0',Cw(3),C(2),'node',ipext(k),'merc_sed aft'
!       stop
       end if 

!	write(*,*) 'Hg2sed' , Hg2sed, 'MeHgsed', MeHgsed
!      write(*,*) ' ' 
!      write(*,*) 'hgit' , hgit, 'mehgt', mehgt
!      write(*,*) ' '    
!	------------------------------------------------------------------
      Hgpsilt= hgit * fsilt1  ! [ug m-3] or [ng(hg) l(w+s)-1]
      Hg2POM = hgit * fPOM1   
      Hg2doc = hgit * fdoc1
      Hg2Cl =  hgit * faq1  

!      write (*,*) '' 
!      write (*,*) 'Hgpsilt',Hgpsilt, 'Hg2POM',Hg2POM
!      write (*,*) ''
!      write (*,*) 'Hg2doc',Hg2doc, 'Hg2Cl',Hg2Cl 
!      write (*,*) ''
!      write (*,*) 'hgit', hgit

        
       MeHgsilt= mehgt * fsilt2  ![ug(hg) m-3(w+s)] or [ng l-1]
       MeHgPOM = mehgt * fPOM2
       MeHgdoc = mehgt * fdoc2
       MeHgCl  = mehgt * faq2  ! mehgD in the tot volume (s+w)   

       Hg2P  = (Hgpsilt + Hg2POM)       ![ug m-3(w+s)] or [ng l-1]  
       MeHgP = (MeHgsilt + MeHgPOM)        

      Hg2D = (Hg2Cl + Hg2doc) 
      Hg2pw = Hg2D/por
     
      MeHgD = (MeHgCl + MeHgdoc) ! [ug m-3(w+s)] or [ng l-1] 
      MeHgpw  = MeHgD/por

      add1 = Hgpsilt/silt    !-- NO, viene molto diverso
      add2 = Hg2POM/POM

      Hg2sed = Hg2P/(silt+POM)
      MeHgsed = MeHgP/(silt+POM)

!     Cdoc  = por*DOC/10.**6. ! g(DOC)/m3(w) * 0.8 m3(w)/m3(s+w) -> mg(DOC)/L(s+w)/10**6 -> kgl-1
!     Csilt = silt/10.**6.  ! from [g m-3] to [g cm-3] or [kg l-1]      
!     Cpom  = POM /10.**6.  ! [kg l-1]   
!     Ctot  = Csilt+Cpom 

!      write (*,*) '' 
!       write (*,*) 'Hg2P', Hg2P
!      write (*,*) ''
!       write (*,*) 'Hg2D ngl', Hg2D
!      write (*,*) ''  
!       write (889,*) (C(m), m=1,nvmerc),'k',k, 'HgSED vars new'
      
!      write (91,*) Hgpsilt, Hg2POM,Hg2doc,Hg2Cl 
!      write (92,*) Shgsil, Shgpom,Rhgsil, Rhgpom,JHgD, bursHg,burpHg
!      write (93,*) Smhgsil,Smhgpom,Rmhgsil,Rmhgpom,JMHgD,bursMHg,burpMHg
!      write (94,*) MeHgsed, MeHgpw 
!      write (95,*) MeHgsilt, MeHgPOM,MeHgdoc,MeHgCl 
!      write (96,*) sksme, sksdem

!------ dissolved hg and mehg in porewater ------------------------------
!     hgd = (faq/por)*hgit    !Thomann and Di Toro, 1983
!       write(*,*) ' '
!       write(*,*) 'Hg2sed', Hg2sed ,'MeHgsed', MeHgsed
!       write (*,*) ''

!     stop
!	Hg0sed=C(1)
!     Hg2sed=C(1)+DepHg2-burHg
!     MeHgsed=C(3)+ DepMeHg-burMeHg

!     write (89,*) time, i,(C(m), m=1,nvmerc) !, ' variables'
!     write (91,*) DepHg2,DepMeHg,burHg,burMeHg,dt !, ' variables'        
! --------------------------------------------------------
!--------------------------------------------------

          if (k .GE. 1) then
           kext=ipext(k)
           fortfilenum=-1
           if(kext==3985)then
               fortfilenum=450
           elseif(kext==3986) then
               fortfilenum=451
           elseif(kext==3982) then
               fortfilenum=452
           elseif(kext==4007) then
               fortfilenum=453
           elseif(kext==3763) then
               fortfilenum=454
           elseif(kext==3765) then
               fortfilenum=455
           elseif(kext==3764)then
               fortfilenum=456
           elseif(kext==3762) then
               fortfilenum=457
           else if(kext==2150)then
               fortfilenum=458
           elseif(kext==2009) then
               fortfilenum=459
           elseif(kext==2359) then
               fortfilenum=460
           elseif(kext==2358) then
               fortfilenum=461
           elseif(kext==2341) then
               fortfilenum=462
           elseif(kext==2408) then
               fortfilenum=463
           elseif(kext==2191)then
               fortfilenum=464
           elseif(kext==2192) then
               fortfilenum=465
           elseif(kext==2654)then
               fortfilenum=466
           elseif(kext==2505) then
                fortfilenum=467
           elseif(kext==2655) then
               fortfilenum=468
          elseif(kext==2653) then
               fortfilenum=469
           else if(kext==1372)then
               fortfilenum=470
           elseif(kext==1375) then
               fortfilenum=471
           elseif(kext==1331) then
               fortfilenum=472
           elseif(kext==1378) then
               fortfilenum=473
           elseif(kext==3216) then
               fortfilenum=475
           elseif(kext==3057)then
               fortfilenum=476
           elseif(kext==2953) then
               fortfilenum=477
           elseif(kext==3217) then
               fortfilenum=478
           elseif(kext==2405) then
               fortfilenum=479
           elseif(kext==2407)then
               fortfilenum=480
           elseif(kext==2284) then
               fortfilenum=481
           elseif(kext==2404) then
               fortfilenum=482
           endif
           if(fortfilenum.ge.0)then
!               if(fortfilenum==450)
!     +             write(*,*) 'stamp to file 350... at iter=',iter,
!     +             ', tday=', tday
               write(fortfilenum,"(7(f20.7,','))") &
     &         Hg2sed,hgit,mehgt,silt,POM,add1,add2
! Hg2sed, MeHgsed,Hg2pw,MeHgpw
           endif
         endif

      end 
!	
