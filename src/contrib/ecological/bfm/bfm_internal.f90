
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016,2018-2020  Georg Umgiesser
!    Copyright (C) 2016,2019  Erik Pascolo
!    Copyright (C) 2016  Leslie Aveytua
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

! routines for BFM module
!
! revision log :
!
! 22.02.2016	ggu&erp	new bfm routines created from newconz
! 22.03.2016	ggu	changed VERS_7_5_6
! 06.06.2016	ggu	initialization from file changed
! 09.06.2016 	laa	modified (Search LAA)
! 10.06.2016	ggu	changed VERS_7_5_13
! 17.06.2016	ggu	changed VERS_7_5_15
! 28.06.2016	ggu	initialize bfmv, new routine bfm_check()
! 09.09.2016	ggu	changed VERS_7_5_17
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 01.05.2019	erp	some routines written for linkage to BFM
! 17.10.2019	ggu	these routines eliminated into own file
! 17.02.2020	ggu	femtime eliminated
! 09.03.2020	ggu	compiler bug fixed (ddt)

!==================================================================
        module mod_bfm_internal
!==================================================================

#include "BFM_module_list.h"

        implicit none

	!integer, save :: ibfm_state = 50	!number of state variables

!==================================================================
	contains
!==================================================================

!==================================================================
        end module mod_bfm_internal
!==================================================================

        subroutine bfm_init_internal

        use levels, ONLY: nlvdi
	implicit none

        call BFM0D_NO_BOXES(nlvdi,1,1,nlvdi,1)
        call BFM0D_INIT_IO_CHANNELS()
        call Initialize()

	write(6,*) 'bfm_init_internal: bfm routines initialized'
        
        end subroutine

!*********************************************************************

      subroutine bfm_reactor_internal
      
      use basin
      use levels
      use mod_sinking
      use mod_bfm

      implicit none

      integer :: i,l,k,bottom,jtr
      double precision :: ddt
      double precision, dimension(nlvdi,14) :: er
      double precision, dimension(4,nlvdi) :: wsink_pft
      logical :: bsur,bot
      double precision, dimension(nlvdi,jptra) :: a_apx
      double precision, dimension(jptra,nlvdi) :: b
      double precision, dimension(jptra_dia,nlvdi) :: d
      double precision, dimension(jptra_dia_2d) :: d2    
      double precision, dimension(nlv,nkn) :: matrix_light
      integer ipext 
      matrix_light(:,:) = 0.
      call light_abs(nkn,nlv,ibfm_state,bfmv,matrix_light)
      call get_ddt(ddt)

      a_apx(:,:) = 1.
      er(:,:)    = 1.
      er(:,14)   = 8.1
       
      DO k=1,nkn
          bottom = ilhkv(k)
!         write(*,*) 'bottom=',bottom
	  a_apx(1:bottom,:) = DBLE(bfmv(1:bottom,k,:))
          !if(k==1) &
!              write(*,*)'before calling BFM0D_Output_EcologyDynamics,', &
!                       ' a_apx(:,11)=',a_apx(1,:)
	  er(1:bottom,1)  = DBLE(tempv(1:bottom,k))!tn (ji,jj,jk)         ! Temperature (Celsius)
	  er(1:bottom,2)  = DBLE(saltv(1:bottom,k)) !(ji,jj,jk)           ! Salinity PSU
	  er(1:bottom,3)  = 1025.0D0 + rhov(1:bottom,k)                   ! Density Kg/m3
	  er(1,       4)  = 0.0D0                                         ! from 0 to 1 adimensional
	  er(1,       5)  = 390.0D0                                       ! CO2 Mixing Ratios (ppm)  390
	  er(1:bottom,6)  = (max(matrix_light(1:bottom,k),0.1))/0.217D0   ! PAR diatoms umoles/m2/s | Watt to umoles photons W2E=1./0.217 (LAA)
	  er(1:bottom,7)  = (max(matrix_light(1:bottom,k),0.1))/0.217D0   ! PAR flagellates umoles/m2/s | Watt to umoles photons W2E=1./0.217 (LAA)
	  er(1:bottom,8)  = (max(matrix_light(1:bottom,k),0.1))/0.217D0   ! PAR picophytoplankton umoles/m2/s | Watt to umoles photons W2E=1./0.217 (LAA)
	  er(1:bottom,9)  = (max(matrix_light(1:bottom,k),0.1))/0.217D0   ! PAR dinoflagellates umoles/m2/s | Watt to umoles photons W2E=1./0.217 (LAA)
	  er(1:bottom,10) = (max(matrix_light(1:bottom,k),0.1))/0.217D0   ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217 (LAA)
	  er(1,       11)  = 24.0D0                                       ! fotoperiod expressed in hours
	  er(1:bottom,12)  = DBLE(hdknv(1:bottom,k))                      ! depth in meters of the given cell
	  er(1       ,13)  = DBLE(sqrt(wxv(k)**2+wxv(k)**2))              ! vatm(ji,jj) * surf_mask(jk) ! wind speed (m/s)
	  er(1:bottom,14) = 8.0D0                                         ! PH
	  
!!!#ifdef BFM_active     
          !write(*,*)'ipext=',ipext(k),' for l=',l 

	  call BFM1D_Input_EcologyDynamics(bottom,a_apx,jptra,er)

	  call BFM1D_reset()

	  call EcologyDynamics()
	  
          call BFM1D_Output_EcologyDynamics(b, wsink_pft, d, d2)
          
!!!#endif         !
!Particolate (default at 3.0 m/day)
         wsink_bfm(1:bottom,k,ppR6c) = 3.0/86400.
         wsink_bfm(1:bottom,k,ppR6n) = 3.0/86400.
         wsink_bfm(1:bottom,k,ppR6p) = 3.0/86400.
         wsink_bfm(1:bottom,k,ppR6s) = 3.0/86400.
! Diatoms
         wsink_bfm(1:bottom,k,ppP1c) = REAL(wsink_pft(1,1:bottom),4)/86400.
         wsink_bfm(1:bottom,k,ppP1n) = REAL(wsink_pft(1,1:bottom),4)/86400.
         wsink_bfm(1:bottom,k,ppP1p) = REAL(wsink_pft(1,1:bottom),4)/86400.
         wsink_bfm(1:bottom,k,ppP1l) = REAL(wsink_pft(1,1:bottom),4)/86400.
         wsink_bfm(1:bottom,k,ppP1s) = REAL(wsink_pft(1,1:bottom),4)/86400.
! Flagellates
         wsink_bfm(1:bottom,k,ppP2c) = REAL(wsink_pft(2,1:bottom),4)/86400.
         wsink_bfm(1:bottom,k,ppP2n) = REAL(wsink_pft(2,1:bottom),4)/86400.
         wsink_bfm(1:bottom,k,ppP2p) = REAL(wsink_pft(2,1:bottom),4)/86400.
         wsink_bfm(1:bottom,k,ppP2l) = REAL(wsink_pft(2,1:bottom),4)/86400.
! Picophytoplankton
         wsink_bfm(1:bottom,k,ppP3c) = REAL(wsink_pft(3,1:bottom),4)/86400.
         wsink_bfm(1:bottom,k,ppP3n) = REAL(wsink_pft(3,1:bottom),4)/86400.
         wsink_bfm(1:bottom,k,ppP3p) = REAL(wsink_pft(3,1:bottom),4)/86400.
         wsink_bfm(1:bottom,k,ppP3l) = REAL(wsink_pft(3,1:bottom),4)/86400.
! Dinoflagellates
         wsink_bfm(1:bottom,k,ppP4c) = REAL(wsink_pft(4,1:bottom),4)/86400.
         wsink_bfm(1:bottom,k,ppP4n) = REAL(wsink_pft(4,1:bottom),4)/86400.
         wsink_bfm(1:bottom,k,ppP4p) = REAL(wsink_pft(4,1:bottom),4)/86400.
         wsink_bfm(1:bottom,k,ppP4l) = REAL(wsink_pft(4,1:bottom),4)/86400.

	  do jtr=1,jptra 
	      bfmv(1:bottom,k,jtr) = REAL(a_apx(1:bottom,jtr) + b(jtr,1:bottom)*ddt,4)
          enddo
!         phosphorus quota mmolP/mgC
	  bfm_d(1:bottom,k,1) = REAL(d(ppqpcPPY_iiP1,1:bottom))            !LAA
	  bfm_d(1:bottom,k,2) = REAL(d(ppqpcPPY_iiP2,1:bottom))
	  bfm_d(1:bottom,k,3) = REAL(d(ppqpcPPY_iiP3,1:bottom))
          bfm_d(1:bottom,k,4) = REAL(d(ppqpcPPY_iiP4,1:bottom))
!         nitrogen quota mmolN/mgC
          bfm_d(1:bottom,k,5) = REAL(d(ppqncPPY_iiP1,1:bottom))
          bfm_d(1:bottom,k,6) = REAL(d(ppqncPPY_iiP2,1:bottom))
          bfm_d(1:bottom,k,7) = REAL(d(ppqncPPY_iiP3,1:bottom))
          bfm_d(1:bottom,k,8) = REAL(d(ppqncPPY_iiP4,1:bottom))
!         silicon quota  mmolS/mgC
          bfm_d(1:bottom,k,9) = REAL(d(ppqscPPY_iiP1,1:bottom))
!         light limitation adimensional
          bfm_d(1:bottom,k,10) = REAL(d(ppeiPPY_iiP1,1:bottom))
          bfm_d(1:bottom,k,11) = REAL(d(ppeiPPY_iiP2,1:bottom))
          bfm_d(1:bottom,k,12) = REAL(d(ppeiPPY_iiP3,1:bottom))
          bfm_d(1:bottom,k,13) = REAL(d(ppeiPPY_iiP4,1:bottom))
!         carbon correction mgC/m3/d
          bfm_d(1:bottom,k,14) = REAL(d(ppuc_iiP1,1:bottom))
          bfm_d(1:bottom,k,15) = REAL(d(ppuc_iiP2,1:bottom))
          bfm_d(1:bottom,k,16) = REAL(d(ppuc_iiP3,1:bottom))
          bfm_d(1:bottom,k,17) = REAL(d(ppuc_iiP4,1:bottom))

          bfm_d(1:bottom,k,18) = REAL(d(ppun_iiP1,1:bottom))
          bfm_d(1:bottom,k,19) = REAL(d(ppun_iiP2,1:bottom))
          bfm_d(1:bottom,k,20) = REAL(d(ppun_iiP3,1:bottom))
          bfm_d(1:bottom,k,21) = REAL(d(ppun_iiP4,1:bottom))

          bfm_d(1:bottom,k,22) = REAL(d(ppup_iiP1,1:bottom))
          bfm_d(1:bottom,k,23) = REAL(d(ppup_iiP2,1:bottom))
          bfm_d(1:bottom,k,24) = REAL(d(ppup_iiP3,1:bottom))
          bfm_d(1:bottom,k,25) = REAL(d(ppup_iiP4,1:bottom))
!         Temperature deg Celsius
          bfm_d(1:bottom,k,26) = REAL(d(ppETW,1:bottom))
!         Salinity
          bfm_d(1:bottom,k,27) = REAL(d(ppESW,1:bottom))
!         Density
          bfm_d(1:bottom,k,28) = REAL(d(ppERHO,1:bottom))
!         PAR
          bfm_d(1:bottom,k,29) = REAL(d(ppPAR_phyto1,1:bottom))
! Flux variables remember to add jptra_var 

! Gross Primary production
          bfm_d(1:bottom,k,30) = REAL(d(jptra_var + ppruPPYc,1:bottom))
! Phytoplankton respiration
          bfm_d(1:bottom,k,31) = REAL(d(jptra_var + ppresPPYc,1:bottom))
! Predation
          bfm_d(1:bottom,k,32) = REAL(d(jptra_var + ppruZOOc,1:bottom))
! Bacterial remineralization
          bfm_d(1:bottom,k,33) = REAL(d(jptra_var + ppremPBAn,1:bottom))
          bfm_d(1:bottom,k,34) = REAL(d(jptra_var + ppremPBAp,1:bottom))
! POC remineralization
          bfm_d(1:bottom,k,35) = REAL(d(jptra_var + ppremPOC,1:bottom))
! pH
          bfm_d(1:bottom,k,36) = REAL(d(pppH,1:bottom))
! uptake and grazing
          bfm_d(1:bottom,k,37) = REAL(d(jptra_var + ppruPBAc,1:bottom))
          bfm_d(1:bottom,k,38) = REAL(d(jptra_var + ppruMEZ1c,1:bottom))
          bfm_d(1:bottom,k,39) = REAL(d(jptra_var + ppruMEZ2c,1:bottom))
          bfm_d(1:bottom,k,40) = REAL(d(jptra_var + ppruMIZ1c,1:bottom))
          bfm_d(1:bottom,k,41) = REAL(d(jptra_var + ppruMIZ2c,1:bottom))
          bfm_d(1:bottom,k,42) = 0.0

!         if(k==1)    &
!              write(*,*)k,' BFM1D_Output_EcologyDynamics, b(1,1)=',b(1,1), &
!                        ', wsink=',wsink_bfm(1,1,ppP1c),', d(1,1)=',d(1,1),   &
!                        ', a_apx(1,1)=',a_apx(1,1)
	ENDDO
      end subroutine
      
!*********************************************************************

!*********************************************************************

      subroutine bfm_benthic_internal

      use basin
      use levels
      use mod_sinking
      use mod_bfm
      use mod_layer_thickness

      implicit none

      integer :: i,l,k,bottom,jtr
      real               :: dt
      real, dimension(ibfm_state_benthic) :: a_apx,b_apx
      real, dimension(ibfm_state_benthic) :: settling_flux, remin, diffusion_flux
      real               :: volnode,areanode
      real               :: vol,volnew,volold,volben,area,benthic_thickness,mixL
      real               :: tau_nodes(nkn)        !bottom stress
      real               :: wsink 
      real               :: fluxD,fluxR
      real               :: tau,Pd,MagR
      real               :: tCDs        ! Input critical shear stress for deposition
      real               :: tCE         ! Input critical shear stress for resuspension
      real               :: rateC, rateN, rateP, rateSi
      real               :: temp
      real               :: DcP_25,DcN_25,DcSi_25,DcC_25
      real               :: DcP_t,DcN_t,DcSi_t,DcC_t
      real               :: por,tor

      benthic_thickness = 5.0! 0.01
      wsink = 3.0
      tCDs = .8
      tCE  = .6
      por  = 0.7
      tor  = 1. - log(por**2.)

      rateC  = 6.*365.*86400.
      rateN  = 6.*365.*86400.
      rateP  = 6.*365.*86400.
      rateSi = 6.*365.*86400.

      DcP_25  = 9.13157E-05/86400.
      DcN_25  = 0.000327073/86400.
      DcSi_25 = 0.000205722/86400.
      DcC_25 =  0.000362645/86400.

      call get_timestep(dt)

      call bottom_stress(tau_nodes)

      DO k=1,nkn

         tau=tau_nodes(k)
         if (tau>1.) then
           tau=1.
         end if
      
         if (tau <= tCDs) then  ! Deposition 
           Pd = (1. - tau/tCDs)  
             else
           Pd = 0.
         end if

        if (tau > tCE) then     
         MagR = (tau/tCE - 1.)**2.       !Magnitude of resuspension and Erosion flux
        else                             
         MagR = 0.
        end if

          settling_flux(:)=0.
          remin(:)=0.
          diffusion_flux(:)=0.
          bottom = ilhkv(k)
          area   = areanode(bottom,k)
          volnew = volnode(bottom,k,+1)
          volold = volnode(bottom,k,-1)
          vol=volnew
          volben = area*benthic_thickness
          mixL   = 0.5*(benthic_thickness+volold/area)

          fluxD=  wsink/86400.*Pd*dt*area
          fluxR=  0.00001*MagR*dt*area ! Check 
!         write(*,*) 'bottom=',bottom


          settling_flux(6)= + bfmv(bottom,k,ppR6c) * fluxD - bfm_benthic(1,k,6) * fluxR
          settling_flux(7)= + bfmv(bottom,k,ppR6n) * fluxD - bfm_benthic(1,k,7) * fluxR
          settling_flux(8)= + bfmv(bottom,k,ppR6p) * fluxD - bfm_benthic(1,k,8) * fluxR
          settling_flux(9)= + bfmv(bottom,k,ppR6s) * fluxD - bfm_benthic(1,k,9) * fluxR

!         Diffusion
!         DeltaConc = Hg2dpw - HgDw
          temp=tempv(bottom,k)
          DcP_t = DcP_25/(1.+0.048*(25.-temp))
          DcN_t = DcN_25/(1.+0.048*(25.-temp))
          DcSi_t = DcSi_25/(1.+0.048*(25.-temp))
          DcC_t = DcC_25/(1.+0.048*(25.-temp))
!         Diff = -((por*DwT)/tor)*DeltaConc/mixL

          diffusion_flux(1)= - ((por*DcP_t)/tor)  * (bfm_benthic(1,k,1)  - bfmv(bottom,k,ppN1p))/mixL*area*dt 
          diffusion_flux(2)= - ((por*DcN_t)/tor)  * (bfm_benthic(1,k,2)  - bfmv(bottom,k,ppN3n))/mixL*area*dt 
          diffusion_flux(3)= - ((por*DcN_t)/tor)  * (bfm_benthic(1,k,3)  - bfmv(bottom,k,ppN4n))/mixL*area*dt 
          diffusion_flux(4)= - ((por*DcSi_t)/tor) * (bfm_benthic(1,k,4) - bfmv(bottom,k,ppN5s))/ mixL*area*dt 
          diffusion_flux(5)= - ((por*DcC_t)/tor)  * (bfm_benthic(1,k,5)  - bfmv(bottom,k,ppO3c))/mixL*area*dt 

          remin(1)        = +1.0/rateP  * bfm_benthic(1,k,8)*dt ! flux to PO4
          remin(2)        =  0.0                             ! flux to NO3
          remin(3)        = +1.0/rateN  * bfm_benthic(1,k,7)*dt ! flux to NH4
          remin(4)        = +1.0/rateSi * bfm_benthic(1,k,9)*dt ! flux to Si
          remin(5)        = +1.0/rateC  * bfm_benthic(1,k,6)*dt ! flux to DIC
          remin(6)        = -1.0/rateC  * bfm_benthic(1,k,6)*dt ! flux from R6c 
          remin(7)        = -1.0/rateN  * bfm_benthic(1,k,7)*dt ! flux from R6n
          remin(8)        = -1.0/rateP  * bfm_benthic(1,k,8)*dt ! flux from R6p
          remin(9)        = -1.0/rateSi * bfm_benthic(1,k,9)*dt ! flux from R6s

          do i=1,ibfm_state_benthic
              bfm_benthic(1,k,i)  = bfm_benthic(1,k,i) + remin(i) + diffusion_flux(i)/volben + settling_flux(i)/volben ! include specific sediment cell thickness
          enddo

          bfmv(bottom,k,ppN1p)= bfmv(bottom,k,ppN1p) - diffusion_flux(1)/vol - settling_flux(1)/vol
          bfmv(bottom,k,ppN3n)= bfmv(bottom,k,ppN3n) - diffusion_flux(2)/vol - settling_flux(2)/vol
          bfmv(bottom,k,ppN4n)= bfmv(bottom,k,ppN4n) - diffusion_flux(3)/vol - settling_flux(3)/vol
          bfmv(bottom,k,ppN5s)= bfmv(bottom,k,ppN5s) - diffusion_flux(4)/vol - settling_flux(4)/vol
          bfmv(bottom,k,ppO3c)= bfmv(bottom,k,ppO3c) - diffusion_flux(5)/vol - settling_flux(5)/vol

          bfmv(bottom,k,ppR6c)= bfmv(bottom,k,ppR6c) - settling_flux(6)/vol
          bfmv(bottom,k,ppR6n)= bfmv(bottom,k,ppR6n) - settling_flux(7)/vol
          bfmv(bottom,k,ppR6p)= bfmv(bottom,k,ppR6p) - settling_flux(8)/vol
          bfmv(bottom,k,ppR6s)= bfmv(bottom,k,ppR6s) - settling_flux(9)/vol


      enddo 

          
      end subroutine
      
!*********************************************************************
      subroutine light_abs(nkn,nlv,nst,bfmv,matrix_light)
      
	use levels, only: ilhkv
	use mod_layer_thickness
        use mod_bfm, ONLY:ppP1l,ppP2l,ppP3l,ppP4l


      implicit none

      integer,intent(in) :: nkn,nlv,nst
      real bfmv(nlv,nkn,nst)
      double precision,dimension(nlv,nkn),intent(inout)::matrix_light

      integer :: k,l
      real :: solrad
      double precision :: coeff,eps1,eps0
      
      eps1= 0.0088
      eps0=0.04

      do k=1,nkn
	  
	  call meteo_get_solar_radiation(k,solrad)
	  matrix_light(1,k) = solrad
      
	  do l=2,ilhkv(k)
	    coeff = -(eps1*(bfmv(l-1,k,ppP1l)+bfmv(l-1,k,ppP2l) +             & 
      	          bfmv(l-1,k,ppP3l)+bfmv(l-1,k,ppP4l))+eps0)*hdknv(l-1,k)      
	    matrix_light(l,k) = matrix_light(l-1,k) * exp(coeff) 
	  enddo
	  
      
	  do l=1,ilhkv(k)
	    coeff = -0.5*(eps1*(bfmv(l,k,ppP1l)+bfmv(l,k,ppP2l)+              &
      	          bfmv(l,k,ppP3l)+bfmv(l,k,ppP4l))+eps0)*hdknv(l,k)
	    matrix_light(l,k) = matrix_light(l,k) * exp(coeff) 
	  enddo
	  
      enddo

      end subroutine

!*********************************************************************

      SUBROUTINE write_matrix(mtw,righe,colonne,percorso)

          IMPLICIT NONE

          INTEGER :: i, j, righe,colonne
          character(LEN=*), INTENT(in) :: percorso
          double precision, INTENT(in) :: mtw(righe,colonne)
    
          OPEN(unit=115, file=percorso)
	  print *,righe,colonne
          Do i=1,righe
          WRITE(115,'(*(F14.7))') (mtw(i,j),j=1,colonne)
          END DO
          CLOSE(0)

       END SUBROUTINE write_matrix
      
!*********************************************************************


	subroutine bfm_init_for_restart()

        use mod_bfm
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
	
	implicit none
	real getpar

        if( mod_bfm_initialized .eqv. .false. ) then 
          ibfm=nint(getpar('ibfm'))

          call mod_bfm_init(ibfm_state,nkn,nlvdi)

	  call bfm_init_internal

          write(6,*) 'bfm initialized for restart: ', &
                     ibfm,ibfm_state,nkn,nlvdi
         endif
         mod_bfm_initialized=.True. 
         end subroutine bfm_init_for_restart

!       ---------------------------------------------

        subroutine write_restart_bfm(iunit)
          use mod_bfm
          use levels, only : nlvdi
          use basin, only : nkndi
          
          implicit none
          integer iunit
          integer l,k,s
            write(iunit) nst_bfm,nlv_bfm,nkn_bfm,-999990
            do s=1,nst_bfm
               write(iunit) ((bfmv(l,k,s),l=1,nlv_bfm),k=1,nkn_bfm)
            enddo
          return
        end subroutine write_restart_bfm

!       ---------------------------------------------

        subroutine skip_restart_bfm(iunit)
          use mod_bfm
          use levels, only : nlvdi
          use basin, only : nkndi
          
          implicit none
          integer iunit
          integer s,r_ns,r_nl,r_nk,r_flag
            read(iunit)  r_ns,r_nl,r_nk,r_flag
            if(r_ns/=nst_bfm .OR. r_nl/= nlv_bfm   .OR. & 
               r_nk/=nkn_bfm .OR. r_flag/=-999990) goto 98
           
            do s=1,nst_bfm
               read(iunit) 
            enddo

          return

   98       continue   
            write(6,*) 'error parameters in skip_restart_bfm '
            write(6,*) r_ns,r_nl,r_nk,r_flag 
            stop 'error stop skip_restart_bfm'

        end subroutine skip_restart_bfm

!       ------------------------------------

        subroutine read_restart_bfm(iunit)
          use mod_bfm
          use levels, only : nlvdi
          use basin, only : nkndi
          
          implicit none
          integer iunit
          integer l,k,s,r_ns,r_nl,r_nk,r_flag
            read(iunit)  r_ns,r_nl,r_nk,r_flag
            if(r_ns/=nst_bfm .OR. r_nl/= nlv_bfm   .OR. &
               r_nk/=nkn_bfm .OR. r_flag/=-999990) goto 98
         
            do s=1,nst_bfm
               read(iunit) ((bfmv(l,k,s),l=1,nlv_bfm),k=1,nkn_bfm)
            enddo

          return

   98       continue   
            write(6,*) 'error parameters in read_restart_bfm '
            write(6,*) r_ns,r_nl,r_nk,r_flag 
            stop 'error stop read_restart_bfm'

        end subroutine read_restart_bfm
