
!--------------------------------------------------------------------------
!
!    Copyright (C) 2014,2016  Donata Melaku Canu
!    Copyright (C) 2016,2019  Georg Umgiesser
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

! weutro_shell -  routines for benthic filter feeding linked to weutro
!
! revision log :
!
! 15.03.2014	dmc	routine for shellfish farming simulation
! 17.06.2016	dmc	last changes integrated
! 27.06.2016	ggu	changed VERS_7_5_16
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 21.05.2019	ggu	changed VERS_7_5_62
!
! notes :
!
! Reference: Melaku Caunu, Aveytua, Camacho-Ibar, Solidoro.
! Trophic properties and effect of upwelling  the San Quintin Bay. in prep.  
!
!********************************************************************

	subroutine wshell(k,t,dt,vol,depth,vel,stp,c,csh)

! EUTRO 0-Dimensional (Sediments)

        implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )
        integer nsstate          !total number of state parameters
        parameter( nsstate = 3 )
	real controllo

	integer k
	real t			!actual time [day]
        real dt                 !time step [day]
        real vol                !volume [m**3]
        real depth              !depth of box [m]
        real vel                !velocity [m/s]
        real stp                !temperature [C]
        real c(nstate)          !state variable [mg/L] == [g/m**3]
        real cold(nstate)       !old state variable [mg/L] == [g/m**3]
        real csh(nsstate)        !state variable (shellfarm) [mg/L] == [g/m**3]

	include 'donata.h'

	logical debug
	logical bshell
	integer i
	real phyto, shell
	real filtration,shellgrowth
	real shellsize
	real shellfarm, decay

	real kfilt, eff, kphy, kdec, kshell

        real cds(nstate)        !source term (right hand side) [g/day]
        real ca(nstate)         !auxiliary state vector [g/m**3]
        real caold(nstate)      !old auxiliary state vector [g/m**3]

        integer icall
        save icall
        data icall / 0 /

	bshell = .false.
	bshell = .true.

	if( .not. bshell ) return


        if( icall .eq. 0 ) then
          call settopseg(.true.)        !marks segment as surface
          call setbotseg(.true.)        !marks segment as bottom
          icall = 1
        end if




!	debug = k .eq. 100
!	debug = k .eq.-1 


!	kfilt= 0.001    !originale
        kfilt= 0.010
!	eff = 1.        !originale
        eff = 3.25
	kphy=0.8  	!semisaturation constant for phytop filtration 
!       kphy=0.5        !originale
	kdec=.001	!decay rate for oyster death
!	kshell=1.0	!semisaturation constant for shellfish growth (originale)
        kshell=0.8      !16/nov/2015

!	call tempcoef(stp,kpt,knt)

	phyto = c(4)

	shellfarm= csh(1)
	shellsize = csh(2)  
	
!	write(6,*) 'shellfarm1'
!	write (6,*) shellfarm,shellsize,phyto

	filtration = kfilt* shellfarm*shellsize*vol*phyto/(kphy+phyto)
	shellgrowth= filtration*eff/(kshell+shellsize)
	decay = shellfarm*kdec*vol
	
!	write(6,*) 'shellfarm2'
!	if(filtration.gt.0) then
!	write (6,*) filtration,shellgrowth, decay
!	end if
	
	ca(1) = phyto 
	ca(2) = shellfarm
	ca(3) = shellsize

	cds(1) = -filtration
	cds(2) = -decay	*1.068**(stp-20)!se move o no? proposta 25 marzo
	cds(3) = shellgrowth

!

	
!
       call euler(1,dt,vol,ca(1),caold,cds(1))     ! euler(1,dt,vol,ca(1),caold(1),cds(1)) 
       call euler(2,dt,vol,ca(2),caold(2),cds(2))  ! euler(1,dt,vol,ca(2),caold(2),cds(2))

	c(4)=ca(1)
	csh(1)=ca(2)
	csh(2)=ca(3)

!	controllo=csh(3)
!	if (controllo.gt.0) then
!	if( debug ) then
!	  write(6,*) 
!	  write(6,*),'ca', (ca(i),i=1,3)
!	  write(88,*),'csh', (csh(i),i=1,3)
!c	  write(6,*) '-----weutro_shell---------------------------'
!	end if
!	end if


	end

!********************************************************************

!	subroutine tempcoef(stemp,kpt,knt)
!
!	implicit none
!
!	real stemp,kpt,knt
!	real stemp20
!
!	include 'weutro.h'
!
!	stemp20 = stemp - 20.
!	kpt = K58T**STemp20
!	knt = K1013T**STemp20
!
!	end
!
!********************************************************************

!********************************************************************

