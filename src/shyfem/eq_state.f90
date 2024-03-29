
!--------------------------------------------------------------------------
!
!    Copyright (C) 1995,1998,2010,2015,2019  Georg Umgiesser
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

! equation of state routine
!
! contents :
!
! subroutine rhotst			to test function sigma
! real function sigma(s,t,p)		computes density
!
! revision log :
!
! 30.08.1995	ggu	$$AUST - austausch coefficient introduced
! 11.10.1995	ggu	$$BCLBND - boundary condition for barocliic runs
! 14.08.1998	ggu	new routine tsmass
! 19.08.1998	ggu	new routines to write NOS file
! 19.08.1998	ggu	call to barcfi changed
! 26.08.1998	ggu	call to bclini changed
! 26.08.1998	ggu	routines deleted: barcpr, chkuvw, vmima
! 26.08.1998	ggu	routines mimari transferred to subssv
! 26.08.1998	ggu	subroutine tsmass transferred to newchk
! 26.08.1998	ggu	subroutine stmima transferred to newut
! 26.08.1998	ggu	subroutines bclini, barcfi, bclbnd deleted
! 23.03.2010	ggu	changed v6.1.1
! 28.09.2010	ggu	changed VERS_6_1_11
! 12.10.2015	ggu	changed VERS_7_3_3
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 17.09.2020	ggu	renamed sigma to sigma_stp
!
!**********************************************************************
!
!	call rhotst
!	end
!
!**********************************************************************
!
	subroutine rhotst
!
! to test function sigma
!
	implicit none
!
	integer itot
	parameter (itot=12)
!
	integer i,is,it,ip
	real sd,td,pd
	real sigma_stp
	real rho,aux
	double precision rhoref(8)
	real s(2),t(2),p(2)
	real sv(itot),tv(itot),pv(itot)
!
	data sv /0.,0.,35.,35.,0.,0.,35.,35.  ,35.,35.,35.,35./
	data tv /0.,0.,5.,5.,0.,0.,5.,5.  ,5.,5.,5.,5./
	data pv /100.,101.,100.,101.,1000.,1010.,1000.,1010. &
     &			,0.,0.1, 100.,100.1/
!
	data s /0.,35./
	data t /5.,25./
	data p /0.,1000./
	data rhoref/ 999.96675,1044.12802, 997.04796,1037.90204 &
     &		   ,1027.67547,1069.48914,1023.34306,1062.53817 &
     &		   /
!
	write(6,1000)
!
	i=0
	do is=1,2
	  do it=1,2
	    do ip=1,2
		i=i+1
		sd=s(is)
		td=t(it)
		pd=p(ip)
		rho = sigma_stp(sd,td,pd)
		write(6,1010) sd,td,pd,rho,real(rhoref(i)-1000.)
	    end do
	  end do
	end do
!
	aux=0.
	do i=1,itot
	  rho = sigma_stp(sv(i),tv(i),pv(i))
	  write(6,1010) sv(i),tv(i),pv(i),rho,rho-aux
	  aux=rho
	  if(mod(i,2).eq.0) aux=0.
	end do
!
	stop
 1000	format('       salinity    temperature       pressure', &
     &         '        density', &
     &         '    ref-density')
 1010	format(5f15.6)
	end
!
!*****************************************************************
!
!---------------------------------------------------------------
!
!         d e n s i t y   o f   s e a - w a t e r
!
!           international unesco equation of state 1980
!                  (double precision version)
!
!            test data:
!
!          s=0.  t=5.  p=0.    rho= 999.96675
!          s=0.  t=5.  p=1000. rho=1044.12802
!          s=0.  t=25. p=0.    rho= 997.04796
!          s=0.  t=25. p=1000. rho=1037.90204
!          s=35. t=5.  p=0.    rho=1027.67547
!          s=35. t=5.  p=1000. rho=1069.48914
!          s=35. t=25. p=0.    rho=1023.34306
!          s=35. t=25. p=1000. rho=1062.53817
!
!          1000 m ~ 1010 dbar
!
!            input units (real) :
!
!         p - must be given as depth in 0.1*meters (bars)
!         s -       "       as salinity in promille
!         t -       "       as degree centigrade
!
!            output units (real) :
!
!         sigma = rho - 1000.   with rho in kg/m**3
!
!------------------------------------------------------------------
!
	real function sigma_stp(s,t,p)
!
! computes density
!
	implicit none
!
	real s,t,p
	double precision t1,t2,t3,t4,t5
	double precision s1,s2,s32
	double precision p1,p2
	double precision rho0,sbmk
!
	double precision one,thousd
	parameter( one = 1.0d+0 , thousd = 1.0d+3 )
!
	double precision r00
	double precision r01,r02,r03,r04,r05,r06,r07
	double precision r08,r09,r10,r11,r12,r13,r14
!
	parameter( r00 = +999.842594d+0  )
	parameter( r01 = +6.793952d-2 , r02 = -9.095290d-3 )
	parameter( r03 = +1.001685d-4 , r04 = -1.120083d-6 )
	parameter( r05 = +6.536332d-9 , r06 = +8.24493d-1  )
	parameter( r07 = -4.0899d-3   , r08 = +7.6438d-5   )
	parameter( r09 = -8.2467d-7   , r10 = +5.3875d-9   )
	parameter( r11 = -5.72466d-3  , r12 = +1.0227d-4   )
	parameter( r13 = -1.6546d-6   , r14 = +4.8314d-4   )
!
	double precision k00
	double precision k01,k02,k03,k04,k05,k06
	double precision k07,k08,k09,k10,k11
	double precision k13,k14,k15,k16,k17,k18,k19
	double precision k20,k21,k22,k23,k24,k25,k26
!
	parameter( k00 = +19652.21d+0 )
	parameter( k01 = +148.4206d+0 , k02 = -2.327105d+0 )
	parameter( k03 = +1.360477d-2 , k04 = -5.155288d-5 )
	parameter( k05 = +3.239908d+0 , k06 = +1.43713d-3  )
	parameter( k07 = +1.16092d-4  , k08 = -5.77905d-7  )
	parameter( k09 = +8.50935d-5  , k10 = -6.12293d-6  )
	parameter( k11 = +5.2787d-8 )
	parameter( k13 = +54.6746d+0  , k14 = -0.603459d+0 )
	parameter( k15 = +1.09987d-2  , k16 = -6.1670d-5   )
	parameter( k17 = +7.944d-2    , k18 = +1.6483d-2   )
	parameter( k19 = -5.3009d-4   , k20 = +2.2838d-3   )
	parameter( k21 = -1.0981d-5   , k22 = -1.6078d-6   )
	parameter( k23 = +1.91075d-4  , k24 = -9.9348d-7   )
	parameter( k25 = +2.0816d-8   , k26 = +9.1697d-10  )
!
	t1=t
	t2=t1*t1
	t3=t2*t1
	t4=t3*t1
	t5=t4*t1
	s1=s
	s2=s1*s1
!	s32=sqrt(real(s2*s1))
	s32=dsqrt(s2*s1)		!ggu
	p1=p
	p2=p1*p1
!
!--------------------> density at standard atmospheric pressure (p=0)
!
	rho0 =   r00 + r01*t1 + r02*t2 + r03*t3 + r04*t4 + r05*t5 &
     &	     + ( r06 + r07*t1 + r08*t2 + r09*t3 + r10*t4 ) * s1 &
     &       + ( r11 + r12*t1 + r13*t2 ) * s32 &
     &	     + r14*s2
!
!--------------------> secant bulk modulus
!
	sbmk = k00 &
     &	     + k01*t1 + k02*t2 + k03*t3 + k04*t4 &
     &       + ( k05 + k06*t1 + k07*t2 + k08*t3 ) * p1 &
     &       + ( k09 + k10*t1 + k11*t2 ) * p2 &
     &       + ( k13 + k14*t1 + k15*t2 + k16*t3 ) * s1 &
     &       + ( k17 + k18*t1 + k19*t2 + k23*p1 ) * s32 &
     &       + ( k20 + k21*t1 + k22*t2 ) * p1 * s1 &
     &       + ( k24 + k25*t1 + k26*t2 ) * p2 * s1
!
!--------------------> effect of sal. and temp. : sigma
!
	sigma_stp = rho0/(one-p/sbmk) - thousd
!
	return
	end
!
!**********************************************************************

