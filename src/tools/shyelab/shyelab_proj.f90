
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Georg Umgiesser
!    Copyright (C) 2018  Debora Bellafiore
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

! projection handling routines
!
! contents :
!
! subroutine shy_proj
!
! revision log :
!
! 20.07.2018	dbf	projection routines for shyelab
! 30.08.2018	ggu	bug fix if mode == 0
! 16.02.2019	ggu	changed VERS_7_5_60
!
!*****************************************************************

        subroutine shy_proj

        use clo
        use basin
        use shympi
        use projection
        use coordinates

        implicit none

	logical, save			:: bwrite = .true.
        integer                         :: mode,i,ioff,j,jj,iproj1
        double precision		:: c_param(9)
        double precision		:: d(9)
        real                            :: xgeov1(nkn),ygeov1(nkn) 
        character*80 sproj,tproj,text

	logical bas_is_spherical
	integer istot,iscand
        
        !mode: +1: cart to geo  -1: geo to cart                
        !iproj1 1=GB 2=UTM 3=CPP 4=UTM non standard
        
        call clo_get_option('proj',sproj)

!	-----------------------------------------------
!	get type of projection
!	-----------------------------------------------

	ioff = 1
	jj = istot(sproj,tproj,ioff)	!in tproj is name of projection
        if( jj == 0 ) return		!no projection requested

	write(6,*) 'setting up projection'

	if( bas_is_spherical() ) then
	  mode = -1
	  text = 'from spherical to cartesian'
	else
	  mode = +1
	  text = 'from cartesian to spherical'
	end if

	call skipwh(sproj,ioff)
	call s2i_proj(tproj,iproj1)	!get numerical code of projection
	if( iproj1 == 0 ) goto 98

	d = 0.
        jj=iscand(sproj(ioff:),d,6)		!mode,iproj,c_param
        if( jj == 0 ) goto 97

	c_param = 0.
	do j=1,jj
          c_param(j)=d(j)
	end do

	if( bwrite ) then
	  write(6,*) 'projecting ',trim(text)
	  write(6,*) 'using projection: ',trim(tproj),iproj1
	  write(6,*) 'parameters: ',(c_param(j),j=1,jj)
	end if

	!call set_projection(iproj1,c_param)
        call init_coords(iproj1,c_param)
        call convert_coords(mode,nkn,xgv,ygv,xgeov1,ygeov1)

        xgv = xgeov1
        ygv = ygeov1

	if( bas_is_spherical() ) then
	  call bas_set_spherical(0)
	else
	  call bas_set_spherical(1)
	end if

	return
   97	continue
	write(6,*) 'no parameters given for projection: ',trim(tproj)
	stop 'error stop shy_proj: no parameters'
   98	continue
	write(6,*) 'unknown projection: ',trim(tproj)
	stop 'error stop shy_proj: unknown projection'
   99	continue
        write(6,*)'coordinates are not lat,lon'
        write(6,*)'  xmin,xmax: ',minval(xgv),maxval(xgv)
        write(6,*)'  ymin,ymax: ',minval(ygv),maxval(ygv)
        stop 'error stop shy_proj: not lat/lon'
        end subroutine shy_proj
       
!*****************************************************************

