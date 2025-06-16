
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016,2018-2019  Georg Umgiesser
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
! 07.06.2016	ggu	changed VERS_7_5_12
! 25.10.2018	ggu	changed VERS_7_5_51
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 17.07.2019	ggu	error in calling nc_output_record()

!**************************************************************************

	subroutine nc_output_init(ncid,title,nvar,ivars,b2d,sncglobal)

	implicit none

	integer ncid			!id of file (return)
	character*(*) title		!name of simulation
	integer nvar			!total number of variables to be written
	integer ivars(nvar)		!variable id of SHYFEM
	logical b2d			!should write 2d fields
	character*(*) sncglobal		!file name of extra global values

	write(6,*) 'Cannot initialize NETCDF module'
	write(6,*) 'NETCDF support has not been enabled'
	write(6,*) 'Please enable NETCDF support in Rules.make'
	write(6,*) '(set "NETCDF=true")'

	stop 'error stop nc_output_init: no netcdf support'

	end

!********************************************************************

	subroutine nc_output_record(ncid,var_id,var_dim,np,cv3)

	use basin
	use levels

	implicit none

	integer ncid
	integer var_id
	integer var_dim
	integer np
	real cv3(nlvdi,np)

	end

!********************************************************************

	subroutine nc_output_hydro(ncid,znv,uprv,vprv)

	use basin
	use levels

	implicit none

	integer ncid
	real znv(nkn)
	real uprv(nlvdi,nkn)
	real vprv(nlvdi,nkn)

	end

!********************************************************************

	subroutine nc_output_time(ncid,dtime)

	implicit none

	integer ncid
	double precision dtime

	end

!********************************************************************

	subroutine nc_output_final(ncid)

	implicit none

	integer ncid

	end

!********************************************************************

        subroutine nc_output_set_vars(br,nx,ny,hc,fm,x,y)

        implicit none

        logical :: br
        integer :: nx
        integer :: ny
        real :: hc(nx*ny)
        real :: fm(4,nx*ny)
        real :: x(nx)
        real :: y(ny)

	end

!********************************************************************

        subroutine nc_output_get_var_id(iv,var_id)

        implicit none

        integer :: iv
        integer :: var_id

        end

!********************************************************************

        subroutine nc_output_get_var_dim(iv,var_dim)

        implicit none

        integer :: iv
        integer :: var_dim

        end

!********************************************************************


        subroutine nc_set_quiet(bquiet)

        implicit none

        logical bquiet

        end

!********************************************************************

