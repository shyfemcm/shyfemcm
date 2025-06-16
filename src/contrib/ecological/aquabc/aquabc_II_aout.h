
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015,2019  Georg Umgiesser
!    
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Parameters and arrays for writing state variables to ASCII files for defined nodes 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccz



!      Output of variables for water column

! revision log :
!
! 26.02.2015	ggu	changed VERS_7_1_6
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61

	      integer NBIOTSMX             !Maximal number of time series(nodes) for WC 
		  parameter(NBIOTSMX=100)    
	      integer NBIOTS               !Actual number of time series(nodes) for WC       
		  integer BIOTSFUN(NBIOTSMX)   !file units for WC state variables ASCII output	
		  integer BIOTSNOD(NBIOTSMX)   !nodes for WC state variables ASCII output

!      Output  variables for sediments
		   integer NBIOTSMX_sed               !Maximal number of time series(nodes) for BS
		   parameter(NBIOTSMX_sed=100)    
	       integer NBIOTS_sed                 !Actual number of time series(nodes) for BS	
		   integer BIOTSFUN_sed(NBIOTSMX_sed) !file units for BS state variables ASCII output 	
		   integer BIOTSNOD_sed(NBIOTSMX_sed) !nodes for BS state variables ASCII output



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Parameters and arrays for writing diagnostics(components of derivatives) to ASCII files for defined nodes 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! WATER COLUMN
		integer   NDIAGVAR           !Maximum number of different intermediate variables in output.  
		parameter (NDIAGVAR=26)      !Number of auxilary variables for output	
		integer   NDGTSMX            !Total number of nodes for intermediates of each state variable 
		parameter (NDGTSMX=100)
		integer NDGTS(nstate)        !Actual number of nodes for intermediates of each state variable
		integer DGTSFUN(NDGTSMX, nstate) !Output file units for EUTRO diagnostics time series(nodes) plots
		integer DGTSNOD(NDGTSMX, nstate) !Nodes for intermediate variables in output
		real    dgar(nstate,NDIAGVAR)    !dgar (Diagnostics data) for each sediments layer


! BOTTOM SEDIMENTS
		integer   NDIAGVAR_sed          !Maximum number of different intermediate variables in output. Chek biotser_write formats
		parameter (NDIAGVAR_sed=25)      !Number of auxilary variables for output	
		integer   NDGTSMX_sed           !Total number of nodes for intermediates of each state variable 
		parameter (NDGTSMX_sed=100)
		integer NDGTS_sed(nsstate)      !Actual number of nodes for intermediates of each state variable
		integer DGTSFUN_sed(NDGTSMX, nsstate) !Output file units for EUTRO diagnostics time series(nodes) plots
		integer DGTSNOD_sed(NDGTSMX, nsstate) !Nodes for intermediate variables in output
	    real    dgar_sed(NOSLAY,nsstate,NDIAGVAR_sed)    ! Diagnostics data, for temporary storage


