
!--------------------------------------------------------------------------
!
!    Copyright (C) 2014,2019  Georg Umgiesser
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

!-------------------------------------------------------------
! header file for ous file format
!-------------------------------------------------------------
!
! ftype		file id
! maxvers	newest version of file
! maxcomp	compatible function calls down to here
!
! ndim		number of possible entries (open files)
! nitdim	number of integer values to be stored
! nchdim	number of string values to be stored
!
! ousitem	number of maximum entries in table
!
! ousvar	integer parameters of open files
! ouschar	string parameters of open files
!
! ousvar(0,n)	iunit
! ousvar(1,n)	nvers
! ousvar(2,n)	nkn
! ousvar(3,n)	nel
! ousvar(4,n)	nlv
! ousvar(5,n)	not used
! ousvar(6,n)	date
! ousvar(7,n)	time
!
! ousreal(1,n)	href
! ousreal(2,n)	hzmin
!
! ouschar(1,n)	title
! ouschar(2,n)	femver
!
!-------------------------------------------------------------
! parameters
!-------------------------------------------------------------

! revision log :
!
! 28.01.2014	ggu	changed VERS_6_1_71
! 16.02.2019	ggu	changed VERS_7_5_60

!-------------------------------------------------------------


        integer ftype,maxvers,maxcomp
        parameter(ftype=27,maxvers=2,maxcomp=1)

        integer ndim,nitdim,nrldim,nchdim
        parameter(ndim=10,nitdim=7,nrldim=2,nchdim=2)

!-------------------------------------------------------------
! common
!-------------------------------------------------------------

        integer ousitem
        common /ousitm/ousitem

        integer ousvar(0:nitdim,ndim)
        common /ousvar/ousvar

        real ousreal(nrldim,ndim)
        common /ousreal/ousreal

        character*80 ouschar(nchdim,ndim)
        common /ouschar/ouschar

!-------------------------------------------------------------
! save
!-------------------------------------------------------------

        save /ousitm/
        save /ousvar/
        save /ousreal/
        save /ouschar/

!-------------------------------------------------------------
! end of header
!-------------------------------------------------------------

