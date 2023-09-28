
!--------------------------------------------------------------------------
!
!    Copyright (C) 2012,2014,2019  Georg Umgiesser
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
! header file for nos file format
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
! nositem	number of maximum entries in table
!
! nosvar	integer parameters of open files
! noschar	string parameters of open files
!
! nosvar(0,n)   iunit
! nosvar(1,n)   nvers
! nosvar(2,n)   nkn
! nosvar(3,n)   nel
! nosvar(4,n)   nlv
! nosvar(5,n)   nvar
! nosvar(6,n)   date
! nosvar(7,n)   time
!
! noschar(1,n)  title
! noschar(2,n)  femver
!
!-------------------------------------------------------------
! parameters
!-------------------------------------------------------------

! revision log :
!
! 05.11.2012	ggu	changed VERS_6_1_60
! 28.01.2014	ggu	changed VERS_6_1_71
! 16.02.2019	ggu	changed VERS_7_5_60
! 31.10.2019	ggu	changed VERS_7_5_65

!-------------------------------------------------------------


        integer ftype,maxvers,maxcomp
        parameter(ftype=161,maxvers=5,maxcomp=3)

        integer ndim,nitdim,nchdim
        parameter(ndim=100,nitdim=7,nchdim=2)

!-------------------------------------------------------------
! common
!-------------------------------------------------------------

        integer nositem
        common /nositm/nositem

        integer nosvar(0:nitdim,ndim)
        common /nosvar/nosvar

        character*80 noschar(nchdim,ndim)
        common /noschar/noschar

!-------------------------------------------------------------
! save
!-------------------------------------------------------------

        save /nositm/
        save /nosvar/
        save /noschar/

!-------------------------------------------------------------
! end of header
!-------------------------------------------------------------

