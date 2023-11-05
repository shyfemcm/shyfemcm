
!--------------------------------------------------------------------------
!
!    Copyright (C) 2013,2019  Georg Umgiesser
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

! maxdimnh = maximum dimension of the hypothetic general matrix row that has to be solved for 
!            non hydrostatic contibution
! kpl(0:maxdimnh) = array with nodes assigned in the position on row of general matrix where 
!                   there is a contribution. 
! lpl(0:maxdimnh) = array with layers assigned in the position on row of general matrix where 
!                   there is a contribution. 
! rownh(0:maxdimnh) = array of mask values for row of general (1 if in the position on row of 
!                     general matrix there is a contribution, 0 otherwise
! kpl1(0:maxdimnh) = array with maximum dimension corresponding to the number of positions on 
!                    row of the general matrix where there is a contribution. The array stores 
!                    the corresponding nodes.
! lpl1(0:maxdimnh) = array with maximum dimension corresponding to the number of positions on row 
!                    of the general matrix where there is a contribution. The array stores the 
!                    corresponding layers.
! ipl(0:nkndim) = auxiliary array that stores for each node k, the sum of lmax for all the nodes 
!       till k. It help in the assignment of pointer position in the final general array for the
!       solution of the general matrix
! iop(nkndim,nlvdim,ngrdim) =  given the row corresponding to (k,l) of the general matrix, 
!                              for each nn that 
!                corresponds to the contribution of the node k1 at layer l, connected with k,
!                iop provides the relative position on the row, respect to the diagonal
! iopup(nkndim,nlvdim) = given the row corresponding to (k,l) of the general matrix, iopup provides
!              the sum of all the positions with non zero values on the previous rows of
!              the general matrix until row (k,l)
! kop(nkndim,nlvdim,ngrdim) = for each row (k,l), array with nodes k1 of column (k1,l) contributing
!                             at non zero values in general matrix
! lop(nkndim,nlvdim,ngrdim) = for each row (k,l), array with layer l of column (k1,l) contributing 
!                             at non zero values in general matrix
! iopp(neldim,9,nlvdim) = pointer that identifies the position in the general array that is built to
!                         solve the non hydrostatic contribution, given the considered element ie, 
!                         the 9 contribution of the combination of its nodes (kn(1:3)=nen3v(ie,1:3)
!                         kn(1:3)<-->kn(1:3)
! ioii((ngrdim3*maxdimnh) = node position k on row of the general matrix corresponding to the pointer
!                          iopp
! iojj((ngrdim3*maxdimnh) = node position k1 on column of the general matrix corresponding to the 
!	                   pointer iopp
! ioll((ngrdim3*maxdimnh) = layer position l on row of the general matrix corresponding to the 
!			   pointer iopp
! ioii1((ngrdim3*maxdimnh) = position on row of the general matrix corresponding to the pointer
!                          iopp
! iojj1((ngrdim3*maxdimnh) = position on the column of the general  matrix corresponding to the 
!	                   pointer iopp


! revision log :
!
! 12.09.2013	ggu	changed VERS_6_1_67
! 16.02.2019	ggu	changed VERS_7_5_60

	integer csrdimnh
	parameter ( csrdimnh = 9 * neldim * nlvdim )

	integer matdimmax
	common /matdimmax/matdimmax
	save /matdimmax/

	integer nnzeronh
	common /nnzeronh/nnzeronh
	save /nnzeronh/

	integer maxdimnh
	parameter( maxdimnh = nlvdim * nkndim ) 

	integer ngrdim3
	parameter (ngrdim3 = ngrdim + 3 )

	integer kpl(0:maxdimnh)  
	integer lpl(0:maxdimnh)
	integer rownh(0:maxdimnh)
	integer kpl1(0:maxdimnh) 
	integer lpl1(0:maxdimnh)
	integer ipl(0:nkndim)
	integer iop(ngrdim,nlvdim,nkndim) 
        integer iopup(0:nlvdim,nkndim)
	integer iopp(neldim,9,nlvdim)	
	integer ikop(ngrdim3,nlvdim,nkndim)
	integer ilop(ngrdim3,nlvdim,nkndim)
	integer ioii(ngrdim3*maxdimnh)
	integer iojj(ngrdim3*maxdimnh)
	integer ioii1(ngrdim3*maxdimnh)
	integer iojj1(ngrdim3*maxdimnh)
	integer ioll(ngrdim3*maxdimnh)
	double precision conh(ngrdim3*maxdimnh)
	double precision vecnot(ngrdim3*maxdimnh)

	common /conh/conh
	common /vecnot/vecnot
	save /conh/
	save /vecnot/

	common /iopp/iopp,/ioii/ioii,/iojj/iojj,/ioll/ioll
	common /ioii1/ioii1,/iojj1/iojj1

	save /iopp/,/ioii/,/iojj/,/ioll/
	save /ioii1/,/iojj1/


	!double precision rvecnh(maxdimnh) 
	real*8 rvecnh(maxdimnh) 
	common /rvecnh/rvecnh
	!double precision rauxnh(maxdimnh)!controllare dimensione
	real*8 rauxnh(maxdimnh)!controllare dimensione
	common /rauxnh/rauxnh

