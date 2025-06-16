
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015,2017-2019  Georg Umgiesser
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

! mpi routines - partition on nodes (nodes are unique)
!
! contents :
!
! revision log :
!
! 24.11.2015	ggu	project started
! 07.12.2017	ggu	changed VERS_7_5_40
! 24.01.2018	ggu	changed VERS_7_5_41
! 22.02.2018	ggu	changed VERS_7_5_42
! 04.04.2018	ggu	new routine shympi_exit
! 04.04.2018	ggu	bug fix on scalar reduction (argument was changed)
! 10.04.2018	ggu	code to exchange arrays
! 19.04.2018	ggu	changed VERS_7_5_45
! 26.04.2018	ggu	changed VERS_7_5_46
! 11.05.2018	ggu	changes in global variables and exchange_arrays()
! 06.07.2018	ggu	changed VERS_7_5_48
! 16.02.2019	ggu	changed VERS_7_5_60
! 07.06.2020	ggu	new routines, 3d exchange array still missing
! 09.04.2021	clr	bug fix in shympi_bcast_array_r() -> real arg
! 17.04.2021	ggu	new shympi_exchange_array_3(), check_external_numbers()
! 22.04.2021	ggu	allocation of some arrays for bounds check
! 22.04.2021	ggu	bug fix in shympi_copy_*()
! 23.04.2021    clr     changes in MODEL PROCEDURE decl. for meson compatibility
! 25.06.2021    ggu     in shympi_init() check if basin has been read
! 10.11.2021    ggu     error fix mixing scalar and array arguments
! 22.11.2021    ggu     error fix when reading grd file
! 01.04.2022    ggu     new routine shympi_set_debug()
! 02.04.2022    ggu     in shympi_check_array() better error output
! 02.04.2022    ggu     new array nlv_domains containing nlv values for domains
! 02.04.2022    ggu     routines shympi_gather_array_3d_*() finally running
! 03.04.2022    ggu     new routine shympi_bcast_array_d()
! 05.04.2022    ggu     new routines to copy from global to local
! 06.04.2022    ggu     new routines for handling double precision
! 10.04.2022    ggu     in shympi_check_array() pass info on h/v dim
! 10.04.2022    ggu     better error reporting (not finished)
! 12.04.2022    ggu     file cleaned
! 01.06.2022    ggu     new routine shympi_gather_root()
! 09.10.2022    ggu     new variable nlv_local
! 11.10.2022    ggu     new routines to deal with fixed first dimension
! 13.10.2022    ggu     bug fix in shympi_g2l_array_3d - wrong indices
! 13.10.2022    ggu     new routine shympi_g2l_array_3d_d()
! 16.10.2022    ggu     shympi_exchange_array_3() eliminated
! 11.11.2022    ggu     in shympi_collect_node_value_3d_r() only copy existing
! 18.03.2023    ggu     id_elem is now (0:3)
! 27.03.2023    ggu     new routines shympi_receive(), more docs
! 27.03.2023    ggu     new shympi_l2g_array_fix_i, shympi_gather_array_fix_i
! 13.04.2023    ggu     flagged potential problems with GGU_NKN_NEL_BUG
! 13.04.2023    ggu     introduced bnode, belem (distinguish calls to node/elem)
! 19.04.2023    ggu     potential bugs fixed with passing belem
! 03.05.2023    ggu     new routine shympi_bdebug()
! 18.05.2023    ggu     in shympi_gather_array() eliminate rectify array
! 09.06.2023    ggu     new routine error_stop()
! 07.03.2024    ggu     double routines for shympi_l2g_array_fix_d/2d_d
! 12.03.2024    ggu     double routines for shympi_g2l_array_fix_ird/2d_d
! 13.03.2024    ggu     bmpi_skip introduced
! 05.04.2024    ggu     shympi_check_3d_elem_r(), nghost_max_global
! 28.08.2024    ggu     new variable bmpi_debug_txt
! 06.09.2024    ggu     better debug info for shympi_check_array()
! 06.09.2024    lrp     nuopc-compliant
! 13.10.2024    ggu     new parameter bextra_exchange to deal with INTEL_BUG
! 13.11.2024    ggu     new routine shympi_find_element()
! 23.11.2024    ggu     new variable pquality
! 04.12.2024    ggu     implement bmpi_allgather for gather routines
! 07.12.2024    ggu     big changes: make sure arrays are same length for reduce
! 13.12.2024    ggu     error_stop() implemented
! 28.01.2025    ggu     compatibility for gfortran 4 inserted
! 06.02.2025    ggu     new routines shympi_exchange_fix_node/elem()

!******************************************************************

! for partioning on nodes the following is true (at least 2 domains)
!
! nkn_global > nkn_local >  nkn_unique == nkn_inner
! nel_global > nel_local >= nel_unique >= nel_inner

!==================================================================
        module shympi
!==================================================================

	implicit none

	public

        logical, parameter :: bcomment = .false. 	!write comment
        integer, parameter :: iuc = 999          	!unit for comment
        character*2, parameter :: indent = '  '		!indent for comment

	logical, parameter :: blocal_shympi_debug = .false. !write debug
	logical, save :: bmpi_debug = .false.		!writes debug messages
	logical, save :: bmpi_debug_txt = .false.	!writes mpi_debug_*.txt
	logical, save :: bmpi_ldebug = .false.		!local debug
	logical, save :: bmpi_unit = .false.		!write debug to my_unit
	integer, save :: my_unit_base = 800		!unit base for output

	logical, save :: bmpi_skip = .true.		!skip if not bmpi
	logical, save :: bextra_exchange = .false.	!deal with INTEL_BUG

	integer, save :: i_code = 0
	integer, save :: i_code_error = 33
	integer, save :: i_code_success = 0

!	---------------------------------
!	next variables are set internally
!	---------------------------------

	logical, save :: bmpi = .false.
	logical, save :: bmpi_master = .false.
	logical, save :: bmpi_support = .true.
	logical, save :: bmpi_alloper = .true.	!do allgather and allreduce

	integer, save :: n_threads = 1
	integer, save :: my_id = 0
	integer, save :: my_unit = 0

	integer, save :: status_size = 0

	integer, save :: ngr_global = 0		!ngr of total basin
	integer, save :: nlv_global = 0		!nlv of total basin

	integer, save :: nlv_local = 0		!nlv of this partition

	integer, save :: nkn_global = 0		!total basin
	integer, save :: nel_global = 0
	integer, save :: nkn_local = 0		!this domain
	integer, save :: nel_local = 0
	integer, save :: nkn_unique = 0		!this domain unique
	integer, save :: nel_unique = 0
	integer, save :: nkn_inner = 0		!only proper, no ghost
	integer, save :: nel_inner = 0

	integer, save :: nk_max = 0		!max of nkn of all domains
	integer, save :: ne_max = 0		!max of nel of all domains
	integer, save :: nn_max = 0		!max of nkn/nel of all domains

! total number of nodes/elems in domains

	integer,save,pointer :: n_domains(:)
	integer,save,target,allocatable :: nkn_domains(:)	!local total
	integer,save,target,allocatable :: nel_domains(:)
	integer,save,target,allocatable :: nkn_domains_u(:)	!local unique
	integer,save,target,allocatable :: nel_domains_u(:)
	integer,save,allocatable :: nkn_cum_domains(:)		!cumulative
	integer,save,allocatable :: nel_cum_domains(:)
	integer,save,allocatable :: nlv_domains(:)

! information on ghost areas (nodes and elements)

	integer, save :: n_ghost_areas = 0
	integer, save :: n_ghost_nodes_max = 0
	integer, save :: n_ghost_elems_max = 0
	integer, save :: n_ghost_max = 0
	integer, save :: n_ghost_max_global = 0
	integer, save :: n_buffer = 0

! arrays for ghost nodes/elems

	integer,save,allocatable :: ghost_areas(:,:)
	integer,save,allocatable :: ghost_nodes_in(:,:)
	integer,save,allocatable :: ghost_nodes_out(:,:)
	integer,save,allocatable :: ghost_elems_in(:,:)
	integer,save,allocatable :: ghost_elems_out(:,:)

! id of nodes and elems
! the id of a node is unique
! a node only belongs to one domain (partitioning on nodes)
! an element can belong to more than one domain (maximum 3)
! the dimension of id_elem is id_elem(0:3,nel)
! id_elem(0,ie) is the the total number of domains the element ie belongs to
! id_elem(1:3,ie) are the domains the element ie belongs to
! in id_elem(1,ie) is the principal domain

	integer,save,allocatable :: id_node(:)		!domain (id) of node
	integer,save,allocatable :: id_elem(:,:)	!domain (id) of elem

	! sorted node/elem numbering as in total domain

	integer,save,allocatable :: ip_sort_node(:)	!sorted internal nodes
	integer,save,allocatable :: ip_sort_elem(:)	!sorted internal elems

	! next are global arrays for external node/elem numbers
	! ip_ext_node(nkn_global)
	! ip_ext_elem(nel_global)

	integer,pointer :: ip_ext(:) !pointer to  external nums
	integer,save,target,allocatable :: ip_ext_node(:) !global external nums
	integer,save,target,allocatable :: ip_ext_elem(:)

	! next are pointers from local to global internal node/elem numbers

	integer,pointer :: ip_int(:) !pointer to  internal nums
	integer,save,target,allocatable :: ip_int_node(:) !global internal nums
	integer,save,target,allocatable :: ip_int_elem(:)

	! next are pointers from local to global internal node/elem numbers

	integer,save,target,allocatable :: ip_int_nodes(:,:) !global int nums
	integer,save,target,allocatable :: ip_int_elems(:,:)

	! nen3v and hlv index of global domain

	integer,save,allocatable :: nen3v_global(:,:)	!global element index
	real,save,allocatable :: hlv_global(:)		!global layer depths

	! quality index of partition

	real, save :: pquality = 0.

	! communication structures

        type communication_info
          integer, public :: numberID
          integer, public :: totalID
          !contains the global ID of the elements belonging to the process
          integer, public, dimension(:), allocatable :: globalID,localID
          !contains the rank ID of the process to which belong the global ID
          !elements with which to communicate
          integer, public, dimension(:,:), allocatable :: rankID
          integer, public, dimension(:,:), allocatable :: neighborID
        end type communication_info

	type (communication_info), SAVE :: univocal_nodes
        integer, allocatable, save, dimension(:) :: allPartAssign

!-------------------------------------------------------
!	receive from given areas
!-------------------------------------------------------

        INTERFACE shympi_receive
        MODULE PROCEDURE   shympi_receive_i &
     &                   , shympi_receive_r
        END INTERFACE

!-------------------------------------------------------
!	exchange ghost node and element information
!-------------------------------------------------------

        INTERFACE shympi_exchange_fix_node
        MODULE PROCEDURE  shympi_exchange_fix_node_i
        END INTERFACE

        INTERFACE shympi_exchange_fix_elem
        MODULE PROCEDURE  shympi_exchange_fix_elem_i
        END INTERFACE

        INTERFACE shympi_exchange_3d_node
        MODULE PROCEDURE  shympi_exchange_3d_node_r &
     &                   ,shympi_exchange_3d_node_d &
     &                   ,shympi_exchange_3d_node_i
        END INTERFACE

        INTERFACE shympi_exchange_3d0_node
       	MODULE PROCEDURE   shympi_exchange_3d0_node_r
!     +                   ,shympi_exchange_3d0_node_d
!     +                   ,shympi_exchange_3d0_node_i
        END INTERFACE

        INTERFACE shympi_exchange_2d_node
        MODULE PROCEDURE  shympi_exchange_2d_node_r &
     &                   ,shympi_exchange_2d_node_d &
     &                   ,shympi_exchange_2d_node_i
        END INTERFACE

        INTERFACE shympi_exchange_3d_elem
        MODULE PROCEDURE   shympi_exchange_3d_elem_r
!     +                   ,shympi_exchange_3d_elem_d
!     +                   ,shympi_exchange_3d_elem_i
        END INTERFACE

        INTERFACE shympi_exchange_2d_elem
        MODULE PROCEDURE  shympi_exchange_2d_elem_r &
     &                   ,shympi_exchange_2d_elem_d &
     &                   ,shympi_exchange_2d_elem_i
        END INTERFACE

!-------------------------------------------------------
!	exchange ghost information and checks if equal
!-------------------------------------------------------

        INTERFACE shympi_check_elem
        MODULE PROCEDURE  shympi_check_2d_elem_r &
     &                   ,shympi_check_2d_elem_d &
     &                   ,shympi_check_2d_elem_i &
     &			 ,shympi_check_3d_elem_r
!     +                   ,shympi_check_3d_elem_d
!     +                   ,shympi_check_3d_elem_i
        END INTERFACE

        INTERFACE shympi_check_node
        MODULE PROCEDURE  shympi_check_2d_node_r &
     &                   ,shympi_check_2d_node_d &
     &                   ,shympi_check_2d_node_i &
     &			 ,shympi_check_3d_node_r
!     +                   ,shympi_check_3d_node_d
!     +                   ,shympi_check_3d_node_i
        END INTERFACE

        INTERFACE shympi_check_2d_node
        MODULE PROCEDURE  shympi_check_2d_node_r &
     &                   ,shympi_check_2d_node_d &
     &                   ,shympi_check_2d_node_i
        END INTERFACE

        INTERFACE shympi_check_2d_elem
        MODULE PROCEDURE  shympi_check_2d_elem_r &
     &                   ,shympi_check_2d_elem_d &
     &                   ,shympi_check_2d_elem_i
        END INTERFACE

        INTERFACE shympi_check_3d_node
        MODULE PROCEDURE  shympi_check_3d_node_r
!     +                   ,shympi_check_3d_node_d
!     +                   ,shympi_check_3d_node_i
        END INTERFACE

        INTERFACE shympi_check_3d0_node
        MODULE PROCEDURE  shympi_check_3d0_node_r
!     +                   ,shympi_check_3d0_node_d
!     +                   ,shympi_check_3d0_node_i
        END INTERFACE

        INTERFACE shympi_check_3d_elem
        MODULE PROCEDURE  shympi_check_3d_elem_r
!     +                   ,shympi_check_3d_elem_d
!     +                   ,shympi_check_3d_elem_i
        END INTERFACE

        INTERFACE shympi_check_array
        MODULE PROCEDURE  shympi_check_array_i &
     &                   ,shympi_check_array_r &
     &                   ,shympi_check_array_d
        END INTERFACE

!-------------------------------------------------------
!	gathers information from domains into one common global array
!-------------------------------------------------------

        INTERFACE shympi_gather
        MODULE PROCEDURE  shympi_gather_scalar_i &
     &                   ,shympi_gather_array_2d_i &
     &                   ,shympi_gather_array_2d_r &
     &                   ,shympi_gather_array_2d_d &
     &                   ,shympi_gather_array_3d_i &
     &                   ,shympi_gather_array_3d_r &
     &                   ,shympi_gather_array_3d_d &
     &			 ,shympi_gather_array_fix_i &
     &			 ,shympi_gather_array_fix_r &
     &			 ,shympi_gather_array_fix_d
        END INTERFACE

        INTERFACE shympi_gather_root
        MODULE PROCEDURE  shympi_gather_root_array_2d_d
        END INTERFACE

!-------------------------------------------------------
!	gathers information from domains and sums it back
!-------------------------------------------------------

        INTERFACE shympi_gather_and_sum
        MODULE PROCEDURE  shympi_gather_and_sum_i &
     &                   ,shympi_gather_and_sum_r &
     &                   ,shympi_gather_and_sum_d
        END INTERFACE

!-------------------------------------------------------
!	broadcats information
!-------------------------------------------------------

        INTERFACE shympi_bcast
        MODULE PROCEDURE  shympi_bcast_scalar_i &
     &                   ,shympi_bcast_array_r &
     &                   ,shympi_bcast_array_d
        END INTERFACE

!-------------------------------------------------------
!	collects nodal values from domains
!-------------------------------------------------------

        INTERFACE shympi_collect_node_value
        MODULE PROCEDURE   shympi_collect_node_value_2d_i &
     &                    ,shympi_collect_node_value_2d_r &
     &                    ,shympi_collect_node_value_3d_r
!     +                    ,shympi_collect_node_value_3d_i
        END INTERFACE

!-------------------------------------------------------
!	general and special reduce routines
!-------------------------------------------------------

        INTERFACE shympi_min
        MODULE PROCEDURE   shympi_min_r &
     &			  ,shympi_min_i &
     &			  ,shympi_min_d &
     &			  ,shympi_min_0_r &
     &			  ,shympi_min_0_i &
     &			  ,shympi_min_0_d
        END INTERFACE

        INTERFACE shympi_max
        MODULE PROCEDURE   shympi_max_r &
     &			  ,shympi_max_i &
     &			  ,shympi_max_d &
     &			  ,shympi_max_0_r        &
     &			  ,shympi_max_0_i	 &
     &			  ,shympi_max_0_d
        END INTERFACE

        INTERFACE shympi_sum
        MODULE PROCEDURE   shympi_sum_r &
     &			  ,shympi_sum_i &
     &			  ,shympi_sum_d &
     &			  ,shympi_sum_0_r &
     &			  ,shympi_sum_0_i &
     &			  ,shympi_sum_0_d
        END INTERFACE

        INTERFACE shympi_array_reduce
        MODULE PROCEDURE   shympi_array_reduce_i &
     &			  ,shympi_array_reduce_r &
     &			  ,shympi_array_reduce_d
        END INTERFACE

        INTERFACE shympi_reduce
        MODULE PROCEDURE shympi_reduce_r
!     +                   ,shympi_reduce_i
        END INTERFACE

!-------------------------------------------------------
!	copies local to global arrays and viceversa
!-------------------------------------------------------

        INTERFACE shympi_l2g_array
        MODULE PROCEDURE   shympi_l2g_array_2d_r &
     &			  ,shympi_l2g_array_2d_i &
     &			  ,shympi_l2g_array_2d_d &
     &			  ,shympi_l2g_array_3d_r &
     &			  ,shympi_l2g_array_3d_i &
     &			  ,shympi_l2g_array_3d_d &
     &			  ,shympi_l2g_array_fix_i &
     &			  ,shympi_l2g_array_fix_r &
     &			  ,shympi_l2g_array_fix_d
        END INTERFACE

        INTERFACE shympi_g2l_array
        MODULE PROCEDURE   shympi_g2l_array_2d_r &
     &			  ,shympi_g2l_array_2d_i &
     &			  ,shympi_g2l_array_2d_d &
     &			  ,shympi_g2l_array_3d_r &
     &			  ,shympi_g2l_array_3d_i &
     &			  ,shympi_g2l_array_3d_d &
     &			  ,shympi_g2l_array_fix_i &
     &			  ,shympi_g2l_array_fix_r &
     &			  ,shympi_g2l_array_fix_d
        END INTERFACE

!-------------------------------------------------------
!	gets array and value (function unclear, not used)
!-------------------------------------------------------

        INTERFACE shympi_get_array
        MODULE PROCEDURE   shympi_get_array_2d_r &
     &			  ,shympi_get_array_2d_i
!     +			  ,shympi_get_array_3d_r
!     +			  ,shympi_get_array_3d_i
        END INTERFACE

        INTERFACE shympi_getvals
        MODULE PROCEDURE   shympi_getvals_2d_node_r &
     &			  ,shympi_getvals_2d_node_i &
     &			  ,shympi_getvals_3d_node_r &
     &			  ,shympi_getvals_3d_node_i
        END INTERFACE

!-------------------------------------------------------
!	routines for element partition (are empty)
!-------------------------------------------------------

        INTERFACE shympi_exchange_and_sum_3d_nodes
        MODULE PROCEDURE   shympi_exchange_and_sum_3d_nodes_r &
     &			  ,shympi_exchange_and_sum_3d_nodes_d
        END INTERFACE

        INTERFACE shympi_exchange_and_sum_2d_nodes
        MODULE PROCEDURE   shympi_exchange_and_sum_2d_nodes_r &
     &			  ,shympi_exchange_and_sum_2d_nodes_d
        END INTERFACE

        INTERFACE shympi_exchange_2d_nodes_min
        MODULE PROCEDURE  shympi_exchange_2d_nodes_min_i &
     &			  ,shympi_exchange_2d_nodes_min_r
        END INTERFACE

        INTERFACE shympi_exchange_2d_nodes_max
        MODULE PROCEDURE   shympi_exchange_2d_nodes_max_i &
     &			  ,shympi_exchange_2d_nodes_max_r
        END INTERFACE

!-------------------------------------------------------
!	error handling
!-------------------------------------------------------

        INTERFACE error_stop
        MODULE PROCEDURE  error_stop_2 &
     &			, error_stop_1 &
     &			, error_stop_0 &
     &			, error_stop_i0
        END INTERFACE

        INTERFACE shympi_bdebug
        MODULE PROCEDURE  shympi_bdebug_0 &
     &			, shympi_bdebug_3 
        END INTERFACE

!-------------------------------------------------------
!	for error reporting
!-------------------------------------------------------

	character*80 :: textd
	character*80, save :: textk,texte,text2

!==================================================================
        contains
!==================================================================

	subroutine shympi_init(b_want_mpi, b_want_mpi_init)

! this is the first call to shympi
!
! it should be called right after having read the basin
! after this shympi_setup should be called
!
! all general data is set up
! data is allocated
!
! still no partitioning of the basin has been done

	use basin
	use levels

	logical b_want_mpi
	logical, optional :: b_want_mpi_init

	logical bstop,binit
	integer ierr,size
	character*10 cunit
	character*80 file

	!-----------------------------------------------------
	! initializing
	!-----------------------------------------------------

	!in some cases you let other initialize MPI for SHYFEM
	!an intialization flag is passed from the interface
	!the default value is true: always initialize MPI

	if (present(b_want_mpi_init)) then
	  binit = b_want_mpi_init
	else
	  binit = .true.
	endif

	call shympi_init_internal(my_id,n_threads,binit)

        if( .not. basin_has_read_basin() ) then
          write(6,*) 'grd file has been read: ',nkn,nel,ngr
          if( nkn == 0 ) then
            stop 'error stop shympi_init: ' // &
     &			'basin has not been initialized'
          end if
        end if

        if( nkn == 0 .or. nel == 0 ) then
         if( shympi_is_master() ) then
          write(6,*) 'nkn = ',nkn,'  nel = ',nel
          write(6,*) '*** basin contains no nodes or elements...'
         end if
         call shympi_stop('error stop shympi_init')
        end if

	bmpi = n_threads > 1
	bmpi_master = my_id == 0

	bmpi_skip = bmpi_skip .and. .not. bmpi

	bstop = .false.
	if( shympi_is_master() ) then
!	 if( b_want_mpi .and. .not. bmpi ) then
!	  write(6,*) 'program wants mpi but only one thread available'
!	  bstop = .true.
!	 end if
	 if( .not. b_want_mpi .and. bmpi ) then
	  write(6,*) 'program does not want mpi '// &
     &			'but program running in mpi mode'
	  bstop = .true.
	 end if
	end if
	if( bstop ) stop 'error stop shympi_init'

	ngr_global = ngr

	nkn_global = nkn
	nel_global = nel
	nkn_local = nkn
	nel_local = nel
	nkn_inner = nkn
	nel_inner = nel
	nkn_unique = nkn
	nel_unique = nel

	!-----------------------------------------------------
	! allocate important arrays
	!-----------------------------------------------------

	call shympi_get_status_size_internal(size)
	status_size = size

        allocate(nkn_domains(n_threads))
        allocate(nel_domains(n_threads))
        allocate(nkn_domains_u(n_threads))
        allocate(nel_domains_u(n_threads))
        allocate(nkn_cum_domains(0:n_threads))
        allocate(nel_cum_domains(0:n_threads))

	nkn_domains(1) = nkn
	nel_domains(1) = nel
	nkn_domains_u(1) = nkn
	nel_domains_u(1) = nel
        nk_max = nkn
        ne_max = nel
        nn_max = max(nkn,nel)

	nkn_cum_domains(0) = 0
	nkn_cum_domains(1) = nkn
	nel_cum_domains(0) = 0
	nel_cum_domains(1) = nel

	call shympi_alloc_global(nkn,nel,nen3v,ipv,ipev)
	call levels_init_2d(nkn,nel)	!needed for bounds check

	!-----------------------------------------------------
	! next is needed if program is not running in mpi mode
	!-----------------------------------------------------

	if( .not. bmpi ) then
	  call shympi_alloc_id(nkn,nel)
          call shympi_alloc_sort(nkn,nel)
          call mpi_sort_index(nkn,nel)
	  call shympi_alloc_ghost(1)	!needed for bounds check
        end if

	!-----------------------------------------------------
	! output to terminal
	!-----------------------------------------------------

	if( bmpi ) then
	  if( bmpi_debug_txt ) then
	    write(cunit,'(i10)') my_id
	    cunit = adjustl(cunit)
	    file = 'mpi_debug_' // trim(cunit) // '.txt'
	    !call shympi_get_new_unit(my_unit)
	    my_unit = my_unit_base + my_id
	    open(unit=my_unit,file=file,status='unknown')
	    write(my_unit,*) '=========================================='
	    write(my_unit,*) 'this is a debug file for domain ',my_id
	    write(my_unit,*) 'set bmpi_debug_txt=.false. to avoid this output'
	    write(my_unit,*) '=========================================='
	    write(my_unit,*) 'shympi initialized: ',my_id,n_threads,my_unit
	  else
	    if( shympi_is_master() ) then
	      write(6,*) 'no mpi_debug_*.txt file opened: ',my_unit
	      write(6,*) 'set bmpi_debug_txt=.true. to write this file'
	    end if
	  end if
	  write(6,*) 'shympi initialized: ',my_id,n_threads,my_unit
	else if( bmpi_debug ) then
	  write(6,*) 'shympi initialized: ',my_id,n_threads
	  write(6,*) 'shympi is not running in mpi mode'
	end if

	flush(6)

	!-----------------------------------------------------
	! finish up
	!-----------------------------------------------------

	call shympi_barrier_internal
	call shympi_syncronize_initial
	call shympi_syncronize_internal

	!-----------------------------------------------------
	! end of routine
	!-----------------------------------------------------

	end subroutine shympi_init

!******************************************************************

	subroutine shympi_alloc_id(nk,ne)

	integer nk,ne	!local domain

	!write(6,*) 'shympi_alloc_id: ',nk,ne

	allocate(id_node(nk))
	allocate(id_elem(0:3,ne))

	id_node = my_id
        id_elem = -1
        id_elem(0,:) = 1		! element just in one domain
        id_elem(1,:) = my_id		! element is in domain my_id

	end subroutine shympi_alloc_id

!******************************************************************

	subroutine shympi_alloc_sort(nk,ne)

	integer nk,ne

        if( allocated(ip_sort_node) ) deallocate(ip_sort_node)
        if( allocated(ip_sort_elem) ) deallocate(ip_sort_elem)

	allocate(ip_sort_node(nk))
	allocate(ip_sort_elem(ne))

	ip_sort_node = 0
	ip_sort_elem = 0

	end subroutine shympi_alloc_sort

!******************************************************************

	subroutine shympi_alloc_global(nk,ne,nen3v,ipv,ipev)

! these are global arrays

	integer nk,ne		!global nkn, nel
	integer nen3v(3,ne)
	integer ipv(nk)
	integer ipev(ne)

	integer i

	allocate(ip_ext_node(nk))
	allocate(ip_ext_elem(ne))
	allocate(ip_int_node(nk))
	allocate(ip_int_elem(ne))
	allocate(ip_int_nodes(nk,1))
	allocate(ip_int_elems(ne,1))
	allocate(nen3v_global(3,ne))

	do i=1,nk
	  ip_int_node(i) = i
	  ip_int_nodes(i,1)= i
	end do

	do i=1,ne
	  ip_int_elem(i) = i
	  ip_int_elems(i,1)= i
	end do

	nen3v_global = nen3v
	ip_ext_node = ipv
	ip_ext_elem = ipev

	end subroutine shympi_alloc_global

!******************************************************************

	subroutine shympi_alloc_ghost(n)

	use basin

	integer n

	if( allocated(ghost_areas) ) then
	  deallocate(ghost_areas)
	  deallocate(ghost_nodes_out,ghost_nodes_in)
	  deallocate(ghost_elems_out,ghost_elems_in)
	end if

	allocate(ghost_areas(5,n_ghost_areas))
        allocate(ghost_nodes_out(n,n_ghost_areas))
        allocate(ghost_nodes_in(n,n_ghost_areas))
        allocate(ghost_elems_out(n,n_ghost_areas))
        allocate(ghost_elems_in(n,n_ghost_areas))

	ghost_areas = 0
        ghost_nodes_out = 0
        ghost_nodes_in = 0
        ghost_elems_out = 0
        ghost_elems_in = 0

	end subroutine shympi_alloc_ghost

!******************************************************************

        subroutine shympi_alloc_buffer(n)

        integer n

        if( n_buffer >= n ) return

!        if( n_buffer > 0 ) then
!          deallocate(i_buffer_in)
!          deallocate(i_buffer_out)
!          deallocate(r_buffer_in)
!          deallocate(r_buffer_out)
!        end if

        n_buffer = n

!        allocate(i_buffer_in(n_buffer,n_ghost_areas))
!        allocate(i_buffer_out(n_buffer,n_ghost_areas))
!
!        allocate(r_buffer_in(n_buffer,n_ghost_areas))
!        allocate(r_buffer_out(n_buffer,n_ghost_areas))

        end subroutine shympi_alloc_buffer

!******************************************************************
!******************************************************************
!******************************************************************

	function shympi_internal_node(k)

	integer shympi_internal_node
	integer k
	integer i

	do i=1,nkn_global
	  if( ip_ext_node(i) == k ) exit
	end do
	if( i > nkn_global ) i = 0

	shympi_internal_node = i

	end function shympi_internal_node

!******************************************************************

	function shympi_internal_elem(ie)

	integer shympi_internal_elem
	integer ie
	integer i

	do i=1,nel_global
	  if( ip_ext_elem(i) == ie ) exit
	end do
	if( i > nel_global ) i = 0

	shympi_internal_elem = i

	end function shympi_internal_elem

!******************************************************************

	function shympi_exist_node(k)

	logical shympi_exist_node
	integer k

	shympi_exist_node = ( shympi_internal_node(k) > 0 )

	end function shympi_exist_node

!******************************************************************

	function shympi_exist_elem(ie)

	logical shympi_exist_elem
	integer ie

	shympi_exist_elem = ( shympi_internal_elem(ie) > 0 )

	end function shympi_exist_elem

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_set_hlv(nlv,hlv)

! sets nlv and hlv globally

        implicit none

        integer nlv
        real hlv(nlv)

        integer ia
        real, parameter :: flag = -1.35472E+10
        real, allocatable :: hlvs(:,:)

        allocate(nlv_domains(n_threads))
	call shympi_gather(nlv,nlv_domains)
        nlv_global = shympi_max(nlv)
        nlv_local = nlv

        allocate(hlv_global(nlv_global))
        allocate(hlvs(nlv_global,n_threads))

        hlv_global = flag
        hlv_global(1:nlv) = hlv
        call shympi_gather(hlv_global,hlvs)

        do ia=1,n_threads
          if( hlvs(nlv_global,ia) /= flag ) exit
        end do
        if( ia > n_threads ) then
          write(6,*) 'error setting global nlv'
          write(6,*) nlv_global,nlv
          write(6,*) hlvs
          stop 'error stop shympi_set_hlv: global nlv'
        end if

        hlv_global = hlvs(:,ia)

        end subroutine shympi_set_hlv

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_get_new_unit(iunit)

	integer iunit

	integer iu,iumin,iumax,ios
        logical opened

	iumin = 20
	iumax = 1000

	do iu=iumin,iumax
          inquire (unit=iu, opened=opened, iostat=ios)
          if (ios.ne.0) cycle
          if (.not.opened) exit
	end do

	if( iu > iumax ) then
	  iu = 0
	  stop 'error stop shympi_get_new_unit: no new unit'
	end if

	iunit = iu

	end subroutine shympi_get_new_unit

!******************************************************************
!******************************************************************
!******************************************************************

	function shympi_partition_on_elements()

	logical shympi_partition_on_elements

	shympi_partition_on_elements = .false.

	end function shympi_partition_on_elements

!******************************************************************

        function shympi_partition_on_nodes()

        logical shympi_partition_on_nodes

        shympi_partition_on_nodes = .true.

        end function shympi_partition_on_nodes

!******************************************************************
!******************************************************************
!******************************************************************

        function shympi_is_inner_node(k)

        logical shympi_is_inner_node
        integer k

        shympi_is_inner_node = ( k <= nkn_inner )

        end function shympi_is_inner_node

!******************************************************************

        function shympi_is_inner_elem(ie)

        logical shympi_is_inner_elem
        integer ie

        shympi_is_inner_elem = ( ie <= nel_inner )

        end function shympi_is_inner_elem

!******************************************************************

        function shympi_is_unique_node(k)

        logical shympi_is_unique_node
        integer k

        shympi_is_unique_node = ( k <= nkn_unique )

        end function shympi_is_unique_node

!******************************************************************

        function shympi_is_unique_elem(ie)

        logical shympi_is_unique_elem
        integer ie

        shympi_is_unique_elem = ( ie <= nel_unique )

        end function shympi_is_unique_elem

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_barrier

	call shympi_barrier_internal

	end subroutine shympi_barrier

!******************************************************************

	subroutine shympi_syncronize

	call shympi_syncronize_internal

	end subroutine shympi_syncronize

!******************************************************************

	subroutine shympi_finalize

! this is called from all processes to finish gracefully

	call shympi_barrier_internal
	call shympi_finalize_internal
	stop 'finalizing...'

	end subroutine shympi_finalize

!******************************************************************

	subroutine shympi_exit(ierr)

! this is called on error - do not call barrier

	integer ierr

	!call shympi_barrier_internal
	call shympi_finalize_internal
	call exit(ierr)
	stop 'exiting...'

	end subroutine shympi_exit

!******************************************************************

	subroutine shympi_stop(text)

! this is called as error stop - do not call barrier

	character*(*) text

	if( shympi_is_master() ) then
	  write(6,*) 'error stop ',trim(text)
	end if
	!call shympi_barrier_internal
	call shympi_abort_internal(i_code_error)
	stop 'aborting...'

	end subroutine shympi_stop

!******************************************************************

	subroutine shympi_abort

	call shympi_abort_internal(i_code)
	stop 'error stop: abort'

	end subroutine shympi_abort

!******************************************************************

	function shympi_wtime()

	double precision shympi_wtime
	double precision shympi_wtime_internal

	shympi_wtime = shympi_wtime_internal()

	end function shympi_wtime

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_receive_i(id_from,id_to &
     &					,n,val_in,val_out)

	integer id_from,id_to,n
        integer val_in(n)
        integer val_out(n)

	call shympi_receive_internal_i(id_from,id_to &
     &					,n,val_in,val_out)

	end subroutine shympi_receive_i

!*******************************

	subroutine shympi_receive_r(id_from,id_to &
     &					,n,val_in,val_out)

	integer id_from,id_to,n
        real val_in(n)
        real val_out(n)

	call shympi_receive_internal_r(id_from,id_to &
     &					,n,val_in,val_out)

	end subroutine shympi_receive_r

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_fix_node_i(nfix,val)

	use basin

	integer nfix
	integer val(:,:)

	logical, parameter :: belem = .false.
	integer ni1,ni2
	integer ilhkv(nkn)

	if( bmpi_skip ) return

	ni1 = size(val,1)
	ni2 = size(val,2)

	if( ni1 /= nfix .or. ni2 /= nkn ) then
	  write(6,*) nfix,nkn,ni1,ni2
	  stop 'error stop shympi_exchange_fix_node_i: incomp dim'
	end if

	ilhkv = nfix
	call shympi_exchange_internal_i(belem,1,nfix,nkn,ilhkv &
     &			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_fix_node_i

!*******************************

	subroutine shympi_exchange_fix_elem_i(nfix,val)

	use basin

	integer nfix
	integer val(:,:)

	logical, parameter :: belem = .true.
	integer ni1,ni2
	integer ilhv(nel)

	if( bmpi_skip ) return

	ni1 = size(val,1)
	ni2 = size(val,2)

	if( ni1 /= nfix .or. ni2 /= nel ) then
	  write(6,*) nfix,nel,ni1,ni2
	  stop 'error stop shympi_exchange_fix_elem_i: incomp dim'
	end if

	ilhv = nfix
	call shympi_exchange_internal_i(belem,1,nfix,nel,ilhv &
     &			,ghost_elems_in,ghost_elems_out,val)

	end subroutine shympi_exchange_fix_elem_i

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_3d_node_i(val)

	use basin
	use levels

	integer val(nlvdi,nkn)
	logical, parameter :: belem = .false.

	if( bmpi_skip ) return

	call shympi_exchange_internal_i(belem,1,nlvdi,nkn,ilhkv &
     &			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_3d_node_i

!*******************************

	subroutine shympi_exchange_3d_node_r(val)

	use basin
	use levels

	real val(nlvdi,nkn)
	logical, parameter :: belem = .false.

	if( bmpi_skip ) return

	call shympi_exchange_internal_r(belem,1,nlvdi,nkn,ilhkv &
     &			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_3d_node_r

!*******************************

	subroutine shympi_exchange_3d_node_d(val)

	use basin
	use levels

	double precision val(nlvdi,nkn)
	logical, parameter :: belem = .false.

	if( bmpi_skip ) return

	call shympi_exchange_internal_d(belem,1,nlvdi,nkn,ilhkv &
     &			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_3d_node_d

!******************************************************************

	subroutine shympi_exchange_3d0_node_r(val)

	use basin
	use levels

	real val(0:nlvdi,nkn)
	logical, parameter :: belem = .false.

	if( bmpi_skip ) return

	call shympi_exchange_internal_r(belem,0,nlvdi,nkn,ilhkv &
     &			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_3d0_node_r

!******************************************************************

	subroutine shympi_exchange_3d_elem_r(val)

	use basin
	use levels

	real val(nlvdi,nel)
	logical, parameter :: belem = .true.

	if( bmpi_skip ) return

	call shympi_exchange_internal_r(belem,1,nlvdi,nel,ilhv &
     &			,ghost_elems_in,ghost_elems_out,val)

	end subroutine shympi_exchange_3d_elem_r

!******************************************************************

	subroutine shympi_exchange_2d_node_i(val)

	use basin
	use levels

	integer val(nkn)
	logical, parameter :: belem = .false.

	if( bmpi_skip ) return

	call shympi_exchange_internal_i(belem,1,1,nkn,ilhkv &
     &			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_2d_node_i

!*******************************

	subroutine shympi_exchange_2d_node_r(val)

	use basin
	use levels

	real val(nkn)
	logical, parameter :: belem = .false.

	if( bmpi_skip ) return

	call shympi_exchange_internal_r(belem,1,1,nkn,ilhkv &
     &			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_2d_node_r

!*******************************

	subroutine shympi_exchange_2d_node_d(val)

	use basin
	use levels

	double precision val(nkn)
	logical, parameter :: belem = .false.

	if( bmpi_skip ) return

	call shympi_exchange_internal_d(belem,1,1,nkn,ilhkv &
     &			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_2d_node_d

!******************************************************************

	subroutine shympi_exchange_2d_elem_i(val)

	use basin
	use levels

	integer val(nel)
	logical, parameter :: belem = .true.

	if( bmpi_skip ) return

	call shympi_exchange_internal_i(belem,1,1,nel,ilhv &
     &			,ghost_elems_in,ghost_elems_out,val)

	end subroutine shympi_exchange_2d_elem_i

!*******************************

	subroutine shympi_exchange_2d_elem_r(val)

	use basin
	use levels

	real val(nel)
	logical, parameter :: belem = .true.

	if( bmpi_skip ) return

	call shympi_exchange_internal_r(belem,1,1,nel,ilhv &
     &			,ghost_elems_in,ghost_elems_out,val)

	end subroutine shympi_exchange_2d_elem_r

!*******************************

	subroutine shympi_exchange_2d_elem_d(val)

	use basin
	use levels

	double precision val(nel)
	logical, parameter :: belem = .true.

	if( bmpi_skip ) return

	call shympi_exchange_internal_d(belem,1,1,nel,ilhv &
     &			,ghost_elems_in,ghost_elems_out,val)

	end subroutine shympi_exchange_2d_elem_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_check_2d_node_i(val,text)

	use basin

	integer val(nkn)
	character*(*) text

	logical, parameter :: belem = .false.
	integer aux(nkn)

	if( bmpi_skip ) return

	aux = val
	call shympi_exchange_2d_node_i(aux)
	call shympi_check_array_i(belem,1,nkn,nkn,val,aux,text)

	end subroutine shympi_check_2d_node_i

!*******************************

	subroutine shympi_check_2d_node_r(val,text)

	use basin

	real val(nkn)
	character*(*) text

	logical, parameter :: belem = .false.
	real aux(nkn)

	if( bmpi_skip ) return

	aux = val
	call shympi_exchange_2d_node_r(aux)
	call shympi_check_array_r(belem,1,nkn,nkn,val,aux,text)

	end subroutine shympi_check_2d_node_r

!*******************************

	subroutine shympi_check_2d_node_d(val,text)

	use basin

	double precision val(nkn)
	character*(*) text

	logical, parameter :: belem = .false.
	double precision aux(nkn)

	if( bmpi_skip ) return

	aux = val
	call shympi_exchange_2d_node_d(aux)
	call shympi_check_array_d(belem,1,nkn,nkn,val,aux,text)

	end subroutine shympi_check_2d_node_d

!******************************************************************

	subroutine shympi_check_3d_node_r(val,text)

	use basin
	use levels

	real val(nlvdi,nkn)
	character*(*) text

	logical, parameter :: belem = .false.
	integer nt
	real aux(nlvdi,nkn)

	if( bmpi_skip ) return

	nt = nlvdi*nkn
	aux = val
	call shympi_exchange_3d_node_r(aux)
	call shympi_check_array_r(belem,nlvdi,nkn,nt,val,aux,text)

	end subroutine shympi_check_3d_node_r

!******************************************************************

	subroutine shympi_check_3d0_node_r(val,text)

	use basin
	use levels

	real val(0:nlvdi,nkn)
	character*(*) text

	logical, parameter :: belem = .false.
	integer nt
	real aux(0:nlvdi,nkn)

	if( bmpi_skip ) return

	nt = (nlvdi+1)*nkn
	aux = val
	call shympi_exchange_3d0_node_r(aux)
	call shympi_check_array_r(belem,nlvdi+1,nkn,nt,val,aux,text)

	end subroutine shympi_check_3d0_node_r

!******************************************************************

	subroutine shympi_check_2d_elem_i(val,text)

	use basin

	integer val(nel)
	character*(*) text

	logical, parameter :: belem = .true.
	integer aux(nel)

	if( bmpi_skip ) return

	aux = val
	call shympi_exchange_2d_elem_i(aux)
	call shympi_check_array_i(belem,1,nel,nel,val,aux,text)

	end subroutine shympi_check_2d_elem_i

!*******************************

	subroutine shympi_check_2d_elem_r(val,text)

	use basin

	real val(nel)
	character*(*) text

	logical, parameter :: belem = .true.
	real aux(nel)

	if( bmpi_skip ) return

	aux = val
	call shympi_exchange_2d_elem_r(aux)
	call shympi_check_array_r(belem,1,nel,nel,val,aux,text)

	end subroutine shympi_check_2d_elem_r

!*******************************

	subroutine shympi_check_2d_elem_d(val,text)

	use basin

	double precision val(nel)
	character*(*) text

	logical, parameter :: belem = .true.
	double precision aux(nel)

	if( bmpi_skip ) return

	aux = val
	call shympi_exchange_2d_elem_d(aux)
	call shympi_check_array_d(belem,1,nel,nel,val,aux,text)

	end subroutine shympi_check_2d_elem_d

!******************************************************************

	subroutine shympi_check_3d_elem_r(val,text)

	use basin
	use levels

	real val(nlvdi,nel)
	character*(*) text

	logical, parameter :: belem = .true.
	integer nt
	real, allocatable :: aux(:,:)

	if( bmpi_skip ) return

	allocate(aux(nlvdi,nel))

	nt = nlvdi*nel
	aux = val
	call shympi_barrier
	call shympi_exchange_3d_elem_r(aux)
	call shympi_check_array_r(belem,nlvdi,nel,nt,val,aux,text)

	end subroutine shympi_check_3d_elem_r

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_check_array_i(belem,nl,nh,n,a1,a2,text)

	logical belem
	integer nl,nh,n
	integer a1(n),a2(n)
	character*(*) text

	integer i,icount,ih,il,ih_ext
	integer, parameter :: imax = 10
	integer maxdif

        if( .not. all( a1 == a2 ) ) then
	  maxdif = maxval( abs(a1-a2) )
          write(6,*) 'arrays are different on ghost items: ' // text
          write(6,*) 'process id: ',my_id
          write(6,*) 'belem: ',belem
          write(6,*) 'total array size: ',n
          write(6,*) 'total differences: ',count(a1/=a2)
          write(6,*) 'max difference: ',maxdif
	  if( count(a1/=a2) > imax ) then
            write(6,*) 'showing only maximum ',imax,' differences'
	  end if
          write(6,*) '      id       i   ihext      ih      l'
	  call shympi_make_debug_text(belem,nh)
          write(6,*) trim(textd)
	  icount = 0
	  do i=1,n
	    if( a1(i) /= a2(i) ) then
	      ih = 1 + (i-1)/nl
	      il = 1 + mod(i-1,nl)
	      ih_ext = make_external_number(belem,ih)
	      icount = icount + 1
	      write(6,1000) my_id,i,ih_ext,ih,il,a1(i),a2(i)
	    end if
	    if( imax > 0 .and. icount >= imax ) exit
	  end do
	  flush(6)
	  call error_stop('error stop shympi_check_array_i')
        end if

	return
 !1000	format(5i8,2f18.6)
 1000	format(7i8)
	end subroutine shympi_check_array_i

!*******************************

	subroutine shympi_check_array_r(belem,nl,nh,n,a1,a2,text)

	logical belem
	integer nl,nh,n
	real a1(n),a2(n)
	character*(*) text

	integer i,icount,ih,il,ih_ext
	integer, parameter :: imax = 10
	real maxdif

	!integer ieext

        if( .not. all( a1 == a2 ) ) then
	  maxdif = maxval( abs(a1-a2) )
          write(6,*) 'arrays are different on ghost items: ' // text
          write(6,*) 'process id: ',my_id
          write(6,*) 'array size: ',n,nh,nl
          write(6,*) 'total differences: ',count(a1/=a2)
          write(6,*) 'max difference: ',maxdif
	  if( count(a1/=a2) > imax ) then
            write(6,*) 'showing only maximum ',imax,' differences'
	  end if
	  call shympi_make_debug_text(belem,nh)
          write(6,*) trim(textd)
	  icount = 0
	  do i=1,n
	    if( a1(i) /= a2(i) ) then
	      ih = 1 + (i-1)/nl
	      il = 1 + mod(i-1,nl)
	      ih_ext = make_external_number(belem,ih)
	      icount = icount + 1
	      write(6,1000) my_id,i,ih_ext,ih,il,a1(i),a2(i)
	    end if
	    if( imax > 0 .and. icount >= imax ) exit
	  end do
	  flush(6)
	  call error_stop('error stop shympi_check_array_r')
        end if

	return
 1000	format(1x,5i8,2f18.6)
	end subroutine shympi_check_array_r

!*******************************

	subroutine shympi_check_array_d(belem,nl,nh,n,a1,a2,text)

	logical belem
	integer nl,nh,n
	double precision a1(n),a2(n)
	character*(*) text

	integer i,icount,ih,il,ih_ext
	integer, parameter :: imax = 10
	double precision maxdif

        if( .not. all( a1 == a2 ) ) then
	  maxdif = maxval( abs(a1-a2) )
          write(6,*) 'arrays are different on ghost items: ' // text
          write(6,*) 'process id: ',my_id
          write(6,*) 'total array size: ',n
          write(6,*) 'total differences: ',count(a1/=a2)
          write(6,*) 'max difference: ',maxdif
	  if( count(a1/=a2) > imax ) then
            write(6,*) 'showing only maximum ',imax,' differences'
	  end if
	  call shympi_make_debug_text(belem,nh)
          write(6,*) trim(textd)
	  icount = 0
	  do i=1,n
	    if( a1(i) /= a2(i) ) then
	      ih = 1 + (i-1)/nl
	      il = 1 + mod(i-1,nl)
	      ih_ext = make_external_number(belem,ih)
	      icount = icount + 1
	      !write(6,*) my_id,i,a1(i),a2(i)
	      write(6,1000) my_id,i,ih_ext,ih,il,a1(i),a2(i)
	    end if
	    if( imax > 0 .and. icount >= imax ) exit
	  end do
	  flush(6)
	  call error_stop('error stop shympi_check_array_d')
        end if

	return
 1000	format(1x,5i8,2f18.6)
	end subroutine shympi_check_array_d

!******************************************************************

	subroutine shympi_make_debug_text(belem,nh)

	logical belem
	integer nh

        textk = '   my_id       i    kext       k       l'
        texte = '   my_id       i   ieext      ie       l'
	text2 = '              val1              val2'
!                1234567890123456789012345678901234567890
	if( belem ) then
	  ip_ext => ip_ext_elem
	  textd = trim(texte) // text2
	else
	  ip_ext => ip_ext_node
	  textd = trim(textk) // text2
	end if

	end

!******************************************************************
	
	function make_external_number(belem,iint)

	use basin

	integer make_external_number
	logical belem
	integer iint

	make_external_number = 0

	if( belem ) then
	  if( iint <= nel ) make_external_number = ipev(iint)
	else
	  if( iint <= nkn ) make_external_number = ipv(iint)
	end if

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_operate_all(balloper)

	logical balloper

	bmpi_alloper = balloper
	call shympi_operate_all_internal(balloper)

	end subroutine shympi_operate_all

!******************************************************************

	subroutine shympi_gather_scalar_i(val,vals)

	integer val
	integer vals(n_threads)

	integer ni,no
	integer valv(1)

	if( bmpi_skip ) then
	  vals(1) = val
	  return
	end if

	valv(1) = val
	ni = 1
	no = 1

	if( bmpi_alloper ) then
	  call shympi_allgather_internal_i(ni,no,valv,vals)
	else
	  call shympi_gather_internal_i(ni,no,valv,vals)
	end if

	end subroutine shympi_gather_scalar_i

!*******************************

	subroutine shympi_gather_array_2d_i(val,vals)

	integer val(:)
	integer vals(:,:)

	integer ni,no
	integer, allocatable :: aux(:)

	if( bmpi_skip ) then
	  vals(:,1) = val(:)
	  return
	end if

	ni = size(val)
	no = size(vals,1)
	if( ni > no ) stop 'error stop shympi_gather_array_2d_i: ni>no'

	allocate(aux(no))
	aux = 0.
	aux(1:ni) = val(1:ni)

	if( bmpi_alloper ) then
	  call shympi_allgather_internal_i(no,no,aux,vals)
	else
	  call shympi_gather_internal_i(no,no,aux,vals)
	end if

	end subroutine shympi_gather_array_2d_i

!*******************************

	subroutine shympi_gather_array_2d_r(val,vals)

	real val(:)
	real vals(:,:)

	integer ni,no
	real, allocatable :: aux(:)

	if( bmpi_skip ) then
	  vals(:,1) = val(:)
	  return
	end if

	ni = size(val)
	no = size(vals,1)
	if( ni > no ) stop 'error stop shympi_gather_array_2d_r: ni>no'

	allocate(aux(no))
	aux(ni+1:) = 0.
	aux(1:ni) = val(1:ni)

	if( bmpi_alloper ) then
	  call shympi_allgather_internal_r(no,no,aux,vals)
	else
	  call shympi_gather_internal_r(no,no,aux,vals)
	end if

	end subroutine shympi_gather_array_2d_r

!*******************************

	subroutine shympi_gather_array_2d_d(val,vals)

	double precision val(:)
	double precision vals(:,:)

	integer ni,no
	double precision, allocatable :: aux(:)

	if( bmpi_skip ) then
	  vals(:,1) = val(:)
	  return
	end if

	ni = size(val)
	no = size(vals,1)
	if( ni > no ) stop 'error stop shympi_gather_array_2d_d: ni>no'

	allocate(aux(no))
	aux = 0.
	aux(1:ni) = val(1:ni)

	if( bmpi_alloper ) then
	  call shympi_allgather_internal_d(no,no,aux,vals)
	else
	  call shympi_gather_internal_d(no,no,aux,vals)
	end if

	end subroutine shympi_gather_array_2d_d

!*******************************

	subroutine shympi_gather_array_3d_i(val,vals)

	integer val(:,:)
	integer vals(:,:,:)

	integer ni1,ni2,no1,no2
	integer ni,no
	integer, allocatable :: aux(:,:)

	if( bmpi_skip ) then
	  vals(:,:,1) = val(:,:)
	  return
	end if

	ni1 = size(val,1)
	ni2 = size(val,2)
	no1 = size(vals,1)
	no2 = size(vals,2)
	if( ni1 > no1 ) stop 'error stop shympi_gather_array_3d_i: ni1>no1'
	if( ni2 > no2 ) stop 'error stop shympi_gather_array_3d_i: ni2>no2'

	allocate(aux(no1,no2))
	aux = 0.
	aux(1:ni1,1:ni2) = val(1:ni1,1:ni2)

	ni = no1 * no2
	no = no1 * no2

	if( bmpi_alloper ) then
	  call shympi_allgather_internal_i(ni,no,aux,vals)
	else
	  call shympi_gather_internal_i(ni,no,aux,vals)
	end if

	end subroutine shympi_gather_array_3d_i

!*******************************

	subroutine shympi_gather_array_3d_r(val,vals)

	real val(:,:)
	real vals(:,:,:)

	integer ni1,ni2,no1,no2
	integer ni,no
	real, allocatable :: aux(:,:)

	if( bmpi_skip ) then
	  vals(:,:,1) = val(:,:)
	  return
	end if

	ni1 = size(val,1)
	ni2 = size(val,2)
	no1 = size(vals,1)
	no2 = size(vals,2)	!this is nn_max
	if( ni1 > no1 ) stop 'error stop shympi_gather_array_3d_r: ni1>no1'
	if( ni2 > no2 ) stop 'error stop shympi_gather_array_3d_r: ni2>no2'

	allocate(aux(no1,no2))
	aux = 0.
	aux(1:ni1,1:ni2) = val(1:ni1,1:ni2)

	ni = no1 * no2
	no = no1 * no2

	if( bmpi_alloper ) then
	  call shympi_allgather_internal_r(ni,no,aux,vals)
	else
	  call shympi_gather_internal_r(ni,no,aux,vals)
	end if

	end subroutine shympi_gather_array_3d_r

!*******************************

	subroutine shympi_gather_array_3d_d(val,vals)

	double precision val(:,:)
	double precision vals(:,:,:)

	integer ni1,ni2,no1,no2
	integer ni,no
	double precision, allocatable :: aux(:,:)

	if( bmpi_skip ) then
	  vals(:,:,1) = val(:,:)
	  return
	end if

	ni1 = size(val,1)
	ni2 = size(val,2)
	no1 = size(vals,1)
	no2 = size(vals,2)
	if( ni1 > no1 ) stop 'error stop shympi_gather_array_3d_d: ni1>no1'
	if( ni2 > no2 ) stop 'error stop shympi_gather_array_3d_d: ni2>no2'

	allocate(aux(no1,no2))
	aux = 0.
	aux(1:ni1,1:ni2) = val(1:ni1,1:ni2)

	ni = no1 * no2
	no = no1 * no2

	if( bmpi_alloper ) then
	  call shympi_allgather_internal_d(ni,no,aux,vals)
	else
	  call shympi_gather_internal_d(ni,no,aux,vals)
	end if

	end subroutine shympi_gather_array_3d_d

!*******************************

	subroutine shympi_gather_array_fix_i(nfix,val,vals)

	integer nfix
	integer val(:,:)
	integer vals(:,:,:)

	integer ni1,ni2,no1,no2
	integer ni,no
	integer, allocatable :: aux(:,:)

	if( bmpi_skip ) then
	  vals(:,:,1) = val(:,:)
	  return
	end if

	ni1 = size(val,1)
	ni2 = size(val,2)
	no1 = size(vals,1)
	no2 = size(vals,2)

	if( ni1 /= nfix .or. no1 /= nfix ) then
	  write(6,*) nfix,ni1,no1
	  stop 'error stop shympi_gather_array_fix_i: incomp first dim'
	end if
	if( ni2 > no2 ) stop 'error stop shympi_gather_array_fix_i: ni2>no2'

	allocate(aux(no1,no2))
	aux = 0.
	aux(1:ni1,1:ni2) = val(1:ni1,1:ni2)

	ni = no1 * no2
	no = no1 * no2

	if( bmpi_alloper ) then
	  call shympi_allgather_internal_i(ni,no,aux,vals)
	else
	  call shympi_gather_internal_i(ni,no,aux,vals)
	end if

	end subroutine shympi_gather_array_fix_i

!*******************************

	subroutine shympi_gather_array_fix_r(nfix,val,vals)

	integer nfix
	real val(:,:)
	real vals(:,:,:)

	integer ni1,ni2,no1,no2
	integer ni,no
	real, allocatable :: aux(:,:)

	if( bmpi_skip ) then
	  vals(:,:,1) = val(:,:)
	  return
	end if

	ni1 = size(val,1)
	ni2 = size(val,2)
	no1 = size(vals,1)
	no2 = size(vals,2)

	if( ni1 /= nfix .or. no1 /= nfix ) then
	  write(6,*) nfix,ni1,no1
	  stop 'error stop shympi_gather_array_fix_r: incomp first dim'
	end if
	if( ni2 > no2 ) stop 'error stop shympi_gather_array_fix_r: ni2>no2'

	allocate(aux(no1,no2))
	aux = 0.
	aux(1:ni1,1:ni2) = val(1:ni1,1:ni2)

	ni = no1 * no2
	no = no1 * no2

	if( bmpi_alloper ) then
	  call shympi_allgather_internal_r(ni,no,aux,vals)
	else
	  call shympi_gather_internal_r(ni,no,aux,vals)
	end if

	end subroutine shympi_gather_array_fix_r

!*******************************

	subroutine shympi_gather_array_fix_d(nfix,val,vals)

	integer nfix
	double precision val(:,:)
	double precision vals(:,:,:)

	integer ni1,ni2,no1,no2
	integer ni,no
	double precision, allocatable :: aux(:,:)

	if( bmpi_skip ) then
	  vals(:,:,1) = val(:,:)
	  return
	end if

	ni1 = size(val,1)
	ni2 = size(val,2)
	no1 = size(vals,1)
	no2 = size(vals,2)

	if( ni1 /= nfix .or. no1 /= nfix ) then
	  write(6,*) nfix,ni1,no1
	  stop 'error stop shympi_gather_array_fix_d: incomp first dim'
	end if
	if( ni2 > no2 ) stop 'error stop shympi_gather_array_fix_d: ni2>no2'

	allocate(aux(no1,no2))
	aux = 0.
	aux(1:ni1,1:ni2) = val(1:ni1,1:ni2)

	ni = no1 * no2
	no = no1 * no2

	if( bmpi_alloper ) then
	  call shympi_allgather_internal_d(ni,no,aux,vals)
	else
	  call shympi_gather_internal_d(ni,no,aux,vals)
	end if

	end subroutine shympi_gather_array_fix_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_gather_root_array_2d_d(val,vals)

	double precision val(:)
	double precision vals(:,:)

	integer ni,no
	double precision, allocatable :: aux(:)

	if( bmpi_skip ) then
	  vals(:,1) = val(:)
	  return
	end if

	ni = size(val)
	no = size(vals,1)
	if( ni > no ) stop 'error stop shympi_gather_array_2d_r: ni>no'

	allocate(aux(no))
	aux = 0.
	aux(1:ni) = val(1:ni)

	call shympi_gather_internal_d(no,no,aux,vals)

	end subroutine shympi_gather_root_array_2d_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_gather_and_sum_i(val)

	integer val(:)

	integer n,no
	integer vals(size(val),n_threads)

	if( bmpi_skip ) return

	n = size(val)
	no = n
	call shympi_allgather_internal_i(n,no,val,vals)
	val(:) = SUM(vals,dim=2)

	end subroutine shympi_gather_and_sum_i

!*******************************

	subroutine shympi_gather_and_sum_r(val)

	real val(:)

	integer n,no
	real vals(size(val),n_threads)

	if( bmpi_skip ) return

	n = size(val)
	no = n
	call shympi_allgather_internal_r(n,no,val,vals)
	val(:) = SUM(vals,dim=2)

	end subroutine shympi_gather_and_sum_r

!*******************************

	subroutine shympi_gather_and_sum_d(val)

	double precision val(:)

	integer n,no
	double precision vals(size(val),n_threads)

	if( bmpi_skip ) return

	n = size(val)
	no = n
	call shympi_allgather_internal_d(n,no,val,vals)
	val(:) = SUM(vals,dim=2)

	end subroutine shympi_gather_and_sum_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_bcast_scalar_i(val)

	integer val

	integer n

	if( bmpi_skip ) return

	n = 1
	call shympi_bcast_internal_i(n,val)

	end subroutine shympi_bcast_scalar_i

!*******************************

	subroutine shympi_bcast_array_r(val)

	real val(:)

	integer n

	if( bmpi_skip ) return

	n = size(val)
	call shympi_bcast_internal_r(n,val)

	end subroutine shympi_bcast_array_r

!*******************************

	subroutine shympi_bcast_array_d(val)

	double precision val(:)

	integer n

	if( bmpi_skip ) return

	n = size(val)
	call shympi_bcast_internal_d(n,val)

	end subroutine shympi_bcast_array_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_collect_node_value_2d_r(k,vals,val)

	integer k
	real vals(:)
	real val

	val = 0
	if( k > 0 .and. k <= nkn_unique ) val = vals(k)
	if( bmpi ) val = shympi_sum(val)

	end subroutine shympi_collect_node_value_2d_r

!*******************************

	subroutine shympi_collect_node_value_2d_i(k,vals,val)

	integer k
	integer vals(:)
	integer val

	val = 0
	if( k > 0 .and. k <= nkn_unique ) val = vals(k)
	if( bmpi ) val = shympi_sum(val)

	end subroutine shympi_collect_node_value_2d_i

!*******************************

	subroutine shympi_collect_node_value_3d_r(k,vals,val)

	integer k
	real vals(:,:)
	real val(:)

	integer ni
	real vaux(size(val),n_threads)

	ni = size(vals,1)

	val = 0
	if( k > 0 .and. k <= nkn_unique ) val(1:ni) = vals(1:ni,k)

	if( bmpi ) then
	  call shympi_gather(val,vaux)
	  val = SUM(vaux,dim=2)
	end if

	end subroutine shympi_collect_node_value_3d_r

!*******************************

!	subroutine shympi_collect_node_value_3d_i(k,vals,val)
!
!	integer k
!	integer vals(:)
!	integer val
!
!	val = 0
!	if( k > 0 .and. k <= nkn_unique ) val = vals(k)
!	if( bmpi ) val = shympi_sum(val)
!
!	end subroutine shympi_collect_node_value_3d_i

!*******************************

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_find_node(ke,ki,id)

! returns internal node number and id of domain given external node number

	use basin
	use mod_compatibility

	integer, intent(in)  ::  ke	! external node number
	integer, intent(out) ::  ki	! internal node number
	integer, intent(out) ::  id	! domain

	integer k,ka,i,j
	integer vals(1,n_threads)
	integer kk(1),kkk(n_threads)

	do k=1,nkn_unique
	  if( ipv(k) == ke ) exit
	end do
	if( k > nkn_unique ) k = 0
	!ka = findloc(ipv(1:nkn_unique),ke,1)
	ka = gfindloc(ipv(1:nkn_unique),ke)
	if( k /= ka ) then
	  write(6,*) k,ka,nkn_unique,my_id
	  stop 'error stop shympi_find_node: internal (1)'
	end if

	if( bmpi_skip ) then
	  ki = k
	  id = 0
	  return
	end if

	kk(1) = k
	call shympi_allgather_internal_i(1,1,kk,vals)
	kkk(:) = vals(1,:)

	call shympi_check_find_item(kkk,ke,k)

	do i=1,n_threads
	  if( kkk(i) /= 0 ) exit
	end do

	j = maxloc(kkk,1)
	if( i /= j ) stop 'error stop shympi_find_node: internal (2)'

	ki = kkk(j)
	id = j-1

	end subroutine shympi_find_node

!*******************************

	subroutine shympi_find_elem(ieext,ieint,id)

! returns internal elem number and id of domain given external elem number

	use basin
	use mod_compatibility

	integer, intent(in)  ::  ieext	! external elem number
	integer, intent(out) ::  ieint	! internal elem number
	integer, intent(out) ::  id	! domain

	integer ie,iee,i,j
	integer vals(1,n_threads)
	integer kk(1),kkk(n_threads)

	do ie=1,nel_unique
	  if( ipev(ie) == ieext ) exit
	end do
	if( ie > nel_unique ) ie = 0
	!iee = findloc(ipev(1:nel_unique),ieext,1)
	iee = gfindloc(ipev(1:nel_unique),ieext)
	if( ie /= iee ) stop 'error stop shympi_find_elem: internal (1)'

	if( bmpi_skip ) then
	  ieint = ie
	  id = 0
	  return
	end if

	kk(1) = ie
	call shympi_allgather_internal_i(1,1,kk,vals)
	kkk(:) = vals(1,:)

	call shympi_check_find_item(kkk,ieext,ie)

	do i=1,n_threads
	  if( kkk(i) /= 0 ) exit
	end do

	j = maxloc(kkk,1)
	if( i /= j ) stop 'error stop shympi_find_elem: internal (2)'

	ieint = kkk(j)
	id = j-1

	end subroutine shympi_find_elem

!*******************************

	subroutine shympi_check_find_item(items,ext,int)

	integer items(n_threads)
	integer ext,int

	integer ic

        ic = count( items /= 0 )

        if( ic /= 1 ) then              ! error - abort
          if( ic > 1 ) then
            write(6,*) 'item found in more than one domain: '
          else
            write(6,*) 'item not found in any domain:'
          end if
          write(6,*) '==========================='
          write(6,*) n_threads,my_id
          write(6,*) 'nkn_unique = ',nkn_unique
          write(6,*) 'nel_unique = ',nel_unique
          write(6,*) 'external number = ',ext
          write(6,*) 'internal number = ',int
          write(6,*) items
          write(6,*) '==========================='
          call shympi_finalize
          stop 'error stop shympi_check_find_item: more than one or no domain'
        end if

	end subroutine shympi_check_find_item

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_get_array_2d_r(n,vals,val_out)

	use basin
	use levels

	integer n
	real vals(:)
	real val_out(n)

	if( bmpi_skip ) then
	  val_out = vals
	  return
	end if

	call shympi_get_array_internal_r(1,n &
     &                                    ,vals,val_out)

	end subroutine shympi_get_array_2d_r

!*******************************

	subroutine shympi_get_array_2d_i(n,vals,val_out)

	use basin
	use levels

	integer n
	integer vals(:)
	integer val_out(n)

	if( bmpi_skip ) then
	  val_out = vals
	  return
	end if

	call shympi_get_array_internal_i(1,n &
     &                                    ,vals,val_out)

	end subroutine shympi_get_array_2d_i

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_l2g_array_2d_r(vals,val_out)

	real vals(:)
	real val_out(:)

	logical bnode,belem
	integer nous
	real, allocatable :: val_domain(:,:)

	if( bmpi_skip ) then
	  val_out = vals
	  return
	end if

	nous = size(val_out,1)

        bnode = ( nous == nkn_global )
        belem = ( nous == nel_global )

	allocate( val_domain(nn_max,n_threads) )
	call shympi_gather_array_2d_r(vals,val_domain)

	if( bnode ) then
	  call shympi_copy_2d_r(val_domain,nous,val_out &
     &				,nkn_domains,nk_max,ip_int_nodes)
	else if( belem ) then
	  call shympi_copy_2d_r(val_domain,nous,val_out &
     &				,nel_domains,ne_max,ip_int_elems)
	else
	  write(6,*) nous,nkn_global,nel_global
	  stop 'error stop shympi_l2g_array_2d_r: (1)'
	end if

	end subroutine shympi_l2g_array_2d_r

!*******************************

	subroutine shympi_l2g_array_2d_i(vals,val_out)

	integer vals(:)
	integer val_out(:)

	logical bnode,belem
	integer nous
	integer, allocatable :: val_domain(:,:)

	if( bmpi_skip ) then
	  val_out = vals
	  return
	end if

	nous = size(val_out,1)

        bnode = ( nous == nkn_global )
        belem = ( nous == nel_global )

	allocate( val_domain(nn_max,n_threads) )
	call shympi_gather_array_2d_i(vals,val_domain)

	if( bnode ) then
	  call shympi_copy_2d_i(val_domain,nous,val_out &
     &				,nkn_domains,nk_max,ip_int_nodes)
	else if( belem ) then
	  call shympi_copy_2d_i(val_domain,nous,val_out &
     &				,nel_domains,ne_max,ip_int_elems)
	else
	  stop 'error stop shympi_l2g_array_2d_i: (1)'
	end if

	end subroutine shympi_l2g_array_2d_i

!*******************************

	subroutine shympi_l2g_array_2d_d(vals,val_out)

	double precision vals(:)
	double precision val_out(:)

	logical bnode,belem
	integer nous
	double precision, allocatable :: val_domain(:,:)

	if( bmpi_skip ) then
	  val_out = vals
	  return
	end if

	nous = size(val_out,1)

        bnode = ( nous == nkn_global )
        belem = ( nous == nel_global )

	allocate( val_domain(nn_max,n_threads) )
	call shympi_gather_array_2d_d(vals,val_domain)

	if( bnode ) then
	  call shympi_copy_2d_d(val_domain,nous,val_out &
     &				,nkn_domains,nk_max,ip_int_nodes)
	else if( belem ) then
	  call shympi_copy_2d_d(val_domain,nous,val_out &
     &				,nel_domains,ne_max,ip_int_elems)
	else
	  stop 'error stop shympi_l2g_array_2d_d: (1)'
	end if

	end subroutine shympi_l2g_array_2d_d

!*******************************

	subroutine shympi_l2g_array_3d_r(vals,val_out)

	real vals(:,:)
	real val_out(:,:)

	logical bnode,belem
	integer noh,nov
	integer nih,niv
	real, allocatable :: val_domain(:,:,:)

	if( bmpi_skip ) then
	  val_out = vals
	  return
	end if

	nih = size(vals,2)
	niv = size(vals,1)
	noh = size(val_out,2)
	nov = size(val_out,1)

        bnode = ( noh == nkn_global )
        belem = ( noh == nel_global )

	if( nih > nn_max ) then
	  stop 'error stop shympi_l2g_array_3d_r: nih>nn_max'
	end if

	allocate(val_domain(nov,nn_max,n_threads))
	val_domain = 0.

	call shympi_gather_array_3d_r(vals,val_domain)

	if( bnode ) then
	  call shympi_copy_3d_r(val_domain,nov,noh,val_out &
     &				,nkn_domains,nk_max,ip_int_nodes)
	else if( belem ) then
	  call shympi_copy_3d_r(val_domain,nov,noh,val_out &
     &				,nel_domains,ne_max,ip_int_elems)
	else
	  write(6,*) noh,nov,nkn_global,nel_global
	  stop 'error stop shympi_l2g_array_3d_r: (1)'
	end if

	end subroutine shympi_l2g_array_3d_r

!*******************************

	subroutine shympi_l2g_array_3d_i(vals,val_out)

	integer vals(:,:)
	integer val_out(:,:)

	logical bnode,belem
	integer noh,nov
	integer nih,niv
	integer, allocatable :: val_domain(:,:,:)

	if( bmpi_skip ) then
	  val_out = vals
	  return
	end if

	nih = size(vals,2)
	niv = size(vals,1)
	noh = size(val_out,2)
	nov = size(val_out,1)

        bnode = ( noh == nkn_global )
        belem = ( noh == nel_global )

	allocate(val_domain(nov,nn_max,n_threads))

	call shympi_gather_array_3d_i(vals,val_domain)

	if( bnode ) then
	  call shympi_copy_3d_i(val_domain,nov,noh,val_out &
     &				,nkn_domains,nk_max,ip_int_nodes)
	else if( belem ) then
	  call shympi_copy_3d_i(val_domain,nov,noh,val_out &
     &				,nel_domains,ne_max,ip_int_elems)
	else
	  write(6,*) noh,nov,nkn_global,nel_global
	  stop 'error stop shympi_l2g_array_3d_i: (1)'
	end if

	end subroutine shympi_l2g_array_3d_i

!*******************************

	subroutine shympi_l2g_array_3d_d(vals,val_out)

	double precision vals(:,:)
	double precision val_out(:,:)

	logical bnode,belem
	integer noh,nov
	integer nih,niv
	double precision, allocatable :: val_domain(:,:,:)

	if( bmpi_skip ) then
	  val_out = vals
	  return
	end if

	nih = size(vals,2)
	niv = size(vals,1)
	noh = size(val_out,2)
	nov = size(val_out,1)

        bnode = ( noh == nkn_global )
        belem = ( noh == nel_global )

	allocate(val_domain(nov,nn_max,n_threads))

	call shympi_gather_array_3d_d(vals,val_domain)

	if( bnode ) then
	  call shympi_copy_3d_d(val_domain,nov,noh,val_out &
     &				,nkn_domains,nk_max,ip_int_nodes)
	else if( belem ) then
	  call shympi_copy_3d_d(val_domain,nov,noh,val_out &
     &				,nel_domains,ne_max,ip_int_elems)
	else
	  write(6,*) noh,nov,nkn_global,nel_global
	  stop 'error stop shympi_l2g_array_3d_d: (1)'
	end if

	end subroutine shympi_l2g_array_3d_d

!*******************************

	subroutine shympi_l2g_array_fix_i(nfix,vals,val_out)

	integer nfix
	integer vals(:,:)
	integer val_out(:,:)

	logical bnode,belem
	integer noh,nov
	integer nih,niv
	integer, allocatable :: val_domain(:,:,:)

	if( bmpi_skip ) then
	  val_out = vals
	  return
	end if

	nih = size(vals,2)
	niv = size(vals,1)
	noh = size(val_out,2)
	nov = size(val_out,1)

        bnode = ( noh == nkn_global )
        belem = ( noh == nel_global )

	if( niv /= nfix .or. nov /= nfix ) then
	  write(6,*) nfix,niv,nov
	  stop 'error stop shympi_l2g_array_fix: incomp first dim'
	end if

	allocate(val_domain(nfix,nn_max,n_threads))
	val_domain = 0.

	call shympi_gather_array_fix_i(nfix,vals,val_domain)

	if( bnode ) then
	  call shympi_copy_3d_i(val_domain,nov,noh,val_out &
     &				,nkn_domains,nk_max,ip_int_nodes)
	else if( belem ) then
	  call shympi_copy_3d_i(val_domain,nov,noh,val_out &
     &				,nel_domains,ne_max,ip_int_elems)
	else
	  write(6,*) noh,nov,nkn_global,nel_global
	  stop 'error stop shympi_l2g_array_fix: (1)'
	end if

	end subroutine shympi_l2g_array_fix_i

!*******************************

	subroutine shympi_l2g_array_fix_r(nfix,vals,val_out)

	integer nfix
	real vals(:,:)
	real val_out(:,:)

	logical bnode,belem
	integer noh,nov
	integer nih,niv
	real, allocatable :: val_domain(:,:,:)

	if( bmpi_skip ) then
	  val_out = vals
	  return
	end if

	nih = size(vals,2)
	niv = size(vals,1)
	noh = size(val_out,2)
	nov = size(val_out,1)

        bnode = ( noh == nkn_global )
        belem = ( noh == nel_global )

	if( niv /= nfix .or. nov /= nfix ) then
	  write(6,*) nfix,niv,nov
	  stop 'error stop shympi_l2g_array_fix: incomp first dim'
	end if

	allocate(val_domain(nfix,nn_max,n_threads))
	val_domain = 0.

	call shympi_gather_array_fix_r(nfix,vals,val_domain)

	if( bnode ) then
	  call shympi_copy_3d_r(val_domain,nov,noh,val_out &
     &				,nkn_domains,nk_max,ip_int_nodes)
	else if( belem ) then
	  call shympi_copy_3d_r(val_domain,nov,noh,val_out &
     &				,nel_domains,ne_max,ip_int_elems)
	else
	  write(6,*) noh,nov,nkn_global,nel_global
	  stop 'error stop shympi_l2g_array_fix: (1)'
	end if

	end subroutine shympi_l2g_array_fix_r

!*******************************

	subroutine shympi_l2g_array_fix_d(nfix,vals,val_out)

	integer nfix
	double precision vals(:,:)
	double precision val_out(:,:)

	logical bnode,belem
	integer noh,nov
	integer nih,niv
	double precision, allocatable :: val_domain(:,:,:)

	if( bmpi_skip ) then
	  val_out = vals
	  return
	end if

	nih = size(vals,2)
	niv = size(vals,1)
	noh = size(val_out,2)
	nov = size(val_out,1)

        bnode = ( noh == nkn_global )
        belem = ( noh == nel_global )

	if( niv /= nfix .or. nov /= nfix ) then
	  write(6,*) nfix,niv,nov
	  stop 'error stop shympi_l2g_array_fix: incomp first dim'
	end if

	allocate(val_domain(nfix,nn_max,n_threads))
	val_domain = 0.

	call shympi_gather_array_fix_d(nfix,vals,val_domain)

	if( bnode ) then
	  call shympi_copy_3d_d(val_domain,nov,noh,val_out &
     &				,nkn_domains,nk_max,ip_int_nodes)
	else if( belem ) then
	  call shympi_copy_3d_d(val_domain,nov,noh,val_out &
     &				,nel_domains,ne_max,ip_int_elems)
	else
	  write(6,*) noh,nov,nkn_global,nel_global
	  stop 'error stop shympi_l2g_array_fix: (1)'
	end if

	call shympi_bdebug('copy_finished')

	end subroutine shympi_l2g_array_fix_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_g2l_array_2d_r(val_g,val_l)

	real val_g(:)
	real val_l(:)

	logical bnode,belem
	integer ng,nl,i,ip

	if( bmpi_skip ) then
	  val_l = val_g
	  return
	end if

	ng = size(val_g,1)
	nl = size(val_l,1)

        bnode = ( ng == nkn_global )
        belem = ( ng == nel_global )

	if( bnode ) then
	  do i=1,nl
	    ip = ip_int_node(i)
	    val_l(i) = val_g(ip)
	  end do
	else if( belem ) then
	  do i=1,nl
	    ip = ip_int_elem(i)
	    val_l(i) = val_g(ip)
	  end do
	else
	  stop 'error stop shympi_g2l_array_2d_r: (1)'
	end if

	end subroutine shympi_g2l_array_2d_r

!*******************************

	subroutine shympi_g2l_array_2d_i(val_g,val_l)

	integer val_g(:)
	integer val_l(:)

	logical bnode,belem
	integer ng,nl,i,ip

	if( bmpi_skip ) then
	  val_l = val_g
	  return
	end if

	ng = size(val_g,1)
	nl = size(val_l,1)

        bnode = ( ng == nkn_global )
        belem = ( ng == nel_global )

	if( bnode ) then
	  do i=1,nl
	    ip = ip_int_node(i)
	    val_l(i) = val_g(ip)
	  end do
	else if( belem ) then
	  do i=1,nl
	    ip = ip_int_elem(i)
	    val_l(i) = val_g(ip)
	  end do
	else
	  stop 'error stop shympi_g2l_array_2d_i: (1)'
	end if

	end subroutine shympi_g2l_array_2d_i

!*******************************

	subroutine shympi_g2l_array_2d_d(val_g,val_l)

	double precision val_g(:)
	double precision val_l(:)

	logical bnode,belem
	integer ng,nl,i,ip

	if( bmpi_skip ) then
	  val_l = val_g
	  return
	end if

	ng = size(val_g,1)
	nl = size(val_l,1)

        bnode = ( ng == nkn_global )
        belem = ( ng == nel_global )

	if( bnode ) then
	  do i=1,nl
	    ip = ip_int_node(i)
	    val_l(i) = val_g(ip)
	  end do
	else if( belem ) then
	  do i=1,nl
	    ip = ip_int_elem(i)
	    val_l(i) = val_g(ip)
	  end do
	else
	  stop 'error stop shympi_g2l_array_2d_d: (1)'
	end if

	end subroutine shympi_g2l_array_2d_d

!*******************************

	subroutine shympi_g2l_array_3d_r(val_g,val_l)

	real val_g(:,:)
	real val_l(:,:)

	logical bnode,belem
	integer ng,nl,nz,i,ip

	if( bmpi_skip ) then
	  val_l = val_g
	  return
	end if

	ng = size(val_g,2)
	nl = size(val_l,2)
	nz = size(val_l,1)

        bnode = ( ng == nkn_global )
        belem = ( ng == nel_global )

	if( bnode ) then
	  do i=1,nl
	    ip = ip_int_node(i)
	    val_l(1:nz,i) = val_g(1:nz,ip)
	  end do
	else if( belem ) then
	  do i=1,nl
	    ip = ip_int_elem(i)
	    val_l(1:nz,i) = val_g(1:nz,ip)
	  end do
	else
	  stop 'error stop shympi_g2l_array_3d_r: (1)'
	end if

	end subroutine shympi_g2l_array_3d_r

!*******************************

	subroutine shympi_g2l_array_3d_i(val_g,val_l)

	integer val_g(:,:)
	integer val_l(:,:)

	logical bnode,belem
	integer ng,nl,nz,i,ip

	if( bmpi_skip ) then
	  val_l = val_g
	  return
	end if

	ng = size(val_g,2)
	nl = size(val_l,2)
	nz = size(val_l,1)

        bnode = ( ng == nkn_global )
        belem = ( ng == nel_global )

	if( bnode ) then
	  do i=1,nl
	    ip = ip_int_node(i)
	    val_l(1:nz,i) = val_g(1:nz,ip)
	  end do
	else if( belem ) then
	  do i=1,nl
	    ip = ip_int_elem(i)
	    val_l(1:nz,i) = val_g(1:nz,ip)
	  end do
	else
	  stop 'error stop shympi_g2l_array_3d_i: (1)'
	end if

	end subroutine shympi_g2l_array_3d_i

!*******************************

	subroutine shympi_g2l_array_3d_d(val_g,val_l)

	double precision val_g(:,:)
	double precision val_l(:,:)

	logical bnode,belem
	integer ng,nl,nz,i,ip

	if( bmpi_skip ) then
	  val_l = val_g
	  return
	end if

	ng = size(val_g,2)
	nl = size(val_l,2)
	nz = size(val_l,1)

        bnode = ( ng == nkn_global )
        belem = ( ng == nel_global )

	if( bnode ) then
	  do i=1,nl
	    ip = ip_int_node(i)
	    val_l(1:nz,i) = val_g(1:nz,ip)
	  end do
	else if( belem ) then
	  do i=1,nl
	    ip = ip_int_elem(i)
	    val_l(1:nz,i) = val_g(1:nz,ip)
	  end do
	else
	  stop 'error stop shympi_g2l_array_3d_d: (1)'
	end if

	end subroutine shympi_g2l_array_3d_d

!*******************************

	subroutine shympi_g2l_array_fix_i(nfix,val_g,val_l)

	integer nfix
	integer val_g(:,:)
	integer val_l(:,:)

	logical bnode,belem
	integer ngh,ngv
	integer nlh,nlv
	integer i,ip,nz

	if( bmpi_skip ) then
	  val_l = val_g
	  return
	end if

	ngh = size(val_g,2)
	ngv = size(val_g,1)
	nlh = size(val_l,2)
	nlv = size(val_l,1)

        bnode = ( ngh == nkn_global )
        belem = ( ngh == nel_global )

	if( ngv /= nfix .or. nlv /= nfix ) then
	  write(6,*) nfix,ngv,nlv
	  stop 'error stop shympi_g2l_array_fix: incomp first dim'
	end if

	nz = nfix

	if( bnode ) then
	  do i=1,nlh
	    ip = ip_int_node(i)
	    val_l(1:nz,i) = val_g(1:nz,ip)
	  end do
	else if( belem ) then
	  do i=1,nlh
	    ip = ip_int_elem(i)
	    val_l(1:nz,i) = val_g(1:nz,ip)
	  end do
	else
	  stop 'error stop shympi_g2l_array_fix: (1)'
	end if

	end subroutine shympi_g2l_array_fix_i

!*******************************

	subroutine shympi_g2l_array_fix_r(nfix,val_g,val_l)

	integer nfix
	real val_g(:,:)
	real val_l(:,:)

	logical bnode,belem
	integer ngh,ngv
	integer nlh,nlv
	integer i,ip,nz

	if( bmpi_skip ) then
	  val_l = val_g
	  return
	end if

	ngh = size(val_g,2)
	ngv = size(val_g,1)
	nlh = size(val_l,2)
	nlv = size(val_l,1)

        bnode = ( ngh == nkn_global )
        belem = ( ngh == nel_global )

	if( ngv /= nfix .or. nlv /= nfix ) then
	  write(6,*) nfix,ngv,nlv
	  stop 'error stop shympi_g2l_array_fix: incomp first dim'
	end if

	nz = nfix

	if( bnode ) then
	  do i=1,nlh
	    ip = ip_int_node(i)
	    val_l(1:nz,i) = val_g(1:nz,ip)
	  end do
	else if( belem ) then
	  do i=1,nlh
	    ip = ip_int_elem(i)
	    val_l(1:nz,i) = val_g(1:nz,ip)
	  end do
	else
	  stop 'error stop shympi_g2l_array_fix: (1)'
	end if

	end subroutine shympi_g2l_array_fix_r

!*******************************

	subroutine shympi_g2l_array_fix_d(nfix,val_g,val_l)

	integer nfix
	double precision val_g(:,:)
	double precision val_l(:,:)

	logical bnode,belem
	integer ngh,ngv
	integer nlh,nlv
	integer i,ip,nz

	if( bmpi_skip ) then
	  val_l = val_g
	  return
	end if

	ngh = size(val_g,2)
	ngv = size(val_g,1)
	nlh = size(val_l,2)
	nlv = size(val_l,1)

        bnode = ( ngh == nkn_global )
        belem = ( ngh == nel_global )

	if( ngv /= nfix .or. nlv /= nfix ) then
	  write(6,*) nfix,ngv,nlv
	  stop 'error stop shympi_g2l_array_fix: incomp first dim'
	end if

	nz = nfix

	if( bnode ) then
	  do i=1,nlh
	    ip = ip_int_node(i)
	    val_l(1:nz,i) = val_g(1:nz,ip)
	  end do
	else if( belem ) then
	  do i=1,nlh
	    ip = ip_int_elem(i)
	    val_l(1:nz,i) = val_g(1:nz,ip)
	  end do
	else
	  stop 'error stop shympi_g2l_array_fix: (1)'
	end if

	end subroutine shympi_g2l_array_fix_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_copy_2d_i(val_domain,nous,val_out &
     &				,ndomains,nmax,ip_ints)

	integer val_domain(nn_max,n_threads)
	integer nous
	integer val_out(nous)
	integer ndomains(n_threads)
	integer nmax
	integer ip_ints(nmax,n_threads)

	integer ia,i,n,ip

        do ia=1,n_threads
          n=ndomains(ia)
          do i=1,n
            ip = ip_ints(i,ia)
	    val_out(ip) = val_domain(i,ia)
          end do
        end do

	end subroutine shympi_copy_2d_i

!*******************************

	subroutine shympi_copy_2d_r(val_domain,nous,val_out &
     &				,ndomains,nmax,ip_ints)

	real val_domain(nn_max,n_threads)
	integer nous
	real val_out(nous)
	integer ndomains(n_threads)
	integer nmax
	integer ip_ints(nmax,n_threads)

	integer ia,i,n,ip

        do ia=1,n_threads
          n=ndomains(ia)
          do i=1,n
            ip = ip_ints(i,ia)
	    val_out(ip) = val_domain(i,ia)
          end do
        end do

	end subroutine shympi_copy_2d_r

!*******************************

	subroutine shympi_copy_2d_d(val_domain,nous,val_out &
     &				,ndomains,nmax,ip_ints)

	double precision val_domain(nn_max,n_threads)
	integer nous
	double precision val_out(nous)
	integer ndomains(n_threads)
	integer nmax
	integer ip_ints(nmax,n_threads)

	integer ia,i,n,ip

        do ia=1,n_threads
          n=ndomains(ia)
          do i=1,n
            ip = ip_ints(i,ia)
	    val_out(ip) = val_domain(i,ia)
          end do
        end do

	end subroutine shympi_copy_2d_d

!*******************************

	subroutine shympi_copy_3d_i(val_domain,nov,nous,val_out &
     &				,ndomains,nmax,ip_ints)

	integer nov,nous
	integer val_domain(nov,nn_max,n_threads)
	integer val_out(nov,nous)
	integer ndomains(n_threads)
	integer nmax
	integer ip_ints(nmax,n_threads)

	integer ia,i,n,ip

        do ia=1,n_threads
          n=ndomains(ia)
          do i=1,n
            ip = ip_ints(i,ia)
	    val_out(:,ip) = val_domain(:,i,ia)
          end do
        end do

	end subroutine shympi_copy_3d_i

!*******************************

	subroutine shympi_copy_3d_r(val_domain,nov,nous,val_out &
     &				,ndomains,nmax,ip_ints)

	integer nov,nous
	real val_domain(nov,nn_max,n_threads)
	real val_out(nov,nous)
	integer ndomains(n_threads)
	integer nmax
	integer ip_ints(nmax,n_threads)

	integer ia,i,n,ip

        do ia=1,n_threads
          n=ndomains(ia)
          do i=1,n
            ip = ip_ints(i,ia)
	    val_out(:,ip) = val_domain(:,i,ia)
          end do
        end do

	end subroutine shympi_copy_3d_r

!*******************************

	subroutine shympi_copy_3d_d(val_domain,nov,nous,val_out &
     &				,ndomains,nmax,ip_ints)

	integer nov,nous
	double precision val_domain(nov,nn_max,n_threads)
	double precision val_out(nov,nous)
	integer ndomains(n_threads)
	integer nmax
	integer ip_ints(nmax,n_threads)

	integer ia,i,n,ip

        do ia=1,n_threads
          n=ndomains(ia)
          do i=1,n
            ip = ip_ints(i,ia)
	    val_out(:,ip) = val_domain(:,i,ia)
          end do
        end do

	end subroutine shympi_copy_3d_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_getvals_2d_node_r(kind,vals,val)

	use basin
	use levels

	integer kind(2)
	real vals(nkn)
	real val(1)

	if( bmpi_skip ) then
	  val = vals(kind(1))
	  return
	end if

	call shympi_getvals_internal_r(kind,1,nkn &
     &                                    ,vals,val)

	end subroutine shympi_getvals_2d_node_r

!*******************************

	subroutine shympi_getvals_2d_node_i(kind,vals,val)

	use basin
	use levels

	integer kind(2)
	integer vals(nkn)
	integer val(1)

	if( bmpi_skip ) then
	  val = vals(kind(1))
	  return
	end if

	call shympi_getvals_internal_i(kind,1,nkn &
     &                                    ,vals,val)

	end subroutine shympi_getvals_2d_node_i

!*******************************

	subroutine shympi_getvals_3d_node_r(kind,vals,val)

	use basin
	use levels

	integer kind(2)
	real vals(nlvdi,nkn)
	real val(nlvdi)

	if( bmpi_skip ) then
	  val(:) = vals(:,kind(1))
	  return
	end if

	call shympi_getvals_internal_r(kind,nlvdi,nkn &
     &                                    ,vals,val)

	end subroutine shympi_getvals_3d_node_r

!*******************************

	subroutine shympi_getvals_3d_node_i(kind,vals,val)

	use basin
	use levels

	integer kind(2)
	integer vals(nlvdi,nkn)
	integer val(nlvdi)

	if( bmpi_skip ) then
	  val(:) = vals(:,kind(1))
	  return
	end if

	call shympi_getvals_internal_i(kind,nlvdi,nkn &
     &                                    ,vals,val)

	end subroutine shympi_getvals_3d_node_i

!******************************************************************
!******************************************************************
!******************************************************************

	function shympi_min_d(vals)

	double precision shympi_min_d
	double precision vals(:)
	double precision val

	val = MINVAL(vals)
	if( bmpi ) call shympi_reduce_internal_d('min',val)

	shympi_min_d = val

	end function shympi_min_d

!******************************************************************

	function shympi_min_r(vals)

	real shympi_min_r
	real vals(:)
	real val

	val = MINVAL(vals)
	if( bmpi ) call shympi_reduce_internal_r('min',val)

	shympi_min_r = val

	end function shympi_min_r

!******************************************************************

	function shympi_min_i(vals)

	integer shympi_min_i
	integer vals(:)
	integer val

	val = MINVAL(vals)
	if( bmpi ) call shympi_reduce_internal_i('min',val)

	shympi_min_i = val

	end function shympi_min_i

!******************************************************************

	function shympi_min_0_d(val0)

	double precision shympi_min_0_d
	double precision val0
	double precision val

	val = val0
	if( bmpi ) call shympi_reduce_internal_d('min',val)

	shympi_min_0_d = val

	end function shympi_min_0_d

!******************************************************************

	function shympi_min_0_r(val0)

	real shympi_min_0_r
	real val0
	real val

	val = val0
	if( bmpi ) call shympi_reduce_internal_r('min',val)

	shympi_min_0_r = val

	end function shympi_min_0_r

!******************************************************************

	function shympi_min_0_i(val0)

	integer shympi_min_0_i
	integer val0
	integer val

	val = val0
	if( bmpi ) call shympi_reduce_internal_i('min',val)

	shympi_min_0_i = val

	end function shympi_min_0_i

!******************************************************************

	function shympi_max_d(vals)

	double precision shympi_max_d
	double precision vals(:)
	double precision val

	val = MAXVAL(vals)
	if( bmpi ) call shympi_reduce_internal_d('max',val)

	shympi_max_d = val

	end function shympi_max_d

!******************************************************************

	function shympi_max_r(vals)

	real shympi_max_r
	real vals(:)
	real val

	val = MAXVAL(vals)
	if( bmpi ) call shympi_reduce_internal_r('max',val)

	shympi_max_r = val

	end function shympi_max_r

!******************************************************************

	function shympi_max_i(vals)

	integer shympi_max_i
	integer vals(:)
	integer val

	val = MAXVAL(vals)
	if( bmpi ) call shympi_reduce_internal_i('max',val)

	shympi_max_i = val

	end function shympi_max_i

!******************************************************************

	function shympi_max_0_i(val0)

! routine for val that is scalar

	integer shympi_max_0_i
	integer val0
	integer val

	val = val0
	if( bmpi ) call shympi_reduce_internal_i('max',val)

	shympi_max_0_i = val

	end function shympi_max_0_i

!******************************************************************

	function shympi_max_0_r(val0)

! routine for val that is scalar

	real shympi_max_0_r
	real val0
	real val

	val = val0
	if( bmpi ) call shympi_reduce_internal_r('max',val)

	shympi_max_0_r = val

	end function shympi_max_0_r

!******************************************************************

	function shympi_max_0_d(val0)

! routine for val that is scalar

	double precision shympi_max_0_d
	double precision val0
	double precision val

	val = val0
	if( bmpi ) call shympi_reduce_internal_d('max',val)

	shympi_max_0_d = val

	end function shympi_max_0_d

!******************************************************************

	function shympi_sum_d(vals)

	double precision shympi_sum_d
	double precision vals(:)
	double precision val

	val = SUM(vals)
	if( bmpi ) call shympi_reduce_internal_d('sum',val)

	shympi_sum_d = val

	end function shympi_sum_d

!******************************************************************

	function shympi_sum_r(vals)

	real shympi_sum_r
	real vals(:)
	real val

	val = SUM(vals)
	if( bmpi ) call shympi_reduce_internal_r('sum',val)

	shympi_sum_r = val

	end function shympi_sum_r

!******************************************************************

	function shympi_sum_i(vals)

	integer shympi_sum_i
	integer vals(:)
	integer val

	val = SUM(vals)
	if( bmpi ) call shympi_reduce_internal_i('sum',val)

	shympi_sum_i = val

	end function shympi_sum_i

!******************************************************************

	function shympi_sum_0_d(val0)

	double precision shympi_sum_0_d
	double precision val0
	double precision val

	val = val0
	if( bmpi ) call shympi_reduce_internal_d('sum',val)

	shympi_sum_0_d = val

	end function shympi_sum_0_d

!******************************************************************

	function shympi_sum_0_r(val0)

	real shympi_sum_0_r
	real val0
	real val

	val = val0
	if( bmpi ) call shympi_reduce_internal_r('sum',val)

	shympi_sum_0_r = val

	end function shympi_sum_0_r

!******************************************************************

	function shympi_sum_0_i(val0)

	integer shympi_sum_0_i
	integer val0
	integer val

	val = val0
	if( bmpi ) call shympi_reduce_internal_i('sum',val)

	shympi_sum_0_i = val

	end function shympi_sum_0_i

!******************************************************************

	subroutine shympi_array_reduce_i(what,vals)
	character*(*) what
	integer vals(:)
	if( bmpi ) then
	  call shympi_reduce_array_internal_i(what,size(vals),vals)
	end if
	end subroutine shympi_array_reduce_i

	subroutine shympi_array_reduce_r(what,vals)
	character*(*) what
	real vals(:)
	if( bmpi ) then
	  call shympi_reduce_array_internal_r(what,size(vals),vals)
	end if
	end subroutine shympi_array_reduce_r

	subroutine shympi_array_reduce_d(what,vals)
	character*(*) what
	double precision vals(:)
	if( bmpi ) then
	  call shympi_reduce_array_internal_d(what,size(vals),vals)
	end if
	end subroutine shympi_array_reduce_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_reduce_r(what,vals,val)

	character*(*) what
	real vals(:)
	real val

	if( what == 'min' ) then
	  val = MINVAL(vals)
	  if( bmpi ) call shympi_reduce_internal_r(what,val)
	else if( what == 'max' ) then
	  val = MAXVAL(vals)
	  if( bmpi ) call shympi_reduce_internal_r(what,val)
	else if( what == 'sum' ) then
	  val = SUM(vals)
	  if( bmpi ) call shympi_reduce_internal_r(what,val)
	else
	  write(6,*) 'what = ',what
	  stop 'error stop shympi_reduce_r: not ready'
	end if

	end subroutine shympi_reduce_r

!******************************************************************
!******************************************************************
!******************************************************************
! next are routines for partition on elements - empty here
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_and_sum_3d_nodes_r(a)
	real a(:,:)
	end subroutine shympi_exchange_and_sum_3d_nodes_r

	subroutine shympi_exchange_and_sum_3d_nodes_d(a)
	double precision a(:,:)
	end subroutine shympi_exchange_and_sum_3d_nodes_d

	subroutine shympi_exchange_and_sum_2d_nodes_r(a)
	real a(:)
	end subroutine shympi_exchange_and_sum_2d_nodes_r

	subroutine shympi_exchange_and_sum_2d_nodes_d(a)
	double precision a(:)
	end subroutine shympi_exchange_and_sum_2d_nodes_d

	subroutine shympi_exchange_2d_nodes_min_i(a)
	integer a(:)
	end subroutine shympi_exchange_2d_nodes_min_i

	subroutine shympi_exchange_2d_nodes_max_i(a)
	integer a(:)
	end subroutine shympi_exchange_2d_nodes_max_i

	subroutine shympi_exchange_2d_nodes_min_r(a)
	real a(:)
	end subroutine shympi_exchange_2d_nodes_min_r

	subroutine shympi_exchange_2d_nodes_max_r(a)
	real a(:)
	end subroutine shympi_exchange_2d_nodes_max_r

!******************************************************************
!******************************************************************
!******************************************************************

	function shympi_output()

	logical shympi_output

	shympi_output = my_id == 0

	end function shympi_output

!******************************************************************

	function shympi_is_master()

	logical shympi_is_master

	shympi_is_master = my_id == 0

	end function shympi_is_master

!******************************************************************

	subroutine shympi_comment(text)

	character*(*), intent(in) :: text

	if( bcomment ) then
	  if( bmpi_master ) then
	    write(iuc,*) indent,trim(text)
	  end if
	end if

	end subroutine shympi_comment

!******************************************************************

	subroutine shympi_comment_main(text)

	character*(*), intent(in) :: text

	if( bcomment ) then
	  if( bmpi_master ) then
	    write(iuc,*) trim(text)
	  end if
	end if

	end subroutine shympi_comment_main

!******************************************************************

        subroutine shympi_parallel_code(text)

        character*(*) text

        text = 'node'

        end subroutine shympi_parallel_code

!******************************************************************

	function shympi_is_parallel()

	logical shympi_is_parallel

	shympi_is_parallel = bmpi

	end function shympi_is_parallel

!******************************************************************

	function shympi_can_parallel()

	logical shympi_can_parallel

	shympi_can_parallel = .true.

	end function shympi_can_parallel

!******************************************************************

        subroutine shympi_set_debug(bdebug)

        logical bdebug

        bmpi_debug = bdebug

        end subroutine shympi_set_debug

!******************************************************************

	subroutine shympi_univocal_nodes

	use basin, only : nkn

	implicit none

	integer ierr

	integer i

	univocal_nodes%numberID = nkn

	allocate(univocal_nodes%localID(nkn))

	do i=1,nkn
	  univocal_nodes%localID(i) = i
	end do

	end subroutine shympi_univocal_nodes

!******************************************************************

        subroutine check_part_basin(what)

        use basin

        implicit none

        character*(5) what
        integer pnkn,pnel,pn_threads,ierr
        integer control
        character*(14) filename
        character*(11) pwhat
        integer i

        call shympi_get_filename(filename,what)

        write(6,*) filename,what

        open(unit=108, file=filename, form="formatted" &
     &   , iostat=control, status="old", action="read")
        
        if(control .ne. 0) then
          if(my_id .eq. 0) then
            write(6,*)'error stop: partitioning file not found'
          end if
          call shympi_barrier
        stop
        end if

        read(unit=108, fmt="(i12,i12,i12,A12)") pnkn,pnel,pn_threads &
     &                  ,pwhat

        if(pnkn .ne. nkndi .or. pnel .ne. neldi .or. pn_threads      &
     &          .ne. n_threads .or. pwhat .ne. what) then
         if(my_id .eq. 0) then
          write(6,*)'basin file does not match'
          write(6,*)'partitioning file:nkn,nel,n_threads,partitioning' &
     &          ,pnkn,pnel,pn_threads,pwhat
          write(6,*)'basin in str file:nkn,nel,n_threads,partitioning' &
     &          ,nkndi,neldi,n_threads,what
         end if
         call shympi_barrier
         stop
        end if

        if(what .eq. 'nodes') then
          allocate(allPartAssign(nkndi))
          read(unit=108,fmt="(i12,i12,i12,i12,i12,i12)") &
     &          (allPartAssign(i),i=1,nkndi)
        else 
          write(6,*)'error partitioning file on nodes'
          stop
        end if

        close(108)

        return

        end subroutine check_part_basin        

!******************************************************************

	subroutine check_external_numbers

	use basin

	implicit none

	integer, allocatable :: ip(:)

	if( bmpi_skip ) return

	allocate(ip(nkn_global))
	call shympi_l2g_array(ipv,ip)
	if( any( ip_ext_node /= ip ) ) then
	  stop 'error stop check_external: node numbers'
	end if
	deallocate(ip)

	allocate(ip(nel_global))
	call shympi_l2g_array(ipev,ip)
	if( any( ip_ext_elem /= ip ) ) then
	  stop 'error stop check_external: elem numbers'
	end if
	deallocate(ip)

	if( shympi_is_master() ) then
	  write(6,*) 'successful check of external numbers'
	end if

	end subroutine check_external_numbers

!******************************************************************

        subroutine shympi_get_filename(filename,what)

          implicit none

          character*(5) what
          character*(14) filename,format_string

          if(n_threads .gt. 999) then
            format_string = "(A5,A5,I4)"
            write(filename,format_string)'part_',what,n_threads
          else if(n_threads .gt. 99) then
            format_string = "(A5,A5,I3)"
            write(filename,format_string)'part_',what,n_threads
          else if(n_threads .gt. 9) then
            format_string = "(A5,A5,I2)"
            write(filename,format_string)'part_',what,n_threads
          else
            format_string = "(A5,A5,I1)"
            write(filename,format_string)'part_',what,n_threads
          end if

          return

        end subroutine shympi_get_filename

!******************************************************************

	subroutine gassert(bassert,text)

	implicit none

	logical bassert
	character*(*) text

	real, save :: r = 0.

	if( bassert ) return

	write(6,*) 'assertion failed: ',trim(text)
	write(6,*) 1./r
	stop 'error stop gassert'

	end subroutine gassert

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine error_stop_2(routine,text)

	implicit none

	character*(*) routine,text

	write(6,*) 'error stop '//trim(routine)//': ',trim(text)
	flush(6)
	i_code = i_code_error
	call shympi_abort

	end subroutine error_stop_2

!******************************************************************

	subroutine error_stop_1(text)

	implicit none

	character*(*) text

	write(6,*) 'error stop: ',trim(text)
	flush(6)
	i_code = i_code_error
	call shympi_abort

	end subroutine error_stop_1

!******************************************************************

	subroutine error_stop_0

	implicit none

	write(6,*) 'error stop'
	flush(6)
	i_code = i_code_error
	call shympi_abort

	end subroutine error_stop_0

!******************************************************************

	subroutine error_stop_i0(ierr)

	implicit none

	integer ierr

	write(6,*) 'error stop'
	flush(6)
	i_code = ierr
	call shympi_abort

	end subroutine error_stop_i0

!******************************************************************

	subroutine success

	implicit none

	i_code = i_code_success
	call shympi_abort

	end subroutine success

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_bdebug_0(text)

	implicit none

	character*(*) text

	if( .not. blocal_shympi_debug ) return

        call shympi_barrier
        write(6,*) trim(text),my_id
        flush(6)
        call shympi_barrier

	end subroutine shympi_bdebug_0

!******************************************************************

	subroutine shympi_bdebug_3(text,i1,i2,i3)

	implicit none

	character*(*) text
	integer i1,i2,i3

	if( .not. blocal_shympi_debug ) return

        call shympi_barrier
        write(6,*) trim(text),i1,i2,i3,my_id
        flush(6)
        call shympi_barrier

	end subroutine shympi_bdebug_3

!==================================================================
        end module shympi
!==================================================================

