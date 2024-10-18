      module mod_bfm_statevars
        type :: bfm_statevars
          private
            integer,public  :: nstates
            integer,public  :: iO2o=0
            integer,public  :: iN1p=0
            integer,public  :: iN3n=0
            integer,public  :: iN4n=0
            integer,public  :: iO4n=0
            integer,public  :: iN5s=0
            integer,public  :: iN6r=0
            integer,public  :: iB1c=0
            integer,public  :: iB1n=0
            integer,public  :: iB1p=0
            integer,public  :: iP1c=0
            integer,public  :: iP1n=0
            integer,public  :: iP1p=0
            integer,public  :: iP1l=0
            integer,public  :: iP1s=0
            integer,public  :: iP2c=0
            integer,public  :: iP2n=0
            integer,public  :: iP2p=0
            integer,public  :: iP2l=0
            integer,public  :: iP3c=0
            integer,public  :: iP3n=0
            integer,public  :: iP3p=0
            integer,public  :: iP3l=0
            integer,public  :: iP4c=0
            integer,public  :: iP4n=0
            integer,public  :: iP4p=0
            integer,public  :: iP4l=0
            integer,public  :: iZ3c=0
            integer,public  :: iZ3n=0
            integer,public  :: iZ3p=0
            integer,public  :: iZ4c=0
            integer,public  :: iZ4n=0
            integer,public  :: iZ4p=0
            integer,public  :: iZ5c=0
            integer,public  :: iZ5n=0
            integer,public  :: iZ5p=0
            integer,public  :: iZ6c=0
            integer,public  :: iZ6n=0
            integer,public  :: iZ6p=0
            integer,public  :: iR1c=0
            integer,public  :: iR1n=0
            integer,public  :: iR1p=0
            integer,public  :: iR1l=0
            integer,public  :: iR2c=0
            integer,public  :: iR2l=0
            integer,public  :: iR3c=0
            integer,public  :: iR3l=0
            integer,public  :: iR6c=0
            integer,public  :: iR6n=0
            integer,public  :: iR6p=0
            integer,public  :: iR6s=0
            integer,public  :: iO3c=0
            integer,public  :: iO3h=0
            integer,public  :: iO5c=0
            integer,allocatable,public  :: pp(:)
            character(len=3),allocatable,public  :: id(:)
          contains
            procedure,public         :: init_bfm_statevars
            procedure,public         :: finalize_bfm_statevars
        end type bfm_statevars
      
        contains
  
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       
        subroutine init_bfm_statevars(self                             &
     &                      , var_list                                 )
#include "BFM_module_list.h"
           class(bfm_statevars),target :: self
           character(len=3), dimension(:), intent(in) :: var_list
           integer::state

           self%nstates = size(var_list)

           allocate(self%pp(self%nstates))
           allocate(self%id(self%nstates))
           do state=1,self%nstates 
             select case( trim(var_list(state)) )
               case('O2o') 
                           self%iO2o=state 
                           self%pp(state)=ppO2o 
                           self%id(state)='O2o'
               case('N1p') 
                           self%iN1p=state 
                           self%pp(state)=ppN1p 
                           self%id(state)='N1p'
               case('N3n') 
                           self%iN3n=state 
                           self%pp(state)=ppN3n 
                           self%id(state)='N3n'
               case('N4n') 
                           self%iN4n=state 
                           self%pp(state)=ppN4n 
                           self%id(state)='N4n'
               case('O4n') 
                           self%iO4n=state 
                           self%pp(state)=ppO4n 
                           self%id(state)='O4n'
               case('N5s') 
                           self%iN5s=state 
                           self%pp(state)=ppN5s 
                           self%id(state)='N5s'
               case('N6r') 
                           self%iN6r=state 
                           self%pp(state)=ppN6r 
                           self%id(state)='N6r'
               case('B1c') 
                           self%iB1c=state 
                           self%pp(state)=ppB1c 
                           self%id(state)='B1c'
               case('B1n') 
                           self%iB1n=state 
                           self%pp(state)=ppB1n 
                           self%id(state)='B1n'
               case('B1p') 
                           self%iB1p=state 
                           self%pp(state)=ppB1p 
                           self%id(state)='B1p'
               case('P1c') 
                           self%iP1c=state 
                           self%pp(state)=ppP1c 
                           self%id(state)='P1c'
               case('P1n') 
                           self%iP1n=state 
                           self%pp(state)=ppP1n 
                           self%id(state)='P1n'
               case('P1p') 
                           self%iP1p=state 
                           self%pp(state)=ppP1p 
                           self%id(state)='P1p'
               case('P1l') 
                           self%iP1l=state 
                           self%pp(state)=ppP1l 
                           self%id(state)='P1l'
               case('P1s') 
                           self%iP1s=state 
                           self%pp(state)=ppP1s 
                           self%id(state)='P1s'
               case('P2c') 
                           self%iP2c=state 
                           self%pp(state)=ppP2c 
                           self%id(state)='P2c'
               case('P2n') 
                           self%iP2n=state 
                           self%pp(state)=ppP2n 
                           self%id(state)='P2n'
               case('P2p') 
                           self%iP2p=state 
                           self%pp(state)=ppP2p 
                           self%id(state)='P2p'
               case('P2l') 
                           self%iP2l=state 
                           self%pp(state)=ppP2l 
                           self%id(state)='P2l'
               case('P3c') 
                           self%iP3c=state 
                           self%pp(state)=ppP3c 
                           self%id(state)='P3c'
               case('P3n') 
                           self%iP3n=state 
                           self%pp(state)=ppP3n 
                           self%id(state)='P3n'
               case('P3p') 
                           self%iP3p=state 
                           self%pp(state)=ppP3p 
                           self%id(state)='P3p'
               case('P3l') 
                           self%iP3l=state 
                           self%pp(state)=ppP3l 
                           self%id(state)='P3l'
               case('P4c') 
                           self%iP4c=state 
                           self%pp(state)=ppP4c 
                           self%id(state)='P4c'
               case('P4n') 
                           self%iP4n=state 
                           self%pp(state)=ppP4n 
                           self%id(state)='P4n'
               case('P4p') 
                           self%iP4p=state 
                           self%pp(state)=ppP4p 
                           self%id(state)='P4p'
               case('P4l') 
                           self%iP4l=state 
                           self%pp(state)=ppP4l 
                           self%id(state)='P4l'
               case('Z3c') 
                           self%iZ3c=state 
                           self%pp(state)=ppZ3c 
                           self%id(state)='Z3c'
               case('Z3n') 
                           self%iZ3n=state 
                           self%pp(state)=ppZ3n 
                           self%id(state)='Z3n'
               case('Z3p') 
                           self%iZ3p=state 
                           self%pp(state)=ppZ3p 
                           self%id(state)='Z3p'
               case('Z4c') 
                           self%iZ4c=state 
                           self%pp(state)=ppZ4c 
                           self%id(state)='Z4c'
               case('Z4n') 
                           self%iZ4n=state 
                           self%pp(state)=ppZ4n 
                           self%id(state)='Z4n'
               case('Z4p') 
                           self%iZ4p=state 
                           self%pp(state)=ppZ4p 
                           self%id(state)='Z4p'
               case('Z5c') 
                           self%iZ5c=state 
                           self%pp(state)=ppZ5c 
                           self%id(state)='Z5c'
               case('Z5n') 
                           self%iZ5n=state 
                           self%pp(state)=ppZ5n 
                           self%id(state)='Z5n'
               case('Z5p') 
                           self%iZ5p=state 
                           self%pp(state)=ppZ5p 
                           self%id(state)='Z5p'
               case('Z6c') 
                           self%iZ6c=state 
                           self%pp(state)=ppZ6c 
                           self%id(state)='Z6c'
               case('Z6n') 
                           self%iZ6n=state 
                           self%pp(state)=ppZ6n 
                           self%id(state)='Z6n'
               case('Z6p') 
                           self%iZ6p=state 
                           self%pp(state)=ppZ6p 
                           self%id(state)='Z6p'
               case('R1c') 
                           self%iR1c=state 
                           self%pp(state)=ppR1c 
                           self%id(state)='R1c'
               case('R1n') 
                           self%iR1n=state 
                           self%pp(state)=ppR1n 
                           self%id(state)='R1n'
               case('R1p') 
                           self%iR1p=state 
                           self%pp(state)=ppR1p 
                           self%id(state)='R1p'
               case('R1l') 
                           self%iR1l=state 
                           self%pp(state)=ppR1l 
                           self%id(state)='R1l'
               case('R2c') 
                           self%iR2c=state 
                           self%pp(state)=ppR2c 
                           self%id(state)='R2c'
               case('R2l') 
                           self%iR2l=state 
                           self%pp(state)=ppR2l 
                           self%id(state)='R2l'
               case('R3c') 
                           self%iR3c=state 
                           self%pp(state)=ppR3c 
                           self%id(state)='R3c'
               case('R3l') 
                           self%iR3l=state 
                           self%pp(state)=ppR3l 
                           self%id(state)='R3l'
               case('R6c') 
                           self%iR6c=state 
                           self%pp(state)=ppR6c 
                           self%id(state)='R6c'
               case('R6n') 
                           self%iR6n=state 
                           self%pp(state)=ppR6n 
                           self%id(state)='R6n'
               case('R6p') 
                           self%iR6p=state 
                           self%pp(state)=ppR6p 
                           self%id(state)='R6p'
               case('R6s') 
                           self%iR6s=state 
                           self%pp(state)=ppR6s 
                           self%id(state)='R6s'
               case('O3c') 
                           self%iO3c=state 
                           self%pp(state)=ppO3c 
                           self%id(state)='O3c'
               case('O3h') 
                           self%iO3h=state 
                           self%pp(state)=ppO3h 
                           self%id(state)='O3h'
               case('O5c') 
                           self%iO5c=state 
                           self%pp(state)=ppO5c 
                           self%id(state)='O5c'
             end select
           enddo

        end subroutine init_bfm_statevars

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        subroutine finalize_bfm_statevars(self)
          implicit none
            class(bfm_statevars),target :: self
            deallocate(self%pp)
            deallocate(self%id)
        end subroutine finalize_bfm_statevars
       
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       
! -----------------------------------------------------------------------------------
      end module mod_bfm_statevars
! -----------------------------------------------------------------------------------

