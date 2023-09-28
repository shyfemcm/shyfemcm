
! 15.11.2018	ccf	call to tide_vuf in do_befor

!----------------------------------------------------------------------
        module befor_after
!----------------------------------------------------------------------
        contains
!----------------------------------------------------------------------
!********************************************************************

	subroutine do_befor

! to do in time loop before time step

	implicit none

	include 'modules.h'

	call modules(M_BEFOR)

	call tide_vuf
        call tideforce       !tidal potential !ccf
	call adjust_chezy

end subroutine do_befor

!********************************************************************

	subroutine do_after

! to do in time loop after time step

	implicit none

	include 'modules.h'

	double precision dtime

	call modules(M_AFTER)

	call get_act_dtime(dtime)

!	call wrouta
	call wrousa
!	call wrexta
	!call wrflxa
	call wrvola(dtime)
	call wrboxa

        call resid
        call rmsvel

        call rst_write_restart

!        call tsmed
	call ts_shell

!	call wrnetcdf		!output in netcdf format - not supported

	call custom(dtime)

end subroutine do_after

!*******************************************************************
!*******************************************************************
!*******************************************************************

!----------------------------------------------------------------------
        end module befor_after
!----------------------------------------------------------------------
