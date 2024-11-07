
!--------------------------------------------------------------
!
! dummy program not used anymore with new connection framework
!
! update_ielt() should be transfered to lnku.f90
!
!--------------------------------------------------------------

!****************************************************************

        subroutine make_links_old(nkn,nel,nen3v)
	implicit none
	integer nkn,nel,nen3v(3,nel)
	end

        subroutine make_links(nkn,nel,nen3v,ibound,kerr)
	implicit none
	integer nkn,nel,nen3v(3,nel),ibound(nkn),kerr
	end

        subroutine checklenk(nkn,ilinkv,lenkv,nen3v)
	implicit none
	integer nkn,ilinkv(*),lenkv(*),nen3v(3,*)
	end

        subroutine checklink(nkn,ilinkv,linkv)
	implicit none
	integer nkn,ilinkv(*),linkv(*)
	end

        subroutine checkkant(nkn,kantv)
	implicit none
	integer nkn,kantv(2,nkn)
	end

        subroutine checkielt(nel,ieltv,nen3v)
	implicit none
	integer nel,ieltv(3,nel),nen3v(3,nel)
	end

	subroutine link_set_write(bset)
	implicit none
	logical bset
	end

	subroutine link_set_stop(bset)
	implicit none
	logical bset
	end

!****************************************************************

        subroutine update_ielt(nkn,nel,ibound,ieltv,nen3v)

! updates vector ieltv with open boundary nodes

        implicit none

        integer nkn
        integer nel
        integer ibound(nkn)       ! >0 => open boundary node
        integer ieltv(3,nel)
        integer nen3v(3,nel)

        integer k,ie,ii,i
        integer ksnext,isbhnd,ksthis

        do ie=1,nel
          do ii=1,3
            k = ksthis(ii,ie,nen3v)
            if( ibound(k) .gt. 0 .and. ibound(ksnext(k,ie,nen3v)) .gt. 0 ) then
              i = isbhnd(k,ie,nen3v)
              ieltv(i,ie) = -1
            end if
          end do
        end do

        end

!****************************************************************

