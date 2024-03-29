
!--------------------------------------------------------------------------
!
!    Copyright (C) 2001,2009  Georg Umgiesser
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

!**************************************************************************

c subroutine axis(xpage,ypage,ibcd,nchar,axlen,angle,firstv,deltav)
c subroutine scale(array,axlen,npts,inc)
c subroutine line(xarray,yarray,npts,inc,lintyp,inteq)
c subroutine number(xpage,ypage,height,fpn,angle,ndec)

c
c revision log :
c
c 30.11.2001	ggu	CW -> avoid compiler warning
c 03.04.2009	ggu	CW -> avoid compiler warning
c 10.11.2021	ggu	avoid compiler warnings for arithmetic if

!**************************************************************************

      subroutine axis(xpage,ypage,ibcd,nchar,axlen,angle,firstv,deltav)

	dimension ibcd(1)
	!data isym/4h*10 /

      character symb*40
      integer isym(40)          !CW -> avoid compiler warning
      equivalence (symb,isym(1))

      do k=1,10
        symb(k:k)=' '
      end do
c
      kn=nchar
      a=1.0
c set constants for annotation on cw or ccw side of axis
	  !      if (kn) 1,2,2
	  continue
	  IF (kn<0) then
	      goto 1
	  ELSE IF (kn==0) then
	      goto 2
	  ELSE IF (kn>0) then
	      goto 2
	  END IF
1     a=-a
      kn=-kn
2     ex=0.0
      adx= abs  (deltav)
	  !      if (adx) 3,7,3
	  continue
	  IF (adx<0) then
	      goto 3
	  ELSE IF (adx==0) then
	      goto 7
	  ELSE IF (adx>0) then
	      goto 3
	  END IF
	  !3     if (adx- 99.0) 6,4,4
 3	  continue
	  IF (adx- 99.0<0) then
	      goto 6
	  ELSE IF (adx- 99.0==0) then
	      goto 4
	  ELSE IF (adx- 99.0>0) then
	      goto 4
	  END IF
4     adx=adx/10.0
      ex=ex+1.0
      go to 3
5     adx=adx*10.0
      ex=ex-1.0
	  !6     if (adx-0.01) 5,7,7
 6	  continue
	  IF (adx-0.01<0) then
	      goto 5
	  ELSE IF (adx-0.01==0) then
	      goto 7
	  ELSE IF (adx-0.01>0) then
	      goto 7
	  END IF
7     xval=firstv*10.0**(-ex)
      adx=2.0*(deltav*10.0**(-ex))
      sth=angle*0.0174533
      cth=cos(sth)
      sth=sin(sth)
      cth2=cth+cth
      sth2=sth+sth
      dxb=-0.254
      dyb=0.38*a-0.127
      xn=xpage+dxb*cth-dyb*sth
      yn=ypage+dyb*cth+dxb*sth
      ntic=axlen+1.0
      ntc=(ntic+1)/2
      nt=ntc/2
      do  i=1,ntc				!GGU
      call number(xn,yn,0.24,xval,angle,2)
	  !      if(i-1) 10,10,105
	  continue
	  IF (i-1<0) then
	      goto 10
	  ELSE IF (i-1==0) then
	      goto 10
	  ELSE IF (i-1>0) then
	      goto 105
	  END IF
10    xn=xpage+2.0*dxb*cth-dyb*sth
      yn=ypage+2.0*dxb*sth+dyb*cth
105   xval=xval+adx
      xn=xn+cth2
      yn=yn+sth2
	  !      if (nt) 20,11,20
	  continue
	  IF (nt<0) then
	      goto 20
	  ELSE IF (nt==0) then
	      goto 11
	  ELSE IF (nt>0) then
	      goto 20
	  END IF
11    z=kn
	  !      if (ex)  12,13,12
	  continue
	  IF (ex<0) then
	      goto 12
	  ELSE IF (ex==0) then
	      goto 13
	  ELSE IF (ex>0) then
	      goto 12
	  END IF
12    z=z+7.0
13    dxb=-.175*z+axlen*0.5
      dyb=0.8*a-0.2
      xt=xpage+dxb*cth-dyb*sth
      yt=ypage+dyb*cth+dxb*sth
      ii = ibcd(1)
      call symbol(xt,yt,0.36,ii,angle,kn)
	  !      if (ex)  14,20,14
	  continue
	  IF (ex<0) then
	      goto 14
	  ELSE IF (ex==0) then
	      goto 20
	  ELSE IF (ex>0) then
	      goto 14
	  END IF
14    z=kn+2
      xt=xt+z*cth*0.36
      yt=yt+z*sth*0.36
      isym1 = isym(1)		!GGU
      call symbol(xt,yt,0.36,isym1,angle,3)
      xt=xt+(3.0*cth-0.6*sth)*0.36
      yt=yt+(3.0*sth+0.6*cth)*0.36
      call number(xt,yt,0.18,ex,angle,-1)
20    nt=nt-1
      end do	!GGU
      call plot(xpage+axlen*cth,ypage+axlen*sth,3)
      dxb=-0.178*a*sth
      dyb=+0.178*a*cth
      a=ntic-1
      xn=xpage+a*cth
      yn=ypage+a*sth
      do  30  i=1 , ntic
      call plot(xn,yn,2)
      call plot(xn+dxb,yn+dyb,2)
      call plot(xn,yn,2)
      xn=xn-cth
      yn=yn-sth
30    continue
      return
      end

!**************************************************************************

      subroutine scale(array,axlen,npts,inc)

      dimension  array(1),save(7)

      save(1)= 1.0
      save(2)= 2.0
      save(3)= 4.0
      save(4)= 5.0
      save(5)= 8.0
      save(6)=10.0
      save(7)=20.0
      fad=0.01
      k=iabs(inc)
      n=npts*k
      y0=array(1)
      yn=y0
      do  25  i=1 ,n,k
      ys=array(i)
	  !      if  (y0-ys)  22,22,21
	  continue
	  IF (y0-ys<0) then
	      goto 22
	  ELSE IF (y0-ys==0) then
	      goto 22
	  ELSE IF (y0-ys>0) then
	      goto 21
	  END IF
21    y0=ys
      go  to  25
	  !22    if  (ys-yn)  25,25,24
 22	  continue
	  IF (ys-yn<0) then
	      goto 25
	  ELSE IF (ys-yn==0) then
	      goto 25
	  ELSE IF (ys-yn>0) then
	      goto 24
	  END IF
24    yn=ys
25    continue
      firstv=y0
	  !      if  (y0)  34,35,35
	  continue
	  IF (y0<0) then
	      goto 34
	  ELSE IF (y0==0) then
	      goto 35
	  ELSE IF (y0>0) then
	      goto 35
	  END IF
34    fad=fad-1.0
35    deltav=(yn-firstv)/axlen
	  !      if (deltav) 70,70,40
	  continue
	  IF (deltav<0) then
	      goto 70
	  ELSE IF (deltav==0) then
	      goto 70
	  ELSE IF (deltav>0) then
	      goto 40
	  END IF
40    i=alog10(deltav) + 1000.0
      p=10.0**(i-1000)
      deltav=deltav/p-0.01
      do  45  i=1,6
      is=i
	  !      if  (save(i)-deltav)  45,50,50
	  continue
	  IF (save(i)-deltav<0) then
	      goto 45
	  ELSE IF (save(i)-deltav==0) then
	      goto 50
	  ELSE IF (save(i)-deltav>0) then
	      goto 50
	  END IF
45    continue
50    deltav=save(is)*p
      firstv=deltav*aint(y0/deltav+fad)
      t=firstv+(axlen+0.01)*deltav
	  !      if (t-yn)  55,57,57
	  continue
	  IF (t-yn<0) then
	      goto 55
	  ELSE IF (t-yn==0) then
	      goto 57
	  ELSE IF (t-yn>0) then
	      goto 57
	  END IF
55    firstv=p*aint(y0/p+fad)
      t=firstv+(axlen+.01)*deltav
	  !      if (t-yn)  56,57,57
	  continue
	  IF (t-yn<0) then
	      goto 56
	  ELSE IF (t-yn==0) then
	      goto 57
	  ELSE IF (t-yn>0) then
	      goto 57
	  END IF
56    is=is+1
      go  to  50
57    firstv=firstv-aint((axlen+(firstv-yn)/deltav)/2.0)*deltav
	  !      if (y0*firstv) 58,58,59
	  continue
	  IF (y0*firstv<0) then
	      goto 58
	  ELSE IF (y0*firstv==0) then
	      goto 58
	  ELSE IF (y0*firstv>0) then
	      goto 59
	  END IF
58    firstv=0.0
	  !59    if  (inc) 61,61,65
 59	  continue
	  IF (inc<0) then
	      goto 61
	  ELSE IF (inc==0) then
	      goto 61
	  ELSE IF (inc>0) then
	      goto 65
	  END IF
61    firstv=firstv+aint(axlen+.5)*deltav
      deltav=-deltav
65    n=n+1
      array(n)=firstv
      n=n+k
      array(n)=deltav
67    return
70    deltav=2.0*firstv
      deltav=abs(deltav/axlen)+1.
      go to 40

      end

!**************************************************************************

      subroutine line(xarray,yarray,npts,inc,lintyp,inteq)

      dimension xarray(1),yarray(1)

      lmin = npts*inc+1
      ldx  = lmin+inc
      nl   = lmin-inc
      firstx = xarray(lmin)
      deltax = xarray(ldx)
      firsty = yarray(lmin)
      deltay = yarray(ldx)
      call where (xn,yn,df)
      df=amax1(abs((xarray( 1)-firstx)/deltax-xn),
     1         abs((yarray( 1)-firsty)/deltay-yn) )
      dl=amax1(abs((xarray(nl)-firstx)/deltax-xn),
     1         abs((yarray(nl)-firsty)/deltay-yn) )

      ipen = 3
      icode = -1
      nt =iabs(lintyp)

      if(lintyp.eq.0) nt = 1
      if(df.gt.dl) then
        nf = nl
        na = ((npts-1)/nt)*nt+nt-(npts-1)
        kk = -inc
      else
        nf = 1
        na = nt
        kk = inc
      end if

      if(lintyp.lt.0) then
        ipena = 3
        icodea = -1
        lsw = 1
      else
        if(lintyp.eq.0) na = ldx
        ipena = 2
        icodea = -2
        lsw=0
      end if

      do i = 1,npts
        xn = (xarray(nf)-firstx)/deltax
        yn = (yarray(nf)-firsty)/deltay
        if(na.eq.nt) then
          call symbol (xn,yn,0.20,inteq,0.0,icode)
          na = 1
        else if(na.lt.nt) then
          if(lsw.eq.0) call plot (xn,yn,ipen)
          na = na + 1
        else
          call plot (xn,yn,ipen)
          na = na + 1
        end if
        nf = nf+kk
        icode = icodea
        ipen = ipena
      end do

      return
      end

!**************************************************************************

      subroutine number(xpage,ypage,height,fpn,angle,ndec)

      character num*20,iar*10,minus*1,ipoint*1,blank*1
      integer inum(20)          !ggu CW -> avoid compiler warning
      equivalence (num,inum(1)) !ggu CW -> avoid compiler warning
      data iar/'0123456789'/
      data minus/'-'/,ipoint/'.'/,blank/' '/

      do k=1,20
        num(k:k)=blank
      end do

      ii=0
      fpv = fpn
      n = ndec
      maxn = 9

      if(n.gt.maxn) n=maxn
      if(-n.gt.maxn) n=-maxn

      if(fpv.lt.0) then
        ii=ii+1
        num(ii:ii)=minus
      end if

      mn = -n
      if(n.lt.0) mn = mn - 1
      fpv = abs(fpv) + (0.5 * 10. ** mn)
      i = alog10(fpv) + 1.0

      ilp = i
      if(n.lt.-1) ilp = ilp + n + 1

      if(ilp.le.0) then
        ii=ii+1
        num(ii:ii)=iar(1:1)
      else
        if(ilp+n.gt.18) then
          n = -1
          if(ilp.gt.19) ilp = 19
        end if

        do j=1,ilp
          k = fpv * 10. ** (j - i)
          if(k.gt.9) k = 9
          ii=ii+1
          isb=k+1
          num(ii:ii)=iar(isb:isb)
          fpv = fpv - (float(k) * 10. ** (i - j))
        end do
      end if

      if(n.ge.0) then
        ii=ii+1
        num(ii:ii)=ipoint
        if(n.gt.0) then
          do j=1,n
            k = fpv * 10.
            if(k.gt.9) k = 9
            ii=ii+1
            isb=k+1
            num(ii:ii)=iar(isb:isb)
            fpv = fpv * 10. - float(k)
          end do
        end if
      end if

c     nchar=ii+1000
      nchar=ii
      inum1 = inum(1)		!GGU

      !call symbol (xpage,ypage,height,num,angle,nchar)
      call symbol (xpage,ypage,height,inum1,angle,nchar) !ggu CW

      return
      end

!**************************************************************************

