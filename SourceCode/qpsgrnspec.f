      subroutine qpsgrnspec(ig)
      implicit none
      integer ig
c
      include 'qpsglobal.h'
c
      integer i,nn,il,istp,ly,lf,ldeg,ldeg0,ldegup,istate
      double precision f,dll,expo
      double precision fac,depst1,depst2,dys2,dxs,rr0a
      double precision disk(0:ldegmax),fdisk(0:ldegmax)
      double complex ca,cb,cag,cs1,cs2,cs3,cs4,ct1,ct2,cll1
      double complex cmuer,cmues,cksir,cksis,ckapr,ckaps
      double complex crrr,crrs,ys3d,yt1d,cgr,csta,cstb
      double complex ys1(4,0:2),ys2(4,0:2),ys3(4,0:2)
      double complex ys4(4,0:2),ys5(4,0:2),ys6(4,0:2)
      double complex yt1(4,0:2),yt2(4,0:2)
      double complex ypsv(6,4),ypsvg(6,4),ysh(2,2)
      double complex ysgr(4),ysgd(4),yspp(4),ystt(4),yttt(4)
c
      double precision expos
      double complex c1,c2,c3,c4
      data expos/12.d0/
      data c1,c2,c3,c4/(1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0),(4.d0,0.d0)/
c
c     Initiation
c
      rr0a=dmin1(rr0,0.1d0*rrup(lys))
c
      if(rr0a.gt.0.d0)then
c
c       tappered disk source:
c       2*[cos(theta)-cos(alpha)]/[1-cos(alpha)]^2
c
        dxs=dcos(rr0a/rrup(lys))
        fdisk(0)=1.d0
        fdisk(1)=(3.d0+dxs)/4.d0
        do ldeg=2,ldegmax
          fdisk(ldeg)=(dble(2*ldeg-1)*dxs*fdisk(ldeg-1)
     &               +dble(4-ldeg)*fdisk(ldeg-2))
     &               /dble(ldeg+3)
        enddo
        do ldeg=0,ldegmax
          disk(ldeg)=fdisk(ldeg)*dble(2*ldeg+1)/(4.d0*PI*rrup(lys)**2)
        enddo
      else
        do ldeg=0,ldegmax
          fdisk(ldeg)=1.d0
          disk(ldeg)=dble(2*ldeg+1)/(4.d0*PI*rrup(lys)**2)
        enddo
      endif
c
      do ldeg=0,ldegmax
        fdisk(ldeg)=dabs(fdisk(ldeg))
        do il=ldeg+1,ldegmax
          fdisk(ldeg)=dmax1(fdisk(ldeg),dabs(fdisk(il)))
        enddo
      enddo
c
      open(21,file=uzgrnfile(ig),
     &     form='unformatted',status='unknown')
      open(22,file=urgrnfile(ig),
     &     form='unformatted',status='unknown')
      open(23,file=utgrnfile(ig),
     &     form='unformatted',status='unknown')
      open(24,file=ppgrnfile(ig),
     &     form='unformatted',status='unknown')
      open(25,file=grgrnfile(ig),
     &     form='unformatted',status='unknown')
      open(26,file=gdgrnfile(ig),
     &     form='unformatted',status='unknown')
      open(27,file=trgrnfile(ig),
     &     form='unformatted',status='unknown')
      open(28,file=ttgrnfile(ig),
     &     form='unformatted',status='unknown')
c
      ldegup=1+ndmax
      do ldegup=1+ndmax,ldegcut
        expo=-dlog(fdisk(ldegup))
        dll=dsqrt(dble(ldegup)*dble(ldegup-1))
        do ly=min0(lys,lyr),max0(lys,lyr)-1
          expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
        enddo
        if(expo.gt.expos)goto 10
      enddo
10    continue
c
      if(ldegup.ge.ldegmax-ndmax-1)then
        print *,' Warning in qpsgrnspec: cutoff degree = ',ldegup,
     &          ' exceeds the limit!'
        ldegup=ldegmax-ndmax-1
      endif
c
      if(ldeggr.gt.0)then
        do ly=1,ly0
          do ldeg=0,ldegup
            nruku(ldeg,ly)=10
          enddo
        enddo
      endif
c
      write(*,*)' ==================================================='
      write(*,'(a,i3,a,f7.2,a)')'  Green functions for ',
     &        ig,'. source at depth ',grndep(ig)/KM2M,' km'
	write(*,*)' cutoff harmonic degree   discrete frequencies'
      write(*,'(i12,i21,a)')ldegup,nf,' + 2'
      write(*,*)' ==================================================='
c
      write(21)nt,dt,nf,df,ldegup
      write(22)nt,dt,nf,df,ldegup
      write(23)nt,dt,nf,df,ldegup
      write(24)nt,dt,nf,df,ldegup
      write(25)nt,dt,nf,df,ldegup
      write(26)nt,dt,nf,df,ldegup
      write(27)nt,dt,nf,df,ldegup
      write(28)nt,dt,nf,df,ldegup
c
      do istp=1,4
        do ldeg=0,1
          yuz(ldeg,istp,0)=(0.d0,0.d0)
          yur(ldeg,istp,0)=(0.d0,0.d0)
          yut(ldeg,istp,0)=(0.d0,0.d0)
          ypp(ldeg,istp,0)=(0.d0,0.d0)
          ygr(ldeg,istp,0)=(0.d0,0.d0)
          ygd(ldeg,istp,0)=(0.d0,0.d0)
          ytr(ldeg,istp,0)=(0.d0,0.d0)
          ytt(ldeg,istp,0)=(0.d0,0.d0)
        enddo
      enddo
      if(.not.nogravity)then
        do ly=1,ly0
          do ldeg=0,ldegup
            nruku(ldeg,ly)=0
          enddo
        enddo
      endif
c
      do ldeg=0,ldegup+1
        dll=dsqrt(dble(ldeg)*dble(ldeg+1))
c
c       determine degree dependent starting layer number
c       of sh solution
c
        lyupt(ldeg)=1
        expo=0.d0
        do ly=lys-1,1,-1
          expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
          if(expo.gt.expos)then
            lyupt(ldeg)=ly
            goto 101
          endif
        enddo
101     continue
c
        lylwt(ldeg)=min0(lycm,ly0)
        expo=0.d0
        do ly=lys,min0(lycm,ly0)-1
          expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
          if(expo.gt.expos)then
            lylwt(ldeg)=ly+1
            goto 102
          endif
        enddo
102     continue
c
c       determine degree dependent starting layer number
c       of psv solution
c
        lyuppsv(ldeg)=1
        expo=0.d0
        do ly=lys-1,1,-1
          expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
          if(expo.gt.expos)then
            lyuppsv(ldeg)=ly
            goto 201
          endif
        enddo
201     continue
c
        lylwpsv(ldeg)=ly0
        expo=0.d0
        do ly=lys,ly0-1
          expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
          if(expo.gt.expos)then
            lylwpsv(ldeg)=ly+1
            goto 202
          endif
        enddo
202     continue
      enddo
c
      ly=lylwpsv(0)
      depst1=(rearth-rrup(ly))/KM2M
      ly=lylwpsv(ldegup+1)
      depst2=(rearth-rrup(ly))/KM2M
c
c     determine layer dependent max. harmonic degree
c     of sh solution
c
      do ly=1,min0(lycm-1,ly0)
        ldegsh(ly)=1
        do ldeg=1,ldegup+1
          if(ly.ge.lyupt(ldeg).and.ly.le.lylwt(ldeg))then
            ldegsh(ly)=ldeg
          endif
        enddo
      enddo
c
c     determine layer dependent max. harmonic degree
c     of psv solution
c
      do ly=1,ly0
        ldegpsv(ly)=0
        do ldeg=0,ldegup+1
          if(ly.ge.min0(lyuppsv(ldeg),lyuppsv(ldeg)).and.
     &       ly.le.max0(lylwpsv(ldeg),lylwpsv(ldeg)))then
            ldegpsv(ly)=ldeg
          endif
        enddo
      enddo
c
      do lf=-1,nf
        f=dmax1(0.d0,dble(lf-1)*df)
c
c       lf = -1: for the fully relaxed case
c          =  0: for the unrelaxed (elastic) case
c
        istate=min0(lf,1)
        call qpsvemodel(f,istate)
c
        do il=0,2
          do istp=1,4
            ys1(istp,il)=(0.d0,0.d0)
            ys2(istp,il)=(0.d0,0.d0)
            ys3(istp,il)=(0.d0,0.d0)
            ys4(istp,il)=(0.d0,0.d0)
            ys5(istp,il)=(0.d0,0.d0)
            ys6(istp,il)=(0.d0,0.d0)
            yt1(istp,il)=(0.d0,0.d0)
            yt2(istp,il)=(0.d0,0.d0)
          enddo
        enddo
c
        do ldeg=0,ldegup+1
          dll=dble(ldeg)*dble(ldeg+1)
          do istp=1,4
            ys1(istp,0)=ys1(istp,1)
            ys2(istp,0)=ys2(istp,1)
            ys3(istp,0)=ys3(istp,1)
            ys4(istp,0)=ys4(istp,1)
            ys5(istp,0)=ys5(istp,1)
            ys6(istp,0)=ys6(istp,1)
c
            ys1(istp,1)=ys1(istp,2)
            ys2(istp,1)=ys2(istp,2)
            ys3(istp,1)=ys3(istp,2)
            ys4(istp,1)=ys4(istp,2)
            ys5(istp,1)=ys5(istp,2)
            ys6(istp,1)=ys6(istp,2)
c
            yt1(istp,0)=yt1(istp,1)
            yt2(istp,0)=yt2(istp,1)
c
            yt1(istp,1)=yt1(istp,2)
            yt2(istp,1)=yt2(istp,2)
          enddo
c
          if(lylwpsv(ldeg).lt.max0(lyr,lys).or.
     &       lyuppsv(ldeg).gt.min0(lyr,lys))then
            do istp=1,4
              do i=1,6
                ypsv(i,istp)=(0.d0,0.d0)
              enddo
            enddo
          else if(nogravity)then
            call qpspsvkern(ldeg,ypsv)
          else
            fac=dble(ldeg)/dble(ldeggr)
            if(fac.le.1.d0)then
              call qpspsvkerng(ldeg,ypsv)
            else if(fac.ge.1.d0+FLTAPER)then
              call qpspsvkern(ldeg,ypsv)
            else
              call qpspsvkerng(ldeg,ypsvg)
              call qpspsvkern(ldeg,ypsv)
              ca=dcmplx(dsin(0.5d0*PI*(fac-1.d0)/FLTAPER)**2,0.d0)
              cb=c1-ca
              do istp=1,4
                do i=1,6
                  ypsv(i,istp)=ca*ypsv(i,istp)+cb*ypsvg(i,istp)
                enddo
              enddo
            endif
          endif
c
          if(lylwt(ldeg).lt.max0(lyr,lys).or.
     &       lyupt(ldeg).gt.min0(lyr,lys))then
            do istp=1,2
              do i=1,2
                ysh(i,istp)=(0.d0,0.d0)
              enddo
            enddo
          else
            call qpsshkern(ldeg,ysh)
          endif
c
          crrr=dcmplx(rrup(lyr),0.d0)
          cgr=dcmplx(grup(lyr),0.d0)
          if(lyr.le.1)then
            cmuer=dcmplx(mueup(lyr),0.d0)*crexup(lyr)
            ckapr=dcmplx(kapup(lyr),0.d0)
            cksir=ckapr+cmuer*c4/c3
            do istp=1,4
              yspp(istp)=-ypsv(2,istp)
            enddo
          else
            cmuer=dcmplx(muelw(lyr-1),0.d0)*crexlw(lyr-1)
            ckapr=dcmplx(kaplw(lyr-1),0.d0)
            cksir=ckapr+cmuer*c4/c3
          endif
c
c         volume strain
c
          cll1=dcmplx(dble(ldeg)*dble(ldeg+1),0.d0)
          do istp=1,4
            yspp(istp)=(ypsv(2,istp)+(c4*ypsv(1,istp)
     &                -c2*cll1*ypsv(3,istp))*cmuer/crrr)/cksir
          enddo
c
c         gravity and geoid
c
          ca=dcmplx(dble(ldeg+1)/rrup(lyr),0.d0)
          cb=dcmplx(freeairgrd,0.d0)
          if(lyr.gt.1)cb=cb-dcmplx(4.d0*PI*BIGG*rholw(lyr-1),0.d0)
          do istp=1,4
            ysgr(istp)=-ypsv(6,istp)+ca*ypsv(5,istp)+cb*ypsv(1,istp)
            ysgd(istp)=ypsv(5,istp)/cgr
          enddo
c
c         spheroidal contribution to tilt
c
          do istp=1,4
            if(cdabs(cmuer).gt.0.d0)then
              ys3d=(-ypsv(1,istp)+ypsv(3,istp))/crrr+ypsv(4,istp)/cmuer
            else
              ys3d=(0.d0,0.d0)
            endif
            ystt(istp)=ypsv(5,istp)/(cgr*crrr)-ys3d
          enddo
c
c         toroidal contribution to tilt
c
          if(cdabs(cmuer).gt.0.d0)then
            do istp=1,2
              yt1d=ysh(1,istp)/crrr+ysh(2,istp)/cmuer
              yttt(istp)=-yt1d
            enddo
          else
            do istp=1,2
              yttt(istp)=(0.d0,0.d0)
            enddo
          endif
c
          crrs=dcmplx(rrup(lys),0.d0)
          cmues=dcmplx(mueup(lys),0.d0)*crexup(lys)
          cksis=dcmplx(kapup(lys),0.d0)+cmues*c4/c3
c
          if(istype.eq.0)then
c
c           for constant dislocation
c
            csta=cksis/dcmplx(kapup(lys)+4.d0*mueup(lys)/3.d0,0.d0)
            cstb=crexup(lys)
          else
c
c           for constant stress drop
c
            cstb=crexup(lys)
     &          /(dcmplx(1.d0-rexfault,0.d0)*crexup(lys)
     &           +dcmplx(rexfault,0.d0))
            csta=cksis/(dcmplx(kapup(lys),0.d0)+cmues*cstb*c4/c3)
          endif
c
c         1. Explosion (M11=M22=M33=1)
c
          cs1=csta*dcmplx( disk(ldeg),0.d0)/cksis
          cs2=csta*dcmplx(-disk(ldeg)*4.d0/rrup(lys),0.d0)*cmues/cksis
          cs4=csta*dcmplx( disk(ldeg)*2.d0/rrup(lys),0.d0)*cmues/cksis
          ys1(1,2)=cs1*ypsv(1,1)+cs2*ypsv(1,2)+cs4*ypsv(1,4)
          ys2(1,2)=cs1*ypsv(3,1)+cs2*ypsv(3,2)+cs4*ypsv(3,4)
          ys3(1,2)=cs1*yspp(1)+cs2*yspp(2)+cs4*yspp(4)
          ys4(1,2)=cs1*ysgr(1)+cs2*ysgr(2)+cs4*ysgr(4)
          ys5(1,2)=cs1*ysgd(1)+cs2*ysgd(2)+cs4*ysgd(4)
          ys6(1,2)=cs1*ystt(1)+cs2*ystt(2)+cs4*ystt(4)
          yt1(1,2)=(0.d0,0.d0)
          yt2(1,2)=(0.d0,0.d0)
c
c         2. Strike-slip (M12=M21=1)
c
          if(ldeg.lt.2)then
            ys1(2,2)=(0.d0,0.d0)
            ys2(2,2)=(0.d0,0.d0)
            ys3(2,2)=(0.d0,0.d0)
            ys4(2,2)=(0.d0,0.d0)
            ys5(2,2)=(0.d0,0.d0)
            ys6(2,2)=(0.d0,0.d0)
            yt1(2,2)=(0.d0,0.d0)
            yt2(2,2)=(0.d0,0.d0)
          else
            ct2=cstb*dcmplx(disk(ldeg)/(dll*rrup(lys)),0.d0)
            cs4=-ct2
            ys1(2,2)=cs4*ypsv(1,4)
            ys2(2,2)=cs4*ypsv(3,4)
            ys3(2,2)=cs4*yspp(4)
            ys4(2,2)=cs4*ysgr(4)
            ys5(2,2)=cs4*ysgd(4)
            ys6(2,2)=cs4*ystt(4)
            yt1(2,2)=ct2*ysh(1,2)
            yt2(2,2)=ct2*yttt(2)
          endif
c
c         3. Dip-slip (M13=M31=1)
c
          if(ldeg.lt.1)then
            ys1(3,2)=(0.d0,0.d0)
            ys2(3,2)=(0.d0,0.d0)
            ys3(3,2)=(0.d0,0.d0)
            ys4(3,2)=(0.d0,0.d0)
            ys5(3,2)=(0.d0,0.d0)
            ys6(3,2)=(0.d0,0.d0)
            yt1(3,2)=(0.d0,0.d0)
            yt2(3,2)=(0.d0,0.d0)
          else
            ct1=cstb*dcmplx(disk(ldeg)/dll,0.d0)/cmues
            cs3=ct1
            ys1(3,2)=cs3*ypsv(1,3)
            ys2(3,2)=cs3*ypsv(3,3)
            ys3(3,2)=cs3*yspp(3)
            ys4(3,2)=cs3*ysgr(3)
            ys5(3,2)=cs3*ysgd(3)
            ys6(3,2)=cs3*ystt(3)
            yt1(3,2)=ct1*ysh(1,1)
            yt2(3,2)=ct1*yttt(1)
          endif
c
c         4. CLVD (M33=1,M11=M22=-0.5)
c
          cs1=cstb*dcmplx( disk(ldeg),0.d0)/cksis
          cs2=cstb*dcmplx( disk(ldeg)/rrup(lys),0.d0)
     &       *(c3-c4*cmues/cksis)
          cs4=-(0.5d0,0.d0)*cs2
          ys1(4,2)=cs1*ypsv(1,1)+cs2*ypsv(1,2)+cs4*ypsv(1,4)
          ys2(4,2)=cs1*ypsv(3,1)+cs2*ypsv(3,2)+cs4*ypsv(3,4)
          ys3(4,2)=cs1*yspp(1)+cs2*yspp(2)+cs4*yspp(4)
          ys4(4,2)=cs1*ysgr(1)+cs2*ysgr(2)+cs4*ysgr(4)
          ys5(4,2)=cs1*ysgd(1)+cs2*ysgd(2)+cs4*ysgd(4)
          ys6(4,2)=cs1*ystt(1)+cs2*ystt(2)+cs4*ystt(4)
          yt1(4,2)=(0.d0,0.d0)
          yt2(4,2)=(0.d0,0.d0)
c
c         convert to (r,t,p,g)-system
c         ===========================
c
c         1. Explosion (M11=M22=M33=1)
c            yr,yg normalized by Plm(l,0,cos(t))
c            yt normalized by Plm(l,1,cos(t)) = -dPlm/dt
c
          yuz(ldeg,1,0)=ys1(1,2)
          yur(ldeg,1,0)=-ys2(1,2)
          ypp(ldeg,1,0)=ys3(1,2)
          ygr(ldeg,1,0)=ys4(1,2)
          ygd(ldeg,1,0)=ys5(1,2)
          ytr(ldeg,1,0)=-ys6(1,2)
c
c         2. Strike-slip (M12=M21=1)
c            yr,yg normalized by Plm(l,2,cos(t))*sin(2p)
c            yt normalized by Plm(l,2,cos(t))*sin(2p)/sin(t)
c            yp normalized by Plm(l,2,cos(t))*cos(2p)/sin(t)
c
          yuz(ldeg,2,0)=ys1(2,2)
          if(ldeg.ge.3)then
            ca=dcmplx(dble(ldeg-2)*dble(ldeg-3)/dble(2*ldeg-3),0.d0)
            cb=dcmplx(dble(ldeg+1)*dble(ldeg+2)/dble(2*ldeg+1),0.d0)
            yur(ldeg-1,2,0)= ca*ys2(2,0)-cb*ys2(2,2)-c2*yt1(2,1)
            yut(ldeg-1,2,0)=-ca*yt1(2,0)+cb*yt1(2,2)+c2*ys2(2,1)
            ytr(ldeg-1,2,0)= ca*ys6(2,0)-cb*ys6(2,2)-c2*yt2(2,1)
            ytt(ldeg-1,2,0)=-ca*yt2(2,0)+cb*yt2(2,2)+c2*ys6(2,1)
          endif
          ypp(ldeg,2,0)=ys3(2,2)
          ygr(ldeg,2,0)=ys4(2,2)
          ygd(ldeg,2,0)=ys5(2,2)
c
c         3. Dip-slip (M13=M31=1)
c            yr,yg normalized by Plm(l,2,cos(t))*cos(p)
c            yt normalized by Plm(l,1,cos(t))*cos(p)/sin(t)
c            yp normalized by Plm(l,1,cos(t))*sin(p)/sin(t)
c
          yuz(ldeg,3,0)=ys1(3,2)
          if(ldeg.ge.2)then
            ca=dcmplx(dble(ldeg-2)*dble(ldeg-2)/dble(2*ldeg-3),0.d0)
            cb=dcmplx(dble(ldeg+1)*dble(ldeg+1)/dble(2*ldeg+1),0.d0)
            yur(ldeg-1,3,0)= ca*ys2(3,0)-cb*ys2(3,2)+yt1(3,1)
            yut(ldeg-1,3,0)=-ca*yt1(3,0)+cb*yt1(3,2)-ys2(3,1)
            ytr(ldeg-1,3,0)= ca*ys6(3,0)-cb*ys6(3,2)+yt2(3,1)
            ytt(ldeg-1,3,0)=-ca*yt2(3,0)+cb*yt2(3,2)-ys6(3,1)
          endif
          ypp(ldeg,3,0)=ys3(3,2)
          ygr(ldeg,3,0)=ys4(3,2)
          ygd(ldeg,3,0)=ys5(3,2)
c
c         4. CLVD (M33=1,M11=M22=-0.5)
c            yr,yg normalized by Plm(l,0,cos(t))
c            yt normalized by Plm(l,0,cos(t)) = -dPlm/dt
c
          yuz(ldeg,4,0)=ys1(4,2)
          yur(ldeg,4,0)=-ys2(4,2)
          ypp(ldeg,4,0)=ys3(4,2)
          ygr(ldeg,4,0)=ys4(4,2)
          ygd(ldeg,4,0)=ys5(4,2)
          ytr(ldeg,4,0)=-ys6(4,2)
        enddo
c
        if(lf.eq.-1)then
          write(*,'(a,i5,a,2(f7.2,a))')
     &       '  fully relaxed response: cut-off degree = ',ldegup,
     &       ', start depth = ',depst1,' - ',depst2,' km'
        else if(lf.eq.0)then
          write(*,'(a,i5,a,2(f7.2,a))')
     &       '     co-seismic response: cut-off degree = ',ldegup,
     &       ', start depth = ',depst1,' - ',depst2,' km'
        else
          write(*,'(i6,a,E14.6,a,i5,a,2(f7.2,a))')lf,'.',f,
     &       ' Hz: cut-off degree = ',ldegup,
     &       ', start depth = ',depst1,' - ',depst2,' km'
        endif
c
        write(21)ldegup
        write(22)ldegup
        write(23)ldegup
        write(24)ldegup
        write(25)ldegup
        write(26)ldegup
        write(27)ldegup
        write(28)ldegup
        write(21)((yuz(ldeg,istp,0),ldeg=0,ldegup),istp=1,4)
        write(22)((yur(ldeg,istp,0),ldeg=0,ldegup),istp=1,4)
        write(23)((yut(ldeg,istp,0),ldeg=0,ldegup),istp=2,3)
        write(24)((ypp(ldeg,istp,0),ldeg=0,ldegup),istp=1,4)
        write(25)((ygr(ldeg,istp,0),ldeg=0,ldegup),istp=1,4)
        write(26)((ygd(ldeg,istp,0),ldeg=0,ldegup),istp=1,4)
        write(27)((ytr(ldeg,istp,0),ldeg=0,ldegup),istp=1,4)
        write(28)((ytt(ldeg,istp,0),ldeg=0,ldegup),istp=2,3)
      enddo
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
      return
      end
