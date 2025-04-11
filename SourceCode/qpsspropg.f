      subroutine qpsspropg(ypsv,ldeg,lyup,lylw)
      implicit none
c
c     calculation of speroidal response (with gravity)
c     ypsv(6,4): solution vector (complex)
c
      integer ldeg,lyup,lylw
      double complex ypsv(6,4)
c
      include 'qpsglobal.h'
c
c     work space
c
      integer i,j,j0,istp,ly,ily,nly,key
      double precision y4max,rr1,rr2,dlnr,h,f
      double complex cdet,alf,bet,cyabs,cyswap,ca,cb
      double complex y0(6,3),c(2)
      double complex yup(6,3),ylw(6,3),yupc(4,2),ylwc(4,2)
      double complex coef6(6,6),b6(6,4),coef4(4,4),b4(4,2)
      external qpsdifmatl,qpsdifmats
c
      if(lyup.eq.1)then
        do j=1,3
          do i=1,6
            yup(i,j)=(0.d0,0.d0)
          enddo
        enddo
        yup(1,1)=(1.d0,0.d0)
        yup(3,2)=(1.d0,0.d0)
        yup(5,3)=(1.d0,0.d0)
        if(ldeg.eq.1)then
          yup(6,3)=(3.d0,0.d0)
        endif
      else
        call qpssmat(ldeg,lyup,lyup)
        do j=1,3
          do i=1,6
            yup(i,j)=mas6x6(i,2*j,lyup)
          enddo
        enddo
      endif
      if(lyr.eq.lyup)call cmemcpy(yup,y0,18)
c
c===============================================================================
c
c     propagation from atmosphere/ocean bottom to source
c
      do ly=lyup,lys-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(0.5d0*h*dble(ldeg)/rrlw(ly))
        dlnr=dlog(rrlw(ly)/rrup(ly))/dble(nly)
        rr2=rrup(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrup(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+yup(i,1)*dconjg(yup(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            yup(i,1)=yup(i,1)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,6
            alf=alf+yup(i,2)*dconjg(yup(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            yup(i,2)=yup(i,2)-alf*yup(i,1)
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+yup(i,2)*dconjg(yup(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            yup(i,2)=yup(i,2)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          bet=(0.d0,0.d0)
          do i=1,6
            alf=alf+yup(i,3)*dconjg(yup(i,1))/cypnorm(i,ly)**2
            bet=bet+yup(i,3)*dconjg(yup(i,2))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            yup(i,3)=yup(i,3)-alf*yup(i,1)-bet*yup(i,2)
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)-alf*y0(i,1)-bet*y0(i,2)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+yup(i,3)*dconjg(yup(i,3))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            yup(i,3)=yup(i,3)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)*cyabs
            enddo
          endif
c
          call ruku(yup,6,3,ly,ldeg,qpsdifmats,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr-1)call cmemcpy(yup,y0,18)
      enddo
c
      do j=1,3
        yup(1,j)=yup(1,j)/dcmplx(rrup(lys),0.d0)
        yup(2,j)=yup(2,j)/dcmplx(rrup(lys)**2,0.d0)
        yup(3,j)=yup(3,j)/dcmplx(rrup(lys),0.d0)
        yup(4,j)=yup(4,j)/dcmplx(rrup(lys)**2,0.d0)
c       yup(5,j)=yup(5,j)
        yup(6,j)=yup(6,j)/dcmplx(rrup(lys),0.d0)
      enddo
c
c===============================================================================
c
c     propagation within inner core
c
      if(lylw.ge.lycc)then
c
c       lowest layer is within inner core
c
        call qpssmat(ldeg,lylw,lylw)
        do j=1,3
          do i=1,6
            ylw(i,j)=mas6x6(i,2*j-1,lylw)
          enddo
        enddo
c
        if(lylw.eq.lyr.and.lylw.gt.lys)then
          call cmemcpy(ylw,y0,18)
        endif
      endif
c
      do ly=lylw-1,lycc,-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(0.5d0*h*dble(ldeg)/rrlw(ly))
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,1)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,1)=ylw(i,1)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,6
            alf=alf+ylw(i,2)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            ylw(i,2)=ylw(i,2)-alf*ylw(i,1)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,2)*dconjg(ylw(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,2)=ylw(i,2)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          bet=(0.d0,0.d0)
          do i=1,6
            alf=alf+ylw(i,3)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
            bet=bet+ylw(i,3)*dconjg(ylw(i,2))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            ylw(i,3)=ylw(i,3)-alf*ylw(i,1)-bet*ylw(i,2)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)-alf*y0(i,1)-bet*y0(i,2)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,3)*dconjg(ylw(i,3))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,3)=ylw(i,3)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)*cyabs
            enddo
          endif
c
          call ruku(ylw,6,3,ly,ldeg,qpsdifmats,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,18)
      enddo
c
c===============================================================================
c
c     propagation within outer core
c
      if(lylw.ge.lycc)then
c
c       interface conditions: solid to liquid
c
        y4max=0.d0
        do j=1,3
          if(y4max.lt.cdabs(ylw(4,j)))then
            y4max=cdabs(ylw(4,j))
            j0=j
          endif
        enddo
        do i=1,6
          cyswap=ylw(i,j0)
          ylw(i,j0)=ylw(i,3)
          ylw(i,3)=cyswap
        enddo
        do j=1,2
          ylwc(1,j)=ylw(4,3)*ylw(1,j)-ylw(4,j)*ylw(1,3)
          ylwc(2,j)=ylw(4,3)*ylw(2,j)-ylw(4,j)*ylw(2,3)
          ylwc(3,j)=ylw(4,3)*ylw(5,j)-ylw(4,j)*ylw(5,3)
          ylwc(4,j)=ylw(4,3)*ylw(6,j)-ylw(4,j)*ylw(6,3)
        enddo
        if(lycc.le.lyr)then
          do i=1,6
            cyswap=y0(i,j0)
            y0(i,j0)=y0(i,3)
            y0(i,3)=cyswap
          enddo
          do j=1,2
            do i=1,6
              y0(i,j)=ylw(4,3)*y0(i,j)-ylw(4,j)*y0(i,3)
            enddo
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
        endif
c
c       y2 = Ut
c
        ca=dcmplx(rholw(lycc-1)*grlw(lycc-1)*rrlw(lycc-1),0.d0)
        cb=dcmplx(rholw(lycc-1)*rrlw(lycc-1)**2,0.d0)
        do j=1,2
          c(j)=ca*ylwc(1,j)-ylwc(2,j)-cb*ylwc(3,j)
        enddo
        do i=1,4
          ylwc(i,1)=c(2)*ylwc(i,1)-c(1)*ylwc(i,2)
        enddo
        ylwc(2,1)=(0.d0,0.d0)
        do i=1,4
          ylwc(i,2)=(0.d0,0.d0)
        enddo
        ylwc(2,2)=c(2)
        if(lycc.le.lyr)then
          do i=1,6
            y0(i,1)=c(2)*y0(i,1)-c(1)*y0(i,2)
            y0(i,2)=(0.d0,0.d0)
          enddo
        endif
      else if(lylw.ge.lycm)then
        call qpssmatc(ldeg,lylw,lylw)
        do j=1,2
          do i=1,4
            ylwc(i,j)=mas4x4(i,2*j-1,lylw)
          enddo
        enddo
        if(lylw.eq.lyr.and.lylw.gt.lys)then
          do j=1,3
            do i=1,6
              y0(i,j)=(0.d0,0.d0)
            enddo
          enddo
          ca=dcmplx(rhoup(lylw)*grup(lylw)*rrup(lylw),0.d0)
          cb=dcmplx(rhoup(lylw)*rrup(lylw)**2,0.d0)
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=ca*ylwc(1,j)-cb*ylwc(3,j)
            y0(3,j)=ylwc(2,j)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
        endif
      endif
c
      do ly=min0(lylw-1,lycc-1),lycm,-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(0.5d0*h*dble(ldeg)/rrlw(ly))
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+ylwc(i,1)*dconjg(ylwc(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            ylwc(i,1)=ylwc(i,1)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,4
            alf=alf+ylwc(i,2)*dconjg(ylwc(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,4
            ylwc(i,2)=ylwc(i,2)-alf*ylwc(i,1)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+ylwc(i,2)*dconjg(ylwc(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            ylwc(i,2)=ylwc(i,2)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          call ruku(ylwc,4,2,ly,ldeg,qpsdifmatl,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr.and.ly.gt.lys)then
          do j=1,3
            do i=1,6
              y0(i,j)=(0.d0,0.d0)
            enddo
          enddo
          ca=dcmplx(rhoup(ly)*grup(ly)*rrup(ly),0.d0)
          cb=dcmplx(rhoup(ly)*rrup(ly)**2,0.d0)
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=ca*ylwc(1,j)-cb*ylwc(3,j)
            y0(3,j)=ylwc(2,j)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
        endif
      enddo
c
c===============================================================================
c
c     propagation from core-mantle boundary to source or ocean bottom
c
      if(lylw.ge.lycm)then
c
c       interface conditions: liquid to solid
c
        ca=dcmplx(rhoup(lycm)*grup(lycm)*rrup(lycm),0.d0)
        cb=dcmplx(rhoup(lycm)*rrup(lycm)**2,0.d0)
        do j=1,2
          ylw(1,j)=ylwc(1,j)
          ylw(2,j)=ca*ylwc(1,j)-cb*ylwc(3,j)
          ylw(3,j)=(0.d0,0.d0)
          ylw(4,j)=(0.d0,0.d0)
          ylw(5,j)=ylwc(3,j)
          ylw(6,j)=ylwc(4,j)
        enddo
        ylw(1,3)=(0.d0,0.d0)
        ylw(2,3)=(0.d0,0.d0)
        ylw(3,3)=(1.d0,0.d0)
        ylw(4,3)=(0.d0,0.d0)
        ylw(5,3)=(0.d0,0.d0)
        ylw(6,3)=(0.d0,0.d0)
      else
        call qpssmat(ldeg,lylw,lylw)
        do j=1,3
          do i=1,6
            ylw(i,j)=mas6x6(i,2*j-1,lylw)
          enddo
        enddo
        if(lylw.eq.lyr.and.lylw.gt.lys)call cmemcpy(ylw,y0,18)
      endif
c
      do ly=min0(lylw-1,lycm-1),lys,-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(0.5d0*h*dble(ldeg)/rrlw(ly))
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,1)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,1)=ylw(i,1)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,6
            alf=alf+ylw(i,2)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            ylw(i,2)=ylw(i,2)-alf*ylw(i,1)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,2)*dconjg(ylw(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,2)=ylw(i,2)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          bet=(0.d0,0.d0)
          do i=1,6
            alf=alf+ylw(i,3)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
            bet=bet+ylw(i,3)*dconjg(ylw(i,2))/cypnorm(i,ly)**2
          enddo
          do i=1,6
            ylw(i,3)=ylw(i,3)-alf*ylw(i,1)-bet*ylw(i,2)
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)-alf*y0(i,1)-bet*y0(i,2)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,6
            cyabs=cyabs+ylw(i,3)*dconjg(ylw(i,3))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,6
            ylw(i,3)=ylw(i,3)*cyabs
          enddo
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,3)=y0(i,3)*cyabs
            enddo
          endif
c
          call ruku(ylw,6,3,ly,ldeg,qpsdifmats,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,18)
      enddo
c
      do j=1,3
        ylw(1,j)=ylw(1,j)/dcmplx(rrup(lys),0.d0)
        ylw(2,j)=ylw(2,j)/dcmplx(rrup(lys)**2,0.d0)
        ylw(3,j)=ylw(3,j)/dcmplx(rrup(lys),0.d0)
        ylw(4,j)=ylw(4,j)/dcmplx(rrup(lys)**2,0.d0)
c        ylw(5,j)=ylw(5,j)
        ylw(6,j)=ylw(6,j)/dcmplx(rrup(lys),0.d0)
      enddo
c
      do j=1,3
        y0(1,j)=y0(1,j)/dcmplx(rrup(lyr),0.d0)
        y0(2,j)=y0(2,j)/dcmplx(rrup(lyr)**2,0.d0)
        y0(3,j)=y0(3,j)/dcmplx(rrup(lyr),0.d0)
        y0(4,j)=y0(4,j)/dcmplx(rrup(lyr)**2,0.d0)
c        y0(5,j)=y0(5,j)
        y0(6,j)=y0(6,j)/dcmplx(rrup(lyr),0.d0)
      enddo
c
c===============================================================================
c     source function
c===============================================================================
c
      do istp=1,4
        do i=1,6
          b6(i,istp)=(0.d0,0.d0)
        enddo
        b6(istp,istp)=(1.d0,0.d0)
      enddo
      do j=1,3
        do i=1,6
          coef6(i,j)=yup(i,j)
          coef6(i,j+3)=-ylw(i,j)
        enddo
      enddo
      key=0
      call cdsvd500(coef6,b6,6,4,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpsspropg: anormal exit from cdsvd500!'
        return
      endif
      if(lyr.le.lys)then
        do istp=1,4
          do i=1,6
            ypsv(i,istp)=(0.d0,0.d0)
            do j=1,3
              ypsv(i,istp)=ypsv(i,istp)+b6(j,istp)*y0(i,j)
            enddo
          enddo
        enddo
      else
        do istp=1,4
          do i=1,6
            ypsv(i,istp)=(0.d0,0.d0)
            do j=1,3
              ypsv(i,istp)=ypsv(i,istp)+b6(j+3,istp)*y0(i,j)
            enddo
          enddo
        enddo
      endif
c
      return
      end