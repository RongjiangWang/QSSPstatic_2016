      subroutine qpssprop(ypsv,ldeg,lyup,lylw)
      implicit none
c
c     calculation of spheroidal response
c     ypsv(6,4): solution vector (complex)
c
      integer ldeg,lyup,lylw
      double complex ypsv(6,4)
c
      include 'qpsglobal.h'
c
c     work space
c
      integer i,j,j0,istp,ly,key
      double precision y4max
      double complex cdet,cyswap,ca,cb
      double complex c0(6,3),c1(6,3),cc0(4,2),cc1(4,2)
      double complex y0(6,3),y1(6,3),yc(6,3)
      double complex yup(6,3),ylw(6,3),yupc(4,2),ylwc(4,2)
      double complex wave(6),orth(3,3),orthc(2,2)
      double complex coef6(6,6),b6(6,4),coef4(4,4),b4(4,2)
c
c===============================================================================
c
c     propagation from surface to atmosphere/ocean bottom
c
      if(ldeg.le.1)return
c
      if(lyup.eq.1)then
        do j=1,3
          do i=1,6
            yup(i,j)=(0.d0,0.d0)
          enddo
        enddo
        yup(1,1)=(1.d0,0.d0)
        yup(3,2)=(1.d0,0.d0)
        if(ldeg.eq.1)then
          yup(1,3)=(1.d0,0.d0)
          yup(3,3)=(1.d0,0.d0)
        endif
        yup(5,3)=-dcmplx(grup(1),0.d0)
      else if(lyup.lt.lycm)then
        do j=1,3
          do i=1,6
            yup(i,j)=mas6x6(i,2*j,lyup)
          enddo
        enddo
      else
        stop ' Errot in qpssprop: wrong source depth!'
      endif
      if(lyr.eq.lyup)call cmemcpy(yup,y0,18)
c
c===============================================================================
c
c     propagation from atmosphere/ocean bottom to source or core-mantle boundary
c
      do ly=lyup,lys-1
        wave(1)=dcmplx((rrlw(ly)/rrup(ly))**ldeg,0.d0)
        wave(2)=dcmplx((rrlw(ly)/rrup(ly))**(ldeg+1),0.d0)
        wave(3)=dcmplx((rrlw(ly)/rrup(ly))**(ldeg+2),0.d0)
        wave(4)=dcmplx((rrlw(ly)/rrup(ly))**(ldeg-1),0.d0)
        wave(5)=dcmplx((rrlw(ly)/rrup(ly))**ldeg,0.d0)
        wave(6)=dcmplx((rrlw(ly)/rrup(ly))**(ldeg+1),0.d0)
c
        call caxcb(mas6x6inv(1,1,ly),yup,6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(2,1)*c0(4,2)*c0(6,3)
     &      +c0(4,1)*c0(6,2)*c0(2,3)
     &      +c0(6,1)*c0(2,2)*c0(4,3)
     &      -c0(6,1)*c0(4,2)*c0(2,3)
     &      -c0(4,1)*c0(2,2)*c0(6,3)
     &      -c0(2,1)*c0(6,2)*c0(4,3)
        orth(1,1)=(c0(4,2)*c0(6,3)-c0(4,3)*c0(6,2))/cdet
        orth(2,1)=(c0(4,3)*c0(6,1)-c0(4,1)*c0(6,3))/cdet
        orth(3,1)=(c0(4,1)*c0(6,2)-c0(4,2)*c0(6,1))/cdet
        orth(1,2)=(c0(2,3)*c0(6,2)-c0(2,2)*c0(6,3))/cdet
        orth(2,2)=(c0(2,1)*c0(6,3)-c0(2,3)*c0(6,1))/cdet
        orth(3,2)=(c0(2,2)*c0(6,1)-c0(2,1)*c0(6,2))/cdet
        orth(1,3)=(c0(2,2)*c0(4,3)-c0(2,3)*c0(4,2))/cdet
        orth(2,3)=(c0(2,3)*c0(4,1)-c0(2,1)*c0(4,3))/cdet
        orth(3,3)=(c0(2,1)*c0(4,2)-c0(2,2)*c0(4,1))/cdet
c
        call caxcb(c0,orth,6,3,3,c1)
        if(ly.ge.lyr)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orth,6,3,3,y1)
          call cmemcpy(y1,y0,18)
          do j=1,3
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j)
            enddo
          enddo
        endif
c
        c1(1,1)=c1(1,1)*wave(1)*wave(2)
        c1(2,1)=(1.d0,0.d0)
        c1(3,1)=c1(3,1)*wave(3)*wave(2)
        c1(4,1)=(0.d0,0.d0)
        c1(5,1)=c1(5,1)*wave(5)*wave(2)
        c1(6,1)=(0.d0,0.d0)
c
        c1(1,2)=c1(1,2)*wave(1)*wave(4)
        c1(2,2)=(0.d0,0.d0)
        c1(3,2)=c1(3,2)*wave(3)*wave(4)
        c1(4,2)=(1.d0,0.d0)
        c1(5,2)=c1(5,2)*wave(5)*wave(4)
        c1(6,2)=(0.d0,0.d0)
c
        c1(1,3)=c1(1,3)*wave(1)*wave(6)
        c1(2,3)=(0.d0,0.d0)
        c1(3,3)=c1(3,3)*wave(3)*wave(6)
        c1(4,3)=(0.d0,0.d0)
        c1(5,3)=c1(5,3)*wave(5)*wave(6)
        c1(6,3)=(1.d0,0.d0)
c
        call caxcb(mas6x6(1,1,ly),c1,6,6,3,yup)
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
        do j=1,3
          do i=1,6
            ylw(i,j)=mas6x6(i,2*j-1,lylw)
          enddo
        enddo
        if(lylw.eq.lyr.and.lylw.gt.lys)call cmemcpy(ylw,y0,18)
      endif
c
      do ly=lylw-1,lycc,-1
        wave(1)=dcmplx((rrlw(ly)/rrup(ly))**ldeg,0.d0)
        wave(2)=dcmplx((rrlw(ly)/rrup(ly))**(ldeg+1),0.d0)
        wave(3)=dcmplx((rrlw(ly)/rrup(ly))**(ldeg+2),0.d0)
        wave(4)=dcmplx((rrlw(ly)/rrup(ly))**(ldeg-1),0.d0)
        wave(5)=dcmplx((rrlw(ly)/rrup(ly))**ldeg,0.d0)
        wave(6)=dcmplx((rrlw(ly)/rrup(ly))**(ldeg+1),0.d0)
c
        call caxcb(mas6x6inv(1,1,ly),ylw,6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(1,1)*c0(3,2)*c0(5,3)
     &      +c0(3,1)*c0(5,2)*c0(1,3)
     &      +c0(5,1)*c0(1,2)*c0(3,3)
     &      -c0(5,1)*c0(3,2)*c0(1,3)
     &      -c0(3,1)*c0(1,2)*c0(5,3)
     &      -c0(1,1)*c0(5,2)*c0(3,3)
        orth(1,1)=(c0(3,2)*c0(5,3)-c0(3,3)*c0(5,2))/cdet
        orth(2,1)=(c0(3,3)*c0(5,1)-c0(3,1)*c0(5,3))/cdet
        orth(3,1)=(c0(3,1)*c0(5,2)-c0(3,2)*c0(5,1))/cdet
        orth(1,2)=(c0(1,3)*c0(5,2)-c0(1,2)*c0(5,3))/cdet
        orth(2,2)=(c0(1,1)*c0(5,3)-c0(1,3)*c0(5,1))/cdet
        orth(3,2)=(c0(1,2)*c0(5,1)-c0(1,1)*c0(5,2))/cdet
        orth(1,3)=(c0(1,2)*c0(3,3)-c0(1,3)*c0(3,2))/cdet
        orth(2,3)=(c0(1,3)*c0(3,1)-c0(1,1)*c0(3,3))/cdet
        orth(3,3)=(c0(1,1)*c0(3,2)-c0(1,2)*c0(3,1))/cdet
c
        call caxcb(c0,orth,6,3,3,c1)
        if(ly.lt.lyr)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orth,6,3,3,y1)
          call cmemcpy(y1,y0,18)
          do j=1,3
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
        c1(1,1)=(1.d0,0.d0)
        c1(2,1)=c1(2,1)*wave(2)*wave(1)
        c1(3,1)=(0.d0,0.d0)
        c1(4,1)=c1(4,1)*wave(4)*wave(1)
        c1(5,1)=(0.d0,0.d0)
        c1(6,1)=c1(6,1)*wave(6)*wave(1)
c
        c1(1,2)=(0.d0,0.d0)
        c1(2,2)=c1(2,2)*wave(2)*wave(3)
        c1(3,2)=(1.d0,0.d0)
        c1(4,2)=c1(4,2)*wave(4)*wave(3)
        c1(5,2)=(0.d0,0.d0)
        c1(6,2)=c1(6,2)*wave(6)*wave(3)
c
        c1(1,3)=(0.d0,0.d0)
        c1(2,3)=c1(2,3)*wave(2)*wave(5)
        c1(3,3)=(0.d0,0.d0)
        c1(4,3)=c1(4,3)*wave(4)*wave(5)
        c1(5,3)=(1.d0,0.d0)
        c1(6,3)=c1(6,3)*wave(6)*wave(5)
c
        call caxcb(mas6x6(1,1,ly),c1,6,6,3,ylw)
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
        y4max=cdabs(ylw(4,3))
        j0=3
        do j=1,2
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
          ylwc(1,j)=ylw(1,j)-ylw(4,j)*ylw(1,3)/ylw(4,3)
          ylwc(2,j)=ylw(2,j)-ylw(4,j)*ylw(2,3)/ylw(4,3)
          ylwc(3,j)=ylw(5,j)-ylw(4,j)*ylw(5,3)/ylw(4,3)
          ylwc(4,j)=ylw(6,j)-ylw(4,j)*ylw(6,3)/ylw(4,3)
        enddo
        if(lycc.le.lyr)then
          do i=1,6
            cyswap=y0(i,j0)
            y0(i,j0)=y0(i,3)
            y0(i,3)=cyswap
          enddo
          do j=1,2
            do i=1,6
              y0(i,j)=y0(i,j)-ylw(4,j)*y0(i,3)/ylw(4,3)
            enddo
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
        endif
        ca=ylwc(2,2)
        cb=ylwc(2,1)
        do i=1,4
          ylwc(i,1)=ca*ylwc(i,1)-cb*ylwc(i,2)
          ylwc(i,2)=(0.d0,0.d0)
        enddo
        ylwc(2,1)=(0.d0,0.d0)
        ylwc(2,2)=(1.d0,0.d0)
        if(lycc.le.lyr.and.lycc.gt.lys)then
          do i=1,6
            y0(i,1)=ca*y0(i,j)-cb*y0(i,2)
            y0(i,2)=(0.d0,0.d0)
          enddo
        endif
      else if(lylw.ge.lycm)then
        do j=1,2
          do i=1,4
            ylwc(i,j)=mas4x4(i,2*j-1,lylw)
          enddo
        enddo
        if(lylw.eq.lyr.and.lylw.gt.lys)then
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=(0.d0,0.d0)
            y0(3,j)=ylwc(2,j)
            y0(4,j)=(0.d0,0.d0)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
        endif
      endif
c
      do ly=min0(lylw,lycc)-1,lycm,-1
        wave(1)=dcmplx((rrlw(ly)/rrup(ly))**ldeg,0.d0)
        wave(2)=dcmplx((rrlw(ly)/rrup(ly))**(ldeg+1),0.d0)
        wave(3)=dcmplx((rrlw(ly)/rrup(ly))**ldeg,0.d0)
        wave(4)=dcmplx((rrlw(ly)/rrup(ly))**(ldeg+1),0.d0)
c
        call caxcb(mas4x4inv(1,1,ly),ylwc,4,4,2,cc0)
c
c       orthonormalization of the p-sv modes
c
        cdet=cc0(1,1)*cc0(3,2)-cc0(3,1)*cc0(1,2)
        orthc(1,1)=cc0(3,2)/cdet
        orthc(1,2)=-cc0(1,2)/cdet
        orthc(2,1)=-cc0(3,1)/cdet
        orthc(2,2)=cc0(1,1)/cdet
c
        call caxcb(cc0,orthc,4,2,2,cc1)
        if(ly.lt.lyr)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orthc,6,2,2,y1)
          call cmemcpy(y1,y0,12)
          do j=1,2
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
c
        cc1(1,1)=(1.d0,0.d0)
        cc1(2,1)=cc1(2,1)*wave(2)*wave(1)
        cc1(3,1)=(0.d0,0.d0)
        cc1(4,1)=cc1(4,1)*wave(4)*wave(1)
c
        cc1(1,2)=(0.d0,0.d0)
        cc1(2,2)=cc1(2,2)*wave(2)*wave(3)
        cc1(3,2)=(1.d0,0.d0)
        cc1(4,2)=cc1(4,2)*wave(4)*wave(3)
c
        call caxcb(mas4x4(1,1,ly),cc1,4,4,2,ylwc)
        if(ly.eq.lyr.and.ly.gt.lys)then
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=(0.d0,0.d0)
            y0(3,j)=ylwc(2,j)
            y0(4,j)=(0.d0,0.d0)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
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
        do j=1,2
          ylw(1,j)=ylwc(1,j)
          ylw(2,j)=(0.d0,0.d0)
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
        do j=1,3
          do i=1,6
            ylw(i,j)=mas6x6(i,2*j-1,lylw)
          enddo
        enddo
        if(lylw.eq.lyr.and.lylw.gt.lys)call cmemcpy(ylw,y0,18)
      endif
c
      do ly=min0(lylw,lycm)-1,lys,-1
        wave(1)=dcmplx((rrlw(ly)/rrup(ly))**ldeg,0.d0)
        wave(2)=dcmplx((rrlw(ly)/rrup(ly))**(ldeg+1),0.d0)
        wave(3)=dcmplx((rrlw(ly)/rrup(ly))**(ldeg+2),0.d0)
        wave(4)=dcmplx((rrlw(ly)/rrup(ly))**(ldeg-1),0.d0)
        wave(5)=dcmplx((rrlw(ly)/rrup(ly))**ldeg,0.d0)
        wave(6)=dcmplx((rrlw(ly)/rrup(ly))**(ldeg+1),0.d0)
c
        call caxcb(mas6x6inv(1,1,ly),ylw,6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(1,1)*c0(3,2)*c0(5,3)
     &      +c0(3,1)*c0(5,2)*c0(1,3)
     &      +c0(5,1)*c0(1,2)*c0(3,3)
     &      -c0(5,1)*c0(3,2)*c0(1,3)
     &      -c0(3,1)*c0(1,2)*c0(5,3)
     &      -c0(1,1)*c0(5,2)*c0(3,3)
        orth(1,1)=(c0(3,2)*c0(5,3)-c0(3,3)*c0(5,2))/cdet
        orth(2,1)=(c0(3,3)*c0(5,1)-c0(3,1)*c0(5,3))/cdet
        orth(3,1)=(c0(3,1)*c0(5,2)-c0(3,2)*c0(5,1))/cdet
        orth(1,2)=(c0(1,3)*c0(5,2)-c0(1,2)*c0(5,3))/cdet
        orth(2,2)=(c0(1,1)*c0(5,3)-c0(1,3)*c0(5,1))/cdet
        orth(3,2)=(c0(1,2)*c0(5,1)-c0(1,1)*c0(5,2))/cdet
        orth(1,3)=(c0(1,2)*c0(3,3)-c0(1,3)*c0(3,2))/cdet
        orth(2,3)=(c0(1,3)*c0(3,1)-c0(1,1)*c0(3,3))/cdet
        orth(3,3)=(c0(1,1)*c0(3,2)-c0(1,2)*c0(3,1))/cdet
c
        call caxcb(c0,orth,6,3,3,c1)
        if(ly.lt.lyr)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orth,6,3,3,y1)
          call cmemcpy(y1,y0,18)
          do j=1,3
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
        c1(1,1)=(1.d0,0.d0)
        c1(2,1)=c1(2,1)*wave(2)*wave(1)
        c1(3,1)=(0.d0,0.d0)
        c1(4,1)=c1(4,1)*wave(4)*wave(1)
        c1(5,1)=(0.d0,0.d0)
        c1(6,1)=c1(6,1)*wave(6)*wave(1)
c
        c1(1,2)=(0.d0,0.d0)
        c1(2,2)=c1(2,2)*wave(2)*wave(3)
        c1(3,2)=(1.d0,0.d0)
        c1(4,2)=c1(4,2)*wave(4)*wave(3)
        c1(5,2)=(0.d0,0.d0)
        c1(6,2)=c1(6,2)*wave(6)*wave(3)
c
        c1(1,3)=(0.d0,0.d0)
        c1(2,3)=c1(2,3)*wave(2)*wave(5)
        c1(3,3)=(0.d0,0.d0)
        c1(4,3)=c1(4,3)*wave(4)*wave(5)
        c1(5,3)=(1.d0,0.d0)
        c1(6,3)=c1(6,3)*wave(6)*wave(5)
c
        call caxcb(mas6x6(1,1,ly),c1,6,6,3,ylw)
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
        print *,' Warning in qpssprop: anormal exit from cdsvd500!'
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
      return
      end