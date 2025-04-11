      subroutine qpsspropg0(ypsv,lyup,lylw)
      implicit none
c
c     calculation of spheroidal response (with gravity) for degree 0
c     ypsv(6,4): solution vector (complex)
c
      integer lyup,lylw
      double complex ypsv(6,4)
c
      include 'qpsglobal.h'
c
c     work space
c
      integer i,istp,ly,ily,nly,key
      double precision f,rr1,rr2,dlnr,h
      double complex y0(3),yup(3),ylw(3),coef(2,2),b(2,2)
      external qpsdifmat0
c
c
c===============================================================================
c
c     propagation from surface to source
c
      if(lyup.eq.1)then
        yup(1)=(1.d0,0.d0)
        yup(2)=(0.d0,0.d0)
        yup(3)=(0.d0,0.d0)
      else
        call qpssmat0(lyup,lyup)
        do i=1,3
          yup(i)=mas3x3(i,2,lyup)
        enddo
      endif
c
      if(lyr.eq.lyup)call cmemcpy(yup,y0,3)
c
      do ly=lyup,lys-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(0.5d0*h/rrlw(ly))
        dlnr=dlog(rrlw(ly)/rrup(ly))/dble(nly)
        rr2=rrup(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrup(ly)*dexp(dble(ily)*dlnr)
          call ruku(yup,3,1,ly,0,qpsdifmat0,rr1,rr2,nruku(0,ly))
        enddo
        if(ly.eq.lyr-1)call cmemcpy(yup,y0,3)
      enddo
      yup(1)=yup(1)/dcmplx(rrup(lys),0.d0)
      yup(2)=yup(2)/dcmplx(rrup(lys)**2,0.d0)
c     yup(3)=yup(3)
c
c===============================================================================
c
c     propagation from bottom to source
c
      call qpssmat0(lylw,lylw)
      do i=1,3
        ylw(i)=mas3x3(i,1,lylw)
      enddo
      if(lylw.eq.lyr.and.lylw.gt.lys)call cmemcpy(ylw,y0,3)
c
      do ly=lylw-1,lys,-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(0.5d0*h/rrlw(ly))
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
          call ruku(ylw,3,1,ly,0,qpsdifmat0,rr1,rr2,nruku(0,ly))
        enddo
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,3)
      enddo
      ylw(1)=ylw(1)/dcmplx(rrup(lys),0.d0)
      ylw(2)=ylw(2)/dcmplx(rrup(lys)**2,0.d0)
c     ylw(3)=ylw(3)
c
      y0(1)=y0(1)/dcmplx(rrup(lyr),0.d0)
      y0(2)=y0(2)/dcmplx(rrup(lyr)**2,0.d0)
c     y0(3)=y0(3)
c
c===============================================================================
c     source function
c===============================================================================
c
      b(1,1)=(1.d0,0.d0)
      b(2,1)=(0.d0,0.d0)
      b(1,2)=(0.d0,0.d0)
      b(2,2)=(1.d0,0.d0)
      do i=1,2
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
      enddo
      key=0
      call cdsvd500(coef,b,2,2,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpsspropg0: anormal exit from cdsvd500!'
        return
      endif
      if(lyr.le.lys)then
        do istp=1,2
          do i=1,2
            ypsv(i,istp)=b(1,istp)*y0(i)
          enddo
          ypsv(5,istp)=b(1,istp)*y0(3)
          ypsv(6,istp)=ypsv(5,istp)/dcmplx(rrup(lyr),0.d0)
        enddo
      else
        do istp=1,2
          do i=1,2
            ypsv(i,istp)=b(2,istp)*y0(i)
          enddo
          ypsv(5,istp)=b(2,istp)*(y0(3)-ylw(3))+b(1,istp)*yup(3)
          ypsv(6,istp)=ypsv(5,istp)/dcmplx(rrup(lyr),0.d0)
        enddo
      endif
      return
      end