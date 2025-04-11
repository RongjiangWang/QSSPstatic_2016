      subroutine qpspsvkerng(ldeg,ypsv)
      implicit none
c
c     calculation of response function in frequency-wavelength domain
c     ldeg: harmonic degree
c     ypsv(6,4): psv solution vector (complex) with gravity effect
c
      integer ldeg
      double complex ypsv(6,4)
c
      include 'qpsglobal.h'
c
      integer i,istp,lyup,lylw
c
      do istp=1,4
        do i=1,6
          ypsv(i,istp)=(0.d0,0.d0)
        enddo
      enddo
c
      lyup=lyuppsv(ldeg)
      lylw=lylwpsv(ldeg)
c
      if(lyr.lt.lyup.or.lyr.gt.lylw)return
c
      if(ldeg.eq.0)then
        call qpsspropg0(ypsv,lyup,lylw)
      else
        call qpsspropg(ypsv,ldeg,lyup,lylw)
      endif
c
      return
      end