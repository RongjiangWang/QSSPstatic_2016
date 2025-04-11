      subroutine qpspsvkern(ldeg,ypsv)
      implicit none
c
c     calculation of response function in frequency-wavelength domain
c     ldeg: harmonic degree
c     ypsv(6,4): psv solution vector (complex)
c
      integer ldeg
      double complex ypsv(6,4)
c
      include 'qpsglobal.h'
c
      integer i,ly,istp,lyup,lylw
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
      if(ldeg.eq.1.or.lyr.lt.lyup.or.lyr.gt.lylw)return
c
      if(ldeg.eq.0)then
        do ly=lyup,lylw
          call qpssmat0(ly,lylw)
        enddo
        call qpssprop0(ypsv,lyup,lylw)
      else
        do ly=lyup,min0(lycm-1,lylw)
          call qpssmat(ldeg,ly,lylw)
        enddo
        do ly=max0(lyup,lycm),min0(lycc-1,lylw)
          call qpssmatc(ldeg,ly,lylw)
        enddo
        do ly=max0(lyup,lycc),min0(ly0,lylw)
          call qpssmat(ldeg,ly,lylw)
        enddo
        call qpssprop(ypsv,ldeg,lyup,lylw)
      endif
      return
      end
