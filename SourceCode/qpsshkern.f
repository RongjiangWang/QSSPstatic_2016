      subroutine qpsshkern(ldeg,ysh)
      implicit none
c
c     calculation of response function in frequency-wavelength domain
c     ldeg: harmonic degree
c     ysh(2,2): sh solution vector (complex)
c
      integer ldeg
      double complex ysh(2,2)
c
      include 'qpsglobal.h'
c
      integer i,ly,istp,lyup,lylw
c
      do istp=1,2
        do i=1,2
          ysh(i,istp)=(0.d0,0.d0)
        enddo
      enddo
c
      lyup=lyupt(ldeg)
      lylw=lylwt(ldeg)
c
      if(ldeg.le.1.or.lyr.lt.lyup.or.lyr.gt.lylw)then
        return
      else if(ldeg.gt.1)then
        do ly=lyup,lylw
          call qpstmat(ldeg,ly,lylw)
        enddo
        call qpstprop(ldeg,ysh,lyup,lylw)
      endif
c
      return
      end

