      subroutine qpstmat(ldeg,ly,lylw)
      implicit none
c
      integer ldeg,ly,lylw
c
      include 'qpsglobal.h'
c
      double complex clm1,clp2,c2lp1
      double complex cmue
c
      double complex c1
      data c1/(1.d0,0.d0)/
c
      clm1=dcmplx(dble(ldeg-1),0.d0)
      clp2=dcmplx(dble(ldeg+2),0.d0)
      c2lp1=dcmplx(dble(2*ldeg+1),0.d0)
      cmue=dcmplx(0.25d0*(mueup(ly)+muelw(ly)),0.d0)
     &    *(crexup(ly)+crexlw(ly))
c
      mat2x2(1,1,ly)=c1
      mat2x2(2,1,ly)=cmue*clm1
      mat2x2(1,2,ly)=-c1
      mat2x2(2,2,ly)=cmue*clp2
c
      if(ly.eq.lylw)return
c
c     calculate inverse matrix at upper radius
c
      mat2x2inv(1,1,ly)=clp2/c2lp1
      mat2x2inv(1,2,ly)=c1/(cmue*c2lp1)
      mat2x2inv(2,1,ly)=-clm1/c2lp1
      mat2x2inv(2,2,ly)=c1/(cmue*c2lp1)
c
      return
      end