      subroutine qpssmatc(ldeg,ly,lylw)
      implicit none
c
c     calculate 4x4 spheroidal layer matrix for a liquid shell
c
      integer ldeg,ly,lylw
c
      include 'qpsglobal.h'
c
      double complex cga,cldeg,clp1,c2lp1,c3lp2
c
      double complex c0,c1,c2,c3
      data c0,c1,c2,c3/(0.d0,0.d0),(1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0)/
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      clp1=cldeg+c1
      c2lp1=c2*cldeg+c1
      c3lp2=c3*cldeg+c2
      cga=dcmplx(2.d0*PI*BIGG*(rhoup(ly)+rholw(ly)),0.d0)
c
c     y1 <- y1 (normal displacement)
c     y2 <- y3 (horizontal displacement)
c     y3 <- y5 (potential)
c     y4 <- y6 (gravity)
c
      mas4x4(1,1,ly)=cldeg
      mas4x4(2,1,ly)=c1
      mas4x4(3,1,ly)=c0
      mas4x4(4,1,ly)=-cga*cldeg
c
      mas4x4(1,2,ly)=-clp1
      mas4x4(2,2,ly)=c1
      mas4x4(3,2,ly)=cga
      mas4x4(4,2,ly)=c0
c
      mas4x4(1,3,ly)=c0
      mas4x4(2,3,ly)=c0
      mas4x4(3,3,ly)=c1
      mas4x4(4,3,ly)=c2lp1
c
      mas4x4(1,4,ly)=c0
      mas4x4(2,4,ly)=c0
      mas4x4(3,4,ly)=c0
      mas4x4(4,4,ly)=c1
c
      if(ly.eq.lylw)return
c
      mas4x4inv(1,1,ly)=c1/c2lp1
      mas4x4inv(1,2,ly)=clp1/c2lp1
      mas4x4inv(1,3,ly)=c0
      mas4x4inv(1,4,ly)=c0
c
      mas4x4inv(2,1,ly)=-c1/c2lp1
      mas4x4inv(2,2,ly)=cldeg/c2lp1
      mas4x4inv(2,3,ly)=c0
      mas4x4inv(2,4,ly)=c0
c
      mas4x4inv(3,1,ly)=cga/c2lp1
      mas4x4inv(3,2,ly)=-cga*cldeg/c2lp1
      mas4x4inv(3,3,ly)=c1
      mas4x4inv(3,4,ly)=c0
c
      mas4x4inv(4,1,ly)=-cga*clp1/c2lp1
      mas4x4inv(4,2,ly)=cga*cldeg*c3lp2/c2lp1
      mas4x4inv(4,3,ly)=-c2lp1
      mas4x4inv(4,4,ly)=c1
c
      return
      end