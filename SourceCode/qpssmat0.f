      subroutine qpssmat0(ly,lylw)
      implicit none
c
c     calculate 3x3 spheroidal layer matrix for a solid shell in case of degree l = 0
c
      integer ly,lylw
c
      include 'qpsglobal.h'
c
      integer i,j,key
      double complex c3ksi,c4mue,cga
      double complex mas(3,3)
c
      double complex c0,c1,c2,c3,c4
      data c0,c1,c2,c3,c4/(0.d0,0.d0),(1.d0,0.d0),
     &                    (2.d0,0.d0),(3.d0,0.d0),(4.d0,0.d0)/
c
      c4mue=dcmplx(mueup(ly)+muelw(ly),0.d0)*(crexup(ly)+crexlw(ly))
      c3ksi=dcmplx(1.5d0*(kapup(ly)+kaplw(ly)),0.d0)
      cga=dcmplx(2.d0*PI*BIGG*(rhoup(ly)+rholw(ly)),0.d0)
c
      mas3x3(1,1,ly)=c1
      mas3x3(2,1,ly)=c3ksi
      mas3x3(3,1,ly)=cga/c2
c
      mas3x3(1,2,ly)=c1
      mas3x3(2,2,ly)=-c4mue
      mas3x3(3,2,ly)=-cga
c
      mas3x3(1,3,ly)=c0
      mas3x3(2,3,ly)=c0
      mas3x3(3,3,ly)=c1
c
      if(ly.eq.lylw)return
c
c     calculate inverse matrix
c
      do j=1,3
        do i=1,3
          mas(i,j)=mas3x3(i,j,ly)
          mas3x3inv(i,j,ly)=(0.d0,0.d0)
        enddo
        mas3x3inv(j,j,ly)=c1
      enddo
      key=0
      call cdsvd500(mas,mas3x3inv(1,1,ly),3,3,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpssmat0: anormal exit from cdsvd500!'
        return
      endif
c
      return
      end