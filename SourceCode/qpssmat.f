      subroutine qpssmat(ldeg,ly,lylw)
      implicit none
c
c     calculate 6x6 spheroidal layer matrix for a solid shell
c
      integer ldeg,ly,lylw
c
      include 'qpsglobal.h'
c
      integer i,j,key
      double complex cldeg,clp1,clp2,clm1,cllm1,cllp1,c2lp1,cllp2
      double complex ca,c2mue,cksi,ceta,cga
      double complex mas(6,6)
c
      double complex c0,c1,c2,c3,c4
      data c0,c1,c2,c3,c4/(0.d0,0.d0),(1.d0,0.d0),(2.d0,0.d0),
     &                    (3.d0,0.d0),(4.d0,0.d0)/
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      c2mue=dcmplx(0.5d0*(mueup(ly)+muelw(ly)),0.d0)
     &     *(crexup(ly)+crexlw(ly))
      cga=dcmplx(2.d0*PI*BIGG*(rhoup(ly)+rholw(ly)),0.d0)
c
      cksi=c2mue/(dcmplx(kapup(ly)+kaplw(ly),0.d0)+c2mue*c4/c3)
      ceta=c1-cksi
c
      clp1=cldeg+c1
      clp2=cldeg+c2
      clm1=cldeg-c1
      cllm1=cldeg*clm1
      cllp1=cldeg*clp1
      cllp2=cldeg*clp2
      c2lp1=c2*cldeg+c1
c
      mas6x6(1,1,ly)=cldeg
      mas6x6(2,1,ly)=c2mue*cllm1
      mas6x6(3,1,ly)=c1
      mas6x6(4,1,ly)=c2mue*clm1
      mas6x6(5,1,ly)=cga
      mas6x6(6,1,ly)=cga*clp1
c
      mas6x6(1,2,ly)=-clp1
      mas6x6(2,2,ly)=c2mue*clp1*clp2
      mas6x6(3,2,ly)=c1
      mas6x6(4,2,ly)=-c2mue*clp2
      mas6x6(5,2,ly)=cga
      mas6x6(6,2,ly)=cga*clp1
c
      mas6x6(1,3,ly)=cldeg*ceta-c2*cksi
      mas6x6(2,3,ly)=c2mue*(cllm1*ceta+c4*cksi-c3)
      mas6x6(3,3,ly)=ceta+c2/clp1
      mas6x6(4,3,ly)=c2mue*(clm1*ceta-c2*cksi+c2lp1/clp1)
      mas6x6(5,3,ly)=-cga*cksi
      mas6x6(6,3,ly)=cga*(clp1*ceta-c2lp1)
c
      mas6x6(1,4,ly)=cllp1*ceta+c2*cldeg*cksi
      mas6x6(2,4,ly)=c2mue*cldeg*(-clp1*clp2*ceta-c4*cksi+c3)
      mas6x6(3,4,ly)=-cldeg*ceta+c2
      mas6x6(4,4,ly)=c2mue*(cllp2*ceta+c2*cldeg*cksi-c2lp1)
      mas6x6(5,4,ly)=cga*cldeg*cksi
      mas6x6(6,4,ly)=-cga*cllp1*ceta
c
      mas6x6(1,5,ly)=c0
      mas6x6(2,5,ly)=c0
      mas6x6(3,5,ly)=c0
      mas6x6(4,5,ly)=c0
      mas6x6(5,5,ly)=c1
      mas6x6(6,5,ly)=c2lp1
c
      mas6x6(1,6,ly)=c0
      mas6x6(2,6,ly)=c0
      mas6x6(3,6,ly)=c0
      mas6x6(4,6,ly)=c0
      mas6x6(5,6,ly)=c1
      mas6x6(6,6,ly)=c0
c
      if(ly.eq.lylw)return
c
      do j=1,6
        do i=1,6
          mas(i,j)=mas6x6(i,j,ly)
          mas6x6inv(i,j,ly)=(0.d0,0.d0)
        enddo
        mas6x6inv(j,j,ly)=c1
      enddo
      key=0
      call cdsvd500(mas,mas6x6inv(1,1,ly),6,6,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpssmat: anormal exit from cdsvd500!'
        return
      endif
c
      return
      end