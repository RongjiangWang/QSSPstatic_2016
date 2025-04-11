      subroutine qpsdifmat0(ly,ldeg,rr,mat)
      implicit none
      integer ly,ldeg
      double precision rr
      double complex mat(6,6)
c
      include 'qpsglobal.h'
c
c     3x3 coefficient matrix for spheroidal mode l = 0
c
      double precision up,lw,xx,x1,r1
      double precision grrr,rhorr,kaprr,muerr,mass,drho,rho1
      double complex crr,crhorr,clamrr,cmuerr,cksirr,cgrrr,cgarr
      double precision sbessj0,sbessj1,sbessy0,sbessy1
c
      double complex c1,c2,c3,c4
      data c1,c2,c3,c4/(1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0),(4.d0,0.d0)/
c
      up=(rr-rrlw(ly))/(rrup(ly)-rrlw(ly))
      lw=1.d0-up
      r1=rrlw(ly)
c
      if(beta(ly).gt.0.d0)then
        xx=beta(ly)*rr
        x1=beta(ly)*r1
        mass=4.d0*PI/beta(ly)
     &      *(rhoa(ly)*(rr**2*sbessj1(xx)
     &                 -r1**2*sbessj1(x1))
     &       +rhob(ly)*(rr**2*sbessy1(xx)
     &                 -r1**2*sbessy1(x1)))
        grrr=(grlw(ly)*r1**2+BIGG*mass)/rr**2
        rhorr=rhoa(ly)*sbessj0(xx)+rhob(ly)*sbessy0(xx)
        kaprr=kapup(ly)*(rhorr/rhoup(ly))**2
      else
        drho=(rhoup(ly)-rholw(ly))/(rrup(ly)-rrlw(ly))
        rho1=rholw(ly)-drho*rrlw(ly) 
        mass=PI*(rr-r1)*((4.d0/3.d0)*rho1*(rr**2+rr*r1+r1**2)
     &      +drho*(rr**3+rr**2*r1+rr*r1**2+r1**3))
        grrr=(grlw(ly)*r1**2+BIGG*mass)/rr**2
        rhorr=up*rhoup(ly)+lw*rholw(ly)
        kaprr=up*kapup(ly)+lw*kaplw(ly)
      endif
      muerr=up*mueup(ly)+lw*muelw(ly)
c
      crr=dcmplx(rr,0.d0)
      crhorr=dcmplx(rhorr,0.d0)
      cmuerr=dcmplx(muerr,0.d0)
     &      *(dcmplx(up,0.d0)*crexup(ly)+dcmplx(lw,0.d0)*crexlw(ly))
      clamrr=dcmplx(kaprr,0.d0)-cmuerr*c2/c3
      cksirr=clamrr+c2*cmuerr
      cgrrr=dcmplx(grrr,0.d0)
      cgarr=dcmplx(4.d0*PI*BIGG*rhorr,0.d0)
c
      mat(1,1)=(c1-c2*clamrr/cksirr)/crr
      mat(1,2)=c1/cksirr/crr
      mat(1,3)=(0.d0,0.d0)
c
      mat(2,1)=c4*(cmuerr*(c1+c2*clamrr/cksirr)/crr-crhorr*cgrrr)
      mat(2,2)=c2*clamrr/cksirr/crr
      mat(2,3)=(0.d0,0.d0)
c
      mat(3,1)=cgarr/crr
      mat(3,2)=(0.d0,0.d0)
      mat(3,3)=(0.d0,0.d0)
c
      return
      end