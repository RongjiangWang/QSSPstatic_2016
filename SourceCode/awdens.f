      subroutine awdens(rho,rup,rlw,grup,grlw,kappa,
     &                  a,b,beta,rhoup,rholw)
      implicit none
      double precision rho,rup,rlw,grup,grlw,kappa,
     &                 a,b,beta,rhoup,rholw
c
      integer n
      double precision rho0,ga,det
      double precision mat(2,2),bat(2)
      double precision sbessj0,sbessj1,sbessy0,sbessy1
c
      integer nmax
      double precision EPS,PI,BIGG
      data nmax/100/
      data EPS,PI,BIGG/1.0d-06,3.14159265358979d0,6.6732d-11/
c
      ga=4.d0*PI*BIGG
      if(rlw.gt.0.d0)then
        n=0
        beta=dsqrt(ga/kappa)*rho
        if(beta.gt.0.5d0*PI/(rup-rlw))then
          beta=0.5d0*PI/(rup-rlw)
          n=nmax
        endif
        rho0=rho
10      n=n+1
        mat(1,1)=sbessj1(beta*rlw)
        mat(1,2)=sbessy1(beta*rlw)
        mat(2,1)=sbessj1(beta*rup)
        mat(2,2)=sbessy1(beta*rup)
        bat(1)=beta*grlw/ga
        bat(2)=beta*grup/ga
        det=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
        a=(bat(1)*mat(2,2)-bat(2)*mat(1,2))/det
        b=(bat(2)*mat(1,1)-bat(1)*mat(2,1))/det
        rhoup=a*sbessj0(beta*rup)+b*sbessy0(beta*rup)
        if(dabs(rhoup-rho0)/rho.gt.EPS.and.n.lt.nmax)then
          rho0=rhoup
          beta=dsqrt(ga/kappa)*rhoup
          goto 10
        endif
        rholw=a*sbessj0(beta*rlw)+b*sbessy0(beta*rlw)
      else
        b=0.d0
        n=0
        beta=dsqrt(ga/kappa)*rho
        if(beta.gt.0.5d0*PI/rup)then
          beta=0.5d0*PI/rup
          n=nmax
        endif
        rho0=rho
20      n=n+1
        a=beta*grup/(ga*sbessj1(beta*rup))
        rhoup=a*sbessj0(beta*rup)
        if(dabs(rhoup-rho0)/rho.gt.EPS.and.n.lt.100)then
          beta=dsqrt(ga/kappa)*rhoup
          goto 20
        endif
        rholw=a*sbessj0(beta*rlw)
      endif
c
      return
      end