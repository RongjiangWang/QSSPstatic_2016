      subroutine qpsvemodel(f,istate)
      implicit none
c
c     calculate q based on the constant q model
c
c     f = frequency
c
      integer istate
      double precision f
c
      include 'qpsglobal.h'
c
      integer ly
      double precision muek
      double complex cs,ctauk,ctaum,ctaumk,cmdk
c
      double complex c1
      data c1/(1.d0,0.d0)/
c
      if(istate.eq.-1)then
        cs=(0.d0,0.d0)
      else
        cs=dcmplx(-2.d0*PI*fi,2.d0*PI*f)
      endif
c
      do ly=1,ly0
        if(eta1up(ly).le.0.d0.and.eta2up(ly).le.0.d0.or.istate.eq.0)then
c
c         Hooke body
c
          crexup(ly)=c1
          crexlw(ly)=c1
        else if(eta1up(ly).gt.0.d0.and.eta2up(ly).le.0.d0.and.
     &          alfaup(ly).gt.0.d0.and.alfaup(ly).lt.1.d0)then
c
c         Standard-Linear-Solid
c
          ctauk=dcmplx(eta1up(ly)/mueup(ly),0.d0)
          cmdk=dcmplx((1.d0-alfaup(ly))/alfaup(ly),0.d0)
          crexup(ly)=(c1+ctauk*cs)/(c1+ctauk*cs+cmdk)
c
          ctauk=dcmplx(eta1lw(ly)/muelw(ly),0.d0)
          cmdk=dcmplx((1.d0-alfalw(ly))/alfalw(ly),0.d0)
          crexlw(ly)=(c1+ctauk*cs)/(c1+ctauk*cs+cmdk)
        else
          if(eta1up(ly).le.0.d0.or.alfaup(ly).ge.1.d0)then
c
c           Maxwell body
c
            ctaum=dcmplx(eta2up(ly)/mueup(ly),0.d0)
            crexup(ly)=ctaum*cs/(c1+ctaum*cs)
c
            ctaum=dcmplx(eta2lw(ly)/muelw(ly),0.d0)
            crexlw(ly)=ctaum*cs/(c1+ctaum*cs)
          else if(eta1up(ly).gt.0.d0.and.eta2up(ly).gt.0.d0.and.
     &            alfaup(ly).gt.0.d0.and.alfaup(ly).lt.1.d0)then
c
c           Burgers body
c
            muek=mueup(ly)*alfaup(ly)/(1.d0-alfaup(ly))
            ctauk=dcmplx(eta1up(ly)/muek,0.d0)
            ctaum=dcmplx(eta2up(ly)/mueup(ly),0.d0)
            ctaumk=dcmplx(eta2up(ly)/muek,0.d0)
            crexup(ly)=(c1+ctauk*cs)*ctaum*cs
     &                /((c1+ctauk*cs)*(c1+ctaum*cs)+ctaumk*cs)
c
            muek=muelw(ly)*alfalw(ly)/(1.d0-alfalw(ly))
            ctauk=dcmplx(eta1lw(ly)/muek,0.d0)
            ctaum=dcmplx(eta2lw(ly)/muelw(ly),0.d0)
            ctaumk=dcmplx(eta2lw(ly)/muek,0.d0)
            crexlw(ly)=(c1+ctauk*cs)*ctaum*cs
     &                /((c1+ctauk*cs)*(c1+ctaum*cs)+ctaumk*cs)
          else
            stop ' Error in qpsvemodel: not identified rheology!'
          endif
c
c         parallel connection with a very small Hooke body
c         for numerical stability
c
          crexup(ly)=dcmplx(1.d0-REXMIN,0.d0)*crexup(ly)
     &              +dcmplx(REXMIN,0.d0)
          crexlw(ly)=dcmplx(1.d0-REXMIN,0.d0)*crexlw(ly)
     &              +dcmplx(REXMIN,0.d0)
        endif
      enddo
      return
      end
