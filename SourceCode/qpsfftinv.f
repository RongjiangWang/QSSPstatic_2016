      subroutine qpsfftinv(ierr)
      implicit none
      integer ierr
c
      include 'qpsglobal.h'
c
      integer i,lf,mf,ir,it,isnap,it1,it2
      double precision f,t,a,b,sr,si,tau,ymax
      double precision ura,urb,uta,utb,uza,uzb,ppa,ppb
      double precision gra,grb,gda,gdb,tra,trb,tta,ttb
      double precision y(nrmax),antia(2*nfmax)
      double complex lpf
      double complex cs(nfmax)
c
      double complex cswap(2*nfmax)
      double precision dswap(4*nfmax)
      equivalence (cswap,dswap)
c
      do lf=1,nf
        cs(lf)=dcmplx(-2.d0*PI*fi,2.d0*PI*dble(lf-1)*df)
      enddo
c
      do it=1,nt
        antia(it)=dexp(-2.d0*PI*fi*dble(it-1)*dt)*df*dt
      enddo
c
      if(seldis)then
c
c       calculate and output north or radial displacement
c
        do ir=1,nr
          ymax=0.d0
          do lf=1,nf
            cswap(lf)=ur(lf,ir)-ur(0,ir)
            ymax=dmax1(ymax,cdabs(cswap(lf)))
          enddo
          if(dt.gt.0.5d0*taumin.and.cdabs(cswap(nf)).gt.0.1d0*ymax)then
            tau=dt*dmin1(2.d0,dsqrt(cdabs(cswap(nf))/ymax))
            do lf=1,nf
              lpf=(1.d0,0.d0)/((1.d0,0.d0)+cs(lf)*dcmplx(tau,0.d0))**2
              cswap(lf)=cswap(lf)*lpf
            enddo
          endif
          mf=1
          do lf=2*nf,nf+2,-1
            mf=mf+1
            cswap(lf)=dconjg(cswap(mf))
          enddo
          cswap(nf+1)=(0.d0,0.d0)
c
c         convention for Fourier transform:
c         f(t)=\int F(f) exp(i2\pi f t) df
c
          call four1(dswap,2*nf,+1)
c
          sr=dreal(ur(0,ir))
          si=sr+dswap(1)*antia(1)
          ur(1,ir)=dcmplx(sr,si)
          i=1
          t=0.d0
          do lf=2,nf
            i=i+2
            t=t+dt
            sr=si+dswap(i)*antia((i+1)/2)
            i=i+2
            t=t+dt
            si=sr+dswap(i)*antia((i+1)/2)
            ur(lf,ir)=dcmplx(sr,si)
          enddo
        enddo
c
        open(20,file=urout,status='unknown')
        write(20,'(a,$)')' Time[day]  '
        do ir=1,nr-1
          write(20,'(a16,$)')rname(ir)
        enddo
        write(20,'(a16)')rname(nr)
        do it=1,ntcut
          t=dble(it-1)*dt/DAY2SEC
          lf=(it+1)/2
          if(it.gt.2*(it/2))then
            do ir=1,nr
              y(ir)=dreal(ur(lf,ir))
            enddo
          else
            do ir=1,nr
              y(ir)=dimag(ur(lf,ir))
            enddo
          endif
          write(20,'(f12.3,$)')t
          do ir=1,nr-1
            write(20,'(E16.8,$)')y(ir)
          enddo
          write(20,'(E16.8)')y(nr)
        enddo
        close(20)
c
c       calculate and output east or transversal displacement
c
        do ir=1,nr
          ymax=0.d0
          do lf=1,nf
            cswap(lf)=ut(lf,ir)-ut(0,ir)
            ymax=dmax1(ymax,cdabs(cswap(lf)))
          enddo
          if(dt.gt.0.5d0*taumin.and.cdabs(cswap(nf)).gt.0.1d0*ymax)then
            tau=dt*dmin1(2.d0,dsqrt(cdabs(cswap(nf))/ymax))
            do lf=1,nf
              lpf=(1.d0,0.d0)/((1.d0,0.d0)+cs(lf)*dcmplx(tau,0.d0))**2
              cswap(lf)=cswap(lf)*lpf
            enddo
          endif
          mf=1
          do lf=2*nf,nf+2,-1
            mf=mf+1
            cswap(lf)=dconjg(cswap(mf))
          enddo
          cswap(nf+1)=(0.d0,0.d0)
c
c         convention for Fourier transform:
c         f(t)=\int F(f) exp(i2\pi f t) df
c
          call four1(dswap,2*nf,+1)
c
          sr=dreal(ut(0,ir))
          si=sr+dswap(1)*antia(1)
          ut(1,ir)=dcmplx(sr,si)
          i=1
          t=0.d0
          do lf=2,nf
            i=i+2
            t=t+dt
            sr=si+dswap(i)*antia((i+1)/2)
            i=i+2
            t=t+dt
            si=sr+dswap(i)*antia((i+1)/2)
            ut(lf,ir)=dcmplx(sr,si)
          enddo
        enddo
c
        open(20,file=utout,status='unknown')
        write(20,'(a,$)')' Time[day]  '
        do ir=1,nr-1
          write(20,'(a16,$)')rname(ir)
        enddo
        write(20,'(a16)')rname(nr)
        do it=1,ntcut
          t=dble(it-1)*dt/DAY2SEC
          lf=(it+1)/2
          if(it.gt.2*(it/2))then
            do ir=1,nr
              y(ir)=dreal(ut(lf,ir))
            enddo
          else
            do ir=1,nr
              y(ir)=dimag(ut(lf,ir))
            enddo
          endif
          write(20,'(f12.3,$)')t
          do ir=1,nr-1
            write(20,'(E16.8,$)')y(ir)
          enddo
          write(20,'(E16.8)')y(nr)
        enddo
        close(20)
c
c       calculate and output vertical displacement
c
        do ir=1,nr
          ymax=0.d0
          do lf=1,nf
            cswap(lf)=uz(lf,ir)-uz(0,ir)
            ymax=dmax1(ymax,cdabs(cswap(lf)))
          enddo
          if(dt.gt.0.5d0*taumin.and.cdabs(cswap(nf)).gt.0.1d0*ymax)then
            tau=dt*dmin1(2.d0,dsqrt(cdabs(cswap(nf))/ymax))
            do lf=1,nf
              lpf=(1.d0,0.d0)/((1.d0,0.d0)+cs(lf)*dcmplx(tau,0.d0))**2
              cswap(lf)=cswap(lf)*lpf
            enddo
          endif
          mf=1
          do lf=2*nf,nf+2,-1
            mf=mf+1
            cswap(lf)=dconjg(cswap(mf))
          enddo
          cswap(nf+1)=(0.d0,0.d0)
c
c         convention for Fourier transform:
c         f(t)=\int F(f) exp(i2\pi f t) df
c
          call four1(dswap,2*nf,+1)
c
          sr=dreal(uz(0,ir))
          si=sr+dswap(1)*antia(1)
          uz(1,ir)=dcmplx(sr,si)
          i=1
          t=0.d0
          do lf=2,nf
            i=i+2
            t=t+dt
            sr=si+dswap(i)*antia((i+1)/2)
            i=i+2
            t=t+dt
            si=sr+dswap(i)*antia((i+1)/2)
            uz(lf,ir)=dcmplx(sr+a,si+b)
          enddo
        enddo
c
        open(20,file=uzout,status='unknown')
        write(20,'(a,$)')' Time[day]  '
        do ir=1,nr-1
          write(20,'(a16,$)')rname(ir)
        enddo
        write(20,'(a16)')rname(nr)
        do it=1,ntcut
          t=dble(it-1)*dt/DAY2SEC
          lf=(it+1)/2
          if(it.gt.2*(it/2))then
            do ir=1,nr
              y(ir)=dreal(uz(lf,ir))
            enddo
          else
            do ir=1,nr
              y(ir)=dimag(uz(lf,ir))
            enddo
          endif
          write(20,'(f12.3,$)')t
          do ir=1,nr-1
            write(20,'(E16.8,$)')y(ir)
          enddo
          write(20,'(E16.8)')y(nr)
        enddo
        close(20)
      endif
c
      if(selvst)then
c
c       calculate and output volume strain
c
        do ir=1,nr
          ymax=0.d0
          do lf=1,nf
            cswap(lf)=upp(lf,ir)-upp(0,ir)
            ymax=dmax1(ymax,cdabs(cswap(lf)))
          enddo
          if(dt.gt.0.5d0*taumin.and.cdabs(cswap(nf)).gt.0.1d0*ymax)then
            tau=dt*dmin1(2.d0,dsqrt(cdabs(cswap(nf))/ymax))
            do lf=1,nf
              lpf=(1.d0,0.d0)/((1.d0,0.d0)+cs(lf)*dcmplx(tau,0.d0))**2
              cswap(lf)=cswap(lf)*lpf
            enddo
          endif
          mf=1
          do lf=2*nf,nf+2,-1
            mf=mf+1
            cswap(lf)=dconjg(cswap(mf))
          enddo
          cswap(nf+1)=(0.d0,0.d0)
c
c         convention for Fourier transform:
c         f(t)=\int F(f) exp(i2\pi f t) df
c
          call four1(dswap,2*nf,+1)
c
          sr=dreal(upp(0,ir))
          si=sr+dswap(1)*antia(1)
          upp(1,ir)=dcmplx(sr,si)
          i=1
          t=0.d0
          do lf=2,nf
            i=i+2
            t=t+dt
            sr=si+dswap(i)*antia((i+1)/2)
            i=i+2
            t=t+dt
            si=sr+dswap(i)*antia((i+1)/2)
            upp(lf,ir)=dcmplx(sr,si)
          enddo
        enddo
c
        open(20,file=ppout,status='unknown')
        write(20,'(a,$)')' Time[day]  '
        do ir=1,nr-1
          write(20,'(a16,$)')rname(ir)
        enddo
        write(20,'(a16)')rname(nr)
        do it=1,ntcut
          t=dble(it-1)*dt/DAY2SEC
          lf=(it+1)/2
          if(it.gt.2*(it/2))then
            do ir=1,nr
              y(ir)=dreal(upp(lf,ir))
            enddo
          else
            do ir=1,nr
              y(ir)=dimag(upp(lf,ir))
            enddo
          endif
          write(20,'(f12.3,$)')t
          do ir=1,nr-1
            write(20,'(E16.8,$)')y(ir)
          enddo
          write(20,'(E16.8)')y(nr)
        enddo
        close(20)
      endif
c
      if(selgra)then
c
c       calculate and output gravity
c
        do ir=1,nr
          ymax=0.d0
          do lf=1,nf
            cswap(lf)=ugr(lf,ir)-ugr(0,ir)
            ymax=dmax1(ymax,cdabs(cswap(lf)))
          enddo
          if(dt.gt.0.5d0*taumin.and.cdabs(cswap(nf)).gt.0.1d0*ymax)then
            tau=dt*dmin1(2.d0,dsqrt(cdabs(cswap(nf))/ymax))
            do lf=1,nf
              lpf=(1.d0,0.d0)/((1.d0,0.d0)+cs(lf)*dcmplx(tau,0.d0))**2
              cswap(lf)=cswap(lf)*lpf
            enddo
          endif
          mf=1
          do lf=2*nf,nf+2,-1
            mf=mf+1
            cswap(lf)=dconjg(cswap(mf))
          enddo
          cswap(nf+1)=(0.d0,0.d0)
c
c         convention for Fourier transform:
c         f(t)=\int F(f) exp(i2\pi f t) df
c
          call four1(dswap,2*nf,+1)
c
          sr=dreal(ugr(0,ir))
          si=sr+dswap(1)*antia(1)
          ugr(1,ir)=dcmplx(sr,si)
          i=1
          t=0.d0
          do lf=2,nf
            i=i+2
            t=t+dt
            sr=si+dswap(i)*antia((i+1)/2)
            i=i+2
            t=t+dt
            si=sr+dswap(i)*antia((i+1)/2)
            ugr(lf,ir)=dcmplx(sr,si)
          enddo
        enddo
c
        open(20,file=grout,status='unknown')
        write(20,'(a,$)')' Time[day]  '
        do ir=1,nr-1
          write(20,'(a16,$)')rname(ir)
        enddo
        write(20,'(a16)')rname(nr)
        do it=1,ntcut
          t=dble(it-1)*dt/DAY2SEC
          lf=(it+1)/2
          if(it.gt.2*(it/2))then
            do ir=1,nr
              y(ir)=dreal(ugr(lf,ir))
            enddo
          else
            do ir=1,nr
              y(ir)=dimag(ugr(lf,ir))
            enddo
          endif
          write(20,'(f12.3,$)')t
          do ir=1,nr-1
            write(20,'(E16.8,$)')y(ir)
          enddo
          write(20,'(E16.8)')y(nr)
        enddo
        close(20)
      endif
c
      if(selgeo)then
c
c       calculate and output geoid
c
        do ir=1,nr
          ymax=0.d0
          do lf=1,nf
            cswap(lf)=ugd(lf,ir)-ugd(0,ir)
            ymax=dmax1(ymax,cdabs(cswap(lf)))
          enddo
          if(dt.gt.0.5d0*taumin.and.cdabs(cswap(nf)).gt.0.1d0*ymax)then
            tau=dt*dmin1(2.d0,dsqrt(cdabs(cswap(nf))/ymax))
            do lf=1,nf
              lpf=(1.d0,0.d0)/((1.d0,0.d0)+cs(lf)*dcmplx(tau,0.d0))**2
              cswap(lf)=cswap(lf)*lpf
            enddo
          endif
          mf=1
          do lf=2*nf,nf+2,-1
            mf=mf+1
            cswap(lf)=dconjg(cswap(mf))
          enddo
          cswap(nf+1)=(0.d0,0.d0)
c
c         convention for Fourier transform:
c         f(t)=\int F(f) exp(i2\pi f t) df
c
          call four1(dswap,2*nf,+1)
c
          sr=dreal(ugd(0,ir))
          si=sr+dswap(1)*antia(1)
          ugd(1,ir)=dcmplx(sr,si)
          i=1
          t=0.d0
          do lf=2,nf
            i=i+2
            t=t+dt
            sr=si+dswap(i)*antia((i+1)/2)
            i=i+2
            t=t+dt
            si=sr+dswap(i)*antia((i+1)/2)
            ugd(lf,ir)=dcmplx(sr,si)
          enddo
        enddo
c
        open(20,file=gdout,status='unknown')
        write(20,'(a,$)')' Time[day]  '
        do ir=1,nr-1
          write(20,'(a16,$)')rname(ir)
        enddo
        write(20,'(a16)')rname(nr)
        do it=1,ntcut
          t=dble(it-1)*dt/DAY2SEC
          lf=(it+1)/2
          if(it.gt.2*(it/2))then
            do ir=1,nr
              y(ir)=dreal(ugd(lf,ir))
            enddo
          else
            do ir=1,nr
              y(ir)=dimag(ugd(lf,ir))
            enddo
          endif
          write(20,'(f12.3,$)')t
          do ir=1,nr-1
            write(20,'(E16.8,$)')y(ir)
          enddo
          write(20,'(E16.8)')y(nr)
        enddo
        close(20)
      endif
c
      if(seltil)then
c
c       calculate and output north or radial tilt
c
        do ir=1,nr
          ymax=0.d0
          do lf=1,nf
            cswap(lf)=utr(lf,ir)-utr(0,ir)
            ymax=dmax1(ymax,cdabs(cswap(lf)))
          enddo
          if(dt.gt.0.5d0*taumin.and.cdabs(cswap(nf)).gt.0.1d0*ymax)then
            tau=dt*dmin1(2.d0,dsqrt(cdabs(cswap(nf))/ymax))
            do lf=1,nf
              lpf=(1.d0,0.d0)/((1.d0,0.d0)+cs(lf)*dcmplx(tau,0.d0))**2
              cswap(lf)=cswap(lf)*lpf
            enddo
          endif
          mf=1
          do lf=2*nf,nf+2,-1
            mf=mf+1
            cswap(lf)=dconjg(cswap(mf))
          enddo
          cswap(nf+1)=(0.d0,0.d0)
c
c         convention for Fourier transform:
c         f(t)=\int F(f) exp(i2\pi f t) df
c
          call four1(dswap,2*nf,+1)
c
          sr=dreal(utr(0,ir))
          si=sr+dswap(1)*antia(1)
          utr(1,ir)=dcmplx(sr,si)
          i=1
          t=0.d0
          do lf=2,nf
            i=i+2
            t=t+dt
            sr=si+dswap(i)*antia((i+1)/2)
            i=i+2
            t=t+dt
            si=sr+dswap(i)*antia((i+1)/2)
            utr(lf,ir)=dcmplx(sr,si)
          enddo
        enddo
c
        open(20,file=trout,status='unknown')
        write(20,'(a,$)')' Time[day]  '
        do ir=1,nr-1
          write(20,'(a16,$)')rname(ir)
        enddo
        write(20,'(a16)')rname(nr)
        do it=1,ntcut
          t=dble(it-1)*dt/DAY2SEC
          lf=(it+1)/2
          if(it.gt.2*(it/2))then
            do ir=1,nr
              y(ir)=dreal(utr(lf,ir))
            enddo
          else
            do ir=1,nr
              y(ir)=dimag(utr(lf,ir))
            enddo
          endif
          write(20,'(f12.3,$)')t
          do ir=1,nr-1
            write(20,'(E16.8,$)')y(ir)
          enddo
          write(20,'(E16.8)')y(nr)
        enddo
        close(20)
c
c       calculate and output east or transversal tilt
c
        do ir=1,nr
          ymax=0.d0
          do lf=1,nf
            cswap(lf)=utt(lf,ir)-utt(0,ir)
            ymax=dmax1(ymax,cdabs(cswap(lf)))
          enddo
          if(dt.gt.0.5d0*taumin.and.cdabs(cswap(nf)).gt.0.1d0*ymax)then
            tau=dt*dmin1(2.d0,dsqrt(cdabs(cswap(nf))/ymax))
            do lf=1,nf
              lpf=(1.d0,0.d0)/((1.d0,0.d0)+cs(lf)*dcmplx(tau,0.d0))**2
              cswap(lf)=cswap(lf)*lpf
            enddo
          endif
          mf=1
          do lf=2*nf,nf+2,-1
            mf=mf+1
            cswap(lf)=dconjg(cswap(mf))
          enddo
          cswap(nf+1)=(0.d0,0.d0)
c
c         convention for Fourier transform:
c         f(t)=\int F(f) exp(i2\pi f t) df
c
          call four1(dswap,2*nf,+1)
c
          sr=dreal(utt(0,ir))
          si=sr+dswap(1)*antia(1)
          utt(1,ir)=dcmplx(sr,si)
          i=1
          t=0.d0
          do lf=2,nf
            i=i+2
            t=t+dt
            sr=si+dswap(i)*antia((i+1)/2)
            i=i+2
            t=t+dt
            si=sr+dswap(i)*antia((i+1)/2)
            utt(lf,ir)=dcmplx(sr,si)
          enddo
        enddo
c
        open(20,file=ttout,status='unknown')
        write(20,'(a,$)')' Time[day]  '
        do ir=1,nr-1
          write(20,'(a16,$)')rname(ir)
        enddo
        write(20,'(a16)')rname(nr)
        do it=1,ntcut
          t=dble(it-1)*dt/DAY2SEC
          lf=(it+1)/2
          if(it.gt.2*(it/2))then
            do ir=1,nr
              y(ir)=dreal(utt(lf,ir))
            enddo
          else
            do ir=1,nr
              y(ir)=dimag(utt(lf,ir))
            enddo
          endif
          write(20,'(f12.3,$)')t
          do ir=1,nr-1
            write(20,'(E16.8,$)')y(ir)
          enddo
          write(20,'(E16.8)')y(nr)
        enddo
        close(20)
      endif
c
c     calculate and output snapshots
c
      do isnap=1,nsnap
        open(20,file=snapshot(isnap),status='unknown')
        write(20,'(a,$)')'    Lat[deg]    Lon[deg]  Station   '
        if(ioutform.eq.1)then
          write(20,'(a)')'      U_n[m]      U_e[m]      U_z[m]'
     &     //'     Vstrain Grav[m/s^2]    Geoid[m] Tilt_n[rad]'
     &     //' Tilt_e[rad]'
        else
          write(20,'(a)')'      U_r[m]      U_t[m]      U_z[m]'
     &     //'     Vstrain Grav[m/s^2]    Geoid[m] Tilt_r[rad]'
     &     //' Tilt_t[rad]'
        endif
c
        if(tsnap(isnap).ge.0.d0)then
          it1=min0(nt,max0(1,1+idint(tsnap(isnap)/dt)))
          it2=min0(nt,1+it1)
          b=dmod(tsnap(isnap),dt)/dt
          a=1.d0-b
        else
          it1=-3
          it2=-3
          a=1.d0
          b=0.d0
        endif
        do ir=1,nr
          write(20,'(2f12.4,a12,$)')latr(ir),lonr(ir),rname(ir)
          lf=(it1+1)/2
          if(it1.ne.2*(it1/2))then
            ura=dreal(ur(lf,ir))
            uta=dreal(ut(lf,ir))
            uza=dreal(uz(lf,ir))
            ppa=dreal(upp(lf,ir))
            gra=dreal(ugr(lf,ir))
            gda=dreal(ugd(lf,ir))
            tra=dreal(utr(lf,ir))
            tta=dreal(utt(lf,ir))
          else
            ura=dimag(ur(lf,ir))
            uta=dimag(ut(lf,ir))
            uza=dimag(uz(lf,ir))
            ppa=dimag(upp(lf,ir))
            gra=dimag(ugr(lf,ir))
            gda=dimag(ugd(lf,ir))
            tra=dimag(utr(lf,ir))
            tta=dimag(utt(lf,ir))
          endif
          lf=(it2+1)/2
          if(it2.ne.2*(it2/2))then
            urb=dreal(ur(lf,ir))
            utb=dreal(ut(lf,ir))
            uzb=dreal(uz(lf,ir))
            ppb=dreal(upp(lf,ir))
            grb=dreal(ugr(lf,ir))
            gdb=dreal(ugd(lf,ir))
            trb=dreal(utr(lf,ir))
            ttb=dreal(utt(lf,ir))
          else
            urb=dimag(ur(lf,ir))
            utb=dimag(ut(lf,ir))
            uzb=dimag(uz(lf,ir))
            ppb=dimag(upp(lf,ir))
            grb=dimag(ugr(lf,ir))
            gdb=dimag(ugd(lf,ir))
            trb=dimag(utr(lf,ir))
            ttb=dimag(utt(lf,ir))
          endif
          write(20,'(8E12.4)')a*ura+b*ura,a*uta+b*utb,a*uza+b*uzb,
     &     a*ppa+b*ppb,a*gra+b*grb,a*gda+b*gdb,a*tra+b*trb,a*tta+b*ttb
        enddo
      enddo
c
      return
      end
