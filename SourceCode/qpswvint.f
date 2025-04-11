      subroutine qpswvint(ierr)
      implicit none
      integer ierr
c
      include 'qpsglobal.h'
c
      integer i,id,is,ir,ig,nd,nt0,nf0,nta0,nfa0,ishift
      integer lf,lf1,istp,ldeg,ldegf
      integer istat,ldegup
      integer ldegtap(4)
      double precision depsarc,anorm
      double precision f,dt0,df0,dta0,dfa0,rn,re,azi,bazi,bazi0
      double precision tap(0:ldegmax)
      double complex cp0,cp1,cp2
      double complex cfac,ca,cb,duz,dur,dut,dpp,dgr,dgd,dtr,dtt
      double complex expl(nsmax),clvd(nsmax),ss12(nsmax)
      double complex ss11(nsmax),ds31(nsmax),ds23(nsmax)
c
      do ldeg=0,ldegmax
        tap(ldeg)=1.d0
      enddo
c
      do is=1,ns
        depsarc=dmax1(deps(is),dpr)/rearth
c
        do ir=1,nr
          call disazi(1.d0,lats(is),lons(is),
     &                     latr(ir),lonr(ir),rn,re)
c
          dis(is,ir)=dsqrt(rn**2+re**2)
c
c         determine order of differential transform
c
          if(dis(is,ir).le.5.d0*depsarc)then
            idr(is,ir)=0
          else
            idr(is,ir)=idnint(dble(ndmax)*(dis(is,ir)-5.d0*depsarc)/PI)
          endif
c
          ssd(is,ir)=dcmplx(dsin(dis(is,ir)),0.d0)
          ssf(is,ir)=dcmplx(2.d0*dsin(0.5d0*dis(is,ir))**2,0.d0)
c
c         azi = receiver azimuth (from south to east)
c
          if(dsqrt(re*re+rn*rn).gt.0.d0)then
            azi=datan2(re,-rn)
          else
c
c           assume southern receiver in case of 0 distance
c
            azi=datan2(0.d0,1.d0)
          endif
          ssa(is,ir)=dcmplx(dsin(azi),0.d0)
          csa(is,ir)=dcmplx(dcos(azi),0.d0)
          ss2a(is,ir)=dcmplx(dsin(2.d0*azi),0.d0)
          cs2a(is,ir)=dcmplx(dcos(2.d0*azi),0.d0)
        enddo
      enddo
c
      if(ioutform.eq.1)then
        do ir=1,nr
          do is=1,ns
            call disazi(1.d0,latr(ir),lonr(ir),
     &                       lats(is),lons(is),rn,re)
            if(dsqrt(re*re+rn*rn).gt.0.d0)then
              bazi=datan2(-re,-rn)
            else
              bazi=datan2(0.d0,-1.d0)
            endif
            ssb(is,ir)=dcmplx(dsin(bazi),0.d0)
            csb(is,ir)=dcmplx(dcos(bazi),0.d0)
          enddo
        enddo
      else
        do ir=1,nr
          call disazi(1.d0,latr(ir),lonr(ir),
     &                     epilat,epilon,rn,re)
          if(dsqrt(re*re+rn*rn).gt.0.d0)then
            bazi0=datan2(-re,-rn)
          else
            bazi0=datan2(0.d0,-1.d0)
          endif
          do is=1,ns
            call disazi(1.d0,latr(ir),lonr(ir),
     &                       lats(is),lons(is),rn,re)
            if(dsqrt(re*re+rn*rn).gt.0.d0)then
              bazi=datan2(-re,-rn)
            else
              bazi=datan2(0.d0,-1.d0)
            endif
            ssb(is,ir)=dcmplx(dsin(bazi-bazi0),0.d0)
            csb(is,ir)=dcmplx(dcos(bazi-bazi0),0.d0)
          enddo
        enddo
      endif
c
      do is=1,ns
        expl(is)=dcmplx((mtt(is)+mpp(is)+mrr(is))/3.d0,0.d0)
        clvd(is)=dcmplx(mrr(is),0.d0)-expl(is)
        ss12(is)=dcmplx(mtp(is),0.d0)
        ss11(is)=dcmplx((mtt(is)-mpp(is))/2.d0,0.d0)
        ds31(is)=dcmplx(mrt(is),0.d0)
        ds23(is)=dcmplx(mpr(is),0.d0)
      enddo
c
      do lf=-1,nf
        do ir=1,nr
          uz(lf,ir)=(0.d0,0.d0)
          ur(lf,ir)=(0.d0,0.d0)
          ut(lf,ir)=(0.d0,0.d0)
          upp(lf,ir)=(0.d0,0.d0)
          ugr(lf,ir)=(0.d0,0.d0)
          ugd(lf,ir)=(0.d0,0.d0)
          utr(lf,ir)=(0.d0,0.d0)
          utt(lf,ir)=(0.d0,0.d0)
        enddo
      enddo
c
      do ig=1,ngrn
        if(nsg(ig).le.0)goto 500
        write(*,'(a)')' '
        write(*,'(a,i4,a,f5.1,a)')' processing ',1+isg2(ig)-isg1(ig),
     &    ' point source(s) at depth ',grndep(ig)/KM2M,' km'
        write(*,'(a)')' open Green function data base: '
     &              //grnfile(ig)(1:40)
        write(*,'(a)')' ... please wait ...'
c
        open(21,file=uzgrnfile(ig),
     &       form='unformatted',status='old')
        open(22,file=urgrnfile(ig),
     &       form='unformatted',status='old')
        open(23,file=utgrnfile(ig),
     &       form='unformatted',status='old')
        open(24,file=ppgrnfile(ig),
     &       form='unformatted',status='old')
        open(25,file=grgrnfile(ig),
     &       form='unformatted',status='old')
        open(26,file=gdgrnfile(ig),
     &       form='unformatted',status='old')
        open(27,file=trgrnfile(ig),
     &       form='unformatted',status='old')
        open(28,file=ttgrnfile(ig),
     &       form='unformatted',status='old')
c
        read(21)nt0,dt0,nf0,df0,ldegup
        read(22)nt0,dt0,nf0,df0,ldegup
        read(23)nt0,dt0,nf0,df0,ldegup
        read(24)nt0,dt0,nf0,df0,ldegup
        read(25)nt0,dt0,nf0,df0,ldegup
        read(26)nt0,dt0,nf0,df0,ldegup
        read(27)nt0,dt0,nf0,df0,ldegup
        read(28)nt0,dt0,nf0,df0,ldegup
c
        if(dt0.ne.dt.or.df0.ne.df.or.
     &     nt0.ne.nt.or.nf0.ne.nf)then
          print *,' Error in qpswvint: t/f sampling'
     &          //' inconsistent with Green functions!'
          write(*,'(a)')'                    nt       nf'
     &                //'       dt       df'
          write(*,'(a,2i5,2f16.8)')' Current input:',nt,nf,dt,df
          write(*,'(a,2i5,2f16.8)')'     Data base:',nt0,nf0,dt0,df0
          stop
        endif
c
        nd=0
        do is=isg1(ig),isg2(ig)
          do ir=1,nr
            nd=max0(nd,idr(is,ir))
          enddo
        enddo
c
        if(ldegup+nd.gt.ldegmax)then
          write(*,'(a,i6,a,i6)')' Error in qpswvint: '
     &     //'max. harmonic degree required = ',ldegup+nd,
     &     ' > ldegmax defined: ',ldegmax
          stop
        endif
c
        do lf=-1,nf
          f=dmax1(0.d0,dble(lf-1)*df)
c
          read(21)ldegf
          read(22)ldegf
          read(23)ldegf
          read(24)ldegf
          read(25)ldegf
          read(26)ldegf
          read(27)ldegf
          read(28)ldegf
c
          read(21)((yuz(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
          read(22)((yur(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
          read(23)((yut(ldeg,istp,0),ldeg=0,ldegf),istp=2,3)
          read(24)((ypp(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
          read(25)((ygr(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
          read(26)((ygd(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
          read(27)((ytr(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
          read(28)((ytt(ldeg,istp,0),ldeg=0,ldegf),istp=2,3)
c
          ldegtap(1)=0
          ldegtap(2)=0
          ldegtap(3)=ldegf*4/5
          ldegtap(4)=ldegf
c
c         use differential filter to suppress spatial aliasing
c
          do id=0,nd
            yur(0,1,id)=(0.d0,0.d0)
            yur(0,4,id)=(0.d0,0.d0)
            ytr(0,1,id)=(0.d0,0.d0)
            ytr(0,4,id)=(0.d0,0.d0)
c
            yuz(0,3,id)=(0.d0,0.d0)
            yur(0,3,id)=(0.d0,0.d0)
            yut(0,3,id)=(0.d0,0.d0)
            ypp(0,3,id)=(0.d0,0.d0)
            ygr(0,3,id)=(0.d0,0.d0)
            ygd(0,3,id)=(0.d0,0.d0)
            ytr(0,3,id)=(0.d0,0.d0)
            ytt(0,3,id)=(0.d0,0.d0)
c
            yuz(0,2,id)=(0.d0,0.d0)
            yur(0,2,id)=(0.d0,0.d0)
            yut(0,2,id)=(0.d0,0.d0)
            ypp(0,2,id)=(0.d0,0.d0)
            ygr(0,2,id)=(0.d0,0.d0)
            ygd(0,2,id)=(0.d0,0.d0)
            ytr(0,2,id)=(0.d0,0.d0)
            ytt(0,2,id)=(0.d0,0.d0)
c
            yuz(1,2,id)=(0.d0,0.d0)
            yur(1,2,id)=(0.d0,0.d0)
            yut(1,2,id)=(0.d0,0.d0)
            ypp(1,2,id)=(0.d0,0.d0)
            ygr(1,2,id)=(0.d0,0.d0)
            ygd(1,2,id)=(0.d0,0.d0)
            ytr(1,2,id)=(0.d0,0.d0)
            ytt(1,2,id)=(0.d0,0.d0)
          enddo
c
          do id=1,nd
c
c           m = 0
c
            ca=dcmplx(1.d0/3.d0,0.d0)
            yuz(0,1,id)=yuz(0,1,id-1)-ca*yuz(1,1,id-1)
            ypp(0,1,id)=ypp(0,1,id-1)-ca*ypp(1,1,id-1)
            ygr(0,1,id)=ygr(0,1,id-1)-ca*ygr(1,1,id-1)
            ygd(0,1,id)=ygd(0,1,id-1)-ca*ygd(1,1,id-1)
            yuz(0,4,id)=yuz(0,4,id-1)-ca*yuz(1,4,id-1)
            ypp(0,4,id)=ypp(0,4,id-1)-ca*ypp(1,4,id-1)
            ygr(0,4,id)=ygr(0,4,id-1)-ca*ygr(1,4,id-1)
            ygd(0,4,id)=ygd(0,4,id-1)-ca*ygd(1,4,id-1)
            do ldeg=1,ldegf-id
              ca=dcmplx(dble(ldeg+1)/dble(2*ldeg+3),0.d0)
              cb=dcmplx(dble(ldeg)/dble(2*ldeg-1),0.d0)
              yuz(ldeg,1,id)=yuz(ldeg,1,id-1)-ca*yuz(ldeg+1,1,id-1)
     &                     -cb*yuz(ldeg-1,1,id-1)
              ypp(ldeg,1,id)=ypp(ldeg,1,id-1)-ca*ypp(ldeg+1,1,id-1)
     &                     -cb*ypp(ldeg-1,1,id-1)
              ygr(ldeg,1,id)=ygr(ldeg,1,id-1)-ca*ygr(ldeg+1,1,id-1)
     &                     -cb*ygr(ldeg-1,1,id-1)
              ygd(ldeg,1,id)=ygd(ldeg,1,id-1)-ca*ygd(ldeg+1,1,id-1)
     &                     -cb*ygd(ldeg-1,1,id-1)
              yuz(ldeg,4,id)=yuz(ldeg,4,id-1)-ca*yuz(ldeg+1,4,id-1)
     &                     -cb*yuz(ldeg-1,4,id-1)
              ypp(ldeg,4,id)=ypp(ldeg,4,id-1)-ca*ypp(ldeg+1,4,id-1)
     &                     -cb*ypp(ldeg-1,4,id-1)
              ygr(ldeg,4,id)=ygr(ldeg,4,id-1)-ca*ygr(ldeg+1,4,id-1)
     &                     -cb*ygr(ldeg-1,4,id-1)
              ygd(ldeg,4,id)=ygd(ldeg,4,id-1)-ca*ygd(ldeg+1,4,id-1)
     &                     -cb*ygd(ldeg-1,4,id-1)
            enddo
c
c           m = 1
c
            do ldeg=1,ldegf-id
              ca=dcmplx(dble(ldeg+2)/dble(2*ldeg+3),0.d0)
              cb=dcmplx(dble(ldeg-1)/dble(2*ldeg-1),0.d0)
              yur(ldeg,1,id)=yur(ldeg,1,id-1)-ca*yur(ldeg+1,1,id-1)
     &                     -cb*yur(ldeg-1,1,id-1)
              yur(ldeg,4,id)=yur(ldeg,4,id-1)-ca*yur(ldeg+1,4,id-1)
     &                     -cb*yur(ldeg-1,4,id-1)
              ytr(ldeg,1,id)=ytr(ldeg,1,id-1)-ca*ytr(ldeg+1,1,id-1)
     &                     -cb*ytr(ldeg-1,1,id-1)
              ytr(ldeg,4,id)=ytr(ldeg,4,id-1)-ca*ytr(ldeg+1,4,id-1)
     &                     -cb*ytr(ldeg-1,4,id-1)
              yuz(ldeg,3,id)=yuz(ldeg,3,id-1)-ca*yuz(ldeg+1,3,id-1)
     &                     -cb*yuz(ldeg-1,3,id-1)
              yur(ldeg,3,id)=yur(ldeg,3,id-1)-ca*yur(ldeg+1,3,id-1)
     &                     -cb*yur(ldeg-1,3,id-1)
              yut(ldeg,3,id)=yut(ldeg,3,id-1)-ca*yut(ldeg+1,3,id-1)
     &                     -cb*yut(ldeg-1,3,id-1)
              ypp(ldeg,3,id)=ypp(ldeg,3,id-1)-ca*ypp(ldeg+1,3,id-1)
     &                     -cb*ypp(ldeg-1,3,id-1)
              ygr(ldeg,3,id)=ygr(ldeg,3,id-1)-ca*ygr(ldeg+1,3,id-1)
     &                     -cb*ygr(ldeg-1,3,id-1)
              ygd(ldeg,3,id)=ygd(ldeg,3,id-1)-ca*ygd(ldeg+1,3,id-1)
     &                     -cb*ygd(ldeg-1,3,id-1)
              ytr(ldeg,3,id)=ytr(ldeg,3,id-1)-ca*ytr(ldeg+1,3,id-1)
     &                     -cb*ytr(ldeg-1,3,id-1)
              ytt(ldeg,3,id)=ytt(ldeg,3,id-1)-ca*ytt(ldeg+1,3,id-1)
     &                     -cb*ytt(ldeg-1,3,id-1)
            enddo
c
c           m = 2
c
            do ldeg=2,ldegf-id
              ca=dcmplx(dble(ldeg+3)/dble(2*ldeg+3),0.d0)
              cb=dcmplx(dble(ldeg-2)/dble(2*ldeg-1),0.d0)
              yuz(ldeg,2,id)=yuz(ldeg,2,id-1)-ca*yuz(ldeg+1,2,id-1)
     &                     -cb*yuz(ldeg-1,2,id-1)
              yur(ldeg,2,id)=yur(ldeg,2,id-1)-ca*yur(ldeg+1,2,id-1)
     &                     -cb*yur(ldeg-1,2,id-1)
              yut(ldeg,2,id)=yut(ldeg,2,id-1)-ca*yut(ldeg+1,2,id-1)
     &                     -cb*yut(ldeg-1,2,id-1)
              ypp(ldeg,2,id)=ypp(ldeg,2,id-1)-ca*ypp(ldeg+1,2,id-1)
     &                     -cb*ypp(ldeg-1,2,id-1)
              ygr(ldeg,2,id)=ygr(ldeg,2,id-1)-ca*ygr(ldeg+1,2,id-1)
     &                     -cb*ygr(ldeg-1,2,id-1)
              ygd(ldeg,2,id)=ygd(ldeg,2,id-1)-ca*ygd(ldeg+1,2,id-1)
     &                     -cb*ygd(ldeg-1,2,id-1)
              ytr(ldeg,2,id)=ytr(ldeg,2,id-1)-ca*ytr(ldeg+1,2,id-1)
     &                     -cb*ytr(ldeg-1,2,id-1)
              ytt(ldeg,2,id)=ytt(ldeg,2,id-1)-ca*ytt(ldeg+1,2,id-1)
     &                     -cb*ytt(ldeg-1,2,id-1)
            enddo
          enddo
c
          do is=isg1(ig),isg2(ig)
            do ir=1,nr
              id=idr(is,ir)
c
              call taper(ldegtap(1),ldegtap(4),tap(0))
c
              call qpslegendre(ldegtap(4),dis(is,ir))
c
              do ldeg=0,ldegf
                cp0=dcmplx(plm(ldeg,0),0.d0)
                cp1=dcmplx(plm(ldeg,1),0.d0)
                cp2=dcmplx(plm(ldeg,2),0.d0)
                duz=(expl(is)*yuz(ldeg,1,id)
     &             +clvd(is)*yuz(ldeg,4,id))*cp0
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &             *yuz(ldeg,2,id)*cp2*ssd(is,ir)**2
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &             *yuz(ldeg,3,id)*cp1*ssd(is,ir)
                dpp=(expl(is)*ypp(ldeg,1,id)
     &             +clvd(is)*ypp(ldeg,4,id))*cp0
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &             *ypp(ldeg,2,id)*cp2*ssd(is,ir)**2
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &             *ypp(ldeg,3,id)*cp1*ssd(is,ir)
                dgr=(expl(is)*ygr(ldeg,1,id)
     &             +clvd(is)*ygr(ldeg,4,id))*cp0
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &             *ygr(ldeg,2,id)*cp2*ssd(is,ir)**2
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &             *ygr(ldeg,3,id)*cp1*ssd(is,ir)
                dgd=(expl(is)*ygd(ldeg,1,id)
     &             +clvd(is)*ygd(ldeg,4,id))*cp0
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &             *ygd(ldeg,2,id)*cp2*ssd(is,ir)**2
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &             *ygd(ldeg,3,id)*cp1*ssd(is,ir)
                dur=(expl(is)*yur(ldeg,1,id)+clvd(is)*yur(ldeg,4,id))
     &             *cp1*ssd(is,ir)
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &             *yur(ldeg,2,id)*cp2*ssd(is,ir)
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &             *yur(ldeg,3,id)*cp1
                dut=(ss12(is)*cs2a(is,ir)-ss11(is)*ss2a(is,ir))
     &             *yut(ldeg,2,id)*cp2*ssd(is,ir)
     &             +(ds31(is)*ssa(is,ir)-ds23(is)*csa(is,ir))
     &             *yut(ldeg,3,id)*cp1
                dtr=(expl(is)*ytr(ldeg,1,id)+clvd(is)*ytr(ldeg,4,id))
     &             *cp1*ssd(is,ir)
     &             +(ss12(is)*ss2a(is,ir)+ss11(is)*cs2a(is,ir))
     &             *ytr(ldeg,2,id)*cp2*ssd(is,ir)
     &             +(ds31(is)*csa(is,ir)+ds23(is)*ssa(is,ir))
     &             *ytr(ldeg,3,id)*cp1
                dtt=(ss12(is)*cs2a(is,ir)-ss11(is)*ss2a(is,ir))
     &             *ytt(ldeg,2,id)*cp2*ssd(is,ir)
     &             +(ds31(is)*ssa(is,ir)-ds23(is)*csa(is,ir))
     &             *ytt(ldeg,3,id)*cp1
                cfac=dcmplx(tap(ldeg),0.d0)
                if(id.gt.0)cfac=cfac/ssf(is,ir)**id
                uz(lf,ir)=uz(lf,ir)+duz*cfac
                upp(lf,ir)=upp(lf,ir)+dpp*cfac
                ugr(lf,ir)=ugr(lf,ir)+dgr*cfac
                ugd(lf,ir)=ugd(lf,ir)+dgd*cfac
                if(ioutform.eq.1)then
                  ur(lf,ir)=ur(lf,ir)
     &                     +(dur*csb(is,ir)+dut*ssb(is,ir))*cfac
                  ut(lf,ir)=ut(lf,ir)
     &                     +(dur*ssb(is,ir)-dut*csb(is,ir))*cfac
                  utr(lf,ir)=utr(lf,ir)
     &                     +(dtr*csb(is,ir)+dtt*ssb(is,ir))*cfac
                  utt(lf,ir)=utt(lf,ir)
     &                     +(dtr*ssb(is,ir)-dtt*csb(is,ir))*cfac
                else
                  ur(lf,ir)=ur(lf,ir)
     &                     +(dur*csb(is,ir)+dut*ssb(is,ir))*cfac
                  ut(lf,ir)=ut(lf,ir)
     &                     -(dur*ssb(is,ir)-dut*csb(is,ir))*cfac
                  utr(lf,ir)=utr(lf,ir)
     &                     +(dtr*csb(is,ir)+dtt*ssb(is,ir))*cfac
                  utt(lf,ir)=utt(lf,ir)
     &                     -(dtr*ssb(is,ir)-dtt*csb(is,ir))*cfac
                endif
              enddo
            enddo
          enddo
          if(lf.eq.-1)then
            write(*,'(a,i5)')'  fully relaxed response: spectra read: ',
     &                       ldegf
          else if(lf.eq.0)then
            write(*,'(a,i5)')'     co-seismic response: spectra read: ',
     &                       ldegf
          else
            write(*,'(i6,a,E14.6,a,i5)')lf,'.',f,
     &                        ' Hz: spectra read: ',ldegf
          endif
        enddo
        close(21)
        close(22)
        close(23)
        close(24)
        close(25)
        close(26)
        close(27)
        close(28)
        write(*,'(i6,a)')lf-1,' spectra read from '
     &                      //grnfile(ig)(1:40)
500     continue
      enddo
      return
      end
