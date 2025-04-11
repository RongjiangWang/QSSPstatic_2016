      subroutine qpsgetinp(unit)
      implicit none
      integer unit
c
      include 'qpsglobal.h'
c
c     work space
c
      integer i,j,k,l,ir,ig,isg,is,is1,flen,sdfsel
      integer icmp(5)
      double precision twindow,suppress,munit,fac,amf
      double precision strike,dip,rake,depdif,dswap(9)
      double precision wup,wlw,deprsmax,ddep,delta,vp0,vs0
      double complex epicenter
      character*80 grndir,outfile(5),fswap
      logical hooke,maxwell,sls,burgers
c
c     uniform receiver depth
c     ======================
c
      lyadd=1
c
      call skipdoc(unit)
      read(unit,*)dpr
      dpr=KM2M*dpr
c
c     time (frequency) sampling
c     =========================
c
      call skipdoc(unit)
      read(unit,*)twindow,dt
      twindow=twindow*DAY2SEC
      dt=dt*DAY2SEC
      ntcut=1+idnint(twindow/dt)
      nt=1
100   nt=2*nt
      if(nt.lt.ntcut)goto 100
      nf=nt/2
      if(nf.gt.nfmax)then
        stop ' Error in qpsgetinp: nfmax (qpsglobal.h) too small!'
      endif
      df=1.d0/(dble(nt)*dt)
c
      call skipdoc(unit)
      read(unit,*)suppress
      if(suppress.le.0.d0.or.suppress.ge.1.d0)then
        suppress=dexp(-1.d0)
      endif
      fi=dlog(suppress)*df/(2.d0*PI)
      call skipdoc(unit)
      read(unit,*)rearth,grsurf
      rearth=rearth*KM2M
c
c     cutoffs of spectra
c     ==================
c
      call skipdoc(unit)
      read(unit,*)ldeggr
      if(ldeggr.lt.0)ldeggr=0
      nogravity=ldeggr.le.0.d0
c
      call skipdoc(unit)
      read(unit,*)ldegcut
      if(ldegcut.gt.ldegmax)then
        stop ' Error in qpsgetinp: max. cutoff degree exceeds limit!'
      endif
      ldegcut=min0(ldegmax-ndmax-1,ldegcut)
c
c     Green's function files
c     ======================
c
      call skipdoc(unit)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     option of stress-drop source deleted on March 21, 2017
c     read(unit,*)ngrn,istype,rexfault,rr0,grndir
      read(unit,*)ngrn,rr0,grndir
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(ngrn.le.0)then
        stop ' Error in qpsgetinp: bad number of source depths!'
      else if(ngrn.gt.ngrnmax)then
        stop ' Error in qpsgetinp: too large number of source depths!'
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     option of stress-drop source deleted on March 21, 2017
c
c
c     if(istype.lt.0.or.istype.gt.1)then
c       stop ' Error in qpsgetinp: bad Green function source type!'
c     endif
c     rexfault=dmax1(0.d0,dmin1(1.d0,rexfault))
      istype=0
      rexfault=1.d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      rr0=rr0*KM2M
      lyadd=lyadd+ngrn
c
      deprsmax=dpr
      do ig=1,ngrn
        call skipdoc(unit)
        read(unit,*)grndep(ig),grnfile(ig),grnsel(ig)
        if(grnsel(ig).lt.0.or.grnsel(ig).gt.1)then
          stop ' Error in qpsgetinp: bad Green function selection!'
        endif
        grndep(ig)=grndep(ig)*KM2M
        deprsmax=dmax1(deprsmax,grndep(ig))
      enddo
c
c     sort green function files by source depth
c
      do i=1,ngrn
        do j=i+1,ngrn
          if(grndep(j).lt.grndep(i))then
            dswap(1)=grndep(i)
            fswap=grnfile(i)
            k=grnsel(i)
c
            grndep(i)=grndep(j)
            grnfile(i)=grnfile(j)
            grnsel(i)=grnsel(j)
c
            grndep(j)=dswap(1)
            grnfile(j)=fswap
            grnsel(j)=k
          endif
        enddo
      enddo
c
      do flen=80,1,-1
        if(grndir(flen:flen).ne.' ')goto 200
      enddo
200   continue
      do ig=1,ngrn
        uzgrnfile(ig)=grndir(1:flen)//'Uz_'//grnfile(ig)
        urgrnfile(ig)=grndir(1:flen)//'Ur_'//grnfile(ig)
        utgrnfile(ig)=grndir(1:flen)//'Ut_'//grnfile(ig)
        ppgrnfile(ig)=grndir(1:flen)//'Pp_'//grnfile(ig)
        grgrnfile(ig)=grndir(1:flen)//'Gr_'//grnfile(ig)
        gdgrnfile(ig)=grndir(1:flen)//'Gd_'//grnfile(ig)
        trgrnfile(ig)=grndir(1:flen)//'Tr_'//grnfile(ig)
        ttgrnfile(ig)=grndir(1:flen)//'Tt_'//grnfile(ig)
      enddo
c
c     multi-event source parameters
c     =============================
c
      call skipdoc(unit)
      read(unit,*)ns,sdfsel
      if(ns.gt.nsmax)then
        stop ' Error in qpsgetinp: too many subevents'
      endif
      if(sdfsel.eq.1)then
        do is=1,ns
          call skipdoc(unit)
c
c         the six moment-tensor elements: Mrr, Mtt, Mpp, Mrt, Mrp, Mtp
c
          read(unit,*)munit,mrr(is),mtt(is),mpp(is),
     &                      mrt(is),mpr(is),mtp(is),
     &                      lats(is),lons(is),deps(is)
          mtt(is)=mtt(is)*munit
          mpp(is)=mpp(is)*munit
          mrr(is)=mrr(is)*munit
          mtp(is)=mtp(is)*munit
          mpr(is)=mpr(is)*munit
          mrt(is)=mrt(is)*munit
          deps(is)=deps(is)*KM2M
        enddo
      else if(sdfsel.eq.2)then
        do is=1,ns
          call skipdoc(unit)
          read(unit,*)munit,strike,dip,rake,
     &                      lats(is),lons(is),deps(is)
          call moments(munit,strike,dip,rake,
     &                 mtt(is),mpp(is),mrr(is),
     &                 mtp(is),mpr(is),mrt(is))
          deps(is)=deps(is)*KM2M
        enddo
      else
        stop ' Error in qpsgetinp: bad selection source data format!'
      endif
c
c     sort sub-events by depth
c
      do i=1,ns
        do j=i+1,ns
          if(deps(j).lt.deps(i))then
            dswap(1)=lats(i)
            dswap(2)=lons(i)
            dswap(3)=deps(i)
            dswap(4)=mtt(i)
            dswap(5)=mpp(i)
            dswap(6)=mrr(i)
            dswap(7)=mtp(i)
            dswap(8)=mpr(i)
            dswap(9)=mrt(i)
c
            lats(i)=lats(j)
            lons(i)=lons(j)
            deps(i)=deps(j)
            mtt(i)=mtt(j)
            mpp(i)=mpp(j)
            mrr(i)=mrr(j)
            mtp(i)=mtp(j)
            mpr(i)=mpr(j)
            mrt(i)=mrt(j)
c
            lats(j)=dswap(1)
            lons(j)=dswap(2)
            deps(j)=dswap(3)
            mtt(j)=dswap(4)
            mpp(j)=dswap(5)
            mrr(j)=dswap(6)
            mtp(j)=dswap(7)
            mpr(j)=dswap(8)
            mrt(j)=dswap(9)
          endif
        enddo
      enddo
c
      isg1(1)=1
      is1=1
      do ig=1,ngrn-1
        depdif=0.5d0*(grndep(ig)+grndep(ig+1))
        isg2(ig)=isg1(ig)-1
        do is=is1,ns
          if(deps(is).lt.depdif)then
            isg2(ig)=is
          endif
        enddo
        isg1(ig+1)=isg2(ig)+1
        is1=isg1(ig+1)
      enddo
      isg2(ngrn)=ns
c
      do ig=1,ngrn
        nsg(ig)=max0(0,1+isg2(ig)-isg1(ig))
      enddo
c
c     receiver parameters
c     ===================
c
      call skipdoc(unit)
      read(unit,*)ioutform,epilat,epilon
      if(ioutform.lt.1.or.ioutform.gt.2)then
        stop ' Error in qpsgetinp: bad selection of output format!'
      endif
c
      call skipdoc(unit)
      read(unit,*)(icmp(i),i=1,5)
      seldis=icmp(1).eq.1
      selvst=icmp(2).eq.1
      seltil=icmp(3).eq.1
      selgra=icmp(4).eq.1
      selgeo=icmp(5).eq.1
c
      call skipdoc(unit)
      read(unit,*)(outfile(i),i=1,5)
c
      do flen=80,1,-1
        if(outfile(1)(flen:flen).ne.' ')goto 301
      enddo
301   continue
      if(ioutform.eq.1)then
        urout=outfile(1)(1:flen)//'.un'
        utout=outfile(1)(1:flen)//'.ue'
        uzout=outfile(1)(1:flen)//'.uz'
      else
        urout=outfile(1)(1:flen)//'.ur'
        utout=outfile(1)(1:flen)//'.ut'
        uzout=outfile(1)(1:flen)//'.uz'
      endif
c
      do flen=80,1,-1
        if(outfile(2)(flen:flen).ne.' ')goto 302
      enddo
302   continue
      ppout=outfile(2)(1:flen)//'.vs'
c
      do flen=80,1,-1
        if(outfile(3)(flen:flen).ne.' ')goto 303
      enddo
303   continue
      if(ioutform.eq.1)then
        trout=outfile(3)(1:flen)//'.tn'
        ttout=outfile(3)(1:flen)//'.te'
      else
        trout=outfile(3)(1:flen)//'.tr'
        ttout=outfile(3)(1:flen)//'.tt'
      endif
c
      do flen=80,1,-1
        if(outfile(4)(flen:flen).ne.' ')goto 304
      enddo
304   continue
      grout=outfile(4)(1:flen)//'.gr'
c
      do flen=80,1,-1
        if(outfile(5)(flen:flen).ne.' ')goto 305
      enddo
305   continue
      gdout=outfile(5)(1:flen)//'.gd'
c
      call skipdoc(unit)
      read(unit,*)nsnap
      if(nsnap.gt.nsnapmax)then
        stop ' Error in qpsgetinp: too many snapshots!'
      endif
      do i=1,nsnap
        call skipdoc(unit)
        read(unit,*)tsnap(i),snapshot(i)
        tsnap(i)=tsnap(i)*DAY2SEC
        if(tsnap(i).gt.twindow)then
          stop ' Error in qpsgetinp: bad snapshot time!'
        endif
      enddo
c
      call skipdoc(unit)
      read(unit,*)nr
      if(nr.gt.nrmax)then
        stop ' Error in qpsgetinp: too many observation stations!'
      endif
c
      do ir=1,nr
        call skipdoc(unit)
        read(unit,*)latr(ir),lonr(ir),rname(ir)
      enddo
c
c     multilayered model parameters
c     =============================
c
      call skipdoc(unit)
      read(unit,*)l,modelmod
      if(l.ge.2*lymax-lyadd)then
        stop ' Error in qpsgetinp: lymax defined too small!'
      endif
c
      taumin=dble(nt)*dt
      do i=1,l
        call skipdoc(unit)
        read(unit,*)j,dp0(i),vp0,vs0,rho0(i),
     &              eta10(i),eta20(i),alfa0(i)
c
c       input units:    -,km,  km/s, km/s, g/cm^3,-,-
c
        if(vp0.le.0.d0.or.vs0.lt.0.d0.or.rho0(i).le.0.d0)then
          stop ' Error in qpsgetinp: bad seismic parameter!'
        endif
c
c       1 = Hooke, 2 = Maxwell, 3 = SLS, 4 = Burgers
c
        hooke=eta10(i).le.0.d0.and.eta20(i).le.0.d0
        maxwell=eta20(i).gt.0.d0.and.
     &          (eta10(i).le.0.d0.or.alfa0(i).eq.1.d0)
        sls=eta10(i).gt.0.d0.and.eta20(i).le.0.d0.and.
     &      alfa0(i).gt.0.d0.and.alfa0(i).lt.1.d0
        burgers=eta10(i).gt.0.d0.and.eta20(i).gt.0.d0.and.
     &          alfa0(i).gt.0.d0.and.alfa0(i).lt.1.d0
        if(.not.(hooke.or.maxwell.or.sls.or.burgers))then
          stop ' Error in qpsgetinp: not identified viscoelastic body!'
        endif
        eta10(i)=dmax1(0.d0,eta10(i))
        eta20(i)=dmax1(0.d0,eta20(i))
c
        dp0(i)=KM2M*dp0(i)
        vp0=KM2M*vp0
        vs0=KM2M*vs0
        rho0(i)=KM2M*rho0(i)
        mue0(i)=rho0(i)*vs0**2
        kap0(i)=rho0(i)*vp0**2-mue0(i)*4.d0/3.d0
        if(i.gt.1)then
          if(dp0(i).lt.dp0(i-1))then
            stop ' Error in qpsgetinp: bad layering of earth model!'
          endif
        endif
        if(kap0(i).le.0.d0)then
          stop ' Error in qpsgetinp: bad Vp/Vs ratio!'
        endif
        if(mue0(i).gt.0.d0)then
          if(maxwell)then
            taumin=dmin1(taumin,eta20(i)/mue0(i))
          else if(burgers)then
            taumin=dmin1(taumin,eta20(i)/mue0(i),
     &                   eta10(i)*(1.d0-alfa0(i))/mue0(i)/alfa0(i))
          endif
        endif
      enddo
c
      if(dp0(1).lt.0.d0.or.dp0(1).gt.0.d0)then
       stop ' Error in qpsgetinp: wrong first interface depth!'
      else if(dp0(l).gt.rearth)then
        stop ' Error in qpsgetinp: too large interface depth!'
      else if(dp0(l).lt.rearth)then
        l=l+1
        if(l.ge.2*lymax-lyadd)then
          stop ' Error in qpsgetinp: lymax defined too small!'
        endif
        dp0(l)=rearth
        kap0(l)=kap0(l-1)
        mue0(l)=mue0(l-1)
        rho0(l)=rho0(l-1)
        eta10(l)=eta10(l-1)
        eta20(l)=eta20(l-1)
        alfa0(l)=alfa0(l-1)
      endif    
c
      if(dpr.lt.0.d0.or.dpr.gt.dp0(l))then
        stop ' Error in qpsgetinp: receiver too shallow or too deep!'
      endif
      do i=1,ngrn
        if(grndep(i).lt.0.d0.or.grndep(i).ge.dp0(l))then
          stop ' Error in qpsgetinp: source too shallow or too deep!'
        endif
      enddo
c
      l0=0
      do i=2,l
        if(dp0(i).gt.dp0(i-1))then
          l0=l0+1
          dp0up(l0)=dp0(i-1)
          kap0up(l0)=kap0(i-1)
          mue0up(l0)=mue0(i-1)
          rho0up(l0)=rho0(i-1)
          eta10up(l0)=eta10(i-1)
          eta20up(l0)=eta20(i-1)
          alfa0up(l0)=alfa0(i-1)
c
          dp0lw(l0)=dp0(i)
          kap0lw(l0)=kap0(i)
          mue0lw(l0)=mue0(i)
          rho0lw(l0)=rho0(i)
          eta10lw(l0)=eta10(i)
          eta20lw(l0)=eta20(i)
          alfa0lw(l0)=alfa0(i)
c
          if(eta1up(l0).le.0.d0.and.eta1lw(l0).gt.0.d0.or.
     &       eta1up(l0).gt.0.d0.and.eta1lw(l0).le.0.d0.or.
     &       eta2up(l0).le.0.d0.and.eta2lw(l0).gt.0.d0.or.
     &       eta2up(l0).gt.0.d0.and.eta2lw(l0).le.0.d0)then
            print *, ' Error in qpsgetinp:'
     &             //' nonunique viscoelastic body within a layer!'
            stop
          endif
        endif
      enddo
c
      ddep=dp0lw(l0)-dp0up(l0)
      if(dp0up(l0).le.deprsmax+THICKMAX)then
        l0=l0+1
        if(l0.ge.lymax-lyadd)then
          stop ' Error in qpsgetinp: lymax defined too small!'
        endif
        if(deprsmax.ge.dp0up(l0-1))then
          delta=deprsmax-dp0up(l0-1)
     &         +dmin1(THICKMAX,0.5d0*(dp0up(l0-1)-deprsmax))
        else
          delta=dmin1(THICKMAX,0.5d0*ddep)
        endif
c
        wlw=delta/ddep
        wup=1.d0-wlw
c
        dp0lw(l0)=dp0lw(l0-1)
        kap0lw(l0)=kap0lw(l0-1)
        mue0lw(l0)=mue0lw(l0-1)
        rho0lw(l0)=rho0lw(l0-1)
        eta10lw(l0)=eta10lw(l0-1)
        eta20lw(l0)=eta20lw(l0-1)
        alfa0lw(l0)=alfa0lw(l0-1)
c
        dp0up(l0)=dp0up(l0-1)+delta
        kap0up(l0)=wup*kap0up(l0-1)+wlw*kap0lw(l0-1)
        mue0up(l0)=wup*mue0up(l0-1)+wlw*mue0lw(l0-1)
        rho0up(l0)=wup*rho0up(l0-1)+wlw*rho0lw(l0-1)
        eta10up(l0)=wup*eta10up(l0-1)+wlw*eta10lw(l0-1)
        eta20up(l0)=wup*eta20up(l0-1)+wlw*eta20lw(l0-1)
        alfa0up(l0)=wup*alfa0up(l0-1)+wlw*alfa0lw(l0-1)
c
        dp0lw(l0-1)=dp0up(l0)
        kap0lw(l0-1)=kap0up(l0)
        mue0lw(l0-1)=mue0up(l0)
        rho0lw(l0-1)=rho0up(l0)
        eta10lw(l0-1)=eta10up(l0)
        eta20lw(l0-1)=eta20up(l0)
        alfa0lw(l0-1)=alfa0up(l0)
      endif
c
c     end of inputs
c     =============
c
      return
      end
