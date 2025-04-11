	subroutine qpssublayer(ierr)
	implicit none
      integer ierr
c
	include 'qpsglobal.h'
c
c	work space
c
	integer i,j,i0,l,ly,lya,lyb,lyp,ig
      double precision f,h,dh,z,zz,slw,wvlcut,up,lw,uplw4
	double precision xup,xlw,rrs,rrr,dkap,dmue,rho1,drho,deta1,deta2
      double precision mass,gr0,gr1,gr2,al,bl,betal,dalfa,rhom,kapm
      double precision rr1,rr2,swapup(7),swaplw(7)
      double precision sbessj0,sbessj1,sbessy0,sbessy1
      logical jump,pseudo(lymax)
c
	ly0=0
      zz=0.d0
      gr2=grsurf
	do l=1,l0-1
c
c       modification using Adam-Williamson condition
c
        gr1=gr2
        rr1=rearth-dp0up(l)
        rr2=rearth-dp0lw(l)
        drho=(rho0up(l)-rho0lw(l))/(rr1-rr2)
        rho1=rho0lw(l)-drho*rr2
        mass=PI*(rr1-rr2)
     &      *((4.d0/3.d0)*rho1*(rr1**2+rr1*rr2+rr2**2)
     &        +drho*(rr1**3+rr1**2*rr2+rr1*rr2**2+rr2**3))
        gr2=(gr1*rr1**2-BIGG*mass)/rr2**2
        if(gr2.le.0.d0)then
          print *,' Error in qpssublayer:'
     &          //' density profile inconsistent with surface gravity!'
          stop
        endif
c
        if(eta20up(l).gt.0.d0)then
c
c         1=up, 2=lw
c
          rhom=rho0up(l)
          kapm=kap0up(l)
          call awdens(rhom,rr1,rr2,gr1,gr2,kapm,
     &                al,bl,betal,rho0up(l),rho0lw(l))
          kap0up(l)=4.d0*PI*BIGG*rho0up(l)**2/betal**2
          kap0lw(l)=4.d0*PI*BIGG*rho0lw(l)**2/betal**2
        else 
          al=0.d0
          bl=0.d0
          betal=0.d0
        endif
c
	  h=dp0lw(l)-dp0up(l)
	  dkap=2.d0*dabs(kap0lw(l)-kap0up(l))/(kap0lw(l)+kap0up(l))
        if(mue0lw(l)+mue0up(l).gt.0.d0)then
	    dmue=2.d0*dabs(mue0lw(l)-mue0up(l))/(mue0lw(l)+mue0up(l))
        else
          dmue=0.d0
        endif
        drho=2.d0*dabs(rho0lw(l)-rho0up(l))/(rho0lw(l)+rho0up(l))
        i0=1+idint(dmax1(dkap/RESOLUT,dmue/RESOLUT,drho/RESOLUT))
        dkap=(kap0lw(l)-kap0up(l))/h
	  dmue=(mue0lw(l)-mue0up(l))/h
        drho=(rho0lw(l)-rho0up(l))/h
	  deta1=(eta10lw(l)-eta10up(l))/h
	  deta2=(eta20lw(l)-eta20up(l))/h
        dalfa=(alfa0lw(l)-alfa0up(l))/h
	  dh=h/dble(i0)
	  do i=1,i0
	    ly0=ly0+1
	    if(ly0.ge.lymax-lyadd)then
	      stop ' Error in qpssublayer: lymax too small!'
	    endif
	    z=dble(i-1)*dh
	    rrup(ly0)=rearth-(zz+z)
          if(i.eq.1)then
            rhoup(ly0)=rho0up(l)
	      kapup(ly0)=kap0up(l)
	      mueup(ly0)=mue0up(l)
	      eta1up(ly0)=eta10up(l)
	      eta2up(ly0)=eta20up(l)
            alfaup(ly0)=alfa0up(l)
          else
            rhoup(ly0)=rholw(ly0-1)
	      kapup(ly0)=kaplw(ly0-1)
	      mueup(ly0)=muelw(ly0-1)
	      eta1up(ly0)=eta1lw(ly0-1)
	      eta2up(ly0)=eta2lw(ly0-1)
            alfaup(ly0)=alfalw(ly0-1)
          endif
          z=z+dh
	    rrlw(ly0)=rearth-(zz+z)
          if(betal.gt.0.d0)then
            rhoa(ly0)=al
            rhob(ly0)=bl
            beta(ly0)=betal
            rholw(ly0)=rhoa(ly0)*sbessj0(beta(ly0)*rrlw(ly0))
     &                +rhob(ly0)*sbessy0(beta(ly0)*rrlw(ly0))
	      kaplw(ly0)=kap0up(l)*(rholw(ly0)/rho0up(l))**2
          else
            rhoa(ly0)=0.d0
            rhob(ly0)=0.d0
            beta(ly0)=0.d0
            rholw(ly0)=rho0up(l)+drho*z
	      kaplw(ly0)=kap0up(l)+dkap*z
          endif
	    muelw(ly0)=mue0up(l)+dmue*z
	    eta1lw(ly0)=eta10up(l)+deta1*z
	    eta2lw(ly0)=eta20up(l)+deta2*z
          alfalw(ly0)=alfa0up(l)+dalfa*z
c
	  enddo
        zz=zz+h
      enddo
c
c     lowest layer assumed to be an imcompressible sphere
c
      ly0=ly0+1
	if(ly0.ge.lymax-lyadd)then
	  stop ' Error in qpssublayer: lymax too small!'
      endif
      rrup(ly0)=rearth-dp0up(l0)
      rhoup(ly0)=gr2/(4.d0*PI*BIGG*rrup(ly0)/3.d0)
	kapup(ly0)=0.5d0*(kap0up(l0)+kap0lw(l0))
	mueup(ly0)=0.5d0*(mue0up(l0)+mue0lw(l0))
	eta1up(ly0)=0.5d0*(eta10up(l0)+eta10lw(l0))
	eta2up(ly0)=0.5d0*(eta20up(l0)+eta20lw(l0))
      alfaup(ly0)=0.5d0*(alfa0up(l0)+alfa0lw(l0))
c
      rrlw(ly0)=0.d0
      rholw(ly0)=rhoup(ly0)
      kaplw(ly0)=kapup(ly0)
      muelw(ly0)=mueup(ly0)
      eta1lw(ly0)=eta1up(ly0)
      eta2lw(ly0)=eta2up(ly0)
      alfalw(ly0)=alfaup(ly0)
c
      rhoa(ly0)=0.d0
      rhob(ly0)=0.d0
      beta(ly0)=0.d0
c
      do ly=1,ly0
        pseudo(ly)=.false.
      enddo
c
c     add source layers
c
      do ig=1,ngrn
        rrs=rearth-grndep(ig)
        if(rrs.lt.0.d0)then
          stop ' Error in qpssublayer: Wrong source depth!'
        endif
        do ly=1,ly0
          if(rrs.ge.rrlw(ly))then
            lys=ly
            goto 100
          endif
        enddo
100     continue
        if(rrs.lt.rrup(lys).and.rrs.gt.rrlw(lys))then
          do ly=ly0,lys,-1
            rrup(ly+1)=rrup(ly)
	      kapup(ly+1)=kapup(ly)
	      mueup(ly+1)=mueup(ly)
	      rhoup(ly+1)=rhoup(ly)
	      eta1up(ly+1)=eta1up(ly)
	      eta2up(ly+1)=eta2up(ly)
            alfaup(ly+1)=alfaup(ly)
c
            rrlw(ly+1)=rrlw(ly)
	      kaplw(ly+1)=kaplw(ly)
	      muelw(ly+1)=muelw(ly)
	      rholw(ly+1)=rholw(ly)
	      eta1lw(ly+1)=eta1lw(ly)
	      eta2lw(ly+1)=eta2lw(ly)
            alfalw(ly+1)=alfalw(ly)
c
            rhoa(ly+1)=rhoa(ly)
            rhob(ly+1)=rhob(ly)
            beta(ly+1)=beta(ly)
            pseudo(ly+1)=pseudo(ly)
          enddo
          lys=lys+1
          up=(rrs-rrlw(lys))/(rrup(lys-1)-rrlw(lys))
          lw=1.d0-up
          rrlw(lys-1)=rrs
          if(beta(lys-1).gt.0.d0)then
            rholw(lys-1)=rhoa(lys-1)*sbessj0(beta(lys-1)*rrlw(lys-1))
     &                  +rhob(lys-1)*sbessy0(beta(lys-1)*rrlw(lys-1))
	      kaplw(lys-1)=kapup(lys-1)*(rholw(lys-1)/rhoup(lys-1))**2
          else
	      rholw(lys-1)=up*rhoup(lys-1)+lw*rholw(lys)
	      kaplw(lys-1)=up*kapup(lys-1)+lw*kaplw(lys)
          endif
	    muelw(lys-1)=up*mueup(lys-1)+lw*muelw(lys)
	    eta1lw(lys-1)=up*eta1up(lys-1)+lw*eta1lw(lys)
	    eta2lw(lys-1)=up*eta2up(lys-1)+lw*eta2lw(lys)
          alfalw(lys-1)=up*alfaup(lys-1)+lw*alfalw(lys)
c
          rrup(lys)=rrs
	    kapup(lys)=kaplw(lys-1)
	    mueup(lys)=muelw(lys-1)
	    rhoup(lys)=rholw(lys-1)
	    eta1up(lys)=eta1lw(lys-1)
	    eta2up(lys)=eta2lw(lys-1)
          alfaup(lys)=alfalw(lys-1)
c
          rhoa(lys)=rhoa(lys-1)
          rhob(lys)=rhob(lys-1)
          beta(lys)=beta(lys-1)
          pseudo(lys)=.true.
c
          ly0=ly0+1       
        endif
      enddo
c
c     add receiver layer
c
      rrr=rearth-dpr
      if(rrr.lt.0.d0)then
        stop ' Error in qpssublayer: Wrong receiver depth!'
      endif
      do ly=1,ly0
        if(rrr.ge.rrlw(ly))then
          lyr=ly
          goto 200
        endif
      enddo
200   continue
      if(rrr.lt.rrup(lyr).and.rrr.gt.rrlw(lyr))then
        do ly=ly0,lyr,-1
          rrup(ly+1)=rrup(ly)
	    kapup(ly+1)=kapup(ly)
	    mueup(ly+1)=mueup(ly)
	    rhoup(ly+1)=rhoup(ly)
	    eta1up(ly+1)=eta1up(ly)
	    eta2up(ly+1)=eta2up(ly)
          alfaup(ly+1)=alfaup(ly)
c
          rrlw(ly+1)=rrlw(ly)
	    kaplw(ly+1)=kaplw(ly)
	    muelw(ly+1)=muelw(ly)
	    rholw(ly+1)=rholw(ly)
	    eta1lw(ly+1)=eta1lw(ly)
	    eta2lw(ly+1)=eta2lw(ly)
          alfalw(ly+1)=alfalw(ly)
c
          rhoa(ly+1)=rhoa(ly)
          rhob(ly+1)=rhob(ly)
          beta(ly+1)=beta(ly)
          pseudo(ly+1)=pseudo(ly)
        enddo
        lyr=lyr+1
        up=(rrr-rrlw(lyr))/(rrup(lyr-1)-rrlw(lyr))
        lw=1.d0-up
        rrlw(lyr-1)=rrr
        if(beta(lyr-1).gt.0.d0)then
          rholw(lyr-1)=rhoa(lyr-1)*sbessj0(beta(lyr-1)*rrlw(lyr-1))
     &                +rhob(lyr-1)*sbessy0(beta(lyr-1)*rrlw(lyr-1))
	    kaplw(lyr-1)=kapup(lyr-1)*(rholw(lyr-1)/rhoup(lyr-1))**2
        else
	    rholw(lyr-1)=up*rhoup(lyr-1)+lw*rholw(lyr)
	    kaplw(lyr-1)=up*kapup(lyr-1)+lw*kaplw(lyr)
        endif
	  muelw(lyr-1)=up*mueup(lyr-1)+lw*muelw(lyr)
	  eta1lw(lyr-1)=up*eta1up(lyr-1)+lw*eta1lw(lyr)
	  eta2lw(lyr-1)=up*eta2up(lyr-1)+lw*eta2lw(lyr)
        alfalw(lyr-1)=up*alfaup(lyr-1)+lw*alfalw(lyr)
c
        rrup(lyr)=rrr
	  kapup(lyr)=kaplw(lyr-1)
	  mueup(lyr)=muelw(lyr-1)
	  rhoup(lyr)=rholw(lyr-1)
	  eta1up(lyr)=eta1lw(lyr-1)
	  eta2up(lyr)=eta2lw(lyr-1)
        alfaup(lyr)=alfalw(lyr-1)
c
        rhoa(lyr)=rhoa(lyr-1)
        rhob(lyr)=rhob(lyr-1)
        beta(lyr)=beta(lyr-1)
        pseudo(lyr)=.true.
        ly0=ly0+1      
      endif
c
c     determine indices of main interfaces
c
      lycm=ly0+1
      do ly=2,ly0
        if(mueup(ly).le.0.d0.and.muelw(ly-1).gt.0.d0)then
          lycm=ly
          goto 300
        endif
      enddo
300   lycc=ly0+1
      do ly=max0(2,lycm+1),ly0
        if(mueup(ly).gt.0.d0.and.muelw(ly-1).le.0.d0)then
          lycc=ly
          goto 400
        endif
      enddo
400   continue
c
c     determine indices of receiver layer
c
      rrr=rearth-dpr
      do ly=1,ly0
        if(rrr.ge.rrup(ly))then
          lyr=ly
          goto 500
        endif
      enddo
500   continue
c
c     determine indices of source layers
c
      do ig=1,ngrn
        rrs=rearth-grndep(ig)
        lygrn(ig)=ly0+1
        do ly=1,ly0
          if(rrs.ge.rrup(ly))then
            lygrn(ig)=ly
            if(lygrn(ig).ge.lycm)then
              stop ' Error in qpssublayer: source below the mantle!'
            endif
            goto 600
          endif
        enddo
600     continue
      enddo
c
      grup(1)=grsurf
      do ly=1,ly0-1
        if(ly.gt.1)then
          grup(ly)=grlw(ly-1)
        endif
        if(beta(ly).le.0.d0)then
          drho=(rhoup(ly)-rholw(ly))/(rrup(ly)-rrlw(ly))
          rho1=rholw(ly)-drho*rrlw(ly) 
          mass=PI*(rrup(ly)-rrlw(ly))*((4.d0/3.d0)*rho1
     &        *(rrup(ly)**2+rrup(ly)*rrlw(ly)+rrlw(ly)**2)
     &        +drho*(rrup(ly)**3+rrup(ly)**2*rrlw(ly)
     &        +rrup(ly)*rrlw(ly)**2+rrlw(ly)**3))
        else
          xup=beta(ly)*rrup(ly)
          xlw=beta(ly)*rrlw(ly)
          mass=4.d0*PI/beta(ly)
     &        *(rhoa(ly)*(rrup(ly)**2*sbessj1(xup)
     &                   -rrlw(ly)**2*sbessj1(xlw))
     &         +rhob(ly)*(rrup(ly)**2*sbessy1(xup)
     &                   -rrlw(ly)**2*sbessy1(xlw)))
        endif
        grlw(ly)=(grup(ly)*rrup(ly)**2-BIGG*mass)/rrlw(ly)**2
      enddo
      if(ly0.gt.1)then
        grup(ly0)=grlw(ly0-1)
      endif
      grlw(ly0)=0.d0
c
      freeairgrd=-grup(lyr)*2.d0/rrup(lyr)
      if(lyr.gt.1)then
        freeairgrd=freeairgrd+4.d0*PI*BIGG*rholw(lyr-1)
      endif
c
      do ly=1,ly0
        crexup(ly)=(1.d0,0.d0)
        crexlw(ly)=(1.d0,0.d0)
      enddo
c
      do ly=1,lycm-1
        cypnorm(1,ly)=(1.d0,0.d0)
        cypnorm(2,ly)=dcmplx(kapup(ly)+mueup(ly)*4.d0/3.d0,0.d0)
        cypnorm(3,ly)=(1.d0,0.d0)
        cypnorm(4,ly)=dcmplx(kapup(ly)+mueup(ly)*4.d0/3.d0,0.d0)
        cypnorm(5,ly)=dcmplx(4.d0*PI*BIGG*Rhoup(ly),0.d0)
        cypnorm(6,ly)=dcmplx(4.d0*PI*BIGG*Rhoup(ly),0.d0)
      enddo
      do ly=lycm,lycc-1
        cypnorm(1,ly)=(1.d0,0.d0)
        cypnorm(2,ly)=(1.d0,0.d0)
        cypnorm(3,ly)=dcmplx(4.d0*PI*BIGG*Rhoup(ly),0.d0)
        cypnorm(4,ly)=dcmplx(4.d0*PI*BIGG*Rhoup(ly),0.d0)
      enddo
      do ly=lycc,ly0
        cypnorm(1,ly)=(1.d0,0.d0)
        cypnorm(2,ly)=dcmplx(kapup(ly)+mueup(ly)*4.d0/3.d0,0.d0)
        cypnorm(3,ly)=(1.d0,0.d0)
        cypnorm(4,ly)=dcmplx(kapup(ly)+mueup(ly)*4.d0/3.d0,0.d0)
        cypnorm(5,ly)=dcmplx(4.d0*PI*BIGG*Rhoup(ly),0.d0)
        cypnorm(6,ly)=dcmplx(4.d0*PI*BIGG*Rhoup(ly),0.d0)
      enddo
c
      open(20,file=modelmod,status='unknown')
	write(*,'(9a)')'  No','       R[km]',' rho[g/cm^3]',
     &    '   kappa[Pa]','     mue[Pa]','  eta1[Pa*s]',
     &    '  eta2[Pa*s]',' alpha[-]',' g[m/s^2]'
      write(20,'(a)')'  no depth[km]  vp[km/s]  vs[km/s] rho[g/cm^3]'
     &             //' eta1[Pa*s] eta2[Pa*s]   alpha[-]'
c
c
      rrlw(ly0)=0.d0
c
      i=0
	do ly=1,ly0
        write(*,1001)ly,rrup(ly)/1.d3,rhoup(ly)/1.d3,
     &               kapup(ly),mueup(ly),eta1up(ly),eta2up(ly),
     &               alfaup(ly),grup(ly)
        j=0
        do ig=1,ngrn
          if(lygrn(ig).eq.ly)then
            j=j+1
            write(*,'(a3,$)')' S '
          endif
        enddo
        if(j.eq.0)write(*,'(a3,$)')'   '
        if(lyr.eq.ly)then
          write(*,'(a3)')' R '
        else
          write(*,'(a3)')'   '
        endif
	  write(*,1002)rrlw(ly)/1.d3,rholw(ly)/1.d3,
     &               kaplw(ly),muelw(ly),eta1lw(ly),eta2lw(ly),
     &               alfalw(ly),grlw(ly)
c
        swapup(1)=(rearth-rrup(ly))/KM2M
        swapup(2)=dsqrt((kapup(ly)+mueup(ly)*4.d0/3.d0)/rhoup(ly))
     &           /KM2M
        swapup(3)=dsqrt(mueup(ly)/rhoup(ly))/KM2M
        swapup(4)=rhoup(ly)/KM2M
        swapup(5)=eta1up(ly)
        swapup(6)=eta2up(ly)
        swapup(7)=alfaup(ly)
        if(ly.eq.1)then
          i=i+1
	    write(20,1000)i,(swapup(j),j=1,7)
        else if(.not.pseudo(ly))then
          jump=.false.
          do j=2,7
            jump=jump.or.swapup(j).ne.swaplw(j)
          enddo
          if(jump)then
            i=i+1
            write(20,1000)i,(swapup(j),j=1,7)
          endif
        endif
        swaplw(1)=(rearth-rrlw(ly))/KM2M
        swaplw(2)=dsqrt((kaplw(ly)+muelw(ly)*4.d0/3.d0)/rholw(ly))
     &           /KM2M
        swaplw(3)=dsqrt(muelw(ly)/rholw(ly))/KM2M
        swaplw(4)=rholw(ly)/KM2M
        swaplw(5)=eta1lw(ly)
        swaplw(6)=eta2lw(ly)
        swaplw(7)=alfalw(ly)
        if(ly.lt.ly0)then
          if(.not.pseudo(ly+1))then
            i=i+1
            write(20,1000)i,(swaplw(j),j=1,7)
          endif
        endif
      enddo
      close(20)
c
1000  format(i4,3f10.3,f12.3,2E11.2,f11.3)
1001	format(i4,2f12.4,4E12.4,2f9.4,$)
1002  format(f16.4,f12.4,4E12.4,2f9.4)
      return
	end
