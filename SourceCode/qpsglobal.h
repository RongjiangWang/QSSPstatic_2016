c     CONSTANTS
c     =========
      double precision PI
      parameter(PI=3.14159265358979d0)
      double precision DEG2RAD,KM2M,DAY2SEC
      parameter(DEG2RAD=1.745329251994328d-02)
	  parameter(KM2M=1.0d+03,DAY2SEC=8.64d+04)
      double precision BIGG
      parameter(BIGG=6.6732d-11)
      double precision RESOLUT
      parameter(RESOLUT=0.01d0)
      double precision FLTAPER
      parameter(FLTAPER=0.2d0)
	  double precision THICKMAX
      parameter(THICKMAX=5.0d+05)
	  double precision REXMIN
      parameter(REXMIN=1.0d-06)
c
c     GLOBAL INDEX PARAMETERS FOR DEFINING ARRAYS
c     ===========================================
c     lymax: max. number of layers
c     nrmax: max. number of receivers
c     nsmax: max. number of sub-events
c     ngrnmax: max. number of source depths (Green's function data bases)
c     nfmax: max. number of frequency samples
c     ldegmax: max. degree of Legendre polynomials
c
      integer lymax,nrmax,nsmax,ngrnmax,nfmax,ldegmax
	  integer ndmax,nsnapmax
      parameter(nrmax=301,nsmax=500,ngrnmax=50)
      parameter(nfmax=16384,ldegmax=40000)
      parameter(lymax=100+ngrnmax)
      parameter(ndmax=5,nsnapmax=50)
c
c     GREEN FUNCTION PARAMETERS
c     =========================
c
c     nt = power of 2 integer nearest to ntcut
c     nf = nt/2
c     lygrn = layer number of Green's function source
c     grnsel = selection of Green's function to be calculated
c     dt = time sampling interval
c     df = frequency sampling interval
c     fi = imaginary frequency derived from anti-aliasing factor
c     ldeggr = critical harmonic degree
c     yr,yt,yp,ys,yg,yv,yw = Green's function spectra
c     (r,t,p = displacement, s = sound, g = gravity, v,w = tilt)
c     grnfile = file name of Green's function spectra
c
	  logical nogravity
	  logical seldis,selvst,seltil,selgra,selgeo
      integer istype,ngrn,nt,ntcut,nf,nsnap
	  integer lyadd,ldeggr,ldegcut
      integer lygrn(ngrnmax),grnsel(ngrnmax)
      integer ldegpsv(lymax),ldegsh(lymax)
      integer nruku(0:ldegmax,lymax)
      double precision dt,df,fi,taumin,rexfault
	  double precision rearth,rr0,grsurf,epilat,epilon
      double precision grndep(ngrnmax)
      complex yuz(0:ldegmax,4,0:ndmax),yur(0:ldegmax,4,0:ndmax)
      complex yut(0:ldegmax,4,0:ndmax),ypp(0:ldegmax,4,0:ndmax)
      complex ygr(0:ldegmax,4,0:ndmax),ygd(0:ldegmax,4,0:ndmax)
      complex ytr(0:ldegmax,4,0:ndmax),ytt(0:ldegmax,4,0:ndmax)
      character*80 modelmod,grnfile(ngrnmax),
     &             uzgrnfile(ngrnmax),urgrnfile(ngrnmax),
     &             utgrnfile(ngrnmax),ppgrnfile(ngrnmax),
     &             grgrnfile(ngrnmax),gdgrnfile(ngrnmax),
     &             trgrnfile(ngrnmax),ttgrnfile(ngrnmax)
      common /lgreen/ nogravity,seldis,selvst,seltil,selgra,selgeo
      common /igreen/ istype,ngrn,nt,ntcut,nf,nsnap,lyadd,ldeggr,
     &                ldegcut,lygrn,grnsel,ldegpsv,ldegsh,nruku
      common /dgreen/ dt,df,fi,taumin,rexfault,rearth,rr0,
     &                grsurf,epilat,epilon,grndep
      common /rgreen/ yuz,yur,yut,ypp,ygr,ygd,ytr,ytt
      common /cgreen/ modelmod,uzgrnfile,urgrnfile,utgrnfile,ppgrnfile,
     &                grgrnfile,gdgrnfile,trgrnfile,ttgrnfile,grnfile
c
c     RECEIVER PARAMETERS
c     ===================
c
c     nr = number of receivers
c     ioutform = output format:
c                1 = vertical(z)/north(n)/east(e)
c                2 = vertical(z)/radial(r)/transversal(t)
c     latr,lonr = geographic receiver coordinates
c     tred = time reduction
c     rname = receiver name
c     dpr = uniform receiver depth
c     uxout, ... = synthetic seismogram outputs
c          (u = displacement,v = velocity,a = acceleration,gr = gravimetric)
c
      integer nr,ioutform
      double precision dpr,freeairgrd
      double precision latr(nrmax),lonr(nrmax),tsnap(nsnapmax)
      double complex uz(-1:nfmax,nrmax),ur(-1:nfmax,nrmax)
      double complex ut(-1:nfmax,nrmax),upp(-1:nfmax,nrmax)
      double complex ugr(-1:nfmax,nrmax),ugd(-1:nfmax,nrmax)
	  double complex utr(-1:nfmax,nrmax),utt(-1:nfmax,nrmax)
      character*10 rname(nrmax)
      character*80 urout,utout,uzout,
     &             ppout,grout,gdout,trout,
     &             ttout,snapshot(nsnapmax)
      common /ireceiver/ nr,ioutform
      common /dreceiver/ dpr,freeairgrd,latr,lonr,
     &                   tsnap,uz,ur,ut,upp,ugr,ugd,utr,utt
      common /creceiver/ rname,urout,utout,uzout,
     &                   ppout,grout,gdout,trout,ttout,snapshot
c
c     SUBEVENT PARAMETERS
c     ===================
c
c     ns = number of sub-events
c     isg1,isg2 = the respective Green's function number of subevents
c     mtt,... = moment tensor
c     lats,lons = geographic coordinates
c
      integer ns
      integer isg1(ngrnmax),isg2(ngrnmax),nsg(ngrnmax)
      double precision mtt(nsmax),mpp(nsmax),mrr(nsmax),
     &                 mtp(nsmax),mpr(nsmax),mrt(nsmax),
     &                 lats(nsmax),lons(nsmax),deps(nsmax)
      common /isubevents/ns,isg1,isg2,nsg
      common /dsubevents/mtt,mpp,mrr,mtp,mpr,mrt,
     &                   lats,lons,deps
c
c     SOURCE-RECEIVER CONFIGURATION PARAMETERS
c     ========================================
c
      integer idr(nsmax,nrmax)
      double precision dis(nsmax,nrmax)
      double complex ssa(nsmax,nrmax),ss2a(nsmax,nrmax)
      double complex csa(nsmax,nrmax),cs2a(nsmax,nrmax)
      double complex ssb(nsmax,nrmax),csb(nsmax,nrmax)
      double complex ssd(nsmax,nrmax),ssf(nsmax,nrmax)
      common /isourcerec/ idr
      common /dsourcerec/ dis,ssa,ss2a,csa,cs2a,ssb,csb,ssd,ssf
c
c     LEGENDRE POLYNOMIALS TABLES
c     ===========================
c
c     plm = associated Legendre polynomials divided by sin(x)**m
c
      double precision plm(0:ldegmax,0:2)
      common /dlegendre/ plm
c
c     ORIGINAL MODEL PARAMETERS
c     =========================
c
c     lys = layer number of source
c     lyr = layer number of receiver
c     lycm = layer number of core-mantle boundary
c     lycc = layer number of inner and outer core boundary
c     ly0 = max. layer number for integration
c     lyuppsv = upper starting layer number for p-sv solution
c     lyupt = upper starting layer number for sh solution
c     lylwpsv = lower starting layer number for p-sv solution
c     lylwt = lower starting layer number for sh solution
c
      integer lys,lyr,lycm,lycc,ly0
      integer lylwpsv(0:ldegmax),lylwt(0:ldegmax)
	  integer lyuppsv(0:ldegmax),lyupt(0:ldegmax)
      common /ilnumber/ lys,lyr,lycm,lycc,ly0,
     &                  lyuppsv,lyupt,lylwpsv,lylwt
c
c     dp = depth (up = top of layer, lw = bottom of layer)
c     vp, vs = p and s wave velocity
c     ro = density
c     eta,alfa = viscoelastic parameters
c
      integer l0
      double precision dp0(2*lymax),dp0up(lymax),dp0lw(lymax)
      double precision kap0(2*lymax),kap0up(lymax),kap0lw(lymax)
      double precision mue0(2*lymax),mue0up(lymax),mue0lw(lymax)
      double precision rho0(2*lymax),rho0up(lymax),rho0lw(lymax)
      double precision eta10(2*lymax),eta10up(lymax),eta10lw(lymax)
      double precision eta20(2*lymax),eta20up(lymax),eta20lw(lymax)
      double precision alfa0(2*lymax),alfa0up(lymax),alfa0lw(lymax)
      common /imodel0/ l0
      common /dmodel0/ dp0,dp0up,dp0lw,kap0,kap0up,kap0lw,
     &                 mue0,mue0up,mue0lw,rho0,rho0up,rho0lw,
     &                 eta10,eta10up,eta10lw,eta20,eta20up,eta20lw,
     &                 alfa0,alfa0up,alfa0lw
c
c     rr = radius (up = top of layer, lw = bottom of layer)
c     else see above.
c
      double precision rrup(lymax),rrlw(lymax)
      double precision kapup(lymax),kaplw(lymax)
      double precision mueup(lymax),muelw(lymax)
      double precision rhoup(lymax),rholw(lymax)
      double precision grup(lymax),grlw(lymax)
      double precision eta1up(lymax),eta1lw(lymax)
	  double precision eta2up(lymax),eta2lw(lymax)
	  double precision alfaup(lymax),alfalw(lymax)
	  double precision rhoa(lymax),rhob(lymax),beta(lymax)
	  double complex crexup(lymax),crexlw(lymax)
      common /dmodel/ rrup,rrlw,kapup,kaplw,
     &                mueup,muelw,rhoup,rholw,grup,grlw,
     &                eta1up,eta1lw,eta2up,eta2lw,
     &                alfaup,alfalw,rhoa,rhob,beta,
     &                crexup,crexlw
c
c     LAYER MATRICES
c     ==============
c
c     mat2x2 = 2x2 toroidal solution matrix
c     mas3x3 = 3x3 spheroidal solution matrix (l = 0)
c     mas4x4 = 4x4 spheroidal solution matrix (l > 0, in liquid)
c     mas6x6 = 6x6 spheroidal solution matrix (l > 0, in solid)
c     mas(t)inv = inverse solution matrix
c
      double complex mat2x2(2,2,lymax),mat2x2inv(2,2,lymax)
	  double complex mas3x3(3,3,lymax),mas3x3inv(3,3,lymax)
	  double complex mas4x4(4,4,lymax),mas4x4inv(4,4,lymax)
	  double complex mas6x6(6,6,lymax),mas6x6inv(6,6,lymax)
      common /matrix/ mat2x2,mat2x2inv,mas3x3,mas3x3inv,
     &                mas4x4,mas4x4inv,mas6x6,mas6x6inv
c
      double complex cypnorm(6,lymax)
      common /normalization/ cypnorm
