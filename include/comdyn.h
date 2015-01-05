c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comdyn.h                                                    
c *** Contents: Common declarations for dynamical part of atmospheric
c               model of ECbilt
c-----------------------------------------------------------------------
c *** COMMON /ec_intpar/ nshm, ll
c     nshm:   contains numbers 22 down to 1 for index 0 to 21
c     ll:     contains total wavenumber n of each spherical harmonic of
c             the corresponding index
c
c *** COMMON /ec_linopr/ rm, rinhel, diss 
c     rm:     contains zonal wavenumber m of each spherical harmonic of
c             the corresponding index for zonal derivative operator
c     rinhel: Laplace and Helmholtz operator for Q-PSI inversion
c     diss:   dissipation coefficients for each spherical harmonic
c             diss(k,1) : hyperviscosity at the three levels 
c                         (tdif sets timescale)
c             diss(k,2) : Ekman friction at lower level 
c                         (tdis sets timescale)
c
c *** COMMON /ec_logpar/ lgdiss,inf
c     lgdiss: if .true. then orography and land-sea mask dependent 
c             friction at the lower level plus Ekman friction, 
c             else only Ekman friction
c     inf:    if .true. then artificial PV forcing read from file
c
c *** COMMON /ec_metras/
c     pp:     Legendre polynomials defined at Gausian latitudes
c     pd:     mu derivative of Legendre polynomials
c     pw:     weights for Legendre integrals
c
c *** COMMON /ec_phypar/ rdiss, ddisdx, ddisdy
c     rdiss:  landsea-mask/orography dependent friction
c     ddisdx: landsea-mask/orography dependent friction
c     ddisdy: landsea-mask/orography dependent friction
c
c *** COMMON /ec_modpar/ rrdef1,rrdef2,rl1,rl2,relt1,relt2 ,
c                      tdis,trel,tdif,addisl,addish,h0,idif
c     rrdef1: Rossby radius of deformation of 200-500 thickness
c     rrdef2: Rossby radius of deformation of 500-800 thickness
c     rl1:    one over Rossby rad. of def. squared of 200-500 thickness
c     rl2:    one over Rossby rad. of def. squared of 500-800 thickness
c     relt1:  nondimensional relaxation coefficient of 200-500 thickness
c     relt2:  nondimensional relaxation coefficient of 500-800 thickness
c     tdis:   Ekman dissipation timescale in days at lower level
c     trel:   relaxation time scale in days of the temperature
c     tdif:   dissipation timescale of scale-selective diffusion in
c             days for wavenumber nm
c     addisl: parameter used in the computation of the dissipation
c             timescale at the lower level over land
c     addish: parameter used in the computation of the dissipation
c             timescale at the lower level as a function of topography
c     h0:     scale factor for the topographically induced upward motion
c             at the lower level
c     idif:   determines scale-selectivity of hyperviscosity; power of
c             laplace operator
c 
c *** COMMON /ec_sfield/ psi, psit, psito, qprime, dqprdt, for, ws
c     psi:    stream function at the nvl levels
c     psit:   thickness at the ntl levels
c     psito:  thickness at the ntl levels at the previous timestep
c     qprime: potential vorticity
c     dqprdt: time derivative of qprime
c     for:    constant potential vorticity forcing at the nvl levels
c     ws:     only used as portable workspace
c
c *** COMMON /ec_corog/ orog, rmount, rmountn
c     orog:   orography in m. divided by h0
c     rmount: mean orography in m.
c     rmountn orography in m for each surface type
c     dorodl: derivative of orog wrt lambda
c     dorodm: derivative of orag wrt sin(fi)
c
c *** COMMON /ec_zotras/ trigd, trigi, wgg
c             arrays used by the nag version of the fft
c
c *** COMMON /ec_cwind/ uv10,utot,vtot,u800,v800,u500,v500,u200,v200
c     uv10:   wind strength at 10 m, extrapolated from 800 hPa at grid
c             with imposed lower bound of uv10m m/s
c     uvw10:  wind strength at 10 m, extrapolated from 800 hPa at grid
c     uv10r:  reduction factor between 800 hPa and surface winds
c     uv10m:  minimum value of 10 m wind
c     utot :   total zonal wind 
c     vtot :   total meridional wind 
c     u800:   800 hPa geostrophic zonal wind
c     v800:   800 hPa geostrophic meridional wind
c     u500:   500 hPa geostrophic zonal wind
c     v500:   500 hPa geostrophic meridional wind
c     u200:   200 hPa geostrophic zonal wind
c     v200:   200 hPa geostrophic meridional wind
c
c *** COMMON /ec_cgpsi/  grpsi1,grpsi2,grpsi3
c     grpsi1: streamfunction on Gaussian grid at 200 mb
c     grpsi2: streamfunction on Gaussian grid at 500 mb
c     grpsi3: streamfunction on Gaussian grid at 800 mb
c
c *** COMMON /ec_cdwin/  udiv,vdiv,udivg,vdivg
c     udiv:   divergent wind east-west direction in spectral form
c     vdiv:   divergent wind north-south direction in spectral form
c     udivg:  divergent wind east-west direction on Gaussian grid
c     vdivg:  divergent wind north-south direction on Gaussian grid
c
c *** COMMON /ec_cdiv/   divs,divg
c     divs:   divergence in spectral form
c     divg:   divergence on Gaussian grid
c
c *** COMMON /ec_forcx/  iartif,ipvf1,ipvf2,ipvf3,ipvf4,ipvf5
C forcing parameter:                                                   
C iartif:     with (1) or without (0) artificial forcing                
C ipvf1 :     with (1) or without (0) diabatic heating                 
C ipvf2 :     with (1) or without (0) advection of f by divergent wind 
C ipvf3 :     with (1) or without (0) stretching term                  
C ipvf4 :     with (1) or without (0) advection of zeta by divergent   
C             wind                                                     
C ipvf5 :     with (1) or without (0) advection of temperature by      
C                                  divergent wind                      
c *** COMMON /ec_clfor/  forcgg1,forcggw1,forcggs1
c     forcgg1: daily variying artificial PV forcing on Gaussian grid
c              at 200 hPa only calculated from forcggw1 and forcggs1
c     forcggw1: artificial PV forcing winter season on Gaussian grid 
c               at 200 hPa
c     forcggs1: artificial PV forcing summer season on Gaussian grid 
c               at 200 hPa
c
c *** COMMON /ec_cvert/  omegs,omegg
c     omegs: vertical velocity omega in spectral form
c     omegg: vertical velocity omega on Gaussian grid
c
c *** COMMON /ec_cdfor/  dfor1,dfor2
c     dfor1: PV forcing of 200-500 thickness due to diabatical terms
c     dfor2: PV forcing of 500-800 thickness due to diabatical terms
c-----------------------------------------------------------------------


      integer nshm(0:nm),ll(nsh),ipert
      real*8  rm(nsh),rinhel(nsh2,0:5),diss(nsh2,2)
      real*8  pp(nlat,nsh),pd(nlat,nsh),pw(nlat,nsh)
      real*8  rdiss(nlat,nlon),ddisdx(nlat,nlon),ddisdy(nlat,nlon)
      real*8  rrdef1,rrdef2,rl1,rl2,relt1,relt2,tdis,trel,tdif
      real*8  addisl,addish,h0
      integer idif
      real*8  psi(nsh2,nvl),psit(nsh2,ntl),psito(nsh2,ntl)
      real*8  qprime(nsh2,nvl),dqprdt(nsh2,nvl),for(nsh2,nvl),ws(nsh2)
      real*8  orog(nsh2),rmount(nlat,nlon)
      real*8  dorodl(nlat,nlon),dorodm(nlat,nlon)
      real*8  trigd(nlon,2),trigi(nlon,2),wgg(nlat,nlon)
      real*8  utot(nlat,nlon,3),vtot(nlat,nlon,3)
      real*8  u800(nlat,nlon),v800(nlat,nlon),u500(nlat,nlon)
      real*8  v500(nlat,nlon),u200(nlat,nlon),v200(nlat,nlon)
      real*8  udiv(nsh2,nvl),vdiv(nsh2,nvl),udivg(nlat,nlon,nvl),
     *        vdivg(nlat,nlon,nvl),chig(nlat,nlon,nvl),chi(nsh2,nvl)
      real*8  psig(nlat,nlon,nvl),qgpv(nlat,nlon,nvl)
      real*8  geopg(nlat,nlon,nvl)
      real*8  divs(nsh2,nvl),divg(nlat,nlon,nvl)
      integer iartif,ipvf1,ipvf2,ipvf3,ipvf4,ipvf5
      real*8  forcgg1(nlat,nlon),forcggw1(nlat,nlon), 
     *        forcggs1(nlat,nlon)
      real*8  omegs(nsh2,nvl),omegg(nlat,nlon,nvl)
      real*8  dfor1(nsh2),dfor2(nsh2)
      real*8  gekdis(nlat,nlon)

      logical lgdiss,inf

      common /ec_intpar/ nshm, ll, ipert
      common /ec_linopr/ rm, rinhel, diss 
      common /ec_logpar/ lgdiss,inf
      common /ec_metras/ pp, pd, pw
      common /ec_phypar/ rdiss, ddisdx, ddisdy
      common /ec_modpar/ rrdef1,rrdef2,rl1,rl2,relt1,relt2 ,
     *                tdis,trel,tdif,addisl,addish,h0,idif
      common /ec_sfield/ psi, psit, psito, qprime, dqprdt, for, ws
      common /ec_corog/  orog, rmount , dorodl , dorodm
      common /ec_zotras/ trigd, trigi, wgg
      common /ec_cwind/  utot, vtot, u800, v800, u500,v500,u200,v200
      common /ec_cgpsi/  psig,qgpv,geopg
      common /ec_cdwin/  udiv,vdiv,udivg,vdivg
      common /ec_cdiv/   divs,divg,chig,chi
      common /ec_forcx/  iartif,ipvf1,ipvf2,ipvf3,ipvf4,ipvf5
      common /ec_clfor/  forcgg1,forcggw1,forcggs1
      common /ec_cvert/  omegs,omegg
      common /ec_cdfor/  dfor1,dfor2
      common /ec_ekman/  gekdis
