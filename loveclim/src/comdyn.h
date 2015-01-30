!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** File:     comdyn.h
! *** Contents: Common declarations for dynamical part of atmospheric
!               model of ECbilt
!-----------------------------------------------------------------------
! *** COMMON /intpar/ nshm, ll
!     nshm:   contains numbers 22 down to 1 for index 0 to 21
!     ll:     contains total wavenumber n of each spherical harmonic of
!             the corresponding index
!
! *** COMMON /linopr/ rm, rinhel, diss
!     rm:     contains zonal wavenumber m of each spherical harmonic of
!             the corresponding index for zonal derivative operator
!     rinhel: Laplace and Helmholtz operator for Q-PSI inversion
!     diss:   dissipation coefficients for each spherical harmonic
!             diss(k,1) : hyperviscosity at the three levels
!                         (tdif sets timescale)
!             diss(k,2) : Ekman friction at lower level
!                         (tdis sets timescale)
!
! *** COMMON /logpar/ lgdiss,inf
!     lgdiss: if .true. then orography and land-sea mask dependent
!             friction at the lower level plus Ekman friction,
!             else only Ekman friction
!     inf:    if .true. then artificial PV forcing read from file
!
! *** COMMON /metras/
!     pp:     Legendre polynomials defined at Gausian latitudes
!     pd:     mu derivative of Legendre polynomials
!     pw:     weights for Legendre integrals
!
! *** COMMON /phypar/ rdiss, ddisdx, ddisdy
!     rdiss:  landsea-mask/orography dependent friction
!     ddisdx: landsea-mask/orography dependent friction
!     ddisdy: landsea-mask/orography dependent friction
!
! *** COMMON /modpar/ rrdef1,rrdef2,rl1,rl2,relt1,relt2 ,
!                      tdis,trel,tdif,addisl,addish,h0,idif
!     rrdef1: Rossby radius of deformation of 200-500 thickness
!     rrdef2: Rossby radius of deformation of 500-800 thickness
!     rl1:    one over Rossby rad. of def. squared of 200-500 thickness
!     rl2:    one over Rossby rad. of def. squared of 500-800 thickness
!     relt1:  nondimensional relaxation coefficient of 200-500 thickness
!     relt2:  nondimensional relaxation coefficient of 500-800 thickness
!     tdis:   Ekman dissipation timescale in days at lower level
!     trel:   relaxation time scale in days of the temperature
!     tdif:   dissipation timescale of scale-selective diffusion in
!             days for wavenumber nm
!     addisl: parameter used in the computation of the dissipation
!             timescale at the lower level over land
!     addish: parameter used in the computation of the dissipation
!             timescale at the lower level as a function of topography
!     h0:     scale factor for the topographically induced upward motion
!             at the lower level
!     idif:   determines scale-selectivity of hyperviscosity; power of
!             laplace operator
!
! *** COMMON /sfield/ psi, psit, psito, qprime, dqprdt, for, ws
!     psi:    stream function at the nvl levels
!     psit:   thickness at the ntl levels
!     psito:  thickness at the ntl levels at the previous timestep
!     qprime: potential vorticity
!     dqprdt: time derivative of qprime
!     for:    constant potential vorticity forcing at the nvl levels
!     ws:     only used as portable workspace
!
! *** COMMON /corog/ orog, rmount, rmountn
!     orog:   orography in m. divided by h0
!     rmount: mean orography in m.
!     rmountn orography in m for each surface type
!     dorodl: derivative of orog wrt lambda
!     dorodm: derivative of orag wrt sin(fi)
!
! *** COMMON /zotras/ trigd, trigi, wgg
!             arrays used by the nag version of the fft
!
! *** COMMON /cwind/ uv10,utot,vtot,u800,v800,u500,v500,u200,v200
!     uv10:   wind strength at 10 m, extrapolated from 800 hPa at grid
!             with imposed lower bound of uv10m m/s
!     uvw10:  wind strength at 10 m, extrapolated from 800 hPa at grid
!     uv10r:  reduction factor between 800 hPa and surface winds
!     uv10m:  minimum value of 10 m wind
!     utot :   total zonal wind
!     vtot :   total meridional wind
!     u800:   800 hPa geostrophic zonal wind
!     v800:   800 hPa geostrophic meridional wind
!     u500:   500 hPa geostrophic zonal wind
!     v500:   500 hPa geostrophic meridional wind
!     u200:   200 hPa geostrophic zonal wind
!     v200:   200 hPa geostrophic meridional wind
!
! *** COMMON /cgpsi/  grpsi1,grpsi2,grpsi3
!     grpsi1: streamfunction on Gaussian grid at 200 mb
!     grpsi2: streamfunction on Gaussian grid at 500 mb
!     grpsi3: streamfunction on Gaussian grid at 800 mb
!
! *** COMMON /cdwin/  udiv,vdiv,udivg,vdivg
!     udiv:   divergent wind east-west direction in spectral form
!     vdiv:   divergent wind north-south direction in spectral form
!     udivg:  divergent wind east-west direction on Gaussian grid
!     vdivg:  divergent wind north-south direction on Gaussian grid
!
! *** COMMON /cdiv/   divs,divg
!     divs:   divergence in spectral form
!     divg:   divergence on Gaussian grid
!
! *** COMMON /forcx/  iartif,ipvf1,ipvf2,ipvf3,ipvf4,ipvf5
! forcing parameter:
! iartif:     with (1) or without (0) artificial forcing
! ipvf1 :     with (1) or without (0) diabatic heating
! ipvf2 :     with (1) or without (0) advection of f by divergent wind
! ipvf3 :     with (1) or without (0) stretching term
! ipvf4 :     with (1) or without (0) advection of zeta by divergent
!             wind
! ipvf5 :     with (1) or without (0) advection of temperature by
!                                  divergent wind
! *** COMMON /clfor/  forcgg1,forcggw1,forcggs1
!     forcgg1: daily variying artificial PV forcing on Gaussian grid
!              at 200 hPa only calculated from forcggw1 and forcggs1
!     forcggw1: artificial PV forcing winter season on Gaussian grid
!               at 200 hPa
!     forcggs1: artificial PV forcing summer season on Gaussian grid
!               at 200 hPa
!
! *** COMMON /cvert/  omegs,omegg
!     omegs: vertical velocity omega in spectral form
!     omegg: vertical velocity omega on Gaussian grid
!
! *** COMMON /cdfor/  dfor1,dfor2
!     dfor1: PV forcing of 200-500 thickness due to diabatical terms
!     dfor2: PV forcing of 500-800 thickness due to diabatical terms
!-----------------------------------------------------------------------


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
      real*8  udiv(nsh2,nvl),vdiv(nsh2,nvl),udivg(nlat,nlon,nvl), &
     &        vdivg(nlat,nlon,nvl),chig(nlat,nlon,nvl),chi(nsh2,nvl)
      real*8  psig(nlat,nlon,nvl),qgpv(nlat,nlon,nvl)
      real*8  geopg(nlat,nlon,nvl)
      real*8  divs(nsh2,nvl),divg(nlat,nlon,nvl)
      integer iartif,ipvf1,ipvf2,ipvf3,ipvf4,ipvf5
      real*8  forcgg1(nlat,nlon),forcggw1(nlat,nlon),forcggs1(nlat,nlon)
      real*8  omegs(nsh2,nvl),omegg(nlat,nlon,nvl)
      real*8  dfor1(nsh2),dfor2(nsh2)
      real*8  gekdis(nlat,nlon)

      logical lgdiss,inf

      common /intpar/ nshm, ll, ipert
      common /linopr/ rm, rinhel, diss
      common /logpar/ lgdiss,inf
      common /metras/ pp, pd, pw
      common /phypar/ rdiss, ddisdx, ddisdy
      common /modpar/ rrdef1,rrdef2,rl1,rl2,relt1,relt2, &
     &                tdis,trel,tdif,addisl,addish,h0,idif
      common /sfield/ psi, psit, psito, qprime, dqprdt, for, ws
      common /corog/  orog, rmount , dorodl , dorodm
      common /zotras/ trigd, trigi, wgg
      common /cwind/  utot, vtot, u800, v800, u500,v500,u200,v200
      common /cgpsi/  psig,qgpv,geopg
      common /cdwin/  udiv,vdiv,udivg,vdivg
      common /cdiv/   divs,divg,chig,chi
      common /forcx/  iartif,ipvf1,ipvf2,ipvf3,ipvf4,ipvf5
      common /clfor/  forcgg1,forcggw1,forcggs1
      common /cvert/  omegs,omegg
      common /cdfor/  dfor1,dfor2
      common /ekman/  gekdis
