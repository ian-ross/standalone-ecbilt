












c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  bloc "ice.com" : 'include' in the routines linked with the ice
c  modif : 25/09/98
c  modif : 18/02/05 AnneM
 

c tfsn      Melting point temperature of the snow
c tfsg      Melting point temperature of the ice
c xkn       Conductivity of the snow
c xkg       Conductivity of the ice
c rcpn      Density times specific heat for the snow
c rcpg      Density times specific heat for the ice
c xlg       Latent heat of fusion of the ice
c xln       Latent heat of fusion of the snow
c rhog      Density of the ice
c rhon      Density of the snow
c emig      Emissivity of the ice
c sglace    Salinity of the ice
c hmelt     Maximum melting at the bottom
c acrit(2)  Minimum fraction for leads
c hgcrit(2) Ice thickness for lateral accretion
c hgmin     Ice thickness corr. to max. energy stored in brine pocket
c hndif     Computation of temp. in snow or not
c hgdif     Computation of temp. in ice or not
c hglim     Minimum ice thickness
c amax      Maximum lead fraction
c uscomi    =1.0/(1.0-amax)
c beta      Numerical caracteritic of the scheme for diffusion in ice
c ddtb      Time step for ice thermodynamics (s)
c swiqst    Energy stored in brine pocket or not
c parlat    Percentage of energy used for lateral ablation
c hakspl    Slope of distr. for Hakkinen-Mellor's lateral melting
c hibspl    Slope of distribution for Hibler's lateral melting
c exld      Exponent for leads-closure rate
c hakdif    Coefficient for diffusions of ice and snow
c hth       Threshold thickness for comp. of eq. thermal conductivity
c hnzst     Thickness of the surf. layer in temp. computation
c parsub    Switch for snow sublimation or not
c alphs     Used to take into account settling during snow-ice formation
c cnscg     ratio  rcpn/rcpg
c nbits     Number of time steps in Newton -Raphson procedure
c stefan    Stefan-Boltzman constant
c xsn       Sublimation heat for the snow
c too       temperature of the triple point for water
c vkarmn    von Karman constant
c cevap     Latent heat of evaporation of water
c zemise    Emissivity of water
c rhoesn    1/rhon
c firg      IR flux over the ice (only used for outputs)
c fcsg      Sensible heat flux over the ice (only used for outputs)
c fleg      Latent heat flux over the ice (only used for outputs)
c dvosbq    Variation of volume at surface (only used for outputs)
c dvobbq    Variation of ice volume at the bottom ice (only used for outputs)
c dvolbq    Total variation of ice volume (only used for outputs)
c dvonbq    Lateral Variation of ice volume (only used for outputs)
c ts        Surface temperature of the ice
c tfu       Melting point temperature of sea water
c hnbq      Snow thickness
c hgbq      Ice thickness
c hgbqp     Ice production/melting
c albq      Leads fraction
c qstobq    Energy stored in the brine pockets
c fbbq      Heat flux at the ice base
c tbq       Temperature inside the ice/snow layer
c dmnbq     Variation of snow mass
c dmgbq     Variation of ice mass
c qlbq      heat balance of the lead (or of the open ocean)
c qcmbq     Energy needed to bring the ocean surface layer until its freezing
c           point (at a factor 2)
c thcm      part of the solar energy used in the lead heat budget
c fstrbq    Solar flux transmitted trough the ice
c ffltbq    Array linked with the max heat contained in brine pockets (?)
c fscmbq    Linked with the solar flux below the ice (?)
c fsbbq     Also linked with the solar flux below the ice (?)
c qfvbq     Array used to store energy in case of toral lateral ablation (?)
c xzo       rugosity of the ice (no more used)
c dmgwi     Variation of the mass of snow ice
c psbq      Surface air pressure
c tabq      Surface air temperature
c qabq      Surface air humidity
c vabq      Surface wind velocity
c hnplbq    Snow precipitation
c fevabq    Evaporation flux
c fsolcn    Solar flux at the ocean surface
c fsolg     Solar flux at the ice surface
c flecn     Latent heat flux at the ocean surface
c fcscn     Sensible heat flux at the ocean surface
c tenagx    Wind stress at the ice surface (x)
c tenagy    Wind stress at the ice surface (y)
c albege    Albedo of the snow or ice (only for outputs)
c tairox    Wind stress at the ocean surface (x)
c tairoy    Wind stress at the ocean surface (y)
c ratbqg    Longwave downward radiation flux over the ice
c ratbqo    Longwave downward radiation flux over the ocean
c cloud     Cloud fraction
c tdew      Air relative humidity
c albecn    Albedo of the ocean (only for outputs)
c tauc      Cloud optical depth
c runoff    river runoff
c sdvt      u*^2/(Stress/density)
c fcm1      Solar flux at the ocean surface or the ocean-ice interface
c fcm2      Non-solar flux at the ocean surface
c fwat      Freshwater flux (change of definition between the routines)
c reslum    Relative absorption of solar radiation in each ocean level
c alct      lead opening/closure rate due to thermodynamics 
c alcd      lead opening/closure rate due to convergence/divergence 
c alcr      lead opening/closure rate due to shear motion if Clds is activated
c hgcol     Collection thickness of sea ice in leads when it is variable (Cvhg enabled)
c ficebergn Heat due to iceberg melting in the Northern Hemisphere
c ficebergs Heat due to iceberg melting in the Southern Hemisphere
c iiceber, jiceber index for location of iceberg melting
c areiceb   surface for iceberg melting
c tmics     Mask for ice shelve melting
c toticesm  total melting due to ice shelve melting
 
c--common blocs :
c
      common / cstbq /
     &  tfsn,tfsg,xkn,xkg,rcpn,rcpg,xlg,xln,rhog,rhon,
     &  emig,sglace,hmelt,acrit(2),hgcrit(2),hgmin,hndif,
     &  hgdif,hglim,amax,uscomi,beta,ddtb,swiqst,parlat,hakspl,
     &  hibspl,exld,hakdif,hth,hnzst,parsub,alphs,cnscg,nbits
c
      common / fluxsf /
     &  stefan,xsn,too,vkarmn,cevap,zemise,rhoesn
c
      common / comdia /
     &  firg(imax,jmax),fcsg(imax,jmax),fleg(imax,jmax),
     &  dvosbq(imax,jmax),dvobbq(imax,jmax),dvolbq(imax,jmax),
     &  dvonbq(imax,jmax)
c
      common / comban /
     &   ts(imax,jmax),tfu(imax,jmax),hnbq(imax,jmax)
     &  ,hgbq(imax,jmax),hgbqp(imax,jmax)
     &  ,qstobq(imax,jmax),fbbq(imax,jmax),tbq(imax,jmax,nkb0)
     &  ,dmnbq(imax,jmax),dmgbq(imax,jmax)
     &  ,qlbq(imax,jmax),qcmbq(imax,jmax),thcm(imax,jmax)
     &  ,fstrbq(imax,jmax),ffltbq(imax,jmax),fscmbq(imax,jmax)
     &  ,fsbbq(imax,jmax),qfvbq(imax,jmax),xzo(imax,jmax)
     &  ,dmgwi(imax,jmax)
     &  ,alct(imax,jmax),alcd(imax,jmax)
     &  ,alcr(imax,jmax)
     &  ,hgcol(imax,jmax)

      common / comban2 /
     &  albq(imax,jmax)

c-uwind et vwind en common pour les icebergs
      common / comfor /
     &  psbq(imax,jmax),tabq(imax,jmax),
     &  qabq(imax,jmax),vabq(imax,jmax),
     &  hnplbq(imax,jmax),fevabq(imax,jmax),fsolcn(imax,jmax),
     &  fsolg(imax,jmax),flecn(imax,jmax),fcscn(imax,jmax),
     &  tenagx(imax,jmax),tenagy(imax,jmax),
     &  albege(imax,jmax),tairox(imax,jmax),tairoy(imax,jmax),
     &  ratbqg(imax,jmax),ratbqo(imax,jmax),cloud(imax,jmax),
     &  tdew(imax,jmax),uwind(imax,jmax),vwind(imax,jmax),
     &  albecn(imax,jmax),tauc(imax,jmax),runoff(imax,jmax),
     &  sdvt(imax,jmax)

      common / comca /
     &  fcm2(imax,jmax),
     &  fwat(imax,jmax),
     &  reslum(imax,jmax,0:kmax+1)
 
      common / comca2 /
     &  fcm1(imax,jmax)

      common / comiceb/
     &  ficebergn,ficebergs,areicebn,areicebs,
     &  tmics(imax,jmax),toticesm,
     &  iicebern1,iicebern2,jicebern1,jicebern2,
     &  iicebers1,iicebers2,jicebers1,jicebers2

      real*4 vwx(57,kmax+1,0:nbsmax)
      common / comstr/
     &  vwx

 
Cage  common / sorage/ agen(imax,jmax),ageg(imax,jmax)
Ccp2  common / fderice / fder(imax,jmax)
c
c--fin du fichier "ice.com"
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
