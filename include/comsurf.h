c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comsurf.h                                                    
c *** Contents: Common declarations for surface dependent variables

      integer noc,nse,nld,ntyps,iclimflux
      parameter (noc=1,nse=2,nld=3,ntyps=3)
      real*8  tsurf(nlat,nlon),tsurfn(nlat,nlon,ntyps)
      real*8  fractoc(nlat,nlon),fractn(nlat,nlon,ntyps)
      real*8  albes(nlat,nlon),albesn(nlat,nlon,ntyps)
      real*8  albesnR(nlat,nlon)
      real*8  alb2esn(nlat,nlon,ntyps),alb2es(nlat,nlon) 
      real*8  emisn(ntyps),albocef
      real*8  hesws(nlat,nlon),heswsn(nlat,nlon,ntyps)
      real*8  hesw0(nlat,nlon),hesw0n(nlat,nlon,ntyps)
      real*8  hesw1(nlat,nlon),hesw1n(nlat,nlon,ntyps)
      real*8  hesw2(nlat,nlon),hesw2n(nlat,nlon,ntyps)
      real*8  ulrads(nlat,nlon),ulradsn(nlat,nlon,ntyps)
      real*8  ulrad0(nlat,nlon),ulrad0n(nlat,nlon,ntyps)
      real*8  ulrad1(nlat,nlon),ulrad1n(nlat,nlon,ntyps)
      real*8  ulrad2(nlat,nlon),ulrad2n(nlat,nlon,ntyps)
      real*8  dlrads(nlat,nlon),dlradsn(nlat,nlon,ntyps)
      real*8  evap(nlat,nlon),evapn(nlat,nlon,ntyps)
      real*8  eflux(nlat,nlon),efluxn(nlat,nlon,ntyps)
      real*8  hflux(nlat,nlon),hfluxn(nlat,nlon,ntyps)
      real*8  hficof(nlat,nlon)
      real*8  cdragv(nlat,nlon),cdragvn(nlat,nlon,ntyps)
      real*8  q10(nlat,nlon),q10n(nlat,nlon,ntyps)
      real*8  qsurf(nlat,nlon),qsurfn(nlat,nlon,ntyps)
      real*8  tempsg(nlat,nlon),tempsgn(nlat,nlon,ntyps)
      real*8  pground(nlat,nlon),pgroundn(nlat,nlon,ntyps)
      real*8  rmountn(nlat,nlon,ntyps),evfacan(nlat,nlon,ntyps)
      real*8  lwrmois(nlat,nlon),epss
      real*8  winstu(nlat,nlon),winstv(nlat,nlon)
      real*8  abmoism(nlat,nlon),adsnow(nlat,nlon),abmoisg(nlat,nlon)
      real*8  arunofl(nlat,nlon),arunofo(nlat,nlon),ahic(nlat,nlon)
      
      common /ec_csurf/  tsurf,tsurfn,albes,albesn,fractoc,fractn,epss,
     &                   alb2esn,alb2es, 
     &                hesws,heswsn,hesw0,hesw0n,hesw1,hesw1n,
     &                hesw2,hesw2n,ulrads,ulradsn,ulrad0,ulrad0n,
     &                ulrad1,ulrad1n,ulrad2,ulrad2n,dlrads,dlradsn,
     &                evap,evapn,eflux,efluxn,hflux,hfluxn,lwrmois,
     &                winstu,winstv,cdragv,cdragvn,q10,q10n,
     &                qsurf,qsurfn,rmountn,tempsg,tempsgn,
     &                pground,pgroundn,evfacan,abmoisg,abmoism,adsnow,
     &                arunofo,arunofl,ahic,hficof,albesnR
      common /ec_isurf/  iclimflux
      common /ec_emis/ emisn, albocef
     
      
