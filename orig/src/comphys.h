c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comphys.h
c *** Contents: Common declarations for physical part of atmospheric
c ***           model of ECbilt
c      COMMON /cgamma/ gamgr,gamean,gamvar,gamsq
c      solarc:      solar constant.
c      common /ctemag/ tempsg,temp2g,temp4g
c      tempsg:      temperature at 10 meters height.
c      temp2g:      temperature at 350 mb.
c      temp4g:      temperature at 650 mb.
c
c      COMMON /ctempm/ tempm,thform,temp2gm,temp4gm,
c                     tempsgm,tsurfm
c      tempm:       global mean temperature.
c      thform:      global mean
c      temp2gm:     global mean temperature at 350 mb.
c      temp4gm:     global mean temperature at 650 mb.
c      tempsgm:     global mean temperature at 10 meters height.
c      tsurfm:      global mean temperature at surface.
c
c      COMMON /cpot/   potm,pot2g,pot4g
c      potm:        global mean potential temperature.
c      pot2g:       potential temperature at 350 mb.
c      pot4g:       potential temperature at 650 mb.
c
c      COMMON /csurf/  tsurf,lsmask,sst,sstday
c      tsurf:       temperature at surface.
c      lsmask:      land sea mask.
c      sst:         sea surface tmeperature.???????????????
c      sstday:      sea surface temperature each day.?????
c
c      COMMON /sunfr/  solarc,Q0
c      solarc;      sun constant.
c
c      COMMON /cswrad/ hesw1,hesw2,hesws,albes,albea,albeaw,albeas,
c      *                abso1,abso2
c      hesw1:       solar radiation heating rate.
c      hesw2:       solar radiation heating rate.
c      hesws:       solar radiation heating rate at the surface.
c      albes:
c      albea:       solar radiation reflective coefficient.
c      albeaw:
c      albeas:
c      abso1:       solar radiation absorbtion coefficient.
c      abso2:       solar radiation absorbtion coefficient.
c
c      COMMON /cpar1/  rowat,roair,cdrag,cpair,cvair,rlatvap,rlatcon,
c      *                sboltz,rlatsub,rlatfus,cwater,gamma,rkappa,ps,
c      *                aa,bb,gamd,tzero,alphad,alphav,alphas,cwdrag
c      rowat:       density of water.
c      roair:       mean air desity at sea level.
c      cdrag:       coefficient in sensible and latent air-sea heat flux.
c      cwdrag:      coefficient in wind stress.
c      cdragv:      coef. in sen. en lat. heat flux depending on roughness
c                   length and stability (Richardson number)
c      richar:      richardson number
c      cpair:       specific heat of dry air at constant pressure.
c      cvair:       specific heat of dry air at constant volum.
c      rlatvap:     latent heat uptake due to evaporation.
c      rlatcon:     latent heat release due to condensition.
c      sboltz:      stefan-boltzmann constant.
c      rlatsub:     latent heat of sublimation.
c      rlatfus:     latent heat of fusion.
c      cwater:      4180
c      gamma:       cpair/cvair
c      rkappa:      =(gamma-1)/gamma
c      ps:          surface pressure.
c      aa:          =(350/1000)**rkappa.
c      bb:          =(650/1000)**rkappa
c      gamad:
c      tzero:       =273.15
c      alphad:      =roair*cdrag*cpair
c      alphav:      =roair*cdrag*rlatvap
c      alphas:      =roair*cdrag*rlatsub
c      albes:       surface albedo.

c      COMMON /clwrad/ ulrad1,ulrad2,ulrads,dlrads
c      ulrad2:      upwards long wave radiation at ?.
c      ulrad1:      upwards long wave radiation at ?.
c      ulrads:      upwards long wave radiation at the surface.
c      dlrads:      downwards long wave radiation at the surface.
c
c      COMMON /clwpar/ pteta,pa,pb,pc,clfrac
c      pteta:
c      pa:
c      pb:
c      pc:
c      clfrac:
c
c      COMMON /csflux/  Hflux(nlat,nlon),Eflux(nlat,nlon),dEflux(nlat,nlon)
c      Hflux:       sensible heat flux at the surface.
c      Eflux:       latent heat flux at the surfac.
c      dEflux:      latent heat flux at the surfac?????????????????.
c
c      COMMON /crain/  corain,dyrain,torain
c      corain;      convective rain.
c      dyrain:      dynamic rain.
c      torain:      total rain.
c
c      COMMON /cmoist/ qsats,qsat4,relhum,qg
c      qsats:
c      qsat4
c      relhum:      relative humidity.
c      qg
c
c      COMMON /moisturew/ relhmax, ihavm, ivavm, imsink
c      relhmax:     maximum relative humidity before saturation.
c      ihavm:       with (1) or without (0) horizontal divergence of moisture.
c      ivavm:       with (1) or without (0) vertical divergence of moisture.
c      imsink:      with (1) or without (0) the source or sink of moisture.
c
c      COMMON /cmois/  rmoiss,rmoisg,precip,evap,dcmoisg
c      rmoisg:      specific humidity.
c      precip:      (not found)
c      evap:        evaporation.
c      dcmoisg:
c      COMMON /conv/   imoisr,imoisc,tetacr,gams,teta,convn
c      imoisr
c      imoisc
c      tetacr
c      gams
c      teta:        0.5*(pot2g - pot4g)
c      convn
c
c      COMMON /cthfor/ thforg1,thforg2,thfor1,thfor2,vhforg1,vhforg2
c      thforg1
c      thforg2
c      thfor1
c      thfor2
c      vhforg1:     heating force for 350 mb.
c      vhforg2:     heating force for 650 mb.
c
c      COMMON /cvorfor/ vfor1,vfor2,vfor3,vforg1,vforg2,vforg3
c      vfor1
c      vfor2
c      vfor3
c      vforg1:      vorticity forcing at 200 mb.
c      vforg2:      vorticity forcing at 500 mb.
c      vforg3:      vorticity forcing at 800 mb.
c
c      COMMON /cfluxpar/ evfac,evfaca
c      evfac:       maximum evaporation factor over land.
c      ecfaca:      actual evaporation factor over land.
c
c      COMMON /calbedo/ albice,albsea,albsnow,albland,albseaw,albseas,
c     *                abstow,abstos
c      albice:      albedo of ice.
c      albsea:      albedo of sea.
c      albsnow:     albedo of snow.
c      albland:     albedo of land.
c
c-----------------------------------------------------------------------

      integer   iqmtab,jqmtab,kqmtab
      parameter (iqmtab=50,jqmtab=20,kqmtab=20)

      real*8  gamgr(nlat,nlon)
      real*8  tempsg(nlat,nlon),temp2g(nlat,nlon),
     *        temp4g(nlat,nlon),pground(nlat,nlon)
      real*8  tempm(2),thform(2),temp2gm,temp4gm
      real*8  tsurf(nlat,nlon)
      integer lsmask(nlat,nlon)
      real*8  solarc,Q0(nlat)
      real*8  hesw1(nlat,nlon),hesw2(nlat,nlon),
     *        hesws(nlat,nlon),albes(nlat,nlon),
     *        albea(nlat,nlon),albeaw(nlat),albeas(nlat),
     *        abso1(nlat,nlon),abso2(nlat,nlon)
      real*8  ulrad1(nlat,nlon),ulrad2(nlat,nlon),
     *        ulrads(nlat,nlon),
     *        dlrads(nlat,nlon)
      real*8  pteta(21),pa(3,21),pb(3,21),pc(3,21),clfrac(nlat)
      real*8  pafac(3,21),pbfac(3,21),pcfac(3,21)
      real*8  Hflux(nlat,nlon),Eflux(nlat,nlon),dEflux(nlat,nlon)
      integer lseaice(nlat,nlon)
      integer landsnow(nlat,nlon)
      real*8  albice(nlat),albsea(nlat),albsnow(nlat),albland(nlat),
     *        albseaw(nlat),albseas(nlat),abstow(nlat),abstos(nlat)
      real*8  rowat,roair,cdrag,cwdrag,cpair,cvair,rlatvap,rlatcon,
     *        sboltz,rlatsub,rlatfus,cwater,gamma,rkappa,
     *        potfac1,potfac2,gamd,tzero,alphad,alphav,alphas
      real*8  cdragv(nlat,nlon),cdragw(nlat,nlon),richar(nlat,nlon)
      real*8  chsea(49),cdsea(49),chlan(49),cdlan(49)
      real*8  tetacr(nlat,nlon),gams(nlat,nlon),teta(nlat,nlon)
      real*8  thforg1(nlat,nlon),thforg2(nlat,nlon),
     *        thfor1(nsh2),thfor2(nsh2),vhforg1(nlat,nlon),
     *        vhforg2(nlat,nlon)
      real*8  vfor1(nsh2),vfor2(nsh2),vfor3(nsh2),
     *        vforg1(nlat,nlon),vforg2(nlat,nlon),
     *        vforg3(nlat,nlon)
      real*8  evfac,evfaca(nlat,nlon)
      real*8  relhum(nlat,nlon),q10(nlat,nlon),qsurf(nlat,nlon)
      real*8  rmoiss(nsh2),rmoisg(nlat,nlon),evap(nlat,nlon)
      real*8  dp2(nlat,nlon),dp1,tdifq,gpm500,relhmax,hmoisr,umoisr,
     *        rainmax
      integer ihavm, ivavm, imsink
      real*8  corain(nlat,nlon),dyrain(nlat,nlon),
     *        torain(nlat,nlon),vemoisg(nlat,nlon)
      real*8  co2(nlat,nlon)
      real*8  cc1,cc2,cc3,tqmimin,tqmi(0:iqmtab),dtqmi,rdtqmi,
     *        tqmjmin,tqmj(0:jqmtab),dtqmj,rdtqmj,
     *        tqmkmin,tqmk(0:kqmtab),dtqmk,rdtqmk
      real*4  qmtabel(0:iqmtab,0:jqmtab,0:kqmtab)
      real*8  tcc(nlat,nlon)
      real*8  relhcrit, relhfac
      real*8  dumt1(nlat,nlon,nvl),dumt2(nlat,nlon,nvl)
      real*8  dumu1(nlat,nlon,nvl),dumu2(nlat,nlon,nvl)

      common /cgamma/ gamgr
      common /ctemag/ tempsg,temp2g,temp4g,pground
      common /ctempm/ tempm,thform,temp2gm,temp4gm
      common /csurf/  tsurf,lsmask
      common /sunfr/  solarc,Q0
      common /cswrad/ hesw1,hesw2,hesws,albes,albea,albeaw,albeas,
     *                abso1,abso2
      common /clwrad/ ulrad1,ulrad2,ulrads,dlrads
      common /clwpar/ pteta,pa,pb,pc,clfrac,pafac,pbfac,pcfac
      common /csflux/ Hflux,Eflux,dEflux
      common /cice/   lseaice
      common /clandp/ landsnow
      common /calbedo/ albice,albsea,albsnow,albland,albseaw,albseas,
     *                abstow,abstos
      common /cpar1/  rowat,roair,cdrag,cwdrag,cpair,cvair,rlatvap,rlatcon,
     *                sboltz,rlatsub,rlatfus,cwater,gamma,rkappa,
     *                potfac1,potfac2,gamd,tzero,alphad,alphav,alphas
      common /cparv/  cdragv,cdragw,richar,chsea,cdsea,chlan,cdlan
      common /conv/   tetacr,gams,teta
      common /cthfor/ thforg1,thforg2,thfor1,thfor2,vhforg1,vhforg2
      common /cvorfor/ vfor1,vfor2,vfor3,vforg1,vforg2,vforg3
      common /cfluxpar/ evfac,evfaca
      common /cmoisa/ relhum,q10,qsurf,rmoiss,rmoisg,evap

      common /cmoisp/ dp2,dp1,tdifq,gpm500,relhmax,hmoisr,umoisr,
     *                rainmax, ihavm, ivavm,imsink
      common /crain/  corain,dyrain,torain,vemoisg
      common /trac/ co2
      common /cqmax/ cc1,cc2,cc3,tqmimin,tqmi,dtqmi,rdtqmi,
     *               tqmjmin,tqmj,dtqmj,rdtqmj,
     *               tqmkmin,tqmk,dtqmk,rdtqmk,qmtabel

      common /ccloud/ tcc, relhcrit, relhfac
      common /dumout/ dumt1,dumt2,dumu1,dumu2
