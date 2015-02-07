!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** File:     comphys.h
! *** Contents: Common declarations for physical part of atmospheric
! ***           model of ECbilt
!      common /cgamma/ gamgr,gamean,gamvar,gamsq
!      solarc:      solar constant.
!      common /ctemag/ tempsg,temp2g,temp4g,tempsgn
!      tempsg:      mean temperature at 10 meters height.
!      tempsgn:     temperature at 10 meters height for each surface
!                   type.
!      temp2g:      temperature at 350 mb.
!      temp4g:      temperature at 650 mb.
!
!      common /ctempm/ tempm,thform,temp2gm,temp4gm,
!                     tempsgm,tsurfm
!      tempm:       global mean temperature.
!      thform:      global mean
!      temp2gm:     global mean temperature at 350 mb.
!      temp4gm:     global mean temperature at 650 mb.
!      tempsgm:     global mean temperature at 10 meters height.
!      tsurfm:      global mean temperature at surface.
!
!      common /cpot/   potm,pot2g,pot4g
!      potm:        global mean potential temperature.
!      pot2g:       potential temperature at 350 mb.
!      pot4g:       potential temperature at 650 mb.
!
!      common /csurf/  tsurf,lsmask,sst,sstday,tsurfn,fractn
!      tsurf:       mean temperature at surface.
!      lsmask:      land sea mask.
!      sst:         sea surface tmeperature.???????????????
!      sstday:      sea surface temperature each day.?????
!      tsurfn       surface temperature for each surface type
!      fractn       fraction of each surface type
!      fracto       fraction of land/ocean
!
!      common /sunfr/  solarc,Q0,solardref
!      solarc;      sun constant.
!
!      common /cswrad/ hesw1,hesw2,hesws,albes,albea,albeaw,albeas,
!      *                abso1,abso2,albesn,heswsn
!      hesw1:       solar radiation heating rate.
!      hesw2:       solar radiation heating rate.
!      hesws:       mean solar radiation heating rate at the surface.
!      albes:       mean albedo of the surface
!      albea:       solar radiation reflective coefficient.
!      albeaw:      albedo of land in winter
!      albeas:      albedo of land in summer
!      abso1:       solar radiation absorbtion coefficient.
!      abso2:       solar radiation absorbtion coefficient.
!      albesn:      albedo of each surface type
!      heswsn:      solar radiation heating rate for each surface type
!
!      common /cpar1/  rowat,roair,cdrag,cpair,cvair,rlatvap,rlatcon,
!      *                sboltz,rlatsub,rlatfus,cwater,gamma,rkappa,ps,
!      *                aa,bb,gamd,tzero,alphad,alphav,alphas,
!      *                dragan,dragla,cdragvn,epss,cwdrag
!      rowat:       density of water.
!      roair:       mean air desity at sea level.
!      cdrag:       coefficient in sensible and latent air-sea heat flux.
!      cwdrag:      coefficient in wind stress.
!      cdragv:      coef. in sen. en lat. heat flux depending on roughness
!                   length and stability (Richardson number)
!      richar:      richardson number
!      dragane:     rotation in the comp. of the wind stress as a function
!                   of latitude
!      cpair:       specific heat of dry air at constant pressure.
!      cvair:       specific heat of dry air at constant volum.
!      rlatvap:     latent heat uptake due to evaporation.
!      rlatcon:     latent heat release due to condensition.
!      sboltz:      stefan-boltzmann constant.
!      rlatsub:     latent heat of sublimation.
!      rlatfus:     latent heat of fusion.
!      cwater:      4180
!      gamma:       cpair/cvair
!      rkappa:      =(gamma-1)/gamma
!      ps:          surface pressure.
!      aa:          =(350/1000)**rkappa.
!      bb:          =(650/1000)**rkappa
!      gamad:
!      tzero:       =273.15
!      epss :       = treshold for surface computation
!      alphad:      =roair*cdrag*cpair
!      alphav:      =roair*cdrag*rlatvap
!      alphas:      =roair*cdrag*rlatsub
!      albes:       surface albedo.
!      dragan:      maximum turning angle in the computation of wind stress
!      dragla:      latitude below which the turning angle is reduced
!      cdragvn:     drag coefficeint for each surface type (see cdragv)

!      common/ fluxcore/ corAN,corPN,corAC,corID,corAS,corPS,corAA
!      corAN       flux correction in the North Atlantic
!      corPN       flux correction in the North Pacific
!      corAC       flux correction in the Arctic
!      corID       flux correction in the Indian
!      corAS       flux correction in the South Atlantic
!      corPS       flux correction in the South Pacific
!      corAA       flux correction in the Southern Ocean


!      common /clwrad/ ulrad1,ulrad2,ulrads,dlrads,dlradsn,ulradsn
!      ulrad2:      upwards long wave radiation at ?.
!      ulrad1:      upwards long wave radiation at ?.
!      ulrads:      mean upwards long wave radiation at the surface.
!      ulradsn:     upwards long wave radiation at the surface.
!                   for each surface type
!      dlrads:      mean downwards long wave radiation at the surface.
!      dlradsn:     downwards long wave radiation at the surface
!                   for each surface type
!
!      common /csflux/  Hflux(nlat,nlon),Eflux(nlat,nlon),dEflux(nlat,nlon),
!                       hfluxn(nlat,nlon,ntyps),eflux(nlat,nlon,ntyps)
!      Hflux:       mean sensible heat flux at the surface.
!      Eflux:       mean latent heat flux at the surfac.
!      dEflux:      latent heat flux at the surfac?????????????????.
!      hfluxn:      sensible heat flux for each surface type
!      efluxn:      latent heat flux for each surface type
!
!      common /crain/  corain,dyrain,torain
!      corain;      convective rain.
!      dyrain:      dynamic rain.
!      torain:      total rain.
!
!      common /cmoist/ qsats,qsat4,relhum,qg
!      qsats:
!      qsat4
!      relhum:      relative humidity.
!      qg
!
!      common /moisturew/ relhmax, ihavm, ivavm, imsink
!      relhmax:     maximum relative humidity before saturation.
!      ihavm:       with (1) or without (0) horizontal divergence of moisture.
!      ivavm:       with (1) or without (0) vertical divergence of moisture.
!      imsink:      with (1) or without (0) the source or sink of moisture.
!
!      common /cmois/  rmoiss,rmoisg,precip,evap,dcmoisg,evapn
!      rmoisg:      specific humidity.
!      precip:      (not found)
!      evap:        evaporation.
!      evapn:       evaporation for each surface type
!      dcmoisg:
!      common /conv/   imoisr,imoisc,tetacr,gams,teta,convn
!      imoisr
!      imoisc
!      tetacr
!      gams
!      teta:        0.5*(pot2g - pot4g)
!      convn
!
!      common /cthfor/ thforg1,thforg2,thfor1,thfor2,vhforg1,vhforg2
!      thforg1
!      thforg2
!      thfor1
!      thfor2
!      vhforg1:     heating force for 350 mb.
!      vhforg2:     heating force for 650 mb.
!
!      common /cvorfor/ vfor1,vfor2,vfor3,vforg1,vforg2,vforg3
!      vfor1
!      vfor2
!      vfor3
!      vforg1:      vorticity forcing at 200 mb.
!      vforg2:      vorticity forcing at 500 mb.
!      vforg3:      vorticity forcing at 800 mb.
!
!      common /cfluxpar/ evfac,evfaca
!      evfac:       maximum evaporation factor over land.
!      ecfaca:      actual evaporation factor over land.
!
!      common /calbedo/ albice,albsea,albsnow,albland,albseaw,albseas,
!     *                abstow,abstos
!      albice:      albedo of ice.
!      albsea:      albedo of sea.
!      albsnow:     albedo of snow.
!      albland:     albedo of land.
!
!      common /cstype/ noc,nse,nld
!      noc:       type number of the ocean
!      nse:       type number of the sea ice
!      nld:       type number of the land
!-----------------------------------------------------------------------


      integer   iqmtab,jqmtab,kqmtab
      parameter (iqmtab=50,jqmtab=20,kqmtab=20)

      real*8  gamgr(nlat,nlon)
      real*8  temp0g(nlat,nlon),temp2g(nlat,nlon),temp4g(nlat,nlon), &
     &        tempm(0:2)
      real*8  solarc,Q0(nlat),solardref
      real*8  ecc,obl,omweb,ecc2,so,perh
      real*8  rowat,roair,cdrag,cwdrag,cpair,cvair,rlatvap,rlatcon, &
     &        sboltz,rlatsub,rlatfus,cwater,gamma,rkappa, &
     &        potfac1,potfac2,gamd,tzero,alphad,alphav,alphas
      real*8  corAN,corPN,corAC,corID,corAS,corPS,corAA
      real*8  cdragw(nlat,nlon),richar(nlat,nlon),dragane(nlat)
      real*8  chsea(49),cdsea(49),chlan(49),cdlan(49)
      real*8  tetacr(nlat,nlon),gams(nlat,nlon),teta(nlat,nlon)
      real*8  thforg1(nlat,nlon),thforg2(nlat,nlon), &
     &        thforg0(nlat,nlon),vhforg0(nlat,nlon),vhforg1(nlat,nlon), &
     &        vhforg2(nlat,nlon)
      real*8  vfor1(nsh2),vfor2(nsh2),vfor3(nsh2), &
     &        vforg1(nlat,nlon),vforg2(nlat,nlon), &
     &        vforg3(nlat,nlon)
      real*8  evfac
      real*8  relhum(nlat,nlon)
      real*8  rmoiss(nsh2),rmoisg(nlat,nlon),qmount(nlat,nlon)
      real*8  dp2,dp1,dp0,tdifq,gpm500,relhmax,hmoisr,umoisr, &
     &        rainmax
      integer ihavm, ivavm, imsink
      real*8  corain(nlat,nlon),dyrain(nlat,nlon),cormois(nlat,nlon), &
     &        torain(nlat,nlon),vemoisg(nlat,nlon), &
     &        cosnow(nlat,nlon),dysnow(nlat,nlon),tosnow(nlat,nlon)
      real*8  co2(nlat,nlon)
      real*8  cc1,cc2,cc3,tqmimin,tqmi(0:iqmtab),dtqmi,rdtqmi, &
     &        tqmjmin,tqmj(0:jqmtab),dtqmj,rdtqmj, &
     &        tqmkmin,tqmk(0:kqmtab),dtqmk,rdtqmk
      real*8  qmtabel(0:iqmtab,0:jqmtab,0:kqmtab)
      real*8  tcc(nlat,nlon),tccd(nlat,nlon)
      integer ndayws,iradcloud,iens,numens
      real*8  bup
      real*8  eccf,oblf,omwebf
      real*8  relhcrit, relhfac, emisoc,emisse,emisld
      real*8 albin,albis,albice,alphd,alphdi,alphs,cgren
      real*8  dumt1(nlat,nlon,nvl),dumt2(nlat,nlon,nvl)
      real*8  dumu1(nlat,nlon,nvl),dumu2(nlat,nlon,nvl)
      real*8  dragan,dragla,uv10rfx,uv10m,uv10rws
      real*8  uv10(nlat,nlon),uvw10(nlat,nlon)


!***  common /swrscheme/costref,salbref,swrref,swrcost,swrsalb,dayfr(nlat),
!***          kosz(nlat),solarf(nlat)
!***  costref contains regional averaged cosine of the zenith angle
!***  salbref contains regional averaged surface albedo from ????????
!***  swrref  contains reference short wave radiation fluxes from
!***          radiation calculated using KRCM with NCEP humidity, ipcc95
!***          greenhousegascontrations, isccpd2 cloud climatology and
!***          ECHAM4 ozone climatology
!***  swrcost linear sensitivity coefficient for SWR changes due to anomalies
!***          from the cosine of the zenith angle wrt the reference value
!***  swrsalb linear sensitivity coefficient for SWR changes due to anomalies
!***          from the surface albedo for clear sky and 3rd order polynomial fit
!***          for unity overcast fluxes.
!*** index 1: flux levels
!*** index 2: regions
!*** index 3: month's
!*** index 4: 0 clear sky, 1 cloudy sky, except for surface albedo: 1, 2,3
!***          correspond to coefficients of 3rd order polynomial fit
!*** region classification and flux level definition same as LWR

      real*4 costref(27,12), salbref(27,12)
      real*4 swrref(8,27,12,0:1)
      real*4 swrcost(8,27,12,0:1)
      real*4 swrsalb(8,27,12,0:3)
      real*8 dayfr(nlat),dso4(nlat,nlon)
      real*8 kosz(nlat)
      real*8 solarf(nlat)
      real*8 tas1(nlat,nlon)
      real*4 sulopt(nlon,nlat)

!***  common lwrscheme/tncep,qancep,ghgipcc,lwrref,lwrt,lwrts,lwrqa,lwrghg,irn

!***  tncep contains reference vertical temperature profiles from the
!***  ncep reanalysis (1982-1994) for 12 month's, 27 regions and 19 levels:
!***  10,20,30,50,70,100,150,200,250,300,400,500,600,700,850,925,1000,p2m,psurf
!***  except when surface pressure lower than 1000, then p2m is inserted
!***  at the appropriate position as given by ipl (see below)

!***  qancep:total precipitable water content below 500 mb from ncep, reduced
!***  in order to tune the LW fluxes (minus 15 %)

!***  ghgref: green house gas concentrations 1990 from ipcc '92:
!           1: pCO2-CO2 conc. (ppmv)
!           2: pCH4-methane conc. (ppbv)
!           3: pN2O-nitrous-oxide (ppbv)
!
!           pCFC-concentrations CFCs (pptv) per CFC type (pptv):
!           4: pCFC(1) - CFC-11
!           5: pCFC(2) - CFC-12
!           6: pCFC(3) - CFC-113
!           7: pCFC(4) - CFC-114
!           8: pCFC(5) - CFC-115
!
!           pHCFC-concentrations HCFCs (pptv) per type (pptv):
!           9: pHCFC(1) - HCFC-22
!          10: pHCFC(2) - HCFC-123
!          11: pHCFC(3) - HCFC-124
!          12: pHCFC(4) - HCFC-125
!          13: pHCFC(5) - HCFC-134A
!          14: pHCFC(6) - HCFC-141B
!          15: pHCFC(7) - HCFC-142B
!          16: pHCFC(8) - HCFC-143A
!          17: pHCFC(9) - HCFC-152A
!
!          18: pCTC-concentration Carbon TetraChloride (CCl4) (pptv)
!          19: pMCF-concentration Methyl Chloroform (CH3CCl3) (pptv)

!***  O3echam4(kg/kg): based on ECHAM4 climatology; pressure levels in file
!***  for each season: djf-mam-jja-son
!***  ccisccp: total cloud cover monthly climatology: derived from isccp, used
!***  in the linearisation procedure of KRCM, tuned
!***  (middle and high clouds -15%) for KRCM to give 'observed'  LWR fluxes.

!***  lwrref: reference long wave radiation fluxes, for four month's, jan,
!***  apr, july, oct, calculated with KRCM with
!***  data from several sources: temperature and humidity from ncep, ground
!***  pressure and cloud cover climatology from ISCCP, ozone from ECHAM4:
!***  values of the first index represent:  1: OLR,  2: upward at 200 mb,
!***  3: upward at 500 hPa, 4: upward at surface, 5: downward at 200 mb,
!***  6: downward at 500 mb, 7: downward at surface
!***  and fourth index corresponds to clearsky (0) and unity overcast (1).
!***  lwrrefc: same as lwrref but corrected for systematic difference of
!***  linearised scheme wrt KRCM for NCEP dataset 1982-1994
!***  lwrt: first partial derivative of lwr fluxes wrt temperature at different
!***  levels (see tncep and lwrref)
!***  lwrts: 4th order polynomial fit of lwr flux dependence on surface temperature
!***  lwrqa: first partial derivative of lwr fluxes wrt total precipitable water
!***  content, distributed according to  ncep mean vertical profiles
!***  lwrghg: fits of lwr flux dependence on several ghg concentrations
!***  irn: index of regions used for definition of the reference profiles and fits
!***  27: Green Land area
!***  26: Rocky Mountain area
!***  25: Himalaya area
!***  24: Andes area
!***  23: Antarctica area
!***  22: zonal mean 75N-90N land
!***  21: zonal mean 75N-90N sea
!***  20: zonal mean 60N-75N land
!***  19: zonal mean 60N-75N sea
!***  18: zonal mean 45N-60N land
!***  17: zonal mean 45N-60N sea
!***  16: zonal mean 30N-45N land
!***  15: zonal mean 30N-45N sea
!***  14: zonal mean 15N-30N land
!***  13: zonal mean 15N-30N sea
!***  12: zonal mean 15S-15N land
!***  11: zonal mean 15S-15N sea
!***  10: zonal mean 30S-15S land
!***   9: zonal mean 30S-15S sea
!***   8: zonal mean 45S-30S land
!***   7: zonal mean 45S-30S sea
!***   6: zonal mean 60S-45S land
!***   5: zonal mean 60S-45S sea
!***   4: zonal mean 75S-60S land
!***   3: zonal mean 75S-60S sea
!***   2: zonal mean 90S-75S land
!***   1: zonal mean 90S-75S sea
!***  ipl: index of pressure level in tncep containing 2mtr
!***       temperature
!***  pisccp: surface pressure anual mean, region averaged

      REAL*4  pisccp(27),pncep(17),z500ncep(27,12)
      REAL*4  tncep(19,27,12),qancep(27,12),ghgipcc(19),ccisccp(32,64,12)
      REAL*4  lwrref(7,27,4,0:1)
      real*8  tsi
      real*8  dtemp(18,nlat,nlon,2),ghg(19),o3, rlogtl(17),rlogts(27)
      real*8  tncep12(2,27,12)
      real*8  lwrflux(7,27,4,0:1,2)
      REAL*4  lwrt(7,18,27,4,0:1),lwrts(7,4,27,4,0:1)
      REAL*4  lwrqa(7,27,4,0:1),lwrqts(7,4,27,4,0:1)
      REAL*4  lwrghg(7,19,27,4,0:1)
      integer irn(nlat,nlon,2),ipl(27)
      real*8  AMPWIR,AMPEQIR,EXPIR,HPROFTROP,HPROFEQ,HPROFW
      real*8  HPROFAN,AMPANIR,HPROFAN2,AMPANIR2

      real*8  solarvol(12,4),solarm,solarcl(nlat)

      common /cgamma/ gamgr
      common /ctemag/ temp2g,temp4g,tempm,temp0g
      common /sunfr/  solarc,q0,omweb,ecc,obl,solarvol, &
     &                   solarm,solarcl,ecc2,so,perh,solardref
      common /irad/  iradcloud,iens,numens, &
     &                  emisoc,emisse,emisld,bup, &
     &                  albin,albis,albice,alphd,alphdi,alphs,cgren, &
     &                  eccf,oblf,omwebf

      common /cpar1/  rowat,roair,cpair,cvair,rlatvap,rlatcon, &
     &                sboltz,rlatsub,rlatfus,cwater,gamma,rkappa, &
     &                potfac1,potfac2,gamd,tzero,alphad,alphav,alphas
      common /cparv/ndayws,cdragw,richar,chsea,cdsea,chlan,cdlan, &
     &          dragane,cdrag,dragan,dragla,uv10rfx,uv10m,uv10rws, &
     &                uv10,uvw10,cwdrag
      common/ fluxcore/ corAN,corPN,corAC,corID,corAS,corPS,corAA
      common /conv/   tetacr,gams,teta
      common /cthfor/ thforg1,thforg2,thforg0,vhforg0,vhforg1,vhforg2
      common /cvorfor/ vfor1,vfor2,vfor3,vforg1,vforg2,vforg3
      common /cfluxpar/ evfac
      common /cmoisa/ relhum,rmoiss,rmoisg,qmount

      common /cmoisp/ dp2,dp1,dp0,tdifq,gpm500,relhmax,hmoisr,umoisr, &
     &                rainmax, ihavm, ivavm,imsink
      common /crain/  corain,dyrain,torain,vemoisg,cormois, &
     &                cosnow,dysnow,tosnow
      common /trac/ co2
      common /cqmax/ cc1,cc2,cc3,tqmimin,tqmi,dtqmi,rdtqmi, &
     &               tqmjmin,tqmj,dtqmj,rdtqmj, &
     &               tqmkmin,tqmk,dtqmk,rdtqmk,qmtabel

      common /ccloud/ tcc, relhcrit, relhfac, tccd
      common /dumout/ dumt1,dumt2,dumu1,dumu2

      common /swrscheme/costref,salbref,swrref,swrcost,swrsalb &
     &          ,dayfr,kosz,solarf,dso4,sulopt,tas1
      common /lwrscheme/tncep,qancep,ghgipcc,tsi, &
           & ccisccp,lwrref,lwrflux,lwrt,lwrts,lwrqts,lwrqa,lwrghg,irn,ipl
      common /modlwr/AMPWIR,AMPEQIR,EXPIR,HPROFTROP,HPROFEQ,HPROFW, &
     &                  HPROFAN,AMPANIR,HPROFAN2,AMPANIR2

      common /cvertint1 /dtemp,ghg,o3
      common /cvertint2 /pisccp,tncep12,pncep,rlogtl,rlogts,z500ncep
