












c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comphys.h                                                    
c *** Contents: Common declarations for physical part of atmospheric 
c ***           model of ECbilt
c      common /ec_cgamma/ gamgr,gamean,gamvar,gamsq
c      solarc:      solar constant.
c      common /ec_ctemag/ tempsg,temp2g,temp4g,tempsgn
c      tempsg:      mean temperature at 10 meters height.
c      tempsgn:     temperature at 10 meters height for each surface 
c                   type.
c      temp2g:      temperature at 350 mb.
c      temp4g:      temperature at 650 mb.
c
c      common /ec_ctempm/ tempm,thform,temp2gm,temp4gm,
c                     tempsgm,tsurfm
c      tempm:       global mean temperature.
c      thform:      global mean 
c      temp2gm:     global mean temperature at 350 mb.
c      temp4gm:     global mean temperature at 650 mb.
c      tempsgm:     global mean temperature at 10 meters height.
c      tsurfm:      global mean temperature at surface.
c
c      common /ec_cpot/   potm,pot2g,pot4g
c      potm:        global mean potential temperature.
c      pot2g:       potential temperature at 350 mb.
c      pot4g:       potential temperature at 650 mb.
c
c      common /ec_csurf/  tsurf,lsmask,sst,sstday,tsurfn,fractn
c      tsurf:       mean temperature at surface.
c      lsmask:      land sea mask.
c      sst:         sea surface tmeperature.???????????????
c      sstday:      sea surface temperature each day.?????
c      tsurfn       surface temperature for each surface type
c      fractn       fraction of each surface type
c      fracto       fraction of land/ocean
c
c      common /ec_sunfr/  solarc,Q0,solardref
c      solarc;      sun constant.
c
c      common /ec_cswrad/ hesw1,hesw2,hesws,albes,albea,albeaw,albeas,
c      *                abso1,abso2,albesn,heswsn
c      hesw1:       solar radiation heating rate.
c      hesw2:       solar radiation heating rate.
c      hesws:       mean solar radiation heating rate at the surface.
c      albes:       mean albedo of the surface
c      albea:       solar radiation reflective coefficient.
c      albeaw:      albedo of land in winter
c      albeas:      albedo of land in summer
c      abso1:       solar radiation absorbtion coefficient.
c      abso2:       solar radiation absorbtion coefficient.
c      albesn:      albedo of each surface type
c      heswsn:      solar radiation heating rate for each surface type
c
c      common /ec_cpar1/  rowat,roair,cdrag,cpair,cvair,rlatvap,rlatcon,
c      *                sboltz,rlatsub,rlatfus,cwater,gamma,rkappa,ps,
c      *                aa,bb,gamd,tzero,alphad,alphav,alphas,
c      *                dragan,dragla,cdragvn,epss,cwdrag
c      rowat:       density of water.
c      roair:       mean air desity at sea level.
c      cdrag:       coefficient in sensible and latent air-sea heat flux.
c      cwdrag:      coefficient in wind stress.
c      cdragv:      coef. in sen. en lat. heat flux depending on roughness
c                   length and stability (Richardson number)
c      richar:      richardson number
c      dragane:     rotation in the comp. of the wind stress as a function
c                   of latitude
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
c      epss :       = treshold for surface computation
c      alphad:      =roair*cdrag*cpair
c      alphav:      =roair*cdrag*rlatvap
c      alphas:      =roair*cdrag*rlatsub
c      albes:       surface albedo.
c      dragan:      maximum turning angle in the computation of wind stress
c      dragla:      latitude below which the turning angle is reduced
c      cdragvn:     drag coefficeint for each surface type (see cdragv)

c      common/ fluxcore/ corAN,corPN,corAC,corID,corAS,corPS,corAA
c      corAN       flux correction in the North Atlantic
c      corPN       flux correction in the North Pacific
c      corAC       flux correction in the Arctic      
c      corID       flux correction in the Indian     
c      corAS       flux correction in the South Atlantic  
c      corPS       flux correction in the South Pacific 
c      corAA       flux correction in the Southern Ocean 


c      common /ec_clwrad/ ulrad1,ulrad2,ulrads,dlrads,dlradsn,ulradsn
c      ulrad2:      upwards long wave radiation at ?.
c      ulrad1:      upwards long wave radiation at ?.
c      ulrads:      mean upwards long wave radiation at the surface.
c      ulradsn:     upwards long wave radiation at the surface.
c                   for each surface type
c      dlrads:      mean downwards long wave radiation at the surface.
c      dlradsn:     downwards long wave radiation at the surface 
c                   for each surface type
c
c      common /ec_csflux/  Hflux(nlat,nlon),Eflux(nlat,nlon),dEflux(nlat,nlon),
c                       hfluxn(nlat,nlon,ntyps),eflux(nlat,nlon,ntyps)
c      Hflux:       mean sensible heat flux at the surface.
c      Eflux:       mean latent heat flux at the surfac.
c      dEflux:      latent heat flux at the surfac?????????????????.
c      hfluxn:      sensible heat flux for each surface type
c      efluxn:      latent heat flux for each surface type
c
c      common /ec_crain/  corain,dyrain,torain
c      corain;      convective rain.
c      dyrain:      dynamic rain.
c      torain:      total rain.
c
c      common /ec_cmoist/ qsats,qsat4,relhum,qg
c      qsats:
c      qsat4
c      relhum:      relative humidity.
c      qg
c
c      common /ec_moisturew/ relhmax, ihavm, ivavm, imsink
c      relhmax:     maximum relative humidity before saturation.
c      ihavm:       with (1) or without (0) horizontal divergence of moisture.
c      ivavm:       with (1) or without (0) vertical divergence of moisture.
c      imsink:      with (1) or without (0) the source or sink of moisture.
c
c      common /ec_cmois/  rmoiss,rmoisg,precip,evap,dcmoisg,evapn
c      rmoisg:      specific humidity.
c      precip:      (not found)
c      evap:        evaporation.
c      evapn:       evaporation for each surface type
c      dcmoisg:
c      common /ec_conv/   imoisr,imoisc,tetacr,gams,teta,convn
c      imoisr
c      imoisc
c      tetacr
c      gams
c      teta:        0.5*(pot2g - pot4g)
c      convn
c
c      common /ec_cthfor/ thforg1,thforg2,thfor1,thfor2,vhforg1,vhforg2
c      thforg1
c      thforg2
c      thfor1
c      thfor2
c      vhforg1:     heating force for 350 mb.
c      vhforg2:     heating force for 650 mb.
c
c      common /ec_cvorfor/ vfor1,vfor2,vfor3,vforg1,vforg2,vforg3
c      vfor1
c      vfor2
c      vfor3
c      vforg1:      vorticity forcing at 200 mb.
c      vforg2:      vorticity forcing at 500 mb.
c      vforg3:      vorticity forcing at 800 mb.
c
c      common /ec_cfluxpar/ evfac,evfaca
c      evfac:       maximum evaporation factor over land.
c      ecfaca:      actual evaporation factor over land.
c
c      common /ec_calbedo/ albice,albsea,albsnow,albland,albseaw,albseas,
c     *                abstow,abstos 
c      albice:      albedo of ice.
c      albsea:      albedo of sea.
c      albsnow:     albedo of snow.
c      albland:     albedo of land.
c      
c      common /ec_cstype/ noc,nse,nld
c      noc:       type number of the ocean
c      nse:       type number of the sea ice
c      nld:       type number of the land
c-----------------------------------------------------------------------


      integer   iqmtab,jqmtab,kqmtab
      parameter (iqmtab=50,jqmtab=20,kqmtab=20)

      real*8  gamgr(nlat,nlon)
      real*8  temp0g(nlat,nlon),temp2g(nlat,nlon),temp4g(nlat,nlon),
     *        tempm(0:2)
      real*8  solarc,Q0(nlat),solardref
      real*8  ecc,obl,omweb,ecc2,so,perh
      real*8  rowat,roair,cdrag,cwdrag,cpair,cvair,rlatvap,rlatcon,
     *        sboltz,rlatsub,rlatfus,cwater,gamma,rkappa,
     *        potfac1,potfac2,gamd,tzero,alphad,alphav,alphas
      real*8  corAN,corPN,corAC,corID,corAS,corPS,corAA
      real*8  cdragw(nlat,nlon),richar(nlat,nlon),
     *        dragane(nlat)
      real*8  chsea(49),cdsea(49),chlan(49),cdlan(49)
      real*8  tetacr(nlat,nlon),gams(nlat,nlon),teta(nlat,nlon)
      real*8  thforg1(nlat,nlon),thforg2(nlat,nlon),
     *        thforg0(nlat,nlon),vhforg0(nlat,nlon),vhforg1(nlat,nlon),
     *        vhforg2(nlat,nlon)
      real*8  vfor1(nsh2),vfor2(nsh2),vfor3(nsh2),
     *        vforg1(nlat,nlon),vforg2(nlat,nlon),
     *        vforg3(nlat,nlon)
      real*8  evfac
      real*8  relhum(nlat,nlon)
      real*8  rmoiss(nsh2),rmoisg(nlat,nlon),qmount(nlat,nlon)
      real*8  dp2,dp1,dp0,tdifq,gpm500,relhmax,hmoisr,umoisr,
     *        rainmax
      integer ihavm, ivavm, imsink
      real*8  corain(nlat,nlon),dyrain(nlat,nlon),cormois(nlat,nlon),
     *        torain(nlat,nlon),vemoisg(nlat,nlon),
     *        cosnow(nlat,nlon),dysnow(nlat,nlon),tosnow(nlat,nlon)
      real*8  co2(nlat,nlon)
      real*8  cc1,cc2,cc3,tqmimin,tqmi(0:iqmtab),dtqmi,rdtqmi,
     *        tqmjmin,tqmj(0:jqmtab),dtqmj,rdtqmj,
     *        tqmkmin,tqmk(0:kqmtab),dtqmk,rdtqmk
      real*4  qmtabel(0:iqmtab,0:jqmtab,0:kqmtab)
      real*8  tcc(nlat,nlon),tccd(nlat,nlon)
      integer ndayws,iradcloud,iscenghg,iscenghg2s,iscentsi,iscenvol,
     *        iscensul,issulstrt,isceno3,
     *        iscencel,iens,numens
      real*8  facttsi,bup
      real*8  eccf,oblf,omwebf
      real*8  relhcrit, relhfac, emisoc,emisse,emisld
      real*8 albin,albis,albice,alphd,alphdi,alphs,cgren
      real*8  dumt1(nlat,nlon,nvl),dumt2(nlat,nlon,nvl)
      real*8  dumu1(nlat,nlon,nvl),dumu2(nlat,nlon,nvl)
      real*8  dragan,dragla,uv10rfx,uv10m,uv10rws
      real*8  uv10(nlat,nlon),uvw10(nlat,nlon)
      

c***  common /ec_swrscheme/costref,salbref,swrref,swrcost,swrsalb,dayfr(nlat),
c***          kosz(nlat),solarf(nlat)
c***  costref contains regional averaged cosine of the zenith angle
c***  salbref contains regional averaged surface albedo from ????????
c***  swrref  contains reference short wave radiation fluxes from 
c***          radiation calculated using KRCM with NCEP humidity, ipcc95 
c***          greenhousegascontrations, isccpd2 cloud climatology and
c***          ECHAM4 ozone climatology
c***  swrcost linear sensitivity coefficient for SWR changes due to anomalies
c***          from the cosine of the zenith angle wrt the reference value
c***  swrsalb linear sensitivity coefficient for SWR changes due to anomalies
c***          from the surface albedo for clear sky and 3rd order polynomial fit
c***          for unity overcast fluxes.
c*** index 1: flux levels
c*** index 2: regions
c*** index 3: month's
c*** index 4: 0 clear sky, 1 cloudy sky, except for surface albedo: 1, 2,3 
c***          correspond to coefficients of 3rd order polynomial fit
c*** region classification and flux level definition same as LWR
     
      real*4 costref(27,12), salbref(27,12)
      real*4 swrref(8,27,12,0:1)
      real*4 swrcost(8,27,12,0:1)
      real*4 swrsalb(8,27,12,0:3)
      real*8 dayfr(nlat),dso4(nlat,nlon)
      real*8 kosz(nlat)
      real*8 solarf(nlat)
      real*8 tas1(nlat,nlon)
      real*4 sulopt(nlon,nlat,6000)
      double precision :: suloptTime(6000)
      integer :: oldmonth

c***  common lwrscheme/tncep,qancep,ghgipcc,lwrref,lwrt,lwrts,lwrqa,lwrghg,irn      

C***  tncep contains reference vertical temperature profiles from the
c***  ncep reanalysis (1982-1994) for 12 month's, 27 regions and 19 levels:
c***  10,20,30,50,70,100,150,200,250,300,400,500,600,700,850,925,1000,p2m,psurf
c***  except when surface pressure lower than 1000, then p2m is inserted
c***  at the appropriate position as given by ipl (see below)

c***  qancep:total precipitable water content below 500 mb from ncep, reduced
c***  in order to tune the LW fluxes (minus 15 %)

c***  ghgref: green house gas concentrations 1990 from ipcc '92: 
C           1: pCO2-CO2 conc. (ppmv)
C           2: pCH4-methane conc. (ppbv)
C           3: pN2O-nitrous-oxide (ppbv)
C
C           pCFC-concentrations CFCs (pptv) per CFC type (pptv):
C           4: pCFC(1) - CFC-11
C           5: pCFC(2) - CFC-12
C           6: pCFC(3) - CFC-113
C           7: pCFC(4) - CFC-114
C           8: pCFC(5) - CFC-115
C
C           pHCFC-concentrations HCFCs (pptv) per type (pptv):
C           9: pHCFC(1) - HCFC-22
C          10: pHCFC(2) - HCFC-123
C          11: pHCFC(3) - HCFC-124
C          12: pHCFC(4) - HCFC-125
C          13: pHCFC(5) - HCFC-134A
C          14: pHCFC(6) - HCFC-141B
C          15: pHCFC(7) - HCFC-142B
C          16: pHCFC(8) - HCFC-143A
C          17: pHCFC(9) - HCFC-152A
C
C          18: pCTC-concentration Carbon TetraChloride (CCl4) (pptv)
C          19: pMCF-concentration Methyl Chloroform (CH3CCl3) (pptv)

c***  O3echam4(kg/kg): based on ECHAM4 climatology; pressure levels in file
c***  for each season: djf-mam-jja-son
c***  ccisccp: total cloud cover monthly climatology: derived from isccp, used
c***  in the linearisation procedure of KRCM, tuned 
c***  (middle and high clouds -15%) for KRCM to give 'observed'  LWR fluxes.

c***  lwrref: reference long wave radiation fluxes, for four month's, jan,
c***  apr, july, oct, calculated with KRCM with
c***  data from several sources: temperature and humidity from ncep, ground
c***  pressure and cloud cover climatology from ISCCP, ozone from ECHAM4:
c***  values of the first index represent:  1: OLR,  2: upward at 200 mb, 
c***  3: upward at 500 hPa, 4: upward at surface, 5: downward at 200 mb, 
c***  6: downward at 500 mb, 7: downward at surface
c***  and fourth index corresponds to clearsky (0) and unity overcast (1).
c***  lwrrefc: same as lwrref but corrected for systematic difference of
c***  linearised scheme wrt KRCM for NCEP dataset 1982-1994
c***  lwrt: first partial derivative of lwr fluxes wrt temperature at different
c***  levels (see tncep and lwrref)
c***  lwrts: 4th order polynomial fit of lwr flux dependence on surface temperature
c***  lwrqa: first partial derivative of lwr fluxes wrt total precipitable water
c***  content, distributed according to  ncep mean vertical profiles 
c***  lwrghg: fits of lwr flux dependence on several ghg concentrations
c***  irn: index of regions used for definition of the reference profiles and fits
c***  27: Green Land area
c***  26: Rocky Mountain area
c***  25: Himalaya area
c***  24: Andes area
c***  23: Antarctica area 
c***  22: zonal mean 75N-90N land
c***  21: zonal mean 75N-90N sea
c***  20: zonal mean 60N-75N land
c***  19: zonal mean 60N-75N sea
c***  18: zonal mean 45N-60N land
c***  17: zonal mean 45N-60N sea
c***  16: zonal mean 30N-45N land
c***  15: zonal mean 30N-45N sea
c***  14: zonal mean 15N-30N land
c***  13: zonal mean 15N-30N sea
c***  12: zonal mean 15S-15N land
c***  11: zonal mean 15S-15N sea
c***  10: zonal mean 30S-15S land
c***   9: zonal mean 30S-15S sea
c***   8: zonal mean 45S-30S land
c***   7: zonal mean 45S-30S sea
c***   6: zonal mean 60S-45S land
c***   5: zonal mean 60S-45S sea
c***   4: zonal mean 75S-60S land
c***   3: zonal mean 75S-60S sea
c***   2: zonal mean 90S-75S land
c***   1: zonal mean 90S-75S sea
c***  ipl: index of pressure level in tncep containing 2mtr
c***       temperature
c***  pisccp: surface pressure anual mean, region averaged
  
      real*4  tncep(19,27,12),qancep(27,12)
      real*4  ghgipcc(19),tsiscen(2,0:12000), ghgscen(20,0:12000), o3scen(2,0:12000)
      integer y1scenghg,nyscenmaxghg,y1scentsi,nyscenmaxtsi,y1sceno3,nyscenmaxo3
      integer y1scenvol,m1scenvol,nyscenmaxvol,nyscenmaxsul
      real*8  dtemp(18,nlat,nlon,2),ghg(19),o3, rlogtl(17),rlogts(27)
      real*4  pisccp(27),pncep(17),z500ncep(27,12)
c     real*4  o3echam4(20,27,12)
      real*4  ccisccp(32,64,12),tncep12(2,27,12)
      real*4  lwrref(7,27,4,0:1),lwrflux(7,27,4,0:1,2)
      real*4  lwrt(7,18,27,4,0:1),lwrts(7,4,27,4,0:1)
      real*4  lwrqa(7,27,4,0:1),lwrqts(7,4,27,4,0:1)
      real*4  lwrghg(7,19,27,4,0:1)
      integer irn(nlat,nlon,2),ipl(27)
      real*8  AMPWIR,AMPEQIR,EXPIR,HPROFTROP,HPROFEQ,HPROFW
      real*8  HPROFAN,AMPANIR,HPROFAN2,AMPANIR2

      real*8  solartsi(2101),solarvol(0:12000,12,4)
      real*8  solarm
      real*8  solarcl(nlat)

      common /ec_cgamma/ gamgr
      common /ec_ctemag/ temp2g,temp4g,tempm,temp0g
      common /ec_sunfr/  solarc,q0,omweb,ecc,obl,solartsi,solarvol,
     *                   solarm,solarcl,ecc2,so,perh,solardref
      common /ec_irad/  iradcloud,iscenghg,iscenghg2s,iscentsi,
     *                  iscenvol,iscensul,issulstrt,isceno3,
     *                  iscencel,iens,numens,
     *                  emisoc,emisse,emisld,facttsi,bup,
     *                  albin,albis,albice,alphd,alphdi,alphs,cgren,
     *                  eccf,oblf,omwebf

      common /ec_cpar1/  rowat,roair,cpair,cvair,rlatvap,rlatcon,
     *                sboltz,rlatsub,rlatfus,cwater,gamma,rkappa,
     *                potfac1,potfac2,gamd,tzero,alphad,alphav,alphas
      common /ec_cparv/ndayws,cdragw,richar,chsea,cdsea,chlan,cdlan,
     *          dragane,cdrag,dragan,dragla,uv10rfx,uv10m,uv10rws,
     *                uv10,uvw10,cwdrag
      common/ fluxcore/ corAN,corPN,corAC,corID,corAS,corPS,corAA
      common /ec_conv/   tetacr,gams,teta
      common /ec_cthfor/ thforg1,thforg2,thforg0,vhforg0,vhforg1,vhforg2
      common /ec_cvorfor/ vfor1,vfor2,vfor3,vforg1,vforg2,vforg3
      common /ec_cfluxpar/ evfac
      common /ec_cmoisa/ relhum,rmoiss,rmoisg,qmount

      common /ec_cmoisp/ dp2,dp1,dp0,tdifq,gpm500,relhmax,hmoisr,umoisr,
     *                rainmax, ihavm, ivavm,imsink
      common /ec_crain/  corain,dyrain,torain,vemoisg,cormois,
     *                cosnow,dysnow,tosnow
      common /ec_trac/ co2
      common /ec_cqmax/ cc1,cc2,cc3,tqmimin,tqmi,dtqmi,rdtqmi,
     *               tqmjmin,tqmj,dtqmj,rdtqmj,
     *               tqmkmin,tqmk,dtqmk,rdtqmk,qmtabel

      common /ec_ccloud/ tcc, relhcrit, relhfac, tccd
      common /ec_dumout/ dumt1,dumt2,dumu1,dumu2
      
      common /ec_swrscheme/costref,salbref,swrref,swrcost,swrsalb
     *          ,dayfr,kosz,solarf,dso4,sulopt,suloptTime,tas1,oldmonth
      common /ec_lwrscheme/tncep,qancep,ghgipcc,ghgscen,o3scen,tsiscen,ccisccp,lwrref,lwrflux
     *        ,lwrt,lwrts,lwrqts,lwrqa,lwrghg,irn,ipl
     *        ,y1scenghg,nyscenmaxghg,y1scentsi,nyscenmaxtsi,y1sceno3,nyscenmaxo3
     *        ,y1scenvol,m1scenvol,nyscenmaxvol,nyscenmaxsul
      common /ec_modlwr/AMPWIR,AMPEQIR,EXPIR,HPROFTROP,HPROFEQ,HPROFW,
     *                  HPROFAN,AMPANIR,HPROFAN2,AMPANIR2
     
      common /ec_cvertint1 /dtemp,ghg,o3
      common /ec_cvertint2 /pisccp,tncep12,pncep,rlogtl,rlogts,z500ncep

