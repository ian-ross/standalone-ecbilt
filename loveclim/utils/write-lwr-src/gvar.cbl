!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!      common block with variables for calculating Green's
!      functions.
!
!     fluxda.....downward fluxes for CLEAR and CLOUDY sky
!     fluxua.....upward fluxes originated at the atmospere
!     fluxus.....upward fluxes originated at the surface
!     atl........heating at pl
!
!      tl.....temperature in Kelvin
!      ts.....surface temperature (K)
!      ql.....specific humidity in Kg/Kg
!      pECHAM..surface pressure field (hPa) on ECHAM grid 48x96
!      o3ECHAM.ozone mixing ratio (KG/KG) on ECHAM grid 48x96, 20 layers
!      pECHAMz..surface pressure field (hPa) zonally averaged
!      o3ECHAMz.ozone mixing ratio (KG/KG) zonally averaged
!      clECHAM.cloud cover (fraction) for 20 layers on ECHAM grid 48x96
!      lwECHAM.liquid water content (Kg/m2)for 20 layers on ECHAM grid 48x96
!      RegMask.ECHAM grid 48x96. Cell value points to region to which cell
!              belongs:
!               1-9: 9 function zonal bands
!               10: Greenland mountains
!               11: Rocky mountains
!               12: Himalaya
!               13: Andes
!               14: Antarctica mountains

!      o3l....ozone mixing ratio in Kg/Kg
!      cl.....cloud cover (fraction) for each layer
!      lwcl...liquid water content (Kg/m2) for each layer
!      co2l...carbon dioxide in ppmv
!      tclc...total cloud cover (fraction) reference profile
!
!      t......temperature of input reference, not interpolated to
!             Dorland structure (K)
!      p......pressure half levels before interpolation to Dorland
!             structure (hPa)
!      q......specific humidity of input reference, not interpolated to
!             Dorland structure (g/g)
!     o3.....ozone concentration of McClatchey input reference, not
!            interpolated to Dorland structure (g/m3)
!    po3....pressure coordinates of McClatchey ozone input reference,
!           not interpolated to Dorland structure (mb)
!                 DT0-
!    pCO2-CO2 conc. (ppmv)
!    pCH4-methane conc. (ppbv)
!    pN2O-nitrous-oxide (ppbv)
!    pCFC-concentrations CFCs (pptv) per CFC type (pptv):
!                      pCFC(1) - CFC-11
!                      pCFC(2) - CFC-12
!                      pCFC(3) - CFC-113
!                      pCFC(4) - CFC-114
!                      pCFC(5) - CFC-115
!   pHCFC-concentrations HCFCs (pptv) per type (pptv):
!                      pHCFC(1) - HCFC-22
!                      pHCFC(2) - HCFC-123
!                      pHCFC(3) - HCFC-124
!                      pHCFC(4) - HCFC-125
!                      pHCFC(5) - HCFC-134A
!                      pHCFC(6) - HCFC-141B
!                      pHCFC(7) - HCFC-142B
!                      pHCFC(8) - HCFC-143A
!                      pHCFC(9) - HCFC-152A
!   pCTC-concentration Carbon TetraChloride (CCl4) (pptv)
!   pMCF-concentration Methyl Chloroform (CH3CCl3) (pptv)
!
!   pname-name of present profile used for calculating
!         Green's functions. name has following format:
!                prof.sss.nnn.mm  sss: southward boundary of zonal band
!                                 nnn: northward boundary
!                boundaries are given as:
!                    hxx with h:  hemisphere "s" or "n"
!                            xx: degrees from equator
!                            mm:  month
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  REAL pl(nlm),pl2(nlmp1),tl(nlm),tl2(nlmp1)
  REAL o3l(nlm),co2l(nlm),cl(nlm),lwcl(nlm),fluxda(nlmp1,2), &
       & fluxua(nlmp1,2),fluxus(nlmp1,2),atl(nlm),ql(nlm), &
       & cldmm(0:ntype),t(profnlm),q(profnlm),p(profnlm+1), &
       & o3(o3profnlm),tclc,cc(3),tp(3),lwp(3), &
       & pCO2,pCH4,pN2O,pCFC(5),pHCFC(9),pCTC,pMCF
  REAL ISCCPopp(MAXBOX),ISCCPlat(MAXBOX),ISCCPlon(MAXBOX)
  REAL psISCCP(MAXBOX),pISCCPz(17,2), &
       & o3ECHAMz(12,20,17),o3ECHAM(20,MAXBOX), &
       & clISCCP(3,MAXBOX),lwISCCP(3,MAXBOX),tpISCCP(3,MAXBOX)
  REAL plz(nlm,17),pl2z(nlmp1,17), &
       & tlz(nlm,17),tl2z(nlmp1,17), &
       & qlz(nlm,17),qcolz(17)
  REAL ts,ps,t2m,q2m
  REAL fluxl(nfluxl)
  REAL THKCOE(7,2,9,12), PLBTAB(12), ZLBTAB(12,72,12)
  INTEGER NPLBTB
  INTEGER RegMask(MAXBOX),LSMask(MAXBOX)
  INTEGER npol,latstepeff,lonstepeff
  COMMON /greenvar/pl,pl2,tl,tl2,o3l,co2l,cl,lwcl,fluxda, &
       & fluxua,fluxus,atl,ql,cldmm,t,q,p,o3,tclc,cc,tp,lwp, &
       & pCO2,pCH4,pN2O,pCFC,pHCFC,pCTC,pMCF,ISCCPopp, &
       & ISCCPlat,ISCCPlon,psISCCP,pISCCPz,o3ECHAMz,o3ECHAM,clISCCP, &
       & lwISCCP,tpISCCP,plz,pl2z,tlz,tl2z,qlz,qcolz, &
       & ts,ps,t2m,q2m,fluxl,THKCOE,PLBTAB,ZLBTAB,NPLBTB, &
       & RegMask,LSMask,npol,latstepeff,lonstepeff
  CHARACTER*15 pname
  CHARACTER*17 pnamels
  CHARACTER(LEN=*), PARAMETER :: pdir='./lw-radiation/profiles/'
  COMMON /greenproc/pname,pnamels
