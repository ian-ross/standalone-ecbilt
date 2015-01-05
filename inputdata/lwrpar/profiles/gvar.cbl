ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      common block with variables for calculating Green's
c      functions.
c
c     fluxda.....downward fluxes for CLEAR and CLOUDY sky
c     fluxua.....upward fluxes originated at the atmospere
c     fluxus.....upward fluxes originated at the surface
c     atl........heating at pl
c
c      tl.....temperature in Kelvin
c      ts.....surface temperature (K)
c      ql.....specific humidity in Kg/Kg
c	 pECHAM..surface pressure field (hPa) on ECHAM grid 48x96
c	 o3ECHAM.ozone mixing ratio (KG/KG) on ECHAM grid 48x96, 20 layers
c	 pECHAMz..surface pressure field (hPa) zonally averaged
c	 o3ECHAMz.ozone mixing ratio (KG/KG) zonally averaged
c      clECHAM.cloud cover (fraction) for 20 layers on ECHAM grid 48x96
c      lwECHAM.liquid water content (Kg/m2)for 20 layers on ECHAM grid 48x96
c      RegMask.ECHAM grid 48x96. Cell value points to region to which cell 
c	         belongs:
C               1-9: 9 function zonal bands
C               10: Greenland mountains
C               11: Rocky mountains
C               12: Himalaya
C               13: Andes
C               14: Antarctica mountains

c      o3l....ozone mixing ratio in Kg/Kg
c	 cl.....cloud cover (fraction) for each layer
C      lwcl...liquid water content (Kg/m2) for each layer
c      co2l...carbon dioxide in ppmv
c	 tclc...total cloud cover (fraction) reference profile
c
c	t......temperature of input reference, not interpolated to 
c		 Dorland structure (K)
c	p......pressure half levels before interpolation to Dorland 
c		 structure (hPa)
c	q......specific humidity of input reference, not interpolated to 
c		 Dorland structure (g/g)
C	o3.....ozone concentration of McClatchey input reference, not interpolated to 
c		 Dorland structure (g/m3)
C	po3....pressure coordinates of McClatchey ozone input reference, not interpolated to 
c		 Dorland structure (mb)
C                 DT0-
C    pCO2-CO2 conc. (ppmv)
C    pCH4-methane conc. (ppbv)
C    pN2O-nitrous-oxide (ppbv)
C    pCFC-concentrations CFCs (pptv) per CFC type (pptv):
C                      pCFC(1) - CFC-11
C                      pCFC(2) - CFC-12
C                      pCFC(3) - CFC-113
C                      pCFC(4) - CFC-114
C                      pCFC(5) - CFC-115
C   pHCFC-concentrations HCFCs (pptv) per type (pptv):
C                      pHCFC(1) - HCFC-22
C                      pHCFC(2) - HCFC-123
C                      pHCFC(3) - HCFC-124
C                      pHCFC(4) - HCFC-125
C                      pHCFC(5) - HCFC-134A
C                      pHCFC(6) - HCFC-141B
C                      pHCFC(7) - HCFC-142B
C                      pHCFC(8) - HCFC-143A
C                      pHCFC(9) - HCFC-152A
C   pCTC-concentration Carbon TetraChloride (CCl4) (pptv)
C   pMCF-concentration Methyl Chloroform (CH3CCl3) (pptv)
C                      
C   pname-name of present profile used for calculating
C	   Green's functions. name has following format:
C	   prof.sss.nnn.mm  sss: southward boundary of zonal band
C				  nnn: northward boundary
C			             boundaries are given as:
C					 hxx with h:  hemisphere "s" or "n"
C						    xx: degrees from equator
C				  mm:  month
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	REAL pl(nlm),pl2(nlmp1),tl(nlm),tl2(nlmp1)
      REAL o3l(nlm),co2l(nlm),cl(nlm),lwcl(nlm),fluxda(nlmp1,2),
     2  fluxua(nlmp1,2),fluxus(nlmp1,2),atl(nlm),ql(nlm),
     3  cldmm(0:ntype),t(profnlm),q(profnlm),p(profnlm+1),
     4  o3(o3profnlm),tclc,cc(3),tp(3),lwp(3),
     5  pCO2,pCH4,pN2O,pCFC(5),pHCFC(9),pCTC,pMCF
      REAL ISCCPopp(MAXBOX),ISCCPlat(MAXBOX),ISCCPlon(MAXBOX)
	REAL psISCCP(MAXBOX),pISCCPz(17,2),
     2  o3ECHAMz(12,20,17),o3ECHAM(20,MAXBOX),
     3  clISCCP(3,MAXBOX),lwISCCP(3,MAXBOX),tpISCCP(3,MAXBOX)
      REAL plz(nlm,17),pl2z(nlmp1,17),
     2  tlz(nlm,17),tl2z(nlmp1,17),
     3  qlz(nlm,17),qcolz(17)
	REAL ts,ps,t2m,q2m
	REAL fluxl(nfluxl)
      REAL THKCOE(7,2,9,12),   PLBTAB(12),
     &       ZLBTAB(12,72,12)
      INTEGER NPLBTB
      INTEGER RegMask(MAXBOX),LSMask(MAXBOX)
	INTEGER npol,latstepeff,lonstepeff
	COMMON /greenvar/pl,pl2,tl,tl2,o3l,co2l,cl,lwcl,fluxda
     2,fluxua,fluxus,atl,ql,cldmm,t,q,p,o3,tclc,cc,tp,lwp
     3,pCO2,pCH4,pN2O,pCFC,pHCFC,pCTC,pMCF,ISCCPopp
     4,ISCCPlat,ISCCPlon,psISCCP,pISCCPz,o3ECHAMz,o3ECHAM,clISCCP
     5,lwISCCP,tpISCCP,plz,pl2z,tlz,tl2z,qlz,qcolz
     6,ts,ps,t2m,q2m,fluxl,THKCOE,PLBTAB,ZLBTAB,NPLBTB
     7,RegMask,LSMask,npol,latstepeff,lonstepeff
      CHARACTER*15 pname
      CHARACTER*17 pnamels
      COMMON /greenproc/pname,pnamels

