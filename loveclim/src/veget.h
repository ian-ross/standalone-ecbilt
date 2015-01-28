












C  modif : 19/05/00 <- size of array replaced by (nlat,nlon)
C        : may 2007 <- different enrichment factors for tree and grass (/riche/)
C        : july 2008 <- deforestation scenario (/forets/, parameter mdfor) and cleaning (suppressed array data_crop)
C
      integer nlon, nlat, nvl, nsh, nm, nsh2,iveg,mdfor
      parameter( nm=21, nlon = 64 , nlat = 32 , nvl=3, 
     *            nsh=((nm+1)*(nm+2))/2, nsh2=2*nsh )
      parameter( iveg = 600 , mdfor = 3100)

      integer numvegvar
      real*8  phi(nlat),veg_fill_value,veg_missing_value,newvegvar(80,20)
      character*60 namevegvar(80,6)
      common /veg_netcdf/newvegvar,veg_fill_value,veg_missing_value,numvegvar,
     *                  namevegvar,phi

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c 2 files "bio.inc" & "buffer.inc" combined in one : "veget.h" (09/12/99)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

        REAL*8 a,bet,gamm,gamm2,fmax,avecube,tmin,npp,nppmax,v1,v2,v3,
     * c1t,c2t,c3t,c1g,c2g,c3g,
     * d1t,d2t,d3t,d1g,d2g,d3g,
     * e1t,e2t,e3t,e1g,e2g,e3g,
     * f1t,f2t,f3t,f1g,f2g,f3g,
     * k1t,k2t,k3t,k1g,k2g,k3g,
     * t1t,t2t,t3t,t4t,t1g,t2g,t3g,t4g,
     * ps1,ps2,ps3,ps4,ps5,soilt,
     * forshare_st,t1tn,t1td, desshare_st,nlshare_st,
     * deng,dentd,dentn,laig,lait,
     * ave_t,ave_pr,ave_pr05,desmin,desmax,
     * ave_pr05_des,
     * ades, acr, k0t, k0g, k4g
 
      integer lon, lat, init_flag
      real*8  b1t,b2t,b3t,b4t,b1g,b2g,b3g,b4g, gdd0,gdd0_min,gdd0_max,
     &  carea,co2ghg,acwd,acwt,acwg,acwn,zrd,zrt,zrg,zrn,rsd,rst,rsg,rsn 
 
        COMMON /BIODAT/
     *  a,bet,gamm,gamm2,fmax,avecube,tmin,npp,nppmax,v1,v2,v3,
     * c1t,c2t,c3t,c1g,c2g,c3g,
     * d1t,d2t,d3t,d1g,d2g,d3g,
     * e1t,e2t,e3t,e1g,e2g,e3g,
     * f1t,f2t,f3t,f1g,f2g,f3g,
     * b1t(nlat,nlon),b2t(nlat,nlon),b3t(nlat,nlon),b4t(nlat,nlon),
     * b1g(nlat,nlon),b2g(nlat,nlon),b3g(nlat,nlon),b4g(nlat,nlon),
     * k1t,k2t,k3t,k1g,k2g,k3g,
     * t1t,t2t,t3t,t4t,t1g,t2g,t3g,t4g,
     * ps1,ps2,ps3,ps4,ps5,soilt,
     * forshare_st,desshare_st,
     * nlshare_st,t1tn,t1td,
     * deng,dentd,dentn,laig,lait,desmin,desmax,
     * gdd0, gdd0_min,gdd0_max,
     * ave_t,ave_pr,ave_pr05,
     * carea(nlat,nlon),
     * ades, acr, k0t, k0g, k4g,co2ghg,
     * lon,lat,init_flag(nlat,nlon),
     * acwd,acwt,acwg,acwn,zrd,zrt,zrg,zrn,rsd,rst,rsg,rsn
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
*********************************************************************
***  BUFFER COMMON: DATA TRANSFER CLIMATE <-> TERRESTR VEG MODEL  ***
*********************************************************************
      integer nstat,ieqveg,ieqvegwr,iscendef,ivegstrt
      real*8 st,sg,sd,snlt,anup,blai,pnpp,b12,b34,b1,b2,b3,b4,
     >      anup_moy,stock,anuploch,stockloch,st_moy
      real*8 stR,sgR,sdR,snltR
Cjmc_D
C       COMMON /CLIMATE/ TATB(nlat,nlon),  PRCB(nlat,nlon),
C    >   PRCB5(nlat,nlon), GDD0B(nlat,nlon)
Cjmc_F
        COMMON /BIOTA/
     >   ST(nlat,nlon), SG(nlat,nlon), SD(nlat,nlon), SNLT(nlat,nlon),
     >   BLAI(nlat,nlon,2), PNPP(nlat,nlon),
     >   B12(nlat,nlon),   B34(nlat,nlon),
     >   B1(nlat,nlon), B2(nlat,nlon), B3(nlat,nlon), B4(nlat,nlon),
     >   ANUP_MOY(nlat,nlon),ANUP(nlat,nlon), STOCK(nlat,nlon),
     >   st_moy(nlat,nlon),
     >   NSTAT(nlat,nlon),STR(nlat,nlon), SGR(nlat,nlon), SDR(nlat,nlon), 
     >   SNLTR(nlat,nlon)	
	
        COMMON /BIOTA2/
     >   ANUPLOCH(nlat,nlon), STOCKLOCH(nlat,nlon)

       common /vegetint/
     >   ieqveg,ieqvegwr,iscendef,ivegstrt
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-AM (2007)
c CO2 enrichment :
c  betat for tree, betag for grass
c  nppt = tree npp, nppg = grass npp
       REAL*8 BETAG,BETAT,NPPG,NPPT
       COMMON /RICHE/BETAG,BETAT,NPPG,NPPT
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-AM (2008)
c Deforestation scenario (data read in VEGET.dat):
c  farea = fraction of the mesh occupied by crops according to R&F (1999)
c  ndfor = nbr of records in deforestation scenario (ndfor =< mdfor)
c  i0dfor = year-1 of first available data in VEGET.dat (date A.D.)
c Reference vegetation distribution against which deforestation is evaluated:
c  sd_const(nlat,nlon) = reference desert distribution.
c  st_const(nlat,nlon) = reference tree distribution.
c Would CO2 flux from vegetation influence atm. CO2 or not : fco2veg
       REAL*8 fco2veg,sd_const,st_const
       REAL*4 VegetTime, farea
       INTEGER i0dfor,ndfor
       COMMON /FORETS/farea(nlon,nlat,mdfor),sd_const(nlat,nlon),
     &                st_const(nlat,nlon),fco2veg,i0dfor,ndfor, VegetTime(mdfor)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
