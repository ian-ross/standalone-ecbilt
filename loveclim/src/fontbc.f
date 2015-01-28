












 
      subroutine fontbc(kideb,kiut)
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c
c This routine determines the time evolution of sea ice
c thickness, concentration and heat content due to the
c vertical and lateral thermodynamic accretion-ablation
c processes.
c
c---
c Ccpl [Ccp0] => ligne specifique a la version avec [sans] couplage .
c---
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'
      include 'thermo.com'
c
      dimension zep(nbpt),zqsat(nbpt),zemin(nbpt),zrchu(nbpt),
     &          ztfs(nbpt),zksdh(nbpt),ztbq(nbpt),zfs(nbpt),zab(nbpt),
     &          zqmax(nbpt),zindn(nbpt),zindg(nbpt),zfsup(nbpt),
     &          zfocea(nbpt),zffs(nbpt)
      dimension qctbqb(nbpt,2)
      dimension zdhg(nbpt),zdhb(nbpt),zalbqb(nbpt)
      dimension ykn(nbpt),ykg(nbpt)
      dimension rcpdt(nbpt),ztst(nbpt)
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  1) Tolerance parameters.                                            |
c-----------------------------------------------------------------------
c
      zeps  = 1.0e-20
      zeps0 = 1.0e-13
      zeps1 = 1.0e-06
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  2) If tbqb(ji,1) > tfsn, tbqb(ji,1) = tfsn.                 |
c     If tbqb(ji,2/3) > tfsg, tbqb(ji,2/3) = tfsg.
c-----------------------------------------------------------------------
c
      do 10 ji=kideb,kiut
 
        zind1        = max(zero,sign(one,hndif-hnbqb(ji)))
        zind23       = max(zero,sign(one,hgdif-hgbqb(ji)))
        qctbqb(ji,1) = max(zero,(tbqb(ji,1)-tfsn)*rcpn*hnbqb(ji))
     &                 *(1.0-zind1)
        qctbqb(ji,2) = (max(zero,(tbqb(ji,2)-tfsg)*rcpg
     &                *(hgbqb(ji)/2.0))+
     &                max(zero,(tbqb(ji,3)-tfsg)*rcpg
     &                *(hgbqb(ji)/2.0)))*
     &                (1.0-zind23)
        tbqb(ji,1)   = min(tfsn,tbqb(ji,1))
        tbqb(ji,2)   = min(tfsg,tbqb(ji,2))
        tbqb(ji,3)   = min(tfsg,tbqb(ji,3))
10    continue
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  3) Calculate some intermediate variables.                           |
c-----------------------------------------------------------------------
c
      do 30 ji=kideb,kiut
        zemin(ji)  = emig
        zignn      = max(zero,sign(one,hnzst-hnbqb(ji)))
        zignm      = max(zero,sign(one,hndif-hnbqb(ji)))
        zignm      = max(zignm,zignn)
        zig        = max(zero,sign(one,-hnbqb(ji)))
        zigm       = max(zero,sign(one,hgdif-hgbqb(ji)))
        ztfs(ji)   = (1.0-zig)*tfsn+zig*tfsg
c
c Effective conductivity
c
        xknw       = xkn/(xkn+xkg)
        xkgw       = xkg/(xkn+xkg)
        he         = xknw*hgbqb(ji)+xkgw*hnbqb(ji)
        zhe        = max(zero,sign(one,2.0*he-hth))
        heshth     = he/hth
c
c       Linear surface temperature profile between he=0 and he=0.5*hth
c
Cold    gdif       = (1.0-zhe)*heshth+zhe*0.5*(1.0+log(2.0*heshth))
C       gdif       = max(zeps,hakspl*log(hgbqb(ji)/(hakspl*hth)))
c
c       quadratic surface temperature profile between he=0 and he=0.5*hth
c
        gdif       = (1.0-zhe)*heshth*(2.0-heshth)
     &                +zhe*0.5*(1.5+log(2.0*heshth))
C       gdif       = min(gdif,1.5*one)
c
        gdif       = 1.0
        ykn(ji)    = gdif*xkn
        ykg(ji)    = gdif*xkg
c
        zk         = 2.0*((1.0-zignm)*ykn(ji)+2.0*zignm*ykg(ji))
        zdh        = (1.0-zignm)*hnbqb(ji)+
     &               zignm*((1.0+3.0*zigm)*hgbqb(ji)
     &               +4.0*ykg(ji)/ykn(ji)*hnbqb(ji))
        zksdh(ji)  = zk/zdh
        ztbq(ji)   = (1.0-zignm)*tbqb(ji,1)
     &               +zignm*(tbqb(ji,2)*(1.0-zigm)
     &               +tfub(ji)*zigm)
30    continue
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  4) Calculate zqmax, fstbqb, qstbqb AND zab.                         |
c-----------------------------------------------------------------------
c
      do 40 ji=kideb,kiut
c
        zhg        = max(zero,sign(one,-hnbqb(ji)))
        zigm       = max(zero,sign(one,hgdif-hgbqb(ji)))
c
c                       GRENFELL AND MAYKUT, 1977.
c
        zst        = zhg*
     &               ((0.18+0.82*max(zero,1.0-(hgbqb(ji)/0.1)))
     &               *(1.0-cldqb(ji))+
     &               (0.35+0.65*max(zero,1.0-(hgbqb(ji)/0.1)))
     &               *cldqb(ji))
        zqmax(ji)  = max(zero,0.5*xlg*(hgbqb(ji)-hgmin))
        zexpo      = min(one,exp(-1.5*(hgbqb(ji)-0.1)))
C       zexpo      = zigm+(1.0-zigm)*min(one,
C    &               exp(-1.5*(hgbqb(ji)-0.1)))
        fstbqb(ji) = zst*fsolgb(ji)*zexpo
        zfsab      = zst*fsolgb(ji)*(1.0-zexpo)
        zhq        = zigm+(1.0-zigm)*max(zero,
     &               sign(one,qstbqb(ji)-zqmax(ji)))
        qstbqb(ji) = (qstbqb(ji)+(1.0-zhq)*zfsab*ddtb)*swiqst
        zab(ji)    = 1.0-zst*
     &               (zhq*zexpo+(1.0-zhq)*(swiqst+(1.0-swiqst)*zexpo))
       if (icoupl .eq. 0) then
        rhoa       = psbqb(ji)/(287.04*tabqb(ji))
C       ce         = ceb(vabqb(ji),tabqb(ji)-tsb(ji))
        ce         = 1.75e-03
C       ce         = 1.30e-03
        zrchu(ji)  = rhoa*ce*vabqb(ji)
       else
        zrchu(ji)  = 0.0
       endif
40    continue
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  5) Calculate surface temperature (Newton-Raphson method).           |
c-----------------------------------------------------------------------
c
cNCAR: Newton-Raphson surface temperature.
c
      surthi = 0.10
      do 45 ji=kideb,kiut
C       es         =  611.0*10.0**(9.5*(tsb(ji)-273.16)
C    &              /(tsb(ji)-7.66))
C       zqsat(ji)  =  (0.622*es)/(psbqb(ji)-(1.0-0.622)*es)
C       fcsb(ji)   =  zrchu(ji)*1004.0*(tsb(ji)-tabqb(ji))
C       fleb(ji)   =  zrchu(ji)*2.834e+06*(zqsat(ji)-qabqb(ji))
        rcpdt(ji)  =  (rcpn*min(hnbqb(ji),surthi)+
     &                 rcpg*max(surthi-hnbqb(ji),zero) )/ddtb
        ztst(ji)   =  tsb(ji)
 45   continue
      do 60 jt=1,nbits
        do 50 ji=kideb,kiut
         if (icoupl .eq. 0) then
          fratsb(ji) =  zemin(ji)*(ratbqb(ji)-
     &                  stefan*tsb(ji)*tsb(ji)*tsb(ji)*tsb(ji))
          es         =  611.0*10.0**(9.5*(tsb(ji)-273.16)
     &                  /(tsb(ji)-7.66))
          zqsat(ji)  =  (0.622*es)/(psbqb(ji)-(1.0-0.622)*es)
          fcsb(ji)   =  zrchu(ji)*1004.0*(tsb(ji)-tabqb(ji))
          fleb(ji)   =  zrchu(ji)*2.834e+06*(zqsat(ji)-qabqb(ji))
          zssdqw     =  (zqsat(ji)*zqsat(ji)*psbqb(ji)/
     &                  (0.622*es))*alog(10.0)*
     &                  9.5*((273.16-7.66)/(tsb(ji)-7.66)**2)
          dzf        =  4.0*zemin(ji)*stefan*tsb(ji)
     &                  *tsb(ji)*tsb(ji)+
     &                  zrchu(ji)*(1004.0+2.834e+06*zssdqw)+
     &                  zksdh(ji)
     &                  +rcpdt(ji)
         else
          fratsb(ji) =  ratbqb(ji)-(zemin(ji)*
     &                  stefan*tsb(ji)*tsb(ji)*tsb(ji)*tsb(ji))
          dzf        =  4.0*zemin(ji)*stefan*tsb(ji)
     &                  *tsb(ji)*tsb(ji)+
     &                  zksdh(ji)
     &                  +rcpdt(ji)+vabqb(ji)
Ccp2 &                  -fderb(ji)
         fcsb(ji)    = 0.0
         fleb(ji)    = 0.0
         endif
          zfs(ji)    =  zksdh(ji)*(ztbq(ji)-tsb(ji))
          zf         = -zab(ji)*fsolgb(ji)-fratsb(ji)
     &                  +fcsb(ji)+fleb(ji)-zfs(ji)
     &                  +rcpdt(ji)*(tsb(ji)-ztst(ji))
          dtsb       = -zf/dzf
C         dtsb       =  dtsb/3.0 
          dtemax     =  10.0
C         dtsb       =  min (max(dtsb/3.0,dtemax*-1.0),dtemax)
          dtsb       =  min (max(dtsb,dtemax*(-1.0)),dtemax)
          tsb(ji)    =  tsb(ji)+dtsb
c
cEND NCAR: Newton-Raphson surface temperature.
c
50      continue
60    continue
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  6) Limitation of surface temperature and update of sensible         |
c     and latent heat fluxes.                                          |
c-----------------------------------------------------------------------
c
      do 70 ji=kideb,kiut
        tsb(ji)    = min(ztfs(ji),tsb(ji))
c
cNCAR: Update surface fluxes.
c
       if (icoupl .eq. 0) then
        fratsb(ji) = zemin(ji)*(ratbqb(ji)-stefan*tsb(ji)*tsb(ji)
     &               *tsb(ji)*tsb(ji))
        es         = 611.0*10.0**(9.5*(tsb(ji)-273.16)/(tsb(ji)-7.66))
        zqsat(ji)  = (0.622*es)/(psbqb(ji)-(1.0-0.622)*es)
        fcsb(ji)   = zrchu(ji)*1004.0*(tsb(ji)-tabqb(ji))
        fleb(ji)   = zrchu(ji)*2.834e+06*(zqsat(ji)-qabqb(ji))
       else
        fratsb(ji) = ratbqb(ji)-(zemin(ji)*stefan*tsb(ji)*tsb(ji)
     &               *tsb(ji)*tsb(ji))
       endif
        zfs(ji)    = zksdh(ji)*(ztbq(ji)-tsb(ji))
c
cNCAR: Update surface fluxes.
c
70    continue
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  7)Calculate available heat for surface ablation.                    |
c-----------------------------------------------------------------------
c
      do 80 ji=kideb,kiut
        zffs(ji)= fratsb(ji)-fleb(ji)-fcsb(ji)
     &            +zfs(ji)+zab(ji)*fsolgb(ji)
        zffs(ji)= max(zero,zffs(ji))
        zffs(ji)= zffs(ji)*max(zero,sign(one,tsb(ji)-ztfs(ji)))
80    continue
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  8) Calculate changes in internal temperature due to                 |
c     vertical diffusion processes.                                    |
c-----------------------------------------------------------------------
c
      do 90 ji=kideb,kiut
c
c--8.1 Calculate intermediate variables.
c----------------------------------------
c
        umb        = 1.0-beta
        zind1      = max(zero,sign(one,hndif-hnbqb(ji)))
        zindn(ji)  = 1.0-zind1
        zind3      = max(zero,sign(one,hgdif-hgbqb(ji)))
        zindg(ji)  = zind1*zind3
        zks        = 2.0*ddtb*ykn(ji)/rcpn
        zkg        = 4.0*ddtb*ykg(ji)*(1.0-zindg(ji))/
     &               max(hgbqb(ji)*hgbqb(ji)*rcpg,zeps)
        zd         = ykn(ji)*hgbqb(ji)+2.0*ykg(ji)*hnbqb(ji)*zindn(ji)
        zaux1      = 2.0*zindn(ji)*zks*ykg(ji)/zd
        zaux2      = 2.0*zindn(ji)*zkg*ykn(ji)*hgbqb(ji)/zd
        zfsb       = ddtb*zfs(ji)
c
c--8.2. Fulfill the linear system matrix.
c-----------------------------------------
c
        zkgb       =  zkg*beta
        za1        = -zaux1*beta
        za2        = -zkgb
        za3        =  0.0
        zc1        =  0.0
        zc2        = -zaux2*beta
        zc3        = -zkgb
        zb1        =  zind1+hnbqb(ji)*zindn(ji)-za1
        zb2        =  1.0+zkgb-zc2
        zb3        =  1.0+3.0*zkgb
c
c--8.3. Fulfill the independent term vector.
c-------------------------------------------
c
        zd1        = hnbqb(ji)*zindn(ji)*tbqb(ji,1)
     &               -zindn(ji)*zfsb/rcpn+
     &               zind1*tsb(ji)+umb*zaux1*(tbqb(ji,2)-tbqb(ji,1))
        zd2        = tbqb(ji,2)-
     &               2.0*zind1*zfsb*(1.0-zindg(ji))/
     &               max(hgbqb(ji)*rcpg,zeps)+
     &               umb*(zaux2*(tbqb(ji,1)-tbqb(ji,2))+
     &                    zkg*(tbqb(ji,3)-tbqb(ji,2)))
        zd3        = tbqb(ji,3)+
     &               zkg*(2.0*tfub(ji)+umb*(tbqb(ji,2)-3.0*tbqb(ji,3)))
c
c--8.4. Solve the system (Gauss elimination method).
c----------------------------------------------------
c
        za1        = -za1/zb1
        zd1        =  zd1/zb1
        za2        = -za2/(zb2+zc2*za1)
        zd2        =  (zd2-zc2*zd1)/(zb2+zc2*za1)
        tbqb(ji,3) =  (zd3-zc3*zd2)/(zb3+zc3*za2)
        tbqb(ji,2) =  za2*tbqb(ji,3)+zd2
        tbqb(ji,1) =  za1*tbqb(ji,2)+zd1
        ztci       =  tbqb(ji,2)*(1.0-zindg(ji))+tfub(ji)*zindg(ji)
        zti        =
     &  (4.0*ykg(ji)*hnbqb(ji)*ztci+ykn(ji)*hgbqb(ji)
     &  *tsb(ji)*(1.0+3.0*zindg(ji)))/
     &  (4.0*ykg(ji)*hnbqb(ji)+ykn(ji)*hgbqb(ji)*(1.0+3.0*zindg(ji)))
        tbqb(ji,2) =  tbqb(ji,2)+zindg(ji)
     &                *((3.0*zti+tfub(ji))/4.0-tbqb(ji,2))
        tbqb(ji,3) =  tbqb(ji,3)+zindg(ji)
     &                *((zti+3.0*tfub(ji))/4.0-tbqb(ji,3))
        tbqb(ji,1) =  tbqb(ji,1)+zind1*(zti-tsb(ji))/2.0
c
90    continue
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  9.Take into account surface ablation and bottom accretion-abalation.|
c----------------------------------------------------------------------
c
      do 100 ji=kideb,kiut
c
c--9.1. Surface ablation and update of snow thickness and qstbqb
c------------------------------------------------------------------
c
c  surface ablation
c
        zdhsm      = -(zffs(ji)*ddtb-qctbqb(ji,1))/xln
        zab(ji)    =  hnbqb(ji)
        zhf        =  hnbqb(ji)+hnpbqb(ji)+zdhsm
        zhn        =  1.0-max(zero,sign(one,-zhf))
        hnbqb(ji)  =  max(zero,zhf)
        dvsbqb(ji) =  (1.0-albqb(ji))*(hnbqb(ji)
     &                -zab(ji)-hnpbqb(ji))
        dvsbqb(ji) =  min(zero,dvsbqb(ji))
        dmnbqb(ji) =  rhon*dvsbqb(ji)
        zqfont     =  max(zero,-zhf)*xln
        zdhg(ji)   = -zqfont/xlg
        zqf        =  qstbqb(ji)+zdhg(ji)*xlg
        ziqf       =  max(zero,sign(one,zqf))
        zhq        =  max(zero,sign(one,qstbqb(ji)-zqmax(ji)))
        zdhg(ji)   =  zdhg(ji)+zhq*(ziqf*zdhg(ji)
     &              -(1.0-ziqf)*qstbqb(ji)/xlg)
        qstbqb(ji) =  qstbqb(ji)+zhq*(ziqf*zqf-qstbqb(ji))
        dvsbqb(ji) =  zdhg(ji)*(1.0-albqb(ji))
c
c If fleb is negative, snow condensates at the surface.
c
        zdhcf      =  hnbqb(ji)-parsub*fleb(ji)/(rhon*xsn)*ddtb
        hnbqb(ji)  =  max(zero,zdhcf)
        zhn        =  1.0-max(zero,sign(one,-hnbqb(ji)))
        zdhg(ji)   =  zdhg(ji)-max(zero,-zdhcf)*rhon/rhog
c
c  Internal temperature and qstbqb.
c
c
        tbqb(ji,1) = zhn*tbqb(ji,1)+(1.0-zhn)*tfub(ji)
        zhnf       = max(zero,sign(one,hnbqb(ji)-zab(ji)))
        tbqb(ji,1) = tbqb(ji,1)+
     &               (1.0-zab(ji)/max(hnbqb(ji),zeps))*
     &               (tsb(ji)-tbqb(ji,1))*zindn(ji)*zhnf
        zqnes      = (tfsg-tbqb(ji,2))*rcpg*(hgbqb(ji)/2.0)
        zqr        = qstbqb(ji)-zqnes
        ziqr       = max(zero,sign(one,zqr))
        tbqb(ji,2) = ziqr*tfsg+(1-ziqr)*(tbqb(ji,2)+
     &               qstbqb(ji)/(rcpg*(hgbqb(ji)/2)))
        qstbqb(ji) = ziqr*zqr*swiqst
c
c--9.2. Calculate bottom accretion-ablation and update qstbqb.
c--------------------------------------------------------------
c
        ztb        = tbqb(ji,3)*(1.0-zindg(ji))+tsb(ji)*zindg(ji)
        zf2        = 4.0*ykg(ji)*(tfub(ji)-ztb)/
     &               (hgbqb(ji)+zindg(ji)*(3.*hgbqb(ji)
     &                +4.0*ykg(ji)/ykn(ji)*hnbqb(ji)))
        zdhb(ji)   = ddtb*(zf2-fbbqb(ji)-qlbbqb(ji))/xlg
        zqf        = qstbqb(ji)+zdhb(ji)*xlg
        ziqf       = max(zero,sign(one,zqf))
        zihq       = max(zero,sign(one,qstbqb(ji)-zqmax(ji)))
        zidhb      = max(zero,sign(one,-zdhb(ji)))
        zdhb(ji)   = zdhb(ji)+zihq*zidhb*(ziqf*zdhb(ji)-
     &               (1.0-ziqf)*qstbqb(ji)/xlg)
     &               -qctbqb(ji,2)/xlg
        qstbqb(ji) = qstbqb(ji)+zihq*zidhb*(ziqf*zqf-qstbqb(ji))
        zdhbf      = max(hmelt,zdhb(ji))
        zfsup(ji)  = xlg*(1.0-albqb(ji))/albqb(ji)
     &               *(zdhbf-zdhb(ji))/ddtb
        zdhb(ji)   = zdhbf
        zhgnew     = hgbqb(ji)+zdhg(ji)+zdhb(ji)
        dvbbqb(ji) = (1.0-albqb(ji))*zdhb(ji)
c
c--9.3. Case of total ablation.
c--------------------------------
c
        zhni       =  hnbqb(ji)
        zihgnew    =  1.0-max(zero,sign(one,-zhgnew))
        zihg       =  max(zero,sign(one,-zhni))
        zdhgm      =  (1.0-zihgnew)*(zhgnew-qstbqb(ji)/xlg)
        zdhnm      =  (1.0-zihg)*zdhgm*rhog/rhon
        zhgnew     =  zihgnew*zhgnew
        dmgbqb(ji) =  dmgbqb(ji)+(1.0-albqb(ji))
     &                *(zhgnew-hgbqb(ji))*rhog
        zhnfi      =  zhni+zdhnm
        hnbqb(ji)  =  max(zero,zhnfi)
        dmnbqb(ji) =  dmnbqb(ji)+(1.0-albqb(ji))*(hnbqb(ji)-zhni)*rhon
        zfocea(ji) = -(zihg*zdhgm*xlg+(zhnfi-hnbqb(ji))*xln)/ddtb
        qstbqb(ji) =  zihgnew*qstbqb(ji)
c
c--9.4. Update internal temperature and ice thickness.
c-------------------------------------------------------
c
        tsb(ji)    =  zihgnew*tsb(ji)+(1.0-zihgnew)*tfub(ji)
        zidhb      =  max(zero,sign(one,-zdhb(ji)))
        zc1        = -zhgnew*0.5
        zpc1       =  min(0.5*one,-hgbqb(ji)*0.5-zdhg(ji))
        zc2        = -zhgnew
        zpc2       =  zidhb*zc2+(1.0-zidhb)*(-hgbqb(ji)-zdhg(ji))
        zp1        =  max(zpc1,zc1)
        zp2        =  max(zpc2,zc1)
        zep(ji)    =  tbqb(ji,2)
        tbqb(ji,2) =
     &            2.0*(-zp1*tbqb(ji,2)+(zp1-zp2)*tbqb(ji,3)
     &            +(zp2-zc1)*tfub(ji))*
     &            zihgnew/max(zhgnew,zeps)+(1.0-zihgnew)*tfub(ji)
        zp1        =  min(zpc1,zc1)
        zp2        =  min(zpc2,zc1)
        zp1        =  max(zc2,zp1)
        tbqb(ji,3) =  2.0*
     &                ((1.0-zidhb)*((zc1-zp2)*tbqb(ji,3)
     &                +(zp2-zc2)*tfub(ji))+
     &                 zidhb*((zc1-zp1)*zep(ji)
     &                 +(zp1-zc2)*tbqb(ji,3)))*
     &                zihgnew/max(zhgnew,zeps)+(1.0-zihgnew)*tfub(ji)
        hgbqb(ji)  =  zhgnew
c
100   continue
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  10) Surface accretion.                                              |
c-----------------------------------------------------------------------
c
      do 110 ji=kideb,kiut
c
c  Archimedes principle:
c       zrmh       = (rhon*hnbqb(ji)+rhog*hgbqb(ji))/rho0
c
c  Lepparanta (1983):
C       dzrmh      = max(zero,(rhon*hnbqb(ji)+(rhog-rho0)
C    &               *hgbqb(ji))/(rhon+rho0-rhog))
        dzrmh      = max(zero,(rhon*hnbqb(ji)+(rhog-rho0)
     &               *hgbqb(ji))/(alphs*rhon+rho0-rhog))
c
c  New ice thickness.
c       zhgnew     = max(hgbqb(ji),zrmh)
c       zhnnew     = (rho0*zrmh-rhog*zhgnew)*rhoesn
c
c Lepparanta (1983):
c
        zhgnew     = max(hgbqb(ji),hgbqb(ji)+dzrmh)
C       zhnnew     = min(hnbqb(ji),hnbqb(ji)-dzrmh)
        zhnnew     = min(hnbqb(ji),hnbqb(ji)-alphs*dzrmh)
c
        zig        = 1.0-max(zero,sign(one,-zhgnew))
c
c  Compute new ice temperatures. snow temperature remains unchanged.
c
        quot       = (1.0-zig)+zig*min(one,hgbqb(ji)/max(zhgnew,zeps))
c       tneq       = cnscg*tbqb(ji,1)
c
c  Lepparanta (1983):
c
C       tneq       = cnscg*tbqb(ji,1)+(1.0-rhon/rhog)*tfub(ji)
        tneq       = alphs*cnscg*tbqb(ji,1)+
     &                   (1.0-alphs*(rhon/rhog))*tfub(ji)
c
        zep(ji)    = tbqb(ji,2)
c
c  Lepparanta (1983) (latent heat released during white ice formation
c  goes to the ocean -for lateral ablation-)
C       qlbqb(ji)  = qlbqb(ji)+dzrmh*(1.0-rhon/rhog)*xlg*(1.0-albqb(ji))
        qlbqb(ji)  = qlbqb(ji)+dzrmh*(1.0-alphs*(rhon/rhog))
     &                                      *xlg*(1.0-albqb(ji))
c
        tbqb(ji,2) = tneq-quot*quot*(tneq-tbqb(ji,2))
        tbqb(ji,3) = 2.0*tneq+quot*(tbqb(ji,3)
     &                   +zep(ji)-2.0*tneq)-tbqb(ji,2)
c
c  Changes in ice volume and ice mass.
c
        dvnbqb(ji) = (1.0-albqb(ji))*(zhgnew-hgbqb(ji))
        dmgwib(ji) = dmgwib(ji)+(1.0-albqb(ji))
     &                *(hnbqb(ji)-zhnnew)*rhon
c
c  Lepparanta (1983):
c
        dmgbqb(ji) = dmgbqb(ji)+(1.0-albqb(ji))
     &              *(zhgnew-hgbqb(ji))*rhog
        dmnbqb(ji) = dmnbqb(ji)
     &              -(1.0-albqb(ji))*(zhgnew-hgbqb(ji))*alphs*rhon
C    &              -(1.0-albqb(ji))*(zhgnew-hgbqb(ji))*rhon
Cfd0    dmgbqb(ji) = dmgbqb(ji)+(1.0-albqb(ji))
Cfd0 &              *(zhgnew-hgbqb(ji))*(rhog-alphs*rhon)
C    &              *(zhgnew-hgbqb(ji))*(rhog-rhon)
c
c  Actualize new snow and ice thickness.
c
        hnbqb(ji)  = zhnnew
        hgbqb(ji)  = zhgnew
110   continue
c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  11. Lateral ablation.                                               |
c-----------------------------------------------------------------------
c
      do 120 ji=kideb,kiut
        zab(ji)    = albqb(ji)
        zia        = 1.0-max(zero,sign(one,-hgbqb(ji)))
        zian       = 1.0-max(zero,sign(one,-hnbqb(ji)))
        albqb(ji)  = (1.0-zia)+zia*zab(ji)
        fscbqb(ji) = (1.0-zab(ji))*(1.0-thcmb(ji))*fstbqb(ji)
        qfvbqb(ji) = (zab(ji)*zfsup(ji)+(1.0-zia)*(1.0-zab(ji))*
     &               zfocea(ji))*ddtb
        qlbqb(ji)  = qlbqb(ji)+qfvbqb(ji)+(1.0-zia)*fscbqb(ji)*ddtb
        zq         = hnbqb(ji)*xln+hgbqb(ji)*xlg
        zqtot      = (1.0-albqb(ji))*zq
        ziq        = 1.0-max(zero,sign(one,qlbqb(ji)-zqtot))
        albqb(ji)  = ziq*(albqb(ji)+max(zero,qlbqb(ji)/max(zeps,zq)))+
     &               (1.0-ziq)
        fltbqb(ji) = ((1.0-zab(ji))*qstbqb(ji)-zqtot)/ddtb
c
c  Opening of leads: HAKKINEN & MELLOR.
c
        zalbqb(ji) = (albqb(ji)+max(zero,(-(zdhg(ji)+zdhb(ji))*hakspl*
     &           (1.0-zab(ji)))/max(zeps0,hgbqb(ji)
     &            +hnbqb(ji)*rhon/rhog)))
c
c  Opening of leads: HIbLER.
c
C       zalbqb(ji) = (albqb(ji)+max(zero,
C    &               ((-zdhg(ji)-zdhb(ji))*hibspl)/
C    &               max(zeps,hgbqb(ji)+hnbqb(ji)*rhon/rhog)))
c
c  Opening of leads: OLD.
c
C       zalbqb(ji) = albqb(ji)
c
        zalbqb(ji) = ziq*min(0.99*one,zalbqb(ji))+(1-ziq)
c
        tsb(ji)    = tsb(ji)+(1.0-ziq)*(tfub(ji)-tsb(ji))
        tbqb(ji,1) = tbqb(ji,1)+(1.0-ziq)*(tfub(ji)-tbqb(ji,1))
        tbqb(ji,2) = tbqb(ji,2)+(1.0-ziq)*(tfub(ji)-tbqb(ji,2))
        tbqb(ji,3) = tbqb(ji,3)+(1.0-ziq)*(tfub(ji)-tbqb(ji,3))
        dvlbqb(ji) = zia*(zab(ji)-albqb(ji))*hgbqb(ji)
        dmgbqb(ji) = dmgbqb(ji)+dvlbqb(ji)*rhog
        dvn        = zian*(zab(ji)-albqb(ji))*hnbqb(ji)
        dmnbqb(ji) = dmnbqb(ji)+dvn*rhon
        hnbqb(ji)  = ziq*hnbqb(ji)
        hnbqb(ji)  = ziq*hnbqb(ji)*(1.0-albqb(ji))
     &               /max(zeps,1.0-zalbqb(ji))
        hgbqb(ji)  = ziq*hgbqb(ji)*(1.0-albqb(ji))
     &               /max(zeps,1.0-zalbqb(ji))
        qstbqb(ji) = ziq*(1.0-zab(ji))*qstbqb(ji)
     &               /max(zeps,1.0-zalbqb(ji))
        albqb(ji)  = zalbqb(ji)
120   continue
c
      return
      end
c
