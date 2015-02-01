!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** File:     comdiag.h
! *** Contents: Common declarations for diagnostics of ECbilt
!-----------------------------------------------------------------------
! *** COMMON  /writectl/ixout,ioutdaily,ioutyearly,itel,minterv,
!                      meantype,meantot,meanyl,ifrendat,instcount
!     ixout:   output frequency for instantanous field in days
!     ioutdaily: output (1) or not (0) instantanous fields.
!     ioutyearly: output (1) or not (0) yearly mean fields.
!     itel:    counts the number of days in the output interval
!     minterv: counter used to compute monthly or seasonal mean.
!     meantype = 1 monthly mean
!     meantype = 2 seasonal mean.
!     meantot = 1 output total integration period monthly or seasonal
!               mean
!     meanyl  = 1 output yearly monthly or seasonal mean.
!     ifrendat: frequency of writing restart data in days.
!     instcount: counter used for output instantaneous fields.
!     irad =1 if calculation of radiative forcings
!-----------------------------------------------------------------------
      integer     ixout,ioutdaily,ioutyearly,itel,minterv,ntotdays, &
     &            meantype,meantot,meanyl,ifrendat,instcount,irad

      common /writectl/ixout,ioutdaily,ioutyearly,itel,minterv, &
     &                 meantype,meantot,meanyl,ifrendat,instcount,irad

      integer numtotvar,thirdd(4)
      real*8  fill_value,missing_value,newtotvar(80,20)
      character*60 nametotvar(80,6)
      common /netcdf/newtotvar,fill_value,missing_value,numtotvar, &
     &                  nametotvar,thirdd

      real*8  s1u200(nlat,nlon),s2u200(nlat,nlon), &
     &        s1u500(nlat,nlon),s2u500(nlat,nlon), &
     &        s1u800(nlat,nlon),s2u800(nlat,nlon), &
     &        s1omeg(nlat,nlon,nvl),s2omeg(nlat,nlon,nvl), &
     &        s1temp4g(nlat,nlon),s2temp4g(nlat,nlon), &
     &        s1temp2g(nlat,nlon),s2temp2g(nlat,nlon), &
     &        s1tstrat(nlat,nlon),s2tstrat(nlat,nlon), &
     &        s1tsurf(nlat,nlon),s2tsurf(nlat,nlon), &
     &        s1tempsg(nlat,nlon),s2tempsg(nlat,nlon), &
     &        s1t2m(nlat,nlon),s2t2m(nlat,nlon), &
     &        s1v200(nlat,nlon),s2v200(nlat,nlon), &
     &        s1v500(nlat,nlon),s2v500(nlat,nlon), &
     &        s1v800(nlat,nlon),s2v800(nlat,nlon), &
     &        s1psi(nlat,nlon,nvl),s2psi(nlat,nlon,nvl), &
     &        s1vhforg1(nlat,nlon),s2vhforg1(nlat,nlon), &
     &        s1vhforg2(nlat,nlon),s2vhforg2(nlat,nlon), &
     &        s1vforg1(nlat,nlon),s2vforg1(nlat,nlon), &
     &        s1pground(nlat,nlon),s2pground(nlat,nlon)

      real*8  s1vforg2(nlat,nlon),s2vforg2(nlat,nlon), &
     &        s1vforg3(nlat,nlon),s2vforg3(nlat,nlon), &
     &        s1udivg(nlat,nlon,nvl),s2udivg(nlat,nlon,nvl), &
     &        s1vdivg(nlat,nlon,nvl),s2vdivg(nlat,nlon,nvl), &
     &        s1dyrain(nlat,nlon),s2dyrain(nlat,nlon), &
     &        s1corain(nlat,nlon),s2corain(nlat,nlon), &
     &        s1snow(nlat,nlon),s2snow(nlat,nlon), &
     &        s1Hflux(nlat,nlon),s2Hflux(nlat,nlon), &
     &        s1Eflux(nlat,nlon),s2Eflux(nlat,nlon), &
     &        s1hesw(nlat,nlon),s2hesw(nlat,nlon), &
     &        s1hesws(nlat,nlon),s2hesws(nlat,nlon), &
     &        s1ulrad1(nlat,nlon),s2ulrad1(nlat,nlon), &
     &        s1nlrads(nlat,nlon),s2nlrads(nlat,nlon), &
     &        s1bmoisg(nlat,nlon),s2bmoisg(nlat,nlon), &
     &        s1rmoisgw3(nlat,nlon),s2rmoisgw3(nlat,nlon), &
     &        s1relhum(nlat,nlon),s2relhum(nlat,nlon), &
     &        s1torain(nlat,nlon),s2torain(nlat,nlon), &
     &        s1evap(nlat,nlon),s2evap(nlat,nlon), &
     &        s1eminp(nlat,nlon),s2eminp(nlat,nlon), &
     &        s1albes(nlat,nlon),s2albes(nlat,nlon), &
     &        s1albep(nlat,nlon),s2albep(nlat,nlon), &
     &        s1winstu1(nlat,nlon),s2winstu1(nlat,nlon), &
     &        s1winstv1(nlat,nlon),s2winstv1(nlat,nlon), &
     &        s1dsnow(nlat,nlon),s2dsnow(nlat,nlon), &
     &        s1hic(nlat,nlon),s2hic(nlat,nlon), &
     &        s1runofo(nlat,nlon),s2runofo(nlat,nlon), &
     &        s1runofl(nlat,nlon),s2runofl(nlat,nlon), &
     &        s1uv10(nlat,nlon),s2uv10(nlat,nlon), &
     &        s1chi(nlat,nlon,nvl),s2chi(nlat,nlon,nvl), &
     &        s1gh(nlat,nlon,nvl),s2gh(nlat,nlon,nvl), &
     &        s1qgpv(nlat,nlon,nvl),s2qgpv(nlat,nlon,nvl), &
     &        s1cdragw(nlat,nlon),s2cdragw(nlat,nlon), &
     &        s1cdragv(nlat,nlon),s2cdragv(nlat,nlon), &
     &        s1tcc(nlat,nlon),s2tcc(nlat,nlon), &
     &        s1dt1(nlat,nlon,nvl),s2dt1(nlat,nlon,nvl), &
     &        s1dt2(nlat,nlon,nvl),s2dt2(nlat,nlon,nvl), &
     &        s1du1(nlat,nlon,nvl),s2du1(nlat,nlon,nvl), &
     &        s1du2(nlat,nlon,nvl),s2du2(nlat,nlon,nvl)

      real*8   sxu200(nlat,nlon,12),syu200(nlat,nlon,12), &
     &        sxu500(nlat,nlon,12),syu500(nlat,nlon,12), &
     &        sxu800(nlat,nlon,12),syu800(nlat,nlon,12), &
     &        sxomeg3(nlat,nlon,12),syomeg3(nlat,nlon,12), &
     &        sxomeg2(nlat,nlon,12),syomeg2(nlat,nlon,12), &
     &        sxomeg1(nlat,nlon,12),syomeg1(nlat,nlon,12), &
     &        sxtemp4g(nlat,nlon,12),sytemp4g(nlat,nlon,12), &
     &        sxtemp2g(nlat,nlon,12),sytemp2g(nlat,nlon,12), &
     &        sxtstrat(nlat,nlon,12),sytstrat(nlat,nlon,12), &
     &        sxtsurf(nlat,nlon,12),sytsurf(nlat,nlon,12), &
     &        sxtempsg(nlat,nlon,12),sytempsg(nlat,nlon,12), &
     &        sxt2m(nlat,nlon,12),syt2m(nlat,nlon,12), &
     &        sxv200(nlat,nlon,12),syv200(nlat,nlon,12), &
     &        sxv500(nlat,nlon,12),syv500(nlat,nlon,12), &
     &        sxv800(nlat,nlon,12),syv800(nlat,nlon,12), &
     &        sxgrpsi3(nlat,nlon,12),sygrpsi3(nlat,nlon,12), &
     &        sxgrpsi2(nlat,nlon,12),sygrpsi2(nlat,nlon,12), &
     &        sxgrpsi1(nlat,nlon,12),sygrpsi1(nlat,nlon,12), &
     &        sxpground(nlat,nlon,12),sypground(nlat,nlon,12)

      real*8 sxvhforg1(nlat,nlon,12),syvhforg1(nlat,nlon,12), &
     &        sxvhforg2(nlat,nlon,12),syvhforg2(nlat,nlon,12), &
     &        sxvforg1(nlat,nlon,12),syvforg1(nlat,nlon,12), &
     &        sxvforg2(nlat,nlon,12),syvforg2(nlat,nlon,12), &
     &        sxvforg3(nlat,nlon,12),syvforg3(nlat,nlon,12), &
     &        sxudivg3(nlat,nlon,12),syudivg3(nlat,nlon,12), &
     &        sxudivg2(nlat,nlon,12),syudivg2(nlat,nlon,12), &
     &        sxudivg1(nlat,nlon,12),syudivg1(nlat,nlon,12), &
     &        sxvdivg3(nlat,nlon,12),syvdivg3(nlat,nlon,12), &
     &        sxvdivg2(nlat,nlon,12),syvdivg2(nlat,nlon,12), &
     &        sxvdivg1(nlat,nlon,12),syvdivg1(nlat,nlon,12), &
     &        sxdyrain(nlat,nlon,12),sydyrain(nlat,nlon,12), &
     &        sxcorain(nlat,nlon,12),sycorain(nlat,nlon,12), &
     &        sxhflux(nlat,nlon,12),syhflux(nlat,nlon,12), &
     &        sxeflux(nlat,nlon,12),syeflux(nlat,nlon,12), &
     &        sxhesw(nlat,nlon,12),syhesw(nlat,nlon,12), &
     &        sxhesws(nlat,nlon,12),syhesws(nlat,nlon,12), &
     &        sxulrad1(nlat,nlon,12),syulrad1(nlat,nlon,12), &
     &        sxnlrads(nlat,nlon,12),synlrads(nlat,nlon,12), &
     &        sxbmoisg(nlat,nlon,12),sybmoisg(nlat,nlon,12), &
     &        sxrmoisgw3(nlat,nlon,12),syrmoisgw3(nlat,nlon,12), &
     &        sxrelhum(nlat,nlon,12),syrelhum(nlat,nlon,12), &
     &        sxtorain(nlat,nlon,12),sytorain(nlat,nlon,12), &
     &        sxsnow(nlat,nlon,12),sysnow(nlat,nlon,12), &
     &        sxevap(nlat,nlon,12),syevap(nlat,nlon,12), &
     &        sxeminp(nlat,nlon,12),syeminp(nlat,nlon,12), &
     &        sxalbes(nlat,nlon,12),syalbes(nlat,nlon,12), &
     &        sxalbep(nlat,nlon,12),syalbep(nlat,nlon,12), &
     &        sxwinstu1(nlat,nlon,12),sywinstu1(nlat,nlon,12), &
     &        sxwinstv1(nlat,nlon,12),sywinstv1(nlat,nlon,12), &
     &        sxdsnow(nlat,nlon,12),sydsnow(nlat,nlon,12), &
     &        sxhic(nlat,nlon,12),syhic(nlat,nlon,12), &
     &        sxrunofo(nlat,nlon,12),syrunofo(nlat,nlon,12), &
     &        sxrunofl(nlat,nlon,12),syrunofl(nlat,nlon,12), &
     &        sxuv10(nlat,nlon,12),syuv10(nlat,nlon,12), &
     &        sxchi3(nlat,nlon,12),sychi3(nlat,nlon,12), &
     &        sxchi2(nlat,nlon,12),sychi2(nlat,nlon,12), &
     &        sxchi1(nlat,nlon,12),sychi1(nlat,nlon,12), &
     &        sxgh3(nlat,nlon,12),sygh3(nlat,nlon,12), &
     &        sxgh2(nlat,nlon,12),sygh2(nlat,nlon,12), &
     &        sxgh1(nlat,nlon,12),sygh1(nlat,nlon,12), &
     &        sxqgpv3(nlat,nlon,12),syqgpv3(nlat,nlon,12), &
     &        sxqgpv2(nlat,nlon,12),syqgpv2(nlat,nlon,12), &
     &        sxqgpv1(nlat,nlon,12),syqgpv1(nlat,nlon,12), &
     &        sxcdragw(nlat,nlon,12),sycdragw(nlat,nlon,12), &
     &        sxcdragv(nlat,nlon,12),sycdragv(nlat,nlon,12), &
     &        sxtcc(nlat,nlon,12),sytcc(nlat,nlon,12), &
     &        sxdt13(nlat,nlon,12),sydt13(nlat,nlon,12), &
     &        sxdt12(nlat,nlon,12),sydt12(nlat,nlon,12), &
     &        sxdt11(nlat,nlon,12),sydt11(nlat,nlon,12), &
     &        sxdt23(nlat,nlon,12),sydt23(nlat,nlon,12), &
     &        sxdt22(nlat,nlon,12),sydt22(nlat,nlon,12), &
     &        sxdt21(nlat,nlon,12),sydt21(nlat,nlon,12), &
     &        sxdu13(nlat,nlon,12),sydu13(nlat,nlon,12), &
     &        sxdu12(nlat,nlon,12),sydu12(nlat,nlon,12), &
     &        sxdu11(nlat,nlon,12),sydu11(nlat,nlon,12), &
     &        sxdu23(nlat,nlon,12),sydu23(nlat,nlon,12), &
     &        sxdu22(nlat,nlon,12),sydu22(nlat,nlon,12), &
     &        sxdu21(nlat,nlon,12),sydu21(nlat,nlon,12)

      common /outxx/ s1u200,s2u200,s1u500,s2u500,s1u800,s2u800, &
     &        s1omeg,s2omeg, &
     &        s1temp4g,s2temp4g,s1temp2g,s2temp2g,s1tempsg,s2tempsg, &
     &        s1tsurf,s2tsurf,s1tstrat,s2tstrat,s1t2m,s2t2m, &
     &        s1v200,s2v200,s1v500,s2v500,s1v800,s2v800, &
     &        s1psi,s2psi, &
     &        s1vhforg1,s2vhforg1,s1vhforg2,s2vhforg2, &
     &        s1vforg1,s2vforg1,s1vforg2,s2vforg2,s1vforg3,s2vforg3, &
     &        s1udivg,s2udivg,s1vdivg,s2vdivg, &
     &        s1dyrain,s2dyrain,s1corain,s2corain, &
     &        s1snow,s2snow, &
     &        s1Hflux,s2Hflux,s1Eflux,s2Eflux, &
     &        s1hesw,s2hesw,s1hesws,s2hesws, &
     &        s1ulrad1,s2ulrad1, &
     &        s1nlrads,s2nlrads, &
     &        s1bmoisg,s2bmoisg, &
     &        s1rmoisgw3,s2rmoisgw3,s1relhum,s2relhum, &
     &        s1torain,s2torain, &
     &        s1evap,s2evap, &
     &        s1eminp,s2eminp, &
     &        s1albes,s2albes, &
     &        s1albep,s2albep, &
     &        s1winstu1,s2winstu1, &
     &        s1winstv1,s2winstv1, &
     &        s1dsnow,s2dsnow,s1hic,s2hic, &
     &        s1runofo,s2runofo, &
     &        s1runofl,s2runofl, &
     &        s1uv10,s2uv10,s1gh,s2gh,s1qgpv,s2qgpv, &
     &        s1chi,s2chi,s1pground,s2pground, &
     &        s1cdragw,s2cdragw,s1cdragv,s2cdragv, &
     &        s1tcc,s2tcc,s1dt1,s2dt1,s1dt2,s2dt2,s1du1,s2du1, &
     &        s1du2,s2du2

      common /outyy/sxu200,syu200,sxu500,syu500,sxu800,syu800, &
     &        sxomeg3,syomeg3,sxomeg2,syomeg2,sxomeg1,syomeg1, &
     &        sxtemp4g,sytemp4g,sxtemp2g,sytemp2g,sxtsurf,sytsurf, &
     &        sxtempsg,sytempsg,sxtstrat,sytstrat,sxt2m,syt2m, &
     &        sxv200,syv200,sxv500,syv500,sxv800,syv800, &
     &        sxgrpsi3,sygrpsi3,sxgrpsi2,sygrpsi2,sxgrpsi1,sygrpsi1, &
     &        sxvhforg1,syvhforg1,sxvhforg2,syvhforg2, &
     &        sxvforg1,syvforg1,sxvforg2,syvforg2,sxvforg3,syvforg3, &
     &        sxudivg3,syudivg3,sxudivg2,syudivg2,sxudivg1,syudivg1, &
     &        sxvdivg3,syvdivg3,sxvdivg2,syvdivg2,sxvdivg1,syvdivg1, &
     &        sxdyrain,sydyrain,sxcorain,sycorain, &
     &        sxhflux,syhflux,sxeflux,syeflux, &
     &        sxhesw,syhesw,sxhesws,syhesws, &
     &        sxulrad1,syulrad1,sxnlrads,synlrads, &
     &        sxbmoisg,sybmoisg, &
     &        sxrmoisgw3,syrmoisgw3, &
     &        sxrelhum,syrelhum, &
     &        sxtorain,sytorain, &
     &        sxsnow,sysnow, &
     &        sxevap,syevap, &
     &        sxeminp,syeminp, &
     &        sxalbes,syalbes, &
     &        sxalbep,syalbep, &
     &        sxwinstu1,sywinstu1, &
     &        sxwinstv1,sywinstv1, &
     &        sxdsnow,sydsnow,sxhic,syhic, &
     &        sxrunofo,syrunofo, &
     &        sxrunofl,syrunofl, &
     &        sxuv10,syuv10, &
     &        sxchi3,sychi3,sxchi2,sychi2,sxchi1,sychi1, &
     &        sxgh3,sygh3,sxgh2,sygh2,sxgh1,sygh1, &
     &        sxqgpv3,syqgpv3,sxqgpv2,syqgpv2,sxqgpv1,syqgpv1, &
     &        sxpground,sypground, &
     &        sxcdragw,sycdragw,sxcdragv,sycdragv, &
     &        sxtcc,sytcc, &
     &        sxdt13,sydt13,sxdt12,sydt12,sxdt11,sydt11, &
     &        sxdt23,sydt23,sxdt22,sydt22,sxdt21,sydt21, &
     &        sxdu13,sydu13,sxdu12,sydu12,sxdu11,sydu11, &
     &        sxdu23,sydu23,sxdu22,sydu22,sxdu21,sydu21
