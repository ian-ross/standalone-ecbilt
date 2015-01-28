












      subroutine thersf(ntnow,ntrmax)


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c   This routine prepares the computation of ice thermodynamics. It also
c   computes surface heat and salt flux at the surface of the ocean
c  modif : 24/12/99

c---
c Ccpl [Ccp0] => ligne specifique a la version avec [sans] couplage .
c---
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'
      include 'thermo.com'
      include 'dynami.com'
Ccfc  include 'trace.com'

      real*8 heaticbism(imax,jmax)
      real*8 heaticbs,heaticbs2
      real*8 phistots,phistotn
      logical flgveg,flgicb,flgisma,flgismg
      real*8 fluxmean
      real   qref,qclim
      real*8  day
      real   dzsdts,DUM
      integer ntstep,nstpyear,nocstpyear
      integer nyears,ndays,irunlabel,irunlabeld,nwrskip_days,end_year,end_day
      common/icbism/ heaticbism
      common /ec_coupl/flgveg,flgicb,flgisma,flgismg
      COMMON/CLIO2A/Qref,Qclim
      COMMON/CLIO2A2/fluxmean
      integer     iyear,imonth,iday,iseason,ntotday,nbclins,nbtrops
      common /ec_timectl/ day,ntstep,nstpyear,nocstpyear,iyear,imonth,
     *                 iday,iseason,ntotday,nbclins,nbtrops
      common /ec_timepar/ nyears,ndays,irunlabel,irunlabeld,iatm,ilan,iice,iobtrop,iobclin,
     *                 nwrskip,nwrskip_days,end_year,end_day

c fdtcn     Transit variable used to compute the vertical heat flux
c           at the ocean surface in ice covered regions (local in thersf)
c
      dimension ain(imax,jmax),zinda(imax,jmax),ifvt(imax,jmax)
      dimension qdtcn(imax,jmax),qlbsbq(imax,jmax),zfwat(imax,jmax)
      dimension zhgbqp(imax,jmax),zhnpbq(imax,jmax)
      dimension fdtcn(imax,jmax)
      common/lowat/ zfwat

c     dimension rrro(imax,jmax),dd1o(imax,jmax),dd2o(imax,jmax)
c

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1) Initialization of diagnostic variables.                          |
c-----------------------------------------------------------------------

C       write(106,*)'albq,hgbq',albq(imax-2,5),hgbq(imax-2,5)
c
Cage    dtj=dts(ks2)/(86400.*365.)
        do 5 j=js1,js2
           do  3 i=is1(j),is2(j)
              dvosbq(i,j) = 0.0
              dvobbq(i,j) = 0.0
              dvolbq(i,j) = 0.0
              dvonbq(i,j) = 0.0
 3         continue
 5       continue
c

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2) Numerical parameters.                                            |
c-----------------------------------------------------------------------
c
      zeps0 = 1.0e-16
      zeps1 = 1.0e-20
      zeps2 = 1.0e-04
      ustmax= 2.0d-02
C     ustmax= 4.0d-02
      ustmin= 5.0d-03
c

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3) Initialization of some arrays.                                   |
c-----------------------------------------------------------------------
c
      if (ntnow.eq.1) then
        do j=js1,js2
          do i=is1(j),is2(j)
            alct(i,j) = 0.0
          enddo
        enddo
      endif
      do 40 j=js1,js2
        do 30 i=is1(j),is2(j)
          fstrbq(i,j) = 0.0
          fscmbq(i,j) = 0.0
          ffltbq(i,j) = 0.0
          qfvbq(i,j)  = 0.0
          hgcol(i,j)  = 0.0
          hnbq(i,j)   = hnbq(i,j)*max(zero,
     &                  sign(one,hnbq(i,j)-zeps2))
          dmnbq(i,j)  = 0.0
          dmgbq(i,j)  = 0.0
          dmgwi(i,j)  = 0.0
Cage      ageg(i,j)  = (1.0-zindg)*ageg(i,j)+zindg*dtj
30      continue
40    continue
c
c
c--division by the number of thermodynamic steps
c--hnplbq and fwat are in kg per day and per m2
      raptime=ddtb/86400.
      do j=js1,js2
        do i=is1(j),is2(j)
           zhnpbq(i,j) = raptime * hnplbq(i,j)
           zfwat(i,j)  = raptime * fwat(i,j)
        enddo
      enddo


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4) Treatment of particular cases.                                   |
c-----------------------------------------------------------------------
c
        do 210 j=js1,js2
          do 200 i=is1(j),is2(j)
            zindb = tms(i,j,ks2)*(1.0-max(zero,
     &              sign(one,-(hnbq(i,j)+hgbq(i,j)))))
c
c-- 4.1. Snow is transformed into ice if the original
c        ice cover disappears.
c-----------------------------------------------------
c
            zindg      = tms(i,j,ks2)*max(zero,
     &                    sign(one,-hgbq(i,j)))
            hgbq(i,j)  = hgbq(i,j)+zindg*rhon*hnbq(i,j)/rho0
            hnbq(i,j)  = (1.0-zindg)*hnbq(i,j)+zindg
     &                    *hgbq(i,j)*(rho0-rhog)/rhon
c           dmgbq(i,j) = zindg*(1.0-albq(i,j))*rhog*hgbq(i,j)
            dmgwi(i,j) = zindg*(1.0-albq(i,j))*rhog*hgbq(i,j)
c
c-- 4.2. The lead fraction, albq, must be little than
c        or equal equal to amax (ice ridging).
c-----------------------------------------------------
c
            za    =   zindb*min(one,(1.0-albq(i,j))*uscomi)
            hnbq(i,j)   = hnbq(i,j)*za
            hgbq(i,j)   = hgbq(i,j)*za
            qstobq(i,j) = qstobq(i,j)*za
            albq(i,j)   = 1.0-zindb*(1.0-albq(i,j))/max(za,zeps1)
Cage        agen(i,j)   = agen(i,j)*vnbq0(i,j)/
Cage &                    max(dvn+vnbq0(i,j),zeps0)+dtj
Cage        ageg(i,j)   = ageg(i,j)*vgbq0(i,j)/
Cage &                    max(dvg+vgbq0(i,j),zeps0)+dtj
Cage        agen(i,j)   = (1.0-zindn)*agen(i,j)
Cage        ageg(i,j)   = (1.0-zindg)*ageg(i,j)
c
c-- 4.3. The in situ ice thickness, hgbq, must be equal to
c        or greater than hglim.
c----------------------------------------------------------
c
            zh          = max(one,zindb*hglim/max(hgbq(i,j),zeps1))
            hnbq(i,j)   = hnbq(i,j)*zh
            hgbq(i,j)   = hgbq(i,j)*zh
            qstobq(i,j) = qstobq(i,j)*zh
            albq(i,j)   = (albq(i,j)+(zh-1.0))/zh
200       continue
210     continue
c
c

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5. Thermodynamics of sea ice.                                       |
c-----------------------------------------------------------------------
c
c--5.1. Partial computation of forcing for the thermo-
c       dynamic sea ice model.
c------------------------------------------------------
c
      do 230 j=js1,js2
        do 220 i=is1(j),is2(j)
c
          zindb      = tms(i,j,ks2)*(one-max(zero
     &                 ,sign(one,-(hnbq(i,j)+hgbq(i,j)))))
          zinda(i,j) = 1.0-max(zero,sign(one,-(1.0-albq(i,j))))
          ain(i,j)   = albq(i,j)
c
c--Solar irradiance transmission at the mixed layer
c--bottom.
c
          thcm(i,j) = 1.-reslum(i,j,ks2)
     &                       *dz(ks2)*(rho0*cpo)
C         thcm(i,j) = 0.0
c
c--Calculate fdtcn and qdtcn. Limitation of mixed
c--layer temperature.
c
CC        fdtcn(i,j) = .5*zindb*rho0*cpo*dz(ks2)*
CC   &                    (scal(i,j,ks2,1)-tfu(i,j))/dts(ks2)
          ustsg = max(min(sqrt(ust2s(i,j)),ustmax),ustmin)
          fdtcn(i,j) = zindb*rho0*cpo*0.006
     &             *ustsg*(scal(i,j,ks2,1)-tfu(i,j))
          qdtcn(i,j)  = zindb*fdtcn(i,j)*albq(i,j)*ddtb
c
c-- Snow accumulation.
c
          zhnpbq(i,j) = zhnpbq(i,j)*rhoesn
c
c-- Partial computation of the lead energy budget and
c-- determination of qcmbq.
c
          zfontn     = zhnpbq(i,j)*xln/ddtb
C         if (icoupl .eq. 0) then
C          zfnsol     = zemise*ratbqo(i,j)-zemise*stefan*
C    &                 scal(i,j,ks2,1)*scal(i,j,ks2,1)*
C    &                 scal(i,j,ks2,1)*scal(i,j,ks2,1)
C    &                 -fcscn(i,j)-flecn(i,j)
C         else
C          zfnsol     = ratbqo(i,j)-zemise*stefan*
C    &                 scal(i,j,ks2,1)*scal(i,j,ks2,1)*
C    &                 scal(i,j,ks2,1)*scal(i,j,ks2,1)
C    &                 -fcscn(i,j)-flecn(i,j)
C         endif
          zfnsol     = ratbqo(i,j)-fcscn(i,j)-flecn(i,j)
          qlbq(i,j)  = tms(i,j,ks2)*(fsolcn(i,j)*(1.0-thcm(i,j))
     &                 +zfnsol+fdtcn(i,j)
     &                 -zfontn+(1.0-zindb)*fsbbq(i,j))*albq(i,j)*ddtb
          fntlat      = 1.0-max(zero,sign(one,-qlbq(i,j)))
          pareff      = 1.0+(parlat-1.0)*zinda(i,j)*fntlat
          qlbsbq(i,j) = qlbq(i,j)*(1.0-pareff)
     &                  /max(((one-albq(i,j))*ddtb),zeps0)
          qlbq(i,j)   = pareff*qlbq(i,j)
          qdtcn(i,j)  = pareff*qdtcn(i,j)
          qcmbq(i,j)  = rho0*cpo*dz(ks2)*(tfu(i,j)-scal(i,j,ks2,1))
     &                 *(1-zinda(i,j))
C    &                 *.5*(1-zinda(i,j))
c
c--Calculate oceanic heat flux.
c
          fbbq(i,j)  = zindb*(fsbbq(i,j)/max((one-albq(i,j))
     &                                      ,zeps1)+fdtcn(i,j))
c
c--Store ice volume per unit area for computation of
c--daily thermodynamic ice production.
c--This variable is only needed for output.
c
          zhgbqp(i,j) = hgbq(i,j)*(1.0-albq(i,j))
c
220     continue
230   continue

c
c-- 5.2. Select icy points and fulfill arrays for the
c        vectorial grid.
c------------------------------------------------------
c
      nbpb = 0
      do 250 j=js1,js2
        do 240 i=is1(j),is2(j)
          if (albq(i,j).ge.1.0) go to 240
          nbpb      = nbpb+1
          npb(nbpb) = (j-1)*imax+i
240     continue
250   continue
c
c--If there is no ice, do nothing.
c
      if (nbpb.eq.0) go to 280
c
      call gather(nbpb,albqb,albq,npb)
      call gather(nbpb,hgbqb,hgbq,npb)
      call gather(nbpb,hnbqb,hnbq,npb)
      call gather(nbpb,hnpbqb,zhnpbq,npb)
      call gather(nbpb,fsolgb,fsolg,npb)
      call gather(nbpb,fbbqb,fbbq,npb)
      call gather(nbpb,thcmb,thcm,npb)
      call gather(nbpb,qlbqb,qlbq,npb)
      call gather(nbpb,qcmbqb,qcmbq,npb)
      call gather(nbpb,qstbqb,qstobq,npb)
      call gather(nbpb,tfub,tfu,npb)
      call gather(nbpb,tsb,ts,npb)
      call gather(nbpb,dmgbqb,dmgbq,npb)
      call gather(nbpb,dmgwib,dmgwi,npb)
      call gather(nbpb,psbqb,psbq,npb)
      call gather(nbpb,tabqb,tabq,npb)
      call gather(nbpb,qabqb,qabq,npb)
      call gather(nbpb,vabqb,vabq,npb)
      call gather(nbpb,ratbqb,ratbqg,npb)
      call gather(nbpb,qlbbqb,qlbsbq,npb)
      call gather(nbpb,cldqb,cloud,npb)
      do 260 k=1,nkb0
        call gather(nbpb,tbqb(1,k),tbq(1,1,k),npb)
260   continue
Ccp2  call gather(nbpb,fderb,fder,npb)
c
c--5.3. Call to ice growth routine.
c-----------------------------------
c
      call fontbc(1,nbpb)
c
c--5.4. Back to the geographic grid.
c-----------------------------------
c
      call scater(nbpb,albq,npb,albqb)
      call scater(nbpb,hnbq,npb,hnbqb)
      call scater(nbpb,hgbq,npb,hgbqb)
      call scater(nbpb,fscmbq,npb,fscbqb)
      call scater(nbpb,ffltbq,npb,fltbqb)
      call scater(nbpb,fstrbq,npb,fstbqb)
      call scater(nbpb,qlbq,npb,qlbqb)
      call scater(nbpb,qfvbq,npb,qfvbqb)
      call scater(nbpb,qstobq,npb,qstbqb)
      call scater(nbpb,ts,npb,tsb)
      call scater(nbpb,dmgbq,npb,dmgbqb)
      call scater(nbpb,dmgwi,npb,dmgwib)
      call scater(nbpb,dmnbq,npb,dmnbqb)
      do 270 k=1,nkb0
        call scater(nbpb,tbq(1,1,k),npb,tbqb(1,k))
270   continue
      call scater(nbpb,dvosbq,npb,dvsbqb)
      call scater(nbpb,dvobbq,npb,dvbbqb)
      call scater(nbpb,dvolbq,npb,dvlbqb)
      call scater(nbpb,dvonbq,npb,dvnbqb)
      call scater(nbpb,firg,npb,fratsb)
      call scater(nbpb,fcsg,npb,fcsb)
      call scater(nbpb,fleg,npb,fleb)
c
280   continue
c
c--5.5. Up-date sea ice thickness.
c---------------------------------
c
c
      do 300 j=js1,js2
        do 290 i=is1(j),is2(j)
          ifvt(i,j) = zinda(i,j)*max(zero,sign(one,-hgbq(i,j)))
          hgbq(i,j) = hgbq(i,j)*(one-max(zero,
     &                          sign(one,-(1.0-albq(i,j)))))
290     continue
300   continue
c
c--Tricky trick: add 2 to albq in the Southern Hemisphere.
c
      do 307 j=js1,jeq-1
         do 305 i=is1(j),is2(j)
            albq(i,j) = albq(i,j)+2.0
 305     continue
 307  continue
c
c-- 5.6. Select points for lateral accretion and fulfill
c--      arrays for the vectorial grid.
c--------------------------------------------------------
c
c collection thickness of frazil ice (with option Cvhg)
c
c hgcol  : collection thickness of consolidated new ice
c hg0    : frazil ice thickness at the leading edge of
c          the lead
c alfagr : fraction of frazil ice in grease ice (the rest
c          is water)
c rhoo   : seawater density
c rhofr  : frazil ice density
c rhogr  : grease ice density
c vafrmn : minimum surface wind velocity for computation
c          of hgcol (in m/s)
c gamafr : ratio between grease ice and wind speeds
c
      alfagr = 0.35
      rhoo   = 1029.0
      rhofr  = 950.0
      rhogr  = rhoo+alfagr*(rhofr-rhoo)
      gred   = gpes*(rhoo-rhogr)/rhoo
      hg0    = 0.1
      fac1hg = sqrt(2.0*hg0/gred)*2.0/acos(-1.0)
      fac2hg = 0.25/gred
      gamafr = 0.06
      do j=js1,js2
        jp1   = j+1
        do i=is1(j),is2(j)
          if (qcmbq(i,j)-qlbq(i,j).gt.0.0) then
            ip1         = i+1
            tenagm      = sqrt( tenagx(i,j)*tenagx(i,j)
     &                         +tenagy(i,j)*tenagy(i,j))
C approx coupled
            vabq(i,j)   = tenagm/0.041
C
            dumfac      = gamafr*vabq(i,j)/max(1.0d-06,tenagm)
            vfrx        = dumfac*tenagx(i,j)
            vfry        = dumfac*tenagy(i,j)
            vgx         = (ug(i,j)+ug(ip1,j)+ug(i,jp1)+ug(ip1,jp1))*0.25
            vgy         = (vg(i,j)+vg(ip1,j)+vg(i,jp1)+vg(ip1,jp1))*0.25
            vrel        = sqrt((vfrx-vgx)**2+(vfry-vgy)**2)
            hgcol(i,j)  = hg0+(fac1hg+fac2hg*vrel)*vrel
          endif
        enddo
      enddo
C     write(99,*) 'collec',int(vabq(50,2)),hgcol(50,2),
C    &             int(vabq(104,59)),hgcol(104,59)
c
      nbpac=0
      do 320 j=js1,js2
        do 310 i=is1(j),is2(j)
          if (qcmbq(i,j)-qlbq(i,j).le.0.0) go to 310
          nbpac       = nbpac+1
          npac(nbpac) = (j-1)*imax+i
310     continue
320   continue
c
c--If nbpac = 0, do nothing.
c
      if (nbpac.eq.0) go to 350
c
      call gather(nbpac,albqb,albq,npac)
      call gather(nbpac,hnbqb,hnbq,npac)
      call gather(nbpac,hgbqb,hgbq,npac)
      call gather(nbpac,qlbqb,qlbq,npac)
      call gather(nbpac,qcmbqb,qcmbq,npac)
      call gather(nbpac,qstbqb,qstobq,npac)
      call gather(nbpac,dmgbqb,dmgbq,npac)
      call gather(nbpac,dmgwib,dmgwi,npac)
      call gather(nbpac,tfub,tfu,npac)
      call gather(nbpac,tsb,ts,npac)
      call gather(nbpac,hgcolb,hgcol,npac)
      do 330 k=1,nkb0
        call gather(nbpac,tbqb(1,k),tbq(1,1,k),npac)
330   continue
      call gather(nbpac,dvlbqb,dvolbq,npac)
c
c-- 5.7. Call lateral accretion routine.
c----------------------------------------
c
      call acrlbq(1,nbpac)
c
c-- 5.8. Back to the geographic grid.
c------------------------------------
c
      call scater(nbpac,albq,npac,albqb)
      call scater(nbpac,hnbq,npac,hnbqb)
      call scater(nbpac,hgbq,npac,hgbqb)
      call scater(nbpac,qstobq,npac,qstbqb)
      call scater(nbpac,ts,npac,tsb)
      call scater(nbpac,dmgbq,npac,dmgbqb)
      call scater(nbpac,dmgwi,npac,dmgwib)
      do 340 k=1,nkb0
        call scater(nbpac,tbq(1,1,k),npac,tbqb(1,k))
340   continue
      call scater(nbpac,dvolbq,npac,dvlbqb)
c
350   continue
c
      do 370 j=js1,js2
        do 360 i=is1(j),is2(j)
c
c--Tricky trick: recover albq values between
c--0 and 1 in the Southern Hemisphere.
c
          albq(i,j) = min(albq(i,j),abs(albq(i,j)-2.0))
c
c rate of lead creation/destruction due to thermodynamics.
c
          alct(i,j) = alct(i,j)+(albq(i,j)-ain(i,j))/ddtb
c
c-- 5.9. Daily thermodynamic ice production.
c--------------------------------------------
c
          hgbqp(i,j) = hgbq(i,j)*(1.0-albq(i,j))-zhgbqp(i,j)
     &                 +hgbqp(i,j)
360     continue
370   continue
c
c-- 5.10. Compute ages.
c-----------------------
c
Cage   do j=1,jmax
Cage    do i=1,imax
Cage      zindn      = max(zero,sign(one,-hnbq(i,j)))
Cage      zindg      = max(zero,sign(one,-hgbq(i,j)))
Cage      dvn        = max(zero,hnbq(i,j)*(1.0-albq(i,j))-vnbq0(i,j))
Cage      dvg        = max(zero,hgbq(i,j)*(1.0-albq(i,j))-vgbq0(i,j))
Cage      agen(i,j)  = agen(i,j)*vnbq0(i,j)/
Cage &                 max(dvn+vnbq0(i,j),zeps0)+dtj
Cage      ageg(i,j)  = ageg(i,j)*vgbq0(i,j)/
Cage &                 max(dvg+vgbq0(i,j),zeps0)+dtj
Cage      agen(i,j)  = (1.0-zindn)*agen(i,j)
Cage      ageg(i,j)  = (1.0-zindg)*ageg(i,j)
Cage    enddo
Cage  enddo
c
c
c-- 5.10. Raccord cyclique
c--------------------------
c
      if (ltest.ge.1) then
c--raccord cyclique pour hgbq,albq,hnbq,ts,tbq,firg,fcsg,fleg,
c                     fsbbq,fdtcn,qstobq,scal:
      call raccord(hgbq(1,1),zero,1,8)
      call raccord(albq(1,1),zero,1,8)
      call raccord(alct(1,1),zero,1,8)
      call raccord(hnbq(1,1),zero,1,8)
      call raccord(ts(1,1),zero,1,8)
      call raccord(tbq(1,1,1),zero,nkb0,8)
      call raccord(firg(1,1),zero,1,8)
      call raccord(fcsg(1,1),zero,1,8)
      call raccord(fleg(1,1),zero,1,8)
      call raccord(fsbbq(1,1),zero,1,8)
      call raccord(qstobq(1,1),zero,1,8)
Cage  call raccord(agen(1,1),zero,1,8)
Cage  call raccord(ageg(1,1),zero,1,8)
      call raccord(scal(1,1,ks2,1),zero,1,8)
      endif

c
c-- 5.11. Sea ice/ocean interface.
c---------------------------------
C
c
c--CALCULATE fcm1, fcm2, fwat, AND fsbbq.
c
      do ns=1,nsmax
        if (scpme(ns).eq.spvr) then
          do j=js1,js2
            do i=is1(j),is2(j)
              scs(i,j,ns) = scal(i,j,ks2,ns)
            enddo
          enddo
        else
          do j=js1,js2
            do i=is1(j),is2(j)
              scs(i,j,ns) = scpme(ns)
            enddo
          enddo
        endif
      enddo

C-AM
C In order to correct for lack of consistency between salinity restoring
C fluxes defined in different routines restoring fluxes are set in this
C routine and transmitted through common array "frapsav"
C
C--borne inf. pour salinite (dans terme rappel sur la sal. pour FW Flx)
      salinf=5.0d0
C
      dzsdts=dz(ks2)/dts(ks2)
C-AM

      do 390 j=js1,js2
        do 380 i=is1(j),is2(j)
          salflux    =  34.7
C         salflux    =  scal(i,j,ks2,2)
          ia         =  1.0-max(zero,sign(one,-(1.0-albq(i,j))))
          iflt       =  zinda(i,j)*(1-ia)*(1-ifvt(i,j))
          ial        =  ifvt(i,j)*ia+(1-ifvt(i,j))*
     &                  (one-max(zero,sign(one,
     &                    albq(i,j)-ain(i,j))))
          fcm1(i,j)  =  ain(i,j)*fsolcn(i,j)+(1.-ain(i,j))*fstrbq(i,j)
          iadv       =  (1-ia)*zinda(i,j)
          zqfv       =  qfvbq(i,j)+qdtcn(i,j)
          fcm2(i,j)  = -fcm1(i,j)*(1-thcm(i,j))+
     &                  iflt*((1-iadv)*ffltbq(i,j)+fscmbq(i,j))+
     &                  (1-ia*(1-ial))*(ial*qcmbq(i,j)+
     &                  (1-ial)*(qlbq(i,j)-iadv*zqfv))/ddtb
CC        zdtm       =  (zqfv+ffltbq(i,j)*ddtb*iflt)/
CC   &                  (rho0*cpo*dz(ks2))
CC        scal(i,j,ks2,1) =  scal(i,j,ks2,1)+iadv*zdtm
          zfdtm      = (zqfv/ddtb+ffltbq(i,j)*iflt)*iadv
          fdtcn(i,j) = fdtcn(i,j)-zfdtm
          fsbbq(i,j) =  (1-(ifvt(i,j)+iflt))*fscmbq(i,j)
Cfd0      prs        =  (zfwat(i,j)-zhnpbq(i,j)*
Cfd0 &                  (1.-ain(i,j))*rhon)*salflux
          fons       =  dmgbq(i,j)*(salflux-sglace)
Cfd0 &                 -dmgwi(i,j)*sglace
Cfd0 &                 +dmnbq(i,j)*salflux
Cfd0      pmess      =  (prs-fons)/ddtb-fevabq(i,j)
Cfd0 &                  *salflux*albq(i,j)
          pme        = -fevabq(i,j)*albq(i,j)
     &                 +(zfwat(i,j)-zhnpbq(i,j)*(1.-ain(i,j))*rhon
     &                  -dmnbq(i,j))
     &                  /ddtb
C
C-AM: frapsav
          frapsav(i,j)=rappes(i,j,0)*
     &              (1.D0-scalr(i,j,ks2,2)/max(scal(i,j,ks2,2),salinf))
C
          DUM = (pme/rho0 + frapsav(i,j) * dzsdts) *ddtb
          phiss(i,j,0) = -DUM
     &                   +phiss(i,j,0)
C
Cdbug     if (abs(scal(i,j,ks2,2)).lt.epsil) then
Cdbug      write(iuo+66,*) 'ARRET : thersf, scal(i,j,ks2,2) too small'
Cdbug      write(iuo+66,*) 'scal(i=',i,',j=',j,')=',scal(i,j,ks2,2)
Cdbug      stop
Cdbug     endif
C
          watsur = -DUM * unsdz(ks2)
C
CC        phiss(i,j,1)  =  -(fcm1(i,j)+fcm2(i,j))
CC   &                      *ddtb*unsdz(ks2)/(rho0*cpo)
          phiss(i,j,1) =  -(fcm2(i,j)-fdtcn(i,j))
     &                      *ddtb*unsdz(ks2)/(rho0*cpo)
     &                      +watsur*scs(i,j,1)
     &                      +phiss(i,j,1)
Cfd0      phiss(i,j,2) =  pmess
Cfd0 &                      *ddtb*unsdz(ks2)/rho0
Cfd0 &                      +phiss(i,j,2)
          phiss(i,j,2) = -fons*unsdz(ks2)/rho0
     &                   +watsur*scs(i,j,2)
     &                   +phiss(i,j,2)
Ccfc      phiss(i,j,3) = watsur*scs(i,j,3)
Ccfc &                   +phiss(i,j,3)
Ccfc      phiss(i,j,4) = watsur*scs(i,j,4)
Ccfc &                   +phiss(i,j,4)
380     continue
390   continue
c
      do 430 k=11,ks2
         do 420 j=js1,js2
            do 410 i=is1(j),is2(j)
               phivs(i,j,k,1)=-reslum(i,j,k)*fcm1(i,j)
     &                         *ddtb
     &                         +phivs(i,j,k,1)
 410        continue
 420     continue
 430  continue

c geothermic flow
      do 450 j=js1,js2
        do 440 i=is1(j),is2(j)
           aflux=max(zero,dfloat(sign(1,11-kfs(i,j))) )
           phivs(i,j,kfs(i,j),1)=-aflux*bheat*
     &          dts(kfs(i,j))*unsdz(kfs(i,j))/(rho0*cpo)
C          if (i.eq.2) then
C             write(99,*) 'bot',j,kfs(i,j),aflux,
C    &              phivs(i,j,kfs(i,j),1),tms(i,j,kfs(i,j))
C          endif
 440     continue
 450  continue

c Iceberg melting
      fnor=(ficebergn/areicebn)
     &           *ddtb*unsdz(ks2)/(rho0*cpo*86400.)
      do i=iicebern1,iicebern2
       do j=jicebern1,jicebern2
          phiss(i,j,1) =  phiss(i,j,1) +fnor*tms(i,j,ks2)
       enddo
      enddo
      fsud=(ficebergs/areicebs)
     &           *ddtb*unsdz(ks2)/(rho0*cpo*86400.)
      do i=iicebers1,iicebers2
       do j=jicebers1,jicebers2
          phiss(i,j,1) =  phiss(i,j,1) +fsud*tms(i,j,ks2)
       enddo
      enddo
C     write(99,*) 'delta t ns', fnor, fsud
c     write(99,*) 'flux iceb', ficebergn/(areicebn*86400),
c    &                   ficebergs/(areicebs*86400)
       write(99,*) 'flux iceb2',fnor,fsud

c mean value of phiss(i,j,2)
        zflux0s=0.0
        do j=js1,js2
         do i=is1(j),is2(j)
           zflux0s=zflux0s+aire(i,j)*tms(i,j,ks2)*phiss(i,j,2)
         enddo
        enddo
        zflux0s=zflux0s/zurfow
        write(99,*) 'zflux salt',zflux0s,phiss(15,15,2)
c mean value of phiss2 set to zero on average
        dum=zfluxms*vcor
        do j=js1,js2
         do i=is1(j),is2(j)
C          phiss(i,j,2)=phiss(i,j,2)-zflux0*tms(i,j,ks2)
           phiss(i,j,2)=phiss(i,j,2)-dum*tms(i,j,ks2)
         enddo
        enddo
c Iceshelf melting

c     gammat is in m/s
      gammat=1E-4

      fluxmean=0.0
      zcums=0.0
      zcums2=0.0
      do j=2,8
       do i=is1(j),is2(j)
        tmean=0.0
        smean=0.0
        dztot=0.0
        do k=9,14
         tmean=tmean+scal(i,j,k,1)*dz(k)*tms(i,j,k)
         smean=smean+scal(i,j,k,2)*dz(k)*tms(i,j,k)
         dztot=dztot+dz(k)*tms(i,j,k)
        enddo
        dztot=max(dztot,1D-5)
        tmean=tmean/dztot*tms(i,j,ks2)
        smean=smean/dztot*tms(i,j,ks2)
        zdepthi=200.0
        tfreezloc= -0.0575*(smean)
     &            + 1.710523e-3*(smean)**1.5
     &            - 2.154996e-4*(smean)
     &                           *(smean)
     &            - 7.53e-4*zdepthi+273.15

C fluxt is the total amount of heat available to melt the ice shelf (W)
        fluxt=rho0*cpo*gammat*(tmean-tfreezloc)*tmics(i,j)
C       if (tmics(i,j).gt.0.0) then
C        write(99,*) i,j,tmean,(tmean-tfreezloc)
C       endif
C       write(99,*) fluxt,i,j,tmics(i,j)
C       write(99,*) tmean,tfreezloc,smean
        phivssup=fluxt*ddtb/(rho0*cpo*dztot*area(i,j))
        phivssups=fluxt*ddtb/(xlg*dztot*area(i,j))*salflux
C       phivssups=0.0
        do k=9,14
          phivs(i,j,k,1)=phivs(i,j,k,1)+phivssup*tms(i,j,k)
          phivs(i,j,k,2)=phivssups*tms(i,j,k)
        enddo
        fluxmean=fluxmean+fluxt
        zcums=zcums+fluxt/xlg*salflux
       enddo
      enddo
      toticesm=fluxmean/xlg*360*86400
c redistribution of the flux
      ficesh=(fluxmean/areicebs)
     &           *ddtb*unsdz(ks2)/(rho0*cpo)
      ficess=(zcums/areicebs)
     &           *ddtb*unsdz(ks2)
C    &           *ddtb*unsdz(ks2)/rho0
C     ficess=0.0

      do i=iicebers1,iicebers2
       do j=jicebers1,jicebers2
          phiss(i,j,1) =  phiss(i,j,1)-ficesh*tms(i,j,ks2)
          phiss(i,j,2) =  phiss(i,j,2)-ficess*tms(i,j,ks2)
       enddo
      enddo
      write(99,*) 'ice-shelf',ficesh,ficess

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C VERIF FLUXES
C     zsurfoz=0.0
C     zsurfzz=0.0
C     znsolo=0.0
C     zsolo=0.0
C     znsolo2=0.0
C     zsolo2=0.0
C     zlph0=0.0
C     zlph2=0.0
C     zvolg=0.0
C     zvoln=0.0
C     zfdt=0.0
C       do j=js1,js2
C        do i=is1(j),is2(j)
C         zsurfoz=zsurfoz+aire(i,j)*ain(i,j)
C         zsurfzz=zsurfzz+aire(i,j)
C         znsolo=znsolo+ratbqo(i,j)*ain(i,j)
C    &                 *aire(i,j)
C         zsolo=zsolo+fsolcn(i,j)*ain(i,j)*aire(i,j)
C         zsolo2=zsolo2+fcm1(i,j)*aire(i,j)
C         znsolo2=znsolo2+fcm2(i,j)*aire(i,j)
C         zvoln=zvoln+(1.0-albq(i,j))*aire(i,j)*hnbq(i,j)
C    &                              *0.33d0*tms(i,j,ks2)
C         zvolg=zvolg+(1.0-albq(i,j))*aire(i,j)
C    &                     *hgbq(i,j)*0.9*tms(i,j,ks2)
C         zlph0=zlph0+phiss(i,j,0)*aire(i,j)
C    &             *360*tms(i,j,ks2)
C    &             *rho0/1000.0
C         zlph2=zlph2+phiss(i,j,2)*aire(i,j)
C    &             *360/(unsdz(ks2)*salflux)*tms(i,j,ks2)
C    &             *rho0/1000.0
C         zfdt=zfdt+fdtcn(i,j)*aire(i,j)
C        enddo
C       enddo
C       write(93,'(4F11.5)') zsolo/zsurfoz,znsolo/zsurfoz,
C    &     zsolo2/zsurfzz,znsolo2/zsurfzz
Ca    &    ,zfdt/zsurfzz

      return

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine thersf -
      end

