












c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  fichier "bloc0.com" : incorpore par instruction 'include' dans 
c                        bloc.com,loch.com
c  creation : 10/02/2004 <- LOCH
 
c--blocs common :
 
      real*8 temp_vent(imax,jmax),temp_press(imax,jmax)
      common /atm2loch/temp_vent,temp_press
 
      real*8 dx,dy,unsdx,unsdy,uns2dx,uns2dy,z,dz,unsdz,zw,dzw,unsdzw,
     &        huy,hux,hu,unshu,tms,tmu
     &       ,xslon,xulon,yslat,yulat,xsedg,xuedg,ysedg,yuedg,
     &        xslonp,xulonp,yslatp,yulatp,hs,angle
      common / domain /
     &  dx, dy, unsdx, unsdy, uns2dx, uns2dy,
     &  z(kmax+1), dz(kmax), unsdz(kmax),
     &  zw(kmax+1), dzw(kmax+1), unsdzw(kmax+1),
     &  huy(imax,jmax), hux(imax,jmax), hu(imax,jmax),
     &  unshu(imax,jmax),
     &  tms(imax,jmax,kmax), tmu(imax,jmax,kmax)
cDFG
     &, xslon(imax,jmax), xulon(imax,jmax),
     &  yslat(imax,jmax), yulat(imax,jmax),
     &  xsedg(imax,jmax,4), xuedg(imax,jmax,4),
     &  ysedg(imax,jmax,4), yuedg(imax,jmax,4),
     &  xslonp(imax+1,jmax+1), xulonp(imax+1,jmax+1),
     &  yslatp(imax+1,jmax+1), yulatp(imax+1,jmax+1),
     &  hs(imax,jmax),angle(imax,jmax)

      integer msks, msku
      common / idomain /
     & msks(imax,jmax), msku(imax,jmax)
cDFG

 
c--WARNING : do not modify "dynam1" without checking vrl() equivalence !
      real*8 tpstot,fss,daeta,eta,ub,vb,u,v,scal,b,bvf,
     &       avsdz,avudz,w,q2turb,fqajc
      common / dynam1 / tpstot,
     &  fss(imax,jmax,0:nsmax), daeta(imax,jmax),
     &  eta(imax,jmax), ub(imax,jmax), vb(imax,jmax),
     &  u(imax,jmax,kmax), v(imax,jmax,kmax),
     &  scal(imax,jmax,kmax,nsmax), b(imax,jmax,kmax),
     &  bvf(imax,jmax,kmax), avsdz(imax,jmax,kmax),
     &  avudz(imax,jmax,kmax), w(imax,jmax,kmax+1),
     &  q2turb(imax,jmax,kmax+1), fqajc(imax,jmax,kmax)
 
      real*8 scal0, spvr, scalr,rappel, rappes, ahrap,
     &  rapint, phimnx,phifu, phifv, phifs,phisu,phisv, phiss
     & ,scpme, scssv, scs,flxus, flxvs,phivs,ust2s,ust2b
      common / forcing /
     &  scal0(kmax,nsmax), spvr, scalr(imax,jmax,kmax,nsmax),
     &  rappel(imax,jmax,kmax), rappes(imax,jmax,0:nsmax), ahrap,
     &  rapint(nrpmax,kmax,nsmax), phimnx(imax,jmax,0:1,0:nsmax),
     &  phifu(imax,jmax), phifv(imax,jmax), phifs(imax,jmax,0:nsmax),
     &  phisu(imax,jmax), phisv(imax,jmax), phiss(imax,jmax,0:nsmax)
     & ,scpme(nsmax), scssv(nsmax), scs(imax,jmax,nsmax)
Cadh & ,flxus(imax,jmax,nsmax), flxvs(imax,jmax,nsmax)
     & ,phivs(imax,jmax,kmax,nsmax)
     & ,ust2s(imax,jmax),ust2b(imax,jmax)
 
      integer nitrun,nsav
      common / lerun2 /nitrun,nsav
 
      integer kniv,kfs, kfu, ks1,ks2, ku1, ku2,is1, is2, ims1,
     &  ims2, js1, js2,iu1, iu2, imu1, imu2, ju1, ju2,isf1, isf2,
     &  iuf1, iuf2,jcl1, jcl2, jeq, jdl1, jdl2, ijsdl, ijudl,
     &  iberp, jberp, ibera, jbera, iberpm, jberpm, iberam, jberam
      common / limites /
     &  kniv(imax,jmax,-1:1),
     &  kfs(imax,jmax), kfu(imax,jmax), ks1, ks2, ku1, ku2,
     &  is1(jmax) , is2(jmax) , ims1, ims2, js1, js2,
     &  iu1(jmax) , iu2(jmax) , imu1, imu2, ju1, ju2,
     &  isf1(jmax), isf2(jmax), iuf1(jmax), iuf2(jmax),
     &  jcl1, jcl2, jeq, jdl1, jdl2, ijsdl, ijudl,
     &  iberp, jberp, ibera, jbera, iberpm, jberpm, iberam, jberam
 
       real*8 dtu, dtb, cdbot,ahu, ahe, alphxu, alphxv, alphyu, alphyv,
     &  slopemax,ahs,ai,
     &  aitd, slopmgm, afilt, ahh, avv,alphgr, algrmn, alphah, alphmi,
     &  alphaz, rifsmx, rifumx,avnu0, avnub, avk0, avkb,qavs, avsn2,
     &  ccfmn, ccfmx, txiadu, txiads, txidfu, txidfs, txeflx,
     &  txifcb, txifcu, bering, ajcmix, xslop
      common / runpara /
     &  dtu, dtb, cdbot,
     &  ahu, ahe, alphxu, alphxv, alphyu, alphyv,
     &  slopemax(kmax),ahs(kmax),ai(kmax),
     &  aitd(kmax), slopmgm(kmax), afilt, ahh, avv,
     &  alphgr(nsmax), algrmn(nsmax), alphah(2), alphmi(kmax),
     &  alphaz(kmax), rifsmx, rifumx,
     &  avnu0(kmax), avnub(kmax), avk0(kmax), 
     &  avkb(kmax),
     &  qavs, avsn2, ccfmn, ccfmx,
     &  txiadu, txiads, txidfu, txidfs, txeflx(nsmax), txifcb, txifcu,
     &  bering, ajcmix, xslop
 
      real*8 dts
      common / runpara2 / dts(kmax)

      real*8 frapsav
      common /fluxsel/frapsav(imax,jmax)
 
      real*8 zurfow,zflux0,zfluxm,vcor,bheat,zflux0s,zfluxms,xfreshsv
      integer ihyster
      common /volcor/
     &  zurfow,zflux0,zfluxm,vcor,bheat,zflux0s,zfluxms,
     &  xfreshsv(6000),ihyster

c--fin du fichier "bloc0.com"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
