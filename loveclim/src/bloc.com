












c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  fichier "bloc.com" : incorpore par instruction 'include' dans les routines :
c      OM, CLASS, FEEDOM, GRADP, STATIS, DIFBIN, TROMNC, UNIBIN,
c      start, flucor, scale, slopez, slopes, scadew, scadns, scali,
c      uve, uvi, barot, uvb0et, uvbfet, uvm,
c      informe, defcst, defgrid, redforc, correct,
c      conti3d, etat, vdiffu, alph2dc, alphdkc, alphdec, raccord,
c      savrunb, redrunb, savrunc, redrunc,
c      ncdfout, moyen, streamv, meridflu, scadhv, checkwst, streamh, stream1h,
c      vague, local, binout, defgl, unigl, foroutp, lisstab, flowsurf.
c  inclus apres "type.com", "para.com".
c  modif : 06/10/99
c  modif : 02/02/04 <- LOCH : 'bloc0.com'
c  modif : 01/03/04 <- LOCH
 
c--blocs common :
 
c--Tableaux de variables en evolution (+ transfert entre routines) :

      integer numit,numspl,ninstb,nclin,nclmoy
      common / iteration /
     &  numit, numspl, ninstb, nclin, nclmoy
 
c--Tableaux de flux et facteurs d'evolution et variables en evolution :

      real*8 deriv,q,fub,fvb,etaspl,ubspl,vbspl,vlturb,avqdz,
     &  tm2tur,avuloc,phizzz,phihhh,umoy,vmoy,alpslp,uvcslp
      common / dynam2 / deriv(nsmax),
     &  q(imax,jmax,kmax), fub(imax,jmax,kmax), fvb(imax,jmax,kmax),
     &  etaspl(imax,jmax), ubspl(imax,jmax), vbspl(imax,jmax),
     &  vlturb(imax,jmax,kmax), avqdz(imax,jmax,kmax),
     &  tm2tur(imax,jmax,kmax), avuloc(imax,jmax,kmax),
     &  phizzz(imax,jmax,kmax+1,2), phihhh(imax,jmax,6),
     &  umoy(imax,jmax), vmoy(imax,jmax),
     &  alpslp(nlpmax), uvcslp(nlpmax)
 
Csai  common / saison /
Csai &  tmens(imax,jmax,nmois),smens(imax,jmax,nseas),
Csai &  txmens(imax,jmax,nmois),tymens(imax,jmax,nmois),
Csai &  d2tmns(imax,jmax,nmois),d2smns(imax,jmax,nseas),
Csai &  d2txms(imax,jmax,nmois),d2tyms(imax,jmax,nmois),
Csai &  flxss(imax,jmax,nsmax,nmois),d2fxss(imax,jmax,nsmax,nmois),
Csai &  flxsur(imax,jmax,nsmax)

      integer  master,icoupl,icoutp,itau_slow
      common / icouplage /
     &  master, icoupl, icoutp, itau_slow
 
      real*8 tfmocn,ttoocn
      common / tbforcing /
     &  tfmocn(imax,jmax,ntocn), ttoocn(imax,jmax,ntatm)

      real*8 cmx,cmy,smx,smy,cmxy,smxy,cmxdy,cmydx,fs2cor,
Cfcc &  fcucor, fcvcor,
     &  aire,covrai,xang1,xang2 
      common / cfmetric /
     &  cmx(imax,jmax,0:3), cmy(imax,jmax,0:3),
     &  smx(imax,jmax,0:3), smy(imax,jmax,0:3),
     &  cmxy(imax,jmax,0:3), smxy(imax,jmax,0:3),
     &  cmxdy(imax,jmax), cmydx(imax,jmax),
     &  fs2cor(imax,jmax), aire(imax,jmax),
Cfcc &  fcucor(imax,jmax), fcvcor(imax,jmax),
     &  covrai(imax,jmax), xang1(imax,jmax),
     &  xang2(imax,jmax)
 
      real*8 ctmi,zsurfs,zsurfo,zsurfv,zvols,zvolo,zvolv,zvolw, 
     &  zsurf,zsurfsla,zsurfola,zsurfsba,zsurfoba,zvolsla, 
     &  zvolola,zvolsba,zvoloba,unsvol
      common / surfvol /
     &  ctmi(imax,jmax,kmax,0:1),
     &  zsurfs(kmax), zsurfo(kmax), zsurfv(kmax),
     &  zvols, zvolo, zvolv, zvolw, zsurf,
     &  zsurfsla(kmax,0:nltmax), zsurfola(kmax,0:nltmax),
     &  zsurfsba(kmax,0:nbsmax), zsurfoba(kmax,0:nbsmax),
     &  zvolsla(0:nltmax), zvolola(0:nltmax),
     &  zvolsba(0:nbsmax), zvoloba(0:nbsmax),
     &  unsvol

      integer ijrap,nrap,i1coin,i2coin,i3coin,i4coin,
     &  n1coin,n2coin,n3coin,n4coin,ijslp,kslp,
     &  lslp,nslp, nslpdw,nxslp,nxyslp 
      common / localise /
     &  ijrap(nrpmax,kmax,nsmax), nrap(kmax,nsmax),
     &  i1coin(ncomax,kmax), i2coin(ncomax,kmax),
     &  i3coin(ncomax,kmax), i4coin(ncomax,kmax),
     &  n1coin(kmax), n2coin(kmax), n3coin(kmax), n4coin(kmax),
     &  ijslp(nlpmax), kslp(nlpmax), lslp(nlpmax),
     &  nslp(nlpmax), nslpdw, nxslp, nxyslp

      integer nstart, nlast, nsplaj, nsplit, ninfo,
     &  nwjl , nwtal , nwm , nwa , nwtest, idyn,
     &  nwtoa, nwtom,
     &  lstab, nsewfr, kstart, kinput, koutpu, nitrap, ntmoy,
     &  kfond, kforc, mdforc, kavsmx, ntout 
      common / lerun /
     &  nstart, nlast, nsplaj, nsplit, ninfo,
     &  nwjl , nwtal , nwm , nwa , nwtest, idyn,
     &  nwtoa, nwtom,
     &  lstab, nsewfr, kstart, kinput, koutpu, nitrap, ntmoy,
     &  kfond, kforc, mdforc, kavsmx, ntout
 
      real*8 vkappa,q2tmin,ghmax,ghmin,zlotur,vlmin,varfor,sqrghm
      integer kajul
      common / coetur/
     &  vkappa,q2tmin,ghmax,ghmin,zlotur,vlmin,varfor,sqrghm,
     &  kajul
 
      character*6 refexp
      common / cerun /
     & refexp
 
      INCLUDE 'bloc0.com'
 
c--fin du fichier "bloc.com"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
