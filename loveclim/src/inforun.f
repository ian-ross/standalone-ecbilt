












      subroutine inforun(nn99)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c    calcule et ecrit  quelques variables globales, sur le fichier "evolu",
c    frequence "ninfo", serie chronologique.
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  modif : 09/08/02
 

      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'reper.com'
      include 'isoslope.com'
      include 'ice.com'
      include 'dynami.com'
 
      real*8  heaticbism(imax,jmax)
      real*8  vicbismn,vicbisms
      real*8  moc,tmc,tmc0,tsurfmean,cland,thex
      logical flgveg,flgicb,flgisma,flgismg

      common /icbism/ heaticbism,vicbismn,vicbisms
      common /ec_coupl/ flgveg,flgicb,flgisma,flgismg
      common/IPCC_out2/moc,tmc,tmc0,tsurfmean,cland,thex 
c--variables locales :
      dimension nnvk(nsmax), vvk(3), zitsum(0:1)
      logical flgout
 
c--variables locales a conserver d'un appel a l'autre -> dans "reper.com" .
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Prepare & debute le remplissage de "vinfor".                    |
c-----------------------------------------------------------------------
 
c-----------------------------------------------------------------------
c  => acces pendant le run, pour calculer et/ou ecrire Var.evolution
c-----------------------------------------------------------------------
 
c- flag = TRUE <=> write on "evolu" file at this iter. :
      flgout = mod(numit,ninfo).eq.0 .or. numit.eq.1
 
c- Initialisation (a zero) de "vinfor":
      do 300 nv=1,nvinfo
        vinfor(nv) = 0.
 300  continue
 
      vinfor(1) = DFLOAT(numit)
      vinfor(2) = tpstot / ( 86400. * yeaday )
      nv = 2
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Calcul de quelques variables globales :                         |
c-----------------------------------------------------------------------
 
      if (flgout) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--2.1 Variables deja cumulees dans une autre routine :
 
c-----
c- evaluation de la variation de eta au cours de l'it. par W(kfond) :
c  -->  transferre le 07/12/93 dans "flucor"
c-----
 
c--Moyenne sur tout le bassin de : 1) l'energie dissipee par Ajust.Conv.
c                                  2) variation de l'elevat. au cours de l'it.
      do 410 j=js1,js2
       do 410 i=is1(j),is2(j)
        vinfor(nv+1) = vinfor(nv+1) + ctmi(i,j,ks2,0) * fqajc(i,j,1)
        vinfor(nv+3) = vinfor(nv+3) + ctmi(i,j,ks2,0) * daeta(i,j)
        daeta(i,j) = 0.
 410  continue
      vinfor(nv+1) = vinfor(nv+1) * zsurf
      vinfor(nv+3) = vinfor(nv+3) * zmdeta
 
c--Integre la Freq.d'AjC sur le volume oceanique global :
      do 430 k=ks1+1,ks2
        sumk = 0.
        kk = k - 1
        do 420 j=js1,js2
         do 420 i=is1(j),is2(j)
           sumk = sumk + ctmi(i,j,kk,0) * fqajc(i,j,k)
 420    continue
        vinfor(nv+2) = vinfor(nv+2) + sumk * dzw(k)
 430  continue
c- fraction du volume total, en o/oo :
      vinfor(nv+2) = vinfor(nv+2) * zvolw * 1000.
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
      nv = nv + 4
 
      if (flgout .or. ntmoy.eq.2) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--2.2 Variables non-cumulees dans une autre routine :
 
c--Moyenne sur tout le bassin de : 3) de l'elevation
c--            (ancienne version : 3) de l'energ. due a l'elevation)
      do 440 j=js1,js2
       do 440 i=is1(j),is2(j)
        vinfor(nv) = vinfor(nv)
     &             + ctmi(i,j,ks2,0) * eta(i,j)
C    &             + ctmi(i,j,ks2,0) * eta(i,j) * eta(i,j)
 440  continue
C     vinfor(nv) = 0.5 * gpes * rho0 * vinfor(nv) * zsurf
      vinfor(nv) = vinfor(nv) * zsurf
c-----
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--2.3 Horizontal Transport (Sverdrup) , entre qq points particuliers :
c- passage de Drake (nhsf=1), Indonesie (nhsf=2)
      do 470 nhsf=1,nvhsf
        nn = nv + nhsf
        if (iehsf(nhsf).eq.0) then
c- Integre (-ub).(-dy) :
          ii = ishsf(nhsf)
          do 450 j=jshsf(nhsf),jehsf(nhsf)
            vinfor(nn) = vinfor(nn) + cmy(ii,j,3) * ub(ii,j)
 450      continue
          vinfor(nn) = vinfor(nn) * dy * svrdrp
        else
c- Integre vb.dx :
          jj = jshsf(nhsf)
          do 460 i=ishsf(nhsf),iehsf(nhsf)
            vinfor(nn) = vinfor(nn) + cmx(i,jj,3) * vb(i,jj)
 460      continue
          vinfor(nn) = vinfor(nn) * dx * svrdrp
        endif
 470  continue
      nv = nv + nvhsf
c-- Florida Strait (i=177, j=79,83)
CL30  ii=177
CL30  jj1=79
CL30  jj2=83
CL30  nv=nv+1
CL30  do 475 j=jj1,jj2
CL30    vinfor(nv) = vinfor(nv) + cmy(ii,j,3) * ub(ii,j)
 475  continue
CL30  vinfor(nv) = vinfor(nv) * dy * svrdrp
 
c--2.4 Horizontal Transport (Sverdrup) , entre qq points particuliers :
c- calcul separe Flux Sud / Flux Nord ou Flux Ouest / Flux Est
c-----
c  4 = Dan.Strait   , 5 = Icel.-Scotl , 6 = Fram Strait , 7 = N.W.Passage
c  8 = Spitz-Norway , 9 = Gibraltar ,  10 = Austr.-NewZe
c-----
Cic0  do 520 nhsf=4,min(6,ndhsf)
      do 520 nhsf=4,min(7,ndhsf)
       if (iehsf(nhsf).eq.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Integre (u.dy.dz) :
        i = ishsf(nhsf)
        do 490 k=ks1,ks2
          vvn = 0.
          vvp = 0.
          do 480 j=jshsf(nhsf),jehsf(nhsf)
            vvn = vvn + tmu(i,j,k) * cmy(i,j,3) * min(zero, u(i,j,k))
            vvp = vvp + tmu(i,j,k) * cmy(i,j,3) * max(zero, u(i,j,k))
 480      continue
          vinfor(nv+1) = vinfor(nv+1) + vvn * dz(k)
          vinfor(nv+2) = vinfor(nv+2) + vvp * dz(k)
 490    continue
        vinfor(nv+1) = vinfor(nv+1) * dy * svrdrp
        vinfor(nv+2) = vinfor(nv+2) * dy * svrdrp
       else
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Integre (v.dx.dz) :
        j = jshsf(nhsf)
        do 510 k=ks1,ks2
          vvn = 0.
          vvp = 0.
          do 500 i=ishsf(nhsf),iehsf(nhsf)
            vvn = vvn + tmu(i,j,k) * cmx(i,j,3) * min(zero, v(i,j,k))
            vvp = vvp + tmu(i,j,k) * cmx(i,j,3) * max(zero, v(i,j,k))
 500      continue
          vinfor(nv+1) = vinfor(nv+1) + vvn * dz(k)
          vinfor(nv+2) = vinfor(nv+2) + vvp * dz(k)
 510    continue
        vinfor(nv+1) = vinfor(nv+1) * dx * svrdrp
        vinfor(nv+2) = vinfor(nv+2) * dx * svrdrp
       endif
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
       nv = nv + 2
 520  continue
 
c-----
 
c--2.5 diagnostics of the THC circulation
c- Atlantique Prod GIN (= le Max entre 68 et 75 N)
      zold=0
      if (zold.eq.1) then
      yy = 68.
      jj1= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      yy = 75.
      jj2= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      do 530 jj=jj1,jj2
       vvpsm = 0.0
       do 530 k=ks1,ks2
         if (z(k).gt.zz) goto 530
         vv = 0.0
         do i=iszon(jj,nbsmax),iezon(jj,nbsmax)
           vv = vv + cmx(i,jj,3)*v(i,jj,k)
         enddo
         vvpsm = vvpsm - vv * dz(k)
         vinfor(nn) = max(vinfor(nn),vvpsm)
 530  continue
      vinfor(nn) = vinfor(nn) * dx * svrdrp
      nv = nv + 1
c- Atlantique Prod (= le Max entre 45 et 65 N)
      yy = 45.
      jj1= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      yy = 65.
      jj2= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      do 540 jj=jj1,jj2
       vvpsm = 0.0
       do 540 k=ks1,ks2
         if (z(k).gt.zz) goto 540
         vv = 0.0
         do i=iszon(jj,nbsmax),iezon(jj,nbsmax)
          vv = vv + cmx(i,jj,3)*v(i,jj,k)
         enddo
         vvpsm = vvpsm - vv * dz(k)
         vinfor(nn) = max(vinfor(nn),vvpsm)
 540  continue
      vinfor(nn) = vinfor(nn) * dx * svrdrp
      nv = nv + 1
c- NADW Outflow and AABW inflow into the Atlantic
c    (= le Max et le min a 20 S)
      yy = -20.0
      jj = 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      vvpsm = 0.0
      do 550 k=ks1,ks2
        if (z(k).gt.zz) goto 550
        vv = 0.0
        do i=iszon(jj,nbsmax),iezon(jj,nbsmax)
          vv = vv + cmx(i,jj,3)*v(i,jj,k)
        enddo
        vvpsm = vvpsm - vv * dz(k)
        vinfor(nn)  = max(vinfor(nn),vvpsm)
        vinfor(nn+3)= min(vinfor(nn+3),vvpsm)
 550  continue
      vinfor(nn) = vinfor(nn) * dx * svrdrp
      vinfor(nn+3) = vinfor(nn+3) * dx * svrdrp
      nv = nv + 1
c- AABW Prod (= le Min Sud de 60S )
      yy = -70.0
      jj1= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      yy = -60.0
      jj2= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      do 560 jj=jj1,jj2
       vvpsm = 0.0
       do 560 k=ks1,ks2
         if (z(k).gt.zz) goto 560
         vv = 0.0
         do i=iu1(jj),iu2(jj)
           vv = vv + cmx(i,jj,3)*v(i,jj,k)
         enddo
         vvpsm = vvpsm - vv * dz(k)
         vinfor(nn) = min(vinfor(nn),vvpsm)
 560  continue
      vinfor(nn) = vinfor(nn) * dx * svrdrp
      nv = nv + 1
c- AABW export (= le Min a 30S)
      yy = -30.0
      jj= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      vvpsm = 0.0
      do 570 k=ks1,ks2
        if (z(k).gt.zz) goto 570
        vv = 0.0
         do i=iu1(jj),iu2(jj)
          vv = vv + cmx(i,jj,3)*v(i,jj,k)
         enddo
        vvpsm = vvpsm - vv * dz(k)
        vinfor(nn)  = min(vinfor(nn),vvpsm)
 570  continue
      vinfor(nn) = vinfor(nn) * dx * svrdrp
      nv = nv + 2
      else
C     write(99,*) start new diag
      nm1n=nv
c- Atlantic Prod GIN (= Max between 68 and 75 N)
      yy = 68
      jj1= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      yy = 75.
      jj2= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      vinfor(nn)=0.0
      do 571 jj=jj1,jj2
        do 571 k=ks1,ks2
         if (z(k).gt.zz) goto 571
         vinfor(nn) = max(vinfor(nn),(vwx(jj,k,3)))
 571  continue
      vinfor(nn) = vinfor(nn) * dx * svrdrp
      nv = nv + 1
c- Atlantique Prod (= Max between 45 and 75 N)
      yy = 45.
      jj1= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      yy = 75.
      jj2= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      vinfor(nn)=0.0
      do 572 jj=jj1,jj2
        do 572 k=ks1,ks2
         if (z(k).gt.zz) goto 572
C        write(99,*) jj,k,vinfor(nn),vwx(jj,k,3)
         vinfor(nn) = max(vinfor(nn),(vwx(jj,k,3)))
 572  continue
      moc=vinfor(nn)
      nv = nv + 1
c- NADW Outflow and AABW inflow into the Atlantic
c    (= Max and  min at 20 S)
      yy = -20.0
      jj = 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      vinfor(nn)=0.0
      vinfor(nn+3)=0.0
      do 573 k=ks1,ks2
        if (z(k).gt.zz) goto 573
        vinfor(nn)  = max(vinfor(nn),(vwx(jj,k,3)))
        vinfor(nn+3)= min(vinfor(nn+3),(vwx(jj,k,3)))
 573  continue
C     write(99,*) 'NADW,AABW',vinfor(nn),vinfor(nn+3)
      nv = nv + 1
c- AABW Prod (= Min southward of 60S )
      yy = -70.0
      jj1= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      yy = -60.0
      jj2= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      vinfor(nn)=0.0
      do 574 jj=jj1,jj2
        do 574 k=ks1,ks2
         if (z(k).gt.zz) goto 574
         vinfor(nn) = min(vinfor(nn),(vwx(jj,k,0)))
 574  continue
      nv = nv + 1
c- AABW export (= Min at 30S)
      yy = -30.0
      jj= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      vinfor(nn)=0.0
      do 575 k=ks1,ks2
       if (z(k).gt.zz) goto 575
       vinfor(nn) = min(vinfor(nn),(vwx(jj,k,0)))
 575  continue
      nv = nv + 2
      endif
      nm2n=nv
c- Heat flux at 30 S
      yy = -30.0
      jj= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      vvpsm = 0.0
      ssc2 = 2.0 * scal0(ks1,1)
C     ssc2 = 2.0 * 273.15
      do 580 k=ks1,ks2
        vv = 0.0
        ccydif = ahs(k) * unsdy
        do i=iszon(jj,nbsmax),iezon(jj,nbsmax)
          vv2 = 0.5 * (v(i,jj,k)+v(i+1,jj,k))
     &        + viso(i,jj,k)
          v2cd2 = 0.5 * cmx(i,jj,2) * vv2
          cny = smy(i,jj,2) * unsdy * dts(k) * vv2
          aalpha = min( one, abs(cny) + alphmi(k) )
          phiy = tms(i,jj,k) * tms(i,jj-1,k) * (
     &         v2cd2 * (scal(i,jj-1,k,1) + scal(i,jj,k,1) - ssc2)
C    &         + ( alphay(i,jj,k) + cmxy(i,jj,2)*ccydif )
     &         + (aalpha*abs(v2cd2) + cmxy(i,jj,2)*ccydif)
     &               * (scal(i,jj-1,k,1) - scal(i,jj,k,1))  )
c-----
C         phi1=tms(i,jj,k) * tms(i,jj-1,k) * (
C    &          v2cd2 * (scal(i,jj-1,k,1) + scal(i,jj,k,1) - ssc2 ))
C         phi2=tms(i,jj,k) * tms(i,jj-1,k) * (
C    &        cmxy(i,jj,2) * ccydif
C    &               * (scal(i,jj-1,k,1) - scal(i,jj,k,1))  )
C         write(95,*) phi1,phi2,phiy
          vv = vv + phiy
        enddo
        vvpsm = vvpsm + vv*dz(k)
580   continue
      nv = nv + 1
      vinfor(nv) = vvpsm * dx * rho0*cpo * 1.D-15
c- Salt flux at 30 S
      yy = -30.0
      jj= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      vvpsm = 0.0
      ssc2 = 2.0 * scal0(ks1,2)
      do 581 k=ks1,ks2
        vv = 0.0
        ccydif = ahs(k) * unsdy
        do i=iszon(jj,nbsmax),iezon(jj,nbsmax)
         vv2 = 0.5 * (v(i,jj,k)+v(i+1,jj,k))
     &         + viso(i,jj,k)
         v2cd2 = 0.5 * cmx(i,jj,2) * vv2
         cny = smy(i,jj,2) * unsdy * dts(k) * vv2
         aalpha = min( one, abs(cny) + alphmi(k) )
         phiy = tms(i,jj,k) * tms(i,jj-1,k) * (
     &       v2cd2 * (scal(i,jj-1,k,2) + scal(i,jj,k,2) - ssc2)
     &       + (aalpha*abs(v2cd2) + cmxy(i,jj,2)*ccydif)
     &       * (scal(i,jj-1,k,2) - scal(i,jj,k,2))  )
         vv = vv + phiy
       enddo
       vvpsm = vvpsm + vv*dz(k)
581   continue
      nv = nv + 1
      vinfor(nv) = vvpsm * dx
c- Salt flux at Bering
      jj=65
      ii=102
      vvpsm = 0.0
      ssc2 = 2.0 * scal0(ks1,2)
      vv = 0.0
      ccydif = ahs(ks2) * unsdy
        vv2 = 0.5 * (v(ii,jj,ks2)+v(ii+1,jj,ks2))
     &         + viso(ii,jj,ks2)
        v2cd2 = 0.5 * cmx(ii,jj,2) * vv2
        cny = smy(ii,jj,2) * unsdy * dts(ks2) * vv2
        aalpha = min( one, abs(cny) + alphmi(ks2) )
         phiy = tms(ii,jj,ks2) * tms(ii,jj-1,ks2) * (
     &       v2cd2 * (scal(ii,jj-1,ks2,2) + scal(ii,jj,ks2,2) - ssc2)
     &       + (aalpha*abs(v2cd2) + cmxy(ii,jj,2)*ccydif)
     &       * (scal(ii,jj-1,ks2,2) - scal(ii,jj,ks2,2))  )
         vv = vv + phiy
       vvpsm = vvpsm + vv*dz(ks2)
       
      nv = nv + 1
      vinfor(nv) = vvpsm * dx

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- downslpoing out of the antarctic shelf ( from j=1 to jmsud ) :
cd     nv = nv + 1
cd    vinfor(nv) = 0.0
cd    ijlim = jmsud*imax
cd    do nnt=1,nslpdw
cd     nl=nslp(nnt)
cd     if ( kslp(nl).ge.ksud .and. ijslp(nl).le.ijlim )
cd   &   vinfor(nv) = vinfor(nv) + abs(uvcslp(nnt))
cd    enddo
cd    vinfor(nv) = vinfor(nv) * dx * svrdrp
c- downslpoing out of the arctic ( from j=jmnor to jmax ) :
cd    nv = nv + 1
cd    vinfor(nv) = 0.0
cd    ijlim = 1+(jmnor-1)*imax
cd    do nnt=1,nslpdw
cd     nl=nslp(nnt)
cd     if ( kslp(nl).ge.knor .and. ijslp(nl).ge.ijlim )
cd   &   vinfor(nv) = vinfor(nv) + abs(uvcslp(nnt))
cd    enddo
cd    vinfor(nv) = vinfor(nv) * dx * svrdrp
cd
      if (flgout .or. ntmoy.eq.2) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--2.6 Integrale globale a faible variabilite temporelle :
 
c--Moyenne sur tout le bassin des Scal., |Scal-Obs| et Scal.-Obs (Niv 1) :
      nnp = 3
      if (ninfo.lt.0) nnp = 2
      nnvk(2) = nocean + ks1
      nnvk(1) = nnvk(2) - ks2
      do 730 ns=1,nsmax
c- 1er niveau (surface) :
        k = ks2
          vvk(1) = 0.
          vvk(2) = 0.
          vvk(3) = 0.
          do 700 j=js1,js2
           do 700 i=is1(j),is2(j)
            ctmobs = min( ctmi(i,j,k,0), (scalr(i,j,k,ns)-spvr) )
            vvk(1) = vvk(1) + ctmi(i,j,k,0) * scal(i,j,k,ns)
            vv2 = ctmobs * (scal(i,j,k,ns) - scalr(i,j,k,ns))
            vvk(2) = vvk(2) + vv2
            vvk(3) = vvk(3) + abs(vv2)
 700      continue
          vinfor(nv+1) = vvk(1) * dz(k)
          vinfor(nv+3) = vvk(3) * dz(k)
          vinfor(nv+2) = vvk(2) * zsurfo(k)
          if (ns.le.2) vinfor(nnvk(ns)-k) = vvk(3) * zsurfo(k)
c- autres niveaux :
        do 720 k=ks2-1,ks1,-1
          vvk(1) = 0.
          vvk(2) = 0.
          vvk(3) = 0.
          do 710 j=js1,js2
           do 710 i=is1(j),is2(j)
            ctmobs = min( ctmi(i,j,k,0), (scalr(i,j,k,ns)-spvr) )
            vvk(1) = vvk(1) + ctmi(i,j,k,0) * scal(i,j,k,ns)
            vv2 = ctmobs * (scal(i,j,k,ns) - scalr(i,j,k,ns))
            vvk(2) = vvk(2) + vv2
            vvk(3) = vvk(3) + abs(vv2)
 710      continue
          vinfor(nv+1) = vinfor(nv+1) + vvk(1) * dz(k)
          vinfor(nv+3) = vinfor(nv+3) + vvk(3) * dz(k)
          if (ns.le.2) vinfor(nnvk(ns)-k) = vvk(nnp) * zsurfo(k)
 720    continue
        vinfor(nv+1) = vinfor(nv+1) * zvols + scalwr(ns)
        if (ns.eq.1) tmc=vinfor(nv+1)*cpo*rho0
        vinfor(nv+3) = vinfor(nv+3) * zvolo
        nv = nv + 3
 730  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--Moyenne sur tout le bassin des |w| (m/s) :
      nv = nv + 1
      do 750 k=ks1+1,ks2
        vvk(1)=0.
        do 740 j=js1,js2
         do 740 i=is1(j),is2(j)
          vvk(1) = vvk(1) + ctmi(i,j,k-1,0) * abs(w(i,j,k))
 740    continue
        vinfor(nv) = vinfor(nv) + vvk(1) * dzw(k)
 750  continue
      vinfor(nv) = vinfor(nv) * zvolw
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--Moyenne sur tout le bassin des Vitesses (m/s) :
C     nvkek = nnvk(1) - ks2
      factke = rho0 * 0.5
      do 790 k=ks1,ks2
        do 760 n=1,3
          vvk(n) = 0.
 760    continue
        do 770 j=ju1,ju2
         do 770 i=iu1(j),iu2(j)
          vvk(1) = vvk(1) + cmxy(i,j,3) * abs(u(i,j,k))
          vvk(2) = vvk(2) + cmxy(i,j,3) * abs(v(i,j,k))
          vvk(3) = vvk(3) + cmxy(i,j,3) *
     &           ( u(i,j,k) * u(i,j,k) + v(i,j,k) * v(i,j,k) )
 770    continue
        vvk(3) = vvk(3) * factke
C       vinfor(nvkek-k) = vvk(3) * zsurfv(k)
        do 780 n=1,3
          vinfor(nv+n) = vinfor(nv+n) + vvk(n) * dz(k)
 780    continue
 790  continue
      do 800 n=1,3
        vinfor(nv+n) = vinfor(nv+n) * zvolv
 800  continue
      nv = nv + 3
c-----
      endif
      nv = nocean
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--2.7 Integrale globale (ou partielle) a forte variabilite temporelle :
 
      do 820 j=jeq,js2
        do 810 i=is1(j),is2(j)
          if (tms(i,j,ks2).eq.1) then
            vinfor(nv+1)=vinfor(nv+1)+(1.0-albq(i,j))*aire(i,j)
            if (albq(i,j).lt.0.15) vinfor(nv+3)=vinfor(nv+3) + aire(i,j)
            if (albq(i,j).lt.0.85) vinfor(nv+5)=vinfor(nv+5) + aire(i,j)
            vinfor(nv+7)=vinfor(nv+7)+albq(i,j)*aire(i,j)*(1-albq(i,j))
     &                        /max((1-albq(i,j)),1e-6*one)
            vinfor(nv+9)=vinfor(nv+9)
     &                   +(1.0-albq(i,j))*aire(i,j)*hgbq(i,j)
            vinfor(nv+11)=vinfor(nv+11)
     &                   +(1.0-albq(i,j))*aire(i,j)*hnbq(i,j)
            vinfor(nv+13)=vinfor(nv+13)
     &                   +(1.0-albq(i,j))*aire(i,j)*hgbq(i,j)*
     &                    (ug(i,j)*ug(i,j)+vg(i,j)*vg(i,j))
          endif
 810    continue
 820  continue
      vinfor(nv+13)=sqrt(vinfor(nv+13)/max(vinfor(nv+9),1e-5*one))
      nv=nv+1
 
      do 840 j=js1,jeq-1
        do 830 i=is1(j),is2(j)
          if (tms(i,j,ks2).eq.1) then
            vinfor(nv+1)=vinfor(nv+1)+(1.0-albq(i,j))*aire(i,j)
            if (albq(i,j).lt.0.15) vinfor(nv+3)=vinfor(nv+3) + aire(i,j)
            if (albq(i,j).lt.0.85) vinfor(nv+5)=vinfor(nv+5) + aire(i,j)
            vinfor(nv+7)=vinfor(nv+7)+albq(i,j)*aire(i,j)*(1-albq(i,j))
     &                        /max((1-albq(i,j)),1e-6*one)
            vinfor(nv+9)=vinfor(nv+9)+(1.0-albq(i,j))
     &                   *aire(i,j)*hgbq(i,j)
            vinfor(nv+11)=vinfor(nv+11)+(1.0-albq(i,j))
     &                   *aire(i,j)*hnbq(i,j)
            vinfor(nv+13)=vinfor(nv+13)
     &                   +(1.0-albq(i,j))*aire(i,j)*hgbq(i,j)*
     &                    (ug(i,j)*ug(i,j)+vg(i,j)*vg(i,j))
          endif
 830    continue
 840  continue
      vinfor(nv+13)=sqrt(vinfor(nv+13)/max(vinfor(nv+9),1e-5*one))
      nv=nv+13
 
C     write(99,*) nv, nvinfo
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Calcul separe  Flux vers Nord / Flux vers Sud .
c-----------------------------------------------------------------------
      nv=nv+1
      j = 56
      ii1 = 107
      ii1 = 106
      ii2 = 108
CL30  j = 109
CL30  ii1 = 209
CL30  ii2 = 215
       vbord = 1.0+(1.0-bound)
       do 850 i=ii1-1,ii2
          ip1    = i+1
C         vfram  = 0.5*(vg(i,j)*tmu(i,j,ks2)+vg(ip1,j)*tmu(ip1,j,ks2))
          vfram  = (vg(i,j)*tmu(i,j,ks2)+vg(ip1,j)*tmu(ip1,j,ks2))/
     &             (max(tmu(i,j,ks2)+tmu(ip1,j,ks2),vbord))
          jj     = j-max(0,int(sign(one,vfram)))
          vinfor(nv) = vinfor(nv)+vfram*dxs1(i,j-1)*1e-6*
     &              (hnbq(i,jj)*0.33+hgbq(i,jj)*0.9)*(1.0-albq(i,jj))
 850   continue
c-----
       nv=nv+1
      j = 55
      ii1 = 109
      ii2 = 114
CL30  j = 108
CL30  ii1 = 216
CL30  ii2 = 225
       do 860 i=ii1-1,ii2
          ip1    = i+1
          vfram  = (vg(i,j)*tmu(i,j,ks2)+vg(ip1,j)*tmu(ip1,j,ks2))/
     &             (max(tmu(i,j,ks2)+tmu(ip1,j,ks2),vbord))
          jj     = j-max(0,int(sign(one,vfram)))
          vinfor(nv) = vinfor(nv)+vfram*dxs1(i,j-1)*1e-6*
     &              (hnbq(i,jj)*0.33+hgbq(i,jj)*0.9)*(1.0-albq(i,jj))
 860   continue
c-----
       nv=nv+1
C      vber=vg(iberp,jberp)*tmu(iberp,jberp,ku2)/2*1e-6
       vber=vg(iberp,jberp)*tmu(iberp,jberp,ku2)/vbord*1e-6
       jj     = jberp-max(0,int(sign(one,vber)))
       vinfor(nv)=vber*dxs1(iberp,jberp-1)*hgbq(iberp,jj)
     &                                *(1.0-albq(iberp,jj))+
     &            vber*dxs1(iberp-1,jberp-1)*hgbq(iberp-1,jj)
     &                                *(1.0-albq(iberp-1,jj))
 
c---Computation of thermal expansion of the ocean
       nv=nv+1
       therma=0.0
       do 870 k=ks1,ks2
        do j=js1,js2
         do i=is1(j),is2(j)
          therma=therma+dz(k)*ctmi(i,j,k,0)*b(i,j,k)
         enddo
       enddo
870    continue
       vinfor(nv)=(-1)*therma*zsurf/gpes
       thex=vinfor(nv)

c---iceshelf and icebergs
c-driess:LOVECLIM includes icebergs from 0, excess of snow and
c-correction linked to the temperature anomalies
       nv=nv+1
       vinfor(nv)=toticesm
       write(99,*) 'informe',ficebergn,ficebergs
       write(99,*) 'informe2',ficebergn/xlg*360.,ficebergs/xlg*360.
       nv=nv+1
       vinfor(nv)=ficebergn/xlg*360.
       nv=nv+1
       vinfor(nv)=ficebergs/xlg*360.
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        if (nv.ne.nvinfo) then
          write(iuo+66,'(2(A,I4))') 'STOP in "informe" : compute nv=', nv,
     &     ' var. <> nvinfo=', nvinfo
          stop
        endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Accumulation before averaging                                   |
c-----------------------------------------------------------------------
 
      do 880 nv=2,nvinfo
        vinfom(nv) = vinfom(nv) + vinfor(nv)
 880  continue
 
      if (.not.flgout) return
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Ecriture sur le fichier type ".evol" :                          |
c-----------------------------------------------------------------------
 
c--Comptabilise le Nb d'iter. cumulees dans "vinfom" :
      nninfo = abs(ninfo)
      nniter = numit - nstart + 1
c- facteur utilise pour Var. Non moyennee : zitsum(0)
      zitsum(0) = 1.
c- facteur utilise pour Var. Moyennee sur tts iter. : zitsum(1)
      if (nniter.eq.0) then
        zitsum(1) = 0.
      elseif (nniter.lt.nninfo) then
c- cas de moyenne sur moins de "ninfo" iter. (en debut de run) :
        zitsum(1) = 1. /  DFLOAT(nniter)
      else
c- cas de moyenne sur "ninfo" iter. (cas standard) :
        zitsum(1) = 1. / DFLOAT(nninfo)
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--Moyenne et transfert : vinfom -> vinfor ; reset "vinfom" :
C     do 890 nv=2,nvinfo
C       vinfor(nv) = vinfom(nv) * zitsum(ktsum(nv))
C       vinfom(nv) = 0.0
C890  continue
      do 891 nv=2,nm1n-1
        vinfor(nv) = vinfom(nv) * zitsum(ktsum(nv))
        vinfom(nv) = 0.0
 891  continue
      do 892 nv=nm2n,nvinfo
        vinfor(nv) = vinfom(nv) * zitsum(ktsum(nv))
        vinfom(nv) = 0.0
 892  continue
      do 893 nv=nm1n,nm2n-1
        vinfom(nv) = 0.0
 893  continue
        
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c- Ecriture sur fichier :
      write(90,fmtw) (titvar(nv),vinfor(nv),nv=1,nvinfo)
 
c- Fermeture du fichier :
      if(numit.eq.nferme) then
        write(90,*)
        close(90)
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  6 ) Traitement pour d'autres sorties et/ou d'autres routines.       |
c-----------------------------------------------------------------------
 
      if (koutpu.le.1 .or. mod(numit,nsav).ne.0) then
c--Remise a zero de fss, Energ. Aj.Conv & Freq.Aj.Conv :
        do 910 k=1,kmax
         do 910 j=1,jmax
          do 910 i=1,imax
            fqajc(i,j,k) = 0.
 910    continue
        do 920 ns=0,nsmax
         do 920 j=1,jmax
          do 920 i=1,imax
           fss(i,j,ns) = 0.
 920    continue
      else
c--Preparation des tableaux avant ecriture sur fichier de resultats :
        do 930 k=ks1,ks2
         do 930 j=js1,js2
          do 930 i=is1(j),is2(j)
           fqajc(i,j,k) = zitsum(1) * fqajc(i,j,k)
 930    continue
c- compute time average surface Fluxes :
        convfx = zitsum(1)
        do 940 ns=0,nsmax
         if (ns.eq.1) convfx = zitsum(1) * ( dz(ks2) / dts(ks2) )
         do 940 j=js1,js2
          do 940 i=is1(j),is2(j)
           fss(i,j,ns) = convfx * fss(i,j,ns)
 940    continue
        call raccord(fqajc(1,1,1), zero, kmax, 4)
        call raccord(fss(1,1,0), zero, 3, 0)
      endif
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine informe -
      end
