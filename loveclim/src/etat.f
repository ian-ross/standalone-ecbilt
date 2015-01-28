












      subroutine etat(mixage, nn99)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--CALCUL DE LA POUSSEE ET CALCUL DE N2 (=bvf) -- Ajust.Conv. selon "mixage"
c   b = gpes * (rho - rhoref) / rho0 ; rhoref = rho(Tref,Sref,z)
c     rho0 = 1030 kg/m3 ; gpes = 9.80 (common static)
c     rho est donne par la formule de C.Eckart (unite : kg/l)
c        (Amerrican Journal of Science, Vol256, April 58, pp225-240)
c Teste de Comparaison : utilise pour g & rho0 & P(z) :
c rho0 = 1.03 kg/l | gpes = 9.8 | P(atm) = z * gpes * rho0 / 101325 .
c -- Calcul de N2 (=Br.Vais.Freq)(bvf) suivant la formule numerique :
c BVF(k-1/2) = [b(T(k-1),S(k-1),z(k-1/2)) - b(T(k),S(k),z(k-1/2))] / dz(k-1/2)
c -- Integration de rho.g.dz par boite pesante.
c  effectue les raccords de grille pour les tableaux w, scal, q et bvf
c -- mixage :  dirige Ajust.Conv. : 0 : rien ;
c  1 : Suppr. TOUTES les instab. ; 2 : melange par paires avec lstab passages
c  3 : Permute et melange partiel (=ajcmix)
c Tient compte de l'Accelerateur de Convergence dans l'Ajust.Convectif.
c  modif : 22/07/97  --  Version Type 7  / Celcius .  --
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'densit.com'
      include 'comunit.h'
 
c--variables locales :
      dimension qintf(imax)
      dimension bup(imax,kmax), bdw(imax,kmax)
      dimension cfb1(imax,kmax), cfb2(imax,kmax)
      dimension rhoref(kmax)
c- pour Ajust. Conv. (tous, 1, 2, 3) :
      dimension scalmi(nsmax), bbajc(kmax), mixtab(kmax)
      dimension ffajc(imax,kmax), ffiajc(kmax)
 
c--Constantes de l'equation d'etat :
      dimension eckart(0:10)
      data eckart /
     &             0.698d0  , 5890.0d0 ,   38.0d0 , 0.375d0 ,  3.0d0  ,
     &             1779.5d0 ,  11.25d0 , 0.0745d0 ,   3.8d0 , 0.01d0  ,
     &           101325.0d0 /
 
      data tkc / 273.15d0 /
Cic0  data tkc / 0.d0 /
 
      data tref, sref / 4.d0 , 34.0d0 /
 
c--Pour l'ecriture en Temp.Potentielle :
      data tavr, savr, cksi / 2.d0 , 34.7d0 , 0.9d-4 /
c- cksi = (T - teta) / |z|, evalue pour S= 35, T=2 , |z|=5km, formule du Gill.
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Mise en place des Coeffs. a la 1ere Iter :                      |
c-----------------------------------------------------------------------
 
      if (numit.le.nstart) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C       write(iuo+66,*) 'etat : version 7  K , I-Boite '
        write(iuo+66,*) 'etat : version 7  C , I-Boite '
     &          //', (lstab,mixage) = ', lstab, mixage
 
c---------
c--1.1 Equation d'etat en Temperature Potentielle :
      do 110 n=0,10
        cstrho(n) = eckart(n)
 110  continue
      gamma = cksi * eckart(10) / ( gpes * rho0 )
      ccpotp = ( eckart(2) - 2.0*eckart(3)*tavr ) * gamma
      unscpp = 1.0 / (1.0 + ccpotp)
      ccpotr = ( eckart(6) - 2.0*eckart(7)*tavr - eckart(9)*savr )
     &         * gamma * unscpp
      cstrho(0)  = cstrho(0) + ccpotr
      cstrho(5)  = cstrho(5) - ccpotr * cstrho(1)
      cstrho(6)  = cstrho(6) - ccpotr * cstrho(2)
      cstrho(7)  = cstrho(7) - ccpotr * cstrho(3)
      cstrho(8)  = cstrho(8) + ccpotr * cstrho(4)
      cstrho(10) = cstrho(10) * unscpp
 
c---------
c--1.2 Mise en place des coefficients dependant de la profondeur.
 
c- pour l'Equ. d'etat (cfb*k*) :
      gravit = gpes * 1000. / rho0
      ccpz = gpes * rho0 / cstrho(10)
      cstrho(1) = cstrho(1) + 1.
      do 140 k=1,kmax
        cfb1z(k)  = cstrho(1) - ccpz * z(k)
 140  continue
      do 150 k=1,kmax+1
        cfb1z4(k) = cstrho(1) - ccpz * zw(k)
 150  continue
 
c---------
c--mise en place de la "densite" de reference ; coeff bilan Ajust.Conv. :
          ccb1 = cstrho(4)*sref
     &         + (cstrho(2)-cstrho(3)*tref)*tref
          ccb2 = cstrho(5)
     &         + (cstrho(6) - cstrho(7)*tref)*tref
     &         - (cstrho(8) + cstrho(9)*tref)*sref
C     ccajc = rho0 / ( DFLOAT(abs(ninfo)) * dts(kmax) )
      ccajc = rho0 / dts(kmax)
      zwtau(1+kmax) = zw(1+kmax)
      do 170 k=kmax,1,-1
        rhoref(k) = 1./(cstrho(0)+ ccb2/(ccb1+cfb1z(k)) )
        bref(k) = gravit * rhoref(k)
        dztau(k) = dz(k) * dts(ks2) / dts(k)
        zwtau(k) = zwtau(k+1) - dztau(k)
        rho0dz(k) = ccajc * dztau(k)
 170  continue
 
c- pour le mixage par paire (m2) :
      do 175 k=2,kmax
        cfm2up(k) = dztau(k)   / (dztau(k-1) + dztau(k))
        cfm2dw(k) = dztau(k-1) / (dztau(k-1) + dztau(k))
        zmix = 0.5 * (zw(k-1) + zw(k+1))
        rhozdz(0,k) = (z(k)   - zmix) * rho0dz(k)
        rhozdz(1,k) = (z(k-1) - zmix) * rho0dz(k-1)
 175  continue
 
c- pour l'Ajust.Conv. par permutation :
      do 180 k2=1,kmax
       do 180 k1=1,kmax
C       dzsdz(k1,k2) = (dz(k1) * unsdz(k2)) * (dts(k2) / dts(k1))
        dzsdz(k1,k2) = dztau(k1) / dztau(k2)
 180  continue
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      if (nn99.eq.2) then
c--Ecriture de controle sur "mouchard" :
        write(99,'(A,3F18.14)') 'tavr,savr, cksi :', tavr, savr, cksi
        write(99,'(A,3F18.14)') 'ccpotp,ccpotr :', ccpotp, ccpotr
        write(99,*) 'Coeff Eq.d''Etat : Eckart | en Temp.Pot :'
        write(99,'(1P2E17.10)') (eckart(k),cstrho(k),k=0,10)
        write(99,'(A,3F18.14)') 'tref,sref, rhoref(tref,sref,k) :',
     &                           tref, sref
        write(99,'(1P5E16.9)') (rhoref(k),k=1,kmax)
        write(99,*)
        write(99,*) 'Coeff dztau(k)=dz(k)*dt(surf)/dt(k) :'
        write(99,'(5F16.9)') (dztau(k),k=1,kmax)
        write(99,*)
      endif
 
      kcrois = 1
      do 190 k=2,kmax
        if (dztau(k).gt.dztau(k-1)) kcrois = k
 190  continue
      if (kcrois.ne.1.and.lstab.eq.-3) then
        write(iuo+66,'(A,I3,A)') 'Ajust.Conv(-3=lstab) & Acc.Conv(k=',
     &                      kcrois, ' ) INCOMPATIBLE !'
        write(iuo+66,'(A)') ' => ARRET , routine "etat" '
        stop
      endif
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
c--Initialisation :
      ninstb = 0
 
c--Traitement differencie selon le type d'ajustement convectif :
 
      if (mixage.eq.1) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Suppression de TOUTES les instabilites (lstab = -1) :           |
c-----------------------------------------------------------------------
 
c--convention : mixtab(k) : plus bas niveau melange avec la boite "k"
c               kdwmix : plus bas  niveau a melanger avec la boite courante.
c               kupmix : plus haut niveau a melanger avec la boite courante.
c               ktopmi : plus haut niveau melange de toute la colonne
 
c--Debut de la boucle externe sur l'indice de latitude j :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,k,kins,ns,kloc,tloc,sloc,ccb1,ccb2,bup,bdw)
C$DIR SPP LOOP_PRIVATE(mixtab,kdwmix,kupmix,ktopmi,zmix,unshmi)
C$DIR SPP LOOP_PRIVATE(scalmi,bdwmix,bbsav,bbajc,kbotmi,qintf,dqd2)
      do 300 j=js1,js2
c-----
 
c--2.1 Calcul de b(i,j,k) dans tout le domaine :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      do 210 k=ks1,ks2
        do 210 i=is1(j),is2(j)
c- conversion Degr.Celcius et tranfert -> variables locales (tloc,sloc) :
          kloc = max(k,kfs(i,j))
          tloc = scal(i,j,kloc,1) - tkc
          sloc = scal(i,j,kloc,2)
          ccb1 = cstrho(4)*sloc
     &         + (cstrho(2)-cstrho(3)*tloc)*tloc
          ccb2 = cstrho(5)
     &         + (cstrho(6) - cstrho(7)*tloc)*tloc
     &         - (cstrho(8) + cstrho(9)*tloc)*sloc
          b(i,j,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z(k)) ) - bref(k)
          bup(i,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z4(k+1)) )
          bdw(i,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z4( k )) )
 210  continue
 
      do 270 i=is1(j),is2(j)
c- boucle NON VECTORISABLE sur l'indice de longitude "i" :
 
c--2.2 Detecte la 1ere instabilite :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
        do 225 kins=kfs(i,j)+1,ks2
c- teste la stabilite de kins / kins-1 :
          if (bdw(i,kins).gt.bup(i,kins-1)) then
c- initialisation :
            do 220 k=ks1,ks2
              mixtab(k) = k
 220        continue
            kdwmix = kins - 1
            mixtab(kins) = kdwmix
            kupmix = kins
            ktopmi = kins
            fqajc(i,j,kins) = fqajc(i,j,kins) + 1.
            goto 230
          endif
 225    continue
        goto 270
 
 230    continue
 
c--2.3 mixing des boites depuis kdwmix (=mixtab(kins)) --> kupmix
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
        zmix   = 0.5 * ( zw(kupmix+1) + zw(kdwmix) )
        unshmi = 1.0 / ( zwtau(kupmix+1) - zwtau(kdwmix) )
        do 245 ns=1,nsmax
          scalmi(ns) = 0.
          do 235 k=kdwmix,kupmix
            scalmi(ns) = scalmi(ns) + dztau(k)*scal(i,j,k,ns)
 235      continue
          scalmi(ns) = scalmi(ns) * unshmi
          do 240 k=kdwmix,kupmix
            scal(i,j,k,ns) = scalmi(ns)
 240      continue
 245    continue
c- fin du mixing .
 
c--2.4 calcule les nouveaux "b" et fait le bilan de l'ajustement convectif :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
          tloc = scal(i,j,kins,1) - tkc
          sloc = scal(i,j,kins,2)
          ccb1 = cstrho(4)*sloc
     &         + (cstrho(2)-cstrho(3)*tloc)*tloc
          ccb2 = cstrho(5)
     &         + (cstrho(6) - cstrho(7)*tloc)*tloc
     &         - (cstrho(8) + cstrho(9)*tloc)*sloc
        bdwmix = gravit / (cstrho(0) + ccb2/(ccb1+cfb1z4(kdwmix)) )
        do 250 k=kdwmix,kupmix
          bbsav = b(i,j,k)
          bup(i,k) = bdwmix
          bdw(i,k) = bdwmix
          b(i,j,k) = gravit /
     &               (cstrho(0)+ ccb2/(ccb1+cfb1z(k)) ) - bref(k)
          bbajc(k) = (bbsav-b(i,j,k)) * (z(k)-zmix) * rho0dz(k)
          fqajc(i,j,1) = fqajc(i,j,1) + bbajc(k)
 250    continue
        bup(i,kupmix) = gravit /
     &                  (cstrho(0) + ccb2/(ccb1+cfb1z4(kupmix+1)) )
        if ( kdwmix.eq.kfs(i,j) .and. kdwmix.gt.ks1 ) then
          do 255 k=ks1,kdwmix-1
            b(i,j,k) = gravit /
     &               (cstrho(0)+ ccb2/(ccb1+cfb1z(k)) ) - bref(k)
 255      continue
        endif
C       if (kupmix.ge.(ks2-1)) hmajc(i,j) = max(hmajc(i,j),-z(kdwmix))
C       ninstb = ninstb + kupmix - kdwmix
 
c--2.5 Test de stab. - determine kdwmix & kupmix :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
        kbotmi = max(kfs(i,j)+1,kdwmix)
        do 265 kins=kbotmi,ks2
c- teste la stabilite de kins / kins-1 :
          if (mixtab(kins).eq.kins .and.
     &        bdw(i,kins).gt.bup(i,kins-1)) then
            kdwmix = mixtab(kins-1)
            mixtab(kins) = kdwmix
            kupmix = kins
            fqajc(i,j,kins) = fqajc(i,j,kins) + 1.
c- Jonction avec une zone (initialement instable) deja melangee :
            do 260 k=kins+1,ktopmi
              if (mixtab(k).eq.kins) then
                kupmix = k
                mixtab(k) = kdwmix
              endif
 260        continue
            ktopmi = kupmix
            goto 230
          endif
 265    continue
 
c--fin de la boucle sur "i" .
 270  continue
 
c--Fin du bloc 'suppression de TOUTES les instabilites'.
 
c--2.6 Calcul de la pression reduite q :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--calcul de la pression (q) (transfere depuis 'uve')
c-   (integration de rho.g par la Methode des Trapezes) :
C
C     do 280 i=is1(j),is2(j)
C       q(i,j,ks2)= 0.5*dz(ks2)*b(i,j,ks2)
C280  continue
C
C     do 285 k=ks2-1,ks1,-1
C      do 285 i=is1(j),is2(j)
C       q(i,j,k)=q(i,j,k+1)+0.5*dz4(k+1)*(b(i,j,k+1)+b(i,j,k))
C285  continue
C
c---------------------
c--calcul de la pression (q) , integration de rho.g "par boite pesante" :
 
      do 280 i=1,imax
        qintf(i) = 0.0
 280  continue
 
      do 285 k=ks2,ks1,-1
       do 285 i=is1(j),is2(j)
         dqd2 = 0.5*dz(k)*b(i,j,k)
         q(i,j,k) = qintf(i) + dqd2
         qintf(i) = q(i,j,k) + dqd2
 285  continue
 
 
c--2.7 calcul de N2 (transfere depuis 'stab') :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      do 290 k=ks1+1,ks2
       do 290 i=is1(j),is2(j)
        bvf(i,j,k) = tms(i,j,k-1) * unsdzw(k) * (bup(i,k-1)-bdw(i,k))
 290  continue
 
c--Fin de la boucle externe sur l'indice de latitude j .
 300  continue
 
      elseif(mixage.eq.2) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Suppression des instabilites : melange par paires (lstab > 0) . |
c-----------------------------------------------------------------------
 
c--Debut de la boucle externe sur l'indice de latitude j :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,k,kloc,tloc,sloc,ccb1,ccb2,bup,bdw)
C$DIR SPP LOOP_PRIVATE(ns,nstab,kk,nk,kk1,kk2,kkd)
C$DIR SPP LOOP_PRIVATE(scalmi,tmix,smix,bbajc,ffajc,qintf,dqd2)
      do 400 j=js1,js2
c-----
 
c--3.1 Calcul de b(i,j,k) dans tout le domaine :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      do 310 k=ks1,ks2
        do 310 i=is1(j),is2(j)
c- conversion Degr.Celcius et tranfert -> variables locales (tloc,sloc) :
          kloc = max(k,kfs(i,j))
          tloc = scal(i,j,kloc,1) - tkc
          sloc = scal(i,j,kloc,2)
          ccb1 = cstrho(4)*sloc
     &              + (cstrho(2)-cstrho(3)*tloc)*tloc
          ccb2 = cstrho(5)
     &              + (cstrho(6) - cstrho(7)*tloc)*tloc
     &              - (cstrho(8) + cstrho(9)*tloc)*sloc
          b(i,j,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z(k)) ) - bref(k)
          bup(i,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z4(k+1)) )
          bdw(i,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z4( k )) )
          ffajc(i,k) = 1.
 310  continue
 
c--3.2 Debut du traitement de l'Ajust.Convectif : boites a melanger ?
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c- le nombre de passages = lstab
      do 370 nstab=1,lstab
       kk2 = ks2 - mod(numit+nstab,2)
c--debut boucle NON VECTORISABLE sur l'indice de longitude "i" :
       do 360 i=is1(j),is2(j)
 
        kk1 = kfs(i,j) + 1
        do 350 kk=kk2,kk1,-2
c--3.2 Detecte les boites a melanger : teste la stabilite de kk / kk-1 :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
          if (bdw(i,kk).gt.bup(i,kk-1)) then
c--3.3 mixing des boites   kk / kk-1 :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
            kkd = kk - 1
            do 320 ns=1,nsmax
              scalmi(ns) = cfm2up(kk) * scal(i,j,kk,ns)
     &                   + cfm2dw(kk) * scal(i,j,kkd,ns)
 320        continue
c--3.4 Mise en place des nouveaux scalaires et calcule des nouveaux "b" :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
            tmix = scalmi(1) - tkc
            smix = scalmi(2)
            ccb1 = cstrho(4)*smix
     &           + (cstrho(2)-cstrho(3)*tmix)*tmix
            ccb2 = cstrho(5)
     &           + (cstrho(6) - cstrho(7)*tmix)*tmix
     &           - (cstrho(8) + cstrho(9)*tmix)*smix
            do 335 nk=0,1
             k = kk - nk
             do 330 ns=1,nsmax
               scal(i,j,k,ns) = scalmi(ns)
 330         continue
             bbajc(k) = b(i,j,k)
             bdw(i,k) = gravit / (cstrho(0) + ccb2/(ccb1+cfb1z4(k)))
             b(i,j,k) = gravit / (cstrho(0)+ ccb2/(ccb1+cfb1z(k)) )
     &                 - bref(k)
             bbajc(k) = rhozdz(nk,kk) * (bbajc(k) - b(i,j,k))
 335        continue
            bup(i,kkd) =  bdw(i,kk)
            bup(i,kk) = gravit / (cstrho(0) + ccb2/(ccb1+cfb1z4(kk+1)))
            if ( kkd.eq.kfs(i,j) .and. kkd.gt.ks1) then
              do 340 k=ks1,kkd-1
               b(i,j,k) = gravit / (cstrho(0)+ ccb2/(ccb1+cfb1z(k)) )
     &                  - bref(k)
 340          continue
            endif
c--3.5 bilan de l'ajustement convectif :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
            fqajc(i,j,1) = fqajc(i,j,1) + bbajc(kkd) + bbajc(kk)
            fqajc(i,j,kk) = fqajc(i,j,kk) + ffajc(i,kk)
            ffajc(i,kk) = 0.
C           hmajc(i,j) = max(hmajc(i,j),-z(kk))
C           if ( hmajc(i,j).ge.-zw(kk+1) )
C    &           hmajc(i,j) = max(hmajc(i,j),-z(kkd))
C           ninstb = ninstb + 1
          endif
c- fin du mixing .
 350    continue
 
 360   continue
c--fin de la boucle sur "i" .
 370  continue
c--fin de la boucle sur "nstab" (<--nombre de passage) .
 
c--Fin du bloc 'suppression des instabilites par paires'.
 
c--3.6 Calcul de la pression reduite q :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--calcul de la pression (q) (transfere depuis 'uve')
c-   (integration de rho.g par la Methode des Trapezes) :
C
C     do 380 i=is1(j),is2(j)
C       q(i,j,ks2)= 0.5*dz(ks2)*b(i,j,ks2)
C380  continue
C
C     do 385 k=ks2-1,ks1,-1
C      do 385 i=is1(j),is2(j)
C       q(i,j,k)=q(i,j,k+1)+0.5*dz4(k+1)*(b(i,j,k+1)+b(i,j,k))
C385  continue
C
c---------------------
c--calcul de la pression (q) , integration de rho.g "par boite pesante" :
 
      do 380 i=1,imax
        qintf(i) = 0.0
 380  continue
 
      do 385 k=ks2,ks1,-1
       do 385 i=is1(j),is2(j)
        dqd2 = 0.5*dz(k)*b(i,j,k)
        q(i,j,k) = qintf(i) + dqd2
        qintf(i) = q(i,j,k) + dqd2
 385  continue
 
 
c--3.7 calcul de N2 (transfere depuis 'stab') :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      do 390 k=ks1+1,ks2
       do 390 i=is1(j),is2(j)
        bvf(i,j,k) = tms(i,j,k-1) * unsdzw(k) * (bup(i,k-1)-bdw(i,k))
 390  continue
 
c--Fin de la boucle externe sur l'indice de latitude j .
 400  continue
 
      elseif (mixage.eq.3) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Suppression des instabilites par permutation (lstab = -3) :     |
c-----------------------------------------------------------------------
 
c--convention : kins : niveau (= interface) instable.
c               kdwm : plus bas atteint par l'Ajust. Conv.
c-----
      ajc0mx = 1.0 - ajcmix
 
c--Debut de la boucle externe sur l'indice de latitude j :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,k,kins,ns,kloc,tloc,sloc,cfb1,cfb2,bup,bdw)
C$DIR SPP LOOP_PRIVATE(kdwm,bdwins,ajchmx,sscdwm,ddssc,scalmi)
C$DIR SPP LOOP_PRIVATE(zmix,bbsav,bbajc,ffiajc,qintf,dqd2)
      do 500 j=js1,js2
c-----
 
c--4.1 Calcul de b(i,j,k) dans tout le domaine :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      do 410 k=ks1,ks2
        do 410 i=is1(j),is2(j)
c- conversion Degr.Celcius et tranfert -> variables locales (tloc,sloc) :
          kloc = max(k,kfs(i,j))
          tloc = scal(i,j,kloc,1) - tkc
          sloc = scal(i,j,kloc,2)
          cfb1(i,k) = cstrho(4)*sloc
     &              + (cstrho(2)-cstrho(3)*tloc)*tloc
          cfb2(i,k) = cstrho(5)
     &              + (cstrho(6) - cstrho(7)*tloc)*tloc
     &              - (cstrho(8) + cstrho(9)*tloc)*sloc
          b(i,j,k) = gravit /
     &      (cstrho(0)+ cfb2(i,k)/(cfb1(i,k)+cfb1z(k)) ) - bref(k)
          bup(i,k) = gravit /
     &      (cstrho(0)+ cfb2(i,k)/(cfb1(i,k)+cfb1z4(k+1)) )
          bdw(i,k) = gravit /
     &      (cstrho(0)+ cfb2(i,k)/(cfb1(i,k)+cfb1z4( k )) )
 410  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      do 470 i=is1(j),is2(j)
c--debut boucle NON VECTORISABLE sur l'indice de longitude "i" :
 
       do 460 kins=kfs(i,j)+1,ks2
         ffiajc(kins) = 1.
       if (bdw(i,kins).gt.bup(i,kins-1))  then
c- debut du traitemnent de l'instabilite de kins / kins-1 :
 
c--bilan de l'Ajust.Conv : Interface (k/k-1) concernees par le melange :
         fqajc(i,j,kins) = fqajc(i,j,kins) + 1.
         ffiajc(kins) = 0.
 
c--4.2 repere les boites a permuter et a melanger partiellement  :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Recherche de "kdwm" tel que bup(kdwm-1) > b(Tkins,Skins,zw(kdwm))
        kdwm = kfs(i,j)
        do 420 k=kins-1,kfs(i,j)+1,-1
          bdwins = gravit /
     &       (cstrho(0)+ cfb2(i,kins)/(cfb1(i,kins)+cfb1z4(k)) )
          if (bdwins.le.bup(i,k-1) ) then
            kdwm = k
            goto 425
          endif
          fqajc(i,j,k) = fqajc(i,j,k) + ffiajc(k)
          ffiajc(k) = 0.
 420    continue
 425    continue
 
c--4.3 Permutation circulaire de H=dztau(kins) + melange kdwm --> kins :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
        ajchmx = ajcmix / ( zwtau(kins+1) - zwtau(kdwm) )
        do 435 ns=1,nsmax
c- Permute & calcule Moy. & New Value :
          sscdwm = scal(i,j,kins,ns)
          scalmi(ns) = dztau(kins) * sscdwm
          scal(i,j,kins,ns) = scal(i,j,kins-1,ns)
          do 430 k=kdwm,kins-1
            ddssc = dzsdz(kins,k) * (sscdwm - scal(i,j,k,ns))
            sscdwm = scal(i,j,k,ns)
            scalmi(ns) = scalmi(ns) + dztau(k) * sscdwm
            scal(i,j,k,ns) = scal(i,j,k,ns) + ddssc
 430      continue
          scalmi(ns) = scalmi(ns) * ajchmx
c- ajoute le melange partiel(=ajcmix) depuis kdwm -> kins :
          do 432 k=kdwm,kins
            scal(i,j,k,ns) = ajc0mx * scal(i,j,k,ns) + scalmi(ns)
 432      continue
 435    continue
c- fin du mixing .
 
c--4.4 calcule les nouveaux "b" et fait le bilan de l'ajustement convectif :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
        zmix   = 0.5 * ( zw(kins+1) + zw(kdwm) )
        do 440 k=kdwm,kins
          tloc = scal(i,j,k,1) - tkc
          sloc = scal(i,j,k,2)
          bbajc(k) = b(i,j,k)
          cfb1(i,k) = cstrho(4)*sloc
     &              + (cstrho(2)-cstrho(3)*tloc)*tloc
          cfb2(i,k) = cstrho(5)
     &              + (cstrho(6) - cstrho(7)*tloc)*tloc
     &              - (cstrho(8) + cstrho(9)*tloc)*sloc
          b(i,j,k) = gravit /
     &      (cstrho(0)+ cfb2(i,k)/(cfb1(i,k)+cfb1z(k)) ) - bref(k)
          bup(i,k) = gravit /
     &      (cstrho(0)+ cfb2(i,k)/(cfb1(i,k)+cfb1z4(k+1)) )
          bdw(i,k) = gravit /
     &      (cstrho(0)+ cfb2(i,k)/(cfb1(i,k)+cfb1z4( k )) )
          bbajc(k) = (bbajc(k)-b(i,j,k)) * (z(k)-zmix) * rho0dz(k)
          fqajc(i,j,1) = fqajc(i,j,1) + bbajc(k)
 440    continue
        if ( kdwm.eq.kfs(i,j) .and. kdwm.gt.ks1 ) then
          do 450 k=ks1,kdwm-1
           b(i,j,k) = gravit /
     &     (cstrho(0)+ cfb2(i,kdwm)/(cfb1(i,kdwm)+cfb1z(k))) - bref(k)
 450     continue
        endif
C       if (kins.ge.(ks2-1)) hmajc(i,j) = max(hmajc(i,j),-z(kdwm))
C       ninstb = ninstb + kins - kdwm
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c- fin du traitemnent de l'instabilite de kins / kins-1 .
       endif
 460   continue
 
c--fin de la boucle sur "i" .
 470  continue
 
c--Fin du bloc 'suppression des instabilites par permutation'.
 
c--4.6 Calcul de la pression reduite q :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--calcul de la pression (q) (transfere depuis 'uve')
c-   (integration de rho.g par la Methode des Trapezes) :
C
C     do 480 i=is1(j),is2(j)
C       q(i,j,ks2)= 0.5*dz(ks2)*b(i,j,ks2)
C480  continue
C
C     do 485 k=ks2-1,ks1,-1
C      do 485 i=is1(j),is2(j)
C       q(i,j,k)=q(i,j,k+1)+0.5*dz4(k+1)*(b(i,j,k+1)+b(i,j,k))
C485  continue
C
c---------------------
c--calcul de la pression (q) , integration de rho.g "par boite pesante" :
 
      do 480 i=1,imax
        qintf(i) = 0.0
 480  continue
      do 485 k=ks2,ks1,-1
       do 485 i=is1(j),is2(j)
        dqd2 = 0.5*dz(k)*b(i,j,k)
        q(i,j,k) = qintf(i) + dqd2
        qintf(i) = q(i,j,k) + dqd2
 485  continue
 
 
c--4.7 calcul de N2 (transfere depuis 'stab') :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      do 490 k=ks1+1,ks2
       do 490 i=is1(j),is2(j)
        bvf(i,j,k) = tms(i,j,k-1) * unsdzw(k) * (bup(i,k-1)-bdw(i,k))
 490  continue
 
c--Fin de la boucle externe sur l'indice de latitude j .
 500  continue
 
      else
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Sans Ajustement Convectif : mixage= 0 (lstab = 0).              |
c-----------------------------------------------------------------------
 
c--Debut de la boucle externe sur l'indice de latitude j :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,k,kloc,tloc,sloc,ccb1,ccb2)
C$DIR SPP LOOP_PRIVATE(bup,bdw,qintf,dqd2)
      do 600 j=js1,js2
c-----
 
c--5.1 Calcul de b(i,j,k) dans tout le domaine :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      do 510 k=ks1,ks2
        do 510 i=is1(j),is2(j)
c- conversion Degr.Celcius et tranfert -> variables locales (tloc,sloc) :
          kloc = max(k,kfs(i,j))
          tloc = scal(i,j,kloc,1) - tkc
          sloc = scal(i,j,kloc,2)
          ccb1 = cstrho(4)*sloc
     &              + (cstrho(2)-cstrho(3)*tloc)*tloc
          ccb2 = cstrho(5)
     &              + (cstrho(6) - cstrho(7)*tloc)*tloc
     &              - (cstrho(8) + cstrho(9)*tloc)*sloc
          b(i,j,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z(k)) ) - bref(k)
          bup(i,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z4(k+1)) )
          bdw(i,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z4( k )) )
 510  continue
 
c--5.6 Calcul de la pression reduite q :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--calcul de la pression (q) (transfere depuis 'uve')
c-   (integration de rho.g par la Methode des Trapezes) :
C
C       do 580 i=is1(j),is2(j)
C        q(i,j,ks2)= 0.5*dz(ks2)*b(i,j,ks2)
C580    continue
C
C       do 585 k=ks2-1,ks1,-1
C        do 585 i=is1(j),is2(j)
C         q(i,j,k)=q(i,j,k+1)+0.5*dz4(k+1)*(b(i,j,k+1)+b(i,j,k))
C585    continue
C
c---------------------
c--calcul de la pression (q) , integration de rho.g "par boite pesante" :
 
        do 580 i=1,imax
         qintf(i) = 0.0
 580    continue
 
        do 585 k=ks2,ks1,-1
         do 585 i=is1(j),is2(j)
          dqd2 = 0.5*dz(k)*b(i,j,k)
          q(i,j,k) = qintf(i) + dqd2
          qintf(i) = q(i,j,k) + dqd2
 585    continue
 
 
c--5.7 calcul de N2 (transfere depuis 'stab') :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      do 590 k=ks1+1,ks2
       do 590 i=is1(j),is2(j)
        bvf(i,j,k) = tms(i,j,k-1) * unsdzw(k) * (bup(i,k-1)-bdw(i,k))
 590  continue
 
c--Fin de la boucle externe sur l'indice de latitude j .
 600  continue
 
c----------------------
 
c--Fin du traitement differencie selon "mixage" .
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  7 ) raccord cyclique + raccord de grille pour w, q, bvf, scal .     |
c-----------------------------------------------------------------------
 
      call raccord(w(1,1,1), zero, kmax, 4)
      call raccord(q(1,1,1), zero, kmax, 0)
      call raccord(b(1,1,1), zero, kmax, 0)
      call raccord(bvf(1,1,1), zero, kmax, 4)
      call raccord(scal(1,1,1,1), zero, kmax*nsmax, 0)
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine etat -
      end
