












      subroutine vdiffu(nn99)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c Calcul des coefficients de DIFFUsion Verticale (divise par dz) : avsdz, avudz
c  Viscosite & Diffusivite Verticale : fonction de la Freq. de Br.Vais. (bvf)
c   suivant la formulation de Pacanowski et Philander (1981)
c Modification : Diffusivite V. Simplifie = avkb + avk0 / (1+5Ri)^3
c   Les seuils portent sur les termes en 1/(1+5Ri)
c--
c Diffusivite Ocean Profond : AvsN2 * (N2)^q  (q = -0.6 environ)
c  Inactif si kavsmx=0 ; Interfaces concernes : k=2,kavsmx
c  bornes (sur Avs) globales (k=1) et par niveau : avsmn & avsmx (lu sur fich)
c  => equivalent sur N2 : bvfmn & bvfmx
c  bvfmn(1) & bvfmx(1) definissent la discretication de la fct Avs(N2)
c   par pas de "dn2" suivant N2 : avsn(nav),nav=0,navmax
c--
c  modif : 11/05/99
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'comunit.h'
 
      parameter ( navmax = 10000 )
c--variables locales conservees d'un appel a l'autre :
      common / vdifloc / epsil2, ref0n2, ref0m2, uv2bnd, cc1zm,
     &  avu0ri(kmax), avu1ri(kmax), avs1ri(kmax),
     &  bvfmix, avs2mx(kmax), unszek(imax,jmax),
     &  ccudm,avs0ri(kmax)
      common / kdifloc / kk0ri
 
      common / vdifbvf / unsdn2,
     &  bvfmn(kmax), bvfmx(kmax), avsn(0:navmax+1)
 
c--variables locales :
      dimension riums(imax,jmax)
      dimension avsmn(kmax), avsmx(kmax), bvfloc(imax,jmax)
 
c--valeurs intervenant dans la formule de Pac&Phi (certaines deja ds le common)
      data alpha / 5.0d0     /
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  0 ) 1ere Iter, Determine 1er niv. a calculer ; Initialisation .     |
c-----------------------------------------------------------------------
 
      if (numit.le.nstart) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Formulation Avs(N2) : mise en place des Coeff.
 
      do 10 k=1,kmax
        avsmn(k) = 1.
        avsmx(k) = 1.
 10   continue
 
c- Lecture des parametres pour la formulation Avs(N2) <- dans defcst
      if (kavsmx.le.0) goto 50
      kavsmx = min(kavsmx,kmax)
 
c- Mise en place des coeff. & Verification :
      unsqav = one / qavs
      do 30 k=1,kavsmx
        avsmn(k) = ccfmn * avsmn(k)
        avsmx(k) = ccfmx * avsmx(k)
        if ( avsmn(k).lt.avsmn(1) .or. avsmn(k).gt.avsmx(k) .or.
     &       avsmx(k).gt.avsmx(1) ) then
          write(iuo+66,*) 'ARRET : vdiffu, coeff. avsmn,avsmx False', k
          stop
        endif
        bvfmn(k) = (avsmx(k)/avsn2)**unsqav
        bvfmx(k) = (avsmn(k)/avsn2)**unsqav
 30   continue
      dn2 = (bvfmx(1) - bvfmn(1)) / DFLOAT(navmax)
      if ( dn2.lt.epsil*epsil ) then
          write(iuo+66,*) 'ARRET : vdiffu, N2 Discr. step too small', dn2
          stop
      endif
      unsdn2 = one / dn2
 
c-  Mise en place de la Fonct. Tabulee : Avs(N2)
      do 40 nav=0,navmax+1
        avsn(nav) = avsn2 * ( bvfmn(1) + dn2 * DFLOAT(nav) )**qavs
 40   continue
 
      if (nn99.eq.2) then
c--Impression sur fichier mouchard :
        write(99,'(2A,I4,1P2E16.8)') 'vdiffu, Avs(N2) :',
     &     ' kavsmx,qavs,avsn2 =', kavsmx, qavs, avsn2
        write(99,*) 'navmax(=N), dn2 =', navmax, dn2
        nn = navmax / 2
        write(99,'(A,1P4E14.6)') ' avsn(0,1,N/2,N) :',
     &       avsn(0), avsn(1), avsn(nn), avsn(navmax)
        write(99,'(A,(1P6E12.5))') 'avsmn :', (avsmn(k),k=1,kavsmx)
        write(99,'(A,(1P6E12.5))') 'avsmx :', (avsmx(k),k=1,kavsmx)
        write(99,'(A,(1P6E12.5))') 'bvfmn :', (bvfmn(k),k=1,kavsmx)
        write(99,'(A,(1P6E12.5))') 'bvfmx :', (bvfmx(k),k=1,kavsmx)
        write(99,*)
      endif
 
 50   continue
c- fin du bloc de definition des Coeff. pour formulation Avs(N2).
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      epsil2 = epsil * epsil
c- valeur typique de N2 : 1.E-5 (s-2)
      ref0n2 = 1.d-5 * epsil
 
      avsmix = avnub(1)
      bvfmix = avnu0(1)
      ref0m2 = avnub(1)
      uvzbnd = avnu0(1)
      zzekm = avk0(1)
      ccudm = avkb(1)
 
      xx = uvzbnd * 100.0
        write(iuo+66,'(2A,F6.2,A1,I4,A1)')
C    &   ' vdiffu : Pac.&Phil. Ris,Riu separes'
C    &   ' vdiffu : Pac.&Phil. avu=f(min{N2}) '
C    &   ' vdiffu : Pac.&Phil. Ris = Moy.Riu  '
     &   ' vdiffu : Pac.&Phil. avu=f(min.N2) & Ris = Moy.Riu'
     &  ,', MinM2 :', ccudm, ',',  nint(zzekm), 'm'
C    &  ,'<-MinM2 :', ccudm, ',',  nint(zzekm), 'm'
C       write(iuo+66,'(2A,1PE8.1,A,I2,2A,0PF4.1,A,I4,A)')
C    &   ' vdiffu : Pac&Phil, avu=f(min.N2), Ris=4.Riu'
C    &   ' vdiffu : Pac&Phil, avu=f(min.N2), Ris=4.RiU'
C    &  ,', MinM2:', ref0m2, ',', nint(xx), 'cm/s'
C    &  ,',', ccudm, ',', nint(zzekm), 'm'
        if (lstab.eq.0) write(iuo+66,'(2(A,1PE10.3))')
     &   '          AVS +', avsmix, ' si bvf <', bvfmix
        if (kavsmx.ne.0) write(iuo+66,'(3X,A,I4,1P2E11.3)')
     &   ' Avs = AvsN2*(N2)^q ; kmx,AvsN2,q =', kavsmx, avsn2, qavs
 
c--1er Niveau a calculer :
      kk0ri = ks2
      avu0ri(ks1) = 0.0
      avs0ri(ks1) = 0.0
      do 70 k=ks1+1,ks2
        if (kk0ri.eq.ks2 .and. (avk0(k).ne.zero .or. avnu0(k).ne.zero))
     &      kk0ri = k - 1
        avu0ri(k) = unsdzw(k) * avnub(k)
        avu1ri(k) = unsdzw(k) * avnu0(k)
        avs0ri(k) = unsdzw(k) * avkb(k)
        avs1ri(k) = unsdzw(k) * avk0(k)
        avs2mx(k) = unsdzw(k) * avsmix * 0.5
        if (k.le.kavsmx) avs0ri(k) = 0.
 70   continue

      uv2bnd = uvzbnd * uvzbnd
      if (ccudm.ge.epsil) then
c--profondeur pour le calcul du cisaillement minimun :
        cc1zm = ccudm * ccudm * 4.
        cc2zm = cc1zm / ( epsil + zzekm * zzekm )
        cc2zm = cc2zm / (omega * omega)
        do 80 j=1,jmax
         do 80 i=1,imax
          unszek(i,j) = cc2zm * fs2cor(i,j) * fs2cor(i,j)
 80     continue
      endif

c- Fin du traitement specifique de la 1ere Iter.
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
c--Debut de la boucle externe sur l'indice de niveau k :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,j,unszz4,uz,vz,unsm2,cczzm,ud,vd,ccm2,ccm2n)
C$DIR SPP LOOP_PRIVATE(riums,ccdz,tm4u,deno,rifac,unsden)
C$DIR SPP LOOP_PRIVATE(nav,bvfloc,bvfavs,xav)
      do 700 k=ks1+1,ks2
c-----
 
      if (k.le.kk0ri) then
c--Niveaux profond : Avu,Avs independent de Ri :
        do 100 j=js1,js2
         do 100 i=ims1,ims2
          avudz(i,j,k) = avu0ri(k)
          avsdz(i,j,k) = avs0ri(k)
 100    continue
 
      else
c--Niveaux peu profond : Calcul de Avu,Avs en fonction de Ri :
        unszz4 = unsdzw(k) * unsdzw(k) * 4.
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Ridchardson Number at velocity point .                          |
c-----------------------------------------------------------------------
c  avudz <- Ri used for momentum,  [=Min(N2)/M2] or [=Moy(N2)/M2]
c  riums <- Ri used for scalar  ,  [=Moy(N2)/M2]
 
c--Initialisation :
      do 110 j=1,jmax
       do 110 i=1,imax
        riums(i,j) = ref0n2
 110  continue


      if (ccudm.lt.epsil) then
c--Evaluation directe du cissaillement :
        do 150 j=ju1,ju2
         do 150 i=iu1(j),iu2(j)
          uz = u(i,j,k) - u(i,j,k-1)
          vz = v(i,j,k) - v(i,j,k-1)
          unsm2 = (tmu(i,j,k-1) + epsil2)
     &          / (unszz4 * (uz*uz + vz*vz) + epsil)
c- Ridchardson Number at velocity point, for scalar   eq. (riums) :
          riums(i,j) = ( (bvf(i-1,j-1,k) + bvf(i,j,k))
     &                 + (bvf(i-1,j,k) + bvf(i,j-1,k)) ) * unsm2
c- Ridchardson Number at velocity point, for momentum eq. (avudz) :
          avudz(i,j,k) = min(bvf(i-1,j-1,k), bvf(i,j,k),
     &                       bvf(i-1,j,k), bvf(i,j-1,k)) * 4. * unsm2
C         avudz(i,j,k) = riums(i,j)
c-
 150    continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      else
c--Evaluation du cissaillement : Minore par (|V| / max(Dekman,Z))
        cczzm = cc1zm / (zw(k) * zw(k))
        do 170 j=ju1,ju2
         do 170 i=iu1(j),iu2(j)
          uz = u(i,j,k) - u(i,j,k-1)
          vz = v(i,j,k) - v(i,j,k-1)
          ud = u(i,j,k) + u(i,j,k-1)
          vd = v(i,j,k) + v(i,j,k-1)
          ccm2  = (uz*uz + vz*vz) * unszz4 + epsil
          ccm2n = (ud*ud + vd*vd) * min(cczzm, unszek(i,j))
          unsm2 = (tmu(i,j,k-1) + epsil2) / max(ccm2, ccm2n)
c- Ridchardson Number at velocity point, for scalar   eq. (riums) :
          riums(i,j) = ( (bvf(i-1,j-1,k) + bvf(i,j,k))
     &               + (bvf(i-1,j,k) + bvf(i,j-1,k)) ) * unsm2
c- Ridchardson Number at velocity point, for momentum eq. (avudz) :
C         avudz(i,j,k) = ( (bvf(i-1,j-1,k) + bvf(i,j,k))
C    &                   + (bvf(i-1,j,k) + bvf(i,j-1,k)) )
          avudz(i,j,k) = min( bvf(i-1,j-1,k), bvf(i,j,k),
     &                        bvf(i-1,j,k), bvf(i,j-1,k) ) * 4.0
C    &                 * (tmu(i,j,k-1) + epsil2) / ccm2
     &                 * unsm2
C         avudz(i,j,k) = riums(i,j)
c-
 170    continue
      endif

 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Ridchardson Number at scalar point :                            |
c-----------------------------------------------------------------------
 
c--calcul (dans avsdz) du nombre de Ridchardson au point scalaire :
C     unszz  = unsdzw(k) * unsdzw(k)
CC    unsz16 = unsdzw(k) * unsdzw(k) * 0.0625
c-  Weat Point Only :
C     do 210 j=js1,js2
C      do 210 i=is1(j),is2(j)
C       uz = ( (u(i,j,k)-u(i,j,k-1)) + (u(i+1,j+1,k)-u(i+1,j+1,k-1)) )
C    &     + ( (u(i,j+1,k)-u(i,j+1,k-1)) + (u(i+1,j,k)-u(i+1,j,k-1)) )
C       vz = ( (v(i,j,k)-v(i,j,k-1)) + (v(i+1,j+1,k)-v(i+1,j+1,k-1)) )
C    &     + ( (v(i,j+1,k)-v(i,j+1,k-1)) + (v(i+1,j,k)-v(i+1,j,k-1)) )
C       tm4u = tmu(i,j,k)+tmu(i+1,j+1,k)+tmu(i,j+1,k)+tmu(i+1,j,k)
C       cm2 = epsil + unszz*(uz*uz + vz*vz)
CC      avsdz(i,j,k) = bvf(i,j,k) / cm2
C       avsdz(i,j,k) = bvf(i,j,k) * tm4u * tm4u / cm2
C210  continue
 
c--raccord cyclique et autre (riums) :
      if (ltest.ge.1) then
        do 250 j=jcl1,jcl2
          riums(ims1-1,j) = riums(ims2,j)
          riums(ims2+1,j) = riums(ims1,j)
 250    continue
      endif
      if (ltest.eq.3) then
        riums(ibera, jbera) = riums(iberp, jberp)
      endif
C     call raccord(riums(1,1), 0., 1, 7)
 
c-  Moyenne des 4 Ri(pt.Vitesse) - Weat Point Only :
      ccdz = unsdzw(k) * unsdzw(k)
      do 260 j=js1,js2
       do 260 i=is1(j),is2(j)
        tm4u = tmu(i,j,k-1) + tmu(i+1,j+1,k-1)
     &       + tmu(i,j+1,k-1) + tmu(i+1,j,k-1) + epsil2
        avsdz(i,j,k) = ( (riums(i+1,j+1) + riums(i,j))
     &                 + (riums(i+1,j) + riums(i,j+1)) ) / tm4u
 260  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Coefficients de diffusion Vert. pour les Vitesses .             |
c-----------------------------------------------------------------------
 
c--Coeff. Diffus. Vert. pour les Scalaires : avudz = Viscosite Vert. / dz
      do 430 j=ju1,ju2
       do 430 i=iu1(j),iu2(j)
        deno =  1.0 + alpha*avudz(i,j,k)
        rifac =  1. / (epsil + deno*deno )
        avudz(i,j,k) = avu0ri(k) + avu1ri(k) * min(rifumx, rifac)
 430  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Coefficients de diffusion Vert. pour les Scalaires .            |
c-----------------------------------------------------------------------
 
c--Coeff. Diffus. Vert. pour les Scalaires : avsdz = Diffusivite Vert. / dz
c- (Pac.Phil. simplifie) :
      do 530 j=js1,js2
       do 530 i=is1(j),is2(j)
        unsden = 1.0 / ( abs(1.0 + alpha*avsdz(i,j,k)) + epsil )
        rifac = unsden * unsden * unsden
        avsdz(i,j,k) = avs0ri(k) + avs1ri(k) * min( rifsmx, rifac )
 530  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Fin du traitement separe selon k, Avu,Avs independant/dependant de Ri
      endif
 
      if (k.le.kavsmx) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  6 ) Diffusion Verticale (Deep Ocean) Fonct. de (N2)^q .
c-----------------------------------------------------------------------
 
c--Initialisation :
      do 610 j=1,jmax
       do 610 i=1,imax
        bvfloc(i,j) = bvfmn(1)
 610  continue
 
c- Applique le minorant/Majorant :
      do 620 j=js1,js2
       do 620 i=is1(j)-1,is2(j)+1
        bvfloc(i,j) = max( bvfmn(k), min( bvfmx(k), bvf(i,j,k) ))
 620  continue
 
c- Raccord cyclique & autres <- inutile.
 
c--Calcul de Avs(N2) et incorpore dans avsdz :
      do 650 j=js1,js2
       do 650 i=is1(j),is2(j)
c- Max de N2 sur 9 points :
        bvfavs = max( bvfloc(i,j),
     & bvfloc(i-1,j-1),bvfloc(i+1,j-1),bvfloc(i-1,j+1),bvfloc(i+1,j+1),
     &   bvfloc(i-1,j), bvfloc(i+1,j), bvfloc(i,j-1), bvfloc(i,j+1) )
c- calcul approche de Avs(N2) :
        xav = (bvfavs - bvfmn(1)) * unsdn2
        nav = int(xav)
        xav = xav - DFLOAT(nav)
        avsdz(i,j,k) = avsdz(i,j,k)
     &    + unsdzw(k)*( (1.-xav)*avsn(nav) + xav*avsn(nav+1) )
 650  continue
 
c--Fin du traitement Niveaux profonds (k=<kavsmx) : Avs Fonct de (N2)^q
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Fin de la boucle externe sur l'indice de niveau k .
 700  continue
 
      if (lstab.eq.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  8 ) Ajust.Conv. par augmentation de la Diffusion Verticale .
c-----------------------------------------------------------------------
 
c- Ajoute avsmix(/dzw) si bvf < bvfmix :
      demi = 0.5
      do 830 k=1+ks1,ks2
       do 830 j=js1,js2
        do 830 i=is1(j),is2(j)
         avsdz(i,j,k) = avsdz(i,j,k)
     &         + ( avs2mx(k) + sign(avs2mx(k), bvfmix-bvf(i,j,k)) )
         fqajc(i,j,k) = fqajc(i,j,k)
     &         + ( demi + sign(demi, bvfmix-bvf(i,j,k)) )
 830  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine vdiffu -
      end
