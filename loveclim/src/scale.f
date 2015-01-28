












      subroutine scale(nsew,nn99)


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c**AVANCEMENT EXPLICITE D'UN PAS DE TEMPS SCALAIRE DES SCALAIRES**
c traitement explicite a "pas de temps fractionnaires", pour les 3 D.
c  et avec taux d'implicitete (a choisir) selon la verticale.
c ordre : 1er=Vert.expli. Si nsew=1 : 2eme=Est-West, 3eme=Nord-Sud, #0 inverse
c--
c Tableaux : scal <- Valeur conserve tel quel jusque fin de "scale" .
c          scalat <- Valeur actualise a chaque "fragment" de  pas de temps.
c          scaldt <- Variation (+correction) totalement explicite.
c--
c En entree : avsdz = Coeff.Diffus.Vert. entier ; en sortie = implicite only
c  advection suivant la verticale :  w(k) * [ S(k-1) + S(k) ] / 2
c  modif : 06/10/99
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'ice.com'
      include 'comunit.h'
 
c--variables locales equivalentes :
      dimension phix(imax,jmax), phiy(imax,jmax)
C     equivalence ( phix(1,1)  , phihhh(1,1,1) )
C     equivalence ( phiy(1,1)  , phihhh(1,1,2) )
 
      dimension alphhk(ixjmax,kmax,2)
C     dimension alphax(imax,jmax,kmax), alphay(imax,jmax,kmax)
C     equivalence ( alphhk(1,1,1) , alphax(1,1,1) )
C     equivalence ( alphhk(1,1,2) , alphay(1,1,1) )
 
c--variables locales :
      dimension scalat(imax,jmax,kmax), scaldt(imax,jmax,kmax)
      dimension difzex(imax,kmax), phiz(imax,kmax+1)
      dimension ccwup(kmax), ccwdn(kmax), ccwabs(kmax)
      dimension cczdt(kmax), cczdte(kmax), ccwdmy(kmax)
      dimension difzmx(kmax)
c-
      dimension kslpdw(nlpmax)
      dimension cdtsxz(kmax), cslpmx(kmax)
 
      common / downsloping / kslpdw
c--Decroissance radioactive du C14 (yr^-1)  [<-> 1/2 vie = 5730.yr ]
      data decay / 1.21d-04 /
 
      if (numit.eq.nstart) then
        write(iuo+66,*) 'scale : Adv. Zex,(X,Y) Alterne, + Corrigee.'
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) preparation des coeff. intervenant dans les flux d'advection .  |
c-----------------------------------------------------------------------
 
c convention :
c    CCW* coeff relatif a l'interface (= place des W).
c    CCZ* coeff relatif au centre des boites (= place de T,S)
c    CC*E : tient compte du taux Explicite (= TXE*) .
c    CCWABS(k) intervient avec |w| ;
c              CCWUP(k) avec scal(k) ; CCWDowN(k) avec scal(k-1).
c  pour le taux de Decentrement calcule a partir du Nb de courant : CCWDUP(k)
c    avec dts(k), CCWDDowN(k) avec dts(k-1), CCWDMY(k) avec la moyenne des 2.
c    CCZDT**(k) intervient pour le bilan des flux V, boite "k" ;
 
      unsyr = one / (yeaday * 86400.)
c--initialisation (indispensable) :
C     do 110 n=1,6
C      do 110 j=1,jmax
C       do 110 i=1,imax
C        phihhh(i,j,n) = 0.0
C110  continue
      do 120 k=1,kmax
       do 120 i=1,imax
        difzex(i,k) = 0.0
 120  continue
      do 125 k=1,kmax+1
       do 125 i=1,imax
        phiz(i,k) = 0.0
 125  continue
 
c--preparation des coeffs dependant de la verticale :
 
      txeadv = 1.0 - txiads
      txedif = 1.0 - txidfs
      ccwe = 0.5 * txeadv
      do 150 k=ks1+1,ks2
C       ccwup(k)  = unsdz4(k) * 0.25 * (1.0 - txiads) *
C    a      ((2.0-alphaz(k))*dz(k-1)+alphaz(k)*dz(k))
C       ccwdn(k)  = unsdz4(k) * 0.25 * (1.0 - txiads) *
C    a      ((2.0-alphaz(k))*dz(k)+alphaz(k)*dz(k-1))
        ccwabs(k) = ccwe * alphaz(k)
C       ccwddn(k) = ccwe * unsdzw(k) * dts(k-1)
        ccwdmy(k) = ccwe * unsdzw(k) * (dts(k-1) + dts(k)) * 0.5
 150  continue
 
c--calcul de la diffusivite/dz explicite maximale (restant stable) :
      do 170 k=ks1,ks2
        difzmx(k) = 0.5 * dz(k) / dts(k)
        cczdt(k)  = unsdz(k) * dts(k)
        cczdte(k) = cczdt(k) * txeadv
        cdtsxz(k) = unsdx * cczdt(k)
        cslpmx(k) = 0.125 * dx / cczdt(k)
 170  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Detecte les courant de Downsloping .
c-----------------------------------------------------------------------
 
      call slopez(cslpmx, kslpdw, nn99)
 
c--Debut de la boucle externe sur tous les scalaires, indices "ns" :
c- N.B : Parallelisable SAUF si txiads ou txidfs different de 1.
CoDIR SPP LOOP_PARALLEL
CoDIR SPP LOOP_PRIVATE(i,j,k,kfd)
CoDIR SPP LOOP_PRIVATE(scaldt,scalat,sscal,ccflxt,fflx)
CoDIR SPP LOOP_PRIVATE(phix,phiy,phiz,phdivz,alphhk)
CoDIR SPP LOOP_PRIVATE(xdecay,ttyear)
      do 800 ns=1,nsmax
c-----
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Conditions Limites : Surface & Fond - Modification de "avsdz" . |
c-----------------------------------------------------------------------
 
      ccflxt = 0.
      if (ns.eq.1) ccflxt = dts(ks2)*unsdz(ks2) / (rho0*cpo)
c- Debut de la 1ere boucle externe sur la latitude, indice "j" :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,k)
C$DIR SPP LOOP_PRIVATE(sscal,fflx,kfd)
      do 360 j=1,jmax
 
c--Initialisation de scalat (<- scal) et scaldt (par le terme de "derive") :
       do 310 k=1,kmax
        sscal = -deriv(ns) * dts(k)
        do 310 i=1,imax
         scalat(i,j,k) = scal(i,j,k,ns)
Cic0     scaldt(i,j,k) = sscal
         scaldt(i,j,k) = sscal - phivs(i,j,k,ns)
 310  continue
 
c--Calcul du Forcage(explicite) en Surface : (phiss <-deja.mult par deltaT/Dz)
      k = ks2
C     do 320 j=js1,js2
       do 320 i=is1(j),is2(j)
         fflx = phiss(i,j,ns)
     &     + rappes(i,j,ns) * ( scal(i,j,k,ns) - scalr(i,j,k,ns) )
Cic0     fflx = max( phimnx(i,j,0,ns), min(phimnx(i,j,1,ns),fflx) )
         fss(i,j,ns) = fss(i,j,ns) + fflx
Cic0 &      + cczdt(ks2) * w(i,j,ks2+1) * (scs(i,j,ns) - scssv(ns))
     &               - ccflxt * fcm1(i,j)
         scaldt(i,j,k) = scaldt(i,j,k) - fflx
Cic0 &      - cczdt(ks2) * w(i,j,ks2+1) * scs(i,j,ns)
C1       fflx0 = txeflx * fflx
C1       scalat(i,j,k) = scalat(i,j,k) - fflx0
C1       scaldt(i,j,k) = scaldt(i,j,k) + fflx0 - fflx
 320  continue
 
c--Expansion par le fond : (<- Not yet multipl. by deltaT)
C     do 350 j=js1,js2
       do 350 i=is1(j),is2(j)
         kfd = kfs(i,j)
C        w(i,j,kfd) = 0.0
         scaldt(i,j,kfd) = scaldt(i,j,kfd)
     &                   + cczdt(kfd) * phifs(i,j,ns)
 350  continue
 
c--fin de la 1ere boucle sur l'indice de latitude "j".
 360  continue
 
Cadh  if (txeflx(ns).gt.epsil) then
c- calcul du Flux derive de l'Advection H. des Obs.
Cadh  do 364 j=js1,js2
Cadh    do 364 i=is1(j),is2(j)+1
Cadh     phix(i,j) = flxus(i,j,ns) * (u(i,j,k) + u(i,j+1,k))
 364   continue
Cadh   do 366 j=js1,js2+1
Cadh    do 366 i=isf1(j),isf2(j)
Cadh     phiy(i,j) = flxvs(i,j,ns) * (v(i,j,k) + v(i+1,j,k))
 366   continue
Cadh   do 368 j=js1,js2
Cadh    do 368 i=is1(j),is2(j)
Cadh     scaldt(i,j,k) = scaldt(i,j,k) - smxy(i,j,0) *
Cadh &    ( (phix(i,j) + phix(i+1,j)) + (phiy(i,j) + phiy(i,j+1)) )
 368   continue
Cadh  endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Forcage thermique : Diffus.Anom.Temp. (<-S.Rahmstorf) :
      if (nitrap.ne.0 .and. ns.eq.1) call rahmflx(ns,nn99,scaldt)
 
      if (nsmax.gt.2) then
c--Terme de viellissement C14 & Age (ns=3,4) :
      if (ns.eq.3) then
        do 380 k=ks1,ks2
         xdecay = decay * dts(k) * unsyr
         do 380 j=js1,js2
          do 380 i=is1(j),is2(j)
            scaldt(i,j,k) = scaldt(i,j,k) - xdecay * scal(i,j,k,ns)
 380    continue
      elseif (ns.eq.4) then
        do 390 k=ks1,ks2
         ttyear = dts(k) * unsyr
         do 390 j=js1,js2
          do 390 i=is1(j),is2(j)
            scaldt(i,j,k) = scaldt(i,j,k) + ttyear
 390    continue
      endif
c-----
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Direction Verticale : Traitement Explicite .                    |
c-----------------------------------------------------------------------
 
      if ( (txiads.ne.1.0).or.(txidfs.ne.1.0) ) then
c-----
 
      do 470 j=js1,js2
c--debut de la 2nd boucle sur l'indice de latitude "j".
 
      if (txidfs.ne.1.0) then
c--4.1 calcul du coefficient de diffusion verticale explicite "difzex" :
      do 410 k=ks1+1,ks2
       do 410 i=is1(j),is2(j)
         difzex(i,k) = txedif * min( avsdz(i,j,k) , difzmx(k) )
 410  continue
      if (ns.eq.nsmax) then
c- transformation du tableaux de Diffusivite Verticale "avsdz" .
        do 420 k=ks1+1,ks2
         do 420 i=is1(j),is2(j)
           avsdz(i,j,k) = avsdz(i,j,k) - difzex(i,k)
 420    continue
      endif
c--
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--4.3 Flux Verticaux advect. & diffus. - Integre les Flux Explicites. |
 
c--calcul des flux verticaux explicite advection et diffusion.
      do 430 k=ks1+1,ks2
       do 430 i=is1(j),is2(j)
         phiz(i,k) = tms(i,j,k-1) * (
C    &      ( difzex(i,k) + ccwdmy(k)*w(i,j,k)*w(i,j,k) )
     &      ( difzex(i,k) + ccwabs(k) * abs(w(i,j,k)) )
     &       * ( scalat(i,j,k-1) - scalat(i,j,k) )
     &    + ccwe * w(i,j,k) * (scalat(i,j,k-1) + scalat(i,j,k))    )
C    a    + (ccwdn(k) * scal(i,j,k-1,ns) + ccwup(k) * scal(i,j,k,ns) )
C    a       * w(i,j,k)                           )
 430  continue
c--fin du calcul des flux verticaux.
 
c--Bilan des flux Verticaux explicites :
      do 450 k=ks1,ks2
       do 450 i=is1(j),is2(j)
        phdivz = cczdte(k) * (w(i,j,k) - w(i,j,k+1)) * scal(i,j,k,ns)
        scalat(i,j,k) = scalat(i,j,k) - phdivz
     &                 + cczdt(k) * (phiz(i,k) - phiz(i,k+1))
        scaldt(i,j,k)  = scaldt(i,j,k)  + phdivz
C       scaldt(i,j,k)  = scaldt(i,j,k)
C    a                 + cczdt(k) * (phiz(i,k) - phiz(i,k+1))
 450  continue
 
c--fin de la 2nd boucle sur l'indice de latitude "j".
 470  continue
 
c--raccord cyclique et autre (scalat) :
C     call raccord(scalat(1,1,1),zero,kmax,0)
 
c--fin distinction tout implicite / en partie explicite .
c----------------------------
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Calcul des Coef. pour le decentrement Horizontal .              |
c-----------------------------------------------------------------------
 
      if (alphgr(ns).lt.-1.0) then
        call alph2dc(alphhk(1,1,1),alphhk(1,1,2),ns)
      else
        call alphdec(alphhk(1,1,1),alphhk(1,1,2),ns)
      endif
 
c- Dowsloping : Modif par Permutation :
 
      call slopes(scalat,alphhk,cdtsxz,kslpdw,ns)
 
c- Isopycn.Diffus.: Compute (fr scal) & Incorporate (to scalat) Diagon.Flx.
      if (ai(kmax).gt.zero) call isodiffu(scalat,ns)
 
c--raccord cyclique et autre (scalat) :
      call raccord(scalat(1,1,1),zero,kmax,0)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  6 ) Flux Horizontaux advect. & diffus - Integre les Flux Explicites |
c-----------------------------------------------------------------------
 
c--bifurcation sur l'ordre de traitement des 2 directions N-S , E-W :
      if(nsew.ne.1) then
 
c--Resolution des Flux advect. & diffus., E-W en 1er, N-S en 2nd :
        call scadew(scalat,scaldt,alphhk(1,1,1),alphhk(1,1,2),ns)
 
      else
 
c--Resolution des Flux advect. & diffus., N-S en 1er, E-W en 2nd :
        call scadns(scalat,scaldt,alphhk(1,1,1),alphhk(1,1,2),ns)
 
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  7 ) Incorpore la partie totalement explicite (stokee dans scaldt).  |
c-----------------------------------------------------------------------
 
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(n,i,j)
      do 760 k=ks1,ks2
 
c--Ajoute le terme de rappel explicite interne (= ailleurs qu'en surf.) :
       do 740 n=1,nrap(k,ns)
         i = 1 + mod(ijrap(n,k,ns),imax)
         j = 1 + ijrap(n,k,ns)/imax
         scaldt(i,j,k) = scaldt(i,j,k)
     &     + rapint(n,k,ns) * ( scalr(i,j,k,ns) - scal(i,j,k,ns) )
 740   continue
c--
       do 750 j=js1,js2
        do 750 i=is1(j),is2(j)
         scal(i,j,k,ns) = scalat(i,j,k) + scaldt(i,j,k) * tms(i,j,k)
 750   continue
c----
 760  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Fin de la boucle externe sur tous les scalaires, indicies "ns" .
 800  continue
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine scale -
      end
