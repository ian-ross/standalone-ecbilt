












      subroutine scali
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c**AVANCEMENT IMPLICITE D'UN PAS DE TEMPS SCALAIRE DES SCALAIRES**
c  Traitement implicite des termes de rappel, diffusion et advection .
c    advect & diffus : taux d'implicitete au choix.
c  advection suivant la verticale :  w(k) * [ S(k-1) + S(k) ] / 2
c- resolution separee ou non, selon "ns"
c  modif : 23/03/98
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'isoslope.com'
      include 'comunit.h'
 
c--variables a conserver d'un appel a l'autre :
      logical bcalc
      common / bscali /
     &  bcalc(nsmax)
 
      common / nscali /
     &  nscom(nsmax)
 
      common / xscali /
     &  z2alph(imax,jmax,2:kmax,nsmax)
 
c--variables locales :
      dimension aa(imax,kmax,nsmax), bb(imax,kmax,nsmax)
      dimension cc(imax,kmax,nsmax), ff(imax,kmax)
C     dimension ccwabs(kmax), ccwcor(kmax), ccwdmy(kmax)
C     dimension ccwddn(kmax), ccwdup(kmax)
      dimension cczdt(kmax)
C     dimension ccvois(-kmax:kmax), alpha(imax,kmax)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  0 ) preparation des coeff. intervenant dans les flux d'advection .  |
c-----------------------------------------------------------------------
 
      if(numit.eq.nstart) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--Resolution commune / separee :
      nnsep = nsmax
      do 25 ns=1,nsmax
        bcalc(ns) = .TRUE.
        nscom(ns) = ns
        do 20 nn=1,ns-1
         if (bcalc(ns).and.(alphaz(1-nn+ks2).eq.alphaz(1-ns+ks2))) then
           bcalc(ns) = .FALSE.
           nscom(ns) = nn
           nnsep = nnsep - 1
         endif
 20     continue
 25   continue
 
C     write(iuo+66,'(A,) ' scali : alphaz = run.par + Dec. near a step'
      write(iuo+66,'(A,I2,A,(6F5.2))') ' scali : cal=', nnsep,
     &    ' ; alphaz = Dec.near_step + F(ns):',
     &    (alphaz(1-ns+ks2),ns=1,nsmax)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- initialisation :
      do 40 ns=1,nsmax
       do 40 k=2,kmax
        do 40 j=1,jmax
         do 40 i=1,imax
          z2alph(i,j,k,ns) = 0.0
 40   continue
 
c--Mise en place des coeff (Decentrement Vertical) :
      cc0w = 0.5 * alphaz(ks1)
      do 55 ns=1,nsmax
       if (bcalc(ns)) then
         cc1w = 0.5 * txiads * alphaz(1-ns+ks2)
         cc2w = 0.5 * (1.0 - alphaz(1-ns+ks2))
         do 50 k=ks1+1,ks2
C         cc1w = 0.5 * txiads * alphaz(k)
C         cc2w = 0.5 * (1.0 - alphaz(k))
          do 50 j=js1,js2
           do 50 i=is1(j),is2(j)
             ccstep = ( tms(i-1,j,k) - tms(i-1,j,k-1) )
     &              + ( tms(i+1,j,k) - tms(i+1,j,k-1) )
     &              + ( tms(i,j-1,k) - tms(i,j-1,k-1) )
     &              + ( tms(i,j+1,k) - tms(i,j+1,k-1) )
             z2alph(i,j,k,ns) = cc1w + min(cc2w , cc0w * ccstep)
 50      continue
       endif
 55   continue
c--Fin du traitement specifique de la 1ere itt.
      endif
 
c convention :
c    CCW* coeff relatif a l'interface (= place des W).
c    CCZ* coeff relatif au centre des boites (= place de T,S)
c    CCWI : tient compte du taux Implicite.
c    CCWABS(k) intervient avec |w| ;
c              CCWUP(k) avec scal(k) ; CCWDowN(k) avec scal(k-1)
c  pour le taux de Decentrement calcule a partir du Nb de courant : CCWDUP(k)
c    avec dts(k), CCWDDowN(k) avec dts(k-1) ey CCWDMY(k) avec la moyenne des 2
c    CCZDT**(k) intervient pour le bilan des flux V, boite "k" ;
 
      ccwi   = 0.5 * txiads
      do 80 k=ks1,ks2
        cczdt(k) = dts(k) * unsdz(k)
C       ccwabs(k) = ccwi * alphaz(k)
 80   continue
C     ccwid2 = 0.25 * txiads
C     do 90 k=ks1+1,ks2
C       ccwcor(k) = ccwid2 * unsdzw(k) * (dz(k) - dz(k-1))
C    &                     * (1. - alphaz(k))
C       ccwdup(k) = ccwi * unsdzw(k) * dts(k)
C       ccwddn(k) = ccwi * unsdzw(k) * dts(k-1)
C       ccwdmy(k) = ccwi * unsdzw(k) * (dts(k) + dts(k-1)) * 0.5
C       cczdif(k) = txidfs * unsdzw(k)
 90   continue
 
c--Debut de la boucle externe sur l'indice de latitude j :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,k,ns,nn,ccdif,ccadv,ccdt,aa,bb,cc,ff)
      do 600 j=js1,js2
 
c-----
c--debut de la boucle sur "ns".
      do 500 ns=1,nsmax
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) construction de la matrice (commune a certains scalaires).      |
c-----------------------------------------------------------------------
 
c--Ininialisation de la diagonale (aa) et Mise en place du 2nd membre (ff) :
c- transfert de scal(*,j,*,ns) dans ff(*,*) & Traite le terme rappel implic.
      do 100 k=ks1,ks2
       do 100 i=is1(j),is2(j)
        aa(i,k,ns) = 1.0 + rappel(i,j,k)
        ff(i,k) = scal(i,j,k,ns) + rappel(i,j,k) * scalr(i,j,k,ns)
 100  continue
 
      if (bcalc(ns)) then
c--Traitement commun a certains scalaires :
 
c-N.B.-: avsdz = txidfs * Dif.Vert. / dzw !!
      do 110 k=ks1+1,ks2
       do 110 i=is1(j),is2(j)
        ccdif = avsdz(i,j,k)
     &        + ai(k)*( c4x(i,j,k)*c4x(i,j,k)
     &                 +c4y(i,j,k)*c4y(i,j,k) )*unsdzw(k)
     &        + z2alph(i,j,k,ns) * abs(w(i,j,k)+wiso(i,j,k))
Cis0 &        + z2alph(i,j,k,ns) * abs(w(i,j,k))
C    &        + ccwabs(k) * abs(w(i,j,k))
C    &        + ccwcor(k) * w(i,j,k)
C    &        + ccwdmy(k) * w(i,j,k) * w(i,j,k)
C    &        + ccwddn(k) * w(i,j,k) * w(i,j,k)
C    &        + ccwdup(k) * w(i,j,k) * w(i,j,k)
        ccadv = ccwi * w(i,j,k)
     &         +ccwi * wiso(i,j,k)
 
c- effet du flux phiz(k) sur S(k) : aa(k)*S(k) + bb(k)*S(k-1)
        ccdt = cczdt(k) * tms(i,j,k-1)
        bb(i,k,ns) = ccdt * (-ccdif - ccadv)
        aa(i,k,ns) = aa(i,k,ns) + ccdt * (ccdif - ccadv)
c- effet du flux phiz(k) sur S(k-1) : aa(k-1)*S(k-1) + cc(k-1)*S(k)
        ccdt = cczdt(k-1) * tms(i,j,k-1)
        aa(i,k-1,ns) = aa(i,k-1,ns) + ccdt * (ccdif + ccadv)
        cc(i,k-1,ns) = ccdt * (-ccdif + ccadv)
 110  continue
c**fin construction de la matrice du systeme commun**
 
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) decomposition LU - (methode et notation : cf Linear Algebra pp165-167)
c-----------------------------------------------------------------------
 
c--calcul de 1/alpha(k) dans aa(i,k) et beta(k) dans bb(i,k)
      do 200 i=is1(j),is2(j)
        aa(i,ks1,ns) = 1.0 / aa(i,ks1,ns)
 200  continue
      do 210 k=ks1+1,ks2
       do 210 i=is1(j),is2(j)
        bb(i,k,ns) = bb(i,k,ns) * aa(i,k-1,ns)
        aa(i,k,ns) = 1.0 / ( aa(i,k,ns) - bb(i,k,ns) * cc(i,k-1,ns) )
 210  continue
c**fin decomposition LU commune**
 
c--Fin de la partie commune a certains scalaires / debut du traitement separe
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) substitutions avant et arriere pour chaque scalaire :           |
c-----------------------------------------------------------------------
 
      nn = nscom(ns)
c--calcul de g(k) dans ff(i,k) :
      do 410 k=ks1+1,ks2
       do 410 i=is1(j),is2(j)
        ff(i,k) = ff(i,k) - bb(i,k,nn) * ff(i,k-1)
 410  continue
 
c--calcul de x(k) dans scal(i,j,k,ns) :
      do 420 i=is1(j),is2(j)
       scal(i,j,ks2,ns) = ff(i,ks2) * aa(i,ks2,nn)
 420  continue
      do 430 k=ks2-1,ks1,-1
       do 430 i=is1(j),is2(j)
        scal(i,j,k,ns) = (ff(i,k) - cc(i,k,nn)*scal(i,j,k+1,ns))
     &                 * aa(i,k,nn)
 430  continue
 
c--fin substitutions avant et arriere pour chaque scalaire.
c-----
 
c--Fin de la boucle sur l'indice du scalaire ns .
 500  continue
 
c--Fin de la boucle externe sur l'indice de latitude j .
 600  continue
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine scali -
      end
