












      subroutine barot
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  prepare l'integration du mode barotrope.
c  remise a zero des variables eta/ub/vb-spl (accumulateur&moyennes).
c  modif : 26/03/99
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
 
c--variables locales equivalentes :
      dimension tm1x2(imax,jmax), tm2x2(imax,jmax)
      equivalence ( tm1x2(1,1)  , phihhh(1,1,1) )
      equivalence ( tm2x2(1,1)  , phihhh(1,1,2) )
      dimension unsdet(imax,jmax), phiypy(imax,jmax)
      dimension phiypx(imax,jmax), phiymx(imax,jmax)
      equivalence ( unsdet(1,1) , phihhh(1,1,3) )
      equivalence ( phiypy(1,1) , phihhh(1,1,4) )
      equivalence ( phiypx(1,1) , phihhh(1,1,5) )
      equivalence ( phiymx(1,1) , phihhh(1,1,6) )
 
c--pour la sortie sur fichier "splout" des variations de eta pendant
c   la derniere iteration :
      include 'reper.com'
C     include 'split1.com'
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Initialisation.                                                 |
c-----------------------------------------------------------------------
 
c--Initialisation :
C     do 100 nn=3,8
C      do 100 j=1,jmax
C       do 100 i=1,imax
C        phihhh(i,j,nn) = 0.0
C100  continue
      do 110 j=1,jmax
       do 110 i=1,imax
        etaspl(i,j) = 0.0
        ubspl(i,j)  = 0.0
        vbspl(i,j)  = 0.0
        tm1x2(i,j) = 0.0
        tm2x2(i,j) = 0.0
        phiypx(i,j) = 0.0
        phiymx(i,j) = 0.0
 110  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Preparation de l'integration du mode barotrope .                |
c-----------------------------------------------------------------------
 
      corfac = 2.*txifcb*dtb

c--Mise en place des masques "interface" :
      j = 1+js2
      do 200 i=ims1,ims2
        tm2x2(i,j) = tms(i,j-1,ks2) * tms(i,j,ks2)
     &             * 2.0 * cmx(i,j,2) * cmy(i,j,2)
 200  continue
 
c--Calcul de la Tension de fond : <- transfere dans "flucor".
 
c--Debut de la boucle externe sur l'indice de latitude j :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,k,afdtb)
      do 300 j=js1,js2
c-----
 
      do 210 i=ims1,ims2+1
        tm1x2(i,j) = tms(i-1,j,ks2) * tms(i,j,ks2)
     &             * 2.0 * cmx(i,j,1) * cmy(i,j,1)
 210  continue
      do 220 i=ims1,ims2
        tm2x2(i,j) = tms(i,j-1,ks2) * tms(i,j,ks2)
     &             * 2.0 * cmx(i,j,2) * cmy(i,j,2)
 220  continue
 
c--Somme les termes de forcage barocline dans fub(ks2) & fvb(ks2) :
c- surface et fond :
      do 230 i=iu1(j),iu2(j)
        fub(i,j,ks2) = fub(i,j,ks2) + phifu(i,j) - phisu(i,j)
        fvb(i,j,ks2) = fvb(i,j,ks2) + phifv(i,j) - phisv(i,j)
C       fub(i,j,ks2) = fub(i,j,ks2) - phisu(i,j)
C       fvb(i,j,ks2) = fvb(i,j,ks2) - phisv(i,j)
 230  continue
      do 240 k=ks1,ks2-1
       do 240 i=iu1(j),iu2(j)
        fub(i,j,ks2) = fub(i,j,ks2) + fub(i,j,k)
        fvb(i,j,ks2) = fvb(i,j,ks2) + fvb(i,j,k)
 240  continue
 
c--determinant, resolution semi Implic. Corriolis :
      do 250 i=iu1(j),iu2(j)
        afdtb = fs2cor(i,j)*corfac
        unsdet(i,j) = tmu(i,j,ku2) / (1.0 + afdtb*afdtb )
 250  continue
 
      if (ahe.ne.zero) then
c--Filtre : calcul des flux diagonaux et N-S (issu de "conti2d"), 1ere ss_iter.
        do 260 i=isf1(j),isf2(j)
          phiypy(i,j) = (eta(i,j-1)-eta(i,j)) * tm2x2(i,j)
 260    continue
        do 270 i=iu1(j),iu2(j)+1
          phiypx(i,j)=(eta(i-1,j-1)-eta(i,j))*tmu(i,j,ks2)*cmxy(i,j,3)
          phiymx(i,j)=(eta(i,j-1)-eta(i-1,j))*tmu(i,j,ks2)*cmxy(i,j,3)
 270    continue
      endif
 
c--Fin de la 1ere boucle externe sur l'indice de latitude j .
 300  continue
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine barot -
      end
