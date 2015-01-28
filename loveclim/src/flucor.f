












      subroutine flucor

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c routine de calcul des corrections de flux (pour les scalaires)
c  qui interviennent au fond de l'ocean.
c ( ATTENTION : Pas encore multiplie par le pas de temps ! )
c  modif : 03/07/98

      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'

c--variables locales :
      real dum,wdum(imax,jmax)


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Transfert dans PhiFond de Wfond , puis mise a  zero .           |
c-----------------------------------------------------------------------

c--Stockage dans phifs(-,-,1) de W(kfd) ; Mise a zero de W(kfd) ;
c- et pour la routine informe : Somme de |W(kfd)| stockee dans daeta .
      do 110 j=js1,js2
       do 110 i=is1(j),is2(j)
        kfd = kfs(i,j)
        phifs(i,j,1) = w(i,j,kfd)
        daeta(i,j) = daeta(i,j) + abs( w(i,j,kfd) )
        w(i,j,kfd) = 0.0
        fss(i,j,0) = fss(i,j,0) + w(i,j,ks2+1)
 110  continue
c
C-AM
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1b ) Modification of sea surface w for mass change .           |
c-----------------------------------------------------------------------
c-- stockage dans wdum de w(ks2+1) modifiee
      dum=zfluxm*vcor
c     dum=zflux0*vcor <- moyenne instantanee
      do j=js1,js2
        do i=is1(j),is2(j)
          wdum(i,j)=w(i,j,ks2+1)-dum*tms(i,j,ks2)
       enddo
      enddo
C-AM

      do 400 ns=nsmax,1,-1
c--boucle sur les scalaires (ns) :


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Calcul des Flux au fond .
c-----------------------------------------------------------------------

c--calcul de la correction sur le flux entrant par le fond :
      do 210 j=js1,js2
       do 210 i=is1(j),is2(j)
         phifs(i,j,ns) = phifs(i,j,1) * scal(i,j,kfs(i,j),ns)
 210  continue

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Calcul de la correction globale (=la derive) .                  |
c-----------------------------------------------------------------------

c--calcul du gain total entrant (=la derive), suite a la correction :
      deriv(ns) = 0.0
C-AM
c--derive associee au Flux de masse en surf. (= Pluie / Evap.) 
c--   -> depend du type de scalaire

      if (scpme(ns).eq.spvr) then
       do 320 j=js1,js2
         do 320 i=is1(j),is2(j)
           deriv(ns) = deriv(ns) + cmxy(i,j,0) *
     &          ( phifs(i,j,ns) - w(i,j,ks2+1)*scal(i,j,ks2,ns) )
 320   continue
      else
       do 330 j=js1,js2
         do 330 i=is1(j),is2(j)
           deriv(ns) = deriv(ns) + cmxy(i,j,0) *
     &          ( phifs(i,j,ns) - wdum(i,j)*scssv(ns) )
 330   continue
       endif
C-AM
       deriv(ns) = deriv(ns) * unsvol

c--fin de la boucle sur les scalaire.
 400  continue

      return

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine flucor -
      end
