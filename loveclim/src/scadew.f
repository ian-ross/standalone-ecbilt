












      subroutine scadew(scalat,scaldt,alphax,alphay,ns)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  resolution de l'equation d'evolution du scalaire "ns",
c  dans les 2 directions Horizontales, pendant 1 pas de temps.
c--------------------------------------------------
c - Direction X (E-W) d'abord, Y (N-S) ensuite -- |
c--------------------------------------------------
c  modif : 31/05/99
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'isoslope.com'
 
c--variables locales equivalentes :
      dimension phix(imax,jmax) , phiy(imax,jmax)
      dimension phidiv(imax,jmax)
C     equivalence ( phix(1,1)  , phihhh(1,1,1) )
C     equivalence ( phiy(1,1)  , phihhh(1,1,2) )
C     equivalence ( phidiv(1,1) , phihhh(1,1,4) )
 
c--dummy variables :
      dimension scalat(imax,jmax,kmax), scaldt(imax,jmax,kmax)
      dimension alphax(imax,jmax,kmax), alphay(imax,jmax,kmax)
 
c--variables locales :
 
C     if(numit.eq.nstart) then
C      if(ns.eq.1) then
C       write(iuo+66,*) 'scadew : Adv.Dec.Bord(X/Y) = alphah(1/2)'
C       write(iuo+66,*) 'scadew : alphah(1/2) SANS EFFET !'
C      endif
C     endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  0 ) initialisation .
c-----------------------------------------------------------------------
 
c--Debut de la boucle externe sur l'indice de niveau k :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,j,divx,ccx,cc2x,ccxdif,u1cd2,phix)
C$DIR SPP LOOP_PRIVATE( phidiv, ccy,cc2y,ccydif,v2cd2,phiy)
      do 800 k=ks1,ks2
c-----
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Calcul des taux de decentrement dans la direction X  (E-W) .    |
c => transfere dans la routine "alphdec" appelee avant "scadew" ou "scadns"
c--raccord cyclique et autre (alphax,alphay) :
C     call raccord (alphax, 0., 1, 12)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Calcul des Flux advect. & diffus., dans la direction X (E-W) .  |
c-----------------------------------------------------------------------
 
      ccx  = unsdx * dts(k)
      cc2x = ccx + ccx
      ccxdif = ahs(k) * unsdx
c--calcul des flux dans la direction x :
      do 330 j=js1,js2
        do 330 i=is1(j),is2(j)+1
          u1cd2 = 0.25 * cmy(i,j,1) * ( u(i,j,k) + u(i,j+1,k) )
     &           +0.5  * cmy(i,j,1) *   uiso(i,j,k)
C         cnx = smx(i,j,1) * cc2x * u1d2
C         alpha = abs(cnx) + alphmi(k)
C    &          + alphah(1) * (2.0 - tmu(i,j,k) - tmu(i,j+1,k))
C         alpha = min(one, max(alpha, alphax(i,j,k), alphax(i-1,j,k)))
          phidiv(i,j) = cc2x * u1cd2
          phix(i,j) = tms(i,j,k) * tms(i-1,j,k) * (
     &          u1cd2 * (scalat(i,j,k) + scalat(i-1,j,k))
     &        + (alphax(i,j,k) + smxy(i,j,1) * ccxdif )
     &               * (scalat(i-1,j,k) - scalat(i,j,k)) )
 330  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Bilan des Flux advect. & diffus., dans la direction X (E-W) .   |
c-----------------------------------------------------------------------
 
c--Prise en compte du bilan des flux selon X : Mise a jour de scalat() :
      do 450 j=js1,js2
        do 450 i=is1(j),is2(j)
          divx = scal(i,j,k,ns) *
     &              smxy(i,j,0) * (phidiv(i,j) - phidiv(i+1,j))
          scaldt(i,j,k) = scaldt(i,j,k)  + divx
          scalat(i,j,k) = scalat(i,j,k)
     &            + smxy(i,j,0) * ccx * (phix(i,j) - phix(i+1,j))
     &            - divx
 450  continue
 
c--Fin de la 1ere boucle externe sur l'indice de niveau k .
C490  continue
 
c--raccord cyclique et autre (scalat) :
      if (ltest.ge.1) then
        do 470 j=jcl1,jcl2
          scalat(ims1-1,j,k) = scalat(ims2,j,k)
          scalat(ims2+1,j,k) = scalat(ims1,j,k)
 470    continue
      endif
      if (ltest.eq.3) then
        scalat(iberpm,jberp,k) = scalat(ibera, jberam,k)
        scalat(iberp, jberp,k) = scalat(iberam,jberam,k)
        scalat(iberam,jbera,k) = scalat(iberp, jberpm,k)
        scalat(ibera, jbera,k) = scalat(iberpm,jberpm,k)
      endif
 
c--raccord cyclique et autre (scalat) :
C     call raccord(scalat(1,1,k), 0., 1, 0)
 
c--Debut de la 2nd  boucle externe sur l'indice de niveau k :
CC$DIR SPP LOOP_PARALLEL
CC$DIR SPP LOOP_PRIVATE(i,j,ccy,cc2y,ccydif,v2cd2,phiy)
C     do 790 k=ks1,ks2
c-----
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Calcul des taux de decentrement dans la direction Y (N-S) .     |
c => transfere dans la routine "alphdec" appelee avant "scadew" ou "scadns"
c--raccord cyclique et autre (alphax,alphay) :
C     call raccord (alphay, 0., 1, 24)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  6 ) Calcul des Flux advect. & diffus., dans la direction Y (N-S) .  |
c-----------------------------------------------------------------------
 
      ccy  = unsdy * dts(k)
      cc2y = ccy + ccy
      ccydif = ahs(k) * unsdy
c--calcul des flux dans la direction y :
      do 630 j=js1,1+js2
        do 630 i=isf1(j),isf2(j)
          v2cd2 = 0.25 * cmx(i,j,2) * ( v(i,j,k) + v(i+1,j,k) )
     &           +0.5  * cmx(i,j,2) *   viso(i,j,k)
C         cny = smy(i,j,2) * cc2y * v2d2
C         alpha = abs(cny) + alphmi(k)
C    &          + alphah(2) * (2.0 - tmu(i,j,k) - tmu(i+1,j,k))
C         alpha = min(one, max(alpha, alphay(i,j,k), alphay(i,j-1,k)))
          phiy(i,j) = tms(i,j,k) * tms(i,j-1,k) * (
     &          v2cd2 * (scalat(i,j-1,k) + scalat(i,j,k))
     &        + (alphay(i,j,k) + cmxy(i,j,2) * ccydif )
     &               * (scalat(i,j-1,k) - scalat(i,j,k))  )
 630  continue
 
c--fin du calcul des flux.
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  7 ) Bilan des Flux advect. & diffus., dans la direction Y (N-S) .   |
c-----------------------------------------------------------------------
 
c--Prise en compte du bilan des flux selon Y : Mise a jour de scalat() :
      do 750 j=js1,js2
        do 750 i=is1(j),is2(j)
          scalat(i,j,k) = scalat(i,j,k)
     &            + smxy(i,j,0) * ccy * (phiy(i,j) - phiy(i,j+1))
 750  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Fin de la boucle externe sur l'indice de niveau k .
 800  continue
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine scadew -
      end
