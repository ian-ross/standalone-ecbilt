












      subroutine alphdec(alphax,alphay,ns)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  calcule les taux de decentrement alphax & alphay dans les 2 directions .
c   en sortie : multiplie par |u1.cmy/2|, |v2.cmx/2|
c  modif : 25/05/99
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'isoslope.com'
 
c--variables locales equivalentes :
C     equivalence ( phiax(1,1)  , phihhh(1,1,5) )
C     equivalence ( phiay(1,1)  , phihhh(1,1,6) )
C     equivalence ( alphxx(1,1)  , phihhh(1,1,7) )
C     equivalence ( alphyy(1,1)  , phihhh(1,1,8) )
 
c--dummy variables :
      dimension alphax(imax,jmax,kmax), alphay(imax,jmax,kmax)
 
c--variables locales :
      dimension phiax(imax,jmax) , phiay(imax,jmax)
      dimension alphxx(imax,jmax) , alphyy(imax,jmax)
 
      if (numit.eq.nstart) then
       if (alphgr(ns).gt.zero) then
Ca2     write(iuo+66,'(A,I3,F10.6,1PE14.6)') ' alphdec : sans alphah, A*A ;'
        write(iuo+66,'(A,I3,F10.6,1PE14.6)') ' alphdec : sans alphah, |A| ;'
     &            //' ns,alphgr,algrmn =', ns, alphgr(ns), algrmn(ns)
       elseif (alphgr(ns).eq.zero .and. ns.eq.1) then
        write(iuo+66,'(A,I3,F10.6,1PE14.6)') ' alphdec : sans alphah ;'
C      else
C       write(iuo+66,'(A,I3,F10.6,1PE14.6)') ' alphdec : sans alphah ;'
C    &            //' ns,alphgr,algrmn =', ns, alphgr(ns), algrmn(ns)
       endif
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  0 ) preparation des coeff. intervenant dans les flux d'advection .  |
c-----------------------------------------------------------------------
 
      if ( alphgr(ns).eq.zero ) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1a) Decentrement Minimum pour stabilite = Nb.de courant :
c-----------------------------------------------------------------------
 
c--Debut de la boucle externe sur l'indice de niveau k :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,j,cc2x,u1pd2,alpha,cc2y,v2pd2)
      do 150 k=1,kmax
c-----
 
      cc2x  = 2. * unsdx * dts(k)
c--calcul du taux de decentrement dans la direction x :
      do 110 j=js1,js2
        do 110 i=is1(j),is2(j)+1
Cis0      u1pd2 = abs( 0.25 * (u(i,j,k) + u(i,j+1,k)) )
          u1pd2 = abs( 0.25 * (u(i,j,k) + u(i,j+1,k))
     &                + 0.5 *  uiso(i,j,k)            )
          alpha = smx(i,j,1) * cc2x * u1pd2 + alphmi(k)
C         cnx = smx(i,j,1) * cc2x * u1d2
C         alpha = abs(cnx) + alphmi(k)
C    &          + alphah(1) * (2.0 - tmu(i,j,k) - tmu(i,j+1,k))
          alphax(i,j,k) = cmy(i,j,1) * u1pd2 * min(one, alpha)
 110  continue
 
c--calcul du taux de decentrement dans la direction y :
      cc2y  = 2. * unsdy * dts(k)
      do 130 j=js1,1+js2
        do 130 i=isf1(j),isf2(j)
Cis0      v2pd2 = abs( 0.25 * (v(i,j,k) + v(i+1,j,k)) )
          v2pd2 = abs( 0.25 * (v(i,j,k) + v(i+1,j,k))
     &                + 0.5 *  viso(i,j,k)            )
          alpha = smy(i,j,2) * cc2y * v2pd2 + alphmi(k)
C         cny = smy(i,j,2) * cc2y * v2d2
C         alpha = abs(cny) + alphmi(k)
C    &          + alphah(2) * (2.0 - tmu(i,j,k) - tmu(i+1,j,k))
          alphay(i,j,k) = cmx(i,j,2) * v2pd2 * min(one, alpha)
 130  continue
 
c--Fin de la boucle externe sur l'indice de niveau k .
 150  continue
 
      return
 
      elseif ( alphgr(ns).eq.2.0d0 ) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1b) Schema Upwind :
c-----------------------------------------------------------------------
 
c--Debut de la boucle externe sur l'indice de niveau k :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,j,cc2x,u1pd2,alpha,cc2y,v2pd2)
      do 180 k=1,kmax
 
c--calcul du taux de decentrement dans la direction x :
      do 160 j=js1,js2
        do 160 i=is1(j),is2(j)+1
Cis0      u1pd2 = abs( 0.25 * (u(i,j,k) + u(i,j+1,k)) )
          u1pd2 = abs( 0.25 * (u(i,j,k) + u(i,j+1,k))
     &                + 0.5 *  uiso(i,j,k)            )
          alphax(i,j,k) = cmy(i,j,1) * u1pd2
 160  continue
 
c--calcul du taux de decentrement dans la direction y :
      do 170 j=js1,1+js2
        do 170 i=isf1(j),isf2(j)
Cis0      v2pd2 = abs( 0.25 * (v(i,j,k) + v(i+1,j,k)) )
          v2pd2 = abs( 0.25 * (v(i,j,k) + v(i+1,j,k))
     &                + 0.5 *  viso(i,j,k)            )
          alphay(i,j,k) = cmx(i,j,2) * v2pd2
 170  continue
 
c--Fin de la boucle externe sur l'indice de niveau k .
 180  continue
c-----
      return
 
      elseif ( alphgr(ns).gt.zero) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Calcul des Coeff. "alpha" : formule "TVM" .                     |
c-----------------------------------------------------------------------
 
      grdmin = max( epsil, algrmn(ns) )
 
c--Debut de la boucle externe sur l'indice de niveau k :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,j,phiax,phiay,alpha1,alpha2,alphxx,alphyy)
C$DIR SPP LOOP_PRIVATE(cc2x,u1pd2,alpha,cc2y,v2pd2)
      do 380 k=ks1,ks2
c----
 
c--calcul des gradients dans les 2 directions :
      do 210 j=js1,js2
       do 210 i=is1(j),1+is2(j)
        phiax(i,j) = abs(scal(i,j,k,ns) - scal(i-1,j,k,ns))
C    &             * tms(i,j,k) * tms(i-1,j,k)
 210  continue
      do 220 j=js1,1+js2
       do 220 i=isf1(j),isf2(j)
        phiay(i,j) = abs(scal(i,j,k,ns) - scal(i,j-1,k,ns))
C    &             * tms(i,j,k) * tms(i,j-1,k)
 220  continue
 
c--calcul des coeff alphax & alphyy :
      do 230 j=js1,js2
       do 230 i=is1(j),is2(j)
        alpha1 = ( phiax(i+1,j) - phiax(i,j) )
     &         / ( phiax(i+1,j) + phiax(i,j) + grdmin )
        alpha2 = ( phiay(i,j+1) - phiay(i,j) )
     &         / ( phiay(i,j+1) + phiay(i,j) + grdmin )
Ca2     alphxx(i,j) = alphgr(ns) * alpha1 * alpha1
        alphxx(i,j) = alphgr(ns) * abs(alpha1)
     &              * tms(i+1,j,k) * tms(i,j,k) * tms(i-1,j,k)
Ca2     alphyy(i,j) = alphgr(ns) * alpha2 * alpha2
        alphyy(i,j) = alphgr(ns) * abs(alpha2)
     &              * tms(i,j+1,k) * tms(i,j,k) * tms(i,j-1,k)
 230  continue
 
C     do 210 j=js1,js2
C      do 210 i=is1(j),is2(j)
C       dp = abs( scal(i+1,j,k,ns) - scal(i,j,k,ns) )
C       dm = abs( scal(i,j,k,ns) - scal(i-1,j,k,ns) )
C       alpha = ( dp - dm ) / ( dp + dm + epsil )
C       alpha = alphgr(ns) * alpha * alpha
CC      alpha = alphgr(ns) * abs( alpha )
C       alphxx(i,j) = alpha * tms(i+1,j,k)*tms(i,j,k)*tms(i-1,j,k)
C210  continue
C     do 230 j=js1,js2
C      do 230 i=is1(j),is2(j)
C       dp = abs( scal(i,j+1,k,ns) - scal(i,j,k,ns) )
C       dm = abs( scal(i,j,k,ns) - scal(i,j-1,k,ns) )
C       alpha = (dp - dm ) / ( dp + dm + epsil )
C       alpha = alphgr(ns) * alpha * alpha
CC      alpha = alphgr(ns) * abs( alpha )
C       alphyy(i,j) = alpha * tms(i,j+1,k)*tms(i,j,k)*tms(i,j-1,k)
C230  continue
 
c--raccord cyclique et autre (alphxx, alphyy) :
      if (ltest.ge.1) then
        do 270 j=jcl1,jcl2
          alphxx(ims1-1,j) = alphxx(ims2,j)
          alphxx(ims2+1,j) = alphxx(ims1,j)
 270    continue
      endif
      if (ltest.eq.3) then
        alphyy(iberpm,jberp) = alphyy(ibera, jberam)
        alphyy(iberp, jberp) = alphyy(iberam,jberam)
        alphyy(iberam,jbera) = alphyy(iberp, jberpm)
        alphyy(ibera, jbera) = alphyy(iberpm,jberpm)
      endif
 
c--raccord cyclique et autre (alphxx, alphyy) :
C     call raccord( alphxx, 0., kmax, 12)
C     call raccord( alphyy, 0., kmax, 24)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Calcul des Taux de decentrement a partir des Coeff "alpha" :    |
c-----------------------------------------------------------------------
 
      cc2x = 2. * unsdx * dts(k)
c--calcul du taux de decentrement dans la direction x :
      do 310 j=js1,js2
        do 310 i=is1(j),is2(j)+1
Cis0      u1pd2 = abs( 0.25 * (u(i,j,k) + u(i,j+1,k)) )
          u1pd2 = abs( 0.25 * (u(i,j,k) + u(i,j+1,k))
     &                +0.5  *  uiso(i,j,k)            )
          alpha = smx(i,j,1) * cc2x * u1pd2 + alphmi(k)
C         cnx = smx(i,j,1) * cc2x * u1d2
C         alpha = abs(cnx) + alphmi(k)
C    &          + alphah(1) * (2.0 - tmu(i,j,k) - tmu(i,j+1,k))
          alphax(i,j,k) = cmy(i,j,1) * u1pd2 *
     &       min(one, max(alpha, alphxx(i,j), alphxx(i-1,j)))
 310  continue
 
c--calcul du taux de decentrement dans la direction y :
      cc2y  = 2. * unsdy * dts(k)
      do 330 j=js1,1+js2
        do 330 i=isf1(j),isf2(j)
Cis0      v2pd2 = abs( 0.25 * (v(i,j,k) + v(i+1,j,k)) )
          v2pd2 = abs( 0.25 * (v(i,j,k) + v(i+1,j,k))
     &                + 0.5 *  viso(i,j,k)            )
          alpha = smy(i,j,2) * cc2y * v2pd2 + alphmi(k)
C         cny = smy(i,j,2) * cc2y * v2d2
C         alpha = abs(cny) + alphmi(k)
C    &          + alphah(2) * (2.0 - tmu(i,j,k) - tmu(i+1,j,k))
          alphay(i,j,k) = cmx(i,j,2) * v2pd2 *
     &       min(one, max(alpha, alphyy(i,j), alphyy(i,j-1)))
 330  continue
 
c--Fin de la boucle externe sur l'indice de niveau k .
 380  continue
 
      else
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Calcul des Coeff. "alpha" visant a supprimer le Depassement .   |
c-----------------------------------------------------------------------
 
      grdmin = max( epsil, algrmn(ns) )
 
c--Debut la de boucle externe sur l'indice de niveau k :
C$DIR SPP LOOP_PARALLEL
      do 480 k=ks1,ks2
c----
        call alphdkc(alphax(1,1,k),alphay(1,1,k),grdmin,k,ns)
c--Fin de la boucle externe sur l'indice de niveau k .
 480  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine alphdec -
      end
