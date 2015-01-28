












      subroutine alph2dc(alphax,alphay,ns)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  calcule les taux de decentrement alphax & alphay dans les 2 directions .
c   en sortie : multiplie par |u1.cmy/2|, |v2.cmx/2|
c  generalisation 2 D du calcul "suppresion du depassement".
c  modif : 30/12/97
c  ppmodif: 19-03-97: gm90
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'isoslope.com'
 
c--variables locales equivalentes :
C     equivalence ( grdx1(1,1)  , phihhh(1,1,5) )
C     equivalence ( grdy2(1,1)  , phihhh(1,1,6) )
 
c--dummy variables :
      dimension alphax(imax,jmax,kmax), alphay(imax,jmax,kmax)
 
c--variables locales :
      dimension grdx1(imax,jmax) , grdy2(imax,jmax)
      dimension u1cd2(imax,jmax), v2cd2(imax,jmax), unsdiv(imax,jmax)
      dimension u0ct(imax,jmax), v0ct(imax,jmax)
      dimension phx0e(imax,jmax), phx0w(imax,jmax)
      dimension phy0s(imax,jmax), phy0n(imax,jmax)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  0 ) preparation des coeff. intervenant dans les flux d'advection .  |
c-----------------------------------------------------------------------
 
      alpgrs = -1. - alphgr(ns)
      grdmin = max( epsil, algrmn(ns) )
      if(numit.eq.nstart) then
       if ( alpgrs.lt.epsil ) then
        write(iuo+66,*) 'ARRET dans "alph2dc" ; param alphgr mauvais !'
        stop
       endif
 
        write(iuo+66,'(A,I3,F10.6,1PE14.6)') ' alph2dc : sans alphah, |A| ;'
Ca2     write(iuo+66,'(A,I3,F10.6,1PE14.6)') ' alph2dc : sans alphah, A*A ;'
     &             //' ns,alphgr,algrmn =', ns, alpgrs, grdmin
      endif
 
c--Debut de la boucle externe sur l'indice de niveau k :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,j,phx0w,phx0e,phy0s,phy0n)
C$DIR SPP LOOP_PRIVATE(u1cd2,v2cd2,grdx1,grdy2,u0ct,v0ct)
C$DIR SPP LOOP_PRIVATE(u0cdiv,v0cdiv,unsdiv,grdx0d,grdy0d,phx0t,phy0t)
C$DIR SPP LOOP_PRIVATE(alpmng,cc2xg,cc2yg,u1pcd2,v2pcd2,alpha)
C$DIR SPP LOOP_PRIVATE(alphyn,alphys,alphxe,alphxw)
      do 600 k=ks1,ks2
c-----
 
c--Initialisation :
      do 30 j=1,jmax
       do 30 i=1,imax
        phx0w(i,j) = 0.
        phx0e(i,j) = 0.
        phy0s(i,j) = 0.
        phy0n(i,j) = 0.
 30   continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Calcul des Coeff. "alpha" visant a supprimer le Depassement .   |
c-----------------------------------------------------------------------
 
c--calcul des vitesses(flux) u1 & v2 et des gradients dans les 2 directions :
c- Attention : suppose dx = dy (si non, doivent etre pris en compte)
      do 410 j=js1,js2
       do 410 i=is1(j),1+is2(j)
        u1cd2(i,j) = 0.25 * cmy(i,j,1) * (u(i,j,k) + u(i,j+1,k))
     &              +0.5  * cmy(i,j,1) *  uiso(i,j,k)
        grdx1(i,j) = (scal(i,j,k,ns) - scal(i-1,j,k,ns))
C    &             * tms(i,j,k) * tms(i-1,j,k)
     &             * tms(i,j,k) * tms(i-1,j,k) * smx(i,j,1)
 410  continue
      do 420 j=js1,1+js2
       do 420 i=isf1(j),isf2(j)
        v2cd2(i,j) = 0.25 * cmx(i,j,2) * (v(i,j,k) + v(i+1,j,k))
     &              +0.5  * cmx(i,j,2) *  viso(i,j,k)
        grdy2(i,j) = (scal(i,j,k,ns) - scal(i,j-1,k,ns))
C    &             * tms(i,j,k) * tms(i,j-1,k)
     &             * tms(i,j,k) * tms(i,j-1,k) * smy(i,j,2)
 420  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--repartition des flux :
      do 450 j=js1,js2
       do 450 i=is1(j),is2(j)
        u0ct(i,j) = min(zero, max(u1cd2(i,j), u1cd2(i+1,j)) )
     &            + max(zero, min(u1cd2(i,j), u1cd2(i+1,j)) )
        v0ct(i,j) = min(zero, max(v2cd2(i,j), v2cd2(i,j+1)) )
     &            + max(zero, min(v2cd2(i,j), v2cd2(i,j+1)) )
        u0cdiv = u1cd2(i,j) - u1cd2(i+1,j)
        v0cdiv = v2cd2(i,j) - v2cd2(i,j+1)
        unsdiv(i,j) = 1. / sign( max(epsil,abs(u0cdiv),abs(v0cdiv)),
     &                           (u0cdiv - v0cdiv) )
c-
C       u0ct(i,j) = 0.
C       v0ct(i,j) = 0.
        phx0w(i,j) = u1cd2( i, j) - u0ct(i,j)
        phx0e(i,j) = u1cd2(i+1,j) - u0ct(i,j)
        phy0s(i,j) = v2cd2(i, j ) - v0ct(i,j)
        phy0n(i,j) = v2cd2(i,j+1) - v0ct(i,j)
 450  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--calcul des coeff phx/y-e/w/n/s :
      do 460 j=js1,js2
       do 460 i=is1(j),is2(j)
c--grad_diag :
        grdx0d = unsdiv(i,j) *
     &    ( phx0w(i,j) * grdx1(i,j) + phx0e(i,j) * grdx1(i+1,j) )
        grdy0d = unsdiv(i,j) *
     &    ( phy0s(i,j) * grdy2(i,j) + phy0n(i,j) * grdy2(i,j+1) )
c--coef alpha transverse, |A| :
        phx0t = abs( u0ct(i,j)*(grdx1(i+1,j) - grdx1(i,j)) )
     &         / ( grdmin + abs(grdx1(i+1,j) + grdx1(i,j)) )
        phy0t = abs( v0ct(i,j)*(grdy2(i,j+1) - grdy2(i,j)) )
     &         / ( grdmin + abs(grdy2(i,j+1) + grdy2(i,j)) )
c--coef alpha transverse, A*A :
Ca2     phx0t = (grdx1(i+1,j) - grdx1(i,j))
Ca2  &   / ( abs(grdx1(i+1,j) + grdx1(i,j)) + grdmin )
Ca2     phx0t = abs(u0ct(i,j)) * phx0t * phx0t
Ca2     phy0t = (grdy2(i,j+1) - grdy2(i,j))
Ca2  &   / ( abs(grdy2(i,j+1) + grdy2(i,j)) + grdmin )
Ca2     phy0t = abs(v0ct(i,j)) * phy0t * phy0t
c--coef alpha complet, |A| :
C       grdy0d = grdx1(i+1,j)          ! Debugg
C       grdy0d = -grdx1(i,j)           ! Debugg
C       grdx0d = grdy2(i,j)            ! Debugg
C       grdx0d = -grdy2(i,j+1)         ! Debugg
        phx0w(i,j) = phx0t + abs( phx0w(i,j)*(grdx1( i, j) - grdy0d) )
     &                       / ( grdmin + abs(grdx1( i, j) + grdy0d) )
        phx0e(i,j) = phx0t + abs( phx0e(i,j)*(grdx1(i+1,j) + grdy0d) )
     &                       / ( grdmin + abs(grdx1(i+1,j) - grdy0d) )
        phy0s(i,j) = phy0t + abs( phy0s(i,j)*(grdy2(i, j ) + grdx0d) )
     &                       / ( grdmin + abs(grdy2(i, j ) - grdx0d) )
        phy0n(i,j) = phy0t + abs( phy0n(i,j)*(grdy2(i,j+1) - grdx0d) )
     &                       / ( grdmin + abs(grdy2(i,j+1) + grdx0d) )
c--coef alpha complet, A*A :
Ca2     alphxw = (grdx1( i, j) - grdy0d)
Ca2  &    / ( abs(grdx1( i, j) + grdy0d) + grdmin )
Ca2     phx0w(i,j) = phx0t + abs(phx0w(i,j)) * alphxw * alphxw
Ca2     alphxe = (grdx1(i+1,j) + grdy0d)
Ca2  &    / ( abs(grdx1(i+1,j) - grdy0d) + grdmin )
Ca2     phx0e(i,j) = phx0t + abs(phx0e(i,j)) * alphxe * alphxe
Ca2     alphys = (grdy2(i, j ) + grdx0d)
Ca2  &    / ( abs(grdy2(i, j ) - grdx0d) + grdmin )
Ca2     phy0s(i,j) = phy0t + abs(phy0s(i,j)) * alphys * alphys
Ca2     alphyn = (grdy2(i,j+1) - grdx0d)
Ca2  &    / ( abs(grdy2(i,j+1) + grdx0d) + grdmin )
Ca2     phy0n(i,j) = phy0t + abs(phy0n(i,j)) * alphyn * alphyn
 460  continue
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--(demi) raccord cyclique pour phx_e/w :
      do 470 j=jcl1,jcl2
        phx0e(ims1-1,j) = phx0e(ims2,j)
        phx0w(ims2+1,j) = phx0w(ims1,j)
 470  continue
 
      if (ltest.eq.3) then
c--Bering : (demi) raccord pour phy_n/s :
        phy0s(iberp-1,jberp) = phy0n( ibera, jbera-1)
        phy0s( iberp, jberp) = phy0n(ibera-1,jbera-1)
        phy0s(ibera-1,jbera) = phy0n( iberp, jberp-1)
        phy0s( ibera, jbera) = phy0n(iberp-1,jberp-1)
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Calcul des Taux de decentrement a partir des Coeff "alpha" :    |
c-----------------------------------------------------------------------
 
      alpmng = alphmi(k) / alpgrs
c--calcul du taux de decentrement dans la direction x :
      cc2xg = 2. * unsdx * dts(k) / alpgrs
      do 510 j=js1,js2
       do 510 i=is1(j),is2(j)+1
         u1pcd2 = abs(u1cd2(i,j))
         alpha =  u1pcd2 *
     &     ( alpmng + smx(i,j,1) * smy(i,j,1) * cc2xg * u1pcd2 )
C    &     ( alphmi(k) + smx(i,j,1) * smy(i,j,1) * cc2x * u1pcd2 )
         alphax(i,j,k) = min( u1pcd2,
     &          alpgrs * max(alpha, phx0e(i-1,j), phx0w(i,j)) )
C    &      max(alpha, alpgrs * max(phx0e(i-1,j), phx0w(i,j)) ))
 510  continue
 
c--calcul du taux de decentrement dans la direction y :
      cc2yg  = 2. * unsdy * dts(k) / alpgrs
      do 530 j=js1,1+js2
       do 530 i=isf1(j),isf2(j)
         v2pcd2 = abs(v2cd2(i,j))
         alpha =  v2pcd2 *
     &      ( alpmng + smx(i,j,2) * smy(i,j,2) * cc2yg * v2pcd2 )
C    &      ( alphmi(k) + smx(i,j,2) * smy(i,j,2) * cc2y * v2pcd2 )
         alphay(i,j,k) = min( v2pcd2,
     &          alpgrs * max(alpha, phy0n(i,j-1), phy0s(i,j)) )
C    &      max(alpha, alpgrs * max(phy0n(i,j-1), phy0s(i,j)) ))
 530  continue
 
c--Fin de la boucle externe sur l'indice de niveau k .
 600  continue
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine alph2dc -
      end
