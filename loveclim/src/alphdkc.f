












      subroutine alphdkc(alphkx, alphky, grdmin, k, ns)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  calcule les taux de decentrement alphax & alphay dans les 2 directions,
c  pour le niveau "k", et les passent par argument (alphkx & alphky).
c -> uniquement cas des Coeff. "alpha" visant a supprimer le Depassement .
c   en sortie : multiplie par |u1.cmy/2|, |v2.cmx/2|
c  modif : 25/05/99
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'isoslope.com'
      include 'comunit.h'
 
c--variables locales equivalentes :
C     equivalence ( phiax(1,1)  , phihhh(1,1,5) )
C     equivalence ( phiay(1,1)  , phihhh(1,1,6) )
C     equivalence ( alphxx(1,1)  , phihhh(1,1,7) )
C     equivalence ( alphyy(1,1)  , phihhh(1,1,8) )
 
c--dummy variables :
      dimension alphkx(imax,jmax), alphky(imax,jmax)
 
c--variables locales :
      dimension phiax(imax,jmax) , phiay(imax,jmax)
      dimension alphxx(imax,jmax) , alphyy(imax,jmax)
      dimension phicx(ixjmax) , phicy(ixjmax)
      equivalence ( phicx(1) , phiax(1,1) )
      equivalence ( phicy(1) , phiay(1,1) )
 
      if (numit.eq.nstart) then
       if (k.eq.ks2) then
Ca2     write(iuo+66,'(A,I3,F10.6,1PE14.6)') ' alphdkc : sans alphah, A*A ;'
        write(iuo+66,'(A,I3,F10.6,1PE14.6)') ' alphdkc : sans alphah, |A| ;'
     &              //' ns,alphgr,algrmn =', ns, -alphgr(ns), grdmin
       endif
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Calcul des Coeff. "alpha" visant a supprimer le Depassement .   |
c-----------------------------------------------------------------------
 
c--calcul des gradients dans les 2 directions :
      do 410 j=js1,js2
       do 410 i=is1(j),1+is2(j)
        phiax(i,j) = (scal(i,j,k,ns) - scal(i-1,j,k,ns))
     &             * tms(i,j,k) * tms(i-1,j,k)
 410  continue
      do 420 j=js1,1+js2
       do 420 i=isf1(j),isf2(j)
        phiay(i,j) = (scal(i,j,k,ns) - scal(i,j-1,k,ns))
     &             * tms(i,j,k) * tms(i,j-1,k)
 420  continue
 
c--virage dans les coins (4 types de coins) :
c- coin 1 : tmu(i,j) = 1.
      do 441 n=1,n1coin(k)
        phicx( i1coin(n,k)+1 )    = -phicy( i1coin(n,k) )
        phicy( i1coin(n,k)+imax ) = -phicx( i1coin(n,k) )
 441  continue
c- coin 2 : tmu(i+1,j) = 1.
      do 442 n=1,n2coin(k)
        phicx( i2coin(n,k) )      =  phicy( i2coin(n,k) )
        phicy( i2coin(n,k)+imax ) =  phicx( i2coin(n,k)+1 )
 442  continue
c- coin 3 : tmu(i,j+1) = 1.
      do 443 n=1,n3coin(k)
        phicx( i3coin(n,k)+1 )    =  phicy( i3coin(n,k)+imax )
        phicy( i3coin(n,k) )      =  phicx( i3coin(n,k) )
 443  continue
c- coin 4 : tmu(i+1,j+1) = 1.
      do 444 n=1,n4coin(k)
        phicx( i4coin(n,k) )      = -phicy( i4coin(n,k)+imax )
        phicy( i4coin(n,k) )      = -phicx( i4coin(n,k)+1 )
 444  continue
 
c--calcul des coeff alphxx & alphyy :
      do 460 j=js1,js2
       do 460 i=is1(j),is2(j)
        alphxx(i,j) = -alphgr(ns)
     &              * abs(phiax(i+1,j) - phiax(i,j))
Ca2     alphx1      =    (phiax(i+1,j) - phiax(i,j))
     &            / ( abs(phiax(i+1,j) + phiax(i,j)) + grdmin )
Ca2     alphxx(i,j) =  -alphgr(ns) * alphx1 * alphx1
        alphyy(i,j) = -alphgr(ns)
     &              * abs(phiay(i,j+1) - phiay(i,j))
Ca2     alphy1      =    (phiay(i,j+1) - phiay(i,j))
     &            / ( abs(phiay(i,j+1) + phiay(i,j)) + grdmin )
Ca2     alphyy(i,j) = -alphgr(ns) * alphy1 * alphy1
 460  continue
 
c--raccord cyclique et autre (alphxx, alphyy) :
      if (ltest.ge.1) then
        do 470 j=jcl1,jcl2
          alphxx(ims1-1,j) = alphxx(ims2,j)
          alphxx(ims2+1,j) = alphxx(ims1,j)
 470    continue
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
      do 510 j=js1,js2
        do 510 i=is1(j),is2(j)+1
Cis0      u1pd2 = abs(0.25 * (u(i,j,k) + u(i,j+1,k)) )
          u1pd2 = abs(0.25 * (u(i,j,k) + u(i,j+1,k))
     &               +0.5  *  uiso(i,j,k)            )
          alpha = smx(i,j,1) * cc2x * u1pd2 + alphmi(k)
C         cnx = smx(i,j,1) * cc2x * u1d2
C         alpha = abs(cnx) + alphmi(k)
C    &          + alphah(1) * (2.0 - tmu(i,j,k) - tmu(i,j+1,k))
          alphkx(i,j) = cmy(i,j,1) * u1pd2 *
     &       min(one, max(alpha, alphxx(i,j), alphxx(i-1,j)))
 510  continue
 
c--calcul du taux de decentrement dans la direction y :
      cc2y  = 2. * unsdy * dts(k)
      do 530 j=js1,1+js2
        do 530 i=isf1(j),isf2(j)
Cis0      v2pd2 = abs( 0.25 * (v(i,j,k) + v(i+1,j,k)) )
          v2pd2 = abs( 0.25 * (v(i,j,k) + v(i+1,j,k))
     &                +0.5  *  viso(i,j,k)            )
          alpha = smy(i,j,2) * cc2y * v2pd2 + alphmi(k)
C         cny = smy(i,j,2) * cc2y * v2d2
C         alpha = abs(cny) + alphmi(k)
C    &          + alphah(2) * (2.0 - tmu(i,j,k) - tmu(i+1,j,k))
          alphky(i,j) = cmx(i,j,2) * v2pd2 *
     &       min(one, max(alpha, alphyy(i,j), alphyy(i,j-1)))
 530  continue
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine alphdkc -
      end
