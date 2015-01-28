












      subroutine uve
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c**AVANCEMENT EXPLICITE D'UN PAS DE TEMPS BAROCLINE DE U et V**
c  avec taux advection verticale implicite. + termes visqueux en 1/Rt.
c    et taux diffusion verticale implicite.
c  advection suivant la verticale :  w(k) * [ V(k-1) + V(k) ] / 2
c  combinee avec version traitement visco.barotrope type "s" (cf notes)
c  Suprime termes visqueux-metriques - Calcule Tension de Fond.
c Cfcc [Cfc0] => option : Avec Coeff Coriolis Reciproque = -2.Omega.cos(phi)
c  modif : 19/08/96
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
 
c--variables locales equivalentes :
      dimension phizu(imax,jmax,kmax+1), phizv(imax,jmax,kmax+1)
      equivalence ( phizu(1,1,1) , phizzz(1,1,1,1) )
      equivalence ( phizv(1,1,1) , phizzz(1,1,1,2) )
 
      dimension phiaxu(imax,jmax), phiaxv(imax,jmax)
      dimension phidxu(imax,jmax), phidxv(imax,jmax)
      dimension phiayu(imax,jmax), phiayv(imax,jmax)
      dimension phidyu(imax,jmax), phidyv(imax,jmax)
      dimension detadx(imax,jmax), detady(imax,jmax)
c--pour raison de place memoire :
C     equivalence ( phiaxu(1,1) , phihhh(1,1,1) )
C     equivalence ( phiaxv(1,1) , phihhh(1,1,2) )
C     equivalence ( phidxu(1,1) , phihhh(1,1,3) )
C     equivalence ( phidxv(1,1) , phihhh(1,1,4) )
C     equivalence ( phiayu(1,1) , phihhh(1,1,5) )
C     equivalence ( phiayv(1,1) , phihhh(1,1,6) )
C     equivalence ( phidyu(1,1) , phihhh(1,1,7) )
C     equivalence ( phidyv(1,1) , phihhh(1,1,8) )
C     equivalence ( detadx(1,1) , phihhh(1,1,9)  )
C     equivalence ( detady(1,1) , phihhh(1,1,10) )
 
      dimension cuhe(imax,jmax) ,  cvhe(imax,jmax)
C     equivalence ( cuhe(1,1)   , phihhh(1,1,9) )
C     equivalence ( cvhe(1,1)   , phihhh(1,1,10))
 
 
c--variables locales :
      dimension wud2dz(imax)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Forcage lie a l'elevation ; Calcul de la Tension de Fond.       |
c-----------------------------------------------------------------------
 
      ccgdx = gpes * uns2dx
      ccgdy = gpes * uns2dy
      do 130 j=ju1,ju2
       do 130 i=iu1(j),iu2(j)
        detadx(i,j) = ccgdx * smx(i,j,3) *
     &              ((eta(i-1,j-1)-eta(i,j))+(eta(i-1,j)-eta(i,j-1)))
        detady(i,j) = ccgdy * smy(i,j,3) *
     &              ((eta(i-1,j-1)-eta(i,j))-(eta(i-1,j)-eta(i,j-1)))
 130  continue
 
      if (cdbot.ne.zero) then
c-----
c--Tension de fond : TauX,Y = Cdrag * |V(bot.)| * u,v(bot.) = Coeff * u,v(bot.)
c  calcul du Coeff. (dans avudz(-,-,1))(pour resolution implicite).
      do 150 j=ju1,ju2
       do 150 i=iu1(j),iu2(j)
         avudz(i,j,1) = sqrt( u(i,j,kfu(i,j))*u(i,j,kfu(i,j))
     &                      + v(i,j,kfu(i,j))*v(i,j,kfu(i,j)) ) * cdbot
         phifu(i,j) = -avudz(i,j,1)*u(i,j,kfu(i,j))
         phifv(i,j) = -avudz(i,j,1)*v(i,j,kfu(i,j))
 150  continue
c-----
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Calcul des Flux Verticaux Explicite ADvection et DiFfusion.    |
c----------------------------------------------------------------------
 
      if (txiadu.eq.one .and. txidfu.eq.one) goto 299
 
      ccexa = (1. - txiadu) * 0.125
      ccexd =  1. - txidfu
c--Debut de la 1ere boucle externe sur l'indice de niveau k :
      do 280 k=ku1+1,ku2
c-----
 
c--calcul des flux verticaux d'advection et de diffusion de qqmvt :
C      ccexak = ccexa * unsdzw(k)
       do 250 j=ju1,ju2
c- construction vitesses verticales (a un fact.mult. pres ) :
        do 210 i=iu1(j),iu2(j)
C         wud2dz(i)= ccexak
          wud2dz(i)= ccexa
     &          *(w(i,j,k)+w(i-1,j,k)+w(i,j-1,k)+w(i-1,j-1,k))
 210    continue
c- calcul des flux :
        do 230 i=iu1(j),iu2(j)
          ccdif =  avudz(i,j,k) * ccexd
          phizu(i,j,k) = tmu(i,j,k-1) * (
C    &            (dz(k)*u(i,j,k-1)+dz(k-1)*u(i,j,k))*wud2dz(i)
     &            ( u(i,j,k-1) + u(i,j,k) ) * wud2dz(i)
     &          + ( u(i,j,k-1) - u(i,j,k) ) * ccdif )
          phizv(i,j,k) = tmu(i,j,k-1) * (
C    &       (dz(k)*v(i,j,k-1)+dz(k-1)*v(i,j,k))*wud2dz(i)
     &       ( v(i,j,k-1) + v(i,j,k) ) * wud2dz(i)
     &          + ( v(i,j,k-1) - v(i,j,k) ) * ccdif )
 230    continue
 250  continue
c--fin du calcul des flux verticaux d'advection et de diffusion.
 
c--Fin de la 1ere boucle externe sur l'indice de niveau k .
 280  continue
 
c--calcul des flux V. Explicites termine.
 299  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      ww3 = 0.
      ccxdif = ahu * unsdx
      ccydif = ahu * unsdy
      corfac = 2.*(1.-txifcu)
c--Debut de la 2nd  boucle externe sur l'indice de niveau k :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,j,u2d2,ccdif,phiaxu,phiaxv,phidxu,phidxv)
C$DIR SPP LOOP_PRIVATE(v1d2,phiayu,phiayv,phidyu,phidyv)
C$DIR SPP LOOP_PRIVATE(cud,cvd,cuhe,cvhe,cu,cv,ww3)
      do 500 k=ku2,ku1,-1
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Resolution Complete des Flux Horizontaux .        |
c---------------------------------------------------------
 
c**place reservee pour le calcul des taux de decentrement**
 
c--calcul des flux (advection "a", diffusion "d") dans les directions x et y :
 
      do 310 j=ju1,ju2
       do 310 i=iu1(j)-1,iu2(j)
        u2d2 = cmy(i,j,2) * 0.25 * (u(i,j,k) + u(i+1,j,k))
        ccdif = smxy(i,j,2) * ccxdif
        phiaxu(i,j) = u2d2  * (u(i,j,k)+u(i+1,j,k))
     &              + alphxu * abs(u2d2) * (u(i,j,k)-u(i+1,j,k))
        phiaxv(i,j) = u2d2  * (v(i,j,k)+v(i+1,j,k))
     &              + alphxv * abs(u2d2) * (v(i,j,k)-v(i+1,j,k))
        phidxu(i,j) = ccdif * (u(i,j,k)-u(i+1,j,k))
        phidxv(i,j) = ccdif * (v(i,j,k)-v(i+1,j,k))
 310  continue
 
      do 320 j=ju1-1,ju2
       do 320 i=iuf1(j),iuf2(j)
        v1d2 = cmx(i,j,1) * 0.25 * (v(i,j,k)+v(i,j+1,k))
        ccdif = cmxy(i,j,1) * ccydif
        phiayu(i,j) = v1d2  * (u(i,j,k)+u(i,j+1,k))
     &              + alphyu * abs(v1d2) * (u(i,j,k)-u(i,j+1,k))
        phiayv(i,j) = v1d2  * (v(i,j+1,k)+v(i,j,k))
     &              + alphyv * abs(v1d2) * (v(i,j,k)-v(i,j+1,k))
        phidyu(i,j) = ccdif * (u(i,j,k)-u(i,j+1,k))
        phidyv(i,j) = ccdif * (v(i,j,k)-v(i,j+1,k))
 320  continue
 
c--fin de calcul des flux horizontaux.
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Bilan de tous les Flux calcules Explicitement .   |
c---------------------------------------------------------
 
      do 420 j=ju1,ju2
       do 420 i=iu1(j),iu2(j)
        cud = smxy(i,j,3)*( unsdx*(phidxu(i-1,j)-phidxu(i,j))
     &                    + unsdy*(phidyu(i,j-1)-phidyu(i,j)) )
        cvd = smxy(i,j,3)*( unsdx*(phidxv(i-1,j)-phidxv(i,j))
     &                    + unsdy*(phidyv(i,j-1)-phidyv(i,j)) )
 
c- test 2nd terme corriolis : du/dt = - 2 . Omega . cos(Phi) . w
Cfcc    ww3 = (w(i-1,j-1,k) + w(i,j,k)) + (w(i-1,j,k) + w(i,j-1,k))
Cfcc &    + (w(i-1,j-1,k+1)+w(i,j,k+1)) + (w(i-1,j,k+1)+w(i,j-1,k+1))
        cuhe(i,j) = tmu(i,j,k) * ( cud
Cfcc &      + fcucor(i,j) * ww3
     &      + smx(i,j,3)*uns2dx
     &          *((q(i-1,j-1,k)-q(i,j,k))+(q(i-1,j,k)-q(i,j-1,k)))
     &      + smxy(i,j,3)*( unsdx*(phiaxu(i-1,j)-phiaxu(i,j))
     &                    + unsdy*(phiayu(i,j-1)-phiayu(i,j))
     &                    - cmxdy(i,j)*u(i,j,k)*v(i,j,k)
     &                    + cmydx(i,j)*v(i,j,k)*v(i,j,k)     ))
        cvhe(i,j) = tmu(i,j,k) * ( cvd
Cfcc &      + fcvcor(i,j) * ww3
     &      + smy(i,j,3)*uns2dy
     &          *((q(i-1,j-1,k)-q(i,j,k))-(q(i-1,j,k)-q(i,j-1,k)))
     &      + smxy(i,j,3)*( unsdx*(phiaxv(i-1,j)-phiaxv(i,j))
     &                    + unsdy*(phiayv(i,j-1)-phiayv(i,j))
     &                    - cmydx(i,j)*u(i,j,k)*v(i,j,k)
     &                    + cmxdy(i,j)*u(i,j,k)*u(i,j,k)     ))
 
        fub(i,j,k) = dz(k)*(cuhe(i,j)-cud)
        fvb(i,j,k) = dz(k)*(cvhe(i,j)-cvd)
 420  continue
 
      do 440 j=ju1,ju2
       do 440 i=iu1(j),iu2(j)
        cu = cuhe(i,j) + detadx(i,j) + corfac*fs2cor(i,j)*v(i,j,k)
     &          + unsdz(k)*(phizu(i,j,k)-phizu(i,j,k+1))
        cv = cvhe(i,j) + detady(i,j) - corfac*fs2cor(i,j)*u(i,j,k)
     &          + unsdz(k)*(phizv(i,j,k)-phizv(i,j,k+1))
 
        u(i,j,k) = u(i,j,k) + dtu * cu
        v(i,j,k) = v(i,j,k) + dtu * cv
 440  continue
 
c--Fin de la 2nd  boucle externe sur l'indice de niveau k .
 500  continue
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine uve -
      end
