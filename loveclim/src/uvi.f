












      subroutine uvi(j)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c**AVANCEMENT IMPLICITE D'UN PAS DE TEMPS BAROCLINE DE U ET V**
c  avec taux advection verticale implicite : txiadu.
c    et taux diffusion verticale implicite : txidfu.
c  advection suivant la verticale :  w(k) * [ V(k-1) + V(k) ] / 2
c  modif : 09/01/95
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
 
c--variables locales equivalentes :
C     dimension detadx(imax,jmax), detady(imax,jmax)
c--pour transfert uvb --> uvi :
C     equivalence ( detadx(1,1), phihhh(1,1,9)  )
C     equivalence ( detady(1,1), phihhh(1,1,10) )
 
c--variables locales :
C     dimension ccwup(kmax), ccwdn(kmax)
      dimension wudtd2(imax,kmax)
      dimension aa(imax,kmax), bb(imax,kmax), cc(imax,kmax)
      dimension ffcor(imax)
c- variables complexes :
      dimension acplex(imax,kmax), bcplex(imax,kmax), vcplex(imax,kmax)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  0 ) Mise en place de coeff dependant de k :
c-------------------
 
      cctdif = dtu * txidfu
      corfac = 2.*txifcu*dtu
 
C     do 10 k=ks1+1,ks2
C       ccwup(k) =  dz(k-1) * unsdzw(k)
C       ccwdn(k) =  dz(k)   * unsdzw(k)
C10   continue
 
      ccz  = dtu * unsdz(ku2)
 
C     do 500 j=ju1,ju2
c**debut boucle sur la latitude** <-- Mise a l'Exterieur.
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Prise en compte des Flux fond & surface .                       |
c-----------------------------------------------------------------------
 
c-  gradient d'elevation  <- deplace dans uve .
 
c--conditions aux limites :
c--fond :
c     do 130 i=iu1(j),iu2(j)
c       u(i,j,kfu(i,j))=u(i,j,kfu(i,j)) +dtu*unsdz(kfu(i,j))*phifu(i,j)
c       v(i,j,kfu(i,j))=v(i,j,kfu(i,j)) +dtu*unsdz(kfu(i,j))*phifv(i,j)
c130  continue
 
c--surface
      do 140 i=iu1(j),iu2(j)
        u(i,j,ku2) = u(i,j,ku2) - ccz * phisu(i,j)
        v(i,j,ku2) = v(i,j,ku2) - ccz * phisv(i,j)
 140  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Traitement Implicite : Calcul et Inversion matricielle .        |
c-----------------------------------------------------------------------
 
c**debut construction vitesses verticales (a un fact.mult. pres)
      cctimp = dtu * 0.125 * txiadu
      do 200 k=ku1+1,ku2
       do 200 i=iu1(j),iu2(j)
        wudtd2(i,k) = cctimp *
     &            ( (w(i-1,j-1,k)+w(i,j,k))+(w(i-1,j,k)+w(i,j-1,k)) )
C    &                (w(i,j,k)+w(i-1,j,k)+w(i,j-1,k)+w(i-1,j-1,k))
 200  continue
c**fin construction vitesses verticales**
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Traitement qqmvt semi-implicite : coriolis, advection, diffusion.
 
c--construction I : partie reelle (advect.+diff.) de la matrice du systeme .
      do 210 k=ku1,ku2
       do 210 i=iu1(j),iu2(j)
        aa(i,k)= 1.0
 210  continue
      do 220 k=ku1+1,ku2
       do 220 i=iu1(j),iu2(j)
        ccdif = cctdif * avudz(i,j,k)
C       wwup  = ccwup(k) * wudtd2(i,k)
C       wwdn  = ccwdn(k) * wudtd2(i,k)
c- effet du flux phiz(k) sur V(k) : a(k)*V(k) + b(k)*V(k-1)
        vt = unsdz(k) * tmu(i,j,k-1)
        bb(i,k) = vt * (-ccdif - wudtd2(i,k))
        aa(i,k) = aa(i,k) +  vt * ( ccdif - wudtd2(i,k))
C       bb(i,k) = vt * (-ccdif - wwdw)
C       aa(i,k) = aa(i,k) +  vt * ( ccdif - wwup)
c- effet du flux phiz(k) sur S(k-1) : a(k-1)*V(k-1) + c(k-1)*V(k)
        vt = unsdz(k-1) * tmu(i,j,k-1)
        aa(i,k-1) = aa(i,k-1) + vt * (ccdif + wudtd2(i,k))
        cc(i,k-1) = vt * (-ccdif + wudtd2(i,k))
C       aa(i,k) = aa(i,k) +  vt * ( ccdif + wwdn)
C       cc(i,k) = vt * (-ccdif + wwup)
 220  continue
 
      if (cdbot.ne.zero) then
c- effet de le tension de fond
      do 230  i=iu1(j),iu2(j)
        aa(i,kfu(i,j)) = aa(i,kfu(i,j)) + avudz(i,j,1)
 230  continue
      endif
 
c**fin construction de la matrice du systeme** --(partie reelle)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--construction II : partie imaginaire (coriolis) + decomposition LU :
 
c**debut decomposition LU** --(cf: Numerical Methods, pp 165-167)
c- calcul de 1/alpha(k) dans acplex(i,k) et de beta(k) dans bcplex(i,k)
      do 240 i=iu1(j),iu2(j)
       ffcor(i) = corfac * fs2cor(i,j)
       acplex(i,ku1) = 1.0 / DCMPLX( aa(i,ku1), ffcor(i) )
 240  continue
      do 250 k=ku1+1,ku2
       do 250 i=iu1(j),iu2(j)
        bcplex(i,k) = bb(i,k) * acplex(i,k-1)
        acplex(i,k) = 1.0 / ( DCMPLX( aa(i,k), ffcor(i) )
     &                      - bcplex(i,k)*cc(i,k-1) )
 250  continue
c**fin decompostion LU **
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Solution du Systeme Complex - Conversion Relle .                |
c-----------------------------------------------------------------------
 
c- charge f(k) dans vcplex(i,k) :
      do 300 k=ku1,ku2
       do 300 i=iu1(j),iu2(j)
        vcplex(i,k) = DCMPLX( u(i,j,k), v(i,j,k) )
 300  continue
 
c**debut substitutions avant et arriere**
c--calcul de g(k) dans vcplex(i,k) :
      do 311 k=ku1+1,ku2
       do 311 i=iu1(j),iu2(j)
        vcplex(i,k) = vcplex(i,k) - bcplex(i,k) * vcplex(i,k-1)
 311  continue
 
c--calcul de x(k) dans vcplex(i,k) :
      do 320 i=iu1(j),iu2(j)
       vcplex(i,ku2) = vcplex(i,ku2) * acplex(i,ku2)
 320  continue
      do 321 k=ku2-1,ku1,-1
       do 321 i=iu1(j),iu2(j)
        vcplex(i,k) = (vcplex(i,k)-cc(i,k)*vcplex(i,k+1)) * acplex(i,k)
 321  continue
 
c--mise en place des valeurs reelles dans u et v :
      do 330 k=ku1,ku2
       do 330 i=iu1(j),iu2(j)
        u(i,j,k) = tmu(i,j,k) * DREAL(vcplex(i,k))
        v(i,j,k) = tmu(i,j,k) * DIMAG(vcplex(i,k))
 330  continue
c**fin substitutions avant et arriere**
 
 500  continue
c**fin boucle sur la latitude**
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine uvi -
      end
