












      subroutine conti3d

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- calcul de la vitesse verticale w depuis la surface jusqu'au fond.
c  modif : 06/10/98

      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
C-AM - pour LOCH: (u,v,w) au meme instant que T & S
C      dynam3 en common avec loch.com
      real*8 uloch,vloch,wloch
      common /dynam3/uloch(imax,jmax,kmax),vloch(imax,jmax,kmax),
     &               wloch(imax,jmax,kmax)
C-AM - pour dilution icesheets
      logical flgveg,flgicb,flgisma,flgismg
      real*8 deltamass,massGold,massAold
      real*8 frwism(imax,jmax)
c
      common /ec_coupl/flgveg,flgicb,flgisma,flgismg
      COMMON/AG2CLIO2/frwism,deltamass,massGold,massAold

c--variables locales equivalentes (pour raison de place memoire) :
      dimension phixj(imax)
C     dimension phix(imax,jmax), phiy(imax,jmax)
C     equivalence ( phix(1,1) , phihhh(1,1,1) )
C     equivalence ( phiy(1,1) , phihhh(1,1,2) )

C-AM - terme rappel sur la sal. pour FW Flx -> voir thersf.f


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      unsdts = 1.0 / dts(ks2)

c--Debut de la boucle externe sur l'indice de latitude j :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,k,phixj,dzd2)
      do 105 j=js1,js2
c-----

c--initialisation :
      do 100 i=is1(j),is2(j)
C       w(i,j,ks1) = 0.0
Cdbug     if (abs(scal(i,j,ks2,2)).lt.epsil) then
Cdbug      write(iuo+66,*) 'ARRET : conti3d, scal(i,j,ks2,2) too small'
Cdbug      write(iuo+66,*) 'scal(i=',i,',j=',j,')=',scal(i,j,ks2,2)
Cdbug      stop
Cdbug     endif
c- 
C-AM
        w(i,j,ks2+1) = tms(i,j,ks2) * unsdts * phiss(i,j,0)
c- 
c- phimnx(0,1) => limits Fresh Water Flux : Unit = m/s !!
Cic0    w(i,j,ks2+1) =
Cic0 &    max( phimnx(i,j,0,0), min(phimnx(i,j,1,0),w(i,j,ks2+1)) )
Cfd0    w(i,j,ks2+1) = 0.0
 100   continue
 105  continue

c mean value of w
        zflux0=0.0
        do j=js1,js2
         do i=is1(j),is2(j)
           zflux0=zflux0+aire(i,j)*tms(i,j,ks2)*w(i,j,ks2+1)
         enddo
        enddo
        zflux0=zflux0/zurfow
        write(99,*) 'zflux',zflux0,w(15,15,ks2+1),zfluxm
c
c ---> Modifications AM ---> see flucor
c

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Bilan de l'Eq. Continuite et Somme directement dans w .         |
c-----------------------------------------------------------------------

      do 300 j=js1,js2
c--boucle sur l'indice de niveau, "k" :
      do 290 k=ks2,ks1,-1

      do 230 i=is1(j),is2(j)+1
        phixj(i) = cmy(i,j,1) * (u(i,j+1,k)+u(i,j,k))
 230  continue
C       phiy(i,j) = cmx(i,j,2) * (v(i+1,j,k)+v(i,j,k))

      dzd2 = 0.5 * dz(k)
      do 250 i=is1(j),is2(j)
        w(i,j,k) = w(i,j,k+1) + dzd2 * smxy(i,j,0) *
C    &           ( unsdx * (phix(i+1,j)-phix(i,j))
C    &           + unsdy * (phiy(i,j+1)-phiy(i,j)) )
     &           ( unsdx * (phixj(i+1)-phixj(i))
     &           + unsdy * ( cmx(i,j+1,2)*(v(i+1,j+1,k)+v(i,j+1,k))
     &                     - cmx(i, j, 2)*(v(i+1, j, k)+v(i, j, k)) ))
C       w(i,j,k+1) = w(i,j,k+1) * tms(i,j,k)
        w(i,j,k) = w(i,j,k) * tms(i,j,k)
 250  continue

c-----
 290  continue
c--fin boucle sur "k".


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Fin de la boucle externe sur l'indice de latitude j .
 300  continue


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Prepare Flux en surface .                                       |
c-----------------------------------------------------------------------

c--Valeur du scalaire associee au Flux de masse en surf. (= Pluie / Evap.)
Cic0  do 330 ns=1,nsmax
Cic0    if (scpme(ns).eq.spvr) then
Cic0      do 320 j=js1,js2
Cic0       do 320 i=is1(j),is2(j)
Cic0         scs(i,j,ns) = scal(i,j,ks2,ns)
 320      continue
Cic0    endif
 330  continue



c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C-AM - pour LOCH: (u,v,w) au meme instant que T & S
c-----------------------------------------------------------------------
      do 400 k=1,kmax
       do 400 j=1,jmax
        do 400 i=1,imax
         uloch(i,j,k)=u(i,j,k)
         vloch(i,j,k)=v(i,j,k)
         wloch(i,j,k)=w(i,j,k)
  400 continue

      return

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine conti3d -
      end
