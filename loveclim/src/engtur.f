












      subroutine engtur
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  evolution of Turbulent Kinetic Energy  (diffusion is completly implicit)
c  modif : 15/04/99
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'ice.com'
      include 'comunit.h'
 
c--variables a conserver d'un appel a l'autre :
      common / xscali /
     &  z2alph(imax,jmax,2:kmax)
      common /aju0/ avs2mx(kmax),bvfmix
 
c--variables locales :
      dimension ff(imax,kmax)
      dimension aa(imax,kmax), bb(imax,kmax), cc(imax,kmax)
      dimension ccwabs(kmax), ccwcor(kmax), ccwdmy(kmax)
C     dimension ccwddn(kmax), ccwdup(kmax)
      dimension cczdt(kmax)
      dimension zaju(imax,jmax,kmax)
C     dimension ccvois(-kmax:kmax), alpha(imax,kmax)
 
c    CCW* coeff relatif a l'interface (= place des W).
c    CCZ* coeff relatif au centre des boites (= place de T,S)
c    CCWI : tient compte du taux Implicite.
c    CCWABS(k) intervient avec |w| ;
c              CCWUP(k) avec scal(k) ; CCWDowN(k) avec scal(k-1)
c  pour le taux de Decentrement calcule a partir du Nb de courant : CCWDUP(k)
c    avec dts(k), CCWDDowN(k) avec dts(k-1) ey CCWDMY(k) avec la moyenne des 2
c    CCZDT**(k) intervient pour le bilan des flux V, boite "k" ;
 
      if (numit.le.nstart) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

        avsmix = avnub(1)
        bvfmix = avnu0(1)
             if (lstab.eq.0) write(iuo+66,'(2(A,1PE10.3))')
     &   '          AVS +', avsmix, ' si bvf <', bvfmix
        do 30 k=ks1+1,ks2
         avs2mx(k) = unsdzw(k) * avsmix * 0.5
 30     continue
      endif

      if (cdbot.ne.zero) then
c-----
c--Tension de fond : TauX,Y = Cdrag * |V(bot.)| * u,v(bot.) = Coeff * u,v(bot.)
c  calcul du Coeff. (dans avudz(-,-,1))(pour resolution implicite).
      do 50 j=ju1,ju2
       do 50 i=iu1(j),iu2(j)
         avudz(i,j,1) = sqrt( u(i,j,kfu(i,j))*u(i,j,kfu(i,j))
     &                      + v(i,j,kfu(i,j))*v(i,j,kfu(i,j)) ) * cdbot
         phifu(i,j) = -avudz(i,j,1)*u(i,j,kfu(i,j))
         phifv(i,j) = -avudz(i,j,1)*v(i,j,kfu(i,j))
 50   continue

      do 60 j=js1,js2
       do 60 i=is1(j),is2(j)
        tm4u = tmu(i,j,ks2) + tmu(i+1,j+1,ks2)
     &       + tmu(i,j+1,ks2) + tmu(i+1,j,ks2) + epsil
        ust2b(i,j)=(sqrt(phifu(i,j)*phifu(i,j)+phifv(i,j)*phifv(i,j))
     &    + sqrt(phifu(i+1,j)*phifu(i+1,j)+phifv(i+1,j)*phifv(i+1,j))
     &    + sqrt(phifu(i,j+1)*phifu(i,j+1)+phifv(i,j+1)*phifv(i,j+1))
     &    + sqrt(phifu(i+1,j+1)*phifu(i+1,j+1)
     &            + phifv(i+1,j+1)*phifv(i+1,j+1))  ) /tm4u
 60   continue
C     write(144,*) ust2b(102,41)
         
c-----
      endif

 
      do 80 k=ks1,ks2
        cczdt(k) = dts(ks2) * unsdzw(k)
 80   continue
 
      if (lstab.eq.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Ajust.Conv. par augmentation de la Diffusion Verticale .
c-----------------------------------------------------------------------

c- Ajoute avsmix(/dzw) si bvf < bvfmix :
      demi = 0.5
Cajo  do k=1+ks1,ks2
Cajo   do j=js1,js2
Cajo    do i=is1(j),is2(j)
Cajo     avsdz(i,j,k) = avsdz(i,j,k)
Cajo &         + ( avs2mx(k) + sign(avs2mx(k), bvfmix-bvf(i,j,k)) )
Cajo     fqajc(i,j,k) = fqajc(i,j,k)
Cajo &         + ( demi + sign(demi, bvfmix-bvf(i,j,k)) )
Cajo    enddo
Cajo   enddo
Cajo  enddo

c- Convective adjustment if instability > kajul levels
      do 85 j=js1,js2
       do 85 i=is1(j),is2(j)
        zinsto = 0.0
	kstab=kajul
        do k=kajul+1,ks2
          zinztz=max(zero,sign(one,bvf(i,j,k) ) )
	  kstab=int(zinztz)*k+(1-int(zinztz))*kstab
          zinsto=min(zinsto+zinztz,one)
        enddo
        do k=ks2,kajul+1,-1
C         avsdz(i,j,k) = avsdz(i,j,k)
          zaju(i,j,k) =
     &         + ( avs2mx(k) + sign(avs2mx(k), bvfmix-bvf(i,j,k)) )
     &         * (1.0-zinsto)
         fqajc(i,j,k) = fqajc(i,j,k)
     &         + ( demi + sign(demi, bvfmix-bvf(i,j,k)) )
     &         * (1.0-zinsto)
        enddo
C       do k=1+ks1,kajul
        do k=1+ks1,kstab
C        avsdz(i,j,k) = avsdz(i,j,k)
         zaju(i,j,k) =
     &         + ( avs2mx(k) + sign(avs2mx(k), bvfmix-bvf(i,j,k)) )
         fqajc(i,j,k) = fqajc(i,j,k)
     &         + ( demi + sign(demi, bvfmix-bvf(i,j,k)) )
     &         * (1.0-zinsto)
C--Note :subsurface conv aju are not taken into account in fqajc
        enddo

 85   continue
      do j=js1,js2
       do i=is1(j),is2(j)
        do k=ks1,ks2
         avqdz(i,j,k)=avqdz(i,j,k)+zaju(i,j,k)
        enddo
       enddo
      enddo

      endif
c--Debut de la boucle externe sur l'indice de latitude j :
      do 500 j=js1,js2
c-----
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Building the matrix                                          |
c-----------------------------------------------------------------------
 
      epsil2=epsil*epsil

c--Boundary condition
      do 90 i=is1(j),is2(j)
       q2turb(i,j,kfs(i,j)) = max(q2tmin,6.51d0*ust2b(i,j))
       q2turb(i,j,ks2+1)    = max(q2tmin,6.51d0*ust2s(i,j))
C      q2turb(i,j,ks2+1)    = max(q2tmin,6.51d0*
C    &                   ust2s(i,j)*(one+sdvt(i,j)))
C    &                   ust2s(i,j)*(one+varfor*.375*sdvt(i,j)))
 90   continue

c--
 
      do 100 k=ks1+1,ks2
       do 100 i=is1(j),is2(j)

C       prod   =  tms(i,j,k-1)* 2 * dts(ks2) * (avudz(i,j,k)*
        prod   =  tms(i,j,k-1)* 2 * dts(ks2) * (avuloc(i,j,k)*
     &            dzw(k)*tm2tur(i,j,k) + 
     &            max(zero,-avsdz(i,j,k)*dzw(k)*bvf(i,j,k)) )

        des    =  2*dts(ks2)*(sqrt(q2turb(i,j,k))/(16.6*vlturb(i,j,k))
     &            + max(zero,avsdz(i,j,k)*dzw(k)*bvf(i,j,k)
     &                     /q2turb(i,j,k)) )

        aa(i,k) = 1.0 + tms(i,j,k-1)*des
        
        ff(i,k) = q2turb(i,j,k) + prod

 100  continue
 
c-N.B.-: avsdz =  Dif.Vert. / dz !!
      do 110 k=ks1+1,ks2
       do 110 i=is1(j),is2(j)
        ccdif = avqdz(i,j,k-1)
 
c- effet du flux phiz(k) sur q2(k) : aa(k)*q2(k) + bb(k)*q2(k-1)
        ccdt = cczdt(k) * tms(i,j,k-1)
        bb(i,k) = ccdt * (-ccdif )
        aa(i,k) = aa(i,k) + ccdt * (ccdif )
c- effet du flux phiz(k) sur q2(k-1) : aa(k-1)*q2(k-1) + cc(k-1)*q2(k)
	km2=max(k-2,ks1)
        ccdt = cczdt(k-1) * tms(i,j,km2)
        aa(i,k-1) = aa(i,k-1) + ccdt * (ccdif )
        cc(i,k-1) = ccdt * (-ccdif)
 110  continue

      do 120 i=is1(j),is2(j)
        ff(i,ks2) = ff(i,ks2) + tms(i,j,ks2)*cczdt(ks2)
     &                         *avqdz(i,j,ks2)*q2turb(i,j,ks2+1)
	aa(i,ks2) = aa(i,ks2) + tms(i,j,ks2)*cczdt(ks2)*avqdz(i,j,ks2)
	aa(i,ks1) = 1.0
 120  continue

 
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) decomposition LU - (methode et notation : cf Linear Algebra pp165-167)
c-----------------------------------------------------------------------
 
c--calcul de 1/alpha(k) dans aa(i,k) et beta(k) dans bb(i,k)
      do 200 i=is1(j),is2(j)
        aa(i,ks1) = 1.0 / aa(i,ks1)
 200  continue
      do 210 k=ks1+1,ks2
       do 210 i=is1(j),is2(j)
        bb(i,k) = bb(i,k) * aa(i,k-1)
        aa(i,k) = 1.0 / ( aa(i,k) - bb(i,k) * cc(i,k-1) )
 210  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) substitutions avant et arriere .                                |
c-----------------------------------------------------------------------
 
c--calcul de g(k) dans ff(i,k) :
      do 410 k=ks1+1,ks2
       do 410 i=is1(j),is2(j)
        ff(i,k) = ff(i,k) - bb(i,k) * ff(i,k-1)
 410  continue
 
c--calcul de x(k) dans scal(i,j,k,ns) :
      do 420 i=is1(j),is2(j)
       q2turb(i,j,ks2) = max(q2turb(i,j,ks2+1),ff(i,ks2) * aa(i,ks2))
 420  continue

      do 430 k=ks2-1,ks1,-1
       do 430 i=is1(j),is2(j)
        q2turb(i,j,k) =  max(q2tmin,(ff(i,k) - cc(i,k)
     &                  *q2turb(i,j,k+1)) * aa(i,k)*tms(i,j,k) )
 430  continue
 
c--Fin de la boucle externe sur l'indice de latitude j .
 500  continue


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Ajust.Conv. par augmentation de la Diffusion Verticale (scal)
c-----------------------------------------------------------------------
      if (lstab.eq.0) then
       do j=js1,js2
        do i=is1(j),is2(j)
         do k=ks1,ks2
          avsdz(i,j,k)=avsdz(i,j,k)+zaju(i,j,k)
         enddo
        enddo
       enddo
      endif
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine engtur -
      end
