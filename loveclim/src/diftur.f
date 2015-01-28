












      subroutine diftur
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c Computation of mixing length : vlturb
c Computation of the coefficients for vertical DIFFUsivities (divided by dz) :
c  avsdz, avudz,avqdz. 
c DIFFUsivities :  avxdz = sxturb * vlturb * vcturb (=sqrt(q2turb) )
c
c--
c  modif : 25/10/96

 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'ice.com'
      include 'comunit.h'

c--variables locales conservees d'un appel a l'autre :
      common / vdifloc / bvfmix,avu0ri(kmax), 
     &           avs0ri(kmax),avq0ri(kmax),
     &           ref0n2,ref0n
      common / kdifloc / kk0ri, kk1ri
Cdn2  common /coefin/ ci1,ci2
C     dimension avuloc(imax,jmax,kmax)
      dimension vcturb(imax,jmax,kmax)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) 1ere Iter, Determine 1er niv. a calculer ; Initialisation .     |
c-----------------------------------------------------------------------
      q2blmi = 1.d-6
      zrimax = 0.7d0
      q2slon = 6.51d0
      if (numit.le.nstart) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
Cdn2   ci1=1.0 - 1.0/2.0
Cdn2   ci2=1.0 - sqrt(2.0)/2.0
        write(iuo+66,'(2A,F6.2,A1,I4,A1)')
     &   ' diftur : mixing length from algebric formula '
 
c--1er Niveau a calculer :
      kk0ri = ks1
      avu0ri(ks1) = 0.0
      avs0ri(ks1) = 0.0
      avq0ri(ks1) = unsdz(ks1) * avkb(ks1+1)
      do 30 k=ks1+1,ks2
        avu0ri(k) = unsdzw(k) * avnub(k)
        avs0ri(k) = unsdzw(k) * avkb(k)
        avq0ri(k) = unsdz(k)  * avkb(k)
 30   continue

      kk1ri = kk0ri + 1
 
       do j=js1,js2
        do i=ims1,ims2
           q2turb(i,j,ks2+1) = max(q2tmin,q2slon*
     &                   ust2s(i,j)*(one+sdvt(i,j)))
        enddo
       enddo

c- Fin du traitement specifique de la 1ere Iter.
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
      epsil2 = epsil * epsil
c- valeur typique de N2 : 1.E-5 (s-2)
      ref0n2 = 1.d-5 * epsil
      ref0n  = sqrt(ref0n2)
 
c--Initialisation :
      do 70 k=ks1+1,kk0ri
       do 70 j=js1,js2
        do 70 i=ims1,ims2
         avudz(i,j,k) = avu0ri(k)
            avsdz(i,j,k) = avs0ri(k)
            avqdz(i,j,k) = avq0ri(k)
 70   continue

      do 80 k=kk0ri,ks2
       do 80 j=js1,js2
        do 80 i=ims1,ims2
         avqdz(i,j,k) = 0.0
 80   continue

 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Mixing Length :                                                 |
c-----------------------------------------------------------------------
c

      do 180 k=kk1ri,ks2
       do 180 j=js1,js2
        do 180 i=is1(j),is2(j)
	  vcturb(i,j,k) = sqrt ( q2turb(i,j,k) )
	  sqbvf  = sqrt (max(zero,bvf(i,j,k))) + ref0n

c Mixing length as an algebric fonction
c
	  ds     = zw(ks2+1) - zw(k)
	  db     = max (epsil, zw(k) - zw(kfs(i,j)))
	  zintr  = vkappa*(ds*db)/(ds+db)
          vlneut = zlotur*zintr/(zintr+zlotur)
	  vlturb(i,j,k) = min (vlneut, sqrghm*vcturb(i,j,k)/sqbvf) 
     &                  + vlmin
 180  continue

      do 190 j=js1,js2
       do 190 i=is1(j),is2(j)
        vlturb(i,j,ks2) = max(vlturb(i,j,ks2),-vkappa*zw(ks2))
 190  continue
	  
c--Debut de la boucle externe sur lindice de niveau k :
      do 500 k=kk1ri,ks2
c-----

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) M2turb                                                          |
c-----------------------------------------------------------------------
      do 400 j=js1,js2
        do 400 i=is1(j),is2(j)

          tm4u = tmu(i,j,k-1) + tmu(i+1,j+1,k-1)
     &         + tmu(i,j+1,k-1) + tmu(i+1,j,k-1) + epsil2

          uz = u(i+1,j+1,k) - u(i+1,j+1,k-1)
          vz = v(i+1,j+1,k) - v(i+1,j+1,k-1)
          riums1 = (uz*uz + vz*vz) * tmu(i+1,j+1,k-1)

          uz = u(i+1,j,k) - u(i+1,j,k-1)
          vz = v(i+1,j,k) - v(i+1,j,k-1)
          riums2 = (uz*uz + vz*vz) * tmu(i+1,j,k-1)

          uz = u(i,j+1,k) - u(i,j+1,k-1)
          vz = v(i,j+1,k) - v(i,j+1,k-1)
          riums3 = (uz*uz + vz*vz) * tmu(i,j+1,k-1)

          uz = u(i,j,k) - u(i,j,k-1)
          vz = v(i,j,k) - v(i,j,k-1)
          riums4 = (uz*uz + vz*vz) * tmu(i,j,k-1)

          tm2tur(i,j,k) = max ( riums1 , riums2 , riums3 , riums4 )
     &                   *max(zero,sign(one,tm4u-one))
     &                         * unsdzw(k) * unsdzw(k)
C    &                   * (1.0+max(zero,sign(one,dfloat(k-kajul))) )
Cvar &    *(1.0+varfor*sdvt(i,j)*max(zero,sign(one,dfloat(k-kajul))) )

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Stabilty functions .                                            |
c-----------------------------------------------------------------------

	  ghturb = - (vlturb(i,j,k)*vlturb(i,j,k))*bvf(i,j,k)
     &                     /max(q2turb(i,j,k),epsil2)
	  ghturb = max(ghmin,min(ghmax,ghturb))
	  sqturb = 0.2d0
	  ssturb = (0.494d0)/(1.-34.7d0*ghturb)
	  suturb = (0.393d0-3.09d0*ghturb)/
     &               (1.0 - 40.8d0*ghturb + 212.*ghturb*ghturb)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Diffusion coefficients .                                        |
c-----------------------------------------------------------------------

	 avcom = vlturb(i,j,k) * vcturb(i,j,k)
	 avuloc(i,j,k)= max (suturb * avcom * unsdzw(k), avu0ri(k))
            avsdz(i,j,k) = max (ssturb * avcom * unsdzw(k), avs0ri(k))

c background diffusivity as a function of N2

c Computation of an averaged N2
Cdn2      zaxmoy = max(epsil,tms(i,j,k)
Cdn2 &             +ci1*(tms(i-1,j,k)+tms(i,j-1,k)
Cdn2 &                  +tms(i+1,j,k)+tms(i,j+1,k))
Cdn2 &             +c12*(tms(i-1,j-1,k)+tms(i-1,j+1,k)
Cdn2 &                  +tms(i+1,j-1,k)+tms(i+1,j+1,k))  )
Cdn2      bvfmoy = (tms(i,j,k)*bvf(i,j,k)
Cdn2 &             +ci1*(tms(i-1,j,k)*bvf(i-1,j,k)
Cdn2 &                  +tms(i,j-1,k)*bvf(i,j-1,k)
Cdn2 &                  +tms(i+1,j,k)*bvf(i+1,j,k)
Cdn2 &                  +tms(i,j+1,k)*bvf(i,j+1,k))
Cdn2 &             +c12*(tms(i-1,j-1,k)*bvf(i-1,j-1,k)
Cdn2 &                  +tms(i-1,j+1,k)*bvf(i-1,j+1,k)
Cdn2 &                  +tms(i+1,j-1,k)*bvf(i+1,j-1,k)
Cdn2 &                  +tms(i+1,j+1,k)*bvf(i+1,j+1,k))  )
Cdn2 &             /zaxmoy
Cdn2      avsmir = max((1.D-7/sqrt(max(4.D-6*one,bvfmoy)) )
Cdn2 &               * unsdzw(k), avs0ri(k) )
Cdn2      avsdz(i,j,k) = max (ssturb * avcom * unsdzw(k), avsmir)

          
c !! location of avqdz(k) : idem scal (k)

	  zintop = (zw(k+1)-z(k))*unsdz(k)
	  zinbot = (z(k-1)-zw(k-1))*unsdz(k-1)
          avqdz(i,j,k) = avqdz(i,j,k) +
     &    zintop * max (sqturb * avcom *unsdz(k), avq0ri(k))
          avqdz(i,j,k-1) = avqdz(i,j,k-1) +
     &    zinbot * max (sqturb * avcom *unsdz(k-1), avq0ri(k))

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  6 ) Diffusion in the highly sheared region below the ML             |
c-----------------------------------------------------------------------

C         ztest = abs (-1.*max(zero,sign(one,q2turb(i,j,k)-q2blmi)) 
C    &         + one-max(zero,sign(one,q2turb(i,j,k+1)-q2blmi)))
C         zri   = bvf(i,j,k)/max(tm2tur(i,j,k),epsil)
C         zri   = max(zero,min(zri,zrimax) )
C         zkhri = (5.0d-3*(1-(zri/zrimax)*(zri/zrimax) )**3 )*unsdzw(k)
C           
C         avuloc(i,j,k)= avuloc(i,j,k)+ (1.0-ztest)*zkhri*tms(i,j,k-1)
C         avsdz(i,j,k) = avsdz(i,j,k) + (1.0-ztest)*zkhri*tms(i,j,k-1)

 400   continue

c-raccord cyclique pour avuloc
       do j=1,jeq
          avuloc(1,j,k)=avuloc(imax-1,j,k)
          avuloc(imax,j,k)=avuloc(2,j,k)
       enddo
       
       do 410 j=ju1,ju2
        do 410 i=iu1(j),iu2(j)
C        avudz(i,j,k)=.25* ( avuloc(i,j,k)+avuloc(i-1,j,k)+
C    &                       avuloc(i,j-1,k)+avuloc(i-1,j-1,k) )
         avudz(i,j,k)=( avuloc(i,j,k)*tms(i,j,k)
     &                 +avuloc(i-1,j,k)*tms(i-1,j,k)
     &                 +avuloc(i,j-1,k)*tms(i,j-1,k)
     &                 +avuloc(i-1,j-1,k)*tms(i-1,j-1,k) )/
     &               max(tms(i,j,k)+tms(i-1,j,k)+tms(i,j-1,k)
     &                  +tms(i-1,j-1,k),q2blmi)

 410   continue
       

 500  continue
 

      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine diftur -
      end
