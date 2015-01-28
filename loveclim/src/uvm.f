












      subroutine uvm
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Resolution des termes verticaux de l'eq.qqmvt par appel a la routine "uvi".
c**AJUSTEMENT DE LA MOYENNE DES VITESSES HORIZONTALES**
c  modif : 08/01/96
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
 
c--variables locales :
      dimension ui(imax), vi(imax)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Calcul des moyennes de eta,ub,vb-spl , remplace eta, ub, vb :   |
c-----------------------------------------------------------------------
 
       j = js1
       do 110 i=is1(j),is2(j)
        eta(i,j) = etaspl(i,j)
 110  continue
 
c--Debut de la boucle externe sur l'indice de latitude j :
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,k,ui,vi)
      do 400 j=ju1,ju2
c-----
 
       do 120 i=is1(j),is2(j)
        eta(i,j) = etaspl(i,j)
 120   continue
 
       do 130 i=iu1(j),iu2(j)
        ub(i,j) = ubspl(i,j)
        vb(i,j) = vbspl(i,j)
 130   continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Resolution Implicite des termes verticaux de l'eq. de qqmvt .   |
c-----------------------------------------------------------------------
 
C     do 280 j=ju1,ju2
        call uvi(j)
 280  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Ajuste la vitesse Horizontale sur sa moyenne sur la Verticale . |
c-----------------------------------------------------------------------
 
      do 300 i=iu1(j),iu2(j)
       ui(i)=0.0
       vi(i)=0.0
 300  continue
 
      do 310 k=ku1,ku2
       do 310 i=iu1(j),iu2(j)
        ui(i) = ui(i) + dz(k)*u(i,j,k)
        vi(i) = vi(i) + dz(k)*v(i,j,k)
 310  continue
      do 320 i=iu1(j),iu2(j)
        ui(i) = unshu(i,j)*(ub(i,j)-ui(i))
        vi(i) = unshu(i,j)*(vb(i,j)-vi(i))
 320  continue
      do 330 k=ku1,ku2
       do 330 i=iu1(j),iu2(j)
        u(i,j,k) = tmu(i,j,k) * ( u(i,j,k) + ui(i) )
        v(i,j,k) = tmu(i,j,k) * ( v(i,j,k) + vi(i) )
 330  continue
 
c-----
c--Fin de la boucle externe sur l'indice de latitude j .
 400  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  6 ) raccord cyclique + raccord de grille pour eta, ub, vb, u, v .   |
c-----------------------------------------------------------------------
 
      if (ltest.eq.3) then
c--vitesse barocline nulle :
        vvber = unshu(iberp,jberp) * vb(iberp,jberp)
        do 650 k=ks1,ks2
          u(iberp,jberp,k) = 0.
          v(iberp,jberp,k) =  tmu(iberp,jberp,k) * vvber
 650    continue
      endif
 
      call raccord( eta, zero, 1, 0)
      call raccord( ub, zero, 1, 23)
      call raccord( vb, zero, 1, 35)
      call raccord( u(1,1,1), zero, kmax, 15)
      call raccord( v(1,1,1), zero, kmax, 27)
 
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  7 ) Computation of the average of the surface velocity              |
c-----------------------------------------------------------------------
 
      unsncl = one / DFLOAT(nclin)
      do j=ju1,ju2
        do i=iu1(j),iu2(j)
          umoy(i,j)=umoy(i,j)+u(i,j,ku2)*unsncl
          vmoy(i,j)=vmoy(i,j)+v(i,j,ku2)*unsncl
        enddo
      enddo
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine uvm -
      end
