












      subroutine start
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  demarrage a partir de l'etat de repos :
c  modif : 24/09/99
 
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'
 
c--variables locales :
 
      write(iuo+66,*) 'debut de start'
 
      tpstot = 0.0
      numit = 0
 
      do 200 n=1,nsmax
       do 200 k=ks1,ks2
        do 200 j=1,jmax
         do 200 i=1,imax
          scal(i,j,k,n) = scal0(k,n)
 200  continue
 
      do 210 k=ks1,ks2
       do 210 j=1,jmax
        do 210 i=1,imax
         u(i,j,k) = 0.0
         v(i,j,k) = 0.0
 210  continue
 
      do 220 j=1,jmax
       do 220 i=1,imax
        eta(i,j) = 0.0
        ub(i,j) = 0.0
        vb(i,j) = 0.0
 220  continue
 
      do 240 j=1,jmax
        do 230 i=1,imax
c
c                        tfu: MELTING POINT OF SEA WATER.
c
            tfu(i,j)    = abs(273.15-0.0575*scal(i,j,ks2,2)+
     &                        1.710523e-03*sqrt(scal(i,j,ks2,2))**3-
     &                        2.154996e-04*scal(i,j,ks2,2)**2)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
            hgbq(i,j)   = 0.0
            albq(i,j)   = 1.0
            hnbq(i,j)   = 0.0
            ts(i,j)     = tfu(i,j)
            tbq(i,j,1)  = tfu(i,j)
            tbq(i,j,2)  = tfu(i,j)
            tbq(i,j,3)  = tfu(i,j)
            firg(i,j)   = 0.0
            fcsg(i,j)   = 0.0
            fleg(i,j)   = 0.0
            fsbbq(i,j)  = 0.0
            qstobq(i,j) = 0.0
            xzo(i,j)    = 0.001
 230    continue
 240  continue
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine start -
      end
