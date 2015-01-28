












      subroutine isodiffu(scalat,ns)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c=============================================================
c  Adaptation of Cox small-slope Isopycnal Tensor
c  References:
c  -Cox, M.D, 1987, Isopycnal diffusion in a z-coordinate
c        model, Ocean Modelling, 74, 1-5
c  -Redi, M.H., 1982, Oceanic isopycnal mixing by coordinate
c        rotation, JPO, 12, 1154-1158
c=============================================================
c  modif : 25/05/99
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'isoslope.com'
      include 'comunit.h'
 
c- dummy variables :
      dimension scalat(imax,jmax,kmax)
 
c--variables locales conservees d'un appel a l'autre :
      common / difisloc / dtsdx(kmax), dtsdy(kmax), dtsdz(kmax)
 
c- local variables :
      dimension dscalx(imax,kmax), dscaly(jmax,kmax)
      dimension fisox(imax,jmax,kmax), fisoy(jmax,kmax)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      if (numit.eq.nstart .and. ns.eq.1) then
c---------------------------------
c ATTENTION : Modif de AHS !!!
c---------------------------------
        write(iuo+66,'(A,3(A,F6.0))') ' *** isodiffu ; Isopycn. *** :',
     &          ' ahs=', ahs(kmax), ' +', ai(kmax), '(=ai)'
        do 10 k=1,kmax
          ahs(k) = ahs(k)+ai(k)
          dtsdz(k) = unsdz(k)*dts(k)
          dtsdx(k) = unsdx*dts(k)
          dtsdy(k) = unsdy*dts(k)
 10     continue
 
c- fin du traitement 1ere iter.
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c ****** GRADIENTS OF SCALAR scal *******
c- NB1: for Gradient or Fluxes: in X => isf2(j)+1 & in Y => js2+1
c- NB2: isf1 is prefered to is1 to take into account
c        the Northern Margin (including Bering).
c- NB3: We take Minus flux * metric coefficient !!
c- NB4: BERING Raccord cmx&cmy coherent ONLY at point 2
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Start  external loop on all Meridional Section, index "j".
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,k,km1,km0,kp1,kp0,dscalx)
      do 300 j=js1,js2
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Merid. Section : Compute Isopycn.Flx in X & Z(part.1) Directions
c-----------------------------------------------------------------------
 
c- Vertical Gradient of scalar "ns" :
      do 210 k=ks1+1,ks2
        do 210 i=is1(j),is2(j)+1
          dscalz(i,j,k) = tms(i,j,k-1) * unsdzw(k) *
     &                  ( scal(i,j,k,ns) - scal(i,j,k-1,ns) )
 210  continue
 
c- Cyclic Conditions :
      if (ltest.ge.1) then
        do 215 k=ks1+1,ks2
          dscalz(1,j,k) = dscalz(ims2,j,k)
          dscalz(imax,j,k) = dscalz(ims1,j,k)
 215    continue
      endif
 
      do 230 k=ks1,ks2
        km1 = max(k-1, ks1)
        km0 = km1 + 1
        kp1 = min(k+1, ks2)
        kp0 = kp1 - 1
        do 230 i=is1(j),is2(j)+1
c- Grad.X of scalar "ns" :
          dscalx(i,k) = ttm1(i,j,k) * unsdx * smx(i,j,1) *
     &                ( scal(i,j,k,ns) - scal(i-1,j,k,ns) )
c- X.Flx of scalar "ns" :
          fisox(i,j,k) = cmy(i,j,1) * ai(k) * c1x(i,j,k) *
     &       ( dscalz(i,j,kp1) + dscalz(i-1,j,kp1)
     &       + dscalz(i,j,km0) + dscalz(i-1,j,km0) ) /
     &       ( tms(i,j,kp0) + tms(i-1,j,kp0)
     &       + tms(i,j,km1) + tms(i-1,j,km1) + epsil )
 230  continue
 
      do 250 k=ks1+1,ks2
c- Z.Flx (1rst part) of scalar "ns" :
        do 250 i=is1(j),is2(j)
          fisoz(i,j,k) = c4x(i,j,k) *
     &       ( dscalx(i,k) + dscalx(i+1,k-1)
     &       + dscalx(i,k-1) + dscalx(i+1,k) ) /
     &       ( ttm1(i,j,k) + ttm1(i+1,j,k-1)
     &       + ttm1(i,j,k-1) + ttm1(i+1,j,k) + epsil )
 250  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- End of external loop on all Meridional Section, index "j".
 300  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Zonal Section : Compute Isopycn.Flx in Y & Z(part.2) Directions
c-----------------------------------------------------------------------
 
c- Raccord a Bering :
      if (ltest.eq.3) then
        do 310 k=ks1+1,ks2
          dscalz(iberpm,jberp,k) = dscalz(ibera, jberam,k)
          dscalz(iberp, jberp,k) = dscalz(iberam,jberam,k)
          dscalz(iberam,jbera,k) = dscalz(iberp, jberpm,k)
          dscalz(ibera, jbera,k) = dscalz(iberpm,jberpm,k)
 310    continue
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Start  external loop on all Zonal Section, index "i".
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(j,k,km1,km0,kp1,kp0,dscaly,fisoy)
      do 500 i=ims1,ims2
 
      do 330 k=ks1,ks2
        km1 = max(k-1, ks1)
        km0 = km1 + 1
        kp1 = min(k+1, ks2)
        kp0 = kp1 - 1
        do 330 j=jsdom1(i),jsdom2(i)+1
c- Grad.Y of scalar "ns" :
          dscaly(j,k) = ttm2(i,j,k) * unsdy * smy(i,j,2) *
     &                ( scal(i,j,k,ns) - scal(i,j-1,k,ns) )
c- Y.Flx of scalar "ns" :
          fisoy(j,k) = cmx(i,j,2) * ai(k) * c2y(i,j,k) *
     &       ( dscalz(i,j,kp1) + dscalz(i,j-1,kp1)
     &       + dscalz(i,j,km0) + dscalz(i,j-1,km0) ) /
     &       ( tms(i,j,kp0) + tms(i,j-1,kp0)
     &       + tms(i,j,km1) + tms(i,j-1,km1) + epsil )
 330  continue
 
      do 350 k=ks1+1,ks2
        do 350 j=jsdom1(i),jsdom2(i)
c- Z.Flx (2nd  part) of scalar "ns" :
          fisoz(i,j,k) = ai(k) * ( fisoz(i,j,k) + c4y(i,j,k) *
     &       ( dscaly(j,k) + dscaly(j+1,k-1)
     &       + dscaly(j,k-1) + dscaly(j+1,k) ) /
     &       ( ttm2(i,j,k) + ttm2(i,j+1,k-1)
     &       + ttm2(i,j,k-1) + ttm2(i,j+1,k) + epsil ) )
 350  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Bilan des Flux Incorpore dans "scalat".
c-----------------------------------------------------------------------
 
c-*******************
c- BALANCE OF FLUXES
c- !! We take Minus the gradient of fluxes we took minus the flux
c-*******************
 
      do 450 k=ks1,ks2
       do 450 j=jsdom1(i),jsdom2(i)
        scalat(i,j,k) = scalat(i,j,k) + smxy(i,j,0) *
     &                     ( dtsdx(k) * (fisox(i,j,k)-fisox(i+1,j,k))
     &                     + dtsdy(k) * (fisoy(j,k)-fisoy(j+1,k))    )
     &                + dtsdz(k) * (fisoz(i,j,k)-fisoz(i,j,k+1))
 450  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- End of external loop on all Zonal Section, index "i".
 500  continue
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine isodiffu -
      end
