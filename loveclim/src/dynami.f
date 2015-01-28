












      subroutine dynami_zh(ih)
c
c determine the velocity field of sea ice
c
c forcing: wind stress, water stress,
c           surface tilt
c
c ice-ice interaction: non-linear viscous-
c                       plastic law and bulk 
c                       rheology
c
c method: Zhang and Hibler, 1997
c
c ih: hemispheric index (=1, NH; =2, SH) 
c iu1: western most demarcation for velocity
c      points along a zonal belt (iu1.ge.2)
c iu2: eastern most demarcation for velocity
c      points along a zonal belt (iu2.le.imax-1)
c
c M.A. Morales Maqueda ?-?-?
c modified M.A. Morales Maqueda 31-12-1993
c modified H. Goosse 08-12-1994
c modified M.A. Morales Maqueda 09-02-2000 
c
c Ccpl [Ccp0] => line pertaining to the
c                coupled (uncoupled) version
c Ccem        => concentric ellipse method 
c                Hibler, 1979
c Ctem        => truncated ellipse method 
c                Geiger et al., 1998
c Clds        => changes in ice concentration
c                associated with shearing
c
      include 'type.com'
      include 'para.com'
      include 'const.com'
      include 'bloc.com'
      include 'ice.com'
      include 'dynami.com'
c
      parameter (ntmax=1000)
c
      dimension presm(imax,jmax),presh(imax,jmax),
     &          gphs1(imax,jmax),gphs2(imax,jmax)
      dimension albqd(imax,jmax)
      dimension zmass(imax,jmax),zcorl(imax,jmax)
      dimension usdtzm(imax,jmax),
     &          uiner(imax,jmax),viner(imax,jmax),
     &          uzcorl(imax,jmax),vzcorl(imax,jmax)
      dimension tagnx(imax,jmax),tagny(imax,jmax)
      dimension u0(imax,jmax),v0(imax,jmax)
      dimension viszeta(imax,jmax),viseta(imax,jmax),
     &          sigm11(imax,jmax),
     &          sigm12(imax,jmax),
     &          sigm22(imax,jmax)
      dimension c1m(imax,jmax),c10(imax,jmax),
     &          c1p(imax,jmax),
     &          r1(imax,jmax),dvsg1(imax,jmax),
     &          c2m(imax,jmax),c20(imax,jmax),
     &          c2p(imax,jmax),
     &          r2(imax,jmax),dvsg2(imax,jmax),
     &          zaw(imax,jmax),zbw(imax,jmax)
      dimension gamy(imax,jmax),bety(ntmax),
     &          gamx(ntmax),betx(ntmax),uu(ntmax)
      dimension up(imax,jmax),vp(imax,jmax),
     &          uaux(imax,jmax),vaux(imax,jmax),
     &          res(imax,jmax),uvmsk(imax,jmax)
      dimension zm_sc(imax,jmax),
     &          za1ct(imax,jmax),za2ct(imax,jmax),
     &          zb1(imax,jmax),zb2(imax,jmax)
      dimension wstrn(imax,jmax)
c
c define limits where to compute dynamics.
c
      if (ih.eq.1) then
        njmin = nlminn
        njmax = nlmaxn
      else
        njmin = nlmins
        njmax = nlmaxs
      endif
      njmm1  = njmin-1
      imaxm1 = imax-1
c
c initialise certain arrays. this could
c be done only once per run by introducing
c a flag and placing the arrays in a common
c
      do j=njmm1,njmax+1
        do i=1,imax
          uvmsk(i,j) = zero
          up(i,j)    = zero
          vp(i,j)    = zero
          zaw(i,j)   = one
          c10(i,j)   = zero
          r1(i,j)    = zero
          dvsg1(i,j) = zero
          c20(i,j)   = zero
          r2(i,j)    = zero
          dvsg2(i,j) = zero
        enddo
      enddo
      uaux(:,:) = 0.0
      vaux(:,:) = 0.0
c
c sign of turning angle for oceanic drag
c
      sang = real(ih)*sangvg
c
c ice mass, ice strength, and wind stress 
c at the center of the grid cells                                                
c
      do j=njmm1,njmax
        do i=1,imax
          zm_sc(i,j) = tms(i,j,ks2)
     &                 *(rhon*hnm(i,j)+rhog*hgm(i,j))
          wstrn(i,j) = exp(-c*albq(i,j))
          presm(i,j) = tms(i,j,ks2)
     &                 *pstarh*hgm(i,j)*wstrn(i,j)
Ccp1     if (icoupl .eq. 0) then
          zb1(i,j)   = tms(i,j,ks2)*(one-albq(i,j))
          zb2(i,j)   = tms(i,j,ks2)*(one-albq(i,j))
Ccp1     else
Ccp1      zb1(i,j)   = tms(i,j,ks2)
Ccp1                   *tenagx(i,j)*(one-albq(i,j))
Ccp1      zb2(i,j)   = tms(i,j,ks2)
Ccp1                   *tenagy(i,j)*(one-albq(i,j))
Ccp1     endif
c
        enddo
      enddo
      
c
c lead area, mass, coriolis coefficients
c and wind stress at the corners of the 
c grid cells + momentum terms that are
c independent of the velocity field
c
      ccgdx = gpes*uns2dx
      ccgdy = gpes*uns2dy
      do 20 j=njmin,njmax
        jm1 = j-1
        do 10 i=iu1(j),iu2(j)
          im1        = i-1
c
          usw        = tmu(i,j,ku2)/
     &                 max( tms(i,j,ks2)*wght(i,j,2,2)
     &                     +tms(im1,j,ks2)*wght(i,j,1,2)
     &                     +tms(i,jm1,ks2)*wght(i,j,2,1)
     &                     +tms(im1,jm1,ks2)*wght(i,j,1,1),
     &                     zepsd2)
c
c leads area
c
         albqd(i,j)  = ( tms(i,j,ks2)*albq(i,j)*wght(i,j,2,2)
     &                  +tms(im1,j,ks2)*albq(im1,j)*wght(i,j,1,2)
     &                  +tms(i,jm1,ks2)*albq(i,jm1)*wght(i,j,2,1)
     &                  +tms(im1,jm1,ks2)*albq(im1,jm1)*wght(i,j,1,1))
     &                 *usw
c
c mass
c
          zmass(i,j) = ( zm_sc(i,j)*wght(i,j,2,2)
     &                  +zm_sc(im1,j)*wght(i,j,1,2)
     &                  +zm_sc(i,jm1)*wght(i,j,2,1)
     &                  +zm_sc(im1,jm1)*wght(i,j,1,1))
     &                 *usw
          uvmsk(i,j) = one-max(zero,sign(one,one-zmass(i,j)))
          zcorl(i,j) = zmass(i,j)*zfn(i,j)
c
c wind stress
c
Ccp1     if (icoupl.eq.0) then
          tagnx(i,j) = (zb1(i,j)*wght(i,j,2,2)
     &                 +zb1(im1,j)*wght(i,j,1,2)
     &                 +zb1(i,jm1)*wght(i,j,2,1)
     &                 +zb1(im1,jm1)*wght(i,j,1,1))
     &                 *usw
     &                 *tenagx(i,j)
C
          tagny(i,j) = (zb2(i,j)*wght(i,j,2,2)
     &                 +zb2(im1,j)*wght(i,j,1,2)
     &                 +zb2(i,jm1)*wght(i,j,2,1)
     &                 +zb2(im1,jm1)*wght(i,j,1,1))
     &                 *usw
     &                 *tenagy(i,j)
Ccp1     else
Ccp1      tagnx(i,j) = (zb1(i,j)*wght(i,j,2,2)
Ccp1 &                 +zb1(im1,j)*wght(i,j,1,2)
Ccp1 &                 +zb1(i,jm1)*wght(i,j,2,1)
Ccp1 &                 +zb1(im1,jm1)*wght(i,j,1,1))
Ccp1 &                 *usw
C
Ccp1      tagny(i,j) = (zb2(i,j)*wght(i,j,2,2)
Ccp1 &                 +zb2(im1,j)*wght(i,j,1,2)
Ccp1 &                 +zb2(i,jm1)*wght(i,j,2,1)
Ccp1 &                 +zb2(im1,jm1)*wght(i,j,1,1))
Ccp1 &                 *usw
Ccp1     endif
c
c terms that are independent of the velocity field
c (surface tilt + wind stress)
c
         detadxi     = zmass(i,j)*ccgdx*smx(i,j,3)
     &                 *( (eta(im1,jm1)-eta(i,j))
     &                   +(eta(im1,j)-eta(i,jm1)))
         detadyi     = zmass(i,j)*ccgdy*smy(i,j,3)
     &                 *( (eta(im1,jm1)-eta(i,j))
     &                   -(eta(im1,j)-eta(i,jm1)))
         za1ct(i,j)  = tagnx(i,j)+detadxi
         za2ct(i,j)  = tagny(i,j)+detadyi
10      continue
20    continue
c
c cyclicity of velocity mask
c
      do j=njmin,njmax
        uvmsk(1,j)    = uvmsk(imaxm1,j)
        uvmsk(imax,j) = uvmsk(2,j)
      enddo
c
c initial conditions
c
      do j=njmin,njmax
        do i=1,imax
          ug(i,j) = uvmsk(i,j)*ug(i,j)
          vg(i,j) = uvmsk(i,j)*vg(i,j)
          u0(i,j) = ug(i,j)
          v0(i,j) = vg(i,j)
        enddo
      enddo
c
      if (ih.eq.1) then
c
c case of bering (compute velocity straight
c                 away!)
c
        zmabe   = ( rhon*hnm(iberp,jberp-1)
     &             +rhog*hgm(iberp,jberp-1)
     &             +rhon*hnm(iberp-1,jberp-1)
     &             +rhog*hgm(iberp-1,jberp-1)
     &             +rhon*hnm(ibera-1,jbera-1)
     &             +rhog*hgm(ibera-1,jbera-1)
     &             +rhon*hnm(ibera,jbera-1)
     &             +rhog*hgm(ibera,jbera-1))/4.0d0
        etapac  =  eta(iberp-1,jberp-3)+eta(iberp,jberp-3)
     &            +eta(iberp-1,jberp-2)+eta(iberp,jberp-2)
        etaarc  =  eta(ibera,jbera-3)+eta(ibera+1,jbera-3)
     &            +eta(ibera,jbera-2)+eta(ibera+1,jbera-2)
        detadyi =  zmabe*ccgdy*smy(iberp,jberp-1,3)
     &            *(etapac-etaarc)*0.25d0
        tabe    =  tenagy(iberp,jberp-1)+detadyi
        pbe     = ( presm(iberp,jberp-1)
     &             +presm(iberp-1,jberp-1)
     &             +presm(ibera,jbera-1)
     &             +presm(ibera-1,jbera-1))/4.0d0
        pslbe   =  pbe/(4.0d0*50.0d03)
        u1      =  max(zero,vo(iberp,jberp)
     &            +sign(one,tabe-pslbe)
     &             *sqrt(abs(tabe-pslbe)/rhoco))
        u2      =  min(zero,vo(iberp,jberp)
     &            +sign(one,tabe+pslbe)
     &             *sqrt(abs(tabe+pslbe)/rhoco))
        ubering =  (u1+u2)*50.0d0/125.0d0*sber
      endif
c
c solution of the momentum equation taking into
c account ice-ice interactions
c
      do 3000 iter=1,2*nbiter
c
        zindu = 0.5d0*mod(iter,2)
        usdtp = (one+2.0d0*zindu)*nbiter*usdt
c
c computation of free drift field for free slip
c boundary conditions
c
c       if (bound.gt.0.5d0) 
c    &    call frdrft(ih,usdtp,u0,v0,zmass,zcorl,tagnx,tagny)
c
c computation of viscosities
c
        do 40 j=njmm1,njmax
          jp1 = j+1
          do 30 i=1,imax
            ip1 = (i+1)-(imax-2)*(i/imax)
c
c rate of strain tensor
c
            ze11         = akappa(i,j,1,1)*
     &                      (ug(ip1,j)+ug(ip1,jp1)-(ug(i,j)+ug(i,jp1)))
     &                    +akappa(i,j,1,2)*
     &                      (vg(ip1,j)+vg(ip1,jp1)+vg(i,j)+vg(i,jp1))
            ze12         = akappa(i,j,2,2)*
     &                      (ug(i,jp1)+ug(ip1,jp1)-(ug(i,j)+ug(ip1,j)))
     &                    -akappa(i,j,2,1)*
     &                      (vg(i,j)+vg(ip1,j)+vg(i,jp1)+vg(ip1,jp1))
            ze22         = akappa(i,j,2,2)*
     &                      (vg(i,jp1)+vg(ip1,jp1)-(vg(i,j)+vg(ip1,j)))
     &                    +akappa(i,j,2,1)*
     &                      (ug(i,j)+ug(ip1,j)+ug(i,jp1)+ug(ip1,jp1))
            ze21         = akappa(i,j,1,1)*
     &                      (vg(ip1,j)+vg(ip1,jp1)-(vg(i,j)+vg(i,jp1)))
     &                    -akappa(i,j,1,2)*
     &                      (ug(ip1,j)+ug(ip1,jp1)+ug(i,j)+ug(i,jp1))
            trace        = ze11+ze22
            trace2       = trace*trace
            deter        = ze11*ze22-0.25d0*(ze12+ze21)**2
            zes2         = max(zero,trace2-4.0d0*deter)
            delta        = max(sqrt(trace2+zes2*usecc2),creepl)
            viszeta(i,j) = max(presm(i,j)/delta,zetamn)
Ccem        viseta(i,j)  = viszeta(i,j)*usecc2
Ccem        presh(i,j)   = presm(i,j)
            presh(i,j)   = delta*viszeta(i,j)
            zes          = max(zepsd2,sqrt(zes2))
            viseta1      = (presh(i,j)-viszeta(i,j)*trace)/zes
            viseta(i,j)  = min(viszeta(i,j)*usecc2,viseta1)
c
c Determination of stress tensor
c
            aa           = viszeta(i,j)*(ze11+ze22)
            sigm11(i,j)  = viseta(i,j)*(ze11-ze22)+aa
            sigm12(i,j)  = viseta(i,j)*(ze12+ze21)
            sigm22(i,j)  = viseta(i,j)*(ze22-ze11)+aa
30        continue
40      continue
c
c gradient of ice strength and some auxiliary arrays
c
        do j=njmin,njmax
          jm1 = j-1
          do i=iu1(j),iu2(j)
            im1        =  i-1
            gphs1(i,j) =  (alambd(i,j,2,2,2,1)-alambd(i,j,2,1,2,1))
     &                    *presh(i,jm1)
     &                   +(alambd(i,j,2,2,2,2)-alambd(i,j,2,1,2,2))
     &                    *presh(i,j)
     &                   -(alambd(i,j,2,2,1,1)+alambd(i,j,2,1,1,1))
     &                    *presh(im1,jm1)
     &                   -(alambd(i,j,2,2,1,2)+alambd(i,j,2,1,1,2))
     &                    *presh(im1,j)
            gphs2(i,j) = -(alambd(i,j,1,1,2,1)+alambd(i,j,1,2,2,1))
     &                    *presh(i,jm1)
     &                   -(alambd(i,j,1,1,1,1)+alambd(i,j,1,2,1,1))
     &                    *presh(im1,jm1)
     &                   +(alambd(i,j,1,1,2,2)-alambd(i,j,1,2,2,2))
     &                    *presh(i,j)
     &                   +(alambd(i,j,1,1,1,2)-alambd(i,j,1,2,1,2))
     &                    *presh(im1,j)
c
          usdtzm(i,j)  = usdtp*zmass(i,j)
          uiner(i,j)   = usdtzm(i,j)*u0(i,j)
          viner(i,j)   = usdtzm(i,j)*v0(i,j)
          uzcorl(i,j)  = ug(i,j)*zcorl(i,j)
          vzcorl(i,j)  = vg(i,j)*zcorl(i,j)
          enddo           
        enddo           

c
        do j=njmin,njmax
          jm1 = j-1
          do i=1,imax
            im1      = (i-1)+(imax-2)*(1/i)
c
c determination of c1m
c
            ze11     = -akappa(im1,jm1,1,1)
            sigm1111 = (viszeta(im1,jm1)+viseta(im1,jm1))*ze11
            sigm2211 = (viszeta(im1,jm1)-viseta(im1,jm1))*ze11
            ze11     = -akappa(im1,j,1,1)
            sigm1112 = (viszeta(im1,j)+viseta(im1,j))*ze11
            sigm2212 = (viszeta(im1,j)-viseta(im1,j))*ze11
            c1m(i,j) =  alambd(i,j,2,2,1,1)*sigm1111
     &                 +alambd(i,j,2,2,1,2)*sigm1112
     &                 +alambd(i,j,2,1,1,1)*sigm2211
     &                 +alambd(i,j,2,1,1,2)*sigm2212
c
c determination of c2m
c
            ze22     = -akappa(im1,jm1,2,2)
            sigm1111 = (viszeta(im1,jm1)-viseta(im1,jm1))*ze22
            sigm2211 = (viszeta(im1,jm1)+viseta(im1,jm1))*ze22
            ze22     = -akappa(i,jm1,2,2)
            sigm1121 = (viszeta(i,jm1)-viseta(i,jm1))*ze22
            sigm2221 = (viszeta(i,jm1)+viseta(i,jm1))*ze22
            c2m(i,j) =  alambd(i,j,1,1,2,1)*sigm2221
     &                 +alambd(i,j,1,1,1,1)*sigm2211
     &                 +alambd(i,j,1,2,1,1)*sigm1111
     &                 +alambd(i,j,1,2,2,1)*sigm1121
          enddo
        enddo
c
c determination of c1p and c2p
c
        do j=njmin,njmax
          jp1 = j+1
          do i=1,imax
            ip1      = (i+1)-(imax-2)*(i/imax)
            c1p(i,j) = c1m(ip1,j)
            c2p(i,j) = c2m(i,jp1)
          enddo
        enddo
c
c special boundary case, c2p(i,njmax)
c
        j = njmax
        do i=1,imax
          im1      = (i-1)+(imax-2)*(1/i)
          ze22     = akappa(im1,j,2,2)
          sigm1112 = (viszeta(im1,j)-viseta(im1,j))*ze22
          sigm2212 = (viszeta(im1,j)+viseta(im1,j))*ze22
          ze22     = akappa(i,j,2,2)
          sigm1122 = (viszeta(i,j)-viseta(i,j))*ze22
          sigm2222 = (viszeta(i,j)+viseta(i,j))*ze22
          c2p(i,j) = -alambd(i,j,1,1,2,2)*sigm2222
     &               -alambd(i,j,1,1,1,2)*sigm2212
     &               +alambd(i,j,1,2,1,2)*sigm1112
     &               +alambd(i,j,1,2,2,2)*sigm1122 
        enddo
c
        do 60 j=njmin,njmax
          jm1 = j-1
          do 50 i=iu1(j),iu2(j)
            im1 = i-1
c
c determination of c10
c
            ze12     =  akappa(im1,jm1,2,2)
            sigm1211 =  viseta(im1,jm1)*ze12
            ze12     =  akappa(i,jm1,2,2)
            sigm1221 =  viseta(i,jm1)*ze12
            ze12     = -akappa(im1,j,2,2)
            sigm1212 =  viseta(im1,j)*ze12
            ze12     = -akappa(i,j,2,2)
            sigm1222 =  viseta(i,j)*ze12
c
            c10(i,j) =  (alambd(i,j,1,1,1,1)-alambd(i,j,1,2,1,1))
     &                  *sigm1211
     &                 +(alambd(i,j,1,1,2,1)-alambd(i,j,1,2,2,1))
     &                  *sigm1221
     &                 -(alambd(i,j,1,1,2,2)+alambd(i,j,1,2,2,2))
     &                  *sigm1222
     &                 -(alambd(i,j,1,1,1,2)+alambd(i,j,1,2,1,2))
     &                  *sigm1212
     &                 -(c1m(i,j)+c1p(i,j))
c
c determination of c20
c 
            ze21     =  akappa(im1,jm1,1,1)
            sigm1211 =  viseta(im1,jm1)*ze21
            ze21     = -akappa(i,jm1,1,1)
            sigm1221 =  viseta(i,jm1)*ze21
            ze21     =  akappa(im1,j,1,1)
            sigm1212 =  viseta(im1,j)*ze21
            ze21     = -akappa(i,j,1,1)
            sigm1222 =  viseta(i,j)*ze21
c
            c20(i,j) = -(alambd(i,j,2,1,2,1)+alambd(i,j,2,2,2,1))
     &                  *sigm1221
     &                 -(alambd(i,j,2,1,2,2)+alambd(i,j,2,2,2,2))
     &                  *sigm1222
     &                 +(alambd(i,j,2,2,1,1)-alambd(i,j,2,1,1,1))
     &                  *sigm1211
     &                 +(alambd(i,j,2,2,1,2)-alambd(i,j,2,1,1,2))
     &                  *sigm1212
     &                 -(c2m(i,j)+c2p(i,j))
c
50        continue
60      continue
c
c mask c1m, c1p, c2m, and c2p
c
        do j=njmin,njmax
          do i=1,imax
            c1m(i,j) = c1m(i,j)*uvmsk(i,j)
            c1p(i,j) = c1p(i,j)*uvmsk(i,j)
            c2m(i,j) = c2m(i,j)*uvmsk(i,j)
            c2p(i,j) = c2p(i,j)*uvmsk(i,j)
          enddo
        enddo
c
c relaxation
c
        do j=njmm1,njmax
          do i=1,imax
            uaux(i,j) = ug(i,j)
            vaux(i,j) = vg(i,j)
          enddo
        enddo
c
        do 1000 jter=1,nbitdr 
c
          do j=njmin,njmax
            jm1 = j-1
            jp1 = j+1
            do i=iu1(j),iu2(j)
              im1        = i-1
              ip1        = i+1 
              up(i,j)    = uaux(i,j)
              vp(i,j)    = vaux(i,j)
c
c determination of right-hand side term of momentum equation
c (it is not possible to linearise the ice-ocean drag term,
c as doing so leads to important inaccuracies in the momentum
c balance when running the routine without pseudo-timesteps...)
c
              ur         = uaux(i,j)-uo(i,j)
              vr         = vaux(i,j)-vo(i,j)
              zmod       = sqrt(ur*ur+vr*vr)*(one-albqd(i,j))
              zcd        = rhoco*zmod
              zaw(i,j)   = one+uvmsk(i,j)*(usdtzm(i,j)+zcd*cangvg-one)
              r1(i,j)    = ( uiner(i,j)
     &                      +vzcorl(i,j)
     &                      +za1ct(i,j)
     &                      +zcd*( cangvg*uo(i,j)
     &                            -sang*(vo(i,j)-vaux(i,j)))
     &                      -gphs1(i,j))
     &                     *uvmsk(i,j)
              r2(i,j)    = ( viner(i,j)
     &                      -uzcorl(i,j)
     &                      +za2ct(i,j)
     &                      +zcd*( cangvg*vo(i,j)
     &                            +sang*(uo(i,j)-uaux(i,j)))
     &                      -gphs2(i,j))
     &                     *uvmsk(i,j)
c
c determination of divergence of stress tensor
c
              dvsg1(i,j) =  alambd(i,j,2,2,2,1)*sigm11(i,jm1)
     &                     +alambd(i,j,2,2,2,2)*sigm11(i,j)
     &                     -alambd(i,j,2,2,1,1)*sigm11(im1,jm1)
     &                     -alambd(i,j,2,2,1,2)*sigm11(im1,j)
     &                     -alambd(i,j,2,1,1,1)*sigm22(im1,jm1)
     &                     -alambd(i,j,2,1,2,1)*sigm22(i,jm1)
     &                     -alambd(i,j,2,1,1,2)*sigm22(im1,j)
     &                     -alambd(i,j,2,1,2,2)*sigm22(i,j)
     &                     +( alambd(i,j,1,2,2,1)
     &                       -alambd(i,j,1,1,2,1))*sigm12(i,jm1)
     &                     +( alambd(i,j,1,2,1,1)
     &                       -alambd(i,j,1,1,1,1))*sigm12(im1,jm1)
     &                     +( alambd(i,j,1,1,2,2)
     &                       +alambd(i,j,1,2,2,2))*sigm12(i,j)
     &                     +( alambd(i,j,1,1,1,2)
     &                       +alambd(i,j,1,2,1,2))*sigm12(im1,j)
              dvsg1(i,j) = ( dvsg1(i,j)
     &                      +c1m(i,j)*uaux(im1,j)
     &                      +c10(i,j)*uaux(i,j)
     &                      +c1p(i,j)*uaux(ip1,j))
     &                     *uvmsk(i,j)
c
              dvsg2(i,j) = -alambd(i,j,1,2,1,1)*sigm11(im1,jm1)
     &                     -alambd(i,j,1,2,2,1)*sigm11(i,jm1)
     &                     -alambd(i,j,1,2,1,2)*sigm11(im1,j)
     &                     -alambd(i,j,1,2,2,2)*sigm11(i,j)
     &                     -alambd(i,j,1,1,2,1)*sigm22(i,jm1)
     &                     -alambd(i,j,1,1,1,1)*sigm22(im1,jm1)
     &                     +alambd(i,j,1,1,2,2)*sigm22(i,j)
     &                     +alambd(i,j,1,1,1,2)*sigm22(im1,j)
     &                     +( alambd(i,j,2,1,1,1)
     &                       -alambd(i,j,2,2,1,1))*sigm12(im1,jm1)
     &                     +( alambd(i,j,2,1,2,1)
     &                       +alambd(i,j,2,2,2,1))*sigm12(i,jm1)
     &                     +( alambd(i,j,2,1,1,2)
     &                       -alambd(i,j,2,2,1,2))*sigm12(im1,j)
     &                     +( alambd(i,j,2,1,2,2)
     &                       +alambd(i,j,2,2,2,2))*sigm12(i,j)
              dvsg2(i,j) = ( dvsg2(i,j)
     &                      +c2m(i,j)*vaux(i,jm1)
     &                      +c20(i,j)*vaux(i,j)
     &                      +c2p(i,j)*vaux(i,jp1))
     &                     *uvmsk(i,j)
            enddo
          enddo
c
          if (ih.eq.1) then
c
c case of bering (ice drift at bering
c strait is used as boundary condition)
c
            i          = ibera
            j          = jbera-1
            dvsg2(i,j) = ( dvsg2(i,j)
     &                    -c2p(i,j)*vaux(i,j+1))
     &                   *uvmsk(i,j)
            i          = iberp
            j          = jberp-1
            dvsg2(i,j) = ( dvsg2(i,j)
     &                    -c2p(i,j)*vaux(i,j+1))
     &                   *uvmsk(i,j)
          endif
c
c solve tridiagonal system in the x-direction.
c
          do 70 j=njmin,njmax
            is         = iu1(j)
            ie         = iu2(j)
c
c cyclic/non-cyclic conditions will be applied 
c depending on whether, at a given latitude, the 
c domain spans a full zonal belt or not. 
c
            if (ie-is.lt.imax-3) then
c
c non-cyclic tridiagonal case
c
              bet        = c10(is,j)+zaw(is,j)
              uaux(is,j) = (r1(is,j)+dvsg1(is,j))/bet
              do i=is+1,ie
                gamx(i)   = c1p(i-1,j)/bet
                bet       = c10(i,j)+zaw(i,j)-c1m(i,j)*gamx(i)
                uaux(i,j) = ( r1(i,j)+dvsg1(i,j)
     &                       -c1m(i,j)*uaux(i-1,j))/bet
              enddo
              do i=ie-1,is,-1
                uaux(i,j) = uaux(i,j)-gamx(i+1)*uaux(i+1,j)
              enddo
            else
c
c cyclic tridiagonal case
c
c first tridiagonal step 
c
              alp        = c1p(ie,j)
              bet        = c1m(is,j)
              gamma      = -(c10(is,j)+zaw(is,j)) 
              betx(is)   = -(gamma+gamma)
              uaux(is,j) = (r1(is,j)+dvsg1(is,j))/betx(is)
              do i=is+1,ie-1
                gamx(i)   = c1p(i-1,j)/betx(i-1)
                betx(i)   = c10(i,j)+zaw(i,j)-c1m(i,j)*gamx(i)
                uaux(i,j) = ( r1(i,j)+dvsg1(i,j)
     &                       -c1m(i,j)*uaux(i-1,j))/betx(i)
              enddo
              gamx(ie)   = c1p(ie-1,j)/betx(ie-1)
              betx(ie)   = c10(ie,j)+zaw(ie,j)-alp*bet/gamma
     &                     -c1m(ie,j)*gamx(ie)
              uaux(ie,j) = ( r1(ie,j)+dvsg1(ie,j)
     &                      -c1m(ie,j)*uaux(ie-1,j))/betx(ie)
              do i=ie-1,is,-1
                uaux(i,j) = uaux(i,j)-gamx(i+1)*uaux(i+1,j)
              enddo
c
c second tridiagonal step
c
              uu(is) = -0.5d0
              do i=is+1,ie-1
                uu(i) = -c1m(i,j)*uu(i-1)/betx(i)
              enddo
              uu(ie) = (alp-c1m(ie,j)*uu(ie-1))/betx(ie)
              do i=ie-1,is,-1
                uu(i) = uu(i)-gamx(i+1)*uu(i+1)
              enddo
              fact = (uaux(is,j)+bet*uaux(ie,j)/gamma)
     &               /(one+uu(is)+bet*uu(ie)/gamma)
              do i=is,ie
                uaux(i,j) = uaux(i,j)-fact*uu(i)
              enddo
            endif
70        continue
c
c solve tridiagonal system in the y-direction
c
          do i=2,imaxm1
            bety(i)       = c20(i,njmin)+zaw(i,njmin)
            vaux(i,njmin) = (r2(i,njmin)+dvsg2(i,njmin))/bety(i)
          enddo
          do j=njmin+1,njmax
            do i=2,imaxm1
              gamy(i,j) = c2p(i,j-1)/bety(i)
              bety(i)   = c20(i,j)+zaw(i,j)-c2m(i,j)*gamy(i,j)
              vaux(i,j) = (r2(i,j)+dvsg2(i,j)-c2m(i,j)*vaux(i,j-1))
     &                    /bety(i)
            enddo
          enddo
          do j=njmax-1,njmin,-1
            do i=2,imaxm1
              vaux(i,j) = vaux(i,j)-gamy(i,j+1)*vaux(i,j+1)
            enddo
          enddo
c
c apply over/under-relaxation
c
          do j=njmin,njmax
            do i=2,imaxm1
              uaux(i,j) = up(i,j)+om*(uaux(i,j)-up(i,j))
              vaux(i,j) = vp(i,j)+om*(vaux(i,j)-vp(i,j))
            enddo
          enddo
c
          if (ih.eq.1) then
c
c case of bering
c
            uaux(ibera,jbera) = zero
            vaux(ibera,jbera) = -ubering
            uaux(iberp,jberp) = zero
            vaux(iberp,jberp) = ubering
          endif
c
c cyclicity
c
          do j=njmin,njmax
            uaux(1,j)    = uaux(imaxm1,j)
            uaux(imax,j) = uaux(2,j)
            vaux(1,j)    = vaux(imaxm1,j)
            vaux(imax,j) = vaux(2,j)
          enddo
c
c convergence test
c
          resm = zero
          do j=njmin,njmax
            do i=iu1(j),iu2(j)
              resm = max(resm,abs(uaux(i,j)-up(i,j)),
     &                        abs(vaux(i,j)-vp(i,j)))
            enddo
          enddo
          if(resm.lt.resl) go to 2000
c
c new stress tensor
c
          do j=njmm1,njmax
            jp1 = j+1
            do i=1,imax
              ip1  = (i+1)-(imax-2)*(i/imax)
c
c rate of strain tensor
c
              ze11 = akappa(i,j,1,1)*
     &               (uaux(ip1,j)+uaux(ip1,jp1)-(uaux(i,j)+uaux(i,jp1)))
     &              +akappa(i,j,1,2)*
     &               (vaux(ip1,j)+vaux(ip1,jp1)+vaux(i,j)+vaux(i,jp1))
              ze12 = akappa(i,j,2,2)*
     &               (uaux(i,jp1)+uaux(ip1,jp1)-(uaux(i,j)+uaux(ip1,j)))
     &              -akappa(i,j,2,1)*
     &               (vaux(i,j)+vaux(ip1,j)+vaux(i,jp1)+vaux(ip1,jp1))
              ze22 = akappa(i,j,2,2)*
     &               (vaux(i,jp1)+vaux(ip1,jp1)-(vaux(i,j)+vaux(ip1,j)))
     &              +akappa(i,j,2,1)*
     &               (uaux(i,j)+uaux(ip1,j)+uaux(i,jp1)+uaux(ip1,jp1))
              ze21 = akappa(i,j,1,1)*
     &               (vaux(ip1,j)+vaux(ip1,jp1)-(vaux(i,j)+vaux(i,jp1)))
     &              -akappa(i,j,1,2)*
     &               (uaux(ip1,j)+uaux(ip1,jp1)+uaux(i,j)+uaux(i,jp1))
c
              aa          = viszeta(i,j)*(ze11+ze22)
c
c is the relaxation of sigm really neccesary?
c the algorithm exhibits sometimes instabilities
c when the stress tensor is not relaxed. is
c this due to the fact that the scheme is a 9-point
c one? 
c
              sigm11(i,j) = sigm11(i,j)+om*( viseta(i,j)*(ze11-ze22)+aa
     &                                      -sigm11(i,j))
              sigm12(i,j) = sigm12(i,j)+om*( viseta(i,j)*(ze12+ze21)
     &                                      -sigm12(i,j))
              sigm22(i,j) = sigm22(i,j)+om*( viseta(i,j)*(ze22-ze11)+aa
     &                                      -sigm22(i,j))
            enddo
          enddo
c
1000    continue
2000    continue
c
c include coriolis term correction
c
        do j=njmin,njmax
          do i=iu1(j),iu2(j)
            ur        = uaux(i,j)-uo(i,j)
            vr        = vaux(i,j)-vo(i,j)
            zmod      = sqrt(ur*ur+vr*vr)*(one-albqd(i,j))
            zcd       = rhoco*zmod
            zbw(i,j)  = alpha*zcorl(i,j)+zcd*sang 
            t1        = zaw(i,j)*uaux(i,j)-zbw(i,j)*vg(i,j)
            t2        = zaw(i,j)*vaux(i,j)+zbw(i,j)*ug(i,j)
            det       = zaw(i,j)*zaw(i,j)+zbw(i,j)*zbw(i,j)
            zden      = sign(one,det)/max(zepsd2,abs(det))
     &                  *uvmsk(i,j)
            zrug      = (zaw(i,j)*t1+zbw(i,j)*t2)*zden
            zrvg      = (zaw(i,j)*t2-zbw(i,j)*t1)*zden
            ug(i,j)   = zindu*(ug(i,j)-zrug)+zrug
            vg(i,j)   = zindu*(vg(i,j)-zrvg)+zrvg
          enddo
        enddo
c
        if (ih.eq.1) then
c
c case of bering
c
          ug(ibera,jbera) = zero
          vg(ibera,jbera) = -ubering
          ug(iberp,jberp) = zero
          vg(iberp,jberp) = ubering
        endif
c
c cyclicity
c
        do j=njmin,njmax
          ug(1,j)    = ug(imaxm1,j)
          ug(imax,j) = ug(2,j)
          vg(1,j)    = vg(imaxm1,j)
          vg(imax,j) = vg(2,j)
        enddo
c
c initial conditions
c
        do j=njmin,njmax
          do i=iu1(j),iu2(j)
            u0(i,j) = 2.0d0*zindu*(u0(i,j)-ug(i,j))+ug(i,j)
            v0(i,j) = 2.0d0*zindu*(v0(i,j)-vg(i,j))+vg(i,j)
          enddo
        enddo
c
c       if (iter.eq.2.or.iter.eq.2*nbiter) write(79,*) ih,iter,resm,jter
c
3000  continue
c
c     normalised rate of energy disipation by 
c     shear and convergence ---> extra open-water production 
c
      do j=njmm1,njmax
        jp1 = j+1
        do i=1,imax
          ip1       = (i+1)-(imax-2)*(i/imax)
          ze11      = akappa(i,j,1,1)*
     &                 (ug(ip1,j)+ug(ip1,jp1)-(ug(i,j)+ug(i,jp1)))
     &               +akappa(i,j,1,2)*
     &                 (vg(ip1,j)+vg(ip1,jp1)+vg(i,j)+vg(i,jp1))
          ze12      = akappa(i,j,2,2)*
     &                 (ug(i,jp1)+ug(ip1,jp1)-(ug(i,j)+ug(ip1,j)))
     &               -akappa(i,j,2,1)*
     &                 (vg(i,j)+vg(ip1,j)+vg(i,jp1)+vg(ip1,jp1))
          ze22      = akappa(i,j,2,2)*
     &                 (vg(i,jp1)+vg(ip1,jp1)-(vg(i,j)+vg(ip1,j)))
     &               +akappa(i,j,2,1)*
     &                 (ug(i,j)+ug(ip1,j)+ug(i,jp1)+ug(ip1,jp1))
          ze21      = akappa(i,j,1,1)*
     &                 (vg(ip1,j)+vg(ip1,jp1)-(vg(i,j)+vg(i,jp1)))
     &               -akappa(i,j,1,2)*
     &                (ug(ip1,j)+ug(ip1,jp1)+ug(i,j)+ug(i,jp1))
          trace     = ze11+ze22
          trace2    = trace*trace
          deter     = ze11*ze22-0.25d0*(ze12+ze21)**2
          delta     = sqrt(trace2+(trace2-4.0d0*deter)*usecc2)
          alcr(i,j) = (1.0d0-max(0.0d0,sign(1.0d0,-(1.0d0-albq(i,j)))))*
     &                0.5*(delta-trace)*wstrn(i,j)
        enddo
      enddo
c
      return
      end
