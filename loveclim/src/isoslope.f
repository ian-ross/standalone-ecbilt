












      subroutine isoslope
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c======================================================================
c  Computation of isopycnal slope : SLOPE CUT-OFF Limiter (Cox, 1987)
c======================================================================
c- Attention au signe !  drho,y,z = - d(rho) / dx,y,z ; idem pour drox,y,z1,2,4
c- modif : "epsil2" est ajoute dans "drhoz"
c  modif : 02/06/99
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'reper.com'
      include 'isoslope.com'
      include 'comunit.h'
 
      parameter ( ntabmx = 100 )
      common / isoslp / drhox(imax,jmax,kmax),
     & drhoy(imax,jmax,kmax), drhoz(imax,jmax,kmax),
     & rossby(imax,jmax),ztabf1(-ntabmx:ntabmx),ztabf2(-1:ntabmx),
C    & slopc, slopd, refcur, radmn, radmx, domdds, ddslop, ddrati,
     & slopc, unsdds, unsddr, ndd1mx, ndd1mn, ndd2mx
c- slopd ... ddrati => facultatif
 
c- n99=2 => ecritures auxiliaires sur fichier "mouchard", unit=99
      common / mchd99 / nn99
 
      epsil2 = epsil * epsil
 
      if (numit.eq.nstart) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- debut du traitement specifique 1ere it.
c-----------------------------------------------------------------------
 
c- initialisation des variables "isopyc + GM" en common :
      do 110 k=1,kmax
       do 110 j=1,jmax
        do 110 i=1,imax
          c1x(i,j,k) = zero
          c2y(i,j,k) = zero
          c4x(i,j,k) = zero
          c4y(i,j,k) = zero
          c4xgm(i,j,k) = zero
          c4ygm(i,j,k) = zero
          drhox(i,j,k) = zero
          drhoy(i,j,k) = zero
          drhoz(i,j,k) = zero
          uiso(i,j,k) = zero
          viso(i,j,k) = zero
          dscalz(i,j,k) = zero
 110  continue
      do 115 k=1,kmax+1
       do 115 j=1,jmax
        do 115 i=1,imax
          wiso(i,j,k) = zero
          fisoz(i,j,k) = zero
 115  continue
 
      if(ai(kmax).eq.zero.and.aitd(kmax).eq.zero) return
 
       write(iuo+66,'(2A,2F6.4)')
     &   ' *** isoslope ; iso,gm90 *** :',
     &   ' SlopMx_ISO,GM =', slopemax(kmax), slopmgm(kmax)
 
c- computation of new masks
      do 125 k=1,kmax
       do 125 j=1,jmax
        do 120 i=1,imax
          if (i.ne.1) ttm1(i,j,k) = tms(i,j,k)*tms(i-1,j,k)
          if (j.ne.1) then
            ttm2(i,j,k) = tms(i,j,k)*tms(i,j-1,k)
          else
            ttm2(i,j,k) = zero
          endif
 120    continue
        ttm1(1,j,k) = ttm1(ims2,j,k)
 125  continue
 
c- setup index jsdom1,jsdom2(i) to cover the computation domain :
      do 130 i=1,imax
        jsdom1(i) = jmax
        jsdom2(i) = 1
 130  continue
      do 135 j=js1,js2
       do 135 i=is1(j),is2(j)
         if (tms(i,j,ks2).eq.one) then
           jsdom1(i) = min(jsdom1(i),j)
           jsdom2(i) = max(jsdom2(i),j)
         endif
 135  continue
c- write index jsdom1,jsdom2 on file "mouchard"
      if (nn99.eq.2) then
        write(99,'(A)') ' isoslope ; i / jsdom1 / jsdom2 :'
        do 140 is=1,imax,18
          ie = min(imax,is+17)
          write(99,'(A,18I4)') '    i = ', (i,i=is,ie)
          write(99,'(A,18I4)') ' jsdom1=', (jsdom1(i),i=is,ie)
          write(99,'(A,18I4)') ' jsdom2=', (jsdom2(i),i=is,ie)
 140    continue
        write(99,*)
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- limiteur F1 & F2 : calcul de "Rossby_Radius" & tabulation de F1 & F2 :
c-----------------------------------------------------------------------
 
c- mise en place des coeffs (cf papier Large et al) :
c  slopc=Sc ; slopd=Sd ; refcur="c"=1rst.Barocl.Wave Speed ;
c  [radmn,radmx] = interval for Rossby radius (R)
c-  discretis.F1 : domaine : -domdds*Sd < S < domdds*Sd ; Pas=ddslop*Sd
c-  discretis.F2 : domaine : 0 < ratio=|z|/R*S < 1 ; Pas=ddrati
c NB: slopc, slopd, ... /Smax in file "run.param"
        slopc = slopmgm(kmax-1)
        slopd = slopmgm(kmax-2)
        refcur = slopmgm(kmax-3)
        radmn  = slopmgm(kmax-4)
        radmx  = slopmgm(kmax-5)
        domdds = slopmgm(kmax-6)
        ddslop = slopmgm(kmax-7)
        ddrati = slopmgm(kmax-8)
c- si Increment(ddslop or ddrati) Nul => pas de limiteur (F1=1 or F2=1)
        if (ddslop.gt.zero) then
          slopd  = max(slopd,epsil)
          domdds = max(domdds,ddslop)
          unsdds = 1.0 / (slopd*ddslop)
          ndd1mx = nint(domdds/ddslop)
          ndd1mn = -ndd1mx-1
        else
          slopd  = 1.
          domdds = 1.
          ddslop = 0.
          unsdds = 0.
          ndd1mx = 0
          ndd1mn = 0
          ztabf1(0) = 1.
        endif
        if (ddrati.gt.zero) then
          unsddr = 1.0 / ddrati
          ndd2mx = nint(unsddr)
        else
          unsddr = 0.
          ndd2mx = -1
          ztabf2(0) = 1.
        endif
c---------
        if (ndd1mx+1.gt.ntabmx .or. ndd2mx+1.gt.ntabmx) then
          write(iuo+66,'(2A)') 'STOP in "isoslope" :',
     &               ' array ztabf1,ztabf2 sous-dimensione !'
          write(iuo+66,'(2(A,2I6))') ' Max Index=', ndd1mx+1, ndd2mx+1,
     &                          ' > Dim=', ntabmx
          stop
        endif
c--------
        zz = 0.5 / tanh(domdds)
        do 160 nn=1+ndd1mn,ndd1mx
          xx = ddslop*nn
          ztabf1(nn) = 0.5 + zz*tanh(-xx)
 160    continue
        ztabf1(-ndd1mx-1) = ztabf1(-ndd1mx)
        ztabf1( ndd1mx+1) = ztabf1( ndd1mx)
c--------
        do 170 nn=0,ndd2mx
          xx = ddrati*nn
          ztabf2(nn) = 0.5 + 0.5*sin( pi*(xx-0.5) )
 170    continue
        ztabf2(-1) = ztabf2(0)
        ztabf2(ndd2mx+1) = ztabf2(ndd2mx)
c--------
        do 180 j=1,jmax
         do 180 i=1,imax
           rossby(i,j) = 0.5*refcur/max(abs(fs2cor(i,j)),epsil)
           rossby(i,j) = min(radmx,max(radmn,rossby(i,j)))
 180    continue
c--------
c- Sortie sur fichier "iso_tab.chk" des fonctions F1 & F2 :
C        include 'iso_slp_tab.inc'
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin du traitement specifique 1ere it.
      endif
 
      if (ai(kmax).eq.zero.and.aitd(kmax).eq.zero) return
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Compute (minus)GRADIENTS OF BUYANCY B .
c-----------------------------------------------------------------------
 
      do 220 k=ks1,ks2
        do 200 j=js1,js2+1
         do 200 i=isf1(j),isf2(j)+1
          drhox(i,j,k) = ttm1(i,j,k) *
     &                 (b(i-1,j,k)-b(i,j,k))*unsdx*smx(i,j,1)
          drhoy(i,j,k) = ttm2(i,j,k) *
     &                 (b(i,j-1,k)-b(i,j,k))*unsdy*smy(i,j,2)
 200    continue
        do 210 j=1,jmax
         do 210 i=1,imax
          drhoz(i,j,k) = max( bvf(i,j,k), epsil2 )
 210    continue
 220  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Start  external loop on all levels, index "k".
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(km1,km0,kp1,kp0,j,i)
C$DIR SPP LOOP_PRIVATE(droz1,droz2,drox4,droy4, slp4x,slp4y,slp4z)
C$DIR SPP LOOP_PRIVATE(nn,ssn,dss,ztap1, rrn,drr,ztap2, ccxy)
      do 800 k=ks1,ks2
c-----
      km1 = max(k-1, ks1)
      km0 = km1 + 1
      kp1 = min(k+1, ks2)
      kp0 = kp1 - 1
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Compute Isopycnal Slope at point 1 (for Diffus. Isopyc).        |
c-----------------------------------------------------------------------
 
      do 310 j=js1,js2
       do 310 i=is1(j),is2(j)+1
         droz1 = ( drhoz(i,j,km0) + drhoz(i-1,j,kp1)
     &           + drhoz(i,j,kp1) + drhoz(i-1,j,km0) ) /
     &           ( tms(i,j,km1) + tms(i-1,j,kp0)
     &           + tms(i,j,kp0) + tms(i-1,j,km1) + epsil )
         c1x(i,j,k) = ttm1(i,j,k) * drhox(i,j,k)
     &              / max( droz1, abs(drhox(i,j,k)/slopemax(k)) )
 310  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Compute Isopycnal Slope at point 2 (for Diffus. Isopyc).        |
c-----------------------------------------------------------------------
 
      do 410 j=js1,js2+1
       do 410 i=isf1(j),isf2(j)
         droz2 = ( drhoz(i,j,km0) + drhoz(i,j-1,kp1)
     &           + drhoz(i,j,kp1) + drhoz(i,j-1,km0) ) /
     &           ( tms(i,j,km1) + tms(i,j-1,kp0)
     &           + tms(i,j,kp0) + tms(i,j-1,km1) + epsil )
         c2y(i,j,k) = ttm2(i,j,k) * drhoy(i,j,k)
     &              / max( droz2, abs(drhoy(i,j,k)/slopemax(k)) )
 410  continue
 
      if (k.eq.ks1) goto 800
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Compute Isopycnal Slope at point 4 (for Diffus. Isopyc & GM scheme).
c-----------------------------------------------------------------------
 
      do 550 j=js1,js2
       do 550 i=is1(j),is2(j)
         drox4 = ( drhox(i,j,k) + drhox(i+1,j,km1)
     &           + drhox(i,j,km1) + drhox(i+1,j,k) ) /
     &           ( ttm1(i,j,k) + ttm1(i+1,j,km1)
     &           + ttm1(i,j,km1) + ttm1(i+1,j,k) + epsil )
         droy4 = ( drhoy(i,j,k) + drhoy(i,j+1,km1)
     &           + drhoy(i,j,km1) + drhoy(i,j+1,k) ) /
     &           ( ttm2(i,j,k) + ttm2(i,j+1,km1)
     &           + ttm2(i,j,km1) + ttm2(i,j+1,k) + epsil )
 
c- slope for Diffus. Isopyc. :
         c4x(i,j,k) = tms(i,j,k-1) * drox4
     &              / max( drhoz(i,j,k), abs(drox4/slopmgm(ks2)) )
         c4y(i,j,k) = tms(i,j,k-1) * droy4
     &              / max( drhoz(i,j,k), abs(droy4/slopmgm(ks2)) )
 
c- slope for GM scheme :
         slp4x = drox4 / max( drhoz(i,j,k), abs(drox4/slopmgm(ks2)) )
         slp4y = droy4 / max( drhoz(i,j,k), abs(droy4/slopmgm(ks2)) )
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c limit slope*Ai <- Large et al, 1997, JPO, Vol.27, p 2445+2446
c test : only applied on A_ITD
c --> ATTENTION : 1) put limit of A_ITD on slope_x,y  --> c4x & c4y
c and 2) use only slopmgm(kmax), for other k<kmax :
c       => slopc,slopd,refcur,Radmn,Radmx, domdds, ddslop, ddrati
c----------------------------------------------------------------------
C     unsdds = 1.0 / (slopd*ddslop)
C     ndd1mx = nint(domdds/ddslop)
C     ndd1mn = -ndd1mx-1
C     unsddr = 1.0 / ddrati
C     ndd2mx = nint(unsddr)
c---------
c- limiteur F1 (slp4z discretise, par pas de "ddslop") :
          slp4z = abs(slp4x) + abs(slp4y) + epsil2
          ssn = (slp4z - slopc) * unsdds
          nn = nint(ssn-0.5)
          dss = ssn - nn
          nn = min(ndd1mx,max(ndd1mn,nn))
          ztap1 = ztabf1(nn) + dss * (ztabf1(nn+1) - ztabf1(nn))
c- limiteur F2 (ratio discretise, par pas de "ddrati") :
C         ratio = -zw(k) / (rossby(i,j) * slp4z)
C         rrn = ratio * unsddr
          rrn = -zw(k)*unsddr / (rossby(i,j) * slp4z)
          rrn=min(rrn,1D9)
          nn = nint(rrn-0.5)
          drr = rrn - nn
          nn = min(ndd2mx,max(-1,nn))
          ztap2 = ztabf2(nn) + drr * (ztabf2(nn+1) - ztabf2(nn))
c---------
          ccxy = tms(i,j,k-1)*ztap1*ztap2
          c4xgm(i,j,k) = ccxy*slp4x
          c4ygm(i,j,k) = ccxy*slp4y
c---------
c- Sortie sur fichier "iso_slp.chk" de la colonne "icheck,jcheck" :
C         include 'iso_slp_out.inc'
c---------
 
 550  continue
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c- CYCLIC Raccord [c4x] for GM90
      if (ltest.ge.1) then
        do 610 j=jcl1,jcl2
          c4xgm(1,j,k) = c4xgm(ims2,j,k)
          c4xgm(imax,j,k) = c4xgm(ims1,j,k)
 610    continue
      endif
c- BERING Raccord [c4y] for GM90 (attention minus sign)
      if (ltest.eq.3) then
        c4ygm(iberpm,jberp,k) = -c4ygm(ibera, jberam,k)
        c4ygm(iberp, jberp,k) = -c4ygm(iberam,jberam,k)
        c4ygm(iberam,jbera,k) = -c4ygm(iberp, jberpm,k)
        c4ygm(ibera, jbera,k) = -c4ygm(iberpm,jberpm,k)
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- End of external loop on all levels, index "k".
 800  continue
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine isoslope -
      end
