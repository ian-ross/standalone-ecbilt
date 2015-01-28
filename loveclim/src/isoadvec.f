












      subroutine isoadvec
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- // BALANCE of DIAGONEL FLUXES //
c ============================================================
c   Implementation of Isopycnal thickness diffusion
c   References :
c   -Danabasoglu, G., McWilliams, J.C., P. Gent, 1994:
c    The role of mesoscale tracer transports in the global
c    ocean circulation. Science, 24, 1123-1128
c   -Gent, P.R., and J.C. McWilliams, 1990: Isopycnal mixing
c    in ocean circulation models. JPO, 20, 150-155
c ============================================================
c  modif : 25/05/99
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'isoslope.com'
      include 'comunit.h'
 
c--variables locales conservees d'un appel a l'autre :
      common / bolusloc / aitds2(kmax)
 
c--variables locales :
      dimension psix(imax), psiy(jmax)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      if(numit.eq.nstart) then
        write(iuo+66,'(2A,F5.0)') ' *** isoadvec ;  G&M.90  *** :',
     &   ' aitd=', aitd(kmax)
 
c- initialisation :
        do 10 k=1,kmax
          aitds2(k) = 0.5 * aitd(k)
 10     continue
 
c- fin du traitement 1ere iter.
      endif
 
c- compute slopes c4xgm & c4ygm (position 4):
Cpp   call isofilter(c4xgm)
Cpp   call isofilter(c4ygm)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Start  external loop on all Meridional Section, index "j".
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(i,k,psix)
      do 300 j=js1,js2
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Merid. Section : Compute "Bolus velocity" uiso & wiso(1rst part).
c-----------------------------------------------------------------------
 
c- initialisation :
      do 210 i=is1(j),is2(j)+1
        uiso(i,j,ks1) = 0.
 210  continue
 
c- Compute Stream Function "psix" (a dy pres) - fill in uiso & wiso :
c- NB : ai*alpha2 is nul at the bottom
 
      do 230 k=ks1+1,ks2
        do 220 i=is1(j),is2(j)+1
          psix(i) = ttm1(i,j,k-1) * aitds2(k) *
     &            ( c4xgm(i,j,k) + c4xgm(i-1,j,k) )
          uiso(i,j,k-1) = uiso(i,j,k-1) + unsdz(k-1)*psix(i)
          uiso(i,j,k) = -unsdz(k)*psix(i)
          psix(i) = cmy(i,j,1)*psix(i)
 220    continue
        do 225 i=is1(j),is2(j)
          wiso(i,j,k) = unsdx * (psix(i+1)-psix(i))
 225    continue
 230  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- End of external loop on all Meridional Section, index "j".
 300  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Start  external loop on all Zonal Section, index "i".
C$DIR SPP LOOP_PARALLEL
C$DIR SPP LOOP_PRIVATE(j,k,psiy)
      do 400 i=ims1,ims2
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Zonal Section : Compute "Bolus velocity" viso & wiso(2nd part).
c-----------------------------------------------------------------------
 
c- initialisation :
      do 310 j=jsdom1(i),jsdom2(i)+1
        viso(i,j,ks1) = 0.
 310  continue
 
c- Compute Stream Function "psiy" (a dx pres) - fill in viso & wiso :
c- NB : ai*alpha2 is nul at the bottom
 
      do 330 k=ks1+1,ks2
        do 320 j=jsdom1(i),jsdom2(i)+1
          psiy(j) = ttm2(i,j,k-1) * aitds2(k) *
     &            ( c4ygm(i,j,k) + c4ygm(i,j-1,k) )
          viso(i,j,k-1) = viso(i,j,k-1) + unsdz(k-1)*psiy(j)
          viso(i,j,k) = -unsdz(k)*psiy(j)
          psiy(j) = cmx(i,j,2)*psiy(j)
 320    continue
        do 325 j=jsdom1(i),jsdom2(i)
          wiso(i,j,k) = -smxy(i,j,0) *
     &      ( wiso(i,j,k) + unsdy * (psiy(j+1)-psiy(j)) )
 325    continue
 330  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- End of external loop on all Zonal Section, index "i".
 400  continue
 
Cpp   call debugrac(uiso,viso,wiso,'GM90_RACCORD')
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine isoadvec -
      end
