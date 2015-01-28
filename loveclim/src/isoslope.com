












c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  fichier "isoslope.com" : incorpore par instruction 'include' dans :
c    isoadvec isodiffu isoslope scali [iso_copy iso_flux isofilter]
c (commun a toutes les routines de diffus. Isopycnale et Genk & McWillliams)
c regroupe isoslope.com + isoadvec.com de P.P.M. (vers. du 23/03/97 & 16/10/96)
c-----
c  modif : 25/05/99
 
c  ppmodif : 23-03-97 : iso,gm90  (<- fichier isoslope.com d'origine)
c------------------------------------------------------------------
c GLOBAL VARIABLES for ISOPYCNAL DIFFUSION
c------------------------------------------------------------------
c              B-GRID faces (V=velocity) (S=scalar)
c------------------------------------------------------------------
c        WEST      SOUTH   MIDDLE   WEST     SOUTH     MIDDLE
c        + + +     + + +   + + +    + + +    + + +     V + V
c        + 1 3     3 2 +   1 0 +    V S V    V S V     + S +
c        + 5 +     + 6 +   + 2 +    + + +    + + +     V + V
c------------------------------------------------------------------
 
c  ai(k) = isopycnal mixing coefficient (m^2/s) (added to ahh)
c  slopemax(k) = maximum slope of isopycnals for ISO  (iroutine=1)
c  slopmgm(k)  = maximum slope of isopycnals for GM90 (iroutine=2)
c  ttm1=mask pt 1, ttm2=mask pt 2
c  c1x=drox/droz,c2y=droy/droz  where (droz=bvf)
      common / isopycnale /
     &  ttm1(imax,jmax,kmax), ttm2(imax,jmax,kmax),
     &  c1x(imax,jmax,kmax), c2y(imax,jmax,kmax),
     &  c4x(imax,jmax,kmax), c4y(imax,jmax,kmax),
     &  fisoz(imax,jmax,kmax+1), dscalz(imax,jmax,kmax)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  ppmodif : 16-10-96 : gm90  (<- fichier isoadvec.com ajoute ici)
c  aitd(k) = isopycnal thickness coefficient (m^2/s)
      common / isobolus /
     &  c4xgm(imax,jmax,kmax),c4ygm(imax,jmax,kmax)
 
      common / isobolus2 /
     &  uiso(imax,jmax,kmax),viso(imax,jmax,kmax),
     &  wiso(imax,jmax,kmax+1)

      common / indexiso /
     &  jsdom1(imax), jsdom2(imax)
 
c--fin du fichier "isoslope.com"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
