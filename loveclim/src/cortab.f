












      subroutine cortab(ytab,ymodif,yspv,ymdspv,
     &                  im,is,ie,js,je,kmodif,kspv)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  Appele par "correct",
c  Modifie les valeurs du tableau ytab, entre les indice [is,ie]x[js,je]
c  suivant kmodif : 0,3 -> substituion, 1,4 -> addition ; 2,5 -> multiplication
c     0,1,2 par val.unique ymodif(1,1) ou 3,4,5 val.correspondante ymodif(i,j)
c          kspv   : 0 sans spv, 1 uniquement spv, 2 sauf spv.
c--------
c  modif : 22/09/97
 
      include 'type.com'
 
c- dummy variables :
      dimension ytab(im,*), ymodif(im,*)
c- variables locale
      logical flmod0, flmod1, flspv1, flspv2
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Initialisation :
c-----------------------------------------------------------------------
 
      flmod0 = mod(kmodif,3).eq.0
      flmod1 = mod(kmodif,3).eq.1
      flspv1 = kspv.eq.1
      flspv2 = kspv.eq.2
 
      if (kmodif.lt.3) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Traitement avec Valeur Unique :
c-----------------------------------------------------------------------
 
      do 210 j=js,je
       do 210 i=is,ie
         if (ytab(i,j).ne.yspv) then
c--Modifie que si = spv :
           if (flspv1) goto 210
         else
c--Ne modifie pas si = spv
           if (flspv2) goto 210
         endif
         if (flmod0) then
c--Substitution :
           ytab(i,j) = ymodif(1,1)
         elseif (flmod1) then
c--Addition :
           ytab(i,j) = ytab(i,j) + ymodif(1,1)
         else
c--Multiplication :
           ytab(i,j) = ytab(i,j) * ymodif(1,1)
         endif
 210  continue
 
      else
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Traitement avec Valeur correspondante :
c-----------------------------------------------------------------------
 
      if (kmodif.eq.0) then
c--Substitution :
        do 310 j=js,je
         do 310 i=is,ie
          if (ymodif(i,j).eq.ymdspv) goto 310
          if (ytab(i,j).ne.yspv) then
c- Modifie que si = spv :
            if(kspv.eq.1) goto 310
          else
c- Ne modifie pas si = spv
            if(kspv.eq.2) goto 310
          endif
          ytab(i,j) = ymodif(i,j)
 310    continue
c-----
 
      elseif (kmodif.eq.1) then
c--Addition :
        do 320 j=js,je
         do 320 i=is,ie
          if (ymodif(i,j).eq.ymdspv) goto 320
          if (ytab(i,j).ne.yspv) then
c- Modifie que si = spv :
            if(kspv.eq.1) goto 320
          else
c- Ne modifie pas si = spv
            if(kspv.eq.2) goto 320
          endif
          ytab(i,j) = ytab(i,j) + ymodif(i,j)
 320    continue
c-----
 
      else
c--Multiplication :
        do 330 j=js,je
         do 330 i=is,ie
          if (ymodif(i,j).eq.ymdspv) goto 330
          if (ytab(i,j).ne.yspv) then
c- Modifie que si = spv :
            if(kspv.eq.1) goto 330
          else
c- Ne modifie pas si = spv
            if(kspv.eq.2) goto 330
          endif
          ytab(i,j) = ytab(i,j) * ymodif(i,j)
 330    continue
c-----
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
      return
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine cortab -
      end
