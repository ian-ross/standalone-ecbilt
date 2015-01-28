












      subroutine sepgl(vinp1, vinp2, voutw, vouta,
     &                 spvbin, spvout, ioutw, joutw, iouta, jouta,
     &                 kmniv, kref, kexcl, ijdl, jeqm, ksgn)
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c   Appele par "ncdfout".
c   Separation des 2 grilles dans 2 tableaux voutw & vouta
c-------
c  modif : 09/08/94
 
      include 'type.com'
C+    include 'const.com'
      include 'para.com'
      include 'reper.com'
 
c--dummy variables :
      real*4 spvbin, spvout, voutw, vouta
      dimension vinp1(imax,jmax), vinp2(imax,jmax)
      dimension kmniv(imax,jmax)
      dimension voutw(ioutw,joutw), vouta(iouta,jouta)
 
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  0 ) Initialisation des tableaux avec "spvout" :
c-----------------------------------------------------------------------
 
      do 10 j=1,joutw
       do 10 i=1,ioutw
        voutw(i,j) = spvout
 10   continue
 
      do 20 j=1,jouta
       do 20 i=1,iouta
        vouta(i,j) = spvout
 20   continue
 
      if (ltest.le.2) then
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  1 ) traitement du tableau grille WW , simple transfert :
c-----------------------------------------------------------------------
 
      if (ksgn.eq.2) then
c--Traitement particulier pour w :
        do 110 j=1,joutw
         do 110 i=1,ioutw
          if(kmniv(i,j).le.kref) then
            voutw(i,j) = vinp1(i,j) - vinp2(i,j)
          elseif(kmniv(i,j).le.kexcl) then
            voutw(i,j) = spvbin
          endif
 110    continue
      else
        do 130 j=1,joutw
         do 130 i=1,ioutw
          if(kmniv(i,j).le.kref) then
            voutw(i,j) = vinp1(i,j)
          elseif(kmniv(i,j).le.kexcl) then
            voutw(i,j) = spvbin
          endif
 130    continue
      endif
 
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  2 ) traitement du tableau grille AA , simple transfert :
c-----------------------------------------------------------------------
 
      if (ksgn.eq.2) then
c--Traitement particulier pour w :
        do 210 ja=1,jouta
         j = ja + jsepar - 1
         do 210 i=1,iouta
          if(kmniv(i,j).le.kref) then
            vouta(i,ja) = vinp1(i,j) - vinp2(i,j)
          elseif(kmniv(i,j).le.kexcl) then
            vouta(i,ja) = spvbin
          endif
 210    continue
      else
        do 230 ja=1,jouta
         j = ja + jsepar - 1
         do 230 i=1,iouta
          if(kmniv(i,j).le.kref) then
            vouta(i,ja) = vinp1(i,j)
          elseif(kmniv(i,j).le.kexcl) then
            vouta(i,ja) = spvbin
          endif
 230    continue
      endif
 
c--
      else
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  3 ) traitement du tableau grille WW , transfert si j =< jsep(i) :
c-----------------------------------------------------------------------
 
      jeq   = jeqm + 2
 
      if (ksgn.eq.2) then
c--Traitement particulier pour w :
        do 310 i=1,ioutw
C        jjm = min(jsep(i), joutw)
         jjm = min(max(jsep(i)-1, jeq), joutw)
         do 310 j=1,jjm
          if(kmniv(i,j).le.kref) then
            voutw(i,j) = vinp1(i,j) - vinp2(i,j)
          elseif(kmniv(i,j).le.kexcl) then
            voutw(i,j) = spvbin
          endif
 310    continue
      else
        do 330 i=1,ioutw
         jjm = min(max(jsep(i)-1, jeq), joutw)
         do 330 j=1,jjm
          if(kmniv(i,j).le.kref) then
            voutw(i,j) = vinp1(i,j)
          elseif(kmniv(i,j).le.kexcl) then
            voutw(i,j) = spvbin
          endif
 330    continue
      endif
 
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  4 ) traitement du tableau grille AA , retournement :
c-----------------------------------------------------------------------
 
      if (ksgn.eq.2) then
c--Traitement particulier pour w :
        do 410 ja=1,jouta
         i = ijdl - ja
         jj = jsep(i) - jeq
         iim = jj + min( 2 , jj + 1 )
         do 410 ia=iim,iouta
           j = jeqm + ia
           if (kmniv(i,j).le.kref) then
             vouta(ia,ja) = vinp1(i,j) - vinp2(i,j)
           elseif(kmniv(i,j).le.kexcl) then
             vouta(ia,ja) = spvbin
           endif
 410    continue
      elseif (ksgn.eq.-1) then
c--Avec changement de signe :
        do 430 ja=1,jouta
         i = ijdl - ja
         jj = jsep(i) - jeq
         iim = jj + min( 2 , jj + 1 )
         do 430 ia=iim,iouta
           j = jeqm + ia
           if (kmniv(i,j).le.kref) then
             vouta(ia,ja) = -vinp2(i,j)
           elseif(kmniv(i,j).le.kexcl) then
             vouta(ia,ja) = spvbin
           endif
 430    continue
      else
c--Sans changement de signe :
        do 450 ja=1,jouta
         i = ijdl - ja
         jj = jsep(i) - jeq
         iim = jj + min( 2 , jj + 1 )
         do 450 ia=iim,iouta
           j = jeqm + ia
           if (kmniv(i,j).le.kref) then
             vouta(ia,ja) = vinp2(i,j)
           elseif(kmniv(i,j).le.kexcl) then
             vouta(ia,ja) = spvbin
           endif
 450    continue
      endif
 
c--
      endif
 
      return
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c- fin de la routine sepgl -
      end
