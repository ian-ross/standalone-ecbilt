












      subroutine correct(nn99,filcor)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  Appelee par "redforc" : Correction du forcage selon le fichier filcor.
c  Convention : ktyp < 3 remplace , 2 < ktyp < 6 additione , 5 < ktyp multiplie
c    mod(ktyp,3) = 0 spv sans effet, = 1 modif si spv, = 2 modif sauf si spv
c  modification locale ([is,ie]x[js,je] -> 2nd ligne) definie par valeur unique
c    ou globale par valeur correspondante lue sur fichier(-> 2nd ligne).
c--------
c  modif : 04/08/98
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
 
c- dummy variables :
      character*(*) filcor
 
c- variables locales equivalentes :
      dimension phisuv(imax,jmax,2)
      equivalence ( phisuv(1,1,1) , phisu(1,1) )
      equivalence ( phisuv(1,1,2) , phisv(1,1) )
 
c- variables locales :
      dimension ycor(imax,jmax)
      character*40 filtab
      character*70 line
 
 1200 format(A,3(2I4,A),1PE11.3)
 1300 format(A,2(2I4,A),I3,2I2,A,1PE11.3)
 1400 format(A,2(2I4,A),2I3,A,2I2,A,1PE11.3)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Initialisation ; Ouverture du fichier "forc.corr" .
c-----------------------------------------------------------------------
 
      nrapmx = 0
 
      if (nn99.eq.2) write(99,'(2A)')
     &  'Forcage Modifie par le fichier ', filcor
 
      open(unit = 29, file = filcor, status = 'old')
      read(29,*)
      read(29,*)
      read(29,*)
      read(29,*)
      read(29,*)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
      read(29,*) ntabc
      do 500 nt=1,ntabc
c--Definition du type de modification :
c  Nbcor = Nb de lignes (a lire) pour definir les modifs ;
c  nvcor = tableau modifie ; kcor = 3eme indice (ex:niveau k) du Tab modifie
c  ktyp : cf en-tete ; spv = Special-Value a prendre en compte.
c-----
      read(29,'(A)') line
      if (nn99.eq.2) write(99,'(A)') line
      read(29,*) nbcor, nvcor, kcor, ktyp, spv
      kspv = mod(ktyp,3)
      kmodif = ktyp / 3
 
      do 400 nc=1,max(nbcor,min(1,-nbcor))
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Modif par portion ou en totalite d'un tableau (2D) de forcage :
c-----------------------------------------------------------------------
 
      if (nbcor.lt.0) then
c- modif de la totalite d'un Tab.2D (Nouvelles valeurs lues sur fichier) :
        read(29,'(A)') filtab
        icor1 = 1
        icor2 = imax
        jcor1 = 1
        jcor2 = jmax
        kmodif = kmodif + 3
c- lecture(1 niveau complet) du fichier "filtab" :
        open(39,file=filtab,status='OLD',form='UNFORMATTED')
        read(39) spvcor
        read(39) ycor
        close(39)
        if (nn99.eq.2) write(99,'(2A)')
     &                ' -> modif globale, fichier = ', filtab
      elseif (nvcor.eq.6) then
c-----
c- ajoute (dans liste specifique) un forcage local :
c- NB : read i1,i2 + j1,j2 + k1,k2 + Y_new
        read(29,*) icor1, icor2, jcor1, jcor2, kcor1, kcor2, ycor(1,1)
        icor2 = max(icor1,icor2,-imax*icor2)
        icor1 = max(icor1,1)
        jcor2 = max(jcor1,jcor2,-jmax*jcor2)
        jcor1 = max(jcor1,1)
        kcor2 = max(kcor1,kcor2,-kmax*kcor2)
        kcor1 = max(kcor1,1)
        spvcor = spv
c-----
      else
c- modif d'un Tab.2D, par portion :
        read(29,*) icor1, icor2, jcor1, jcor2, ycor(1,1)
        icor2 = max(icor1,icor2,-imax*icor2)
        icor1 = max(icor1,1)
        jcor2 = max(jcor1,jcor2,-jmax*jcor2)
        jcor1 = max(jcor1,1)
        spvcor = spv
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Traitements des differents cas de tableaux a modifier :
c-----------------------------------------------------------------------
 
      if (nvcor.eq.1) then
          if (kcor.lt.1 .or. kcor.gt.2) goto 910
          call cortab(phisuv(1,1,kcor),ycor,spv,spvcor,
     &                imax,icor1,icor2,jcor1,jcor2,kmodif,kspv)
          if (nn99.eq.2) write(99,1200) ' modif phiss(ns) : i',
     &      icor1,icor2,' ; j', jcor1,jcor2, ' ; u/v,ktyp', kcor,ktyp,
     &      ' ; V.cor=', ycor(1,1)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      elseif (nvcor.eq.2) then
          if (kcor.lt.0 .or. kcor.gt.nsmax) goto 910
          call cortab(phiss(1,1,kcor),ycor,spv,spvcor,
     &                imax,icor1,icor2,jcor1,jcor2,kmodif,kspv)
          if (nn99.eq.2) write(99,1200) ' modif phiss(ns) : i',
     &      icor1,icor2,' ; j', jcor1,jcor2, ' ;  ns,ktyp', kcor,ktyp,
     &      ' ; V.cor=', ycor(1,1)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      elseif (nvcor.eq.3) then
c-- kcor = ns * k
          if (kcor.lt.1 .or. kcor.gt.nsmax*kmax) goto 910
          k = kcor - 1
          ns = k / kmax + 1
          k = mod(k,kmax) + 1
          call cortab(scalr(1,1,k,ns),ycor,spv,spvcor,
     &                imax,icor1,icor2,jcor1,jcor2,kmodif,kspv)
          if (nn99.eq.2) write(99,1300) ' modif scalr(ns) : i',
     &      icor1,icor2,' ; j', jcor1,jcor2, ' ; k,ns,ktyp', k,ns,ktyp,
     &      ' ; V.cor=', ycor(1,1)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      elseif (nvcor.eq.4) then
          if (kcor.lt.1 .or. kcor.gt.kmax) goto 910
          call cortab(rappel(1,1,kcor),ycor,spv,spvcor,
     &                imax,icor1,icor2,jcor1,jcor2,kmodif,kspv)
          if (nn99.eq.2) write(99,1200) ' modif rappel(k) : i',
     &      icor1,icor2,' ; j', jcor1,jcor2, ' ;  k, ktyp', kcor,ktyp,
     &      ' ; V.cor=', ycor(1,1)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      elseif (nvcor.eq.5) then
          if (kcor.lt.0 .or. kcor.gt.nsmax) goto 910
          call cortab(rappes(1,1,kcor),ycor,spv,spvcor,
     &                imax,icor1,icor2,jcor1,jcor2,kmodif,kspv)
          if (nn99.eq.2) write(99,1200) ' modif rappes(ns): i',
     &      icor1,icor2,' ; j', jcor1,jcor2, ' ;  ns,ktyp', kcor,ktyp,
     &      ' ; V.cor=', ycor(1,1)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      elseif (nvcor.eq.6) then
c---Rappel (explic.) autre qu'en surface (--> tableau specifique) :
          if (kcor.lt.1 .or. kcor.gt.nsmax) goto 910
          ns = kcor
          nnadd = (icor2 - icor1 + 1) * (jcor2 - jcor1 + 1)
          do 365 k=kcor1,kcor2
            nn = nrap(k,ns)
            nrap(k,ns) = nrap(k,ns) + nnadd
            nrapmx = max(nrapmx, nrap(k,ns))
            if (nrapmx.le.nrpmax) then
              do 360 j=jcor1,jcor2
               do 360 i=icor1,icor2
                 nn = nn + 1
                 ijrap(nn,k,ns) = (j - 1)*imax + i - 1
                 rapint(nn,k,ns) = ycor(1,1)
 360          continue
            endif
 365      continue
          if (nn99.eq.2) write(99,1400) ' rapint(k,ns): i',
     &      icor1,icor2,' ; j', jcor1,jcor2, ' ; k', kcor1,kcor2,
     &      ' ; ns,typ', ns,ktyp, ' ; V.cor=', ycor(1,1)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      elseif (nvcor.eq.7) then
c-- kcor = 0,1 min/Max Frsh-W-Flx ; 2,3 min/Max Heat-Flx ; 4,5 min/Max Salt-Flx
          if (kcor.lt.0 .or. kcor.ge.3*nsmax) goto 910
          ns = kcor/2
          nx = mod(kcor,2)
          call cortab(phimnx(1,1,nx,ns),ycor,spv,spvcor,
     &                imax,icor1,icor2,jcor1,jcor2,kmodif,kspv)
          if (nn99.eq.2) write(99,1300) ' modif phimnx(ns): i',
     &      icor1,icor2,' ; j', jcor1,jcor2, ' ; n,ns,ktyp',nx,ns,ktyp,
     &      ' ; V.cor=', ycor(1,1)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      else
          goto 900
      endif
 
c--fin du traitement d'un tableau 2 D.
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 400  continue
 
c--fin du traitement de la modificcation.
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 500  continue
 
      close(29)
 
      if (nrapmx.le.nrpmax) then
        write(iuo+66,'(2A)') 'Forcage Modifie par le fichier ', filcor
        return
      else
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  9 ) Traitements des cas d'erreurs .
c-----------------------------------------------------------------------
 
        write(iuo+66,'(A,I10,A)')'STOP in "correct", dimension nrpmax=',
     &              nrpmax, '= TOO SMALL !'
        write(iuo+66,'(2A)')     '=> Change file : ', filcor
        write(iuo+66,'(A,I10,A)')'or set nrpmax (in "para.Com") to at least'
     &            , nrapmx, ' and compile the code again.'
        stop
      endif
 
 900  continue
      write(iuo+66,'(2A)') 'STOP in "correct", ERROR in : ',filcor
      write(iuo+66,'(2A)') ' modif : ', line
      write(iuo+66,'(A,I8,A)') ' var. nvcor=', nvcor, '  <-- Not found !'
      stop
 
 910  continue
      write(iuo+66,'(2A)') 'STOP in "correct", ERROR in : ',filcor
      write(iuo+66,'(2A)') ' modif : ', line
      write(iuo+66,'(A,I3,A,I8,A)') ' var. nvcor=', nvcor,
     &   ' index kcor=', kcor, '  <-- out of range !'
      stop
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine correct -
      end
