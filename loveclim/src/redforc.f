












      subroutine redforc(nn99)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c mise en place des tableaux servant au forcage du modele (ocean uniquement)
c cas nn99=2 : impression de controle sur fichier "mouchard" ;
c Preparation du rappel en surface ou en profondeur.
c kforc -> quels fichiers de donnees a lire.
c  modif : 30/09/99
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'reper.com'
      include 'comunit.h'
 
c--variables locales :
 
      character*50 ccfile
      character*70 line
      dimension zrap(imax,jmax)
 
      line = ' '
      nnline = 1
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Forcage Constant (ex: moyenne annuelle) .                       |
c-----------------------------------------------------------------------
 
c---------------
c Options : kforc(for all Scalar) ;
c kk=mod(abs(kforc),100) -> Constant Field ; kforc/100 -> Time Dependant
c     kforc >= 0   : read file "tsobs.om" (spvr+T+S for all levels )
c 1+mod(kk-1,nsmax)  : Nb of Scalar to read ; and for each of theses :
c 1+(kk-1/nsmax) = k : Nb of level, from the surface to the bottom
c  File_name : "ts1up.om" (k=1), "ts2up.om" (k=2), ...
c-
c concerning each scalar (N), yforc(N) :
c   non.zero : read file "flux"N".om", and multiply surface flux by yflux.
c---------------
 
 
c--2.1 Lecture du forcage Dynamique :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--tension de vent en surface (Wind Stress) :
Cic0  open(unit=22,file='wsxy.om',status='old',form='UNFORMATTED')
Cic0  read (22, end=900) phisu
Cic0  read (22, end=900) phisv
Cic0  close(22)
Cic0  line = line(:nnline)//'wsxy.om '
Cic0  nnline = nnline + 8
 
c--2.2 Lecture des temperatures et salinites observees :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c- spvr (=Special Value) doit verifier, pour tout scalaire ns :
c  spvr < mini { Inf(Scal)ns - kmax , 2*Inf(Scal)ns - Sup(Scal)ns }
      spvr = -100.
 
      kkyf = mod(abs(kforc),100)
      if (kforc.ge.0) then
c- Lecture de spvr & kmax*2 niveaux :
        kknb = 2 * ijkmax
        open(unit=24,file='tsobs.om',status='old',form='UNFORMATTED')
        read (24, end=900) spvr
        call redtab( scalr(1,1,1,1), kknb, 24 , kkerr)
        if (kkerr.ne.1) goto 900
        close(24)
        line = line(:nnline)//'tsobs.om '
        nnline = nnline + 9
c- precaution : valeur speciale si rappel non nul :
        do 215 ns=3,nsmax
         do 215 k=1,kmax
          if (rapp1(k).ne.0.) then
            do 210 j=1,jmax
             do 210 i=1,imax
               scalr(i,j,k,ns) = spvr
 210        continue
          endif
 215    continue
      endif
 
      if (kkyf.eq.0 .or. kkyf.gt.kmax*nsmax) then
        k2yf = 0
        n2yf = 0
        kk1 = kmax
      else
c- n2yf = Nb de scalaire a lire ;
c- k2yf = Nb de niveaux(pour chaque Scal) a partir de la surface ;
c  et donc  kkyf = (k2yf-1)*nsmax + n2yf
        k2yf = 1 + (kkyf-1) / nsmax
        n2yf = 1 + mod(kkyf-1,nsmax)
 
c- Lecture de spvr & (k2yf*n2yf) niveaux :
        kk1 = kmax - k2yf + 1
        write(ccfile,'(A,I1,A)') 'ts', k2yf, 'up.om'
 
        open(unit=24, file=ccfile, status='old', form='UNFORMATTED')
        read (24) spvr2
        if (kforc.ge.0 .and. spvr.ne.spvr2) then
          write(iuo+66,*) 'ARRET dans routine Redforc :'
          write(iuo+66,*) 'SPVR Differentes entre tsobs.om et '//ccfile
          close(24)
          stop
        endif
        spvr = spvr2
        do 240 ns=1,n2yf
         do 240 k=kk1,kmax
          call redtab( scalr(1,1,k,ns), ixjmax, 24 , kkerr)
          if (kkerr.ne.1) goto 900
 240    continue
        close(24)
        line = line(:nnline)//ccfile
        nnline = nnline + 9
      endif
 
      if (kforc.ge.0) then
        kk1 = 1
      else
        kk2 = kmax - k2yf
        do 250 ns=1,nsmax
         if (ns.gt.n2yf) kk2 = kmax
         do 250 k=1,kk2
          do 250 j=1,jmax
           do 250 i=1,imax
             scalr(i,j,k,ns) = spvr
 250    continue
      endif
 
c correction celcius-kelvin
      do 255 k=ks1,ks2
       do 255 j=1,jmax
        do 255 i=1,imax
          if (scalr(i,j,k,1).ne.spvr) scalr(i,j,k,1) =
     &                                scalr(i,j,k,1) + 273.15d0
 255  continue
 
c--2.3 Mise en place des tableaux "rappes & rappel" (A Delta-t*unstyr pres) :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--Rappel explicite differencie T / S : Rappes .
        do 260 ns=0,nsmax
         do 260 j=1,jmax
          do 260 i=1,imax
           rappes(i,j,ns) = rapp0(ns)
 260    continue
 
c--Rappel implicite commun a T & S : Rappel .
      do 265 k=kk1,kmax
       do 265 j=1,jmax
        do 265 i=1,imax
          if(j.eq.js1.and.kbath1(i).ge.(ks2+ks1-k)) then
c- Rappel sur le mur S :
            rappel(i,j,k) = max(rapp1(k), rapp0(nsmax+2))
          elseif(j.eq.js2.and.kbath2(i).ge.(ks2+ks1-k)) then
c- Rappel sur le mur N :
            rappel(i,j,k) = max(rapp1(k), rapp0(nsmax+1))
          else
c- Rappel a l'interieur du bassin :
            rappel(i,j,k) = rapp1(k)
          endif
 265  continue
 
c--2.4 Flux en surface pour chaque scalaire :
c------------------------------------------------------------------
 
c- Convention Flux : + vers le Haut ; unites : Scal() x L(m) / Temps(s)
c-  dans modele : phiss = (Flux->haut) x (Delta_T) / Dz(1er_niv)
 
      do 280 ns=1,nsmax
        if (yforc(ns).ne.zero) then
c- lecture fichier flux :
          write(ccfile,'(A,I1,A)') 'flux', ns, '.om'
          open(unit=28,file=ccfile,status='old',form='UNFORMATTED')
          call redtab( phiss(1,1,ns), ixjmax, 28 , kkerr)
          close(28)
          line = line(:nnline)//ccfile
          nnline = nnline + 9
        endif
 280  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Forcage Saisonnier / Mensuel .                                  |
c-----------------------------------------------------------------------
 
c--LECTURE des donnes saisoniers de T, S et tension du vent
c--tension de vent mensuelle en surface (Wind Stress) :
Csai  open(unit=23,file='wsxy.mens.3x3.om',
Csai &    status='old',form='UNFORMATTED')
Csai   read(23) txmens
Csai   read(23) tymens
Csai  close(23)
 
c--temperature mensuelle et salinite saisoniere observees :
Csai  open(unit=25,file='tlev.mens.3x3.om',
Csai &    status='old',form='UNFORMATTED')
Csai   read(25) tmens
Csai  close(25)
 
Csai  open(unit=26,file='slev.sais.3x3.om',
Csai &    status='old',form='UNFORMATTED')
Csai   read(26) smens
Csai  close(26)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  7 ) Modification du Forcage, Conversion pour utilisation directe .  |
c-----------------------------------------------------------------------
 
c- modification localisee par appel a "correct" :
      if (mdforc.eq.1) call correct(nn99,filcor)
 
c- modification du rappel explicite par la lecture de filcor
      if (mdforc.eq.2) then
        write(iuo+66,*) 'Restoring not applied in some regions'
        write(iuo+66,*) 'file = ', filcor
        open(unit = 29, file = filcor, status = 'old')
        read(29,*) zrap
        do ns=0,nsmax
         do j=1,jmax
          do i=1,imax
           rappes(i,j,ns)=rappes(i,j,ns)*zrap(i,j)
          enddo
         enddo
        enddo
        write(iuo+66,*) 'restoring 10,10 = ', rappes(10,10,0),rappes(10,10,1)
     &               ,rappes(10,10,2)
       endif
 
c- prepare le calcul du Flux derive de l'Advec.H.Obs :
Cadh  nflag = 1
Cadh  call initflx(nflag, nn99)
 
c--Mise en place definitive des tableaux "rappel & rappes" (= Delta-t / tau) :
      cctsr = unstyr / (yeaday * 86400.)
      do 730 k=1,kmax
        cctsrk = dts(k) * cctsr
        do 730 j=1,jmax
         do 730 i=1,imax
          rappel(i,j,k) = tms(i,j,k) * cctsrk * rappel(i,j,k)
 730  continue
      cctsrk = dts(ks2) * cctsr
      do 735 ns=0,nsmax
       do 735 j=1,jmax
        do 735 i=1,imax
         rappes(i,j,ns) = cctsrk * rappes(i,j,ns)
 735  continue
 
c--Mise en place du tableau "rapint" (=Delta-t/tau) ; si pas d'Obs -> rapint=0
      do 740 ns=1,nsmax
       do 740 k=1,kmax
        cctsrk = dts(k) * cctsr
        do 740 n=1,nrap(k,ns)
         i = 1 + mod(ijrap(n,k,ns),imax)
         j = 1 + ijrap(n,k,ns)/imax
         rapint(n,k,ns) = cctsrk * rapint(n,k,ns)
     &                  * min(tms(i,j,k), (scalr(i,j,k,ns)-spvr))
 740  continue
 
c--Terre ou Pas d'Obs -> Rappes = 0
      if (nn99.eq.2) then
       write(99,'(A)') 'Points sans Observations en surface :'
       do 755 j=1,jmax
        do 755 i=1,imax
          xxmsq = nsmax * tms(i,j,ks2)
          do 750 ns=1,nsmax
            xxns = min(tms(i,j,ks2), (scalr(i,j,ks2,ns)-spvr))
            rappes(i,j,ns) = rappes(i,j,ns) * xxns
            xxmsq = xxmsq - xxns
 750      continue
          xxns = min(tms(i,j,ks2), (scalr(i,j,ks2,2)-spvr))
          rappes(i,j,0) = rappes(i,j,0) * xxns
          if (xxmsq.gt.epsil) write(99,'(2I4,3x,1P10E13.5)')
     &        i,j,(scalr(i,j,ks2,ns),ns=1,nsmax)
 755   continue
      else
        do 760 ns=0,nsmax
         nsrp = max(ns,2-ns)
         do 760 j=1,jmax
          do 760 i=1,imax
           rappes(i,j,ns) = rappes(i,j,ns)
     &              * min(tms(i,j,ks2), (scalr(i,j,ks2,nsrp)-spvr))
 760    continue
      endif
c--Pas d'Obs -> Rappel = 0
      do 770 ns=1,nsmax
       do 770 k=1,kmax
        do 770 j=1,jmax
         do 770 i=1,imax
          rappel(i,j,k) = min(rappel(i,j,k), (scalr(i,j,k,ns)-spvr))
 770  continue
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c--conversion & changements d'unites (mais pas de signe !), limiteurs de flux :
c- unitfx(0) = Year : Fx. en m/y ; ATTENTION : phiss(0) en m et phimnx en m/s.
c- unitfx(1) = rho.Cp : Fx. en W/m2 ; unitfx(2) = Year / 34.7 g/l : Fx. en m/y
c- conversion tableau phiss <- utilise directement ds equation des scalaires.
      do 790 ns=0,nsmax
        if (yforc(ns).ne.zero) then
          ccmult = yforc(ns) * dts(ks2) * unsdz(ks2)
          if (ns.eq.0) ccmult = yforc(ns) * dts(ks2)
          do 780 j=1,jmax
           do 780 i=1,imax
            phiss(i,j,ns) = ccmult * phiss(i,j,ns)
 780      continue
        endif
Cic0      ccmult = dts(ks2) * unsdz(ks2) / abs(unitfx(ns))
Cic0      if (ns.eq.0) ccmult = 1.0 / abs(unitfx(ns))
Cic0      do 785 j=1,jmax
Cic0       do 785 i=1,imax
Cic0        phimnx(i,j,0,ns) = ccmult * phimnx(i,j,0,ns)
Cic0        phimnx(i,j,1,ns) = ccmult * phimnx(i,j,1,ns)
 785      continue
 790  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  8 ) Impression de controle sur le fichier "mouchard" .              |
c-----------------------------------------------------------------------
 
      if(nn99.ne.2) goto 899
c--ecriture de controle :
      write(99,*) 'DeltaT(Surf)/tau , Coeff rappel Expl. FW, T, S,'
     &          //' rappel N, S :'
      write(99,*) cctsrk, (rapp0(k),k=0,nsmax+2)
      write(99,*) 'coeff rappel de 1 a kmax :'
      write(99,*) rapp1
 
 899  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  9 ) Sortie de la routine .
c-----------------------------------------------------------------------
 
      write(iuo+66,'(A)') 'Files read :'//line(:nnline)
 
      return
 
 900  continue
 
      write(iuo+66,*) 'Arret routine "redforc" :'
      write(iuo+66,*) 'Probleme de lecture apres les fichiers suivant :'
      write(iuo+66,*) line(:nnline)
 
      stop
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine redforc -
      end
