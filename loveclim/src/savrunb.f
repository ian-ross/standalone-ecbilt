












      subroutine savrunb(nnt,nn99,ccfile)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  sortie des resultats sur fichier binaire
c     1ere partie utilisable pour faire redemarrer le programme ;
c  fichier de sortie 'res/n/.om' avec "n" decroissant jusqu'a zero
c-----
c  modif : 02/07/98
 
      include 'type.com'
      include 'const.com'
      include 'para.com'
      include 'bloc.com'
      include 'varno.com'
      include 'vareq.com'
      include 'ice.com'
      include 'dynami.com'
      include 'moment.com'
      include 'comunit.h'
 
c- dummy variables :
      character*(*) ccfile
c- variables locales equivalentes :
 
c- variables locales :
      dimension nnvv(nvmax)
      character*120 ccline
 
c--instructions "data" :
      include 'datadc.com'
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Preparation de l'Ecriture ; definition et ecriture de l'entete  |
c-----------------------------------------------------------------------
 
c- Unite(fortran) du Fichier "ccfile" a ecrire :
      nnfbin = 50
      nnct = len(ccline)
      ccline(:6) = 'Wbin :'
      nnc = 6
 
c- calcul du nb de tableaux a ecrire :
      nbvar = 0
      do 110 nv=1,nvmax
C       if (nvrl(nv).ne.0) write(iuo+66,'(A,2I4)') titcv(nv), nv, nvrl(nv)
        if (mod(nvrl(nv),2).eq.1) then
          nbvar = nbvar + 1
          nnvv(nbvar) = nv
        endif
 110  continue
 
      if (nbvar.eq.0) then
        write(iuo+66,*) 'savrunb : aucun tableau a ecrire sur fichier !'
        return
      endif
 
c--Ouverture et ecriture de l'en-tete
 
      open(nnfbin,file=ccfile,status='unknown',form='UNFORMATTED')
      write(nnfbin) numit, tpstot, refexp
      write(nnfbin) nbvar
 
      do 300 n=1,nbvar
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Ecriture sur fichier variable par variable .
c-----------------------------------------------------------------------
 
      nv = nnvv(n)
C     write(iuo+66,'(A,2I4)') titcv(nv), nv, nvrl(nv)
 
c- taille du tableau (=Nb d'elements) a ecrire :
      kksize = krlm(nv)
      if (ltyp(nv).ge.0) kksize = kksize * ixjmax
 
c- Ecriture du No="nv", de la taille et du nom (3cc) de la variable "nv"
      write(nnfbin) nv, kksize, titcv(nv)
 
c- Recherche du vrai tableau et ecriture par appel a la routine "savtab" :
      if (nv.le.nsmax) then
           call savtab(scal(1,1,1,nv),kksize,nnfbin,kkerr)
      elseif (nv.ge.nvrfw .and. nv.le.nvrfw+nsmax) then
           ns = nv-nvrfw
           call savtab(fss(1,1,ns),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvret) then
           call savtab(eta(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrub) then
           call savtab(ub(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrvb) then
           call savtab(vb(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrajc) then
           call savtab(fqajc(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrb) then
           call savtab(b(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvru) then
           call savtab(u(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrv) then
           call savtab(v(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrtke) then
           call savtab(q2turb(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrw) then
           call savtab(w(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrn2) then
           call savtab(bvf(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvras) then
           call savtab(avsdz(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrau) then
           call savtab(avudz(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrusl) then
           call savtab(uslpfx(1,1,-1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrvsl) then
           call savtab(vslpfx(1,1,-1),kksize,nnfbin,kkerr)
c-----
      elseif (nv.eq.nvrum) then
           call savtab(umoy(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrvm) then
           call savtab(vmoy(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrhg) then
           call savtab(hgbq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrfq) then
           call savtab(fsbbq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrqs) then
           call savtab(qstobq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvral) then
           call savtab(albq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrhn) then
           call savtab(hnbq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrts) then
           call savtab(ts(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrug) then
           call savtab(ug(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrvg) then
           call savtab(vg(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrtbq) then
           call savtab(tbq(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrxzo) then
           call savtab(xzo(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrtgx) then
           call savtab(tenagx(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrtgy) then
           call savtab(tenagy(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrmom) then
           call savtab(vicmom(1,1,1),kksize,nnfbin,kkerr)
c----
      else
        write(iuo+66,'(3A,I3,A)') 'ARRET, savrunb : Ecriture de ',
     &         titcv(nv), ' nv=', nv, ' Pas Prevue !'
        goto 900
      endif
      if (kkerr.eq.-1) then
        write(iuo+66,'(3A,I3,A)') 'ARRET, savrunb : Error writing ',
     &         titcv(nv), ' nv=', nv
        goto 900
      endif
 
c--Consigne les noms des variables ecrites sur "ccline" :
      nnc1 = 1 + nnc
      nnc = min(3+nnc1,nnct)
      if(nnc1.le.nnc) ccline(nnc1:nnc) = ' '//titcv(nv)
 
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--fin de l'ecriture de la nariable "nv".
 300  continue
 
      close(nnfbin)
 
      write(iuo+66,'(A,I11,A,I3,2A)') 'Exp '//refexp//', Iter', numit,
     &  ' :', nbvar, ' variables written  on  File ', ccfile
      if (nn99.eq.2) then
        write(99,'(A)') ccline(:nnc)
        write(99,'(A,I11,A,I3,2A)') 'Exp '//refexp//', Iter', numit,
     &  ' :', nbvar, ' variables written  on  File ', ccfile
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  7 ) Remise a zero des variables moyennnes et extremes .
c-----------------------------------------------------------------------
 
      if (nnt.ge.2 .and. numit.ge.nstart) then
        do 710 k=1,kmax
         do 710 j=1,jmax
          do 710 i=1,imax
           fqajc(i,j,k) = 0.0
 710    continue
        do 715 ns=0,nsmax
         do 715 j=1,jmax
          do 715 i=1,imax
           fss(i,j,ns) = 0.0
 715    continue
      endif
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  9 ) Traitement des cas "problematiques" .                           |
c-----------------------------------------------------------------------
 
 900  continue
      write(iuo+66,'(A,I11,A,I3,2A)') 'Exp '//refexp//', Iter', numit,
     &  ' :', n-1, ' var. written & STOP,   File ', ccfile
      write(iuo+66,'(A)') ccline(:nnc)
 
      if (nn99.eq.2) close(99)
      close(nnfbin)
 
      stop
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine savrunb -
      end
