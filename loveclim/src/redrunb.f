












      subroutine redrunb(nnt,nn99,ccfile)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  redemarrage a partir de l'etat definit par le fichier binaire "rest.om" .
c Entree : nnt = 0,1  => lecture du fichier de resultat Nouveau(>05/96) format
c          nnt > 3    => lecture du fichier de resultat Ancien( <05/96) format
c              = 4,5  => lit  ancien fichier de resultat simple.
c              = 6,7  => lit  ancien fichier de resultat complet (sans w).
c              = 8,9  => lit +ancien(<01/96) fich. res.  complet (avec w).
c----- nnt(Sortie) = nnt(Entree)
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
 
      parameter (nszmax = ijkmax)
 
c- dummy variables :
      character*(*) ccfile
c- variables locales equivalentes :
      dimension vloc(nszmax)
      equivalence ( q(1,1,1), vloc(1) )
c- pour ancien fichier :
      dimension wloc(imax,jmax,kmax)
      equivalence ( w(1,1,1) ,  wloc(1,1,1) )
C     dimension egajc(imax,jmax)
C     equivalence ( fqajc(1,1,1), egajc(1,1) )
C     dimension h1ajc(imax,jmax)
C     equivalence ( fqajc(1,1,kmax), h1ajc(1,1) )
 
c- local varaibles :
      character*3 cc3
      character*6 cc6exp
      character*120 ccline
 
c--instructions "data" :
      include 'datadc.com'
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Initialisation & Ouverture du fichier binaire a lire .          |
c-----------------------------------------------------------------------
 
c- Unite(fortran) du Fichier "ccfile" a lire :
      nnfbin = 60
 
      nnct = len(ccline)
      ccline(:6) = 'Rbin :'
      nnc = 6
 
      open(nnfbin,file=ccfile,status='old',form='UNFORMATTED')
 
      if (nnt.ge.4) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Lecture, fichier binaire "ccfile", Ancien format                |
c-----------------------------------------------------------------------
 
c--lecture de la premiere partie :
      read (nnfbin, end=281) numit
      read (nnfbin, end=281) tpstot
      read (nnfbin, end=281) eta
      read (nnfbin, end=281) ub
      read (nnfbin, end=281) vb
      read (nnfbin, end=281) u
      read (nnfbin, end=281) v
      read (nnfbin, end=281) scal
 
      read (nnfbin, end=281) umoy
      read (nnfbin, end=281) vmoy
      read (nnfbin, end=281) hgbq
      read (nnfbin, end=281) fsbbq
      read (nnfbin, end=281) qstobq
      read (nnfbin, end=281) albq
      read (nnfbin, end=281) hnbq
      read (nnfbin, end=281) ts
      read (nnfbin, end=281) ug
      read (nnfbin, end=281) vg
      read (nnfbin, end=281) tbq
      read (nnfbin, end=281) xzo
      read (nnfbin, end=281) tenagx
      read (nnfbin, end=281) tenagy
      read (nnfbin, end=281) sxg
      read (nnfbin, end=281) syg
      read (nnfbin, end=281) sxxg
      read (nnfbin, end=281) syyg
      read (nnfbin, end=281) sxyg
      read (nnfbin, end=281) sxn
      read (nnfbin, end=281) syn
      read (nnfbin, end=281) sxxn
      read (nnfbin, end=281) syyn
      read (nnfbin, end=281) sxyn
      read (nnfbin, end=281) sxa
      read (nnfbin, end=281) sya
      read (nnfbin, end=281) sxxa
      read (nnfbin, end=281) syya
      read (nnfbin, end=281) sxya
      read (nnfbin, end=281) sxc0
      read (nnfbin, end=281) syc0
      read (nnfbin, end=281) sxxc0
      read (nnfbin, end=281) syyc0
      read (nnfbin, end=281) sxyc0
      read (nnfbin, end=281) sxc1
      read (nnfbin, end=281) syc1
      read (nnfbin, end=281) sxxc1
      read (nnfbin, end=281) syyc1
      read (nnfbin, end=281) sxyc1
      read (nnfbin, end=281) sxc2
      read (nnfbin, end=281) syc2
      read (nnfbin, end=281) sxxc2
      read (nnfbin, end=281) syyc2
      read (nnfbin, end=281) sxyc2
      read (nnfbin, end=281) sxst
      read (nnfbin, end=281) syst
      read (nnfbin, end=281) sxxst
      read (nnfbin, end=281) syyst
      read (nnfbin, end=281) sxyst
 
 
c-----
      nvrl(nvret) = mod(nvrl(nvret),2) + 2
      nvrl(nvrub) = mod(nvrl(nvrub),2) + 2
      nvrl(nvrvb) = mod(nvrl(nvrvb),2) + 2
      nvrl(nvru)  = mod(nvrl(nvru),2)  + 2
      nvrl(nvrv)  = mod(nvrl(nvrv),2)  + 2
      nvrl(nvrt)  = mod(nvrl(nvrt),2)  + 2
      nvrl(nvrs)  = mod(nvrl(nvrs),2)  + 2
 
      nvrl(nvrum)  = mod(nvrl(nvrum),2)  + 2
      nvrl(nvrvm)  = mod(nvrl(nvrvm),2)  + 2
      nvrl(nvrhg)  = mod(nvrl(nvrhg),2)  + 2
      nvrl(nvrfq)  = mod(nvrl(nvrfq),2)  + 2
      nvrl(nvrqs)  = mod(nvrl(nvrqs),2)  + 2
      nvrl(nvrhn)  = mod(nvrl(nvrhn),2)  + 2
      nvrl(nvral)  = mod(nvrl(nvral),2)  + 2
      nvrl(nvrts)  = mod(nvrl(nvrts),2)  + 2
      nvrl(nvrug)  = mod(nvrl(nvrug),2)  + 2
      nvrl(nvrvg)  = mod(nvrl(nvrvg),2)  + 2
      nvrl(nvrtbq) = mod(nvrl(nvrtbq),2) + 2
      nvrl(nvrxzo) = mod(nvrl(nvrxzo),2) + 2
      nvrl(nvrtgx) = mod(nvrl(nvrtgx),2) + 2
      nvrl(nvrtgy) = mod(nvrl(nvrtgy),2) + 2
      nvrl(nvrmom) = mod(nvrl(nvrmom),2) + 2
 
      if (nnt.eq.6 .or.nnt.eq.7) then
c--lecture de la seconde partie :
        read (nnfbin, end=282) cc6exp
        read (nnfbin, end=282) b
        read (nnfbin, end=282) bvf
        read (nnfbin, end=282) avsdz
        read (nnfbin, end=282) avudz
        read (nnfbin, end=282) fqajc
        read (nnfbin, end=282)
     &     (((fss(i,j,ns),i=1,imax),j=1,jmax),ns=1,nsmax)
c-----
        nvrl(nvrb)   = mod(nvrl(nvrb),2)   + 2
        nvrl(nvrn2)  = mod(nvrl(nvrn2),2)  + 2
        nvrl(nvras)  = mod(nvrl(nvras),2)  + 2
        nvrl(nvrau)  = mod(nvrl(nvrau),2)  + 2
        nvrl(nvrajc) = mod(nvrl(nvrajc),2) + 2
        nvrl(nvrfc)  = mod(nvrl(nvrfc),2)  + 2
        nvrl(nvrfs)  = mod(nvrl(nvrfs),2)  + 2
c--fin de la seconde partie .
        write(iuo+66,'(A,I11,3A)') 'Exp '//cc6exp//', Iter', numit,
     &      ' , File ', ccfile, ' : FULL RESULTS read - OK.'
 
      elseif (nnt.ge.8) then
c--lecture de la seconde partie, ancien fichier :
        read (nnfbin, end=282) cc6exp
        read (nnfbin, end=282) b
        read (nnfbin, end=282) bvf
        read (nnfbin, end=282) avsdz
        read (nnfbin, end=282) avudz
        read (nnfbin, end=282) wloc
c-----
        nvrl(nvrb)   = mod(nvrl(nvrb),2)   + 2
        nvrl(nvrn2)  = mod(nvrl(nvrn2),2)  + 2
        nvrl(nvras)  = mod(nvrl(nvras),2)  + 2
        nvrl(nvrau)  = mod(nvrl(nvrau),2)  + 2
        nvrl(nvrw)   = mod(nvrl(nvrw),2)   + 2
        read (nnfbin, end=286) egajc
        nvrl(nvreac) = mod(nvrl(nvreac),2) + 2
        read (nnfbin, end=286) hmajc
        nvrl(nvrhac) = mod(nvrl(nvrhac),2) + 2
        read (nnfbin, end=287)
     &     (((fss(i,j,ns),i=1,imax),j=1,jmax),ns=1,nsmax)
        nvrl(nvrfc)  = mod(nvrl(nvrfc),2)  + 2
        nvrl(nvrfs)  = mod(nvrl(nvrfs),2)  + 2
c--fin de la seconde partie .
        write(iuo+66,'(A,I11,3A)') 'Exp '//cc6exp//', Iter', numit,
     &      ' , File ', ccfile, ' : OLD FILE read I & II - OK.'
 
      else
        read (nnfbin,end=285) cc6exp
        write(iuo+66,'(A,I11,3A)') 'Exp '//cc6exp//', Iter', numit,
     &      ' , File ', ccfile, ' : Restart (I)  read - OK.'
      endif
 
      goto 290
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Cas de lecture Incomplete :
 281  continue
        write(iuo+66,'(A)') ' 1ere partie du fichier Incomplete ! =>  '
        write(iuo+66,'(3A)') 'ARRET dans "redrunb", fichier : ', ccfile
        STOP
 
 282  continue
        write(iuo+66,'(A)') ' 2nd  partie du fichier Incomplete ! =>  '
        write(iuo+66,'(3A)') 'ARRET dans "redrunb", fichier : ', ccfile
        STOP
 
 285  continue
        write(iuo+66,'(3A,I11,A)') 'Lecture de ', ccfile,
     &       ' terminee ( sans Ref., Iter', numit, ' ) .'
        goto 500
 
 286  continue
        write(iuo+66,'(A,I11,3A)') 'Exp '//cc6exp//', Iter', numit,
     &      ' , File ', ccfile, ' : read I & II except E/Hajc & fss!'
        goto 290
 
 287  continue
        write(iuo+66,'(A,I11,3A)') 'Exp '//cc6exp//', Iter', numit,
     &      ' , File ', ccfile, ' : read I & II except fss !'
 
 290  continue
 
      else
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Lecture, fichier binaire "ccfile", Nouveau format               |
c-----------------------------------------------------------------------
 
c--Lecture de l'en-tete
      read(nnfbin) numit, tpstot, cc6exp
      read(nnfbin) nbvar
 
      if (nbvar.gt.nvmax) then
        write(iuo+66,'(2(A,I4))') 'ARRET, redrunb : Nb.var=', nbvar,
     &     ' > nvmax=', nvmax
        goto 950
      endif
 
      llost = 0
      nnvrd = 0
      do 350 n=1,nbvar
c--Lecture variable par variable .
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--3.1 Lecture du No,taille & nom de la variable & verifications :
c-------------
      read(nnfbin,end=910,err=920) nv, kkr, cc3
 
c--Verification :
      if (nv.lt.1 .or. nv.gt.nvmax) then
        write(iuo+66,'(2(A,I3))') 'ARRET, redrunb : Var. nv= ', nv,
     &     ' out of range 1,nvmax=', nvmax
        goto 950
      elseif (titcv(nv).ne.cc3) then
        write(iuo+66,'(A,I3,4A)') 'ARRET, redrunb : nv=', nv,
     &        ', Var.= ',titcv(nv), ' <-> on File= ', cc3
        goto 950
      endif
c- taille (supposee) du tableau (=Nb d'elements) :
      kksize = krlm(nv)
      if (ltyp(nv).ge.0) kksize = kksize * ixjmax
      if (kksize.gt.kkr) then
        write(iuo+66,'(A,I3,2A,2(A,I8))') 'WARNING, redrunb : nv=', nv,
     &        ', Var.= ', titcv(nv), ', Nb.Elm. to read', kksize,
     &        ' > on file', kkr
        kksize = kkr
      elseif (kksize.eq.0) then
        write(iuo+66,'(3A,I3,A)') 'WARNING, redrunb : Lecture de ',
     &         titcv(nv), ' nv=', nv, ' Pas Prevue (Tab.Vide) !'
        llost = llost + 1
        kksize = 1
        call redtab(vloc(1),kksize,nnfbin,kkerr)
        vloc(1) = 0.
        nnr = mod(nvrl(nv),2)
        goto 320
      elseif (kksize.lt.kkr) then
        write(iuo+66,'(A,I3,2A,2(A,I8))') 'WARNING, redrunb : nv=', nv,
     &        ', Var.= ', titcv(nv), ', Nb.Elm.lu=', kksize,
     &        ' < on file=', kkr
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--3.2 Recherche du vrai tableau ; lecture par appel a "redtab" :
c-------------
      nnr = mod(nvrl(nv),2) + 2
      if (nv.le.nsmax) then
           call redtab(scal(1,1,1,nv),kksize,nnfbin,kkerr)
      elseif (nv.ge.nvrfw .and. nv.le.nvrfw+nsmax) then
           ns = nv-nvrfw
           call redtab(fss(1,1,ns),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvret) then
           call redtab(eta(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrub) then
           call redtab(ub(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrvb) then
           call redtab(vb(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrajc) then
           call redtab(fqajc(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrb) then
           call redtab(b(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvru) then
           call redtab(u(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrv) then
           call redtab(v(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrtke) then
           call redtab(q2turb(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrw) then
           call redtab(w(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrn2) then
           call redtab(bvf(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvras) then
           call redtab(avsdz(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrau) then
           call redtab(avudz(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrusl) then
           call redtab(uslpfx(1,1,-1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrvsl) then
           call redtab(vslpfx(1,1,-1),kksize,nnfbin,kkerr)
c-----
      elseif (nv.eq.nvrum) then
           call redtab(umoy(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrvm) then
           call redtab(vmoy(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrhg) then
           call redtab(hgbq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrfq) then
           call redtab(fsbbq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrqs) then
           call redtab(qstobq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvral) then
           call redtab(albq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrhn) then
           call redtab(hnbq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrts) then
           call redtab(ts(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrug) then
           call redtab(ug(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrvg) then
           call redtab(vg(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrtbq) then
           call redtab(tbq(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrxzo) then
           call redtab(xzo(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrtgx) then
           call redtab(tenagx(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrtgy) then
           call redtab(tenagy(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrmom) then
           call redtab(vicmom(1,1,1),kksize,nnfbin,kkerr)
c-----
      else
c- cas d'une variable dont la lecture n'est pas prevue :
        kksize = 1
        call redtab(vloc(1),kksize,nnfbin,kkerr)
c- reinitialise a zero "vloc" :
        vloc(1) = 0.
        write(iuo+66,'(3A,I3,A)') 'WARNING, redrunb : Lecture de ',
     &         titcv(nv), ' nv=', nv, ' Pas Prevue !'
        nnr = nnr - 2
        llost = llost + 1
c-----
      endif
 320  continue
 
      if (kkerr.eq.-1) then
        write(iuo+66,'(3A,I3,A)') 'ARRET, redrunb : Error reading ',
     &         titcv(nv), ' nv=', nv
        goto 950
      elseif (kkerr.eq.0) then
        write(iuo+66,'(3A,I3,A)') 'ARRET, redrunb : End_of_File reading ',
     &         titcv(nv), ' nv=', nv
        goto 950
      endif
c- fin de la lecture de la "n"ieme variable .
      if (ltyp(nv).ne.99) nvrl(nv) = nnr
 
      if (nnr.ge.2) then
c--Consigne les noms des variables luees sur "ccline" :
        nnvrd = nnvrd + 1
        nnc1 = 1 + nnc
        nnc = min(3+nnc1,nnct)
        if(nnc1.le.nnc) ccline(nnc1:nnc) = ' '//titcv(nv)
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 350  continue
c--fin de la lecture des resultats.
 
      write(iuo+66,'(A,I11,3A,I3,A)') 'Exp '//cc6exp//', Iter', numit,
     &    ' , File ', ccfile, ' :',  nnvrd, ' variables read - OK.'
      if (nn99.eq.2) then
        write(99,'(A)') ccline(:nnc)
        write(99,'(A,I11,3A,I3,A)') 'Exp '//cc6exp//', Iter', numit,
     &    ' , File ', ccfile, ' :',  nnvrd, ' variables read - OK.'
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
      if (refexp.eq.'      ') refexp = cc6exp
 
 500  continue
c--Fin de la lecture.
      close(nnfbin)
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  7 ) Traitement des cas d'un changement de stockage / bathymetrie    |
c-----------------------------------------------------------------------
 
C     if(nnt.lt.0) then
C       do 710 j=1,jmax
C        do 710 i=1,imax
C         ub(i,j) = ub(i,j) * hu(i,j)
C         vb(i,j) = vb(i,j) * hu(i,j)
C710    continue
C     endif
 
c--Tester si la bathymetrie a ete lue :
      if (unsvol.gt.zero) then
c--precaution : vitesse nulle en dehors du domaine .
        do 730 k=1,kmax
         do 730 j=1,jmax
          do 730 i=1,imax
            u(i,j,k) = u(i,j,k) * tmu(i,j,k)
            v(i,j,k) = v(i,j,k) * tmu(i,j,k)
 730    continue
        do 740 j=1,jmax
         do 740 i=1,imax
           ub(i,j) = ub(i,j) * tmu(i,j,ks2)
           vb(i,j) = vb(i,j) * tmu(i,j,ks2)
 740    continue
      endif
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  9 ) Traitement des cas "problematiques" .                           |
c-----------------------------------------------------------------------
 
 910  continue
      write(iuo+66,'(A)') 'ARRET, redrunb : End of File !'
      goto 950
 
 920  continue
      write(iuo+66,'(A)') 'ARRET, redrunb : Read Error  !'
 
 950  continue
      write(iuo+66,'(A,I11,3A,I3,A)') 'Exp '//cc6exp//', Iter', numit,
     &    ' , File ', ccfile, ' :',  nnvrd, ' variables read & STOP'
      write(iuo+66,'(A)') ccline(:nnc)
 
      if (nn99.eq.2) close(99)
      close(nnfbin)
 
      stop
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine redrunb -
      end
