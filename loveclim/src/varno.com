












c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c fichier "varno.com" : incorpore par instruction 'include' dans les routines :
c  CLASS, GRADP, STATIS, TROMNC, FLTROM,
c    defcst, savrunb, redrunb, savrunc, redrunc, foroutp,
c    ncdfout, rdclpar, process, scadhv, meridflu, moyen, local, binout
c  (commun a toutes les routines de sortie (output) de resultats)
c   inclus apres "type.com", "para.com", "bloc.com".
c------------------
c Chaque variable est reperee par un numero "nv" (nv > 0),
c  nv = -1,0 pour fichiers (et autre) regroupant les moyennes longitudinales
c  nv = -2 pour le fichier regroupant les bilans globaux et par niveaux.
c  nv de -10-nsmax a -10 pour transports meridiens de masse (-10) et de scalaire
c-----
c La localisation de la variable sur la maille est consignee
c  dans le tableau "ltyp", croissant dans le sens trigo., =99 => non defini
c     0 a 3 : verticalement centre ; 4 a 7 interface entre 2 niveaux ;
c     8 a 11 : sans localisation verticale ( eta , ub , vb ...)
c     + 12 si 1ere composantE vectorrielle, + 24 si 2nd composante.
c  titres et formats : lus dans le fichier "class.par"
c------------------
c  modif : 30/01/98
 
c--common definissant les correspondances (variables, niveaux associes) :
      common / indeks/
     & krl1(0:nvmax), krlm(0:nvmax), kvar2d(0:nvmax), ltyp(0:nvmax),
     & nvrl(0:nvmax)
 
      character*3 titcv
      common / cvname / titcv(-2:nvmax)
 
c--fin du fichier "varno.com"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
