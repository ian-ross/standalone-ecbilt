












c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  fichier "para.com" : incorpore par instruction 'include' dans les programmes
c   (et les routines des programmes) :
c      OM, CLASS, FEEDOM, GRADP, STATIS, DIFBIN, TROMNC, UNIBIN, TRSBATH.
c contient la plupart des parametres definissant les tableaux courrants :
c--ocean mondial, 3x3 deg, version 2 Sous-grilles reunies en une seule.
c-----
c Ctke [Ctk0] => ligne specifique a la version avec [sans] TKE .
c Cice [Cic0] => ligne specifique a la version avec [sans] glace marine .
c-----
c  modif : 01/05/98
c
c Reorganisation : 1/03/2004; A. Mouchet pour couplage avec Loch
C -> cree para0.com inclus dans para.com et dans declarations Loch
c
C-parametres lies a la taille du domaine; freq. des donnees; etc.
C


      INCLUDE 'para0.com'
 
c--parametres donnant le nombre de niveaux dans la glace
      integer nkb0
      parameter ( nkb0 = 3 )
 
c--parametres fixant le nombre de zones (= nb bassins), nb de detroits :
c [indispensable pour inclusion du fichier  reper.com]
      integer nbsmax,nhsfmx,nltmax
      parameter ( nbsmax = 3 , nhsfmx = 10 )
      parameter ( nltmax = 4 )
c--parametres pour les sorties sur fichier "evolu" par la routine "informe" :
c [indispensable pour inclusion du fichier  reper.com]
      integer ninfmx,nchinf,nchsep
      parameter ( ninfmx = 30 + nhsfmx + 2*kmax + 20)
      parameter ( nchinf = 5 , nchsep = nchinf+2 )
 
c--parametres pour couplage : Nombre de tableaux envoyes et recus par l'ocean :
      integer ntocn,ntatm
      parameter ( ntocn = 4 , ntatm = 10 )
 
c--parametre indiquant le rang "k" (ds le tableau general),
c   du 1er niveau occupe par la variable:
c [indispensable pour inclusion du fichier  vareq.com]
      integer krlu,krlfw,krlfc,krlfs,krlfs3,krlfs4,krlps,krlet,krlub,
     &        krlvb,krlv,krlt,krls,krls3,krls4,krlb,krln2,krlas,krlau,
     &        krlw,krltke,krlajc
      parameter( krlu  =  0 )
      parameter( krlfw = krlu - nsmax - 5 )
      parameter( krlfc = krlfw + 1 )
      parameter( krlfs = krlfw + 2 )
      parameter( krlfs3= krlfw + 3 )
      parameter( krlfs4= krlfw + 4 )
      parameter( krlps = krlu - 4 )
      parameter( krlet = krlu - 3 )
      parameter( krlub = krlu - 2 )
      parameter( krlvb = krlu - 1 )
      parameter( krlv  = krlu + kmax   )
      parameter( krlt  = krlu + kmax*2 )
      parameter( krls  = krlu + kmax*3 )
      parameter( krls3 = krlu + kmax*4 )
      parameter( krls4 = krlu + kmax*5 )
      parameter( krlb  = krlu + kmax*(2+nsmax) )
      parameter( krln2 = krlu + kmax*(3+nsmax) )
      parameter( krlas = krlu + kmax*(4+nsmax) )
      parameter( krlau = krlu + kmax*(5+nsmax) )
      parameter( krlw  = krlu + kmax*(6+nsmax) )
      parameter( krltke= krlu + kmax*(7+nsmax) + 1 )
      parameter( krlajc= krlu + kmax*(8+nsmax) + 2 )
 
c  Avu=>Avi S=>Psx T=>Psy U=>Ugm V=>Vgm W=>Wgm B=>Slx N2=>Sly
      integer krlavi,krlavs,krlslx,krlsly,krlpsx,krlpsy,krlugm,krlvgm,
     &   krlwgm,krlusl,krlvsl,krlhac,krleac,krlhg,krlfq,krlqs,krlal,
     &   krlhn,krlts,krlug,krlvg,kvsice,krlvaf,krlvdf,krlaxt,krlayt,
     &   krlhat,krlhdt,krlvat,krlvdt
      parameter( krlavi=krlau)
      parameter( krlavs=krlas)
      parameter( krlslx=krlb)
      parameter( krlsly=krln2)
      parameter( krlpsx=krlt)
      parameter( krlpsy=krls)
      parameter( krlugm=krlu)
      parameter( krlvgm=krlv)
      parameter( krlwgm=krlw)
 
      parameter( krlusl= krlajc + kmax )
      parameter( krlvsl= krlusl + nsmax   + 2 )
      parameter( krlhac= krlusl + nsmax*2 + 4 )
      parameter( krleac= krlhac + 1 )
 
Cic0  parameter( krlhg = krlajc )
      parameter( krlhg = krleac + 1 )
      parameter( krlfq = krlhg  + 1 )
      parameter( krlqs = krlhg  + 2 )
      parameter( krlal = krlhg  + 3 )
      parameter( krlhn = krlhg  + 4 )
      parameter( krlts = krlhg  + 5 )
      parameter( krlug = krlhg  + 6 )
      parameter( krlvg = krlhg  + 7 )
Cic0  parameter( kvsice = 0 )
      parameter( kvsice = 8 )
 
      parameter( krlvaf= krlt )
      parameter( krlvdf= krlvaf + kmax )
      parameter( krlaxt= krlvaf + kmax*2 )
      parameter( krlayt= krlvaf + kmax*3 )
      parameter( krlhat= krlvaf + kmax*4 )
      parameter( krlhdt= krlvaf + kmax*5 )
      parameter( krlvat= krlvaf + kmax*6 )
      parameter( krlvdt= krlvaf + kmax*7 + 1)
 
c--parametres lies a la definition du tableau general utilise pour les sorties:
c [indispensable pour inclusion du fichier  var??.com]
      integer ltymax,krlmin,nvmax,nv3dmx,nv2dmx,kv2dmx,krlmax
      parameter( ltymax = 11 )
      parameter( krlmin = -5 - nsmax)
      parameter( nvmax = 99 , nv3dmx = 9+nsmax )
c- nv2dmx = Nb. Var. 2D (fixe+Ns+Ice) ; kv2dmx = Nb. Niv. reserves pour var. 2D
      parameter( nv2dmx = 9+nsmax+15 , kv2dmx = 11 + 3*nsmax + kvsice )
      parameter( krlmax = krlmin + 1 + nv3dmx*kmax + kv2dmx )
 
c--fin du fichier "para.com"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
