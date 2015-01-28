












c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c fichier "datadc.com" (anciennement dans fichier "varno.com")
c  regroupant les instructions "DATA" ; inclus apres TOUTES les declarations !
c incorpore par instruction 'include' dans les routines :
c  CLASS, STATIS, FLTROM, defcst, savrunb, redrunb, savrunc, redrunc, foroutp,
c     ncdfout, rdclpar, process, meridflu, local, binout
c------------------
c  modif : 30/06/98
 
c--definition du numero de chaque variable (pour + de transparence) :
 
      data  nvrt   /  1 /
      data  nvrs   /  2 /
      data  nvrs3  /  3 /
      data  nvrs4  /  4 /
      data  nvrfw  / 10 /
      data  nvrfc  / 11 /
      data  nvrfs  / 12 /
      data  nvrfs3 / 13 /
      data  nvrfs4 / 14 /
 
      data  nvreac / 20 /
      data  nvret  / 21 /
      data  nvrub  / 22 /
      data  nvrvb  / 23 /
      data  nvrps  / 24 /
      data  nvrhac / 25 /
 
      data  nvrajc / 30 /
      data  nvrb   / 31 /
      data  nvru   / 32 /
      data  nvrv   / 33 /
      data  nvrtke / 34 /
      data  nvrw   / 35 /
      data  nvrn2  / 36 /
      data  nvras  / 37 /
      data  nvrau  / 38 /
      data  nvrusl / 40 /
      data  nvrvsl / 41 /
 
c- variables pour Sea Ice :
      data  nvrhg  / 50 /
      data  nvrfq  / 51 /
      data  nvrqs  / 52 /
      data  nvral  / 53 /
      data  nvrhn  / 54 /
      data  nvrts  / 55 /
 
c- variables pour Bilan detaille :
      data  nvraxt / 60 /
      data  nvrayt / 61 /
      data  nvrhat / 62 /
      data  nvrhdt / 63 /
      data  nvrvat / 64 /
      data  nvrvdt / 65 /
      data  nvrvaf / 66 /
      data  nvrvdf / 67 /
 
c- variables pour Sea Ice :
      data  nvrum  / 70 /
      data  nvrvm  / 71 /
      data  nvrug  / 72 /
      data  nvrvg  / 73 /
      data  nvrtbq / 74 /
      data  nvrxzo / 75 /
      data  nvrtgx / 76 /
      data  nvrtgy / 77 /
      data  nvrmom / 80 /
 
c- variables pour diffus. Isopyc. ou G.&M.W. :
      data  nvravi / 90 /
      data  nvravs / 91 /
      data  nvrslx / 92 /
      data  nvrsly / 93 /
      data  nvrpsx / 94 /
      data  nvrpsy / 95 /
      data  nvrugm / 96 /
      data  nvrvgm / 97 /
      data  nvrwgm / 98 /
 
c--fin du fichier "datadc.com"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
