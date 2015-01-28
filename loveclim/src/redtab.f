












      subroutine redtab(ytab,nmax,numfil,kkr)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c Appele par "redforc", "redrunb".
c Lecture sur fich. numfil, en un seul bloc, des nmax elements du tableau ytab.
c  modif : 16/15/96
 
      include 'type.com'
 
c- dummy variables :
      dimension ytab(nmax)
 
      kkr = 1
      read(numfil,end=920,err=910) ytab
 
      return
 910  continue
      kkr = -1
      return
 920  continue
      kkr = 0
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine redtab -
      end
