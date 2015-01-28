












      subroutine savtab(ytab,nmax,numfil,kkw)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c Appele par "savrunb".
c Ecriture sur fich. numfil, en un seul bloc, des nmax elements du tableau ytab.
c  modif : 16/15/96
 
      include 'type.com'
 
c- dummy variables :
      dimension ytab(nmax)
 
      kkw = 1
      write(numfil,err=910) ytab
 
      return
 
 910  continue
      kkw = -1
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine savtab -
      end
