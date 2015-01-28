












c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  fichier "const.com", incorpore par instruction 'include' dans les programmes
c   (et les routines des programmes) :
c      OM, CLASS, FEEDOM, GRADP, STATIS, DIFBIN, TROMNC, UNIBIN, TRSBATH.
c  modif : 06/02/98
c  modif : 01/03/04 <- declarations explicites pour inclusion dans LOCH
 
      DOUBLE PRECISION cpo, cstmin, cstmax, epsil, gpes,
     &                 omega, one, pi, radian, rho0, rterre, svrdrp,
     &                 unsrt, untour, yeaday, zero

c--blocs common :
 
      common / cstfix / cpo, cstmin, cstmax, epsil, gpes,
     &                  omega, one, pi, radian, rho0, rterre, svrdrp,
     &                  unsrt, untour, yeaday, zero

 
c--fin du fichier "const.com"
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
