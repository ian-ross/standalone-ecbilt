












c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  bloc "trace.com" : incorpore par instruction 'include' dans les routines
c      liees aux traceurs : cfc_flx.f forcng.f
c  modif : 12/04/99
 
c--blocs common :
c   Constants needed for solubility computations of CFCs in seawater
c       from Warner and Weiss (1985, Deep Sea Res., equation 6).
c   Note that F is in mol.kg-1.atm-1 following table 5 of WW(1985).
c==========================================================================

 
      parameter (a1cfc11=-232.0411, a2cfc11=322.5546, a3cfc11=120.4956,
     &           a4cfc11=-1.39165,  b1cfc11=-.146531, b2cfc11=0.093621,
     &           b3cfc11=-.0160693)
 
      parameter (a1cfc12=-220.2120, a2cfc12=301.8695, a3cfc12=114.8533,
     &           a4cfc12=-1.39165,  b1cfc12=-.147718, b2cfc12=0.093175,
     &           b3cfc12=-.0157340)
 
c       cfc11atm = variable holding time-interpolated SH 0-11
c       cfc12atm = variable holding time-interpolated SH 0-12
 
        common /cfcmhe1/ cfc11(66,2), cfc12(66,2),
     &                   atmcfc11(2), atmcfc12(2)
 
c--fin du fichier "trace.com"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
