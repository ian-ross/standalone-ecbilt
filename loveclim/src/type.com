












c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  fichier "type.com" : incorpore par instruction 'include' dans les programmes
c   (et les routines des programmes) :
c      OM, CLASS, FEEDOM, GRADP, STATIS, DIFBIN, TROMNC, UNIBIN, TRSBATH.
c contient les definitions de type (sauf les chaines de caracteres).
c  modif : 08/01/95
 
c--declaration implicite de type (standard fortran) :
      implicit double precision (a-h,o-z)
Cray  implicit real (a-h,o-z)
 
c--declaration explicite de type :
      complex*16 acplex, bcplex, vcplex
Cray  complex*8  acplex, bcplex, vcplex
 
c--declaration explicite PVM de type :
c      integer PVMDEFAULT/0/
c      integer INTEGER4/3/
c      integer REAL8/6/
 
c--fin du fichier "type.com"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
