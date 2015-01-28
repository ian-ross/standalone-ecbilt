












c
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  Thermo.com is incorparated by an instruction include in thersf.f ,
c      fontbc.f and acrlbq.f. It comprises the commons associated to
c      thermodynamic ice computation
c
C Correspondance between the variables
c qlbqb   qlbq
c qcmbqb  qcmbq
c thcmb   thcm
c fstbqb  fstrbq
c fltbqb  ffltbq
c fscbqb  fscmbq
c fsolgb  fsolg
c ratbqb  ratbqg
c psbqb   psbq
c tabqb   tabq
c qabqb   qabq
c vabqb   vabq
c qfvbqb  qfvbq
c tsb     ts
c tfub    tfu
c hnpbqb  zhnpbq
c hnbqb   hnbq
c hgbqb   hgbq
c albqb   albq
c qstbqb  qstobq
c fbbqb   fbbq
c tbqb    tbq
c dmgbqb  dmgbq
c dmnbqb  dmnbq
c qlbbqb  qlbsbq
c cldqb   cloud
c dmgwib  dmgwi
c npb     number of points where computations has to be done
c npac    correspondance between the points
c fratsb  firg
c fcsb    fcsg
c fleb    fleg
c dvsbqb  dvosbq
c dvbbqb  dvobbq
c dvlbqb  dvolbq
c dvnbqb  dvonbq
c hgcolb  hgcol (applies only when Cvhg is activated)
 
      common/combq/qlbqb(nbpt),qcmbqb(nbpt),thcmb(nbpt),fstbqb(nbpt),
     &             fltbqb(nbpt),fscbqb(nbpt),fsolgb(nbpt),
     &             ratbqb(nbpt),psbqb(nbpt),tabqb(nbpt),
     &             qabqb(nbpt),vabqb(nbpt),qfvbqb(nbpt),tsb(nbpt),
     &             tfub(nbpt),hnpbqb(nbpt),hnbqb(nbpt),hgbqb(nbpt),
     &             albqb(nbpt),qstbqb(nbpt),fbbqb(nbpt),tbqb(nbpt,nkb0),
     &             dmgbqb(nbpt),dmnbqb(nbpt),qlbbqb(nbpt),cldqb(nbpt),
     &             dmgwib(nbpt),
     &             hgcolb(nbpt),
     &             npb(nbpt),npac(nbpt)
c
      common/comdbq/fratsb(nbpt),fcsb(nbpt),
     &              fleb(nbpt),dvsbqb(nbpt),dvbbqb(nbpt),dvlbqb(nbpt),
     &              dvnbqb(nbpt)
c
Ccp2  common/fdericeb/fderb(nbpt)
c-- End of file 'thermo.com'
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
