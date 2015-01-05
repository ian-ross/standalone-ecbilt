c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** open statements of atmosphere input files:
c *** units 90-99 are reserved for initial states
c *** units 20 and higher are reserved for other parts of ecbilt
c-----------------------------------------------------------------------

c *** initialisation data

      open(iuo+1,file='inputdata/coef.dat',
     &     status='old',form='unformatted')
c      open(iuo+3,file='inputdata/berg.dat',
c     &     status='old',form='unformatted')
      open(iuo+7,file='inputdata/gauss.asc',
     &     status='old',form='formatted')
      open(iuo+11,file='inputdata/win.dat',
     &                             status='old',form='formatted')
      open(iuo+12,file='inputdata/sum.dat',
     &                   status='old',form='formatted')

      open(iuo+14,file='inputdata/ocbas.dat',
     &                   status='old',form='formatted')

      open(iuo+16,file='inputdata/lwrref.dat'
     &      ,form='unformatted')

      open(iuo+17,file='inputdata/lwrcoef.dat'
     &      ,form='unformatted')

      open(iuo+18,file='inputdata/swrref.dat'
     &      ,form='unformatted')

      open(iuo+19,file='inputdata/swrcoef.dat'
     &      ,form='unformatted')

      open(iuo+15,file = 'namelistecbilt',status='old',form='formatted')

      open(iuo+33,file='inputdata/GHG.dat')
      open(iuo+34,file='inputdata/TSI.dat')
      open(iuo+35,file='inputdata/VOLC.dat')
      open(iuo+36,file='inputdata/SUL.dat',form='unformatted')
      open(iuo+37,file='inputdata/OZONE.dat')

      open(iuo+38,file='inputdata/mbcs2_cor')
      open(iuo+39,file='inputdata/scenario2Xco2.dat')
