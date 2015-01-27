c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** open statements of atmosphere input files:
c *** units 90-99 are reserved for initial states
c *** units 20 and higher are reserved for other parts of ecbilt
c-----------------------------------------------------------------------

c *** initialisation data

      open(unit=1,file='inputdata/atmos/coef.dat',
     &     status='old',form='unformatted')
      open(unit=3,file='inputdata/atmos/berg.dat',
     &     status='old',form='unformatted')
      open(unit=4,file='inputdata/atmos/albedo.dat',
     &     status='old',form='formatted')
      open(unit=7,file='inputdata/atmos/gauss.dat',
     &     status='old',form='formatted')
      open(unit=8,file='inputdata/atmos/cloudpar.dat',
     &     status='old',form='formatted')
      open(unit=9,file='inputdata/atmos/lwavepar.dat',
     &     status='old',form='formatted')

      open(unit=10,file='inputdata/land/labas.dat')

      open(13,file='inputdata/atmos/cdrag.dat',
     &     status='old',form='formatted')


      open(15,file = 'namelist',status='old',form='formatted')
