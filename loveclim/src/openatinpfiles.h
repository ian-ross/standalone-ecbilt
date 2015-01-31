!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** open statements of atmosphere input files:
! *** units 90-99 are reserved for initial states
! *** units 20 and higher are reserved for other parts of ecbilt
!-----------------------------------------------------------------------

! *** initialisation data

      open(iuo+1,file='inputdata/coef.dat',status='old',form='unformatted')
      open(iuo+7,file='inputdata/gauss.asc',status='old',form='formatted')

      open(iuo+14,file='inputdata/ocbas.dat',status='old',form='formatted')

      open(iuo+16,file='inputdata/lwrref.dat',form='unformatted')
      open(iuo+17,file='inputdata/lwrcoef.dat',form='unformatted')

      open(iuo+18,file='inputdata/swrref.dat',form='unformatted')
      open(iuo+19,file='inputdata/swrcoef.dat',form='unformatted')

      open(iuo+15,file = 'namelistecbilt',status='old',form='formatted')

      open(iuo+38,file='inputdata/mbcs2_cor')
