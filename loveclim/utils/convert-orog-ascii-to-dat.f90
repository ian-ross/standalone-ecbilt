PROGRAM convert_orog_to_dat
  IMPLICIT NONE
  INCLUDE 'comatm.h'
  INCLUDE 'comdyn.h'

  REAL*8 agg1(nlat,nlon), agg2(nlat,nlon)

  OPEN(unit=3,file='berg.dat',form='unformatted')
  OPEN(unit=13,file='berg.asc',form='formatted')
  READ (13,*) agg1
  READ (13,*) agg2
  WRITE (3) agg1
  WRITE (3) agg2
END PROGRAM convert_orog_to_dat
