PROGRAM gencoef
  IMPLICIT NONE
  INCLUDE 'comatm.h'
  INCLUDE 'comdyn.h'

  OPEN(unit=1,file='../fixed-data/coef.dat',form='unformatted')
  OPEN(unit=11,file='../fixed-data/coef.asc',form='formatted')
  READ (11,*) nshm, ll
  READ (11,*) pp
  READ (11,*) pd
  READ (11,*) pw
  WRITE (1) nshm, ll
  WRITE (1) pp
  WRITE (1) pd
  WRITE (1) pw
END PROGRAM gencoef
