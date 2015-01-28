












c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  creation : 18/06/2007 <-- Mask of Precip, Lev Tarasov forcing !

c--blocs common :

      INTEGER, DIMENSION(nlat, nlon):: icemask 
      REAL(KIND=4) :: totalprecipcases
      common /lev_tarasov_forcing/icemask,totalprecipcases

