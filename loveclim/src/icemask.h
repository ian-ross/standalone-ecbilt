!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  creation : 18/06/2007 <-- Mask of Precip, Lev Tarasov forcing !

!--blocs common :

      INTEGER, DIMENSION(nlat, nlon):: icemask
      REAL(KIND=4) :: totalprecipcases
      common /lev_tarasov_forcing/icemask,totalprecipcases
