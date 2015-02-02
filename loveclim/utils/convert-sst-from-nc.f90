PROGRAM convert_sst

  ! Convert ECBILT SST forcing files from NetCDF to binary.

  USE netcdf

  IMPLICIT NONE

  INTEGER, PARAMETER :: NLON = 64, NLAT = 32
  INTEGER :: ncid, sstvarid
  INTEGER :: iday, ilat, ilon
  REAL*4, DIMENSION(NLON, NLAT, 360) :: sst

  OPEN (14, FILE='sst_daily.dat', FORM='unformatted')
  CALL chk(NF90_OPEN('ecbilt-sst.nc', NF90_NOWRITE, ncid))
  CALL chk(NF90_INQ_VARID(ncid, 'sst', sstvarid))
  CALL chk(NF90_GET_VAR(ncid, sstvarid, sst))
  DO iday = 1, 360
     WRITE (14) ((sst(ilon, ilat, iday), ilon = 1, NLON), ilat = 1, NLAT)
  END DO
  CALL chk(NF90_CLOSE(ncid))

CONTAINS
  SUBROUTINE chk(ncstatus)
    USE netcdf
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ncstatus

    IF (ncstatus /= NF90_NOERR) THEN
       WRITE (*, *) 'NetCDF error: ', NF90_STRERROR(ncstatus)
       STOP
    END IF
  END SUBROUTINE chk
END PROGRAM convert_sst
