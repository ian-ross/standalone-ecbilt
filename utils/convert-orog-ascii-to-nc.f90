PROGRAM convert_orog

  ! Convert ECBILT ASCII orography to NetCDF.

  USE netcdf

  IMPLICIT NONE

  INTEGER, PARAMETER :: NLON = 64, NLAT = 32
  DOUBLE PRECISION, PARAMETER :: LON0 = 0.0, DLON = 5.625
  DOUBLE PRECISION, PARAMETER :: LAT0 = -85.76, DLAT = 5.54

  INTEGER :: ncid
  INTEGER :: latdimid, londimid
  INTEGER :: latvarid, lonvarid, orogvarid, sdfricvarid
  INTEGER :: ncstatus

  DOUBLE PRECISION, DIMENSION(NLON) :: lons
  DOUBLE PRECISION, DIMENSION(NLAT) :: lats

  INTEGER :: i, j, k
  REAL*4, DIMENSION(NLAT, NLON) :: orog, sdfric
  INTEGER, DIMENSION(0:NLON, 0:NLAT) :: indu
  INTEGER, DIMENSION(0:NLON-1, 0:NLAT-1) :: indt
  INTEGER, DIMENSION(NLAT, NLON) :: lsmask


  OPEN (13, FILE='inputdata/atmos/berg.asc', FORM='formatted')
  OPEN (40, FILE='inputdata/ocean/mask.dat')

  ! Create NetCDF file:
  !  Dimensions: lat, lon
  !  Variables: lat, lon, field variables

  CALL chk(NF90_CREATE('orog.nc', NF90_CLOBBER, ncid))

  CALL chk(NF90_DEF_DIM(ncid, 'lat', NLAT, latdimid))
  CALL chk(NF90_DEF_DIM(ncid, 'lon', NLON, londimid))

  CALL chk(NF90_DEF_VAR(ncid, 'lat', NF90_DOUBLE, (/ latdimid /), latvarid))
  CALL chk(NF90_PUT_ATT(ncid, latvarid, 'units', 'degrees_north'))
  CALL chk(NF90_PUT_ATT(ncid, latvarid, 'standard_name', 'latitude'))
  CALL chk(NF90_PUT_ATT(ncid, latvarid, 'long_name', 'Latitude'))
  CALL chk(NF90_PUT_ATT(ncid, latvarid, 'axis', 'Y'))

  CALL chk(NF90_DEF_VAR(ncid, 'lon', NF90_DOUBLE, (/ londimid /), lonvarid))
  CALL chk(NF90_PUT_ATT(ncid, lonvarid, 'units', 'degrees_east'))
  CALL chk(NF90_PUT_ATT(ncid, lonvarid, 'standard_name', 'longitude'))
  CALL chk(NF90_PUT_ATT(ncid, lonvarid, 'long_name', 'Longitude'))
  CALL chk(NF90_PUT_ATT(ncid, lonvarid, 'axis', 'X'))

  CALL chk(NF90_DEF_VAR(ncid, 'orog', NF90_DOUBLE, &
       (/ londimid, latdimid /), orogvarid))
  CALL chk(NF90_PUT_ATT(ncid, orogvarid, 'units', 'm'))
  CALL chk(NF90_PUT_ATT(ncid, orogvarid, 'missing_value', 9.99E10))

  CALL chk(NF90_DEF_VAR(ncid, 'sdfric', NF90_DOUBLE, &
       (/ londimid, latdimid /), sdfricvarid))
  CALL chk(NF90_PUT_ATT(ncid, sdfricvarid, 'missing_value', 9.99E10))

  CALL chk(NF90_ENDDEF(ncid))

  DO i = 1, NLAT
     lats(i) = LAT0 + (i - 1) * DLAT
  END DO
  CALL chk(NF90_PUT_VAR(ncid, latvarid, lats))
  DO i = 1, NLON
     lons(i) = LON0 + (i - 1) * DLON
  END DO
  CALL chk(NF90_PUT_VAR(ncid, lonvarid, lons))

  READ (13, *) orog
  READ (13, *) sdfric
  READ (40,120)
  DO j = NLAT, 0, -1
     READ (40, 120) (indu(i, j), i = 0, NLON)
  END DO
  DO j = NLAT - 1, 0, -1
     READ (40, 310) k, (indt(i, j), i = 0, NLON - 1)
  END DO
  DO i = 1, NLAT
     DO j = 1, NLON
        lsmask(i, j) = indt(j - 1, i - 1)
     END DO
  END DO

310 FORMAT (i4, i2, 90i1)
120 FORMAT (65i1)

  WHERE (lsmask == 1)
     orog = 9.99E10
     sdfric = 9.99E10
  END WHERE

  CALL chk(NF90_PUT_VAR(ncid, orogvarid, TRANSPOSE(orog), &
          (/ 1, 1 /), (/ NLON, NLAT /)))
  CALL chk(NF90_PUT_VAR(ncid, sdfricvarid, TRANSPOSE(sdfric), &
          (/ 1, 1 /), (/ NLON, NLAT /)))

200 CONTINUE

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
END PROGRAM convert_orog
