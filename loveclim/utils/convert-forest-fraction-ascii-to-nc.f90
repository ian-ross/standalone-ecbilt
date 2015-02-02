PROGRAM convert_forfrac

  ! Convert LOVECLIM ASCII albedo to NetCDF for further processing.

  USE netcdf

  IMPLICIT NONE

  INTEGER, PARAMETER :: NLON = 64, NLAT = 32
  DOUBLE PRECISION, PARAMETER :: LON0 = 0.0, DLON = 5.625
  DOUBLE PRECISION, PARAMETER :: LAT0 = -85.76, DLAT = 5.54

  INTEGER :: ncid
  INTEGER :: latdimid, londimid
  INTEGER :: latvarid, lonvarid, forfrvarid
  INTEGER :: maskncid, maskvarid
  INTEGER :: ncstatus

  DOUBLE PRECISION, DIMENSION(NLON) :: lons
  DOUBLE PRECISION, DIMENSION(NLAT) :: lats

  INTEGER :: i, j, is
  REAL*4, DIMENSION(NLON, NLAT) :: forestfr, mask


  OPEN (13, FILE='../fixed-data/forest-fraction.dat', FORM='formatted')

  CALL chk(NF90_OPEN('../fixed-data/ecbilt-orig-mask.nc', 0, maskncid))
  CALL chk(NF90_INQ_VARID(maskncid, 'mask', maskvarid))
  CALL chk(NF90_GET_VAR(maskncid, maskvarid, mask))
  CALL chk(NF90_CLOSE(maskncid))

  CALL chk(NF90_CREATE('../fixed-data/forest-fraction.nc', NF90_CLOBBER, ncid))

  CALL chk(NF90_DEF_DIM(ncid, 'lat',  NLAT, latdimid))
  CALL chk(NF90_DEF_DIM(ncid, 'lon',  NLON, londimid))

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

  CALL chk(NF90_DEF_VAR(ncid, 'forestfr', NF90_DOUBLE, &
       (/ londimid, latdimid /), forfrvarid))
  CALL chk(NF90_PUT_ATT(ncid, forfrvarid, 'missing_value', -1.0D0))

  CALL chk(NF90_ENDDEF(ncid))

  DO i = 1, NLAT
     lats(i) = LAT0 + (i - 1) * DLAT
  END DO
  CALL chk(NF90_PUT_VAR(ncid, latvarid, lats))
  DO i = 1, NLON
     lons(i) = LON0 + (i - 1) * DLON
  END DO
  CALL chk(NF90_PUT_VAR(ncid, lonvarid, lons))

  DO i = 1, NLAT
     DO j = 1, NLON
        READ (13,110) forestfr(j,i)
     END DO
  END DO

110 FORMAT(F8.4)

  WHERE (mask == 0.0) forestfr = -1.0D0

  CALL chk(NF90_PUT_VAR(ncid, forfrvarid, forestfr, &
       & (/ 1, 1 /), (/ NLON, NLAT /)))

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
END PROGRAM convert_forfrac
