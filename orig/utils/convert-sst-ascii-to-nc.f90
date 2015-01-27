PROGRAM convert_sst

  ! Convert ECBILT ASCII SST forcing files to NetCDF.

  USE netcdf

  IMPLICIT NONE

  INTEGER, PARAMETER :: NLON = 64, NLAT = 32
  DOUBLE PRECISION, PARAMETER :: LON0 = 0.0, DLON = 5.625
  DOUBLE PRECISION, PARAMETER :: LAT0 = -85.76, DLAT = 5.54

  INTEGER :: iday

  INTEGER :: ncid
  INTEGER :: latdimid, londimid, timdimid
  INTEGER :: latvarid, lonvarid, timvarid, sstvarid
  INTEGER :: ncstatus

  DOUBLE PRECISION, DIMENSION(NLON) :: lons
  DOUBLE PRECISION, DIMENSION(NLAT) :: lats

  INTEGER :: i, j
  REAL*4, DIMENSION(NLAT, NLON) :: sstday


  OPEN (14, FILE='sst_daily.asc', FORM='formatted')

  ! Create NetCDF file:
  !  Dimensions: lat, lon, lev, lev2, time
  !  Variables: lat, lon, lev, lev2, time, field variables

  CALL chk(NF90_CREATE('sst_daily.nc', NF90_CLOBBER, ncid))

  CALL chk(NF90_DEF_DIM(ncid, 'lat', NLAT, latdimid))
  CALL chk(NF90_DEF_DIM(ncid, 'lon', NLON, londimid))
  CALL chk(NF90_DEF_DIM(ncid, 'time', NF90_UNLIMITED, timdimid))

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

  CALL chk(NF90_DEF_VAR(ncid, 'time', NF90_DOUBLE, (/ timdimid /), timvarid))
  CALL chk(NF90_PUT_ATT(ncid, timvarid, 'units', 'days since 1-1-1 00:00:0.0'))
  CALL chk(NF90_PUT_ATT(ncid, timvarid, 'standard_name', 'time'))
  CALL chk(NF90_PUT_ATT(ncid, timvarid, 'long_name', 'Time'))
  CALL chk(NF90_PUT_ATT(ncid, timvarid, 'axis', 'T'))
  CALL chk(NF90_PUT_ATT(ncid, timvarid, 'calendar', '360_day'))

  CALL chk(NF90_DEF_VAR(ncid, 'sst', NF90_DOUBLE, &
       (/ londimid, latdimid, timdimid /), sstvarid))
  CALL chk(NF90_PUT_ATT(ncid, sstvarid, 'units', 'degC'))
  CALL chk(NF90_PUT_ATT(ncid, sstvarid, 'missing_value', 9.99E10))

  CALL chk(NF90_ENDDEF(ncid))

  DO i = 1, NLAT
     lats(i) = LAT0 + (i - 1) * DLAT
  END DO
  CALL chk(NF90_PUT_VAR(ncid, latvarid, lats))
  DO i = 1, NLON
     lons(i) = LON0 + (i - 1) * DLON
  END DO
  CALL chk(NF90_PUT_VAR(ncid, lonvarid, lons))

  DO iday = 1, 360
     READ (14, *) ((sstday(i, j), j = 1, NLON), i = 1, NLAT)

     CALL chk(NF90_PUT_VAR(ncid, timvarid, (/ iday - 1 /), (/ iday /), (/ 1 /)))
     CALL chk(NF90_PUT_VAR(ncid, sstvarid, TRANSPOSE(sstday), &
          (/ 1, 1, iday /), (/ NLON, NLAT, 1 /)))
  END DO

200 CONTINUE

    CALL chk(NF90_CLOSE(ncid))

CONTAINS
  FUNCTION convtime(iyear, imonth, iday)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iyear, imonth, iday
    DOUBLE PRECISION :: convtime

    convtime = REAL((iyear - 1) * 360 + (imonth - 1) * 30 + iday - 1)
  END FUNCTION convtime

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
