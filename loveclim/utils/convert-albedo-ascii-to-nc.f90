PROGRAM convert_albedo

  ! Convert LOVECLIM ASCII albedo to NetCDF for further processing.

  USE netcdf

  IMPLICIT NONE

  INTEGER, PARAMETER :: NLON = 64, NLAT = 32
  DOUBLE PRECISION, PARAMETER :: LON0 = 0.0, DLON = 5.625
  DOUBLE PRECISION, PARAMETER :: LAT0 = -85.76, DLAT = 5.54

  INTEGER :: ncid
  INTEGER :: timdimid, latdimid, londimid
  INTEGER :: timvarid, latvarid, lonvarid, albvarid
  INTEGER :: ncstatus

  DOUBLE PRECISION, DIMENSION(NLON) :: lons
  DOUBLE PRECISION, DIMENSION(NLAT) :: lats

  INTEGER :: i, j, is
  REAL*4, DIMENSION(NLON, NLAT, 4) :: albedo


  OPEN (13, FILE='../fixed-data/land_albedo.dat', FORM='formatted')

  ! Create NetCDF file:
  !  Dimensions: lat, lon
  !  Variables: lat, lon, field variables

  CALL chk(NF90_CREATE('../fixed-data/land-albedo.nc', NF90_CLOBBER, ncid))

  CALL chk(NF90_DEF_DIM(ncid, 'time', 4,    timdimid))
  CALL chk(NF90_DEF_DIM(ncid, 'lat',  NLAT, latdimid))
  CALL chk(NF90_DEF_DIM(ncid, 'lon',  NLON, londimid))

  CALL chk(NF90_DEF_VAR(ncid, 'time', NF90_INT, (/ timdimid /), timvarid))
  CALL chk(NF90_PUT_ATT(ncid, timvarid, 'units', 'season'))
  CALL chk(NF90_PUT_ATT(ncid, timvarid, 'standard_name', 'time'))
  CALL chk(NF90_PUT_ATT(ncid, timvarid, 'long_name', 'Time'))
  CALL chk(NF90_PUT_ATT(ncid, timvarid, 'axis', 'T'))

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

  CALL chk(NF90_DEF_VAR(ncid, 'albedo', NF90_DOUBLE, &
       (/ londimid, latdimid, timdimid /), albvarid))
  CALL chk(NF90_PUT_ATT(ncid, albvarid, 'missing_value', -1.0D0))

  CALL chk(NF90_ENDDEF(ncid))

  DO i = 1, NLAT
     lats(i) = LAT0 + (i - 1) * DLAT
  END DO
  CALL chk(NF90_PUT_VAR(ncid, latvarid, lats))
  DO i = 1, NLON
     lons(i) = LON0 + (i - 1) * DLON
  END DO
  CALL chk(NF90_PUT_VAR(ncid, lonvarid, lons))
  CALL chk(NF90_PUT_VAR(ncid, timvarid, (/ 1, 2, 3, 4 /)))

  READ (13,*)
  DO i = 1, NLAT
     DO j = 1, NLON
        READ(13,100) (albedo(j,i,is), is = 1, 4)
     END DO
  END DO

100 FORMAT(4(2X,F8.4))

  WHERE (albedo == 0.0) albedo = -1.0D0

  CALL chk(NF90_PUT_VAR(ncid, albvarid, albedo, &
          (/ 1, 1, 1 /), (/ NLON, NLAT, 4 /)))

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
END PROGRAM convert_albedo
