PROGRAM convert_restart

  ! Convert LOVECLIM binary restart files to NetCDF for debugging.

  USE netcdf

  IMPLICIT NONE

  INTEGER, PARAMETER :: NLON = 64, NLAT = 32, NTYP=3
  DOUBLE PRECISION, PARAMETER :: LON0 = 0.0, DLON = 5.625
  DOUBLE PRECISION, PARAMETER :: LAT0 = -85.76, DLAT = 5.54

  CHARACTER(LEN=10) :: runid

  INTEGER :: ncid
  INTEGER :: latdimid, londimid, typdimid
  INTEGER :: latvarid, lonvarid, typvarid
  INTEGER :: tsurfnvarid, tempmvarid, temp0gvarid
  INTEGER :: rmoisgvarid, torainvarid, tosnowvarid
  INTEGER :: ncstatus

  DOUBLE PRECISION, DIMENSION(NLON) :: lons
  DOUBLE PRECISION, DIMENSION(NLAT) :: lats

  INTEGER :: i
  REAL*8 tsurfn(NLAT, NLON, NTYP), tempm(0:2), temp0g(NLAT, NLON)
  REAL*8 rmoisg(NLAT, NLON), torain(NLAT, NLON), tosnow(NLAT, NLON)

  CALL GETARG(1, runid)
  OPEN(95, FILE='inatphy'//runid//'.dat', FORM='unformatted')
  READ (95) tsurfn, tempm, temp0g
  READ (95) rmoisg, torain, tosnow
  CLOSE(95)

  CALL chk(NF90_CREATE('inatphy'//runid//'.nc', NF90_CLOBBER, ncid))

  CALL chk(NF90_DEF_DIM(ncid, 'lat',  NLAT, latdimid))
  CALL chk(NF90_DEF_DIM(ncid, 'lon',  NLON, londimid))
  CALL chk(NF90_DEF_DIM(ncid, 'typ',  NTYP, typdimid))

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

  CALL chk(NF90_DEF_VAR(ncid, 'typ', NF90_INT, (/ typdimid /), typvarid))

  CALL chk(NF90_DEF_VAR(ncid, 'tsurfn', NF90_DOUBLE, &
       (/ londimid, latdimid, typdimid /), tsurfnvarid))
  CALL chk(NF90_DEF_VAR(ncid, 'temp0g', NF90_DOUBLE, &
       (/ londimid, latdimid /), temp0gvarid))
  CALL chk(NF90_DEF_VAR(ncid, 'rmoisg', NF90_DOUBLE, &
       (/ londimid, latdimid /), rmoisgvarid))
  CALL chk(NF90_DEF_VAR(ncid, 'torain', NF90_DOUBLE, &
       (/ londimid, latdimid /), torainvarid))
  CALL chk(NF90_DEF_VAR(ncid, 'tosnow', NF90_DOUBLE, &
       (/ londimid, latdimid /), tosnowvarid))

  CALL chk(NF90_ENDDEF(ncid))

  DO i = 1, NLAT
     lats(i) = LAT0 + (i - 1) * DLAT
  END DO
  CALL chk(NF90_PUT_VAR(ncid, latvarid, lats))
  DO i = 1, NLON
     lons(i) = LON0 + (i - 1) * DLON
  END DO
  CALL chk(NF90_PUT_VAR(ncid, lonvarid, lons))

  DO i = 1, NTYP
     CALL chk(NF90_PUT_VAR(ncid, tsurfnvarid, TRANSPOSE(tsurfn(:,:,i)), &
          & (/ 1, 1, i /), (/ NLON, NLAT, 1 /)))
  END DO
  CALL chk(NF90_PUT_VAR(ncid, temp0gvarid, TRANSPOSE(temp0g), &
       & (/ 1, 1 /), (/ NLON, NLAT /)))
  CALL chk(NF90_PUT_VAR(ncid, rmoisgvarid, TRANSPOSE(rmoisg), &
       & (/ 1, 1 /), (/ NLON, NLAT /)))
  CALL chk(NF90_PUT_VAR(ncid, torainvarid, TRANSPOSE(torain), &
       & (/ 1, 1 /), (/ NLON, NLAT /)))
  CALL chk(NF90_PUT_VAR(ncid, tosnowvarid, TRANSPOSE(tosnow), &
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
END PROGRAM convert_restart
