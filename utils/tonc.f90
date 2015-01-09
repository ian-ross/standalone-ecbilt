PROGRAM tonc

  ! Convert ECBILT output files to NetCDF.

  USE netcdf

  IMPLICIT NONE

  INTEGER, PARAMETER :: MAXCODE = 1000
  INTEGER, PARAMETER :: NLON = 64, NLAT = 32
  DOUBLE PRECISION, PARAMETER :: LON0 = 0.0, DLON = 5.625
  DOUBLE PRECISION, PARAMETER :: LAT0 = -85.76, DLAT = 5.54

  LOGICAL, DIMENSION(MAXCODE) :: codepresent = .FALSE.
  INTEGER, DIMENSION(MAXCODE) :: codencid
  INTEGER :: ntime, itime
  INTEGER :: idate, idateold, idate0
  INTEGER :: icode, icodeold, ilevel, iyear, imonth, iday, ix, iy, ifield
  INTEGER :: ilev

  INTEGER :: ncid
  INTEGER :: latdimid, londimid, lvdimid, lv2dimid, timdimid
  INTEGER :: latvarid, lonvarid, lvvarid, lv2varid, timvarid
  INTEGER :: ncstatus
  INTEGER, DIMENSION(:), ALLOCATABLE :: vdims
  CHARACTER(LEN=10) :: vunits

  INTEGER, DIMENSION(3) :: tlevels = (/ 1000, 650, 350 /)
  INTEGER, DIMENSION(3) :: ulevels = (/ 800, 500, 200 /)
  DOUBLE PRECISION, DIMENSION(NLON) :: lons
  DOUBLE PRECISION, DIMENSION(NLAT) :: lats

  INTEGER :: i, j
  INTEGER :: reclen
  REAL*4 :: x(NLON, NLAT)
  CHARACTER(LEN=256) :: infn, outfn

  CALL GETARG(1, infn)
  IF (infn .EQ. '') THEN
     WRITE (*,*) 'Usage: tonc <filename>'
     STOP
  END IF
  outfn = TRIM(infn) // '.nc'


  ! First pass through input file: determine start time, variables for
  ! output.

  OPEN (999, FILE=infn, FORM='unformatted')
  INQUIRE (IOLENGTH=reclen) x
  READ (999, END=999) icode, ilevel, iyear, imonth, iday, ix, iy, ifield
  READ (999, END=999) ((x(i, j), j = 1, NLAT), i = 1, NLON)
  codepresent(icode) = .TRUE.
  idate = 20000000 + iyear * 10000 + imonth * 100 + iday
  idate0 = idate
  ntime = 1
  idateold = idate
  DO
     READ (999, END=100) icode, ilevel, iyear, imonth, iday, ix, iy, ifield
     READ (999, END=100) ((x(i, j), j = 1, NLAT), i = 1, NLON)
     codepresent(icode) = .TRUE.
     idate = 20000000 + iyear * 10000 + imonth * 100 + iday
     IF (idate .NE. idateold) THEN
        ntime = ntime + 1
        idateold = idate
     END IF
  END DO

999 WRITE (*, *) 'No records in input file!'
  STOP

100 CONTINUE


  ! Create NetCDF file:
  !  Dimensions: lat, lon, lev, lev2, time
  !  Variables: lat, lon, lev, lev2, time, field variables

  CALL chk(NF90_CREATE(outfn, NF90_CLOBBER, ncid))

  CALL chk(NF90_DEF_DIM(ncid, 'lat', NLAT, latdimid))
  CALL chk(NF90_DEF_DIM(ncid, 'lon', NLON, londimid))
  CALL chk(NF90_DEF_DIM(ncid, 'lev', 3, lvdimid))
  CALL chk(NF90_DEF_DIM(ncid, 'lev2', 3, lv2dimid))
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

  CALL chk(NF90_DEF_VAR(ncid, 'lev', NF90_DOUBLE, (/ lvdimid /), lvvarid))
  CALL chk(NF90_PUT_ATT(ncid, lvvarid, 'units', 'millibars'))
  CALL chk(NF90_PUT_ATT(ncid, lvvarid, 'standard_name', 'level'))
  CALL chk(NF90_PUT_ATT(ncid, lvvarid, 'axis', 'Z'))

  CALL chk(NF90_DEF_VAR(ncid, 'lev2', NF90_DOUBLE, (/ lv2dimid /), lv2varid))
  CALL chk(NF90_PUT_ATT(ncid, lv2varid, 'units', 'millibars'))
  CALL chk(NF90_PUT_ATT(ncid, lv2varid, 'standard_name', 'level'))
  CALL chk(NF90_PUT_ATT(ncid, lv2varid, 'axis', 'Z'))

  CALL chk(NF90_DEF_VAR(ncid, 'time', NF90_DOUBLE, (/ timdimid /), timvarid))
  CALL chk(NF90_PUT_ATT(ncid, timvarid, 'units', 'days since 1-1-1 00:00:0.0'))
  CALL chk(NF90_PUT_ATT(ncid, timvarid, 'standard_name', 'time'))
  CALL chk(NF90_PUT_ATT(ncid, timvarid, 'long_name', 'Time'))
  CALL chk(NF90_PUT_ATT(ncid, timvarid, 'axis', 'T'))
  CALL chk(NF90_PUT_ATT(ncid, timvarid, 'calendar', '360_day'))

  DO icode = 1, MAXCODE
     IF (codepresent(icode)) THEN
        vdims = vardims(icode, &
             latdimid, londimid, timdimid, lvdimid, lv2dimid)
        CALL chk(NF90_DEF_VAR(ncid, varname(icode, ifield), NF90_DOUBLE, &
             vdims, codencid(icode)))
        DEALLOCATE(vdims)
        vunits = varunits(icode)
        IF (TRIM(vunits) /= '') THEN
           CALL chk(NF90_PUT_ATT(ncid, codencid(icode), 'units', vunits))
        END IF
     END IF
  END DO

  CALL chk(NF90_ENDDEF(ncid))

  DO i = 1, NLAT
     lats(i) = LAT0 + (i - 1) * DLAT
  END DO
  CALL chk(NF90_PUT_VAR(ncid, latvarid, lats))
  DO i = 1, NLON
     lons(i) = LON0 + (i - 1) * DLON
  END DO
  CALL chk(NF90_PUT_VAR(ncid, lonvarid, lons))
  CALL chk(NF90_PUT_VAR(ncid, lvvarid, tlevels))
  CALL chk(NF90_PUT_VAR(ncid, lv2varid, ulevels))


  ! Second pass through input file: read field slices and write to
  ! output.

  REWIND (999)
  idateold = -1
  itime = 0
  icodeold = -1
  ilev = 1
  DO
     READ (999, END=200) icode, ilevel, iyear, imonth, iday, ix, iy, ifield
     READ (999, END=200) ((x(i, j), j = 1, NLAT), i = 1, NLON)
     idate = 20000000 + iyear * 10000 + imonth * 100 + iday
     IF (idate .NE. idateold) THEN
        itime = itime + 1
        CALL chk(NF90_PUT_VAR(ncid, timvarid, &
             (/ convtime(iyear, imonth, iday) /), (/ itime /), (/ 1 /)))
        idateold = idate
     END IF
     IF (icode .NE. icodeold) THEN
        ilev = 1
        icodeold = icode
     END IF
     IF (varlevels(icode) == 1) THEN
        CALL chk(NF90_PUT_VAR(ncid, codencid(icode), TRANSPOSE(x), &
             (/ 1, 1, itime /), (/ NLON, NLAT, 1 /)))
     ELSE
        CALL chk(NF90_PUT_VAR(ncid, codencid(icode), TRANSPOSE(x), &
             (/ 1, 1, ilev, itime /), (/ NLON, NLAT, 1, 1 /)))
        ilev = ilev + 1
     END IF
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

  FUNCTION vardims(icode, latdimid, londimid, timdimid, lvdimid, lv2dimid)
    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: icode
    INTEGER, INTENT(IN)   :: latdimid, londimid, timdimid
    INTEGER, INTENT(IN)   :: lvdimid, lv2dimid
    INTEGER, DIMENSION(:), ALLOCATABLE :: vardims

    SELECT CASE (icode)
    CASE (134, 139, 140, 308, 141, 142, 143, 146, 147, 133, 157, 159, &
          160, 161, 164, 182, 174, 175, 176, 177, 178, 179, 180, 181, &
          260, 309, 307, 305, 306)
       ALLOCATE(vardims(3))
       vardims = (/ londimid, latdimid, timdimid /)
    CASE (130, 135, 301, 996, 997)
       ALLOCATE(vardims(4))
       vardims = (/ londimid, latdimid, lvdimid, timdimid /)
    CASE (131, 132, 148, 149, 156, 302, 303, 304, 310, 998, 999)
       ALLOCATE(vardims(4))
       vardims = (/ londimid, latdimid, lv2dimid, timdimid /)
    CASE DEFAULT
       WRITE (*, *) 'Unknown code for level determination: ', icode
       STOP
    END SELECT
  END FUNCTION vardims

  FUNCTION varlevels(icode)
    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: icode
    INTEGER :: varlevels

    SELECT CASE (icode)
    CASE (134, 139, 140, 308, 141, 142, 143, 146, 147, 133, 157, 159, &
          160, 161, 164, 182, 174, 175, 176, 177, 178, 179, 180, 181, &
          260, 309, 307, 305, 306)
       varlevels = 1
    CASE (130, 135, 301, 996, 997, 131, 132, 148, 149, 156, 302, 303, &
          304, 310, 998, 999)
       varlevels = 3
    CASE DEFAULT
       WRITE (*, *) 'Unknown code for level determination: ', icode
       STOP
    END SELECT
  END FUNCTION varlevels

  FUNCTION varname(icode, ifield)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: icode, ifield
    CHARACTER(LEN=10)    :: varname

    SELECT CASE (icode)
    CASE (134)
       varname = 'sp'
    CASE (139)
       varname = 'ts'
    CASE (140)
       varname = 'bmu'
    CASE (308)
       varname = 'bml'
    CASE (141)
       varname = 'sdl'
    CASE (142)
       varname = 'lsp'
    CASE (143)
       varname = 'cp'
    CASE (146)
       varname = 'shf'
    CASE (147)
       varname = 'lhf'
    CASE (133)
       varname = 'q'
    CASE (157)
       varname = 'r'
    CASE (159)
       varname = 'uv10'
    CASE (160)
       varname = 'runoffo'
    CASE (161)
       varname = 'runoffl'
    CASE (164)
       varname = 'tcc'
    CASE (182)
       varname = 'evap'
    CASE (174)
       varname = 'palb'
    CASE (175)
       varname = 'alb'
    CASE (176)
       varname = 'ssr'
    CASE (177)
       varname = 'str'
    CASE (178)
       varname = 'tsr'
    CASE (179)
       varname = 'ttr'
    CASE (180)
       varname = 'ustress'
    CASE (181)
       varname = 'vstress'
    CASE (260)
       varname = 'pp'
    CASE (309)
       varname = 'eminp'
    CASE (130)
       varname = 't'
    CASE (135)
       varname = 'omega'
    CASE (301)
       varname = 'hforc'
    CASE (307)
       varname = 'richarson'
    CASE (305)
       varname = 'cdragw'
    CASE (306)
       varname = 'cdragv'
    CASE (996)
       varname = 'dumt1'
    CASE (997)
       varname = 'dumt2'
    CASE (131)
       varname = 'u'
    CASE (132)
       varname = 'v'
    CASE (148)
       varname = 'psi'
    CASE (149)
       varname = 'chi'
    CASE (156)
       varname = 'gh'
    CASE (302)
       varname = 'vforc'
    CASE (303)
       varname = 'ageu'
    CASE (304)
       varname = 'agev'
    CASE (310)
       varname = 'pv'
    CASE (998)
       varname = 'dumu1'
    CASE (999)
       varname = 'dumu2'
    END SELECT
    IF (ifield .EQ. 2) varname = varname // 'd'
  END FUNCTION varname

  FUNCTION varunits(icode)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: icode
    CHARACTER(LEN=10)    :: varunits

    SELECT CASE (icode)
    CASE (139, 130) ! ts, t
       varunits = 'degC'
    CASE (141)
       varunits = 'sdl'
    CASE (146, 147, 176, 177, 178, 179) ! shf, lhf, ssr, str, tsr, ttr
       varunits = 'W m-2'
    CASE (133) ! q
       varunits = 'g m-2'
    CASE (260, 182) ! pp, evap
       varunits = 'cm yr-1'
    CASE (131, 132) ! u, v
       varunits = 'm s-1'
    CASE DEFAULT
       varunits = ''
    END SELECT
  END FUNCTION varunits

  SUBROUTINE chk(ncstatus)
    USE netcdf
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ncstatus

    IF (ncstatus /= NF90_NOERR) THEN
       WRITE (*, *) 'NetCDF error: ', NF90_STRERROR(ncstatus)
       STOP
    END IF
  END SUBROUTINE chk
END PROGRAM tonc
