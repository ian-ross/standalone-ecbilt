load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Read UM SST forcing.
sstinfp = addfile("$SRCDIR/um-sst.nc", "r")
umlat = sstinfp->latitude
umlon = sstinfp->longitude
sstin = sstinfp->temp(:, 0, :, :)

;; Read ECBILT land/sea mask.
ecmaskfp = addfile("$TMPDIR/ecbilt-mask.nc", "r")
eclat = ecmaskfp->lat
nlat = dimsizes(eclat)
eclon = ecmaskfp->lon
nlon = dimsizes(eclon)
ecmask = ecmaskfp->mask

;; Fill in missing values in SST fields.
guess = 1
nscan = 1500
eps = 1.0E-2
relc = 0.6
poisson_grid_fill(sstin, False, guess, nscan, eps, relc, 0)

;; Regrid SST and mask by ECBILT land/sea mask.
ecsst = area_hi2lores_Wrap(umlon, umlat, sstin, False, 1, eclon, eclat, False)
;;ecsst = mask(ecsst, ecmask, 0)

;; Middle of month days, time points to interpolate to and values to
;; use for output days.
mdays = ispan(15, 15 + 11 * 30, 30)
idays = ispan(1, 360, 1) - 0.5
odays = ispan(0, 359, 1)

;; Create output data field.
sstout = new((/ 360, nlat, nlon /), double)
sstout@units = "K"
sstout!0 = "time"
sstout!1 = "lat"
sstout!2 = "lon"
odays@standard_name = "time"
odays@long_name = "Time"
odays@calendar = "360_day"
odays@units = "days since 0000-01-01 00:00:00"
sstout&time = odays
sstout&lat = eclat
sstout&lon = eclon

;; Do periodic interpolation in time (from monthly to daily data) for
;; each non-missing data point.
do i = 0, nlat - 1
  do j = 0, nlon - 1
    if (ismissing(ecsst(0, i, j))) then
      sstout(:, i, j) = sstout@_FillValue
    else
      sstout(:, i, j) = ftcurvp(mdays, ecsst(:, i, j), 360, idays)
    end if
  end do
end do
sstout = where(sstout .lt. 272.15, 272.20, sstout)
sstout = where(sstout .ge. 310.00, 309.90, sstout)

;; Write NetCDF output.
system("/bin/rm -f $TMPDIR/ecbilt-sst.nc")
ofile = addfile("$TMPDIR/ecbilt-sst.nc", "c")
ofile->sst = sstout
