load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Read original forest-fraction.
forfrinfp = addfile("$FIXDIR/forest-fraction.nc", "r")
forfrinlat = forfrinfp->lat
forfrinlon = forfrinfp->lon
forfrin = forfrinfp->forestfr

;; Read ECBILT land/sea masks.
ecmaskfp = addfile("$TMPDIR/ecbilt-mask.nc", "r")
eclat = ecmaskfp->lat
nlat = dimsizes(eclat)
eclon = ecmaskfp->lon
nlon = dimsizes(eclon)
ecmask = ecmaskfp->mask

;; Fill in missing values in forest fraction field.
guess = 1
nscan = 1500
eps = 1.0E-2
relc = 0.6
poisson_grid_fill(forfrin, False, guess, nscan, eps, relc, 0)

;; Regrid forest fraction and mask by ECBILT land/sea mask.
ecforfr = area_hi2lores_Wrap(forfrinlon, forfrinlat, forfrin, \
                             False, 1, eclon, eclat, False)
ecforfr = mask(ecforfr, ecmask, 1)

;; Write NetCDF output.
system("/bin/rm -f $TMPDIR/forest-fraction.nc")
ofile = addfile("$TMPDIR/forest-fraction.nc", "c")
ofile->forestfr = ecforfr

;; Write ASCII output.
ecforfr = where(ecmask .eq. 1, ecforfr, 0.0)
out = new(nlat * nlon, string)
iout = 0
do ilat = 0, nlat - 1
  do ilon = 0, nlon - 1
      out(iout) = sprintf("%8.4f", ecforfr(ilat, ilon))
      iout = iout + 1
   end do
end do
asciiwrite("$DSTDIR/inputdata/forfr.dat", out)
