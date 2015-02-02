load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Read original albedo.
albinfp = addfile("$FIXDIR/land-albedo.nc", "r")
albinlat = albinfp->lat
albinlon = albinfp->lon
albedoin = albinfp->albedo

;; Read ECBILT land/sea masks.
ecmaskfp = addfile("$TMPDIR/ecbilt-mask.nc", "r")
eclat = ecmaskfp->lat
nlat = dimsizes(eclat)
eclon = ecmaskfp->lon
nlon = dimsizes(eclon)
ecmask = ecmaskfp->mask

;; Fill in missing values in albedo fields.
guess = 1
nscan = 1500
eps = 1.0E-2
relc = 0.6
poisson_grid_fill(albedoin, False, guess, nscan, eps, relc, 0)

;; Regrid albedo and mask by ECBILT land/sea mask.
ecalbedo = new((/ 4, nlat, nlon /), double)
do is = 0, 3
  ecalbedo(is, :, :) = \
    area_hi2lores_Wrap(albinlon, albinlat, albedoin(is, :, :), \
                       False, 1, eclon, eclat, False)
  ecalbedo(is, :, :) = mask(ecalbedo(is, :, :), ecmask, 1)
end do

;; Write NetCDF output.
system("/bin/rm -f $TMPDIR/land-albedo.nc")
ofile = addfile("$TMPDIR/land-albedo.nc", "c")
ofile->albedo = ecalbedo

;; Write ASCII output.
do is = 0, 3
  ecalbedo(is, :, :) = where(ecmask .eq. 1, ecalbedo(is, :, :), 0.0)
end do
out = new(nlat * nlon + 1, string)
out(0) = " ***   DJF       MAM       JJA       SON"
iout = 1
do ilat = 0, nlat - 1
  do ilon = 0, nlon - 1
      of1 = sprintf("%10.4f", ecalbedo(0, ilat, ilon))
      of2 = sprintf("%10.4f", ecalbedo(1, ilat, ilon))
      of3 = sprintf("%10.4f", ecalbedo(2, ilat, ilon))
      of4 = sprintf("%10.4f", ecalbedo(3, ilat, ilon))
      out(iout) = str_concat((/ of1, of2, of3, of4 /))
      iout = iout + 1
   end do
end do
asciiwrite("$DSTDIR/inputdata/land_albedo.dat", out)
