load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Read UM orography.
oroginfp = addfile("$SRCDIR/um-orog.nc", "r")
umlat = oroginfp->latitude
umlon = oroginfp->longitude
orogin = oroginfp->ht(0, 0, :, :)

;; Read ETOPO2 topography (used for land fraction calculation).
topoinfp = addfile("/big/data/ETOPO2v2c_f4.nc", "r")
topolat = topoinfp->lat
topolon = topoinfp->lon
topolon = lonFlip(topolon)
topolon = mod(topolon + 360.0, 360.0)
topomask = topoinfp->orog
topomask = lonFlip(topomask)
topomask = where(topomask .gt. 0, 1, 0)

;; Read ECBILT land/sea mask.
ecmaskfp = addfile("$TMPDIR/ecbilt-mask.nc", "r")
eclat = ecmaskfp->lat
nlat = dimsizes(eclat)
eclon = ecmaskfp->lon
nlon = dimsizes(eclon)
ecmask = ecmaskfp->mask

;; Fill in missing values in orography fields.
guess = 1
nscan = 1500
eps = 1.0E-2
relc = 0.6
poisson_grid_fill(orogin, False, guess, nscan, eps, relc, 0)

;; Regrid orography and mask by ECBILT land/sea mask.
ecorog = area_hi2lores_Wrap(umlon, umlat, orogin, False, 1, eclon, eclat, False)

;; Calculate land fraction and mask by ECBILT land/sea mask.
landfrac = area_hi2lores_Wrap(topolon, topolat, topomask, False, 1, \
                              eclon, eclat, False)
landfrac = where(ecmask, landfrac, 0.0)
oceanfrac = 1.0 - landfrac

;; Write NetCDF output.
system("/bin/rm -f $DSTDIR/inputdata/berg.nc")
ofile = addfile("$DSTDIR/inputdata/berg.nc", "c")
ofile->h = ecorog
ofile->sfric = landfrac

;; Write ASCII output.
bergorog = ecorog
bergorog = where(ismissing(bergorog), 0.0, bergorog)
berglfrac = landfrac
replace_ieeenan(berglfrac, 0.0, 0)
ooberg = new(nlat * nlon, string)
ofberg = new(nlat * nlon, string)
of = new(nlat * nlon, string)
iout = 0
do ilon = 0, nlon - 1
   do ilat = 0, nlat - 1
      ooberg(iout) = sprintf("%24.16f", bergorog(ilat, ilon))
      ofberg(iout) = sprintf("%24.16f", berglfrac(ilat, ilon))
      iout = iout + 1
   end do
end do
iout = 0
do ilat = 0, nlat-1
  do ilon = 0, nlon - 1
      of(iout) = sprintf("%24.16f", oceanfrac(ilat, ilon))
      iout = iout + 1
   end do
end do
noutrows = (nlat * nlon - 1) / 3 + 1
out = new(noutrows + 1, string)
out(0) = "fraction d ocean"
do i = 0, noutrows - 1
   if (i .lt. (noutrows - 1)) then
      out(i + 1) = str_concat((/ of(3*i), of(3*i+1), of(3*i+2) /))
   else
      out(i + 1) = str_concat((/ of(3*i), of(3*i+1) /))
   end if
end do
asciiwrite("$DSTDIR/inputdata/fractoc.dat", out)
noutrowsberg = (nlat * nlon - 1) / 3 + 1
outberg = new(2 * noutrowsberg, string)
do i = 0, noutrowsberg - 1
   if (i .lt. (noutrowsberg - 1)) then
      outberg(i) = str_concat((/ ooberg(3*i), ooberg(3*i+1), ooberg(3*i+2) /))
      outberg(i+noutrowsberg) = \
          str_concat((/ ofberg(3*i), ofberg(3*i+1), ofberg(3*i+2) /))
   else
      outberg(i) = str_concat((/ ooberg(3*i), ooberg(3*i+1) /))
      outberg(i+noutrowsberg) = \
          str_concat((/ ofberg(3*i), ofberg(3*i+1) /))
   end if
end do
asciiwrite("$TMPDIR/berg.asc", outberg)
