load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Read UM orography.
oroginfp = addfile("../../um-inidata/tcszd.qrparm.orog.nc", "r")
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
ecmaskfp = addfile("../inidata/ecbilt-um-mask.nc", "r")
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
ecorog = mask(ecorog, .not. ecmask, 0)

;; Calculate land fraction and mask by ECBILT land/sea mask.
landfrac = area_hi2lores_Wrap(topolon, topolat, topomask, False, 1, \
                              eclon, eclat, False)
eclfrac = mask(landfrac, .not. ecmask, 0)

;; Write NetCDF output.
system("/bin/rm -f ../inidata/ecbilt-orog.nc")
ofile = addfile("../inidata/ecbilt-orog.nc", "c")
ofile->orog = ecorog
ofile->lfrac = eclfrac

;; Write ASCII output.
ecorog = where(ismissing(ecorog), 0.0, ecorog)
replace_ieeenan(eclfrac, 0.0, 0)
oo = new(nlat * nlon, string)
of = new(nlat * nlon, string)
iout = 0
do ilon = 0, nlon - 1
   do ilat = 0, nlat - 1
      oo(iout) = sprintf("%24.16f", ecorog(ilat, ilon))
      of(iout) = sprintf("%24.16f", eclfrac(ilat, ilon))
      iout = iout + 1
   end do
end do
noutrows = (nlat * nlon - 1) / 3 + 1
out = new(2 * noutrows, string)
do i = 0, noutrows - 1
   if (i .lt. (noutrows - 1)) then
      out(i) = str_concat((/ oo(3*i), oo(3*i+1), oo(3*i+2) /))
      out(i+noutrows) = str_concat((/ of(3*i), of(3*i+1), of(3*i+2) /))
   else
      out(i) = str_concat((/ oo(3*i), oo(3*i+1) /))
      out(i+noutrows) = str_concat((/ of(3*i), of(3*i+1) /))
   end if
end do
asciiwrite("../inidata/berg.asc", out)
