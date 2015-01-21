load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Read UM orography.
oroginfp = addfile("../../um-inidata/tcszd.qrparm.orog.nc", "r")
umlat = oroginfp->latitude
umlon = oroginfp->longitude
orogin = oroginfp->ht(0, 0, :, :)

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

;; Write NetCDF output.
system("/bin/rm -f ../inidata/ecbilt-orog.nc")
ofile = addfile("../inidata/ecbilt-orog.nc", "c")
ofile->orog = ecorog
