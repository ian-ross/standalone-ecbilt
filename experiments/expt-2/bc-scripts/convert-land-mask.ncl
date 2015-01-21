umfp = addfile("../../um-inidata/tcszd.qrparm.mask.nc", "r")
ecfp = addfile("../inidata/ecbilt-orig-mask.nc", "r")
umlat = umfp->latitude
umlon = umfp->longitude
ummask = umfp->lsm
eclat = ecfp->lat
nlat = dimsizes(eclat)
eclon = ecfp->lon
nlon = dimsizes(eclon)
wgt = cos(umlat)
ecmask1 = area_hi2lores(umlon, umlat, ummask, False, 1, eclon, eclat, False)
ecmask = reshape(ecmask1, (/ nlat, nlon /))
ecmask!0 = "lat"
ecmask!1 = "lon"
ecmask&lat = eclat
ecmask&lon = eclon
msk = where(ecmask .ge. 0.5, 1, 0)
msk!0 = "lat"
msk!1 = "lon"
msk&lat = eclat
msk&lon = eclon
system("/bin/rm -f ../inidata/ecbilt-um-mask.nc")
ofile = addfile("../inidata/ecbilt-um-mask.nc", "c")
ofile->mask = msk

array0 = tocharacter("0")
char0 = array0(0)
array1 = tocharacter("1")
char1 = array1(0)
chmask = new((/ nlat /), string)
chmask = where(msk, char1, char0)
chmaskneg = new((/ nlat /), string)
chmaskneg = where(msk, char0, char1)

out = new((/ 2 * nlat + 1 /), string)
out(0) = "12345678901234567890123456789012345678901234567890123456789012345"
do i = 0, nlat - 1
  mskrow = chmask(nlat - 1 - i)
  tmp = stringtocharacter(mskrow)
  rpt = charactertostring(tmp(nlon - 1))
  out(i + 1) = str_concat((/ mskrow, rpt /))
  row = sprinti("%2d", nlat - 1 - i)
  mskrowneg = chmaskneg(nlat - 1 - i)
  out(i + 1 + nlat) = str_concat((/ "  ", row, " ", mskrowneg /))
end do

asciiwrite("../inidata/mask.dat", out)
