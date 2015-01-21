ecfp = addfile("../inidata/ecbilt-orig-mask.nc", "r")
eclat = ecfp->lat
nlat = dimsizes(eclat)

out = new((/ 2 * nlat /), string)
mskrow = "0000000000000000000000000000000000000000000000000000000000000000"
do i = 0, nlat - 1
  row = sprinti("%2d", nlat - 1 - i)
  out(i)        = str_concat((/ "  ", row, " ", mskrow /))
  out(i + nlat) = str_concat((/ "  ", row, " ", mskrow /))
end do

asciiwrite("../inidata/lakemask.dat", out)
