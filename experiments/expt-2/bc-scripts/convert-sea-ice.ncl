load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Read UM sea-ice forcing.
iceinfp = addfile("../../um-inidata/tcszd.qrclim.ice.nc", "r")
umlat = iceinfp->latitude
umlon = iceinfp->longitude
iceconcin = iceinfp->iceconc(:, 0, :, :)
icedptin = iceinfp->field92(:, 0, :, :)

;; Read ECBILT land/sea mask.
ecmaskfp = addfile("../inidata/ecbilt-um-mask.nc", "r")
eclat = ecmaskfp->lat
nlat = dimsizes(eclat)
eclon = ecmaskfp->lon
nlon = dimsizes(eclon)
ecmask = ecmaskfp->mask

;; Fill in missing values in sea-ice fields.
guess = 1
nscan = 1500
eps = 1.0E-2
relc = 0.6
poisson_grid_fill(iceconcin, False, guess, nscan, eps, relc, 0)
poisson_grid_fill(icedptin, False, guess, nscan, eps, relc, 0)

;; Regrid sea-ice concentration and mask by ECBILT land/sea mask.
eciceconc = area_hi2lores_Wrap(umlon, umlat, iceconcin, \
                               False, 1, eclon, eclat, False)
eciceconc = mask(eciceconc, ecmask, 0)
ecicedpt = area_hi2lores_Wrap(umlon, umlat, icedptin, \
                              False, 1, eclon, eclat, False)
ecicedpt = mask(ecicedpt, ecmask, 0)

;; Middle of month days, time points to interpolate to and values to
;; use for output days.
mdays = ispan(15, 15 + 11 * 30, 30)

;; Write NetCDF output.
system("/bin/rm -f ../inidata/ecbilt-ice.nc")
ofile = addfile("../inidata/ecbilt-ice.nc", "c")
ofile->iceconc = eciceconc
ofile->icedepth = ecicedpt

;; Set up ice birth and death month data.
birth = new((/ nlat /), string)
death = new((/ nlat /), string)
w_array = tocharacter("W")
w_char = w_array(0)
dot_array = tocharacter(".")
dot_char = dot_array(0)
birth = where(ecmask, dot_char, w_char)
death = where(ecmask, dot_char, w_char)
m_array = tocharacter("123456789ABCD")

;; Determine ice birth and death month for each grid point and encode.
icethreshold = 0.01
do ilat = 0, nlat - 1
  do ilon = 0, nlon - 1
    if (.not. ismissing(eciceconc(0, ilat, ilon))) then
      hasice = eciceconc(:, ilat, ilon) .gt. icethreshold
      if (all(hasice)) then
        tmp = stringtocharacter(birth(ilat))
        tmp(ilon) = m_array(0)
        birth(ilat) = charactertostring(tmp)
        tmp = stringtocharacter(death(ilat))
        tmp(ilon) = m_array(12)
        death(ilat) = charactertostring(tmp)
      else
        if (any(hasice)) then
          if (hasice(0)) then
            if (isdefined("notices")) then
              delete(notices)
            end if
            notices = ind(.not. hasice)
            tmp = stringtocharacter(birth(ilat))
            tmp(ilon) = m_array(notices(dimsizes(notices) - 1))
            birth(ilat) = charactertostring(tmp)
            tmp = stringtocharacter(death(ilat))
            tmp(ilon) = m_array(notices(0) - 1)
            death(ilat) = charactertostring(tmp)
          else
            if (isdefined("ices")) then
              delete(ices)
            end if
            ices = ind(hasice)
            tmp = stringtocharacter(birth(ilat))
            tmp(ilon) = m_array(ices(0) - 1)
            birth(ilat) = charactertostring(tmp)
            tmp = stringtocharacter(death(ilat))
            tmp(ilon) = m_array(ices(dimsizes(ices) - 1) - 1)
            death(ilat) = charactertostring(tmp)
          end if
        end if
      end if
      btmp = stringtocharacter(birth(ilat))
      b = btmp(ilon)
      dtmp = stringtocharacter(death(ilat))
      d = dtmp(ilon)
      if (b .eq. d) then
        btmp(ilon) = w_char
        birth(ilat) = charactertostring(btmp)
        dtmp(ilon) = w_char
        death(ilat) = charactertostring(dtmp)
      end if
    end if
  end do
end do

out = new((/ 2 * nlat /), string)
do i = 0, nlat - 1
  row = sprinti("%2d", nlat - 1 - i)
  out(i) = str_concat((/ "  ", row, " ", birth(nlat - 1 - i) /))
  out(i + nlat) = str_concat((/ "  ", row, " ", death(nlat - 1 - i) /))
end do

asciiwrite("../inidata/icemask.dat", out)
