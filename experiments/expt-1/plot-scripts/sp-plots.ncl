load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

infiles = (/ "ecbilt.nc", "um-tracer.nc" /)

begin
  wks = gsn_open_wks("ps", "sp-plots")
  gsn_define_colormap(wks, "BlueRed")
  NhlSetColor(wks, 255, 192/255.0, 192/255.0, 192/255.0)

  plots = new(8, graphic)
  panel_labels = new(8, string)
  plotidx = 0

  row_labels = (/ "DJF", "MAM", "JJA", "SON" /)
  col_labels = (/ "ECBILT", "UM:tcszd" /)

  res = True

  res@tiMainString = ""
  res@tiYAxisString = ""
  res@tiXAxisString = ""
  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLineLabelsOn = True
  res@cnLineLabelPlacementMode = "Computed"
  res@cnInfoLabelOn = False
  res@cnFillOn = False
  res@cnLineThicknessF = 0.125

  res@mpLandFillColor = 255

  res@lbLabelBarOn = False

  res@gsnDraw = False
  res@gsnFrame = False

  res@tmXTOn = False
  res@tmXBOn = False
  res@tmYLOn = False
  res@tmYROn = False

  do row = 0, 3
    do col = 0, 1
      fp = addfile(infiles(col), "r")
      sp = fp->sp
      sp = sp / 100.0

      if (col .eq. 0) then
        res@cnLevelSpacingF = 20
        res@cnMinLevelValF = 500
        res@cnMaxLevelValF = 1100
      else
        res@cnLevelSpacingF = 5
        res@cnMinLevelValF = 500
        res@cnMaxLevelValF = 1100
      end if

      plots(plotidx) = gsn_csm_contour_map_ce(wks, sp(row, :, :), res)
      panel_labels(plotidx) = col_labels(col) + " (" + row_labels(row) + ")"
      plotidx = plotidx + 1
      delete(sp)
    end do
  end do

  pw = 0.406212
  ph = 0.203106
  pgapx = 0.015
  pgapy = 0.015
  px1 = 0.086288
  py1 = 0.98
  panelx = new(8, float)
  panely = new(8, float)
  pidx = 0
  do r = 0, 3
    do c = 0, 1
      panelx(pidx) = px1 + (pw + pgapx) * c
      panely(pidx) = py1 - (ph + pgapy) * r
      pidx = pidx + 1
    end do
  end do

  pres = True
;  pres@gsnPanelDebug = True
  pres@gsnPanelXF = panelx
  pres@gsnPanelYF = panely
  pres@gsnPanelFigureStrings = panel_labels
  pres@gsnPanelFigureStringsPerimOn = False
  pres@gsnPanelFigureStringsPerimSpaceF = 0.125
  pres@gsnPanelFigureStringsFontHeightF = 0.01
  pres@amJust = "BottomLeft"
  gsn_panel(wks, plots, (/ 4, 2 /), pres)
end
