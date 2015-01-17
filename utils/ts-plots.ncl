load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

infiles = (/ "ecbilt.nc", "um-tracer.nc" /)

begin
  wks = gsn_open_wks("ps", "ts-plots")
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
  res@gsnSpreadColors = True
  res@gsnSpreadColorEnd = 253
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 10
  res@cnMinLevelValF = 210
  res@cnMaxLevelValF = 330
  res@cnLineLabelsOn = False
  res@cnFillOn = True
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
      ts = fp->ts

      plots(plotidx) = gsn_csm_contour_map_ce(wks, ts(row, :, :), res)
      panel_labels(plotidx) = col_labels(col) + " (" + row_labels(row) + ")"
      plotidx = plotidx + 1
      delete(ts)
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
  pres@gsnFrame = False
  pres@gsnPanelFigureStrings = panel_labels
  pres@gsnPanelFigureStringsPerimOn = False
  pres@gsnPanelFigureStringsPerimSpaceF = 0.125
  pres@gsnPanelFigureStringsFontHeightF = 0.01
  pres@amJust = "BottomLeft"
  gsn_panel(wks, plots, (/ 4, 2 /), pres)

  contour = get_overlay_plot(plots(0), "contourPlotClass", 0)
  getvalues contour
    "cnLevels": lblevels
    "cnFillColors": lbcolors
  end getvalues

  do idx = 0, dimsizes(lblevels) - 1
    if (abs(lblevels(idx)) .lt. 0.001) then
      lblevels(idx) = 0.0
    end if
  end do

  lbres = True
  lbw = 0.7
  lbh = 0.05
  lbres@vpWidthF = lbw
  lbres@vpHeightF = lbh
  lbres@lbMonoFillPattern = True
  lbres@lbLabelFontHeightF = 0.015
  lbres@lbLabelAlignment = "InteriorEdges"
  lbres@lbPerimOn = False
  lbres@lbOrientation = "Horizontal"
  lbres@lbFillColors = lbcolors
  lbres@lbTitleString = "Surface temperature (K)"
  lbres@lbTitleFontHeightF = 0.015
  lbres@lbTitleOffsetF = 0.2
  lbres@lbTitlePosition = "Bottom"
  lbres@lbLabelStride = 2
  lblabels = new(dimsizes(lblevels), string)
  lblabels = lblevels

  lbx = (1.0 - lbw) / 2.0
  lby = 0.08
  gsn_labelbar_ndc(wks, dimsizes(lbcolors), lblabels, lbx, lby, lbres)

  frame(wks)
end
