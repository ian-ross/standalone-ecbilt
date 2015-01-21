load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

infiles = (/ "ecbilt-orog.nc", "um-orog.nc" /)

begin
  wks = gsn_open_wks("ps", "orog-plots")
  gsn_define_colormap(wks, "OceanLakeLandSnow")

  plots = new(2, graphic)
  panel_labels = new(2, string)
  plotidx = 0

  col_labels = (/ "ECBILT", "UM:tcszd" /)

  res = True

  res@tiMainString = ""
  res@tiYAxisString = ""
  res@tiXAxisString = ""
  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@gsnSpreadColors = True
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 500
  res@cnMinLevelValF = 0
  res@cnMaxLevelValF = 5500
  res@cnLineLabelsOn = False
  res@cnLinesOn = False
  res@cnFillOn = True
  res@cnFillMode = "CellFill"

  res@lbLabelBarOn = False

  res@gsnDraw = False
  res@gsnFrame = False

  res@tmXTOn = False
  res@tmXBOn = False
  res@tmYLOn = False
  res@tmYROn = False
  res@mpOutlineBoundarySets = "NoBoundaries"
  res@mpFillBoundarySets = "NoBoundaries"

  do col = 0, 1
    fp = addfile(infiles(col), "r")
    orog = fp->orog

    plots(plotidx) = gsn_csm_contour_map_ce(wks, orog, res)
    panel_labels(plotidx) = col_labels(col)
    plotidx = plotidx + 1
    delete(orog)
  end do

  pres = True
  pres@gsnFrame = False
  pres@gsnPanelFigureStrings = panel_labels
  pres@gsnPanelFigureStringsPerimOn = False
  pres@gsnPanelFigureStringsPerimSpaceF = 0.125
  pres@gsnPanelFigureStringsFontHeightF = 0.01
  pres@amJust = "BottomLeft"
  gsn_panel(wks, plots, (/ 1, 2 /), pres)

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
  lbres@lbTitleString = "Orography (m)"
  lbres@lbTitleFontHeightF = 0.015
  lbres@lbTitleOffsetF = 0.2
  lbres@lbTitlePosition = "Bottom"
  lbres@lbLabelStride = 2
  lblabels = new(dimsizes(lblevels), string)
  lblabels = lblevels

  lbx = (1.0 - lbw) / 2.0
  lby = 0.34
  gsn_labelbar_ndc(wks, dimsizes(lbcolors), lblabels, lbx, lby, lbres)

  frame(wks)
end
