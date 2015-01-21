load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

winfiles = (/ "ecbilt.nc", "um-tracer.nc" /)
oinfiles = (/ "ecbilt.nc", "um-offset.nc" /)

begin
  wks = gsn_open_wks("ps", "wind-plots")
  gsn_define_colormap(wks, "BlueWhiteOrangeRed")
  gsn_reverse_colormap(wks)
  NhlSetColor(wks, 253,   0/255.0, 224/255.0,   0/255.0)
  NhlSetColor(wks, 254, 224/255.0,   0/255.0, 224/255.0)
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
  res@gsnSpreadColorEnd = 252
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 0.01
  res@cnMinLevelValF = -0.1
  res@cnMaxLevelValF = 0.1
  res@cnLineLabelsOn = False
  res@cnFillOn = True
  res@cnLinesOn = False

  res@mpLandFillColor = 255

  res@lbLabelBarOn = False

  res@gsnDraw = False
  res@gsnFrame = False

  res@tmXTOn = False
  res@tmXBOn = False
  res@tmYLOn = False
  res@tmYROn = False

  vres1 = True
  vres1@tiYAxisString = ""
  vres1@tiXAxisString = ""
  vres1@gsnAddCyclic = False
  vres1@gsnLeftString = ""
  vres1@gsnRightString = ""
  vres1@gsnFrame = False
  vres1@gsnDraw = False
  vres1@vcVectorDrawOrder = "PostDraw"
  vres1@vcGlyphStyle = "FillArrow"
  vres1@vcFillArrowFillColor = 253
  vres1@vcFillArrowEdgeColor = "black"
  vres1@vcRefLengthF = 0.02
  vres1@vcRefMagnitudeF = 10
  vres1@vcRefAnnoString1 = "10 m/s"
  vres1@vcRefAnnoString2 = "Lower level wind"
  vres1@vcRefAnnoOn = False
  vres1@vcRefAnnoPerimOn = False
  vres1@vcRefAnnoFontHeightF = 0.015
  vres1@vcRefAnnoOrthogonalPosF = 0.15
  vres1@vcRefAnnoParallelPosF = 0.975

  vres2 = True
  vres2@tiYAxisString = ""
  vres2@tiXAxisString = ""
  vres2@gsnAddCyclic = False
  vres2@gsnLeftString = ""
  vres2@gsnRightString = ""
  vres2@gsnFrame = False
  vres2@gsnDraw = False
  vres2@vcVectorDrawOrder = "PostDraw"
  vres2@vcGlyphStyle = "FillArrow"
  vres2@vcFillArrowFillColor = 254
  vres2@vcFillArrowEdgeColor = "black"
  vres2@vcRefLengthF = 0.02
  vres2@vcRefMagnitudeF = 20
  vres2@vcRefAnnoString1 = "20 m/s"
  vres2@vcRefAnnoString2 = "Upper level wind"
  vres2@vcRefAnnoOn = False
  vres2@vcRefAnnoPerimOn = False
  vres2@vcRefAnnoFontHeightF = 0.015
  vres2@vcRefAnnoOrthogonalPosF = 0.15
  vres2@vcRefAnnoParallelPosF = 0.7

  do row = 0, 3
    do col = 0, 1
      wfp = addfile(winfiles(col), "r")
      ofp = addfile(oinfiles(col), "r")
      omega = ofp->omega
      uup = wfp->uupper
      vup = wfp->vupper
      ulo = wfp->ulower
      vlo = wfp->vlower

      if (col .eq. 0) then
        ws = 4
      else
        ws = 7
      end if
      if (row .eq. 3 .and. col .eq. 1) then
        vres1@vcRefAnnoOn = True
        vres2@vcRefAnnoOn = True
      end if

      plots(plotidx) = gsn_csm_contour_map_ce(wks, omega(row, :, :), res)
      v1 = gsn_csm_vector(wks, ulo(row, ::ws,::ws), \
                               vlo(row, ::ws,::ws), vres1)
      v2 = gsn_csm_vector(wks, uup(row, 1::ws,1::ws), \
                               vup(row, 1::ws,1::ws), vres2)
      overlay(plots(plotidx), v1)
      overlay(plots(plotidx), v2)

      panel_labels(plotidx) = col_labels(col) + " (" + row_labels(row) + ")"
      plotidx = plotidx + 1
      delete(omega)
      delete(uup)
      delete(vup)
      delete(ulo)
      delete(vlo)
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
  lbw = 0.55
  lbh = 0.05
  lbres@vpWidthF = lbw
  lbres@vpHeightF = lbh
  lbres@lbMonoFillPattern = True
  lbres@lbLabelFontHeightF = 0.015
  lbres@lbLabelAlignment = "InteriorEdges"
  lbres@lbPerimOn = False
  lbres@lbOrientation = "Horizontal"
  lbres@lbFillColors = lbcolors
  lbres@lbTitleString = "Vertical pressure velocity (Pa/s)"
  lbres@lbTitleFontHeightF = 0.015
  lbres@lbTitleOffsetF = 0.2
  lbres@lbTitlePosition = "Bottom"
  lbres@lbLabelStride = 2
  lblabels = new(dimsizes(lblevels), string)
  lblabels = lblevels

  lbx = (0.75 - lbw) / 2.0
  lby = 0.08
  gsn_labelbar_ndc(wks, dimsizes(lbcolors), lblabels, lbx, lby, lbres)

  frame(wks)
end
