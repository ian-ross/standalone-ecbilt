!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       original Green's functions  for T(both Ta and Ts), q, c (cloud
!    cover) and pt (cloud top).
!
!    greenfda-.....Green's functions for downward fluxes
!    greenfua-.....Green's functions for upward fluxes originated at
!                  the atmosphere
!    greenfus-.....Green's functions for upward fluxes originated at
!                  the surface
!    -t............temperature (both atmospheric and surface)
!    -q............moisture
!    -c............cloud cover
!    -pt...........cloud top
!    -ghg..........greenhouse gases
!    greenfdat(x,y,z), where
!             x....flux layers
!             y....perturbation layers
!             z....cloud cover type
!               definitions see Chou and Neelin 1995, and note: the number
!               is different from the reference. Here is an example wich
!               was used in the calculation of Chou and Neelin 1995.
!               ntype=0, for clear sky
!               ntype=1, for ECHAM cloudy sky
!
!     greenfghg(x,g,z) where x&z as above, g is number of GHG, see 'gpar.cbl'
!
  REAL reffu(nfluxl,0:ntype),reffd(nfluxl,0:ntype)
  REAL greenfdat(nfluxl,nlm,0:ntype), &
       & greenfuat(nfluxl,nlm,0:ntype), &
       & greenfuts(nfluxl,0:ntype,4),greenfdts(0:ntype,4), &
       & greenfdaq(nfluxl,nlm,0:ntype), &
       & greenfuaq(nfluxl,nlm,0:ntype), &
       & greenfuaqcol(nfluxl,0:ntype), &
       & greenfdaqcol(nfluxl,0:ntype)
  REAL greenfughg(nfluxl,nghg,0:ntype), &
       & greenfdghg(nfluxl,nghg,0:ntype), &
       & greenfuco2(nfluxl,0:ntype,2),greenfdco2(nfluxl,0:ntype,2), &
       & greenfuch4(nfluxl,0:ntype,2),greenfdch4(nfluxl,0:ntype,2), &
       & greenfun2o(nfluxl,0:ntype,2),greenfdn2o(nfluxl,0:ntype,2), &
       & greenfuo3(nfluxl,nlm,0:ntype),greenfdo3(nfluxl,nlm,0:ntype)
  common /green0/ reffu,reffd,greenfdat,greenfuat,greenfuts, &
       & greenfdts,greenfdaq,greenfuaq,greenfuaqcol,greenfdaqcol, &
       & greenfughg,greenfdghg,greenfuco2,greenfdco2,greenfuch4, &
       & greenfdch4,greenfun2o,greenfdn2o,greenfuo3,greenfdo3
