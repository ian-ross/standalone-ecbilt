ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       original Green's functions  for T(both Ta and Ts), q, c (cloud 
c    cover) and pt (cloud top).
c
c    greenfda-.....Green's functions for downward fluxes
c    greenfua-.....Green's functions for upward fluxes originated at
c                  the atmosphere
c    greenfus-.....Green's functions for upward fluxes originated at
c                  the surface
c    -t............temperature (both atmospheric and surface)
c    -q............moisture
c    -c............cloud cover
c    -pt...........cloud top
c    -ghg..........greenhouse gases
c    greenfdat(x,y,z), where
c             x....flux layers
c             y....perturbation layers
c             z....cloud cover type
c               definitions see Chou and Neelin 1995, and note: the number
c               is different from the reference. Here is an example wich
c               was used in the calculation of Chou and Neelin 1995. 
c               ntype=0, for clear sky
c               ntype=1, for ECHAM cloudy sky
c
c     greenfghg(x,g,z) where x&z as above, g is number of GHG, see 'gpar.cbl'
c
	REAL reffu(nfluxl,0:ntype),reffd(nfluxl,0:ntype)
      REAL greenfdat(nfluxl,nlm,0:ntype)
     &,  greenfuat(nfluxl,nlm,0:ntype)
     &,  greenfuts(nfluxl,0:ntype,4),greenfdts(0:ntype,4)
     &,  greenfdaq(nfluxl,nlm,0:ntype)
     &,  greenfuaq(nfluxl,nlm,0:ntype)
     &,  greenfuaqcol(nfluxl,0:ntype)
     &,  greenfdaqcol(nfluxl,0:ntype)
	REAL greenfughg(nfluxl,nghg,0:ntype)
     &,     greenfdghg(nfluxl,nghg,0:ntype)
     &,     greenfuco2(nfluxl,0:ntype,2),greenfdco2(nfluxl,0:ntype,2)
     &,     greenfuch4(nfluxl,0:ntype,2),greenfdch4(nfluxl,0:ntype,2)
     &,     greenfun2o(nfluxl,0:ntype,2),greenfdn2o(nfluxl,0:ntype,2)
     &,     greenfuo3(nfluxl,nlm,0:ntype),greenfdo3(nfluxl,nlm,0:ntype)
      common /green0/ reffu,reffd,greenfdat,greenfuat,greenfuts
     &,  greenfdts,greenfdaq,greenfuaq,greenfuaqcol,greenfdaqcol
     &,  greenfughg,greenfdghg
     &,  greenfuco2,greenfdco2
     &,  greenfuch4,greenfdch4
     &,  greenfun2o,greenfdn2o
     &,  greenfuo3,greenfdo3
c
