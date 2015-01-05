ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Parameters for linear and weakly nonlnear shcemes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      nlm......number of model layers
c               nlm=1 is the highest layer in the model
c               nlm=20 is the lowest layer in the model
c	 profnlm..number of layers of original reference profile,
c		    to be interpolated
c      nfluxl...number of flux levels for which reference fluxes and
c	          Greens' functions are output
c      fluxl....for each flux level the pressure level in hPa
c		    value -1 indicates surface level.
c	 o3profnlm..number of layers of original reference profile for ozone,
c		    to be interpolated
c      ntype....cloud type
c               definitions see Chou and Neelin 1995, and note: the number
c               is different from the reference.
c               ntype=0, for clear sky
c               ntype=1, for ECHAM4 cloudy sky
c      latstep..step taken in latitude loops when calculating zonal mean
c               responses
c      lonstep..step in longitude loops
c      nvar.....number of variables, 
c               nvar=1, for T (atmospheric temperature) (K)
c               nvar=2, for q (Kg/Kg)
c               nvar=3, for cloud cover (%)
c               nvar=4, for cloud top (mb)
c               nvar=5, for Ts (K)
c               nvar=6, for O3 (kg/kg)
c               nvar=7, for CO2 (ppmv)
c               nvar=8, for CH4 (ppbv)
c               nvar=9, for N2O(ppbv)
c               nvar=10, for CFC-11 (pptv)
c               nvar=11, for CFC-12 (pptv)
c               nvar=12, for CFC-113 (pptv)
c               nvar=13, for CFC-114 (pptv)
c               nvar=14, for CFC-115 (pptv)
c               nvar=15, for HCFC-22 (pptv)
c               nvar=16, for HCFC-123 (pptv)
c               nvar=17, for HCFC-124 (pptv)
c               nvar=18, for HCFC-125 (pptv)
c               nvar=19, for HCFC-134a (pptv)
c               nvar=20, for HCFC-141b (pptv)
c               nvar=21, for HCFC-142b (pptv)
c               nvar=22, for HCFC-143a (pptv)
c               nvar=23, for HCFC-152a (pptv)
c               nvar=24, for CTC (pptv)
c               nvar=25, for MCF (pptv)
c	 nghg.....number of greenhouse gases, excluding H2O and O3
c		    nghg=1, for CO2 (ppmv)
c               nghg=2, for CH4 (ppbv)
c               nghg=3, for N2O(ppbv)
c               nghg=4, for CFC-11 (pptv)
c               nghg=5, for CFC-12 (pptv)
c               nghg=6, for CFC-113 (pptv)
c               nghg=7, for CFC-114 (pptv)
c               nghg=8, for CFC-115 (pptv)
c               nghg=9, for HCFC-22 (pptv)
c               nghg=10, for HCFC-123 (pptv)
c               nghg=11, for HCFC-124 (pptv)
c               nghg=12, for HCFC-125 (pptv)
c               nghg=13, for HCFC-134a (pptv)
c               nghg=14, for HCFC-141b (pptv)
c               nghg=15, for HCFC-142b (pptv)
c               nghg=16, for HCFC-143a (pptv)
c               nghg=17, for HCFC-152a (pptv)
c               nghg=18, for CTC (pptv)
c               nghg=19, for MCF (pptv)
c      g........acceleration due to gravity at surface of earth (m/s**2)   
c      Cp.......specific heat at constant pressure (J/K/Kg)
c      day.....one day (second)
c      dp.......transferring pressure from mb to Pa.
c      pl.......at the model layer (mb)
c      pl2......at the interface between two layers (mb)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	INTEGER nlm,ntype,nfluxl,nvar,nghg,profnlm,o3profnlm,NLMp1,NLMP2,
     &        NLMP3,MAXBOX
	REAL dp,day,sfcem,QFAC
	CHARACTER*56 wdir
      parameter(nlm=20,ntype=1,nfluxl=4,nvar=25,profnlm=17,o3profnlm=31)
      PARAMETER ( NLMP1=NLM+1,    NLMP2=NLM+2,    NLMP3=NLM+3)
      PARAMETER(MAXBOX=6596)
      INTEGER latstep,lonstep
      parameter(latstep=1,lonstep=2)
      parameter(dp=100., day=86400., sfcem=0.96, QFAC=0.85)
	PARAMETER(nghg=19)
	PARAMETER(wdir=
     &'                            /rivm/users/imagerun/greens/')
c      12345678901234567890123456789012345678901234567890123456
