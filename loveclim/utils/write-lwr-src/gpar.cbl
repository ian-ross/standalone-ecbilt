!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Parameters for linear and weakly nonlnear shcemes
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!      nlm......number of model layers
!               nlm=1 is the highest layer in the model
!               nlm=20 is the lowest layer in the model
!      profnlm..number of layers of original reference profile,
!               to be interpolated
!      nfluxl...number of flux levels for which reference fluxes and
!               Greens' functions are output
!      fluxl....for each flux level the pressure level in hPa
!               value -1 indicates surface level.
!      o3profnlm..number of layers of original reference profile for ozone,
!                 to be interpolated
!      ntype....cloud type
!               definitions see Chou and Neelin 1995, and note: the number
!               is different from the reference.
!               ntype=0, for clear sky
!               ntype=1, for ECHAM4 cloudy sky
!      latstep..step taken in latitude loops when calculating zonal mean
!               responses
!      lonstep..step in longitude loops
!      nvar.....number of variables,
!               nvar=1, for T (atmospheric temperature) (K)
!               nvar=2, for q (Kg/Kg)
!               nvar=3, for cloud cover (%)
!               nvar=4, for cloud top (mb)
!               nvar=5, for Ts (K)
!               nvar=6, for O3 (kg/kg)
!               nvar=7, for CO2 (ppmv)
!               nvar=8, for CH4 (ppbv)
!               nvar=9, for N2O(ppbv)
!               nvar=10, for CFC-11 (pptv)
!               nvar=11, for CFC-12 (pptv)
!               nvar=12, for CFC-113 (pptv)
!               nvar=13, for CFC-114 (pptv)
!               nvar=14, for CFC-115 (pptv)
!               nvar=15, for HCFC-22 (pptv)
!               nvar=16, for HCFC-123 (pptv)
!               nvar=17, for HCFC-124 (pptv)
!               nvar=18, for HCFC-125 (pptv)
!               nvar=19, for HCFC-134a (pptv)
!               nvar=20, for HCFC-141b (pptv)
!               nvar=21, for HCFC-142b (pptv)
!               nvar=22, for HCFC-143a (pptv)
!               nvar=23, for HCFC-152a (pptv)
!               nvar=24, for CTC (pptv)
!               nvar=25, for MCF (pptv)
!      nghg.....number of greenhouse gases, excluding H2O and O3
!               nghg=1, for CO2 (ppmv)
!               nghg=2, for CH4 (ppbv)
!               nghg=3, for N2O(ppbv)
!               nghg=4, for CFC-11 (pptv)
!               nghg=5, for CFC-12 (pptv)
!               nghg=6, for CFC-113 (pptv)
!               nghg=7, for CFC-114 (pptv)
!               nghg=8, for CFC-115 (pptv)
!               nghg=9, for HCFC-22 (pptv)
!               nghg=10, for HCFC-123 (pptv)
!               nghg=11, for HCFC-124 (pptv)
!               nghg=12, for HCFC-125 (pptv)
!               nghg=13, for HCFC-134a (pptv)
!               nghg=14, for HCFC-141b (pptv)
!               nghg=15, for HCFC-142b (pptv)
!               nghg=16, for HCFC-143a (pptv)
!               nghg=17, for HCFC-152a (pptv)
!               nghg=18, for CTC (pptv)
!               nghg=19, for MCF (pptv)
!      g........acceleration due to gravity at surface of earth (m/s**2)
!      Cp.......specific heat at constant pressure (J/K/Kg)
!      day.....one day (second)
!      dp.......transferring pressure from mb to Pa.
!      pl.......at the model layer (mb)
!      pl2......at the interface between two layers (mb)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  INTEGER, PARAMETER :: nlm=20, ntype=1, nfluxl=4, nvar=25
  INTEGER, PARAMETER :: profnlm=17, o3profnlm=31
  INTEGER, PARAMETER :: NLMP1=NLM+1, NLMP2=NLM+2, NLMP3=NLM+3
  INTEGER, PARAMETER :: MAXBOX=6596
  INTEGER, PARAMETER :: latstep=1, lonstep=2
  REAL, PARAMETER :: dp=100.0, day=86400.0, sfcem=0.96, QFAC=0.85
  INTEGER, PARAMETER :: nghg=19
  CHARACTER*56, PARAMETER :: wdir= &
       & '                            /rivm/users/imagerun/greens/'
!         12345678901234567890123456789012345678901234567890123456
