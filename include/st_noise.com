c-- parameters needed for introducing a stochastic noise to LOVECLIM :
c--  nmodesMax_st is the maximum value of number of modes used in EOF 
c--  ntMax_st is the maximum value of number of time steps used in EOF
c--  rhocp_st is \rho_cp in the heat flux equation
c--  h_st is h in the heat flux equation
      integer nmodesMax_st,ntMax_st
      parameter ( nmodesMax_st = 2000 , ntMax_st = 2000 )
      real*8 rhocp_st,h_st
      parameter ( rhocp_st=4*1.0e+6 , h_st =50 )

