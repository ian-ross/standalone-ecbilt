












c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comatm.h                                                    
c *** Contents: General Parameter and Common declarations for ECbilt
c-----------------------------------------------------------------------
c *** PARAMETERS
c     nm  :   the truncation is of type T(riangular) nm. 
c     nlon:   number of longitude points of the Gaussian grid
c     nlat:   number of latitude  points of the Gaussian grid
c     nvl :   number of vorticity levels in the vertical 
c             (should be set to 3)
c     ntl :   number of temperature levels in the vertical 
c             (equal to nvl-1)
c     nsh :   half of nsh2
c     nsh2:   number of coefficients needed to define one level of the 
c             T nm model
c     ngp:    number of grid points of the Gaussian grid
c 
c *** COMMON  /rvari/ pi,dp,om,rgas,grav,radius
c     pi :     value of pi
c     fzero:   value of f at 45 degrees
c     dp :     layer thicknesses [Pa]
c     om :     angular velocity of earth [rad/s]
c     rgas :   gas constant
c     grav :   gravity acceleration [m/s^2]
c     radius:  radius of earth [m]
c
c *** COMMON  /ipar/  iadyn,iaphys,initdate,initfield
c     iadyn :  with (1) atmospheric dynamics or without (0)
c     iaphys:  with (1) atmospheric physics or without (0)
c     iaphys:  with (1) atmospheric physics or without (0)
c     initfield:  initiqlisation form rest(0) or from file (1)
c     initdate: date of the beginning of the WHOLE experiment)
c               (actually it is the difference between the starting point
c                of the experiment and initdate)
c
c *** COMMON  /ggrid/ phi,dphi,dlab,darea,tarea,cosfi,cosfid,sinfi,
c                     sinfid,tanfi
c     phi :    Gauss points in radians
c     darea:   area [m**2] of gridboxes as a function of latitude
c     tarea:   total area of earth [m**2]
c     cosfi:   cosine of phi
c     sinfi:   sine of phi
c     tanfi:   tangens of phi
c
c *** COMMON  /ctstep/ dt,dtime,dtt
c     dt :     timestep in fraction of one day
c     dtime :  timestep in seconds
c     dtt :    timestep in non-dimensional units
c
c *** COMMON  /plev/  plevel,tlevel,rlogtl12
c     plevel:   contains value of the pressure at the nvl levels [Pa]
c     tlevel:   contains value of the pressure at the ntl levels [Pa]
c     rlogtl12: contains log(tlevel(1)/tlevel(2))
c-----------------------------------------------------------------------

      integer   nm,nsh,nlon,nlat,nsh2,ngp,nvl,ntl
      parameter ( nm=21, nlon=64, nlat=32, nvl=3, ntl=nvl-1,
     *            nsh=((nm+1)*(nm+2))/2, nsh2=2*nsh, ngp=nlon*nlat)

      integer iadyn,iaphys,initdate,initfield
      real*8  pi,fzero,dp,om,rgas,grav,radius
      real*8  phi(nlat)
      real*8  cosfi(nlat),sinfi(nlat),tanfi(nlat)
      real*8  dt,dtime,dtt,rdtime
      real*8  plevel(nvl),tlevel(ntl),rlogtl12,p0,tarea,tareas
      real*8  alogtl12,alogtl1pl2,alogpl2tl2,darea(nlat),dareas(nlat)


      common /ec_rvari/ pi,fzero,dp,om,rgas,grav,radius
      common /ec_ipar/  iadyn,iaphys,initdate,initfield  
      common /ec_ggrid/ phi,cosfi,sinfi,tanfi,dareas,darea,tarea,tareas
      common /ec_ctstep/ dt,dtime,dtt,rdtime
      common /ec_plev/  plevel,tlevel,rlogtl12,p0,alogtl12,alogtl1pl2,
     *               alogpl2tl2
