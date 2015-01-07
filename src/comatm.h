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
c *** COMMON  /ipar/  iadyn,iaphys
c     iadyn :  with (1) atmospheric dynamics or without (0)
c     iaphys:  with (1) atmospheric physics or without (0)
c
c *** COMMON  /ggrid/ phi,dphi,dlab,darea,tarea,cosfi,cosfid,sinfi,
c                     sinfid,tanfi
c     phi :    Gauss points
c     dphi:    approximate distance between Gauss points
c     dlab:    distance between grid points in zonal direction
c     darea:   dphi*dlab
c     tarea:   total area of earth
c     cosfi:   cosine of phi
c     cosfid:  cosine of phi times darea
c     sinfi:   sine of phi
c     sinfid:  sine of phi times darea
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

      integer iadyn,iaphys
      real*8  pi,fzero,dp,om,rgas,grav,radius
      real*8  phi(nlat),dphi,dlab,darea,tarea
      real*8  cosfi(nlat),cosfid(nlat),sinfi(nlat),sinfid(nlat),
     *        tanfi(nlat)
      real*8  dt,dtime,dtt
      real*8  plevel(nvl),tlevel(ntl),rlogtl12,p0
      real*8  alogtl12,alogtl1pl2,alogpl2tl2


      common /rvari/ pi,fzero,dp,om,rgas,grav,radius
      common /ipar/  iadyn,iaphys
      common /ggrid/ phi,dphi,dlab,darea,tarea,cosfi,cosfid,sinfi,
     *               sinfid,tanfi
      common /ctstep/ dt,dtime,dtt
      common /plev/  plevel,tlevel,rlogtl12,p0,alogtl12,alogtl1pl2,
     *               alogpl2tl2
