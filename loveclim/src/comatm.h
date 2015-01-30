!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** File:     comatm.h
! *** Contents: General Parameter and Common declarations for ECbilt
!-----------------------------------------------------------------------
! *** PARAMETERS
!     nm  :   the truncation is of type T(riangular) nm.
!     nlon:   number of longitude points of the Gaussian grid
!     nlat:   number of latitude  points of the Gaussian grid
!     nvl :   number of vorticity levels in the vertical
!             (should be set to 3)
!     ntl :   number of temperature levels in the vertical
!             (equal to nvl-1)
!     nsh :   half of nsh2
!     nsh2:   number of coefficients needed to define one level of the
!             T nm model
!     ngp:    number of grid points of the Gaussian grid
!
! *** COMMON  /rvari/ pi,dp,om,rgas,grav,radius
!     pi :     value of pi
!     fzero:   value of f at 45 degrees
!     dp :     layer thicknesses [Pa]
!     om :     angular velocity of earth [rad/s]
!     rgas :   gas constant
!     grav :   gravity acceleration [m/s^2]
!     radius:  radius of earth [m]
!
! *** COMMON  /ipar/  iadyn,iaphys,initdate,initfield
!     iadyn :  with (1) atmospheric dynamics or without (0)
!     iaphys:  with (1) atmospheric physics or without (0)
!     iaphys:  with (1) atmospheric physics or without (0)
!     initfield:  initiqlisation form rest(0) or from file (1)
!     initdate: date of the beginning of the WHOLE experiment)
!               (actually it is the difference between the starting point
!                of the experiment and initdate)
!
! *** COMMON  /ggrid/ phi,dphi,dlab,darea,tarea,cosfi,cosfid,sinfi,
!                     sinfid,tanfi
!     phi :    Gauss points in radians
!     darea:   area [m**2] of gridboxes as a function of latitude
!     tarea:   total area of earth [m**2]
!     cosfi:   cosine of phi
!     sinfi:   sine of phi
!     tanfi:   tangens of phi
!
! *** COMMON  /ctstep/ dt,dtime,dtt
!     dt :     timestep in fraction of one day
!     dtime :  timestep in seconds
!     dtt :    timestep in non-dimensional units
!
! *** COMMON  /plev/  plevel,tlevel,rlogtl12
!     plevel:   contains value of the pressure at the nvl levels [Pa]
!     tlevel:   contains value of the pressure at the ntl levels [Pa]
!     rlogtl12: contains log(tlevel(1)/tlevel(2))
!-----------------------------------------------------------------------

      integer   nm,nsh,nlon,nlat,nsh2,ngp,nvl,ntl
      parameter ( nm=21, nlon=64, nlat=32, nvl=3, ntl=nvl-1, &
     &            nsh=((nm+1)*(nm+2))/2, nsh2=2*nsh, ngp=nlon*nlat)

      integer iadyn,iaphys,initdate,initfield
      real*8  pi,fzero,dp,om,rgas,grav,radius
      real*8  phi(nlat)
      real*8  cosfi(nlat),sinfi(nlat),tanfi(nlat)
      real*8  dt,dtime,dtt,rdtime
      real*8  plevel(nvl),tlevel(ntl),rlogtl12,p0,tarea,tareas
      real*8  alogtl12,alogtl1pl2,alogpl2tl2,darea(nlat),dareas(nlat)


      common /rvari/ pi,fzero,dp,om,rgas,grav,radius
      common /ipar/  iadyn,iaphys,initdate,initfield
      common /ggrid/ phi,cosfi,sinfi,tanfi,dareas,darea,tarea,tareas
      common /ctstep/ dt,dtime,dtt,rdtime
      common /plev/  plevel,tlevel,rlogtl12,p0,alogtl12,alogtl1pl2, &
     &               alogpl2tl2
