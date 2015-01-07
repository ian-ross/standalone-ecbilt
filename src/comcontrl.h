c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comcontrl.h
c *** contents: common declarations for control routines
c-----------------------------------------------------------------------

      real*8    bmoisgo(nlat,nlon),dsnowo(nlat,nlon)
      real*8    toraino(nlat,nlon),evapo(nlat,nlon)
      real*8    rmoisgo(nlat,nlon),runoflo(nlat,nlon)
      real*8    runofoo(nlat,nlon)

      common /mocntrl/ bmoisgo,dsnowo,toraino,evapo,
     *                 rmoisgo,runoflo,runofoo
