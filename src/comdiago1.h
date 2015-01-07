c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comdiago1.h                                                   
c *** Contents: Common declarations for diagnostics of ocean
c-----------------------------------------------------------------------
      real*8 s1hefx(nlat,nlon),s2hefx(nlat,nlon),
     &       s1safx(nlat,nlon),s2safx(nlat,nlon),
     &       s1brine(nlat,nlon),s2brine(nlat,nlon),
     &       s1fbtfx(nlat,nlon),s2fbtfx(nlat,nlon),
     &       s1hic(nlat,nlon),s2hic(nlat,nlon),
     &       s1hsn(nlat,nlon),s2hsn(nlat,nlon),
     &       s1tijs(nlat,nlon),s2tijs(nlat,nlon)


      common /meanocx1/s1hefx,s2hefx,
     &       s1safx,s2safx, s1brine,
     &       s2brine, s1fbtfx,s2fbtfx,
     &       s1hic,s2hic, s1hsn,s2hsn, 
     &       s1tijs,s2tijs

