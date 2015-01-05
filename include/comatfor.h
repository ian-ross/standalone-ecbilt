c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comato.h
c *** Contents: Common declarations for driving the ocean of ECbilt
c ***           with atmospheric data during spin-up phase
c-----------------------------------------------------------------------

      integer  nocfile
      integer  nbsatfor,nafyear
      real*8   otempsg(nlat,nlon),odlrads(nlat,nlon)
      real*8   ouv10(nlat,nlon),owinstu(nlat,nlon),owinstv(nlat,nlon)
      real*8   oq10(nlat,nlon),otorain(nlat,nlon),orunofo(nlat,nlon)
      real*8   oevap(nlat,nlon)

      common /ec_ato/  otempsg,odlrads,ouv10,owinstu,owinstv,
     &              oq10,otorain,orunofo,oevap
      common /ec_iato/ nocfile,nbsatfor,nafyear
