!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** File:     comato.h
! *** Contents: Common declarations for driving the ocean of ECbilt
! ***           with atmospheric data during spin-up phase
!-----------------------------------------------------------------------

      integer  nocfile
      integer  nbsatfor,nafyear
      real*4   otempsg(nlat,nlon),odlrads(nlat,nlon)
      real*4   ouv10(nlat,nlon),owinstu(nlat,nlon),owinstv(nlat,nlon)
      real*4   oq10(nlat,nlon),otorain(nlat,nlon),orunofo(nlat,nlon)
      real*4   oevap(nlat,nlon)

      common /ato/  otempsg,odlrads,ouv10,owinstu,owinstv, &
     &              oq10,otorain,orunofo,oevap
      common /iato/ nocfile,nbsatfor,nafyear
