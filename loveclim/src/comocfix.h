!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** File:     comocfix.h
! *** Contents: Common declarations for fixed ocean of ECbilt
!-----------------------------------------------------------------------
      real*8  sstday(nlat,nlon,360)
      integer lsicebirth(nlat,nlon),lsicedeath(nlat,nlon)
      integer lseaice(nlat, nlon)

      common /rocef/ sstday
      common /iocef/ lsicebirth,lsicedeath,lseaice
