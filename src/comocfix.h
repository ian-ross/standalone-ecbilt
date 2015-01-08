c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comocfix.h
c *** Contents: Common declarations for fixed ocean of ECbilt
c-----------------------------------------------------------------------
        real*8   sstday(nlat,nlon,360)
        integer  lsicebirth(nlat,nlon),lsicedeath(nlat,nlon)

	common /rocef/ sstday
	common /iocef/ lsicebirth,lsicedeath
