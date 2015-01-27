c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comland.h                                                    
c *** Contents: Common declarations land parameterizations of ECbilt
c-----------------------------------------------------------------------

        integer   nbasins
        parameter (nbasins=24)

	real*8    bmoisg(nlat,nlon),bmoism	
	real*8    dsnow(nlat,nlon),dsnm,tland(nlat,nlon)
	integer   ilabas(nlat,nlon),iocbas(nlat,nlon)
	real*8    runofl(nlat,nlon),runofo(nlat,nlon),arocbas(nbasins)

	common /bmois/ bmoisg,bmoism,tland	
	common /csnow/ dsnow,dsnm
	common /iruno/ ilabas,iocbas
	common /runol/ runofl,runofo,arocbas

