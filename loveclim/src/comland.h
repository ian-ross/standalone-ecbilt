!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** File:     comland.h
! *** Contents: Common declarations land parameterizations of ECbilt
!-----------------------------------------------------------------------

        integer   nbasins,nlat,nlon,mbasins
        real*8    epsl
        parameter (mbasins=26,nlat=32,nlon=64,epsl=1e-10)

        real*8    bmoisg(nlat,nlon),bmoism(nlat,nlon)
        real*8    lhcap,rlhcap,dtland,rdtland
        real*8    pi,radius,tareas,bmoismfix
        real*8    tzero,nethfxland(nlat,nlon),landheat(nlat,nlon)
        real*8    meltheat(nlat,nlon),evapl(nlat,nlon),evapoc(nlat,nlon)
        real*8    dsnow(nlat,nlon),dsnm,tland(nlat,nlon)
        integer   ilabas(nlat,nlon),iocbas(nlat,nlon)
        real*8    runofl(nlat,nlon),runofo(nlat,nlon),arocbas(mbasins)
        real*8    albsnow(nlat),albland(nlat,nlon,4)
        real*8    alblandR(nlat,nlon,4),forestfrR(nlat,nlon),alblbmR(nlat,nlon)
        real*8    alblandismR(nlat,nlon,4,3)
        real*8    forestfr(nlat,nlon),fractl(nlat,nlon)
        real*8    rlatfus,rlatsub,rlatvap,rowat,dareas(nlat)
        real*8    rainf(nlat,nlon),snowf(nlat,nlon),alblbm(nlat,nlon)
        real*8    heatsnown,heatsnows
        real*8    runo_yn,runo_ys

        common /lbmbmois/ bmoisg,bmoism,rainf,snowf
        common /lbmheat/ tland,lhcap,rlhcap,dtland,rdtland,tzero
        common /lbmflux/ nethfxland,landheat,meltheat,evapl
        common /lbmwater/rlatfus,rlatsub,rlatvap,rowat
        common /lbmcsnow/ dsnow,dsnm
        common /lbmiruno/ ilabas,iocbas,nbasins
        common /lbmrunol/ runofl,runofo,arocbas,runo_yn,runo_ys
        common /lbmcalbedo/albsnow,albland,forestfr,alblbm
        common /lbmcalbedo2/alblandR,forestfrR,alblbmR,alblandismR
        common /lbmgrid/fractl,dareas,pi,radius,tareas
        common /lbmhsnow/ heatsnown,heatsnows
