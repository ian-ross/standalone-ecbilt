        integer ivlevel(nvl),itlevel(0:nvl)
        real*8  tsurf1(nlat,nlon),temp4g1(nlat,nlon),temp2g1(nlat,nlon),
     *          tempsg1(nlat,nlon),albep(nlat,nlon),temp0g1(nlat,nlon)
        real*8  dyrain1(nlat,nlon),corain1(nlat,nlon),torain1(nlat,nlon),
     *          evap1(nlat,nlon),eminp1(nlat,nlon),hesw(nlat,nlon),
     *          nlrads(nlat,nlon),runofl1(nlat,nlon),runofo1(nlat,nlon),
     *          winstu1(nlat,nlon),winstv1(nlat,nlon),snow1(nlat,nlon)

        common /ec_local1/  ivlevel,itlevel
        common /ec_local2/  tsurf1, temp4g1,temp2g1,albep, temp0g1,
     *         tempsg1,dyrain1,corain1,torain1,snow1,
     *         evap1,  eminp1, hesw,   nlrads,
     *         runofl1, runofo1,winstu1,winstv1







