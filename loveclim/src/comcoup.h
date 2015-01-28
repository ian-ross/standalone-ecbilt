












c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comcoup.h                                                    
c *** Contents: Common declarations for coupling module of ECbilt
c-----------------------------------------------------------------------

      integer    kamax,komax,ijatm,ijocn
      parameter  (ijatm=64*32,ijocn=122*65,komax=17,kamax=14)
      integer iobclint,iobtropt,indo2a(ijatm,kamax),inda2o(ijocn,komax)
      real*8  cohesws(nlat,nlon),clhesws(nlat,nlon)
      real*8  cohesw0(nlat,nlon),clhesw0(nlat,nlon)
      real*8  cohesw1(nlat,nlon),clhesw1(nlat,nlon)
      real*8  cohesw2(nlat,nlon),clhesw2(nlat,nlon)
      real*8  coulrads(nlat,nlon),clulrads(nlat,nlon)
      real*8  coulrad0(nlat,nlon),clulrad0(nlat,nlon)
      real*8  coulrad1(nlat,nlon),clulrad1(nlat,nlon)
      real*8  coulrad2(nlat,nlon),clulrad2(nlat,nlon)
      real*8  codlrads(nlat,nlon),cldlrads(nlat,nlon)
      real*8  coeflux(nlat,nlon),cleflux(nlat,nlon)
      real*8  cohflux(nlat,nlon),clhflux(nlat,nlon)
      real*8  coevap(nlat,nlon),clevap(nlat,nlon)
      real*8  winstua(nlat,nlon),winstva(nlat,nlon)
      real*8  sumohfx(nlat,nlon),sumoswr(nlat,nlon)
      real*8  sumihfx(nlat,nlon),sumiswr(nlat,nlon) 
      real*8  sumofwf(nlat,nlon),sumisno(nlat,nlon) 
      real*8  sumty(nlat,nlon),sumtx(nlat,nlon)
      real*8  sumuv10(nlat,nlon)
      real*8  couptcc(nlat,nlon)
      real*8  couprf(nlat,nlon),coupsf(nlat,nlon)
      real*8  samix(nlat,nlon)
      real*8  sumrl(nlat,nlon),sumro(nlat,nlon)
      real*8  sumhsn,sumhss,sumohsn,sumohss
      real*8  couprunl(nlat,nlon),coupruno(nlat,nlon)
      real*8  couphsnn,couphsns
      real*8  sumicof(nlat,nlon),sumpress(nlat,nlon)
      real*8  wo2a(ijatm,kamax),wa2o(ijocn,komax)

      common /ec_rcoup/winstua,winstva,sumtx,sumty,
     &      couprf,coupsf,couptcc,samix,
     &      cohesws,clhesws,cohesw0,clhesw0,cohesw1,clhesw1,
     &      cohesw2,clhesw2,coulrads,clulrads,coulrad0,clulrad0,
     &      coulrad1,clulrad1,coulrad2,clulrad2,codlrads,cldlrads,
     &      coeflux,cleflux,cohflux,clhflux,coevap,clevap,
     &      sumrl,sumro,couprunl,coupruno,couphsnn,couphsns,
     &      sumohfx,sumoswr,sumihfx,sumiswr,sumofwf,sumisno,
     &      sumhsn,sumhss,sumohsn,sumohss,sumicof,sumuv10,sumpress
     
      common /ec_icoup/ iobtropt,iobclint,indo2a,inda2o
      
      common /ec_gridinfo/ wo2a,wa2o

