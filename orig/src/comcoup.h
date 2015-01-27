c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comcoup.h                                                    
c *** Contents: Common declarations for coupling module of ECbilt
c-----------------------------------------------------------------------

      integer icoupleh,icouples,icouplew,ndayws  
      real*8  samix(nlat,nlon)
      real*8  hefx(nlat,nlon),stmix(nlat,nlon),safx(nlat,nlon)
      real*8  winstu(nlat,nlon),winstv(nlat,nlon),hefxo(nlat,nlon)
      real*8  winstua(nlat,nlon),winstva(nlat,nlon)
      real*8  safxo(nlat,nlon),winstuo(nlat,nlon),winstvo(nlat,nlon)
      real*8  urun(nlat,nlon),vrun(nlat,nlon)
      real*8  corsafx,yearsum
        
      common /rcoup/ hefx,stmix,samix,safx,winstu,winstv,hefxo,safxo,
     &      winstuo,winstvo,corsafx,yearsum,winstua,winstva,urun,vrun
      common /icoup/ icoupleh,icouples,icouplew,ndayws

