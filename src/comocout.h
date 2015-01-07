c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File: comocout.h contains common declarations for ocean output
c *** routines in file oceandiag0.f
c *** rein addition of variables for oamflux
c-----------------------------------------------------------------------

      integer   noutbas,ioutmst,ioutfst,ioutzmst,ioutto
      integer   ioutc,iouterm,ioutlakm,ioutova
      parameter (noutbas=9)
      integer   nmh,nlh,nzh
      parameter (nmh=31,nlh=63,nzh=12)
      integer   iocvbas(0:nlh,0:nmh+1),ioctbas(0:nlh,0:nmh)
      integer   izatbas(0:nlh,0:nmh)


      real*8    arocvbas(noutbas),aroctbas(noutbas)
      real*8    temeant(nzh,noutbas),comeant(nzh,noutbas)
      real*8    sameant(nzh,noutbas)
      real*8    tmefxt(noutbas),tmhfxt(noutbas),tmdlst(noutbas)
      real*8    tmulst(noutbas),tmsst(noutbas),tmhet(noutbas)
      real*8    tmppt(noutbas),tmrunt(noutbas),tmbrit(noutbas)
      real*8    tmsat(noutbas),tmicet(noutbas),taicet(noutbas)

      real*8    ovmaxt(2),vtemsut(2)


      common /comoutmr/  temeant,sameant,comeant
      common /comoutmi/  ioutmst,ioutfst,ioutzmst,ioutto
      common /comoutbai/ iocvbas,ioctbas,izatbas
      common /comoutbar/ arocvbas,aroctbas
      common /comoutfr/  tmefxt,tmhfxt,tmdlst,tmulst,tmsst,tmhet,
     &                   tmppt,tmrunt,tmbrit,tmsat,tmicet,taicet

      common /comer/ ovmaxt,vtemsut
      common /icoupop/ ioutc
      common /ioutte/ iouterm,ioutova
      common /ioutla/ ioutlakm
