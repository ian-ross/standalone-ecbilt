c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File: comlakeout.h contains common declarations for lake output
c *** routines in file oceandiag0.f
c-----------------------------------------------------------------------

      integer   nlakes,ioutml,ioutfl
      parameter (nlakes=8)
      real*8    arlake(nlakes)
      real*8    temlt(nzl,nlakes),comlt(nzl,nlakes)
      real*8    samlt(nzl,nlakes)
      real*8    tmefxlt(nlakes),tmhfxlt(nlakes),tmdlslt(nlakes)
      real*8    tmulslt(nlakes),tmsslt(nlakes),tmhelt(nlakes)
      real*8    tmpplt(nlakes),tmrunlt(nlakes),tmbrilt(nlakes)
      real*8    tmsalt(nlakes),tmicelt(nlakes),taicelt(nlakes)



      common /comoutlmr/  temlt,samlt,comlt
      common /comoutlmi/  ioutml,ioutfl
      common /comoutlbar/ arlake
      common /comoutlfr/  tmefxlt,tmhfxlt,tmdlslt,tmulslt,tmsslt,
     &                    tmhelt,tmpplt,tmrunlt,tmbrilt,tmsalt,
     &                    tmicelt,taicelt
