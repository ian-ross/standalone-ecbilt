












c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comclio.h                                                    
c *** Contents: Common declarations for diagnostics of clio
c-----------------------------------------------------------------------
      integer numtotvaro,maxmrecs,maxarecs
      real*8  fill_valo,missing_valo,newtotvaro(40,20)
      character*60 dimtotvaro(40,6)
      character*80 nametotvaro(40,6)

      common /ec_netcdf_o/ newtotvaro,fill_valo,missing_valo,
     *                     numtotvaro,nametotvaro,dimtotvaro,
     *                     maxmrecs,maxarecs

      logical     initializo,flagmon,flagyear
      character*120 globalatto(26,2)

      common /globalo/ globalatto,initializo,flagmon,flagyear
