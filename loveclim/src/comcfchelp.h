












c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comcfchelp.h
c *** Contents: common declaration for time (from ecbilt to clio)
c-----------------------------------------------------------------------
c *** COMMON /ec_timectl/ day,ntstep,nyears,nstpyear,nocstpyear,iatm,
c                      iyear,imonth,iday
c     day:       keeps track of time in days during the integration
c     ntstep:    total number of timesteps of the integration
c     nyears:    total integration period in years
c     nstpyear:  number of atmospheric timesteps that fit in one year
c     nocstpyear:number of ocean timesteps that fit in one year
c     iatm:      number of atmospheric timesteps that fit in one day
c     nwrskip    number of years between writing the model state to disk
c     iyear:     counts the years during the integration
c     imonth:    counts the months during the integration
c     iday:      counts the days during the integration
c     iseason:   counts the seasons during the integration
c
c-----------------------------------------------------------------------

      real*8      day
      integer     nyears,irunlabel,iatm,ilan,iice,iobtrop,iobclin,
     *            nwrskip,nbclins,nbtrops
      integer     ntstep,nstpyear,nocstpyear,ntotday
      integer     iyear,imonth,iday,iseason

      common /ec_timepar/ nyears,irunlabel,iatm,ilan,iice,iobtrop,
     *                 iobclin,nwrskip
      common /ec_timectl/ day,ntstep,nstpyear,nocstpyear,iyear,imonth,
     *                    iday,iseason,ntotday,nbclins,nbtrops

