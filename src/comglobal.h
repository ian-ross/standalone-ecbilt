c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comglobal.h                                                    
c *** Contents: Common declarations global control variables of ECbilt
c-----------------------------------------------------------------------
c *** COMMON  /timectl/ day,ntstep,nyears,nstpyear,nocstpyear,iatm,
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
c *** COMMON  /startctl/ irunlabel,fini,fend
c     irunlabel:index of startfile containing initial conditions
c               if zero, then initial conditions are prescribed
c               read in iatmdyn,iatmphys, iniland,
c               iniocean, iniseaice
c     fini:     character string which contains irunlabel
c     fend:     character string which contains irunlabel+nyears
c
c *** COMMON  /coupctl/ idtbclin,idtbtrop,ibclins,nbclins,ibtrops,nbtrops
c     idtbclin: baroclinic timestep of ocean in days
c     idtbtrop: barotropic timestep of ocean in days
c     ibclins:  counts the number of atmospheric steps in one baroclinic
c               ocean step
c     nbclins:  total number of atmospheric timesteps in one barocinic 
c               ocean step
c     ibtrops:  counts the number of atmospheric steps in one barotropic
c               ocean step
c     nbtrops:  total number of atmospheric timesteps in one barotropic 
c               ocean step
c
c *** COMMON  /writectl/ixout,ioutdaily,itel,minterv,
c                      meantype,meantot,meanyl,ifrendat,instcount
c     ixout:   output frequency for instantanous field in days
c     ioutdaily: output (1) or not (0) instantanous fields.
c     itel:    counts the number of days in the output interval
c     minterv: counter used to compute monthly or seasonal mean.
c     meantype = 1 monthly mean
c     meantype = 2 seasonal mean.
c     meantot = 1 output total integration period monthly or seasonal
c               mean
c     meanyl  = 1 output yearly monthly or seasonal mean.
c     ifrendat: frequency of writing restart data in days.
c     instcount: counter used for output instantaneous fields.
c-----------------------------------------------------------------------

      real*8      day,undef
      integer     ntstep,nyears,nstpyear,nocstpyear,iatm,irunatm,isatfor
      integer     nwrskip,iyear,imonth,iday,iseason
      integer     irunlabel
      character*4 fini,fend
      integer     ixout,ioutdaily,itel,minterv,
     *            meantype,meantot,meanyl,ifrendat,instcount 
      integer     idtbclin,idtbtrop,ibclins,nbclins,ibtrops,nbtrops

      common /timectl/ day,ntstep,nyears,nstpyear,nocstpyear,iatm,
     *                 nwrskip,iyear,imonth,iday,iseason
      common /startctl/irunlabel,fini,fend
      common /coupctl/ idtbclin,idtbtrop,ibclins,nbclins,ibtrops,
     *                 nbtrops,irunatm,isatfor
      common /writectl/undef,ixout,ioutdaily,itel,minterv,
     *                 meantype,meantot,meanyl,ifrendat,instcount                







