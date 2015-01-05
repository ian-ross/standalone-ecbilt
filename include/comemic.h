












c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comemic.h
c *** Contents: Common declarations global control variables of EMIC
c-----------------------------------------------------------------------
c *** COMMON /ec_timectl/ day,ntstep,nyears,nstpyear,iatm,
c                      iyear,imonth,iday
c     day:       keeps track of time in days during the integration
c     ntstep:    total number of timesteps of the integration
c     nyears:    total integration period in years
c     ndays:     total integration period in days mod 360
c     => total integration period in days = nyears*360 + ndays
c     nstpyear:  number of atmospheric timesteps that fit in one year
c     iatm:      number of atmospheric timesteps that fit in one day
c     nwrskip    number of years between writing the model state to disk
c     nwrskip_days: number of days mod 360 between writing the model state to disk
c     => model state writen every nwrskip*360 + nwrskip_days days
c     iyear:     counts the years during the integration
c     imonth:    counts the months during the integration
c     iday:      counts the days during the integration
c     iseason:   counts the seasons during the integration
c     end_year   end year of the run
c     end_day    end day of the run
c
c *** COMMON /ec_startctl/ irunlabel,fini,fend
c     irunlabel:index of startfile containing initial conditions
c               if zero, then initial conditions are prescribed
c               read in iatmdyn,iatmphys, iniland,
c               iniocean, iniseaice
c     irunlabeld: day index of startfile
c     fini:     character string which contains irunlabel
c     fend:     character string which contains irunlabel+nyears
c
c-----------------------------------------------------------------------

      real*8      day,undef,fracto(nlat,nlon),dareafac(nlat)
      real*8      PCO2ref,PGACO2
      integer     nyears,ndays,irunlabel,irunlabeld,iatm,ilan,iice,
     *            nwrskip,nwrskip_days,end_year,end_day
      integer     ntstep,nstpyear,ntotday
      integer     iyear,imonth,iday,iseason
      character*10 fini
      logical     flgtsi,flgvol,flgghg,flgsul
      logical     initialization
      character*120 globalatt(26,2)

      common /ec_timepar/ nyears,ndays,irunlabel,irunlabeld,iatm,ilan,iice,
     *                 nwrskip,nwrskip_days,end_year,end_day
      common /ec_timectl/ day,ntstep,nstpyear,iyear,imonth,
     *                    iday,iseason,ntotday,initialization
      common /ec_startctl/fini

      common /ec_globaldef/undef,fracto,dareafac

      common /ec_trans/ flgtsi,flgvol,flgghg,flgsul

      common /globala/ globalatt

      common /lo2atm/ PGACO2,PCO2ref
