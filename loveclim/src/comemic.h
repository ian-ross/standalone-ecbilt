












c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comemic.h                                                    
c *** Contents: Common declarations global control variables of EMIC
c-----------------------------------------------------------------------
c *** COMMON /ec_timectl/ day,ntstep,nyears,nstpyear,nocstpyear,iatm,
c                      iyear,imonth,iday
c     day:       keeps track of time in days during the integration 
c     ntstep:    total number of timesteps of the integration
c     nyears:    total integration period in years
c     ndays:     total integration period in days mod 360
c     => total integration period in days = nyears*360 + ndays
c     nstpyear:  number of atmospheric timesteps that fit in one year
c     nocstpyear:number of ocean timesteps that fit in one year
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
c *** COMMON /ec_coupctl/ idtbclin,idtbtrop,ibclins,nbclins,ibtrops,nbtrops
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
c     lferCO2: fertilization effect of CO2 if T
c     lradCO2: radiative forcing of CO2 if T
c
c-----------------------------------------------------------------------

      real*8      day,undef,fracto(nlat,nlon),dareafac(nlat)
      real*8      PCO2ref,PGACO2
      integer     nyears,ndays,irunlabel,irunlabeld,iatm,ilan,iice,iobtrop,iobclin,
     *            nwrskip,nwrskip_days,nbclins,nbtrops,end_year,end_day
      integer     ntstep,nstpyear,nocstpyear,isatfor,ntotday
      integer     iyear,imonth,iday,iseason
      integer     is_ism,if_ism,kism,tstartism,flgloch
      character*10 fini
      logical     flgveg,flgicb,flgisma,flgismg
      logical     flgtsi,flgvol,flgghg,flgsul
      logical     lferCO2,lradCO2
      logical     initialization
      character*120 globalatt(26,2)

      common /ec_timepar/ nyears,ndays,irunlabel,irunlabeld,iatm,ilan,iice,iobtrop,
     *                 iobclin,nwrskip,nwrskip_days,end_year,end_day
      common /ec_timectl/ day,ntstep,nstpyear,nocstpyear,iyear,imonth,
     *                    iday,iseason,ntotday,nbclins,nbtrops,initialization
      common /ec_startctl/fini

      common /ec_globaldef/undef,fracto,dareafac

      common /ec_coupl/ flgveg,flgicb,flgisma,flgismg

      common /ec_couplint/ is_ism,if_ism,kism,tstartism
      common /ec_trans/ flgtsi,flgvol,flgghg,flgsul

      common /ec_dioxid/lferCO2,lradCO2

      common /globala/ globalatt

      common /lo2atm/ PGACO2,PCO2ref
