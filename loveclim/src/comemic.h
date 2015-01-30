!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** File:     comemic.h
! *** Contents: Common declarations global control variables of EMIC
!-----------------------------------------------------------------------
! *** COMMON /timectl/ day,ntstep,nyears,nstpyear,nocstpyear,iatm,
!                      iyear,imonth,iday
!     day:       keeps track of time in days during the integration
!     ntstep:    total number of timesteps of the integration
!     nyears:    total integration period in years
!     ndays:     total integration period in days mod 360
!     => total integration period in days = nyears*360 + ndays
!     nstpyear:  number of atmospheric timesteps that fit in one year
!     nocstpyear:number of ocean timesteps that fit in one year
!     iatm:      number of atmospheric timesteps that fit in one day
!     nwrskip    number of years between writing the model state to disk
!     nwrskip_days: number of days mod 360 between writing model state to disk
!     => model state writen every nwrskip*360 + nwrskip_days days
!     iyear:     counts the years during the integration
!     imonth:    counts the months during the integration
!     iday:      counts the days during the integration
!     iseason:   counts the seasons during the integration
!     end_year   end year of the run
!     end_day    end day of the run
!
! *** COMMON /startctl/ irunlabel,fini,fend
!     irunlabel:index of startfile containing initial conditions
!               if zero, then initial conditions are prescribed
!               read in iatmdyn,iatmphys, iniland,
!               iniocean, iniseaice
!     irunlabeld: day index of startfile
!     fini:     character string which contains irunlabel
!     fend:     character string which contains irunlabel+nyears
!
! *** COMMON /coupctl/ idtbclin,idtbtrop,ibclins,nbclins,ibtrops,nbtrops
!     idtbclin: baroclinic timestep of ocean in days
!     idtbtrop: barotropic timestep of ocean in days
!     ibclins:  counts the number of atmospheric steps in one baroclinic
!               ocean step
!     nbclins:  total number of atmospheric timesteps in one barocinic
!               ocean step
!     ibtrops:  counts the number of atmospheric steps in one barotropic
!               ocean step
!     nbtrops:  total number of atmospheric timesteps in one barotropic
!               ocean step
!
!     lferCO2: fertilization effect of CO2 if T
!     lradCO2: radiative forcing of CO2 if T
!
!-----------------------------------------------------------------------

      real*8      day,undef,fracto(nlat,nlon),dareafac(nlat)
      real*8      PCO2ref,PGACO2
      integer     nyears,ndays,irunlabel,irunlabeld,iatm, &
           & iobtrop,iobclin,nwrskip,nwrskip_days,nbclins,nbtrops, &
           & end_year,end_day
      integer     ntstep,nstpyear,nocstpyear,isatfor,ntotday
      integer     iyear,imonth,iday,iseason
      integer     is_ism,if_ism,kism,tstartism
      character*10 fini
      logical     flgveg
      logical     flgtsi,flgvol,flgghg,flgsul
      logical     lferCO2,lradCO2
      logical     initialization
      character*120 globalatt(26,2)

      common /timepar/ nyears,ndays,irunlabel,irunlabeld,iatm, &
           & iobtrop,iobclin,nwrskip,nwrskip_days,end_year,end_day
      common /timectl/ day,ntstep,nstpyear,nocstpyear,iyear,imonth, &
           & iday,iseason,ntotday,nbclins,nbtrops,initialization
      common /startctl/fini

      common /globaldef/undef,fracto,dareafac

      common /coupl/ flgveg

      common /couplint/ is_ism,if_ism,kism,tstartism
      common /trans/ flgtsi,flgvol,flgghg,flgsul

      common /dioxid/lferCO2,lradCO2

      common /globala/ globalatt

      common /lo2atm/ PGACO2,PCO2ref
