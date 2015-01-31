!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** File:     comemic.h
! *** Contents: Common declarations global control variables of EMIC
!-----------------------------------------------------------------------
! *** COMMON /timectl/ day,ntstep,nyears,nstpyear,iatm,
!                      iyear,imonth,iday
!     day:       keeps track of time in days during the integration
!     ntstep:    total number of timesteps of the integration
!     nyears:    total integration period in years
!     ndays:     total integration period in days mod 360
!     => total integration period in days = nyears*360 + ndays
!     nstpyear:  number of atmospheric timesteps that fit in one year
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
! *** COMMON /coupctl/ idtbclin,idtbtrop
!     idtbclin: baroclinic timestep of ocean in days
!     idtbtrop: barotropic timestep of ocean in days
!
!     lferCO2: fertilization effect of CO2 if T
!     lradCO2: radiative forcing of CO2 if T
!
!-----------------------------------------------------------------------

      real*8      day,undef,fracto(nlat,nlon),dareafac(nlat)
      real*8      PCO2ref,PGACO2
      integer     nyears,ndays,irunlabel,irunlabeld,iatm, &
           & nwrskip,nwrskip_days,end_year,end_day
      integer     ntstep,nstpyear,isatfor,ntotday
      integer     iyear,imonth,iday,iseason
      integer     is_ism,if_ism,kism,tstartism
      character*10 fini
      logical     flgtsi,flgvol,flgghg,flgsul
      logical     lferCO2,lradCO2
      logical     initialization
      character*6 num_startyear

      common /timepar/ nyears,ndays,irunlabel,irunlabeld,iatm, &
           & nwrskip,nwrskip_days,end_year,end_day,num_startyear
      common /timectl/ day,ntstep,nstpyear,iyear,imonth, &
           & iday,iseason,ntotday,initialization
      common /startctl/fini

      common /globaldef/undef,fracto,dareafac

      common /couplint/ is_ism,if_ism,kism,tstartism
      common /trans/ flgtsi,flgvol,flgghg,flgsul

      common /dioxid/lferCO2,lradCO2

      common /lo2atm/ PGACO2,PCO2ref
