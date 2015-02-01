!23456789012345678901234567890123456789012345678901234567890123456789012
      program emic
!-----------------------------------------------------------------------
! *** coupled coupled atmosphere-ocean-seaice-land-carbon model
! *** atmosphere      : ecbilt
! *** ocean-seaice    : clio
! *** land            : lbm (land bucket model)
! *** carbon          - ocean & atmosphere : loch
!                     - continents         : vecode
! *** Ice-sheets      : agism
!-----------------------------------------------------------------------
!
! CO2 atmospheric concentration:
!
!  ECBilt reads PGACO2 in a forcing file
!
!-----------------------------------------------------------------------
! ... The reference value PCO2ref is read by EcBilt
!                                 is modified only in case of equilibrium run
!                                                            (iscenghg=0 or 3)
!     It allows sensitivity runs with/without CO2 radiative and/or
!     fertilisation effects
!-----------------------------------------------------------------------
! ...... Vecode takes patmCO2
!
!-----------------------------------------------------------------------


      implicit none

      include 'comatm.h'
      include 'comsurf.h'
      include 'comemic.h'
      integer ittt
      double precision patmCO2
      integer istep,k


      PCO2ref=277.4D0
      PGACO2=PCO2ref

      call initdriver
      call initecbilt
      call inioceanfixed
      call initlbm
      call initcoup

      initialization=.false.
      patmCO2=PGACO2

      do istep = 1, ntstep
         call update(istep)
         call co2at
         call at2co
         call fluxes(noc)
         call fluxes(nse)
         ! *** integrate atmosphere
         call ecbilt(istep)

         call la2co
         call fluxes(nld)
         call co2la
         call lbm(istep)
         call lae2co
         call sumfluxland

         call sumfluxocean(istep)

         call oceanfixed(int((day+0.5*dt)/(iatm*dt)) + 1)

         !!! CHECK!
         call writestate(iday)
      enddo

      call writestate(ntotday)
      call error(999)
!-AM

      close(57)
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine initdriver
!-----------------------------------------------------------------------
! *** initialisation of climate model
!-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comemic.h'
      include 'comunit.h'

      integer   ijatm,ija,i,j,k,ismfile
      parameter (ijatm=nlat*nlon)
      real*8    fractocn(ijatm)
      parameter (ismfile = 400)
      character*3 num_startday
      logical dummy

      NAMELIST /tstepctl/nyears,ndays,irunlabel,irunlabeld,iatm, &
           & nwrskip,nwrskip_days

      include 'openfiles.h'

      write(iuo+99,*) 'Initialize'

! *** open emic.param

      lradCO2 = .TRUE.
      lferCO2 = .TRUE.

! *** open namelist

      nyears=10
      ndays=0
      irunlabel=000000
      irunlabeld=0
      iatm=6
      nwrskip=50
      nwrskip_days=0

      read(iuo+46, NML = tstepctl)

      if(irunlabeld.lt.360) then
        write(num_startyear,'(i6.6)')irunlabel
        write(num_startday,'(i3.3)')irunlabeld+1
      else
        write(num_startyear,'(i6.6)')irunlabel+1
        write(num_startday,'(i3.3)')1
      endif

      fini=num_startyear//'_'//num_startday

      kism=1
      if_ism=15
      is_ism=15

      write(iuo+30, 900) 'nyears       =', nyears
      write(iuo+30, 900) 'ndays        =', ndays
      write(iuo+30, 900) 'irunlabel    =', irunlabel
      write(iuo+30, 900) 'irunlabeld   =', irunlabeld
      write(iuo+30, 900) 'iatm         =', iatm
      write(iuo+30, 900) 'nwrskip      =', nwrskip
      write(iuo+30, 900) 'nwrskip_days =', nwrskip_days

      undef = 9.99E10

! *** nstpyear is number of atmospheric timesteps per year
! *** ntstep is total number of timesteps

      nstpyear   = iatm*360
      ntstep     = nstpyear*nyears
      ntotday    = nyears*360+ndays

      read(iuo+48,*)
      read(iuo+48,*) (fractocn(ija),ija=1,ijatm)
      rewind(iuo+48)
      do ija=1,ijatm
        j=int((ija-1)/nlon)+1
        i=ija-(j-1)*nlon
        fracto(j,i)=fractocn(ija)
        if (fracto(j,i).gt.0.990) fracto(j,i)=1.0d0
      enddo


      do i=1,nlat
        read(iuo+50,*) dareafac(i)
      enddo

900   format(a14,1x,i6)
901   format(a14,1x,l6)

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine inimdldate
!-----------------------------------------------------------------------------
! ***
! *** This routine initialises the day, month, year of the model run
! ***
!-----------------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comrunlabel.h'

      day     = 0
      iyear   = 0
      imonth  = int(((irunlabeld - 1) - mod(irunlabeld - 1, 30)) / 30) + 1
      iday    = mod(irunlabeld - 1, 30) + 1
      iseason = mod(irunlabeld, 90)
      initialization = .true.

      write (iuo+99,*) 'Init Date', irunlabelf + iyear, imonth, iday

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine mdldate(istep)
!-----------------------------------------------------------------------------
! ***
! *** This routine calculates the day, month, year of the model run from
! *** integration steps.
! *** Written by xueli wang April 1995.
! ***
!-----------------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comdyn.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comrunlabel.h'

      integer iy,im,idate,istep

      day = (mod(irunlabeld*iatm+istep-1,nstpyear)) * dt

      if (mod(istep-1,iatm).eq.0) then
        iday = iday + 1
        if (iday .gt. 30) then
          iday = 1
          imonth = imonth + 1
          if (imonth.gt.12) then
            imonth = 1
            iyear = iyear + 1
          endif
          call progress(irunlabelf+iyear,imonth-1,12-1)
        endif
!        write(iuo+99,'(A,I,A1,I2,A1,I2)') '>>>Update Date', irunlabelf+iyear,"/",imonth,"/", iday
      endif

      iy = iyear * 10000
      im = imonth * 100

      idate = 20000000 + iy + im +iday

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012

       subroutine writestate(ist)
!-----------------------------------------------------------------------
!*** this routine writes the current state of ecbilt
!    to datafiles for each day
!-----------------------------------------------------------------------
       implicit none

       include 'comatm.h'
       include 'comemic.h'
       include 'comunit.h'

       integer kday
       integer kyear
       integer kInDays
       integer ist
       integer nwrskip_totdays
       character*6 numyear
       character*3 numday

       nwrskip_totdays = nwrskip*360+nwrskip_days


       if (mod(ist,nwrskip_totdays).eq.0.or.ist.eq.ntotday) then
          write(*,*) " "
          kInDays=irunlabeld+ist+irunlabel*360
          kday=mod(kInDays-1,360)+1
          kyear=(kInDays-kday)/360
!           kday=mod(irunlabeld+ist,360)
!           kyear=irunlabel+int((irunlabeld+ist)/360)

          write(numyear,'(i6.6)') kyear
          write(numday,'(i3.3)') kday
          open(iuo+95,file='startdata/inatdyn'//numyear//'_'//numday//'.dat' &
     &           ,form='unformatted')
          call wrenddyn
          close(iuo+95)
          open(iuo+95,file='startdata/inatphy'//numyear//'_'//numday//'.dat' &
     &           ,form='unformatted')
          call wrendphy
          close(iuo+95)
          open(iuo+95,file='startdata/inland'//numyear//'_'//numday//'.dat' &
     &           ,form='unformatted')
          call wrendland
          close(iuo+95)
          open(iuo+95,file='startdata/incoup'//numyear//'_'//numday//'.dat' &
     &           ,form='unformatted')
          call wrendcoup
          close(iuo+95)

       endif

       return

       end

      subroutine progress(date,ndone,ntotal)
        implicit none
        integer :: date
        character*255 prog,oldprog
        double precision oldtime,hires_time,tl
        integer ndone,ntotal,i
        save oldprog,oldtime

        write(prog,'(I6,''['')') date
        do i=1,40
          prog(7+i:7+i)=' '
        enddo
        write(prog(23:31),'(f7.1,''%'')') 100.0*ndone/ntotal
        do i=1,40
          if ((1.0*ndone/ntotal).gt.(1.0*i/40)) then
            if (prog(7+i:7+i).eq.' ') prog(7+i:7+i)='#'
          endif
        enddo
        prog(47:47)=']'
        if (prog.ne.oldprog) write(0,'(a,a,$)') prog(1:77),char(13)
        oldprog=prog
        if (ndone.eq.ntotal) write(0,*)
        return
      end
