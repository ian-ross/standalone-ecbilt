c23456789012345678901234567890123456789012345678901234567890123456789012
      program emic
c-----------------------------------------------------------------------
! STANDALONE VERSION OF ECBILT ATMOSPHERE
c-----------------------------------------------------------------------

      implicit none

      include 'comatm.h'
      include 'comsurf.h'
      include 'comemic.h'
      include 'para.com'

      double precision patmCO2
      integer i, j

      PCO2ref = 277.4D0
      PGACO2 = PCO2ref

      call ec_initdriver
      call ec_initecbilt
      call ec_initcoup

      initialization = .false.
      patmCO2 = PGACO2

      do i = 1, ntotday
        do j = 1, iatm
          call ec_update(i, j)
          call ec_co2at
          call ec_fluxes(noc)
          call ec_fluxes(nse)
          call ec_ecbilt(i, j)

          if (j .eq. iatm) then
            call ec_writestate(i)
          endif
          call IPCC_output(i, j)
        enddo
      enddo

      call ec_writestate(ntotday)
      call ec_error(999)
c-AM

      close(57)
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_initdriver
c-----------------------------------------------------------------------
c *** initialisation of climate model
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comemic.h'
      include 'comunit.h'

      integer   ijatm,ija,i,j,ismfile
      parameter (ijatm=nlat*nlon)
      real*8    fractocn(ijatm)
      parameter (ismfile = 400)
      character*6 num_startyear
      character*3 num_startday

      NAMELIST /tstepctl/ nyears,ndays,irunlabel,irunlabeld,iatm,
     &                    nwrskip,nwrskip_days


      include 'openatoutfiles.h'

      write(iuo+99,*) 'Initialize'

c *** open namelist

      include 'openemicinfiles.h'
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
      globalatt(19,1)="branch_time"
      globalatt(19,2)=""//num_startyear


      include 'openemicoutfiles.h'

      write(iuo+30, 900) 'nyears       =', nyears
      write(iuo+30, 900) 'ndays        =', ndays
      write(iuo+30, 900) 'irunlabel    =', irunlabel
      write(iuo+30, 900) 'irunlabeld   =', irunlabeld
      write(iuo+30, 900) 'iatm         =', iatm
      write(iuo+30, 900) 'nwrskip      =', nwrskip
      write(iuo+30, 900) 'nwrskip_days =', nwrskip_days


      undef = 9.99E10

c *** nstpyear is number of atmospheric timesteps per year
c *** ntstep is total number of timesteps


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

      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_inimdldate
C--------------------------------------------------------------------------------
C ***
C *** This routine initialises the day, month, year of the model run
C ***
C-------------------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comrunlabel.h'

      day    = 0
      iyear  = 0
      imonth = int(((irunlabeld-1)-mod(irunlabeld-1,30))/30)+1
      iday   = mod(irunlabeld-1,30)+1
      iseason= mod(irunlabeld,90)
      initialization=.true.

      write(iuo+99,*) 'Init Date', irunlabelf+iyear, imonth, iday

      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_mdldate(istep)
C--------------------------------------------------------------------------------
C ***
C *** This routine calculates the day, month, year of the model run from
C *** integration steps.
C *** Written by xueli wang April 1995.
C ***
C-------------------------------------------------------------------------------
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


c23456789012345678901234567890123456789012345678901234567890123456789012

       subroutine ec_writestate(ist)
c-----------------------------------------------------------------------
c*** this routine writes the current state of ecbilt
c    to datafiles for each day
c-----------------------------------------------------------------------
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
          open(iuo+95,file='startdata/inatdyn'//numyear//'_'//numday//'.dat'
     *           ,form='unformatted')
          call ec_wrenddyn
          close(iuo+95)
          open(iuo+95,file='startdata/inatphy'//numyear//'_'//numday//'.dat'
     *           ,form='unformatted')
          call ec_wrendphy
          close(iuo+95)
       endif

       return

       end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine IPCC_output(ist,jst)
c-----------------------------------------------------------------------
c *** this routine writes the current state of ecbilt to datafiles
c-----------------------------------------------------------------------
      implicit none

      include 'comatm.h'
      include 'comemic.h'
      include 'comunit.h'
      include 'comrunlabel.h'
      include 'comdiag.h'

      DOUBLE PRECISION CARSTOK
      integer istep,ist,jst
      real*8  tmc_ipcc,ts_ipcc,fr_ipcc,cland_ipcc,coc_ipcc,ctot_ipcc
      real*8  moc_ipcc,th_ipcc
      real*8  tmc0,tmc,tsurfmean,moc,thex,cland,oh_ipcc
      real*8  co2_ipcc
      real*8  ulrad0nT,ulrad1nT,fswutoaGA,fswdtoaGA
      real*8  ulrad0nz(nlat,nlon),ulrad1nz(nlat,nlon)
      real*8  fswutoa0(nlat,nlon),fswdtoa0(nlat,nlon)

      common/IPCC_out/tmc_ipcc,ts_ipcc,fr_ipcc,cland_ipcc,coc_ipcc,
     &      moc_ipcc,th_ipcc,oh_ipcc,co2_ipcc
      common/IPCC_out2/moc,tmc,tmc0,tsurfmean,cland,thex
      common /rad031/ulrad0nz,ulrad1nz,ulrad0nT,ulrad1nT
      common /rad_sul1 /fswutoa0,fswutoaGA,fswdtoa0,fswdtoaGA
      COMMON /SORTIE/CARSTOK(3)

      istep=(ist-1)*iatm+jst

      if((mod(nint(day*real(iatm)),nstpyear).eq.0).or.(initialization.eqv..true.)) then
      !if (mod(istep,nstpyear).eq.1) then
        ts_ipcc=0.
        th_ipcc=0.
        moc_ipcc=0.
        oh_ipcc=0.
        co2_ipcc=PGACO2
      endif
      if (mod(istep,iatm).eq.0) then
        ts_ipcc=ts_ipcc+(tsurfmean/360.)
      endif

      if ((initialization.eqv..false.).and.(nint(day*real(iatm)).ne.0)) then
       if(mod(nint(day*real(iatm)),30*iatm).eq.0) then
       !if (mod(istep,(iatm*30)).eq.0) then
        th_ipcc=th_ipcc+(thex/12.)
        moc_ipcc=moc_ipcc+(moc/12.)
        oh_ipcc=oh_ipcc+((tmc-tmc0)*1.4E+18)
        tmc0=tmc
       endif
      endif

      if(mod(int(day*iatm)+1,nstpyear).eq.0) then
      !if (mod(istep,nstpyear).eq.0) then
      if (irad.ne.1) ulrad1nT=0.
      cland_ipcc=Carstok(3)
      coc_ipcc=Carstok(2)
      ctot_ipcc=Carstok(1)+Carstok(2)+Carstok(3)
      write(iuo+74,'(I4,1X,F8.3,1X,F8.3,1X,E15.5,1X,F9.4,1X,F6.2,
     &               1X,F8.3,1X,3(E17.7,1X))')
     &       iyear+irunlabelf,(ulrad1nT-fswutoaGA+fswdtoaGA),
     &       ts_ipcc,oh_ipcc,th_ipcc,moc_ipcc,co2_ipcc,
     &       ctot_ipcc,cland_ipcc,coc_ipcc
      ulrad1nT=0.
      fswutoaGA=0.
      fswdtoaGA=0.

      endif

      !attention ce test de marchera pas en mensuel
      if ((mod(istep,nstpyear).eq.0).and.(iyear.eq.nyears)) then
        close (iuo+74)
        write(42,*) tmc
        close(36)
      endif


      end

      subroutine progress(date,ndone,ntotal)
        implicit none
        integer :: date
        character*255 prog,oldprog
        integer ndone,ntotal,i
        save oldprog

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
