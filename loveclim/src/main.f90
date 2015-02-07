!23456789012345678901234567890123456789012345678901234567890123456789012
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
      PROGRAM emic
      IMPLICIT NONE
      INCLUDE 'comatm.h'
      INCLUDE 'comsurf.h'
      INCLUDE 'comemic.h'
      DOUBLE PRECISION patmCO2
      INTEGER istep

      PCO2ref = 277.4D0
      PGACO2 = PCO2ref

      CALL initdriver
      CALL initecbilt
      CALL initlbm
      CALL inioceanfixed
      CALL initcoup

      initialization = .FALSE.
      patmCO2 = PGACO2

      DO istep = 1, ntstep
         IF (MOD(istep - 1, iatm) == 0) THEN
            WRITE (*,100) iyear, imonth, iday
         END IF
100      FORMAT (I6.6, '/', I2.2, '/', I2.2)

         CALL update(istep)
         CALL co2at
         CALL at2co
         CALL fluxes(noc)
         CALL fluxes(nse)
         ! *** integrate atmosphere
         CALL ecbilt(istep)

         CALL la2co
         CALL fluxes(nld)
         CALL co2la
         CALL lbm(istep)
         CALL lae2co
         CALL sumfluxland

         CALL sumfluxocean(istep)

         CALL oceanfixed(INT((day + 0.5 * dt) / (iatm * dt)) + 1)

         CALL atmout(istep)
         CALL writestate(istep)
      END DO

      CALL writestate(ntotday)
      CALL error(999)

      END PROGRAM emic

!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** initialisation of climate model
!-----------------------------------------------------------------------
      SUBROUTINE initdriver
      IMPLICIT NONE
      INCLUDE 'comatm.h'
      INCLUDE 'comemic.h'
      INCLUDE 'comunit.h'

      INTEGER, PARAMETER :: ijatm = nlat * nlon
      INTEGER, PARAMETER :: ismfile = 400
      INTEGER ija, i, j
      REAL*8 fractocn(ijatm)
      CHARACTER*3 num_startday

      NAMELIST /tstepctl/ nyears,ndays,irunlabel,irunlabeld,iatm, &
           & nwrskip,nwrskip_days

      INCLUDE 'openfiles.h'

! *** open emic.param
      lradCO2 = .TRUE.
      lferCO2 = .TRUE.

! *** open namelist
      nyears = 10
      ndays = 0
      irunlabel = 000000
      irunlabeld = 0
      iatm = 6
      nwrskip = 50
      nwrskip_days =0
      READ(iuo+46, NML=tstepctl)

      IF (irunlabeld < 360) THEN
         WRITE(num_startyear,'(i6.6)') irunlabel
         WRITE(num_startday,'(i3.3)') irunlabeld + 1
      ELSE
         WRITE(num_startyear,'(i6.6)') irunlabel + 1
         WRITE(num_startday,'(i3.3)') 1
      END IF

      fini = num_startyear // '_' // num_startday
      OPEN(iuo+20, FILE='book'//fini, FORM='formatted')
      OPEN(iuo+99, FILE='info'//fini, FORM='formatted')
      OPEN(iuo+29, FILE='error'//fini, FORM='formatted')

      kism = 1
      if_ism = 15
      is_ism = 15

      WRITE(iuo+30,900) 'nyears       =', nyears
      WRITE(iuo+30,900) 'ndays        =', ndays
      WRITE(iuo+30,900) 'irunlabel    =', irunlabel
      WRITE(iuo+30,900) 'irunlabeld   =', irunlabeld
      WRITE(iuo+30,900) 'iatm         =', iatm
      WRITE(iuo+30,900) 'nwrskip      =', nwrskip
      WRITE(iuo+30,900) 'nwrskip_days =', nwrskip_days

      undef = 9.99E10

! *** nstpyear is number of atmospheric timesteps per year
! *** ntstep is total number of timesteps
      nstpyear = iatm * 360
      ntstep = nstpyear * nyears
      ntotday = nyears * 360 + ndays

      READ(iuo+48,*)
      READ(iuo+48,*) (fractocn(ija), ija=1, ijatm)
      REWIND(iuo+48)
      DO ija = 1, ijatm
         j = INT((ija - 1) / nlon) + 1
         i = ija - (j - 1) * nlon
         fracto(j,i) = fractocn(ija)
         IF (fracto(j,i) > 0.990) fracto(j,i) = 1.0d0
      END DO

      DO i = 1, nlat
         READ(iuo+50,*) dareafac(i)
      END DO

900   FORMAT(a14,1x,i6)

      RETURN
      END SUBROUTINE initdriver

!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------------
! ***
! *** This routine initialises the day, month, year of the model run
! ***
!-----------------------------------------------------------------------------
      SUBROUTINE inimdldate
      IMPLICIT NONE
      INCLUDE 'comatm.h'
      INCLUDE 'comemic.h'
      INCLUDE 'comunit.h'
      INCLUDE 'comrunlabel.h'

      day     = 0
      iyear   = 0
      imonth  = INT(((irunlabeld - 1) - MOD(irunlabeld - 1, 30)) / 30) + 1
      iday    = MOD(irunlabeld - 1, 30) + 1
      iseason = MOD(irunlabeld, 90)
      initialization = .TRUE.

      WRITE (iuo+99,*) 'Init Date', irunlabelf + iyear, imonth, iday

      RETURN
      END SUBROUTINE inimdldate


!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------------
! ***
! *** This routine calculates the day, month, year of the model run from
! *** integration steps.
! *** Written by xueli wang April 1995.
! ***
!-----------------------------------------------------------------------------
      SUBROUTINE mdldate(istep)
      IMPLICIT NONE
      INCLUDE 'comatm.h'
      INCLUDE 'comdyn.h'
      INCLUDE 'comemic.h'
      INCLUDE 'comunit.h'
      INCLUDE 'comrunlabel.h'

      INTEGER iy, im, idate, istep

      day = (MOD(irunlabeld * iatm + istep - 1, nstpyear)) * dt

      IF (MOD(istep - 1, iatm) == 0) THEN
         iday = iday + 1
         IF (iday > 30) THEN
            iday = 1
            imonth = imonth + 1
            IF (imonth > 12) THEN
               imonth = 1
               iyear = iyear + 1
            END IF
         END IF
      END IF

      iy = iyear * 10000
      im = imonth * 100
      idate = 20000000 + iy + im +iday

      RETURN
      END SUBROUTINE mdldate

!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
!*** this routine writes the current state of ecbilt
!    to datafiles for each day
!-----------------------------------------------------------------------
       SUBROUTINE writestate(istep)
       IMPLICIT NONE
       INCLUDE 'comatm.h'
       INCLUDE 'comemic.h'
       INCLUDE 'comunit.h'

       INTEGER kday, kyear, kInDays, istep, nwrskip_totdays
       CHARACTER*6 numyear
       CHARACTER*3 numday

       IF (MOD(istep,nstpyear) == 0) THEN
          IF (MOD(iyear,nwrskip) == 0 .OR. iyear == nyears) THEN
             nwrskip_totdays = nwrskip * 360 + nwrskip_days

             kInDays = irunlabeld + istep + irunlabel * 360
             kday = MOD(kInDays - 1, 360) + 1
             kyear = (kInDays - kday) / 360

             WRITE(numyear,'(i6.6)') kyear
             WRITE(numday,'(i3.3)') kday
             OPEN(iuo+95, FILE='startdata/inatdyn'//numyear//'_'// &
                  & numday//'.dat', FORM='unformatted')
             CALL wrenddyn
             CLOSE(iuo+95)
             OPEN(iuo+95, FILE='startdata/inatphy'//numyear//'_'// &
                  & numday//'.dat', FORM='unformatted')
             CALL wrendphy
             CLOSE(iuo+95)
             OPEN(iuo+95, FILE='startdata/inland'//numyear//'_'// &
                  & numday//'.dat', FORM='unformatted')
             CALL wrendland
             CLOSE(iuo+95)
             OPEN(iuo+95, FILE='startdata/incoup'//numyear//'_'// &
                  & numday//'.dat', FORM='unformatted')
             CALL wrendcoup
             CLOSE(iuo+95)
          END IF
       END IF

       RETURN
       END SUBROUTINE writestate
